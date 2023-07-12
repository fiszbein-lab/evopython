import glob
import os
import re

from collections import defaultdict

from Bio.SeqIO import SeqRecord
from Bio.AlignIO import MafIO

from .feature import Feature


BASE_REGEX = re.compile('[ACGTRYSWKMBDHVN]')
COMPLEMENT = {
    'A': "T",
    'C': "G",
    'G': "C",
    'T': "A",
    'R': "Y",  # A|G to T|C
    'Y': "R",  # T|C to A|G
    'S': "S",  # G|C to C|G
    'W': "W",  # A|T to T|A
    'K': "M",  # G|T to C|A
    'M': "K",  # C|A to G|T
    'B': "V",  # C|G|T to G|C|A
    'D': "H",  # A|G|T to T|C|A
    'H': "D",  # T|C|A to A|G|T
    'V': "B",  # G|C|A to C|G|T
    'N': "N",  # A|C|G|T to T|G|C|A
    '-': "-"
}


class MAF:
    def __init__(self, maf_dir: str, aligned_on: str) -> None:
        """Inits MAF.

        Args:
            maf_dir: The path to a directory of MAF files, each following the
                naming scheme "<chromosome>.maf."
            aligned_on: The species that the chromosome names map to.
        """
        self._maf_dir = maf_dir
        self._aligned_on = aligned_on

        self._maf_index = dict()
        self._chrom_size = dict()

        for path in glob.iglob(f"{self._maf_dir}/*.maf"):
            chrom, _ = os.path.basename(path).split(".", maxsplit=1)

            basename = f"{self._maf_dir}/{chrom}"
            self._maf_index[chrom] = MafIO.MafIndex(
                f"{basename}.mafindex",
                f"{basename}.maf",
                f"{self._aligned_on}.{chrom}"
            )

            # I peek at the MAF files to avoid asking for a genome index.
            with open(path, 'r') as f:
                for line in f:
                    if not line.startswith("s"):
                        continue

                    _, source, *_, source_size, _ = line.split()
                    species, _ = source.split(".", maxsplit=1)

                    if species == self._aligned_on:
                        self._chrom_size[chrom] = int(source_size)
                        break

    def get(self, feat: Feature, match_strand: bool = True) -> dict[int: dict]:
        """Finds an alignment for a feature.

        Args:
            feat: The feature to get an alignment for.
            match_strand: A bool expressing whether the alignment should match
                the feature's strand; if False, the alignment is mapped to the
                forward strand.

        Returns:
            A list of dicts mapping species to tuple, where tuple[0] is the
            chromosome name; tuple[1] the 0-based, inclusive starting
            coordinate; tuple[2] the 0-based, exclusive ending coordinate;
            tuple[3] the strand, plus or minus for forward or reverse; and
            tuple[4] the alignment.
        """
        alignments = list()
        search_result = self._search(feat=feat)

        for b, (derived_feat, overlap) in enumerate(search_result):
            alignment = dict()

            on_species_overlap = overlap[self._aligned_on]
            on_chrom, on_start, on_end, on_strand, on_seq = on_species_overlap

            M, N = _resolve_indexes(derived_feat, on_species_overlap)

            # M and N are the indexes of the start and ending bases of the
            # feature on the forward-mapped sequence (it's easier), so we may
            # need to reverse complement each alignment.
            if on_strand == "+":
                rev_comp = False
            else:
                rev_comp = True

            for species in overlap:
                chrom, start, end, strand, seq = overlap[species]

                if rev_comp:
                    seq = _reverse_complement(seq)
                    strand = "+" if strand == "-" else "+"

                if strand == "+":
                    # The number of bases preceding the feature alignment is
                    # our start coordinate adjustment.
                    start += _count_bases(seq[:M])
                else:
                    # The number of bases after the feature alignment, e.g.,
                    # the number of bases preceding the true start coordinate
                    # on the forward strand, is how much we need to shift the
                    # starting coordinate.
                    start += _count_bases(seq[N+1:])

                seq = seq[M:N+1]
                end = start + _count_bases(seq)

                alignment[species] = chrom, start, end, strand, seq

            # The default is to forward-map the alignment, so we update only
            # when both `match_strand` is True and the feature is has a reverse
            # orientation.
            if match_strand and feat.is_reverse:
                for species, seq in alignment.items():
                    alignment[species] = _reverse_complement(seq)

            alignments.append(alignment)

        return alignments

    def _search(self, feat: Feature) -> dict[str: dict] | None:
        """Finds alignments overlapping a feature.

        Most features and reference genomes are mapped to the forward-strand
        by convention, such that we reverse-complement reverse-stranded
        features, where the former operation corrects for the coordinates'
        obligate forward orientation.

        But MAF records on the reverse strand have reverse-stranded
        coordinates [1], and because feature-overlapping alignment blocks can
        be on both strands, we need to check both possible mapping; see [2] for
        the intuition behind between-strand mapping.

        Args:
             feat: The feature to search for.

        Yields:
            A tuple, where tuple[0] is the feature and tuple[1] a dict mapping
            species to alignment; the feature is returned back because in cases
            where the original does not have a contiguous alignment, it's
            adjusted to achieve contiguity.

        Raises:
            KeyError: The feature lies on a chromosome with no matching MAF
                file in the specified directory.
            ValueError: MafIO.MafIndex.search found an overlapping record with
                incorrect offset.

        Refs:
            1. https://genome.ucsc.edu/FAQ/FAQformat.html#format5
            2. http://genomewiki.ucsc.edu/index.php/Coordinate_Transforms
        """
        overlap = defaultdict(dict)

        try:
            maf_index = self._maf_index[feat.chrom]
        except KeyError as err:
            message = (
                f"{feat.locus()} is on {feat.chrom}, but {feat.chrom}.maf was "
                f"not found in {self._maf_dir}.")
            raise Exception(message) from err

        for ref_strand in "+-":
            if ref_strand == "+":
                search_start = feat.start
                search_end = feat.end
            else:
                # Map the coordinates to the reverse strand.
                search_start = self._chrom_size[feat.chrom] - feat.end
                search_end = search_start + len(feat)

            try:
                search = list(maf_index.search([search_start], [search_end]))
            except ValueError as err:
                message = "Found a record with incorrect offset."
                raise Exception(message) from err

            b = 0
            for alignment in search:
                does_overlap_feat = True
                block = dict()

                for record in alignment:
                    species, *pos, chrom_size, seq = _unpack_record(record)
                    chrom, start, end, on_strand, seq_size = pos

                    if species == "ancestral_sequences":
                        continue

                    if (species == self._aligned_on and
                            on_strand != ref_strand):
                        does_overlap_feat = False
                        break

                    if on_strand == "-":
                        # Map the coordinates to the forward strand.
                        start = chrom_size - end
                        end = start + seq_size

                    block[species] = chrom, start, end, on_strand, seq

                if does_overlap_feat:
                    overlap[b] = block
                    b += 1

        N = len(overlap)
        if N > 0:
            # In cases where we get 1 alignment block, the feature must be in-
            # full or in-part inside the block; in the former case, the derived
            # feature is identical to the queried feature, whereas in the
            # latter case, we use the alignment's start or end to generate a
            # contiguous sub-feature. In cases where we get >1 alignment block,
            # the feature is broken up in an analogous manner.
            for b in overlap:
                chrom, start, end, strand, _ = overlap[b][self._aligned_on]

                if b in (0, N-1):
                    if feat.start > start:
                        start = feat.start
                    else:
                        # In these cases, the feature began before the
                        # alignment, so we use the overlapping alignment's
                        # start, i.e., no update is needed.
                        pass

                    if feat.end < end:
                        end = feat.end
                    else:
                        # In these cases, the feature ended after the
                        # alignment, so we use the overlapping alignment's end,
                        # i.e., no update is needed.
                        pass

                yield Feature(chrom, start, end, feat.strand), overlap[b]


def _count_bases(alignment: str) -> int:
    """Counts the number of bases in an alignment.

    Notes:
        The valid bases are the IUPAC codes.

    Args:
        alignment: An aligned sequence.

    Returns:
        The number of bases in the alignment.
    """
    return len(re.findall(BASE_REGEX, alignment))


def _get_base_index(base_number: int, alignment: str) -> int:
    """Finds the index of the Nth base in an alignment.

    Notes:
        The valid bases are the IUPAC codes.

    Args:
        base_number: The base number to index.

    Returns:
        The index of the base at the given number.

    Raises:
        IndexError: The base number was not found in the alignment.
    """
    match_iter = re.finditer(BASE_REGEX, alignment)

    for bases_covered, match in enumerate(match_iter):
        if bases_covered == base_number:
            return match.start()

    raise IndexError(f"{base_number}th base was not found in the alignment.")


def _resolve_indexes(feat: Feature, overlap: tuple) -> tuple[int, int]:
    """Resolves a feature's indexes from within a alignment.

    Args:
        feat: The feature to resolve.
        overlap: A tuple representing an alignment overlapping the feature,
            where tuple[0] is the feature-overlapping alignment's chromosome;
            tuple[1] is the 0-based, inclusive starting coordinate; tuple[2]
            the 0-based, exclusive ending coordinate; tuple[3] the strand, plus
            or minus for forward or reverse; and tuple[4] the alignment.

    Returns:
        A tuple. where tuple[0] is the feature's starting index in the
        alignment and tuple[1] the ending index.
    """
    chrom, start, end, strand, seq = overlap
    if strand == "-":
        seq = _reverse_complement(seq)

    start_base = feat.start - start
    end_base = start_base + len(feat) - 1

    M = _get_base_index(start_base, alignment=seq)
    N = _get_base_index(end_base, alignment=seq)

    return M, N


def _reverse_complement(seq: str) -> str:
    """Reverse complements a sequence.

    Notes:
        The valid bases are the IUPAC codes, as well as dashes.

    Args:
        seq: The sequence to reverse complement.

    Returns:
        The input sequence's reverse complement.
    """
    return "".join([COMPLEMENT[char] for char in seq][::-1])


def _unpack_record(record: SeqRecord) -> tuple:
    """Unpacks needed Bio.SeqIO.SeqRecord information.

    Args:
        record: The Bio.SeqIO.SeqRecord instance to unpack.

    Returns:
        A tuple, where
            tuple[0] is the species name;
            tuple[1] is the chromosome name;
            tuple[2] is the 0-based, inclusive starting coordinate;
            tuple[3] is the 0-based, exclusive ending coordinate;
            tuple[4] is the strnd, plus or minus for forward or reverse;
            tuple[4] is the size of the sequence;
            tuple[5] is the source chromosome size; and
            tuple[6] is the alignment.
    """
    species, chrom = record.id.split(".", maxsplit=1)
    chrom_size = record.annotations['srcSize']

    start = record.annotations['start']
    seq_size = record.annotations['size']
    end = start + seq_size
    strand = "+" if record.annotations['strand'] > 0 else "-"

    seq = str(record.seq).upper()

    return (
        species,
        chrom,
        start,
        end,
        strand,
        seq_size,
        chrom_size,
        seq
    )
