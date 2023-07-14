import collections
import csv
import glob
import os
import re
import unittest

import requests

from evopython import Feature, MAF


class TestResolution(unittest.TestCase):
    comp_wga = dict()
    comp_feat = collections.defaultdict(list)

    alignment_info = dict()
    with open("tests/data/meta_data.csv", 'r') as f:
        reader = csv.reader(f)
        next(reader)

        for row in reader:
            name, *info = row
            alignment_info[name] = info

    for maf_dir in glob.iglob("tests/data/*/*"):
        comp = os.path.basename(maf_dir)

        aligned_on, aligned_against, method = alignment_info[comp]
        wga = MAF(maf_dir, aligned_on)

        with open(f"{maf_dir}/features.txt", 'r') as f:
            for line in f:
                locus = line.strip("\n")

                chrom, start, end, strand = re.split('[:-]', locus)
                strand = "+" if int(strand) == 1 else "-"

                feat = Feature(chrom, int(start), int(end), strand)
                comp_feat[comp].append(feat)

        comp_wga[comp] = {
            'info': (aligned_on, aligned_against, method),
            'wga': wga
        }

    def test_resolution(self):
        for comp in self.comp_wga:
            print(comp)

            on, against, method = self.comp_wga[comp]['info']
            wga = self.comp_wga[comp]['wga']

            for feat in self.comp_feat[comp]:
                our_alignment = wga.get(feat, match_strand=False)
                ens_alignment = _fetch(feat, on, against, method)

                if not ens_alignment:
                    continue

                for block in our_alignment:
                    for species in block:
                        *pos, our_seq = block[species]
                        chrom, start, end, strand = pos

                        if set(our_seq) == {"-"}:
                            # The API won't return species if the whole
                            # sequence is gapped but we do, so the `assertIn`
                            # below fails; these cases don't reflect failed
                            # resolution, so I skip them.
                            continue

                        if species in ens_alignment:
                            locus = f"{chrom}:{start}-{end}:{strand}"
                            self.assertIn(locus, ens_alignment[species])

                            ens_seq = ens_alignment[species][locus]
                            self.assertEqual(our_seq, ens_seq)
                        else:
                            # Podarcis muralis is missing from the
                            # 24:4996761-4996769 alignment in 17 sauropsids
                            # fetched from Ensembl, but I checked it in the
                            # browser, and our resolution is right — so not
                            # sure why the API doesn't retun that species.
                            continue


def _fetch(
    feat: Feature,
    on: str,
    against: str,
    method: str = "LASTZ_NET"
) -> dict[str: dict]:
    """Fetches an alignment from the Ensembl REST API.

    Args:
        feat: The feature to fetch an alignment for.
        on: The species the feature belongs to.
        against: The comparison set to get; see [1].
        method: The alignment method; see [1].

    Returns:
        A nested dict mapping species to 0-based, stranded locus, represented
        in generic genome browser, to aligned sequence.

    Refs:
        1. https://rest.ensembl.org/documentation/info/compara_species_sets
    """
    alignments = collections.defaultdict(dict)

    locus = feat.locus(base=1, strand=False)
    url = (
        f"https://rest.ensembl.org/"
        f"alignment/region/"
        f"{on}/"
        f"{locus}?"
        f"method={method};"
    )

    if method in ("EPO", "EPO_EXTENDED", "PECAN"):
        url += f"species_set_group={against}"
    else:
        url += f"species_set={on};species_set={against}"

    response = requests.get(url, headers={"Content-Type": "application/json"})
    if response.ok:
        for b, overlap in enumerate(response.json()):
            for alignment in overlap['alignments']:
                species, *pos, seq = _unpack_alignment(alignment)
                chrom, start, end, strand = pos

                locus = f"{chrom}:{start}-{end}:{strand}"
                alignments[species][locus] = seq

    return alignments


def _unpack_alignment(alignment: dict[str: str | int]) -> tuple:
    """Unpacks the alignment from an Ensembl Compara REST API respone.

    Args:
        alignment: A dict representing the alignment.

    Returns:
        A tuple, where tuple[0] is the species name; tuple[1] the chromosome;
        tuple[2] the 0-based, inclusive starting coordinate; tuple[3] the
        0-based, exclusive ending coordinate; tuple[4] is the strand, plus or
        minus for forward or reverse; and tuple[5] the sequence.
    """
    return (
        alignment['species'],
        alignment['seq_region'],
        alignment['start'] - 1,
        alignment['end'],
        "+" if alignment['strand'] == 1 else "-",
        alignment['seq']
    )


if __name__ == '__main__':
    unittest.main()
