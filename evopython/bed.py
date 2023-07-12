import collections
import csv

from evopython.feature import Feature


class BED(collections.UserDict):
    def __init__(self, bed: str, on_name: bool = False):
        """Inits. BED.

        Args:
            bed: The BED file path.
            on_name: A bool expressing whether name field values should be
                used to hash the features.
        """
        featurome = _read_bed(bed, on_name=on_name)
        super().__init__(featurome)


def _read_bed(bed: str, on_name: bool = False) -> dict[tuple | str: Feature]:
    """Reads BED file.

    Args:
        bed: The BED path.
        on_name: A bool expressing whether name field values should be used as
            keys to the features.

    Returns:
        A dict mapping locus tuple or name value to `Feature` instance; in the
        former case, tuples have the form `(seqname, start, end, strand)`.
    """
    featurome = dict()

    with open(bed, 'r') as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            chrom, start, end, name, _, strand = row

            if on_name:
                on = name
            else:
                on = chrom, start, end, strand

            featurome[on] = Feature(chrom, start, end, strand)

    return featurome
