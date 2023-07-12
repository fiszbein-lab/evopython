import argparse
import os
import random


COUNT = 100


def generate_features(maf: str, aligned_on: str) -> None:
    """Makes features within alignment bounds.

    Notes:
        I assume the 1st species is the "aligned on" species that the feature
        coordinates should map to.

    Args:
        maf: The MAF file path.
        aligned_on: The species that the MAF file's chromosome name corresponds
            to.

    Writes:
        A text file with a feature on each line in the same directory as the
        MAF file, where each feature has the form <chrom:start-end:strand> with
        1 indicating the forward strand and 0 indicating the reverse strand.

    Raises:
        ValueError: The chromosome could not be resolved.
    """
    chrom_name = ""

    chrom_min = -1
    chrom_max = -1

    with open(maf, 'r') as f:
        initial_record = True
        for line in f:
            if not line.startswith("s"):
                continue

            _, source, start, seq_size, strand, chrom_size, _ = line.split()
            species, chrom = source.split(".", maxsplit=1)

            if species != aligned_on:
                continue

            chrom_name = chrom

            seq_size = int(seq_size)
            start = int(start)
            end = start + seq_size

            if strand == "-":
                # Map coordinates to the forward strand.
                start = int(chrom_size) - end
                end = start + seq_size

            if initial_record:
                chrom_min = start
                chrom_max = end
                initial_record = False
            else:
                if start < chrom_min:
                    chrom_min = start
                if end > chrom_max:
                    chrom_max = end

    if not chrom_name or chrom_min < 0 or chrom_max < 0:
        raise ValueError("Unable to resolve chromosome.")

    diranme = os.path.dirname(maf)
    with open(f"{diranme}/features.txt", 'w') as f:
        for _ in range(COUNT):
            start = random.randint(chrom_min, chrom_max)
            strand = random.randint(0, 1)
            size = random.randint(1, 1000)

            if start + size > chrom_max:
                end = chrom_size
            else:
                end = start + size

            f.write(f"{chrom_name}:{start}-{end}:{strand}\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--maf', type=str)
    parser.add_argument('--aligned-on', type=str)

    args = parser.parse_args()
    generate_features(args.maf, args.aligned_on)
