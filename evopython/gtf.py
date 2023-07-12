import re
import collections

from evopython.feature import Feature


class GTF(collections.UserDict):
    def __init__(self, gtf: str, types: list):
        """Inits. GTF.

        Args:
            gtf: The GTF file path.
            types: The tuple of feature types to parse.
        """
        featurome = _read_gtf(gtf, types=types)
        super().__init__(featurome)


def _read_gtf(gtf: str, types: list) -> dict[str: dict]:
    """Reads GTF.

    Args:
        gtf: The GTF file path.
        types: The feature type/s to parse.

    Returns:
        A nested dict mapping gene name to feature name to a list of `Feature`
        instances; each high-level gene dict two additional keys, "attr" and
        "feat," with the former mapping to a dict with info. such as
        "gene_biotype" indicating whether the gene is protein-coding and the
        latter mapping to the gene's `Feature` instance.
    """
    featurome = dict()

    if "gene" not in types:
        types.append("gene")

    with open(gtf, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue

            row = line.split("\t")
            chrom, _, type, *pos, _, strand, _, attr = row

            if type in types:
                start, end = map(int, pos)
                feat = Feature(chrom, start, end, strand)

                gene_attr = _unpack_attr(attr)
                if 'gene_name' in gene_attr:
                    gene = gene_attr['gene_name']
                else:
                    gene = gene_attr['gene_id']

                if gene not in featurome:
                    featurome[gene] = {
                        type: list() for type in types if type != "gene"}

                    featurome[gene]['attr'] = gene_attr
                    featurome[gene]['feat'] = feat

                if type == "gene":
                    featurome[gene][feat] = feat
                else:
                    featurome[gene][feat].append(feat)

    return featurome


def _unpack_attr(attr_str: str) -> dict[str: str | int]:
    """Unpacks GTF attributes.

    Args:
        attr_str: The attribute field of a GTF file line.

    Returns:
        A dict mapping attribute tag to value.
    """
    attr_dict = dict()

    for pair in attr_str.split(";"):
        pair = pair.strip()

        if not pair or "(" in pair or re.match('^tag', pair):
            continue

        tag, attr = pair.split()
        attr = attr.strip('"')

        attr_dict[tag] = int(attr) if attr.isdigit() else attr

    return attr_dict
