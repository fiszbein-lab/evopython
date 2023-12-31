# evopython
[![pypi version](https://img.shields.io/pypi/v/evopython)](https://pypi.org/project/evopython/0.0.1/)

`evopython` is an object-oriented Python package designed for genome-scale
feature resolution from whole-genome alignment data.
- [Installation](#installation)
- [Usage](#usage)
- [Documenation](#documentation)
- [Testing](#testing)
---

## Installation
`evoython` depends on just 
[`biopython`](https://github.com/biopython/biopython) and can be installed with
```commandline
pip install evopython
```

## Usage
`evopython` parses genome annotation data and whole-genome alignment data, 
providing an interface for accessing the former in the context of the latter.
*Ensembl* is a great resource for both; whole-genome alignment 
[MAF](https://genome.ucsc.edu/FAQ/FAQformat.html#format5) files and gene
annotation [GTF](https://genome.ucsc.edu/FAQ/FAQformat.html#format4) files can
be downloaded from their FTP site, indexed 
[here](https://useast.ensembl.org/info/data/ftp/index.html).

`evopython` was designed to support a linear, progressive form of analysis, 
where
1. `GTF` or `BED` instances are initialized with their respective files;
2. these dictionary-like data structures are queried for instances of the
`Feature` class; and
3. these `Feature` instances are resolved from a pairwise or multiple 
whole-genome alignment, represented with the `MAF` class.

In general, we have analyses of the form:
```python
from evopython import GTF, MAF

genes = GTF("path/to/genes")
wga = MAF("path/to/wga", aligned_on="species_name")

for gene_name in genes:
    feat = genes[gene_name]['feat']
    alignments = wga.get(feat)
    
    if len(alignments) == 1:
        # The alignment is contiguous; do something.
        pass
    else:
        # The alignment is discontiguous; do something else.
        pass
```
We can parse any feature type in the GTF (see the `GTF` description 
below) and further generate derived, secondary features with the `Feature` 
class's pad method (see the `Feature` description below).

For specific usage examples, see the Jupyter notebooks in the **examples** 
directory.

## Documentation
### class `evopython.GTF(gtf: str, types: tuple = tuple())`
> A nested `dict` mapping gene name to feature name to a list of `Feature`
> instances; each high-level gene `dict` has two additional keys, `attr` and 
> `feat`, with the former mapping to a `dict` with annotated information such 
> as "gene_biotype" indicating whether the gene is protein-coding and the 
> latter mapping to the gene's `Feature` instance.
> 
> *Arguments*:
> - `gtf`: The GTF file path.
> - `types`: The feature types to parse; gene features are parsed by default.
----
### class `evopython.BED(bed: str, on_name: bool = False)`
> A `dict` mapping locus `tuple` or name value to `Feature` instance; in the 
> former case, loci have the form `(seqname, start, end, strand)`.
> 
> *Arguments*:
> - `bed`: The BED file path.
> - `on_name`: A `bool` expressing whether name field values should be 
used as keys to the features.
----
### class `evopython.Feature`
> A stranded, genomic feature.
>
> ### Attributes:
> - `chrom`: The chromosome name.
> - `start`: The forward-mapped, 0-based, inclusive starting coordinate.
> - `end`: The forward-mapped, 0-based, exclusive ending coordinate.
> - `strand`: The strand, plus or minus for forward or reverse.
> ----
> ### Instance properties:
> - `is_forward`: A `bool` expressing forward strand orientation.
> - `is_reverse`: A `bool` expressing reverse strand orientation.
> ----
> ### Methods:
> 
> `locus(self, base: int = 0, strand: bool = False)`
> 
> Returns the locus in a generic genome browser format.
> 
> *Arguments*:
> - `base`: The coordinate system to use, 0 or 1, where the former is
> half-open on the end and the latter fully closed.
> - `strand`: A `bool` expressing whether to include the strand at the end
> of the locus.
>
> *Raises*:
> - `ValueError`: An invalid base was given.
> ----
> `pad(self, pad5: int, pad3: int, center: int = 0)`
> 
> Pads the feature.
> 
> Positive padding is tanatamount to feature extension and negative 
> padding to feature shrinkage; with centering, both can be used to 
> derive features that do not overlap the source feature.
> 
> *Arguments*:
> - `pad5`: The number of bases to add to the 5'-end of the feature.
> - `pad3`: The number of bases to add to the 3'-end of the feature.
> - `center`: 5, 3, or 0, indicating how to, or to not, center the padding: 
passing 5 prompts 5'-centering, such that padding is applied on the 5' 
coordinate; 3 likewise prompts 3'-centering; and 0, the default, prompts no 
centering, such that the whole feature is padded.
> 
> *Returns:*
> - A new, padded `Feature` instance.
----
### class `evopython.MAF`
> A resolver for multiple alignment formatted whole-genome alignment data.
>
> *Arguments*:
> - `maf_dir`: The path to a directory of MAF files, each following the 
naming scheme *chromosome_name.maf*.
> - `aligned_on`: The species that the chromosome names correspond to.
>
> ### Methods:
> 
> `get(self, feat: Feature, match_strand: bool = True)`
> 
> Finds an alignment for a feature.
> 
> *Arguments:*
> - `feat`: The feature to get an alignment for.
> - `match_strand`: A bool expressing whether the alignment should match the 
> feature's strand; if `False`, the alignment is mapped to the forward strand.
>
> *Returns:*
> - A `list` of `dicts` mapping species to `tuple`, where `tuple[0]` is a 
> `Feature` instance describing the alignment's position and `tuple[1]` the 
> aligned sequence.

## Testing

To test feature resolution,
1. clone the repository with
`git clone https://github.com/fiszbein-lab/evopython`,
2. download the MAF files into their respective directories using the 
provided FTP links;
3. generate random test features using the command-line script, 
`python features.py --maf path/to/maf --aligned-on species_name`, where 
`aligned_on` is the name of the species the file is indexed on 
(see `tests/data/meta_data.csv`); and
4. run the test from the command line with 
`python -m unittest tests/test_resolution.py`.
