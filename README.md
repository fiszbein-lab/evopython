# evopython
`evopython` is an object-oriented Python package designed for genome-scale
feature resolution from whole-genome alignment data.

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
3. these `Feature` instances are resolved from whole-genome alignment MAF
files, represented with the `MAF` class.

In general, we have analyses of the form:
```python
from evopython.gtf import GTF
from evopython.maf import MAF

genes = GTF("path/to/genes")
wga = MAF("path/to/wga", aligned_on="homo_sapiens")

for gene_name in genes:
    feat = genes[gene_name]['feat']
    alignment = wga.get(feat)
    
    if len(alignment) == 1:
        # The alignment is contiguous; do something.
        pass
    else:
        # The alignment is not contiguous; do something else.
        pass
```
We can parse any feature type in the GTF (see the `GTF` class description 
below) and further generate derived, secondary features with the `Feature` 
class's pad method (see the `Feature` class description below), placing minimal
constraints on the features we can query alignments for.

For specific usage examples, see the Jupyter notebooks in the **examples** 
directory.

# Documentation
### class `evopython.gtf.GTF(gtf: str, types: list)`
> A nested `dict` mapping gene name to feature name to a list of `Feature`
> instances; each high-level gene `dict` has two additional keys, "attr" and 
> "feat," with the former mapping to a `dict` with information such as 
> "gene_biotype" indicating whether the gene is protein-coding and the latter 
> mapping to the gene's `Feature` instance.
> 
> *Arguments*:
> - `gtf`: The GTF file path.
> - `types`: The feature types to parse.
----
### class `evopython.bed.BED(bed: str, on_name: bool = False)`
> A `dict` mapping locus `tuple` or name value to `Feature` instance; in the 
> former case, loci have the form `(seqname, start, end, strand)`.
> 
> *Arguments*:
> - `bed`: The BED file path.
> - `on_name`: A `bool` expressing whether name field values should be 
used as keys to the features.
----
### class `evopython.feature.Feature`
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
> of the locus; 1 is used for forward and 0 for reverse.
>
> *Raises*:
> - `ValueError`: An invalid base was given.
> ----
> `pad(self, pad5: int, pad3: int, center: int = 0)`
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
### class `evopython.maf.MAF`
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
> - A list of dictionaries mapping species to `tuple`, where `tuple[0]` is the 
> chromosome name; `tuple[1]` the 0-based, inclusive starting coordinate; 
> `tuple[2]` the 0-based, exclusive ending coordinate; `tuple[3`] the strand, 
> plus or minus for forward or reverse; and `tuple[4]` the alignment.

# Testing

To test feature resolution,
1. clone the repository with
`git clone https://github.com/fiszbein-lab/evopython`,
2. download the needed MAF files into their respective directories using the 
provided FTP links; and
3. run the test from the command line with `python -m unittest tests/test.py`.

The test features have been pre-generated, but new features can be created
with the `features.py` command-line script
```pycon
python features.py --maf path/to/maf --aligned-on species_name
```
where `--aligned-on` is the name of the species that the files are indexed on; 
atext file, **features.txt**, will be written in the same directory as the MAF file.
