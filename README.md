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
Ensembl is a great resource for both; whole-genome alignment 
[MAF](https://genome.ucsc.edu/FAQ/FAQformat.html#format5) files and gene
annotation [GTF](https://genome.ucsc.edu/FAQ/FAQformat.html#format4) files can
be downloaded from their FTP site, indexed 
[here](https://useast.ensembl.org/info/data/ftp/index.html).

`evopython` was designed to support a linear, progressive form of analysis:
1. `GTF` or `BED` classes are initialized with their respective files;
2. instances of the `Feature` class, a data container for stranded, genomic 
interval representation, are gathered from the above parser class instance; 
and
3. these instances are passed to the `get()` method of an instance of the `MAF` 
class, a representation of the whole-genome alignment.

For example, we could resolve the HES3 transcription factor's core promoter 
from a multiple whole-genome alignment comprising 10 primate species:
```python
from evopython.gtf import GTF
from evopython.maf import MAF

genes = GTF("path/to/genes")
wga = MAF("path/to/wga", aligned_on="homo_sapiens")

HES3_gene = genes['HES3']['feature']
HES3_prom = HES3_gene.pad(center=5, pad5=50, pad3=0)

alignments = wga.get(HES3_prom)

if len(alignments) == 1:  # The alignment is contiguous.
    alignment, = alignments
    
    for species in alignment:
        *pos, seq = alignment[species]
        print("{: <20}".format(species), seq)
```
```
homo_sapiens         ACATGTAAACGAGGGT-CCCTATAAAGGCGGCGCTGCCTCGC-ACTGGTGCA
pan_paniscus         ACATGTAAACGAGGGT-CCCTATAAAGGCGGCGCTGCCTCGC-ACTGGTGCA
pan_troglodytes      ACATGTAAACGAGGGT-CCCTATAAAGGCGGCGCTGCCTCGC-ACTGGTGCA
gorilla_gorilla      ACATGTAAACAAGGGT-CCCTATAAAGGCGGCGCTGCCTCGC-ACTGGTGCA
pongo_abelii         ACATGTAAACGAGGGT-CCCTATAAAGGCGGCGCTGCCTCGC-ACTGGTGCA
nomascus_leucogenys  ACATGTAAACGAGGGT-CCCTATAAAGGCGGCGCTGCCTCGC-ACTGGCACA
chlorocebus_sabaeus  ACATGTAAACGAGGGT-CCCTATAAAGGCAGCACTGCCTCAC-ACTGGTGCA
macaca_fascicularis  ACATGTAAACGAGGGT-CCCTATAAAGGCAGCACTGCCTCAC-ACTGGTGCA
macaca_mulatta       ACATGTAAACGAGGGT-CCCTATAAAGGCAGCACTGCCTCAC-ACTGGTGCA
microcebus_murinus   ACATGTAAACGAGGGGCCCCAATAAAGGCGGCACTGACTTGCTGCTTGCGCA
```
For more detailed examples, see the Jupyter notebooks in the *examples*
directory.

# Index
### *class* `evopython.gtf.GTF(gtf: str, types: list)`
> A dictionary-like data structure for feature representation.
> 
> *Arguments*:
> - `gtf`: The GTF file path.
> - `types`: The feature types to parse.
----
### *class* `evopython.bed.BED(bed: str, on_name: bool = False)`
> A dictionary-like data structure for feature representation.
> 
> *Arguments*:
> - `bed`: The BED file path.
> - `on_name`: A bool expressing whether name field values should be 
used as keys to the features.
----
### *class* `evopython.feature.Feature`
> A stranded, genomic feature.
>
> ### Attributes:
> - `chrom`: The chromosome name.
> - `start`: The forward-mapped, 0-based, inclusive starting coordinate.
> - `end`: The forward-mapped, 0-based, exclusive ending coordinate.
> - `strand`: The strand, plus or minus for forward or reverse.
> ----
> ### Instance properties:
> - `is_forward`: A bool expressing forward strand orientation.
> - `is_reverse`: A bool expressing reverse strand orientation.
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
> - `strand`: A bool expressing whether to include the strand at the end
> of the locus; 1 is used for forward and 0 for reverse.
>
> *Raises*:
> - ValueError: An invalid base was given.
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
passing 5 prompts 5'-centering, such that padding is  applied on the 5' 
coordinate; 3 likewise prompts 3'-centering; and 0, the default, prompts no 
centering, such that the whole feature is padded.
> 
> *Returns:*
> - A new, padded `Feature` instance.
----
### *class* `evopython.maf.MAF`
> A resolver for multiple alignment formatted whole-genome alignment data.
>
> *Arguments*:
> - `maf_dir`: The path to a directory of MAF files, each following the 
naming scheme `<chromosome>.maf`.
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
> - A list of dicts mapping species to tuple, where `tuple[0]` is the 
> chromosome name; `tuple[1]` the 0-based, inclusive starting coordinate; 
> `tuple[2]` the 0-based, exclusive ending coordinate; `tuple[3`] the strand, 
> plus or minus for forward or reverse; and `tuple[4]` the alignment.
> 
