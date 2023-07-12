from __future__ import annotations
from dataclasses import dataclass


@dataclass(eq=True, frozen=True)
class Feature:
    """A genomic feature.

    Attr:
        chrom: The chromosome name.
        start: The forward-mapped, 0-based, inclusive starting coordinate.
        end: The forward-mapped, 0-based, exclusive ending coordinate.
        strand: The strand, plus or minus for forward or reverse.
    """
    chrom: str
    start: int
    end: int
    strand: str

    @property
    def is_forward(self) -> bool:
        """A bool expressing forward-strand orientation."""
        return self.strand == "+"

    @property
    def is_reverse(self) -> bool:
        """A bool expressing reverse-strand orientation."""
        return self.strand == "-"

    def locus(self, base: int = 0, strand: bool = False) -> str:
        """Returns the locus in a generic genome browser format.

        Args:
            base: The coordinate system to use, 0 or 1, where the former is
                half-open on the end and the latter fully closed.
            strand: A bool expressing whether to include the strand at the end
                of the locus; 1 is used for forward and -1 for reverse.

        Raises:
            ValueError: An invalid base was given.
        """
        match base:
            case 0:
                start = self.start
            case 1:
                start = self.start + 1
            case _:
                raise ValueError(f"`base` must be in 0 or 1 — got {base}.")

        locus = f"{self.chrom}:{start}-{self.end}"

        if strand:
            locus += ":1" if self.is_forward else ":0"

        return locus

    def pad(self, pad5: int, pad3: int, center: int = 0) -> Feature:
        """Pads the feature.

        Positive padding is tanatamount to feature extension and negative
        padding to feature shrinkage; with centering, both can be used to
        derive features that do not overlap the original feature.

        Args:
            pad5: The number of bases to add to the 5'-end of the feature.
            pad3: The number of bases to add to the 3'-end of the feature.
            center: 5, 3, or 0, indicating how to, or to not, center the
                padding: passing 5 prompts 5'-centering, such that padding is
                applied on the 5' coordinate; 3 likewise prompts 3'-centering;
                and 0, the default, prompts no centering, such that the whole
                feature is padded.

        Returns:
            A new, padded `Feature` instance.

        Raises:
            ValueError: The `center` argument was not 5, 3, or 0.
        """
        if center in (5, 3):
            if self.is_forward:
                if center == 5:
                    start = self.start - pad5
                    end = self.start + pad3
                else:
                    start = self.end - pad5
                    end = self.end + pad3
            else:
                if center == 5:
                    start = self.end - pad3
                    end = self.end + pad5
                else:
                    start = self.start - pad3
                    end = self.start + pad5
        elif center == 0:
            if self.is_forward:
                start = self.start - pad5
                end = self.end + pad3
            else:
                start = self.end - pad3
                end = self.end + pad5
        else:
            raise ValueError(f"`center` must be 5, 3, or 0 — got {center}.")

        return Feature(self.chrom, start, end, self.strand)

    def __len__(self) -> int:
        """Returns the feature's length."""
        return self.end - self.start
