"""LD (linkage disequilibrium) input: per-chromosome R / R^2 lookup."""

import pandas as pd

_EMPTY = ["SNP_A", "SNP_B", "R"]


class LDIndex:
    """Per-chromosome LD long tables keyed by chromosome number."""

    def __init__(self, frames):
        self.frames = frames

    @classmethod
    def from_frames(cls, frames):
        return cls({int(k): v for k, v in frames.items()})

    @property
    def chromosomes(self):
        return sorted(self.frames)

    def get(self, chrom):
        return self.frames.get(int(chrom), pd.DataFrame(columns=_EMPTY))

    def restrict(self, rsids):
        """Return a new index keeping only pairs whose both SNPs are in ``rsids``."""
        keep = set(rsids)
        frames = {}
        for chrom, df in self.frames.items():
            mask = df["SNP_A"].isin(keep) & df["SNP_B"].isin(keep)
            frames[chrom] = df.loc[mask]
        return LDIndex(frames)
