"""Gene coordinate blocks (``chrom, start, end, label``) input."""

from pathlib import Path

import numpy as np
import pandas as pd

_COLUMNS = ["chrom", "start", "end", "label"]
_CHROM_MAP = {"x": 23, "y": 24, "m": 26, "mt": 26}


def _norm_chrom(c):
    if isinstance(c, str):
        c = c.lower().replace("chr", "")
        if c in _CHROM_MAP:
            return _CHROM_MAP[c]
    try:
        return int(c)
    except (TypeError, ValueError):
        return np.nan


def load_block_annotation(path):
    """Load gene coordinate blocks from parquet/tsv/csv.

    Returns a frame with ``chrom`` (int), ``start``, ``end`` and ``label``.
    Export the bundled R data with ``scripts/export_r_assets.R``.
    """
    path = Path(path)
    if path.suffix == ".parquet":
        df = pd.read_parquet(path)
    elif path.suffix == ".gz":
        df = pd.read_csv(path, sep="\t")
    elif path.suffix == ".csv":
        df = pd.read_csv(path)
    else:
        df = pd.read_csv(path, sep=None, engine="python")
    missing = [c for c in _COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"block annotation missing columns: {missing}")
    out = df.loc[:, _COLUMNS].copy()
    out["chrom"] = out["chrom"].map(_norm_chrom)
    out = out.dropna(subset=["chrom"])
    out["chrom"] = out["chrom"].astype(int)
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)
    return out
