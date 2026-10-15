"""GWAS summary-statistics input and cleaning."""

import numpy as np
import pandas as pd

REQUIRED = ["chrom", "pos", "rsid", "beta", "se", "maf"]

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


def read_gwas(path):
    """Read a whitespace/tab/comma-delimited GWAS file and clean it."""
    df = pd.read_csv(path, sep=None, engine="python")
    return clean_gwas(df, None)


def clean_gwas(gwas, cfg):
    """Filter and normalize GWAS summary statistics.

    Normalizes chromosome to int 1..22, flips MAF > 0.5, drops duplicates,
    filters sex chromosomes and the MHC region, and requires ``abs(beta) < 1``.
    """
    missing = [c for c in REQUIRED if c not in gwas.columns]
    if missing:
        raise ValueError(f"GWAS missing columns: {missing}")
    df = gwas.loc[:, REQUIRED].copy()
    df["chrom"] = df["chrom"].map(_norm_chrom)
    df = df.dropna(subset=["chrom"])
    df["chrom"] = df["chrom"].astype(int)
    df["pos"] = df["pos"].astype(int)
    for c in ("beta", "se", "maf"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["beta", "se", "maf"])
    df["maf"] = np.where(df["maf"] > 0.5, 1 - df["maf"], df["maf"])
    df = df[df["maf"] > (cfg.maf_filter if cfg else 0.01)]
    df = df[~df["rsid"].duplicated()]
    if cfg and cfg.sex_filter:
        df = df[df["chrom"].between(1, 22)]
    if cfg and cfg.mhc_filter:
        df = df[~((df["chrom"] == 6) & df["pos"].between(25_000_000, 34_000_000))]
    if cfg and cfg.gwas_z_filter is not None:
        df = df[(df["beta"] / df["se"]).abs() < cfg.gwas_z_filter]
    df = df[df["beta"].abs() < 1]
    return df.reset_index(drop=True)
