"""Access the built-in assets ported from the scPagwas R package.

Assets live in ``scpagwas/data/`` and are produced by
``scripts/port_r_assets.R`` (pathways, gene blocks, LD, example GWAS).
"""

from pathlib import Path

import pandas as pd

from .annotation import load_block_annotation
from .ld import LDIndex
from .pathways import read_gmt


def asset_dir():
    """Directory holding the bundled data assets."""
    return Path(__file__).resolve().parent.parent / "data"


def load_builtin_pathways(name="kegg"):
    """Load a bundled pathway gene-set GMT (``kegg``, ``gobp``, ``reactome``)."""
    path = asset_dir() / f"pathways_{name}.gmt"
    if not path.exists():
        available = sorted(p.stem.replace("pathways_", "")
                           for p in asset_dir().glob("pathways_*.gmt"))
        raise FileNotFoundError(f"no builtin pathways '{name}'; have {available}")
    return read_gmt(path)


def load_builtin_block_annotation(genome="hg38"):
    """Load bundled gene coordinate blocks (``hg38`` or ``hg37``)."""
    suffix = "" if genome in ("hg38", "GRCh38") else f"_{genome}"
    return load_block_annotation(asset_dir() / f"block_annotation{suffix}.tsv.gz")


def load_builtin_ld(chromosomes=None):
    """Load bundled LD tables as an :class:`LDIndex`.

    ``chromosomes`` optionally restricts which chromosomes to read.
    """
    frames = {}
    for path in sorted((asset_dir() / "chrom_ld").glob("chr*.tsv.gz")):
        chrom = int(path.name[3:].split(".")[0])
        if chromosomes is not None and chrom not in set(chromosomes):
            continue
        frames[chrom] = pd.read_csv(path, sep="\t")
    return LDIndex.from_frames(frames)


def example_gwas_path():
    """Path to the bundled example GWAS summary statistics file."""
    return asset_dir() / "gwas_example.txt"
