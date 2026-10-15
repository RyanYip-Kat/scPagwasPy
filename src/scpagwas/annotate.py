"""Map GWAS SNPs to genes and build per-pathway annotation blocks.

``map_snps_to_genes`` reproduces ``SnpToGene``: a SNP is assigned to every
gene whose TSS window ``[start - marg, end + marg]`` contains it.
"""

import numpy as np
import pandas as pd


def map_snps_to_genes(gwas, block_annotation, marg):
    """Return a tidy frame of ``rsid, label, pos`` SNP-to-gene assignments."""
    rows = []
    for chrom, gw in gwas.groupby("chrom"):
        anno = block_annotation[block_annotation["chrom"] == chrom]
        if anno.empty:
            continue
        pos = gw["pos"].to_numpy()
        rsids = gw["rsid"].to_numpy()
        order = np.argsort(pos)
        pos_sorted = pos[order]
        rsid_sorted = rsids[order]
        for start, end, label in zip(
            anno["start"].to_numpy(), anno["end"].to_numpy(), anno["label"].to_numpy()
        ):
            lo = np.searchsorted(pos_sorted, start - marg, side="left")
            hi = np.searchsorted(pos_sorted, end + marg, side="right")
            for k in range(lo, hi):
                rows.append((rsid_sorted[k], label, pos_sorted[k]))
    return pd.DataFrame(rows, columns=["rsid", "label", "pos"])


def build_pathway_blocks(pathways, block_annotation, snp_gene, min_blocks=10):
    """Map each pathway to gene coordinate blocks, grouped by chromosome."""
    anno = block_annotation.drop_duplicates(subset=["label"])
    blocks = {}
    for name, genes in pathways.items():
        sub = anno[anno["label"].isin(genes)]
        if len(sub) < min_blocks:
            continue
        sub = sub.sort_values(["chrom", "start"])
        blocks[name] = {int(c): g for c, g in sub.groupby("chrom")}
    return blocks
