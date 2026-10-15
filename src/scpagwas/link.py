"""Link pathways to GWAS SNPs through gene annotation and LD."""

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class PathwayLink:
    """Per-pathway SNP data ready for regression.

    ``y`` is ``beta^2``; ``ldsq`` is the squared LD matrix used as weights.
    """

    pathway: str
    block_info: pd.DataFrame
    snps: pd.DataFrame
    y: np.ndarray
    ld_matrix: np.ndarray
    ldsq: np.ndarray


def make_ld_matrix(snps, ld_rows):
    """Build a symmetric LD matrix over ``snps``; missing pairs default to 0."""
    pos = {s: i for i, s in enumerate(snps)}
    M = np.eye(len(snps))
    for a, b, r in zip(ld_rows["SNP_A"], ld_rows["SNP_B"], ld_rows["R"]):
        if a in pos and b in pos:
            i, j = pos[a], pos[b]
            M[i, j] = M[j, i] = float(r)
    return M


def link_pathway(pathway, blocks_by_chrom, gwas_by_chrom, snp_gene, ld_index):
    """Collect a pathway's SNPs across chromosomes and build its LD weights."""
    frames, ys = [], []
    for chrom, block in blocks_by_chrom.items():
        gw = gwas_by_chrom.get(int(chrom))
        if gw is None or gw.empty:
            continue
        sg = snp_gene[snp_gene["label"].isin(block["label"])]
        merged = sg[["rsid", "label"]].merge(
            gw[["rsid", "chrom", "pos", "beta", "se"]], on="rsid", how="inner"
        )
        if merged.empty:
            continue
        frames.append(merged)
        ys.append(merged["beta"].to_numpy() ** 2)
    if not frames:
        return None
    snps = pd.concat(frames, ignore_index=True)
    y = np.concatenate(ys)

    ld_frames = [ld_index.get(int(c)) for c in blocks_by_chrom]
    ld_rows = (pd.concat(ld_frames, ignore_index=True)
               if ld_frames else pd.DataFrame(columns=["SNP_A", "SNP_B", "R"]))
    snp_ids = set(snps["rsid"])
    if not ld_rows.empty:
        ld_rows = ld_rows[ld_rows["SNP_A"].isin(snp_ids)
                          & ld_rows["SNP_B"].isin(snp_ids)]
    M = make_ld_matrix(snps["rsid"].tolist(), ld_rows)

    block_info = pd.concat(blocks_by_chrom.values(), ignore_index=True)
    block_info["pathway"] = pathway
    return PathwayLink(pathway, block_info, snps, y, M, M * M)
