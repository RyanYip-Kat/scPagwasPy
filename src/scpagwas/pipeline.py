"""End-to-end scpagwas pipeline orchestration."""

import numpy as np
import pandas as pd

from .annotate import build_pathway_blocks, map_snps_to_genes
from .config import PagwasConfig
from .genes import gene_pcc
from .io.gwas import clean_gwas
from .link import link_pathway
from .pathway_score import score_pathways
from .preprocess import prepare_singlecell
from .regress import (bootstrap_celltype, build_celltype_design,
                      celltype_regression, singlecell_regression)
from .result import PagwasResult
from .score import compute_gpas
from .significance import background_correction, merge_celltype_p


def run(adata, gwas, pathways, block_annotation, ld, config=None, layer=None):
    """Run the full scpagwas pipeline.

    Parameters
    ----------
    adata : AnnData
        Single-cell data; ``obs["cell_type"]`` holds the cell-type labels and
        ``X`` (or ``layers[layer]``) the log-normalized expression.
    gwas : DataFrame
        Summary statistics with ``chrom, pos, rsid, beta, se, maf``.
    pathways : dict
        ``{pathway: [gene, ...]}``.
    block_annotation : DataFrame
        Gene coordinates with ``chrom, start, end, label``.
    ld : LDIndex
        Per-chromosome LD long tables.
    config : PagwasConfig, optional
    layer : str, optional
        ``adata.layers`` key to use instead of ``X``.
    """
    cfg = config or PagwasConfig()
    sc = prepare_singlecell(adata, pathways, cfg, layer=layer)
    print(f"[scpagwas] prepared {sc.data_mat.shape[0]} genes x "
          f"{sc.data_mat.shape[1]} cells, {len(sc.celltypes)} celltypes",
          flush=True)
    ps = score_pathways(sc, pathways, cfg)
    print(f"[scpagwas] scored {len(ps.pathways)} pathways", flush=True)

    gwas_clean = clean_gwas(gwas, cfg)
    snp_gene = map_snps_to_genes(gwas_clean, block_annotation, cfg.marg)
    blocks = build_pathway_blocks(ps.pathways, block_annotation, snp_gene)
    gwas_by_chrom = {int(c): g for c, g in gwas_clean.groupby("chrom")}
    ld = ld.restrict(set(gwas_clean["rsid"]))
    print(f"[scpagwas] mapped {len(snp_gene)} snp-gene pairs, "
          f"{len(blocks)} pathways with blocks", flush=True)

    links = {}
    for name, bc in blocks.items():
        pl = link_pathway(name, bc, gwas_by_chrom, snp_gene, ld)
        if pl is not None and pl.y.size > 0:
            links[name] = pl
    if not links:
        raise ValueError("No pathway could be linked to GWAS SNPs")
    print(f"[scpagwas] linked {len(links)} pathways to SNPs", flush=True)

    celltypes = list(ps.pca_cell_df.columns)
    designs = {n: build_celltype_design(pl, ps.pca_cell_df, sc.merge_expr, ps.pathways[n])
               for n, pl in links.items()}
    base = celltype_regression(list(designs.values()), columns=celltypes)
    boot = bootstrap_celltype(list(designs.values()), celltypes, base.to_numpy(), cfg)
    print(f"[scpagwas] celltype regression + {cfg.iters_celltype} bootstraps done",
          flush=True)

    dense = np.asarray(sc.data_mat.todense())
    coef = {n: singlecell_regression(pl, ps.pca_cell_mat, dense, list(sc.genes),
                                     ps.pathways[n], list(sc.cells))
            for n, pl in links.items()}
    coef = pd.DataFrame(coef, index=sc.cells).T
    print(f"[scpagwas] single-cell regression done ({coef.shape})", flush=True)

    cs = compute_gpas(coef, ps.pca_cell_mat, dense, list(sc.genes), ps.pathways, cfg)
    pcc = gene_pcc(cs.gpas, dense, sc.genes, sc.cells)
    top = list(pcc.sort_values(ascending=False).index[:cfg.n_topgenes])
    bg = background_correction(sc, top, cfg)
    print(f"[scpagwas] gPAS + gene PCC + background correction done", flush=True)

    cell_scores = pd.DataFrame({
        "scPagwas.TRS.Score": bg["raw_score"].to_numpy(),
        "scPagwas.gPAS.score": cs.gpas.reindex(sc.cells).to_numpy(),
        "Random_Correct_BG_p": bg["pooled_p"].to_numpy(),
        "Random_Correct_BG_adjp": bg["adj_p"].to_numpy(),
        "Random_Correct_BG_z": bg["pooled_z"].to_numpy(),
    }, index=sc.cells)

    adata = adata.copy()
    for col in cell_scores.columns:
        adata.obs[col] = cell_scores[col].reindex(adata.obs_names).to_numpy()
    adata.obsm["scPagwas_pa_pca"] = (
        ps.pca_cell_mat.T.reindex(adata.obs_names).fillna(0.0).to_numpy()
    )

    merged = merge_celltype_p(bg["pooled_p"].to_numpy(), sc.cell_type.to_numpy())
    return PagwasResult(
        adata=adata,
        config=cfg,
        celltype_bootstrap=boot,
        gene_pcc=pcc,
        pathway_cell_lm=coef,
        cell_scores=cell_scores,
        merged_celltype_p=merged,
        pathway_list=ps.pathways,
        pca_cell_df=ps.pca_cell_df,
        lm_results={"parameters": base},
    )
