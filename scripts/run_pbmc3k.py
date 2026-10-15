#!/usr/bin/env python
"""Run the scpagwas pipeline on testdata/pbmc3k_data.h5ad with ported assets.

Result files are written by ``PagwasResult.to_files`` (unchanged). On top of
that this script saves the figure set mirroring the R scPagwas plotting
functions:

  scPagwas_Visualization      -> scPagwas.gPAS.score / scPagwas.TRS.Score UMAPs
                                 and the significant-cell UMAP
  Bootstrap_estimate_Plot     -> bootstrap forest plot
  Bootstrap_P_Barplot         -> bootstrap -log2(p) bar plot
  heritability_cor_scatterplot-> gene PCC scatter
  plot_bar_positie_nagtive    -> positive-cell fraction bar chart
  plot_vln_Corgenes           -> top-gene violin grid
"""

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")

import matplotlib  # noqa: E402

matplotlib.use("Agg")

import time  # noqa: E402

import anndata as ad  # noqa: E402

from scpagwas.config import PagwasConfig  # noqa: E402
from scpagwas.io import assets  # noqa: E402
from scpagwas.io.gwas import read_gwas  # noqa: E402
from scpagwas.pipeline import run  # noqa: E402


def make_figures(result, outdir, p_thre=0.05):
    """Save the R-equivalent scPagwas figure set; return the written paths."""
    from scpagwas.plot.bar import plot_bootstrap_bar, plot_positive_cells_bar
    from scpagwas.plot.forest import plot_bootstrap_forest
    from scpagwas.plot.gene import plot_gene_pcc
    from scpagwas.plot.umap import plot_significant_cells, plot_umap_scores
    from scpagwas.plot.violin import plot_top_gene_violin

    figdir = Path(outdir) / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    adata = result.adata
    pcc = result.gene_pcc.sort_values(ascending=False)
    top_genes = list(pcc.index[:6])

    figures = {
        "scPagwas.gPAS.score_umap.png": lambda p: plot_umap_scores(
            adata, key="scPagwas.gPAS.score", save_path=p
        ),
        "scPagwas.TRS.Score_umap.png": lambda p: plot_umap_scores(
            adata, key="scPagwas.TRS.Score", save_path=p
        ),
        "scPagwas_Correct_BG_adjp0.05_umap.png": lambda p: plot_significant_cells(
            adata, p_thre=p_thre, save_path=p
        ),
        "Bootstrap_estimate_forest.png": lambda p: plot_bootstrap_forest(
            result.celltype_bootstrap, save_path=p
        ),
        "Bootstrap_P_barplot.png": lambda p: plot_bootstrap_bar(
            result.celltype_bootstrap, save_path=p
        ),
        "gene_PCC.png": lambda p: plot_gene_pcc(pcc, topn=15, save_path=p),
        "positive_cells_bar.png": lambda p: plot_positive_cells_bar(
            adata, cell_type_key="cell_type", p_thre=p_thre, save_path=p
        ),
        "top_gene_violin.png": lambda p: plot_top_gene_violin(
            adata, top_genes, cell_type_key="cell_type", save_path=p
        ),
    }

    written = []
    for name, make in figures.items():
        path = figdir / name
        make(path)
        written.append(path)
        print(f"[scpagwas] figure: {path}", flush=True)
    return written


def main():
    t0 = time.time()
    adata = ad.read_h5ad("testdata/pbmc3k_data.h5ad")
    adata.obs["cell_type"] = adata.obs["leiden"].astype(str)
    print(f"loaded adata {adata.shape} in {time.time() - t0:.1f}s", flush=True)

    gwas = read_gwas(assets.example_gwas_path())
    pathways = assets.load_builtin_pathways("kegg")
    anno = assets.load_builtin_block_annotation("hg38")
    print(f"gwas={len(gwas)} pathways={len(pathways)} blocks={len(anno)}", flush=True)

    t1 = time.time()
    ld = assets.load_builtin_ld()
    print(f"loaded LD {sum(len(ld.get(c)) for c in ld.chromosomes)} rows "
          f"in {time.time() - t1:.1f}s", flush=True)

    cfg = PagwasConfig(iters_celltype=100, iters_singlecell=50, n_topgenes=500)
    t2 = time.time()
    result = run(adata, gwas, pathways, anno, ld, config=cfg)
    print(f"pipeline finished in {time.time() - t2:.1f}s", flush=True)

    files = result.to_files("pbmc3k", "outputs/pbmc3k")
    print("wrote:", [str(f) for f in files], flush=True)

    make_figures(result, "outputs/pbmc3k")

    print(result.celltype_bootstrap.round(4).to_string(), flush=True)
    print(result.merged_celltype_p.round(4).to_string(), flush=True)


if __name__ == "__main__":
    main()
