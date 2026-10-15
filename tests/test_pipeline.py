import anndata as ad
import numpy as np
import pandas as pd

from scpagwas.config import PagwasConfig
from scpagwas.io.ld import LDIndex
from scpagwas.pipeline import run


def _adata():
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(120)]
    cells = [f"c{i}" for i in range(120)]
    X = rng.random((120, 120))
    obs = pd.DataFrame({"cell_type": ["A"] * 40 + ["B"] * 40 + ["C"] * 40},
                       index=cells)
    return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes)), genes


def _gwas(genes):
    rows = [(1, 1000 + i * 100, f"rs{i}", 0.01 * (i % 5), 0.1, 0.2)
            for i in range(len(genes))]
    return pd.DataFrame(rows, columns=["chrom", "pos", "rsid", "beta", "se", "maf"])


def _annotation(genes):
    return pd.DataFrame({
        "chrom": [1] * len(genes),
        "start": [1000 + i * 100 for i in range(len(genes))],
        "end": [1000 + i * 100 + 50 for i in range(len(genes))],
        "label": genes,
    })


def test_pipeline_end_to_end(tmp_path):
    a, genes = _adata()
    pathways = {"P1": genes[:40], "P2": genes[40:80], "P3": genes[80:]}
    ld = LDIndex.from_frames({1: pd.DataFrame(columns=["SNP_A", "SNP_B", "R"])})
    cfg = PagwasConfig(min_clustercells=5, min_pathway_size=5, max_pathway_size=100,
                       iters_celltype=5, iters_singlecell=3, n_topgenes=10)
    res = run(a, _gwas(genes), pathways, _annotation(genes), ld, config=cfg)

    assert res.cell_scores.shape[0] == a.n_obs
    assert {"scPagwas.TRS.Score", "scPagwas.gPAS.score",
            "Random_Correct_BG_p"} <= set(res.cell_scores.columns)
    assert "scPagwas.gPAS.score" in res.adata.obs
    assert res.adata.obsm["scPagwas_pa_pca"].shape[0] == a.n_obs
    assert not res.gene_pcc.empty
    assert len(res.to_files("t", tmp_path)) == 5


def test_pipeline_records_celltype_bootstrap():
    a, genes = _adata()
    pathways = {"P1": genes[:40], "P2": genes[40:80], "P3": genes[80:]}
    ld = LDIndex.from_frames({1: pd.DataFrame(columns=["SNP_A", "SNP_B", "R"])})
    cfg = PagwasConfig(min_clustercells=5, min_pathway_size=5, max_pathway_size=100,
                       iters_celltype=5, iters_singlecell=3, n_topgenes=10)
    res = run(a, _gwas(genes), pathways, _annotation(genes), ld, config=cfg)
    assert list(res.celltype_bootstrap.index) == ["Intercept", "A", "B", "C"]
    assert set(res.merged_celltype_p["celltype"]) == {"A", "B", "C"}
