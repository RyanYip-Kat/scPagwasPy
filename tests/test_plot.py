import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")
os.environ.setdefault("MPLBACKEND", "Agg")

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd

from scpagwas.plot.bar import plot_positive_cells_bar
from scpagwas.plot.bar import plot_bootstrap_bar
from scpagwas.plot.forest import plot_bootstrap_forest
from scpagwas.plot.gene import plot_gene_pcc
from scpagwas.plot.umap import plot_significant_cells, plot_umap_scores
from scpagwas.plot.violin import plot_top_gene_violin


def _adata(n=20):
    rng = np.random.default_rng(0)
    obs = pd.DataFrame(
        {"cell_type": ["A"] * (n // 2) + ["B"] * (n - n // 2),
         "scPagwas.gPAS.score": rng.random(n),
         "Random_Correct_BG_adjp": rng.random(n)},
        index=[f"c{i}" for i in range(n)],
    )
    a = ad.AnnData(X=rng.random((n, 5)), obs=obs,
                   var=pd.DataFrame(index=[f"g{i}" for i in range(5)]))
    a.obsm["X_umap"] = rng.random((n, 2))
    return a


def test_plot_umap(tmp_path):
    plot_umap_scores(_adata(), save_path=tmp_path / "u.png")
    assert (tmp_path / "u.png").exists()


def test_plot_significant_cells(tmp_path):
    plot_significant_cells(_adata(), save_path=tmp_path / "s.png")
    assert (tmp_path / "s.png").exists()


def test_plot_forest(tmp_path):
    df = pd.DataFrame({"bias_corrected_estimate": [0.1, -0.2],
                       "CI_lo": [0.0, -0.3], "CI_hi": [0.2, -0.1],
                       "bp_value": [0.01, 0.2]}, index=["A", "B"])
    plot_bootstrap_forest(df, save_path=tmp_path / "f.png")
    assert (tmp_path / "f.png").exists()


def test_plot_bootstrap_bar(tmp_path):
    df = pd.DataFrame({"bp_value": [0.001, 0.02, 0.5]}, index=["A", "B", "C"])
    plot_bootstrap_bar(df, save_path=tmp_path / "bp.png")
    assert (tmp_path / "bp.png").exists()


def test_plot_positive_cells_bar(tmp_path):
    plot_positive_cells_bar(_adata(), save_path=tmp_path / "b.png")
    assert (tmp_path / "b.png").exists()


def test_plot_gene_pcc(tmp_path):
    pcc = pd.Series(np.linspace(-0.5, 0.9, 30), index=[f"g{i}" for i in range(30)])
    plot_gene_pcc(pcc, save_path=tmp_path / "p.png")
    assert (tmp_path / "p.png").exists()


def test_plot_violin(tmp_path):
    plot_top_gene_violin(_adata(), ["g0", "g1"], save_path=tmp_path / "v.png")
    assert (tmp_path / "v.png").exists()
