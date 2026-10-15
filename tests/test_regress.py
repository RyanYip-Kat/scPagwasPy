import numpy as np
import pandas as pd

from scpagwas.config import PagwasConfig
from scpagwas.link import PathwayLink
from scpagwas.regress import (bootstrap_celltype, build_celltype_design,
                              celltype_regression, singlecell_regression)


def _link():
    snps = pd.DataFrame({
        "chrom": [1, 1, 1],
        "pos": [1, 2, 3],
        "rsid": ["a", "b", "c"],
        "beta": [0.2, 0.3, 0.1],
        "se": [0.1, 0.1, 0.1],
        "maf": [0.2, 0.2, 0.2],
        "label": ["G1", "G1", "G2"],
    })
    M = np.eye(3)
    return PathwayLink("P", snps.assign(pathway="P"), snps,
                       snps["beta"].to_numpy() ** 2, M, M * M)


def test_singlecell_regression_shape():
    link = _link()
    cells = pd.Index(["c1", "c2", "c3", "c4"])
    pca_mat = pd.DataFrame(np.ones((1, 4)), index=["P"], columns=cells)
    coefs = singlecell_regression(link, pca_mat, np.ones((2, 4)),
                                  ["G1", "G2"], ["G1", "G2"], cells)
    assert coefs.shape == (4,)


def test_celltype_regression_index():
    link = _link()
    pca_df = pd.DataFrame({"A": [1.0], "B": [0.5]}, index=["P"])
    merge = pd.DataFrame({"A": [1.0, 2.0], "B": [2.0, 1.0]}, index=["G1", "G2"])
    x, y, noise = build_celltype_design(link, pca_df, merge, ["G1", "G2"])
    params = celltype_regression([(x, y, noise)], columns=["A", "B"])
    assert list(params.index) == ["Intercept", "A", "B"]


def test_build_celltype_design_handles_duplicate_genes():
    link = _link()
    pca_df = pd.DataFrame({"A": [1.0], "B": [0.5]}, index=["P"])
    merge = pd.DataFrame({"A": [1.0, 2.0], "B": [2.0, 1.0]}, index=["G1", "G2"])
    x, y, noise = build_celltype_design(link, pca_df, merge, ["G1", "G1", "G2"])
    assert x.shape == (3, 2)


def test_bootstrap_celltype_columns():
    link = _link()
    pca_df = pd.DataFrame({"A": [1.0], "B": [0.5]}, index=["P"])
    merge = pd.DataFrame({"A": [1.0, 2.0], "B": [2.0, 1.0]}, index=["G1", "G2"])
    x, y, noise = build_celltype_design(link, pca_df, merge, ["G1", "G2"])
    base = celltype_regression([(x, y, noise)], columns=["A", "B"])
    out = bootstrap_celltype([(x, y, noise)], ["A", "B"], base.to_numpy(),
                             PagwasConfig(iters_celltype=5, seed=1))
    expected = {"bootstrap_estimate", "bootstrap_error", "bt_value", "bp_value",
                "bias_corrected_estimate", "CI_lo", "CI_hi"}
    assert expected <= set(out.columns)
    assert list(out.index) == ["Intercept", "A", "B"]
