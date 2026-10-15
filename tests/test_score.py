import numpy as np
import pandas as pd

from scpagwas.config import PagwasConfig
from scpagwas.score import compute_gpas


def test_compute_gpas_sum_over_pathways():
    cells = pd.Index(["c1", "c2", "c3"])
    coef = pd.DataFrame({"c1": [1.0], "c2": [2.0], "c3": [3.0]}, index=["P"])
    pca = pd.DataFrame({"c1": [1.0], "c2": [1.0], "c3": [1.0]}, index=["P"])
    cs = compute_gpas(coef, pca, np.ones((2, 3)), ["G1", "G2"],
                      {"P": ["G1", "G2"]}, PagwasConfig(remove_outlier=False))
    assert cs.gpas.shape == (3,)
    assert list(cs.gpas.index) == ["c1", "c2", "c3"]
    assert np.allclose(cs.gpas.to_numpy(), [1.0, 2.0, 3.0])


def test_compute_gpas_winsorizes_by_default():
    cells = pd.Index([f"c{i}" for i in range(100)])
    coef = pd.DataFrame([np.arange(100.0)], index=["P"], columns=cells)
    pca = pd.DataFrame([np.ones(100)], index=["P"], columns=cells)
    cs = compute_gpas(coef, pca, np.ones((1, 100)), ["G1"],
                      {"P": ["G1"]}, PagwasConfig(remove_outlier=True))
    assert cs.gpas.max() < 99.0
