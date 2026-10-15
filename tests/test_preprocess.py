import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import warnings

from scpagwas.config import PagwasConfig
from scpagwas.preprocess import prepare_singlecell


def _adata():
    rng = np.random.default_rng(0)
    X = sp.csr_matrix(rng.random((30, 10)))
    obs = pd.DataFrame(
        {"cell_type": ["A"] * 20 + ["B"] * 10},
        index=[f"c{i}" for i in range(30)],
    )
    var = pd.DataFrame(index=[f"G{i}" for i in range(10)])
    return ad.AnnData(X=X, obs=obs, var=var)


def _pw():
    return {f"p{i}": [f"G{j}" for j in range(10)] for i in range(12)}


def test_prepare_shapes():
    cfg = PagwasConfig(min_clustercells=5, min_shared_genes=5)
    sc = prepare_singlecell(_adata(), _pw(), cfg)
    assert sc.data_mat.shape == (10, 30)
    assert list(sc.merge_expr.columns) == ["A", "B"]


def test_prepare_drops_small_clusters():
    cfg = PagwasConfig(min_clustercells=15, min_shared_genes=5)
    sc = prepare_singlecell(_adata(), _pw(), cfg)
    assert sc.celltypes == ["A"]


def test_prepare_requires_overlap():
    try:
        prepare_singlecell(_adata(), {"p": ["ZZZ"]}, PagwasConfig(min_shared_genes=5))
        assert False
    except ValueError:
        pass


def test_prepare_dedupes_gene_names():
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(10)] + ["G0"]  # duplicate symbol
    obs = pd.DataFrame({"cell_type": ["A"] * 20 + ["B"] * 10},
                       index=[f"c{i}" for i in range(30)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = ad.AnnData(X=sp.csr_matrix(rng.random((30, 11))), obs=obs,
                       var=pd.DataFrame(index=genes))
    pw = {f"p{i}": [f"G{j}" for j in range(10)] for i in range(12)}
    sc = prepare_singlecell(a, pw, PagwasConfig(min_clustercells=5, min_shared_genes=5))
    assert sc.genes.is_unique
    assert sc.data_mat.shape[0] == 10
