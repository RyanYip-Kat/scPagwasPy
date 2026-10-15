import anndata as ad
import numpy as np
import pandas as pd

from scpagwas.config import PagwasConfig
from scpagwas.result import PagwasResult


def _result():
    return PagwasResult(
        adata=ad.AnnData(X=np.array([[1.0]])),
        config=PagwasConfig(),
        celltype_bootstrap=pd.DataFrame({"bp_value": [0.1]}),
        gene_pcc=pd.Series({"g": 0.5}),
        pathway_cell_lm=pd.DataFrame({"c": [1.0]}, index=["P"]),
        cell_scores=pd.DataFrame({"scPagwas.gPAS.score": [1.0]}),
        merged_celltype_p=pd.DataFrame({"celltype": ["A"], "pvalue": [0.2]}),
        pathway_list={"P": ["g"]},
        pca_cell_df=pd.DataFrame({"A": [1.0]}, index=["P"]),
        lm_results={"parameters": pd.Series({"Intercept": 1.0})},
    )


def test_to_files_writes_five(tmp_path):
    files = _result().to_files("t", tmp_path)
    assert len(files) == 5
    assert all(f.exists() for f in files)
    names = sorted(f.name for f in files)
    assert "t_gene_PCC.csv" in names
    assert "t_Merged_celltype_pvalue.csv" in names
    assert "t_Pathway_singlecell_lm_results.txt" in names


def test_to_files_creates_outdir(tmp_path):
    outdir = tmp_path / "nested" / "out"
    _result().to_files("x", outdir)
    assert outdir.is_dir()
