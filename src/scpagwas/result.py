"""Result container returned by :func:`scpagwas.pipeline.run`."""

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd


@dataclass
class PagwasResult:
    """Bundle of all scpagwas outputs.

    ``adata`` carries the per-cell results in ``obs``/``obsm``; the DataFrames
    are the tabular results also written by :meth:`to_files`.
    """

    adata: object
    config: object
    celltype_bootstrap: pd.DataFrame
    gene_pcc: pd.Series
    pathway_cell_lm: pd.DataFrame
    cell_scores: pd.DataFrame
    merged_celltype_p: pd.DataFrame
    pathway_list: dict
    pca_cell_df: pd.DataFrame
    lm_results: dict = field(default_factory=dict)

    def to_files(self, prefix, outdir):
        """Write the five scPagwas result files; return their paths."""
        out = Path(outdir)
        out.mkdir(parents=True, exist_ok=True)
        written = []
        pairs = [
            (f"{prefix}_celltypes_bootstrap_results.csv", self.celltype_bootstrap, {}),
            (f"{prefix}_gene_PCC.csv", self.gene_pcc.rename("PCC").to_frame(), {}),
            (f"{prefix}_Pathway_singlecell_lm_results.txt", self.pathway_cell_lm,
             {"sep": "\t"}),
            (f"{prefix}_singlecell_scPagwas_score_pvalue.Result.csv", self.cell_scores,
             {}),
            (f"{prefix}_Merged_celltype_pvalue.csv", self.merged_celltype_p,
             {"index": False}),
        ]
        for name, df, kw in pairs:
            p = out / name
            df.to_csv(p, **kw)
            written.append(p)
        return written
