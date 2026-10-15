"""Configuration for a scpagwas run."""

from dataclasses import dataclass


@dataclass
class PagwasConfig:
    """Parameters controlling the scpagwas pipeline.

    Defaults mirror the R package ``scPagwas`` (``scPagwas_main``).
    """

    min_clustercells: int = 10
    min_pathway_size: int = 5
    max_pathway_size: int = 300
    min_shared_genes: int = 100
    marg: int = 10_000
    maf_filter: float = 0.01
    sex_filter: bool = True
    mhc_filter: bool = True
    gwas_z_filter: float | None = None
    n_topgenes: int = 1000
    iters_celltype: int = 200
    iters_singlecell: int = 100
    remove_outlier: bool = True
    seed: int = 1234
    n_cores: int = 1
