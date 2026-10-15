# scPagwasPy

Python reimplementation of the R package **scPagwas** — a pathway-based
polygenic integration method that links GWAS summary statistics to scRNA-seq
data in order to prioritize disease-relevant cell subpopulations and genes.

The implementation is **functionally aligned** with the R package (same inputs,
same pipeline stages, same result files), but it is not bit-identical:
R-specific internals (`RMTstat::qWishartMax`, `Seurat::AddModuleScore`,
`bigstatsr::big_univLinReg`, Storey `qvalue`) are replaced by equivalent
Python implementations. See [Deviations](#deviations-from-the-r-package).

- Input: `AnnData` + GWAS summary statistics + pathway gene sets + gene
  coordinates + LD
- Output: per-cell scores in `adata.obs`, five result files, and eight figures

---

## Contents

- [Installation](#installation)
- [Inputs](#inputs)
- [Quick start](#quick-start)
- [Configuration](#configuration)
- [Bundled assets](#bundled-assets)
- [Outputs](#outputs)
- [Plots](#plots)
- [End-to-end example: PBMC3K](#end-to-end-example-pbmc3k)
- [How the pipeline works](#how-the-pipeline-works)
- [Deviations from the R package](#deviations-from-the-r-package)
- [Development](#development)

---

## Installation

Requires Python >= 3.11.

```bash
git clone <this-repo> scPagwasPy
cd scPagwasPy
uv venv
uv pip install -e .
```

Runtime dependencies: `numpy`, `scipy`, `pandas`, `anndata`, `statsmodels`,
`matplotlib`, `seaborn`, `pyarrow`. `pytest` is an optional extra
(`uv pip install -e ".[test]"`).

The bundled genomic assets (pathways, gene blocks, LD) make the installed
package roughly 70 MB.

The package uses a `src/` layout. To run without installing, set
`PYTHONPATH=src` (this is what `pytest.ini` does for the tests).

---

## Inputs

| Input | Type | Required columns / shape | Notes |
|---|---|---|---|
| Single-cell data | `anndata.AnnData` | `obs["cell_type"]`, log-normalized expression in `X` (or a `layers[...]` key) | Cell-type labels drive the pseudobulk step; use `Idents`/cluster labels |
| GWAS summary statistics | `pandas.DataFrame` | `chrom, pos, rsid, beta, se, maf` | Tab/space/comma-delimited files can be read with `read_gwas` |
| Pathway gene sets | `dict[str, list[str]]` | `{pathway: [gene, ...]}` | Gene **symbols**; read from GMT with `read_gmt` |
| Gene coordinates | `pandas.DataFrame` | `chrom, start, end, label` | `label` is the gene symbol; hg38/hg19 bundled |
| LD | `LDIndex` | per chromosome `SNP_A, SNP_B, R` | R² long table from PLINK; bundled for 1,000 Genomes |

Notes on the data:

- Expression should be **normalized** (e.g. `normalize_total` + `log1p`).
  Batch correction is recommended, and the matrix may be sparse.
- Chromosomes are normalized internally to integers `1..22`; `X/Y/M`
  become `23/24/26`, and any other contig is dropped.
- Duplicate gene symbols and duplicated genes inside a pathway list are
  deduplicated automatically.

---

## Quick start

```python
import anndata as ad

from scpagwas.config import PagwasConfig
from scpagwas.io import assets
from scpagwas.io.gwas import read_gwas
from scpagwas.pipeline import run

adata = ad.read_h5ad("my_data.h5ad")
adata.obs["cell_type"] = adata.obs["leiden"].astype(str)   # or your annotation

gwas = read_gwas("gwas_summary.txt")
pathways = assets.load_builtin_pathways("kegg")
blocks = assets.load_builtin_block_annotation("hg38")
ld = assets.load_builtin_ld()

cfg = PagwasConfig(n_topgenes=1000, iters_celltype=200, iters_singlecell=100)
result = run(adata, gwas, pathways, blocks, ld, config=cfg)

result.to_files("run1", "outputs/run1")      # five result files
result.adata                                 # AnnData with per-cell columns
```

`run(adata, gwas, pathways, block_annotation, ld, config=None, layer=None)` is
the single entry point. Pass `layer="counts"` to use a named layer instead of
`X`.

---

## Configuration

All behaviour is controlled by the `PagwasConfig` dataclass
(`scpagwas/config.py`). Defaults mirror `scPagwas::scPagwas_main`.

| Parameter | Default | Meaning |
|---|---|---|
| `min_clustercells` | `10` | Drop cell types with fewer cells than this |
| `min_pathway_size` | `5` | Drop pathways with fewer genes than this |
| `max_pathway_size` | `300` | Drop pathways with more genes than this |
| `min_shared_genes` | `100` | Error if the expression matrix and pathway sets share fewer genes |
| `marg` | `10000` | SNP-to-gene window: TSS ± `marg` bp (±10 kb by default) |
| `maf_filter` | `0.01` | Drop SNPs with MAF below this (applied after flipping MAF > 0.5) |
| `sex_filter` | `True` | Drop SNPs on sex chromosomes |
| `mhc_filter` | `True` | Drop SNPs in the MHC region (chr6:25–34 Mb) |
| `gwas_z_filter` | `None` | If set, drop SNPs with `abs(beta/se)` above this value |
| `n_topgenes` | `1000` | Number of top-correlated genes used for the TRS module score |
| `iters_celltype` | `200` | Bootstrap iterations for the cell-type regression |
| `iters_singlecell` | `100` | Iterations for the per-cell empirical background p-value; `0` skips it (much faster, no cell p-values) |
| `remove_outlier` | `True` | Winsorize the per-cell gPAS score at the 1st/99th percentile |
| `seed` | `1234` | Random seed for the bootstrap and the background null |
| `n_cores` | `1` | Reserved for parallel regression |

Typical tuning:

```python
# Fast exploratory run: no per-cell p-values, fewer bootstraps
cfg = PagwasConfig(iters_singlecell=0, iters_celltype=50, n_topgenes=500)

# Full run
cfg = PagwasConfig(iters_singlecell=100, iters_celltype=200, n_topgenes=1000)
```

---

## Bundled assets

The configuration data from the R package is bundled under
`src/scpagwas/data/` and exposed through `scpagwas.io.assets`:

| Loader | Returns | Size |
|---|---|---|
| `load_builtin_pathways("kegg")` | 345 KEGG pathways | 0.2 MB |
| `load_builtin_pathways("gobp")` | 7,658 GO biological-process sets | 4.2 MB |
| `load_builtin_pathways("reactome")` | 1,615 Reactome sets | 0.6 MB |
| `load_builtin_block_annotation("hg38")` | 62,700 gene blocks | 0.8 MB |
| `load_builtin_block_annotation("hg37")` | 64,486 gene blocks (GRCh37) | 0.8 MB |
| `load_builtin_ld()` | 6,789,678 LD pairs, 22 chromosomes | 62 MB |
| `example_gwas_path()` | example GWAS summary file | 0.6 MB |

`load_builtin_ld(chromosomes=[1, 2, ...])` restricts which chromosomes are read.

These were exported from the R package with
[`scripts/port_r_assets.R`](scripts/port_r_assets.R):

```bash
Rscript scripts/port_r_assets.R /path/to/scPagwas/data src/scpagwas/data
```

User-supplied data uses the same loaders, e.g.
`read_gmt("my_pathways.gmt")`, `load_block_annotation("blocks.parquet")`, and
`LDIndex.from_frames({1: df_chr1, 2: df_chr2, ...})`.

---

## Outputs

### Return value

`run()` returns a `PagwasResult`:

| Field | Content |
|---|---|
| `adata` | Copy of the input `AnnData` with results added (see below) |
| `config` | The `PagwasConfig` used |
| `celltype_bootstrap` | Per-cell-type bootstrap estimates, errors, p-values, 95% CIs |
| `gene_pcc` | Pearson correlation of each gene with the gPAS score (a `Series`) |
| `pathway_cell_lm` | Single-cell regression coefficients, pathways × cells |
| `cell_scores` | Per-cell score/p-value table |
| `merged_celltype_p` | Cell-type p-values merged from the per-cell p-values |
| `pathway_list` | The filtered pathway gene sets actually used |
| `pca_cell_df` | Pathway activity scores per cell type |
| `lm_results` | Cell-type level regression parameters |

Added to `adata`:

- `obs["scPagwas.TRS.Score"]` — enrichment score of the top-correlated genes
- `obs["scPagwas.gPAS.score"]` — per-cell genetic pathway association score
- `obs["Random_Correct_BG_p"]`, `obs["Random_Correct_BG_adjp"]`,
  `obs["Random_Correct_BG_z"]` — empirical p-value, FDR and z-score
- `obsm["scPagwas_pa_pca"]` — pathway activity matrix (cells × pathways)

### Result files (`PagwasResult.to_files(prefix, outdir)`)

| File | Content |
|---|---|
| `{prefix}_celltypes_bootstrap_results.csv` | `annotation, bootstrap_estimate, bootstrap_error, bt_value, bp_value, bias_corrected_estimate, CI_lo, CI_hi` |
| `{prefix}_gene_PCC.csv` | `PCC` per gene, sorted-able |
| `{prefix}_Pathway_singlecell_lm_results.txt` | Pathways × cells coefficient matrix |
| `{prefix}_singlecell_scPagwas_score_pvalue.Result.csv` | `scPagwas.TRS.Score, scPagwas.gPAS.score, Random_Correct_BG_p, _adjp, _z` |
| `{prefix}_Merged_celltype_pvalue.csv` | `celltype, pvalue` |

---

## Plots

Each function returns a matplotlib `Figure` and writes a file when `save_path`
is given. They mirror the R plotting functions one-to-one.

| R function | Python | Signature |
|---|---|---|
| `scPagwas_Visualization` | `plot_umap_scores` | `(adata, key="scPagwas.gPAS.score", obsm="X_umap", save_path=None)` |
| `scPagwas_Visualization` (significant cells) | `plot_significant_cells` | `(adata, p_thre=0.05, obsm="X_umap", save_path=None, key="Random_Correct_BG_adjp")` |
| `Bootstrap_estimate_Plot` | `plot_bootstrap_forest` | `(bootstrap_df, save_path=None)` |
| `Bootstrap_P_Barplot` | `plot_bootstrap_bar` | `(bootstrap_df, save_path=None)` |
| `heritability_cor_scatterplot` | `plot_gene_pcc` | `(gene_pcc, topn=10, save_path=None)` |
| `plot_bar_positie_nagtive` | `plot_positive_cells_bar` | `(adata, cell_type_key="cell_type", p_thre=0.05, save_path=None)` |
| `plot_vln_Corgenes` | `plot_top_gene_violin` | `(adata, genes, cell_type_key="cell_type", save_path=None)` |

```python
from scpagwas.plot.umap import plot_umap_scores, plot_significant_cells
from scpagwas.plot.forest import plot_bootstrap_forest
from scpagwas.plot.bar import plot_bootstrap_bar, plot_positive_cells_bar
from scpagwas.plot.gene import plot_gene_pcc
from scpagwas.plot.violin import plot_top_gene_violin
```

---

## End-to-end example: PBMC3K

[`scripts/run_pbmc3k.py`](scripts/run_pbmc3k.py) runs the whole pipeline on
`testdata/pbmc3k_data.h5ad` with the bundled KEGG pathways, hg38 blocks and LD:

```bash
PYTHONPATH=src .venv/bin/python scripts/run_pbmc3k.py
```

Everything below is the actual output of that script.

### Dataset and run

| Item | Value |
|---|---|
| Cells | 2,700 |
| Genes | 32,738 |
| Cell types | 9 (Leiden clusters used as `cell_type`) |
| GWAS SNPs (after cleaning) | 9,981 |
| KEGG pathways | 345 → 321 after the size filter → 313 scored |
| Pathways linked to SNPs | 306 |
| Single-cell coefficient matrix | 306 pathways × 2,700 cells |
| Significant cells (FDR < 0.05) | 685 (25.4%) |
| Wall-clock time | ~162 s (48 cores available, single-threaded steps) |

### Per-cell-type summary

| Cluster | Cells | Mean gPAS | Mean TRS | Significant cells | % significant |
|---:|---:|---:|---:|---:|---:|
| 0 | 629 | −1.50 | −0.17 | 0 | 0.0% |
| 1 | 544 | −1.45 | −0.14 | 0 | 0.0% |
| **2** | 480 | −0.25 | 0.48 | **480** | **100.0%** |
| 3 | 350 | −2.17 | −0.15 | 0 | 0.0% |
| 4 | 315 | −1.55 | −0.17 | 0 | 0.0% |
| **5** | 175 | −0.82 | 0.44 | **174** | **99.4%** |
| 6 | 161 | −2.31 | −0.16 | 0 | 0.0% |
| **7** | 33 | 0.16 | 0.28 | **30** | **90.9%** |
| 8 | 13 | −4.83 | −0.10 | 1 | 7.7% |
| **Total** | **2,700** | −1.36 | — | **685** | **25.4%** |

### Cell-type bootstrap (`*_celltypes_bootstrap_results.csv`)

`bp_value` is the one-sided bootstrap p-value
(`pnorm(estimate/error, lower.tail = FALSE)`).

| Cell type | Bias-corrected estimate | 95% CI | bt_value | bp_value |
|---|---:|---|---:|---:|
| Intercept | −0.0035 | [−0.0084, 0.0071] | −0.35 | 0.637 |
| 0 | 0.0030 | [−0.0042, 0.0069] | 0.49 | 0.313 |
| 1 | −0.0049 | [−0.0097, −0.0000] | −2.01 | 0.978 |
| **2** | **0.0059** | [−0.0002, 0.0107] | 2.07 | **0.019** |
| 3 | −0.0021 | [−0.0069, 0.0011] | −1.22 | 0.888 |
| **4** | **0.0043** | [0.0004, 0.0069] | 2.13 | **0.017** |
| 5 | −0.0011 | [−0.0061, 0.0027] | −0.63 | 0.735 |
| 6 | −0.0018 | [−0.0046, −0.0006] | −2.17 | 0.985 |
| **7** | **0.0260** | [0.0198, 0.0375] | 5.99 | **~1e−9** |
| 8 | −0.0023 | [−0.0052, −0.0012] | −2.58 | 0.995 |

### Merged cell-type p-values (`*_Merged_celltype_pvalue.csv`)

Stouffer-Z merge of the per-cell empirical p-values within each cell type.

| Cluster | p-value | Cluster | p-value | Cluster | p-value |
|---|---:|---|---:|---|---:|
| 0 | 0.498 | 3 | 0.495 | 6 | 0.494 |
| 1 | 0.496 | **4** | 0.499 | **7** | **1.7e−8** |
| **2** | **1.6e−10** | **5** | **1.9e−10** | 8 | 0.312 |

The two cell-type strategies mostly agree (clusters **2**, **7** stand out in
both). They can differ on borderline clusters — bootstrap flags cluster 4
(p = 0.017) while the single-cell merge flags cluster 5 (p = 1.9e−10) — because
bootstrap works on pseudobulk data while the merge aggregates single-cell
results. This mirrors the R package's own note about the two strategies.

### Top genes correlated with gPAS (`*_gene_PCC.csv`)

| Rank | Gene | PCC | Rank | Gene | PCC |
|---:|---|---:|---:|---|---:|
| 1 | ALDH2 | 0.379 | 7 | FCN1 | 0.191 |
| 2 | CFP | 0.207 | 8 | CFD | 0.189 |
| 3 | LYZ | 0.200 | 9 | CST3 | 0.189 |
| 4 | AIF1 | 0.197 | 10 | LST1 | 0.187 |
| 5 | LGALS2 | 0.194 | 11 | S100A9 | 0.186 |
| 6 | CALM2 | 0.194 | 12 | CD14 | 0.185 |

These are canonical **monocyte/myeloid** markers (LYZ, AIF1, FCN1, CFD, CFP,
CST3, S100A9, CD14) — the expected biology for PBMC3K.

### Figures

**gPAS score and TRS score on the UMAP** — the myeloid/monocyte cluster at the
bottom right carries the highest score.

![scPagwas gPAS score UMAP](docs/figures/scPagwas.gPAS.score_umap.png)

![scPagwas TRS score UMAP](docs/figures/scPagwas.TRS.Score_umap.png)

**Significant cells** (empirical adjusted p <= 0.05) — the monocyte cluster and
the small platelet cluster are highlighted.

![significant cells UMAP](docs/figures/scPagwas_Correct_BG_adjp0.05_umap.png)

**Cell-type bootstrap forest plot** — cluster 7 stands out with a tight CI.

![bootstrap forest plot](docs/figures/Bootstrap_estimate_forest.png)

**Cell-type bootstrap -log2(p) bar plot**

![bootstrap p-value bar plot](docs/figures/Bootstrap_P_barplot.png)

**Gene PCC scatter** (top-ranked genes labelled)

![gene PCC scatter](docs/figures/gene_PCC.png)

**Positive-cell fraction per cluster**

![positive cell fraction](docs/figures/positive_cells_bar.png)

**Top correlated genes per cluster** (violin grid)

![top gene violins](docs/figures/top_gene_violin.png)

---

## How the pipeline works

1. **Prepare single cells** — drop small cell types, keep the normalized
   expression matrix (genes × cells) and build a pseudobulk mean per cell type.
2. **Score pathway activity** — for each pathway and cell type, take PC1 of the
   centred expression submatrix and standardize its variance by the analytic
   Marchenko–Pastur upper edge (signal / noise).
3. **Clean GWAS** — MAF flip/filter, deduplicate, drop sex chromosomes and MHC.
4. **Map SNPs to genes** — each SNP is assigned to genes whose TSS window
   (TSS ± `marg`) contains it.
5. **Build pathway blocks** — pathway genes → coordinate blocks, grouped by
   chromosome.
6. **Link pathway blocks to GWAS** and **regress** `beta^2 ~ offset(se^2) + x`,
   where `x` is LD-weighted (R⁴) pathway-activity × expression. This yields a
   per-cell-type contribution and, at single-cell level, a coefficient per cell.
7. **gPAS score** — combine the per-pathway coefficients with pathway activity
   and mean expression, summed over pathways for each cell.
8. **Gene PCC** — Pearson correlation of each gene with the gPAS score.
9. **Significance** — module score of the top `n_topgenes` genes, corrected
   against a randomized expression-matched background, then merged per cell
   type with a Stouffer-Z test.

---

## Deviations from the R package

Functionally equivalent, numerically different:

- **Null standardization** uses the analytic Marchenko–Pastur edge instead of
  `RMTstat::qWishartMax`.
- **Module score** uses expression-binned control genes (24 bins, 100 controls)
  instead of `Seurat::AddModuleScore`.
- **Single-cell regression** is vectorized (`beta = cov(x, y - offset) / var(x)`)
  instead of `bigstatsr::big_univLinReg`.
- **LD lookup** is restricted to SNPs present in the GWAS before building the
  matrix (same pairs, faster).

Known gaps:

- `bp_value` in the cell-type bootstrap is **not** FDR-adjusted (the R package
  applies `p.adjust(method = "fdr")`).
- The v2 additions are not implemented: the `FindAllMarkers` cell-type pathway
  p-value and the randomized weighted gene PCC (`scGet_PCC2`).
- Large-data helpers (`run_split` / `Merge_gPas` / `Corr_Random`) are not ported.

---

## Development

```bash
uv venv
uv pip install -p .venv numpy scipy pandas anndata statsmodels \
    matplotlib seaborn pyarrow pytest pyreadr
MPLCONFIGDIR=/tmp/mplconfig .venv/bin/python -m pytest -p no:cacheprovider -q
```

`MPLCONFIGDIR` is needed because the plotting tests use the headless `Agg`
backend. Design and plan documents live in
`docs/superpowers/specs/` and `docs/superpowers/plans/`.
