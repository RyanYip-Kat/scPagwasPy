#!/usr/bin/env Rscript
# Export scPagwas R package assets into formats the Python package reads.
#
# Run this once in an environment that has the R `scPagwas` package installed
# (plus `arrow`, `Seurat`, `SeuratDisk`). It writes into ./scpagwas_assets/.
#
#   Rscript scripts/export_r_assets.R
#
# Outputs:
#   block_annotation.parquet        gene coordinates (GRCh38)
#   block_annotation_hg37.parquet   gene coordinates (GRCh37)
#   chrom_ld_<chrom>.parquet        per-chromosome LD (SNP_A, SNP_B, R)
#   pathways_<name>.gmt             pathway gene sets
#   scRNAexample.h5ad               example single-cell object

suppressMessages({
  library(scPagwas)
  library(arrow)
})

outdir <- "scpagwas_assets"
dir.create(outdir, showWarnings = FALSE)

# --- gene coordinate blocks -------------------------------------------------
data(block_annotation)
write_parquet(block_annotation[, c("chrom", "start", "end", "label")],
              file.path(outdir, "block_annotation.parquet"))
data(block_annotation_hg37)
write_parquet(block_annotation_hg37[, c("chrom", "start", "end", "label")],
              file.path(outdir, "block_annotation_hg37.parquet"))

# --- LD per chromosome ------------------------------------------------------
data(chrom_ld)
for (i in seq_along(chrom_ld)) {
  write_parquet(as.data.frame(chrom_ld[[i]]),
                file.path(outdir, paste0("chrom_ld_", i, ".parquet")))
}

# --- pathway gene sets (GMT) ------------------------------------------------
write_gmt <- function(pathways, file) {
  lines <- vapply(names(pathways), function(nm) {
    paste(c(nm, "na", pathways[[nm]]), collapse = "\t")
  }, character(1))
  writeLines(lines, file)
}
data(Genes_by_pathway_kegg)
write_gmt(Genes_by_pathway_kegg, file.path(outdir, "pathways_kegg.gmt"))
data(genes.by.gobp.pathway)
write_gmt(genes.by.gobp.pathway, file.path(outdir, "pathways_gobp.gmt"))
data(genes.by.reactome.pathway)
write_gmt(genes.by.reactome.pathway, file.path(outdir, "pathways_reactome.gmt"))

# --- example single-cell object --------------------------------------------
Single_data <- readRDS(system.file("extdata", "scRNAexample.rds",
                                   package = "scPagwas"))
if (requireNamespace("SeuratDisk", quietly = TRUE)) {
  h5seurat <- file.path(outdir, "scRNAexample.h5seurat")
  SeuratDisk::SaveH5Seurat(Single_data, filename = h5seurat, overwrite = TRUE)
  SeuratDisk::Convert(h5seurat, dest = "h5ad", overwrite = TRUE)
}

message("Assets written to ", normalizePath(outdir))
