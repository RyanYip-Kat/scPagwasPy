#!/usr/bin/env Rscript
# Port scPagwas R package data assets into the Python package's data dir.
#
# Usage:
#   Rscript scripts/port_r_assets.R [rdata_dir] [out_dir]
#
# Defaults assume the R package lives at ../scPagwas relative to this repo.

args <- commandArgs(trailingOnly = TRUE)
rdata_dir <- if (length(args) >= 1) args[1] else "../scPagwas/data"
extdata_dir <- file.path(dirname(rdata_dir), "inst", "extdata")
out_dir <- if (length(args) >= 2) args[2] else file.path("src", "scpagwas", "data")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ld_dir <- file.path(out_dir, "chrom_ld")
dir.create(ld_dir, recursive = TRUE, showWarnings = FALSE)

load_first <- function(path) {
  e <- new.env()
  load(path, envir = e)
  e[[ls(e)[1]]]
}

write_gmt <- function(pathways, file) {
  lines <- vapply(names(pathways), function(nm) {
    genes <- pathways[[nm]]
    paste(c(nm, "na", genes), collapse = "\t")
  }, character(1))
  writeLines(lines, file)
}

message("* block annotation")
anno <- load_first(file.path(rdata_dir, "block_annotation.RData"))
data.table::fwrite(anno, file.path(out_dir, "block_annotation.tsv.gz"),
                   sep = "\t", compress = "gzip")
anno37 <- load_first(file.path(rdata_dir, "block_annotation_hg37.RData"))
data.table::fwrite(anno37, file.path(out_dir, "block_annotation_hg37.tsv.gz"),
                   sep = "\t", compress = "gzip")

message("* pathway gene sets")
write_gmt(load_first(file.path(rdata_dir, "Genes_by_pathway_kegg.RData")),
          file.path(out_dir, "pathways_kegg.gmt"))
if (file.exists(file.path(rdata_dir, "genes.by.gobp.pathway.RData"))) {
  write_gmt(load_first(file.path(rdata_dir, "genes.by.gobp.pathway.RData")),
            file.path(out_dir, "pathways_gobp.gmt"))
}
if (file.exists(file.path(rdata_dir, "genes.by.reactome.pathway.RData"))) {
  write_gmt(load_first(file.path(rdata_dir, "genes.by.reactome.pathway.RData")),
            file.path(out_dir, "pathways_reactome.gmt"))
}

message("* LD per chromosome")
ld <- load_first(file.path(rdata_dir, "chrom_ld.RData"))
for (i in seq_along(ld)) {
  data.table::fwrite(as.data.frame(ld[[i]]),
                     file.path(ld_dir, sprintf("chr%d.tsv.gz", i)),
                     sep = "\t", compress = "gzip")
}

message("* example GWAS")
gwas_example <- file.path(extdata_dir, "GWAS_summ_example.txt")
if (file.exists(gwas_example)) {
  file.copy(gwas_example, file.path(out_dir, "gwas_example.txt"), overwrite = TRUE)
}

message("Assets written to ", normalizePath(out_dir))
