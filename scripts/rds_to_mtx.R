# Convert a GEO count-matrix RDS into MatrixMarket + row/column name files.
# Usage: Rscript scripts/rds_to_mtx.R <input.RDS> <output_dir>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript rds_to_mtx.R <input.RDS> <output_dir>")
suppressPackageStartupMessages(library(Matrix))

x <- readRDS(args[1])
cat("Object class:", class(x), "\n")

if (inherits(x, "Seurat")) {
  # Requires install.packages("SeuratObject")
  suppressPackageStartupMessages(library(SeuratObject))
  x <- SeuratObject::LayerData(x, assay = SeuratObject::DefaultAssay(x), layer = "counts")
}
if (is.data.frame(x)) x <- as.matrix(x)
if (!inherits(x, "dgCMatrix")) x <- as(x, "CsparseMatrix")
if (is.null(rownames(x)) || is.null(colnames(x))) stop("Matrix is missing row or column names")

dir.create(args[2], showWarnings = FALSE, recursive = TRUE)
writeMM(x, file.path(args[2], "matrix.mtx"))
writeLines(rownames(x), file.path(args[2], "rows.txt"))
writeLines(colnames(x), file.path(args[2], "cols.txt"))
cat("Dimensions (rows x cols):", dim(x), "\n")
cat("First row names:", head(rownames(x), 3), "\n")
cat("First col names:", head(colnames(x), 3), "\n")
