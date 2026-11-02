# Run with Rscript or source(); locate the repository independently of the working directory.
repo_root <- local({
  candidates <- c(Sys.getenv("ELNMF_ROOT", unset = ""), getwd())
  for (frame in sys.frames()) {
    if (!is.null(frame$ofile)) candidates <- c(dirname(frame$ofile), candidates)
  }
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_arg)) candidates <- c(dirname(sub("^--file=", "", script_arg[[1]])), candidates)
  found <- NULL
  for (candidate in candidates[nzchar(candidates)]) {
    current <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
    repeat {
      if (dir.exists(file.path(current, "main", "shared")) && dir.exists(file.path(current, "data"))) {
        found <- current; break
      }
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
    if (!is.null(found)) break
  }
  if (is.null(found)) stop("Cannot locate ELNMF repository. Run from its root or set ELNMF_ROOT.")
  found
})
setwd(repo_root)
# Outputs are created only when this script is run, outside main/.

suppressPackageStartupMessages({
  library(Seurat)
})

data_dir <- "data/Human_Breast_Cancer"
output_dir <- Sys.getenv("HBC_BENCH_OUTPUT", unset = "results/HBC_baseline_computational_benchmark/labels")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Input preparation and spot alignment are outside the timed region.
hbc <- Load10X_Spatial(
  data.dir = data_dir,
  filename = "filtered_feature_bc_matrix.h5",
  slice = "HBC"
)
meta <- read.delim(file.path(data_dir, "metadata.tsv"), row.names = 1,
                   stringsAsFactors = FALSE, check.names = FALSE)
coord <- read.csv(file.path(data_dir, "spatial/tissue_positions_list.csv"),
                  header = FALSE, stringsAsFactors = FALSE)
colnames(coord) <- c("barcode", "in_tissue", "array_row", "array_col", "x", "y")
valid_ids <- rownames(meta)[!is.na(meta$ground_truth) & nzchar(meta$ground_truth)]
tissue_ids <- coord$barcode[coord$in_tissue == 1]
final_spots <- Cells(hbc)[Cells(hbc) %in% valid_ids & Cells(hbc) %in% tissue_ids]
if (length(final_spots) != 3798L) stop("Expected 3798 HBC spots, found ", length(final_spots))
hbc <- subset(hbc, cells = final_spots)
hbc@meta.data <- cbind(hbc@meta.data, meta[Cells(hbc), , drop = FALSE])

set.seed(1234)
cat("BENCHMARK_START|Seurat|spots=3798|target_clusters=20\n")
flush.console()
timer_start <- proc.time()[["elapsed"]]

options(future.globals.maxSize = 800 * 1024^2)
hbc <- SCTransform(hbc, assay = "Spatial", verbose = FALSE)
hbc <- RunPCA(hbc, assay = "SCT", verbose = FALSE)
hbc <- FindNeighbors(hbc, reduction = "pca", dims = 1:40, verbose = FALSE)

# The supplied HBC setting gives 18 clusters at resolution 1.21.  A nearby,
# predeclared resolution is used to approach the requested target of 20 without
# running a performance-based search during the benchmark.
seurat_resolution <- 1.35
hbc <- FindClusters(hbc, resolution = seurat_resolution,
                    random.seed = 1234, verbose = FALSE)
labels <- hbc$seurat_clusters

runtime_sec <- proc.time()[["elapsed"]] - timer_start
cat(sprintf("BENCHMARK_END|Seurat|runtime_sec=%.6f|clusters=%d|resolution=%.2f\n",
            runtime_sec, length(unique(labels)), seurat_resolution))
flush.console()

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(Barcode = Cells(hbc), Cluster = as.character(labels)),
          file.path(output_dir, "HBC_Seurat_labels.csv"), row.names = FALSE)

