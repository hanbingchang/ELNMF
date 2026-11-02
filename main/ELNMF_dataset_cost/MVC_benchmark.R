# Historical cost run: q = 0.20 (not the Table S1 value 0.30).
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

getwd()

# 加载必要的R包 ------------------------------------------------------------
library(Seurat)
library(hdf5r)
library(ggplot2)
library(patchwork)
library(Matrix)  
library(mclust)    # 提供 NMI 计算（adjusted = TRUE）
library(aricode)   # 提供高效 ARI/NMI（推荐）

# 设置随机种子以保证结果可重复
set.seed(1234)

#------------------------------读取数据----------------------------------------
cat("========== 数据加载 ==========\n")
data_path <- "data/MVC/filtered_feature_bc_matrix.h5"
st_data <- Read10X_h5(data_path)
seurat_obj <- CreateSeuratObject(
  counts = st_data,
  project = "MVC",
  min.cells = 0,
  min.features = 0
)
cat("原始数据维度（genes × spots）:", dim(seurat_obj), "\n")

#------------------------------低质量过滤----------------------------------------
cat("========== 低质量过滤 ==========\n")
counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
n_spots <- ncol(counts_mat)
zero_ratio <- Matrix::rowSums(counts_mat == 0) / n_spots
tau <- 0.86
genes_keep <- names(zero_ratio)[zero_ratio < tau]
seurat_obj <- subset(seurat_obj, features = genes_keep)
cat("零表达比例过滤后基因数:", nrow(seurat_obj), "\n")

#------------------------------归一化-------------------------------------------
cat("========== 数据归一化 ==========\n")
counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
spot_totals <- Matrix::colSums(counts_mat)
spot_totals[spot_totals == 0] <- NA
scale_factor <- 1e4
norm_mat <- t(t(counts_mat) / spot_totals) * scale_factor

#------------------------------对数变换-------------------------------------------
log_norm_mat <- log1p(norm_mat)
seurat_obj <- SetAssayData(
  object = seurat_obj,
  assay = "RNA",
  layer = "data",
  new.data = log_norm_mat
)

#------------------------------筛选高变异基因-------------------------------------------
cat("========== 筛选高变异基因 ==========\n")
seurat_obj <- FindVariableFeatures(seurat_obj,
                                   selection.method = "vst",
                                   nfeatures = 2000)
X <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")[VariableFeatures(seurat_obj), ]
cat("高变异基因表达矩阵维度:", dim(X), "\n")
#------------------------------数据对齐-------------------------------------------
#---------------- 关键：统一斑点过滤（组织内+有标签）
# ============================================

# 1. 读取空间坐标文件，确定哪些斑点在组织内
coord_file <- "data/MVC/tissue_positions_list.csv"
coord_df <- read.csv(coord_file, header = FALSE, stringsAsFactors = FALSE)
colnames(coord_df) <- c("barcode", "in_tissue", "array_row", "array_col", "imagerow", "imagecol")

tissue_spots <- coord_df$barcode[coord_df$in_tissue == 1]
cat(sprintf("组织内斑点数: %d\n", length(tissue_spots)))

# 2. 读取metadata，文件无表头，共5列：细胞ID、Total_couX、X、Y、label
meta_path <- "data/MVC/metadata.csv"
metadata <- read.csv(meta_path, header = TRUE, stringsAsFactors = FALSE)

# 预览前几行确认
head(metadata)
# 查看列名，确保正确
colnames(metadata)  # 应该显示 "ID", "Total_couX", "X", "Y", "label"
# 创建标签映射（ID -> label）
label_map <- setNames(metadata$label, metadata$ID)

# 处理空白或各种NA形式
na_patterns <- c("", "NA", "na", "Na", "nA")
for (pattern in na_patterns) {
  label_map[label_map == pattern] <- NA
}

# 只保留非NA的标签
valid_labels <- label_map[!is.na(label_map)]
labeled_spots <- names(valid_labels)
cat(sprintf("有标签（非NA）的斑点数: %d\n", length(labeled_spots)))
# 3. 获取当前表达矩阵X中的所有斑点
all_spots <- colnames(X)
cat(sprintf("表达矩阵X中的斑点数: %d\n", length(all_spots)))

# 4. 取三个集合的交集：既要存在于X中，又要在组织内，还要有标签
keep_spots <- colnames(X) %in% tissue_spots & colnames(X) %in% labeled_spots


# 5. 统计信息
cat(sprintf("\n=== 斑点过滤统计 ===\n"))
cat(sprintf("表达矩阵原始斑点数: %d\n", length(all_spots)))
cat(sprintf("组织内斑点数: %d\n", length(tissue_spots)))
cat(sprintf("有标签（非NA）斑点数: %d\n", length(labeled_spots)))
cat(sprintf("最终保留的斑点数（组织内+有标签）: %d\n", length(keep_spots)))
cat(sprintf("过滤掉的斑点数: %d (%.1f%%)\n", 
            length(all_spots) - length(keep_spots),
            (length(all_spots) - length(keep_spots)) / length(all_spots) * 100))

# 6. 过滤表达矩阵X
X <- X[, keep_spots, drop = FALSE]
cat(sprintf("\n过滤后X维度: %d × %d\n", nrow(X), ncol(X)))


# 用 colnames(X) 对齐标签（关键！！！）
true_labels <- valid_labels[colnames(X)]
true_labels <- factor(true_labels)

cat(sprintf("真实类别: %s\n", paste(levels(true_labels), collapse=", ")))
cat(sprintf("真实类别簇数: %d\n", nlevels(true_labels)))



#------------------------------混合图构建-------------------------------------------
cat("========== 构建混合图 ==========\n")
source("main/shared/hybrid graph.R") 

hybrid_df <- hybrid_graph(
  coord_header = "legacy", # Preserve the reader used for this historical experiment.
  X = X,
  coord_file = "data/MVC/tissue_positions_list.csv",
  model = "KNN",
  k_cutoff = 27,
  alpha = 0.8,
  sigma = 1,
  verbose = TRUE
)

cat("========== 重构邻接矩阵 (基于 Barcode 严格对齐) ==========\n")

# 1. 获取主脚本中当前 X 的列名 (这是我们的金标准)
current_barcodes <- colnames(X)
n_spots <- length(current_barcodes)

# 2. 检查 hybrid_df 中的条码是否都在当前 X 中
all_barcodes_in_df <- unique(c(hybrid_df$Cell1_barcode, hybrid_df$Cell2_barcode))
missing_barcodes <- setdiff(all_barcodes_in_df, current_barcodes)

if (length(missing_barcodes) > 0) {
  warning(sprintf("发现 %d 个条码在混合图中存在但在当前表达矩阵 X 中缺失。将丢弃相关边。\n", length(missing_barcodes)))
  # 过滤掉包含缺失条码的边
  keep_rows <- !(hybrid_df$Cell1_barcode %in% missing_barcodes | hybrid_df$Cell2_barcode %in% missing_barcodes)
  hybrid_df <- hybrid_df[keep_rows, ]
}

# 3. 【关键步骤】使用 match 将 Barcode 映射为当前 X 的索引
# match 返回的是 barcode 在 current_barcodes 向量中的位置 (1 到 n_spots)
# 这保证了无论函数内部顺序如何，这里生成的索引都严格对应主脚本的 X
idx_1 <- match(hybrid_df$Cell1_barcode, current_barcodes)
idx_2 <- match(hybrid_df$Cell2_barcode, current_barcodes)

# 双重检查：确保没有 NA (理论上上面过滤后不应该有)
if (any(is.na(idx_1)) || any(is.na(idx_2))) {
  stop("致命错误：仍有条码无法匹配。请检查 barcode 字符串是否完全一致（注意空格、大小写）。")
}

# 4. 构建稀疏矩阵
# 维度强制设定为当前 X 的维度
W <- sparseMatrix(
  i = idx_1,
  j = idx_2,
  x = hybrid_df$Similarity,
  dims = c(n_spots, n_spots)
)

# 5. 对称化 (因为之前去重时可能只保留了一个方向，或者为了保险)
W <- (W + t(W)) / 2

# 6. 最终验证
cat(sprintf("邻接矩阵 W 维度: %d x %d\n", nrow(W), ncol(W)))
cat(sprintf("表达矩阵 X 维度: %d x %d (spots)\n", nrow(X), ncol(X)))

# 这里现在应该能顺利通过了
stopifnot(nrow(W) == ncol(X))
stopifnot(ncol(W) == ncol(X))

cat(">>> 矩阵维度对齐成功！\n")


#------------------------------设置参数-------------------------------------------
cat("========== 参数设置 ==========\n")

k_values <- seq(10, 28, by = 2)
lambda_values <- c(180, 185, 190, 195, 200)

max_iter <- 500
tol <- 1e-4
n_clusters <- length(unique(na.omit(true_labels)))

cat("k_values:", k_values, "\n")
cat("lambda_values:", lambda_values, "\n")

source("main/shared/NMF.R")
#------------------------------保存V矩阵-------------------------------------------
cat("\n========== 生成并保存所有V矩阵 ==========\n")

# 创建保存目录
save_dir <- "results/ELNMF_computational_cost/work/MVC"
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("创建目录: %s\n", save_dir))
}

# 准备存储结构
V_list <- list()
V_metadata <- data.frame()
counter <- 1

# 记录开始时间
start_time <- Sys.time()

# 循环生成所有V矩阵
for (lambda1 in lambda_values) {
  for (k in k_values) {
    cat(sprintf("生成 V[%d/%d]: λ=%d, k=%d\n", 
                counter, length(k_values) * length(lambda_values),
                lambda1, k))
    
    # 运行NMF
    nmf_res <- regularized_nmf(
      X = X,
      W = W,
      k = k,
      lambda1 = lambda1,
      max_iter = max_iter,
      tol = tol,
      verbose = FALSE
    )
    
    # 提取V矩阵
    V <- nmf_res$V
    
    # 创建唯一标识符
    v_name <- sprintf("V_lambda%d_k%d", lambda1, k)
    
    # 存储到列表
    V_list[[v_name]] <- V
    
    # 保存元数据
    V_metadata <- rbind(V_metadata, data.frame(
      id = v_name,
      lambda = lambda1,
      k = k,
      n_rows = nrow(V),
      n_cols = ncol(V),
      file_name = paste0(v_name, ".rds"),
      stringsAsFactors = FALSE
    ))
    
    # 单独保存每个V矩阵为RDS文件
    saveRDS(V, file.path(save_dir, paste0(v_name, ".rds")))
    
    counter <- counter + 1
  }
}

# 保存完整的V列表和元数据
saveRDS(V_list, file.path(save_dir, "all_V_matrices.rds"))
write.csv(V_metadata, 
          file.path(save_dir, "V_metadata.csv"),
          row.names = FALSE)

# 记录结束时间
end_time <- Sys.time()
time_diff <- difftime(end_time, start_time, units = "mins")

cat(sprintf("\n========== V矩阵保存完成 ==========\n"))
cat(sprintf("• 总共生成 %d 个V矩阵\n", length(V_list)))
cat(sprintf("• 保存目录: %s\n", save_dir))
cat(sprintf("• 主要文件: \n"))
cat(sprintf("  - all_V_matrices.rds (所有V矩阵的列表)\n"))
cat(sprintf("  - V_metadata.csv (元数据)\n"))
cat(sprintf("  - V_lambda*_k*.rds (每个V矩阵单独文件)\n"))
cat(sprintf("• 生成耗时: %.2f 分钟\n", time_diff))




# ==================== 加载方法函数 ====================

source("main/shared/ensemble_clustering_filtered.R")
source("main/shared/aggregative_indicator.R")
source("main/shared/filter_top.R")


# ==================  质量控制集成聚类 ==================
cat("========== 集成聚类（带筛选） ==========\n")

ensemble_res <- ensemble_clustering_filtered(
  X = X,
  V_dir = save_dir,
  k_values = k_values,
  lambda_values = lambda_values,
  n_clusters = n_clusters,
  clustering_methods = c("kmeans", "pam"),
  q = 0.20,          # ⭐ 你只需要改这个
  verbose = TRUE
)

ensemble_res$summary_before
ensemble_res$summary_after
kept_S_matrices <- ensemble_res$all_S_matrices  # 或直接命名为 all_S_matrices
if (length(kept_S_matrices) > 0) {
  S_consensus <- Reduce("+", kept_S_matrices) / length(kept_S_matrices)
  S_consensus <- as.matrix(S_consensus)
  cat("保留的高质量聚类个数:", ensemble_res$n_kept, "\n")
} else {
  stop("没有聚类结果通过质量筛选！请降低分位数 q")
  
}

# ==================== 最终层次聚类 ====================
cat("\n========== 层次聚类 ==========\n")

# 层次聚类函数
hierarchical_clustering <- function(S, n_clusters, method = "average") {
  # S: 共识矩阵（相似度矩阵）
  # n_clusters: 目标簇数
  # method: 链接方法
  
  # 转换为距离矩阵（相似度越高，距离越小）
  dist_matrix <- as.dist(1 - S)
  
  # 层次聚类
  hc <- hclust(dist_matrix, method = method)
  
  # 切割树得到指定数量的簇
  pred_labels <- cutree(hc, k = n_clusters)
  
  return(list(
    labels = pred_labels,
    hclust_obj = hc,
    method = method
  ))
}

# 执行层次聚类
hc_result <- hierarchical_clustering(
  S = S_consensus,
  n_clusters = n_clusters,
  method = "average"  # 可以使用"ward.D", "complete", "average"等
)

pred_labels <- hc_result$labels
names(pred_labels) <- colnames(X)  # 给预测标签命名

# ==================== 评估聚类结果 ====================
cat("\n========== 聚类评估 ==========\n")

# 确保标签长度一致
stopifnot(length(true_labels) == length(pred_labels))
stopifnot(length(true_labels) == ncol(X))

cat(sprintf("评估斑点数: %d\n", length(true_labels)))

# 计算评估指标
ari <- ARI(true_labels, pred_labels)   
nmi <- NMI(true_labels, pred_labels)  

cat(sprintf("聚类性能: ARI=%.4f, NMI=%.4f\n", ari, nmi))



