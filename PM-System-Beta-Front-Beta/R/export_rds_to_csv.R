# 1. 加载必要的库
cat("正在加载库...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
})

# 2. 定义输入和输出路径
# ---------------------------------
# 包含 Seurat 对象的 RDS 文件。
rds_file_path <- "D:/xjh/laj/PM-System-Beta-Front-Beta/R/10DPA_L7.rds"

# CSV 文件将被保存的目录。
output_base_dir <- "D:/xjh/laj/PM-System-Beta-Front-Beta/R/csv_exports"

# =================================================================
# 主要处理函数
# =================================================================
process_seurat_object <- function(rds_path, output_dir) {
  
  # 检查 RDS 文件是否存在
  if (!file.exists(rds_path)) {
    stop(paste("错误: 未在以下路径找到 RDS 文件", rds_path))
  }
  
  # 如果输出目录不存在，则创建它
  if (!dir.exists(output_dir)) {
    cat(paste("正在创建输出目录:", output_dir, "\n"))
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 3. 读取 Seurat 对象
  cat(paste("正在从以下路径读取 Seurat 对象:", rds_path, "\n"))
  seurat_object <- readRDS(rds_path)
  cat("Seurat 对象加载成功。\n")
  
  # 4. 提取坐标
  cat("正在提取坐标...\n")
  coords <- NULL
  if (length(seurat_object@images) > 0) {
    cat("找到图像数据，尝试提取空间坐标。\n")
    coords <- seurat_object@images[[1]]@coordinates
    if ("tissue.x" %in% colnames(coords) && "tissue.y" %in% colnames(coords)) {
      coords <- coords %>% rename(x = tissue.x, y = tissue.y)
    }
  } else if (length(seurat_object@reductions) > 0) {
    cat("未找到空间图像数据。尝试从降维结果中提取坐标 (UMAP/tSNE)。\n")
    reduction_name <- if ("umap" %in% names(seurat_object@reductions)) "umap" else if ("tsne" %in% names(seurat_object@reductions)) "tsne" else NULL
    if (!is.null(reduction_name)) {
      cat(paste("使用 '", reduction_name, "' 降维结果作为坐标。\n", sep=""))
      coords <- as.data.frame(seurat_object@reductions[[reduction_name]]@cell.embeddings)
      colnames(coords) <- c("x", "y")
    }
  }
  if (is.null(coords)) stop("错误: 无法从 Seurat 对象中提取任何坐标 (空间或降维)。")
  coords <- coords %>% tibble::rownames_to_column("barcode")
  
  # 5. 导出聚类信息
  cat("正在导出聚类信息...\n")
  metadata <- seurat_object@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("barcode")
  cluster_col_name <- "seurat_clusters"
  if (!cluster_col_name %in% names(metadata)) {
    potential_cols <- names(metadata)[grepl("cluster", names(metadata), ignore.case = TRUE)]
    cluster_col_name <- if (length(potential_cols) > 0) potential_cols[1] else if ("ident" %in% names(metadata)) "ident" else "seurat_clusters"
    if (cluster_col_name == "seurat_clusters") {
      warning("无法自动找到聚类列。请检查元数据。")
      metadata$seurat_clusters <- 0
    }
  }
  cat(paste("使用 '", cluster_col_name, "' 列作为聚类 ID。\n", sep=""))
  cluster_data <- metadata %>%
    select(barcode, cluster = !!sym(cluster_col_name)) %>%
    left_join(coords, by = "barcode") %>%
    select(barcode, x, y, cluster)
  cluster_csv_path <- file.path(output_dir, "clusters.csv")
  write_csv(cluster_data, cluster_csv_path)
  cat(paste("聚类数据已保存至:", cluster_csv_path, "\n"))
  
  # 6. 高效导出基因表达数据
  cat("正在高效导出基因表达数据...\n")
  expression_data <- GetAssayData(seurat_object, slot = "data")
  
  # 将表达矩阵转换为长格式数据框
  cat("  - 正在重塑表达矩阵...\n")
  long_expression_df <- as.data.frame(as.matrix(expression_data)) %>%
    tibble::rownames_to_column("gene") %>%
    pivot_longer(
      cols = -gene,
      names_to = "barcode",
      values_to = "expression"
    )
  
  # 与坐标进行一次性连接
  cat("  - 正在与坐标数据连接...\n")
  full_data <- long_expression_df %>%
    left_join(coords, by = "barcode") %>%
    select(barcode, x, y, gene, expression)
  
  # 创建基因数据子目录
  genes_output_dir <- file.path(output_dir, "genes")
  if (!dir.exists(genes_output_dir)) {
    dir.create(genes_output_dir)
  }
  
  # 分组并写入文件
  cat("  - 正在分组并将数据写入 CSV 文件...\n")
  full_data %>%
    group_by(gene) %>%
    group_walk(~ write_csv(.x, file.path(genes_output_dir, paste0(.y$gene, ".csv"))))
  
  cat(paste("已在以下目录保存", length(unique(full_data$gene)), "个基因的表达数据:", genes_output_dir, "\n"))
  cat("处理完成。\n")
}

# =================================================================
# 脚本执行
# =================================================================
output_dir_name <- tools::file_path_sans_ext(basename(rds_file_path))
final_output_dir <- file.path(output_base_dir, output_dir_name)

tryCatch({
  process_seurat_object(rds_file_path, final_output_dir)
}, error = function(e) {
  cat(paste("发生错误:", e$message, "\n"))
})
