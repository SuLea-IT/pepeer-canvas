# 1. 加载必要的库
cat("正在加载库...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
})

# 2. 定义输入和输出路径
# ---------------------------------
# 包含 Seurat 对象的 RDS 文件路径。
rds_file_path <- "D:/path/to/your/seurat_object.rds" # <--- 请将此路径更改为您的 RDS 文件路径

# CSV 文件将被保存的目录。
output_base_dir <- "D:/path/to/your/output_directory" # <--- 请将此路径更改为您希望保存CSV的目录

# =================================================================
# 主要处理函数
# =================================================================
process_seurat_clusters <- function(rds_path, output_dir) {
  
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
  
  # 4. 导出聚类信息
  cat("正在导出聚类信息...\n")
  metadata <- seurat_object@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("barcode")
  
  # 自动检测聚类列
  cluster_col_name <- "seurat_clusters"
  if (!cluster_col_name %in% names(metadata)) {
    potential_cols <- names(metadata)[grepl("cluster", names(metadata), ignore.case = TRUE)]
    if (length(potential_cols) > 0) {
      cluster_col_name <- potential_cols[1]
    } else if ("ident" %in% names(metadata)) {
      cluster_col_name <- "ident"
    } else {
      stop("错误: 无法在元数据中自动找到聚类列 (如 'seurat_clusters').")
    }
  }
  
  cat(paste("使用 '", cluster_col_name, "' 列作为聚类 ID。\n", sep=""))
  
  # 5. 提取坐标信息 (例如 UMAP)
  cat("正在提取 UMAP 坐标...\n")
  
  # 检查是否存在 UMAP reduction
  if (!"umap" %in% names(seurat_object@reductions)) {
    stop("错误: 在 Seurat 对象中未找到 'umap' reduction。请确保已运行 UMAP。")
  }
  
  # 提取坐标并转换为数据框
  coords <- as.data.frame(seurat_object@reductions$umap@cell.embeddings)
  colnames(coords) <- c("x", "y") # 通常是 UMAP_1 和 UMAP_2
  coords <- coords %>%
    tibble::rownames_to_column("barcode")
  
  # 6. 合并聚类和坐标数据
  cat("正在合并聚类和坐标数据...\n")
  
  # 准备聚类数据
  cluster_data <- metadata %>%
    select(barcode, cluster = !!sym(cluster_col_name))
  
  # 将聚类数据与坐标数据合并
  combined_data <- dplyr::inner_join(cluster_data, coords, by = "barcode")
  
  # 定义输出CSV文件路径
  cluster_csv_path <- file.path(output_dir, "clusters_with_coords.csv")
  
  # 写入CSV文件
  write_csv(combined_data, cluster_csv_path)
  cat(paste("包含坐标的聚类数据已成功保存至:", cluster_csv_path, "\n"))
  
  cat("处理完成。\n")
}

# =================================================================
# 脚本执行
# =================================================================
# 基于RDS文件名创建特定的输出子目录
output_dir_name <- tools::file_path_sans_ext(basename(rds_file_path))
final_output_dir <- file.path(output_base_dir, output_dir_name)

# 使用 tryCatch 块来捕获和报告错误
tryCatch({
  process_seurat_clusters(rds_file_path, final_output_dir)
}, error = function(e) {
  cat(paste("发生错误:", e$message, "\n"))
})