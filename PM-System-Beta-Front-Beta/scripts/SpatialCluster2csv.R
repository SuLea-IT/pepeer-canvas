# 加载所需的库
library(Seurat)
library(dplyr)
library(readr)
library(tibble)

# --- 用户需修改的路径 ---

# 1. RDS 文件路径
# 这个文件包含 Seurat 对象，其中有细胞的 cluster 分组信息。
rds_file_path <- "path/to/your/seurat_object.rds"

# 2. 细胞条码与空间坐标文件路径
# 这是包含细胞条码及其空间坐标的文件。
# 用户提供的示例: "D:/xjh/leaf_rds/leaf_barcodes_pos.tsv.gz"
barcode_pos_path <- "path/to/your/leaf_barcodes_pos.tsv.gz"

# 3. 输出 CSV 文件路径
# 这是最终生成的包含合并后信息的 CSV 文件的保存路径。
csv_output_path <- "path/to/your/spatial_clusters_output.csv"

# --- 脚本执行逻辑 ---

# 检查输入文件是否存在
if (!file.exists(rds_file_path)) {
  stop("错误: RDS 文件不存在: ", rds_file_path)
}
if (!file.exists(barcode_pos_path)) {
  stop("错误: 细胞条码坐标文件不存在: ", barcode_pos_path)
}

# 读取 Seurat 对象
seurat_object <- readRDS(rds_file_path)

# 读取空间坐标文件
# 假设文件是 .tsv.gz 格式，并且没有列名。
# 根据 10x Visium 的标准格式，列通常是: barcode, in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres
barcode_positions <- read_tsv(barcode_pos_path, col_names = FALSE)
# 为列指定更清晰的名称
colnames(barcode_positions) <- c("barcode", "in_tissue", "array_row", "array_col", "y", "x")


# 提取 cluster 信息
# cluster 信息通常存储在 Seurat 对象的 meta.data 中
# 'seurat_clusters' 是常见的列名，如果您的列名不同，请在此处修改
if (!"seurat_clusters" %in% colnames(seurat_object@meta.data)) {
    warning("警告: 在 meta.data 中未找到 'seurat_clusters' 列。请检查您的 Seurat 对象并相应地修改脚本。")
}

cluster_data <- seurat_object@meta.data %>%
  as.data.frame() %>%
  select(any_of("seurat_clusters")) %>% # 使用 any_of 避免列不存在时出错
  rownames_to_column("barcode")

# 将空间坐标数据与 cluster 数据通过 'barcode' 列进行合并
# inner_join 会保留两个数据集中都存在的条码
merged_data <- inner_join(barcode_positions, cluster_data, by = "barcode")

# 将合并后的数据写入 CSV 文件
write.csv(merged_data, file = csv_output_path, row.names = FALSE)

# 打印成功信息
cat("成功将空间坐标和 cluster 数据合并并导出到: ", csv_output_path, "\n")
