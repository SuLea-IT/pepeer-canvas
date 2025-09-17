###############################################################################
# run_BASS_local.R —— 读取单个 Seurat 对象 → BASS 建模 → 写回 meta.data
#                   → 绘制空间聚类图（支持 Seurat v3–v5 & 多切片）
###############################################################################

## ---------- 可调参数 ---------- ##
workdir  <- "D:/xjh/jsq"         # *.rds 所在目录
rds_file <- "Merged_object2.rds"  # Seurat 对象
level    <- "level4"             # 仅用于输出文件名
nct      <- 25                   # 期望细胞类型数 (cluster 上限)
nsd      <- 5                    # 空间域数
n_cores  <- max(1, parallel::detectCores() - 2) # 并行核数 (自动检测并保留2个核心)
## --------------------------------##

setwd(workdir)

## 1. 加载基础包 ------------------------------------------------------------
need_pkgs <- c("Seurat", "SeuratObject", "ggplot2", "dplyr", "parallel")
for (p in need_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}
options(mc.cores = n_cores) # 为 R 会话设置全局并行核心数

## 2. 加载 / 回退 BASS -------------------------------------------------------
bass_ok <- FALSE
if (requireNamespace("BASS", quietly = TRUE)) {
  library(BASS); bass_ok <- TRUE
} else if (file.exists("createBASSObject.R")) {
  message("⚠️  BASS 包未安装，使用本地源码 …")
  src_files <- list.files(pattern = "\\.R$", full.names = TRUE)
  invisible(lapply(src_files, source))
  bass_ok <- exists("createBASSObject") &&
    exists("BASS.preprocess")  &&
    exists("BASS.run")         &&
    exists("BASS.postprocess")
}
if (!bass_ok) stop("❌ 无法加载 BASS。请安装包或确认源码完整。")

## 3. 读取 Seurat 对象 -------------------------------------------------------
seu <- readRDS(rds_file)
assay_name <- DefaultAssay(seu)
assay_obj  <- seu[[assay_name]]

## 4. 提取表达矩阵 —— 合并全部 counts layer (兼容 v3~v5) --------------------

assay_name  <- DefaultAssay(seu)
layer_names <- tryCatch(Layers(seu[[assay_name]]), error = function(e) character(0))

if (length(layer_names) > 0) {                      # ---- v5: 多 layer ----
  counts_layers <- grep("^counts", layer_names, value = TRUE)
  if (length(counts_layers) == 0)
    stop("❌ 没有任何 counts layer。当前 layer：",
         paste(layer_names, collapse = ", "))
  
  message("ℹ️  检测到 ", length(counts_layers),
          " 个 counts layer，将全部合并：\n  ",
          paste(counts_layers, collapse = ", "))
  
  ## 1) 逐层读取
  expr_list <- lapply(counts_layers, function(lay) {
    GetAssayData(seu, assay = assay_name, layer = lay)
  })
  
  ## 2) 找所有 layer 共有的基因
  common_genes <- Reduce(intersect, lapply(expr_list, rownames))
  message("ℹ️  共同基因数：", length(common_genes))
  
  if (length(common_genes) < 200)
    warning("⚠️  共同基因 < 200，可能过少，建议检查 layer 基因一致性。")
  
  ## 3) 对齐行顺序并裁剪
  expr_list <- lapply(expr_list, function(mat) mat[common_genes, , drop = FALSE])
  
  ## 4) 合并并去重列
  expr_all <- do.call(cbind, expr_list)
  expr_all <- expr_all[, !duplicated(colnames(expr_all))]
  
} else {                                            # ---- v3/v4 ----
  expr_all <- GetAssayData(seu, assay = assay_name, slot = "counts")
}


## 5. 准备每张切片的 expr & 坐标列表 ----------------------------------------
#   BASS 需要 X = list(expr1, expr2, …), xy = list(xy1, xy2, …)
X_list  <- list()
xy_list <- list()
cell_order <- c()   # 保存所有细胞顺序，稍后给 names()

coord_candidates <- list(
  c("x", "y"),
  c("imagecol", "imagerow"),
  c("pxl_col_in_fullres", "pxl_row_in_fullres"),
  c("col", "row")
)

message("\nℹ️  总共检测到 ", length(names(seu@images)), " 个图像 (切片/样本)，开始逐一处理...")

for (img in names(seu@images)) {
  message("\n--- 开始处理切片: '", img, "' ---")
  coords <- seu@images[[img]]@coordinates
  
  message("  - 在 '", img, "' 的坐标中找到 ", nrow(coords), " 个细胞。")
  
  # 找 (x,y)
  found <- NULL
  for (pair in coord_candidates) {
    if (all(pair %in% colnames(coords))) { found <- pair; break }
  }
  if (is.null(found))
    stop("❌ 图像 ", img, " 的 coordinates 不含 (x,y) 列，列名：",
         paste(colnames(coords), collapse = ", "))
  
  cells_img <- rownames(coords)
  cells_in_expr <- colnames(expr_all)
  
  message("  - 在合并后的总表达矩阵中找到 ", length(cells_in_expr), " 个细胞。")
  
  cells_img_intersect <- intersect(cells_img, cells_in_expr)
  
  message("  - '", img, "' 的坐标与总表达矩阵的【共同细胞数】: ", length(cells_img_intersect))
  
  if (length(cells_img_intersect) == 0) {
    message("  ⚠️  [跳过]：由于共同细胞为 0，此切片无法处理。")
    next
  }
  
  message("  ✅ [成功]：此切片将被纳入 BASS 分析。")
  
  expr_mat <- expr_all[, cells_img_intersect, drop = FALSE]
  xy_mat   <- as.matrix(coords[cells_img_intersect, found])
  colnames(xy_mat) <- c("x", "y")
  
  X_list[[img]]  <- expr_mat
  xy_list[[img]] <- xy_mat
  cell_order <- c(cell_order, cells_img_intersect)
}

if (length(X_list) == 0)
  stop("❌ 没找到任何有效切片 / 细胞可用于 BASS。")

## 6. 运行 BASS -------------------------------------------------------------
bass <- createBASSObject(
  X  = X_list,
  xy = xy_list,
  C  = nct,
  R  = nsd,
  beta_method = "SW"
) |>
  BASS.preprocess() |>
  BASS.run() |>
  BASS.postprocess()

## 7. 写回 Seurat 元数据 -----------------------------------------------------
# BASS@results$c 是一个包含各切片聚类结果的列表。
# unlist 会将其合并成一个向量，向量的 names 就是细胞 ID，值是聚类编号。
# 这样可以确保细胞名和聚类结果的长度和顺序完全一致。
cluster_vec <- unlist(bass@results$c)

# 将聚类结果安全地添加回 Seurat 对象的 meta.data
seu$bass_clusters <- NA  # 初始化列，所有细胞设为 NA
# 仅为有聚类结果的细胞赋值
seu$bass_clusters[names(cluster_vec)] <- factor(cluster_vec)

## 8. 绘图 -------------------------------------------------------------------
palette_vec <- c(
  "#F56867","#FEB915","#59BE86","#7495D3","#6D1A9C","#15821E","#3A84E6",
  "#DB4C6C","#99db27","#e07233","#ce2523","#f7aa5d","#9E7A7A","#AF5F3C",
  "#F9BD3F","#DAB370","#268785","#ed1299","#09f9f5","#246b93","#cc8e12",
  "#d561dd","#c93f00","#ddd53e","#4aef7b","#e86502","#9ed84e","#39ba30",
  "#6ad157","#cebb10","#03827f","#931635","#373bbf","#a1ce4c","#ef3bb6",
  "#d66551","#1a918f","#ff66fc","#2927c4","#7149af","#57e559","#8e3af4",
  "#f9a270","#22547f","#db5e92","#edd05e","#6f25e8","#0dbc21","#280f7a",
  "#6373ed","#5b910f","#7b34c1","#0cf29a","#997273","#d80fc1","#877F6C",
  "#93796C","#C798EE","#554236","#8249aa"
)

uniq_clu <- sort(na.omit(unique(seu$bass_clusters)))
if (length(uniq_clu) > length(palette_vec))
  stop("❌ Cluster 数超过调色板长度，请自行扩充 palette_vec。")

cluster_cols <- setNames(palette_vec[seq_along(uniq_clu)], uniq_clu)

outfile <- sprintf("c1.nct%s.nsd%s.%s.pdf", nct, nsd, level)
pdf(outfile, 7, 7)
p <- SpatialDimPlot(
  seu, group.by = "bass_clusters",
  label = FALSE, pt.size.factor = 3.5,
  cols  = cluster_cols,
  image.alpha = 0,
  crop = FALSE
)
print(p)
dev.off()

message("✅ BASS 聚类完成：", length(uniq_clu),
        " 个 cluster；PDF 已保存至：", outfile)
