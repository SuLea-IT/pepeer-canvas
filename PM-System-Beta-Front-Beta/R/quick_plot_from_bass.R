###############################################################################
# quick_plot_from_bass.R —— 从已保存的 BASS 对象快速重新生成聚类图
#
# 使用方法：
# 1. 确保已完整运行过 `run_BASS_local_fixed.R`，并生成了 `bass_object.rds`。
# 2. 根据需要，调整此脚本中的绘图参数（如 pt.size.factor, palette_vec 等）。
# 3. 运行此脚本，即可在数秒内生成新的 PDF 图。
###############################################################################

## ---------- 可调参数 ---------- ##
workdir      <- "D:/xjh/jsq"         # *.rds 所在目录
rds_file     <- "Merged_object.rds"  # 原始 Seurat 对象
bass_rds     <- "bass_object.rds"    # 已保存的 BASS 结果
level        <- "level4_rerun"       # 输出文件名，可修改以防覆盖
pt_size      <- 3.5                  # 点的大小
## --------------------------------##

setwd(workdir)

## 1. 加载基础包 ------------------------------------------------------------
need_pkgs <- c("Seurat", "SeuratObject", "ggplot2", "dplyr")
for (p in need_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

## 2. 加载 BASS 包或源码 (仅为加载对象所需) --------------------------------
if (!requireNamespace("BASS", quietly = TRUE)) {
  if (file.exists("createBASSObject.R")) {
    message("⚠️  BASS 包未安装，使用本地源码 …")
    src_files <- list.files(pattern = "\\.R$", full.names = TRUE)
    invisible(lapply(src_files, source))
  } else {
    stop("❌ 无法加载 BASS。请安装包或确认源码完整，以便读取 BASS 对象。")
  }
} else {
  library(BASS)
}

## 3. 加载数据对象 ---------------------------------------------------------
message("ℹ️  正在加载 Seurat 对象: ", rds_file)
seu <- readRDS(rds_file)

message("ℹ️  正在加载已计算的 BASS 对象: ", bass_rds)
if (!file.exists(bass_rds)) {
  stop("❌ 未找到 ", bass_rds, "！请先完整运行 `run_BASS_local_fixed.R`。")
}
bass <- readRDS(bass_rds)

## 4. 写回 Seurat 元数据 (同主脚本第 7 步) ---------------------------------
final_cell_names <- unlist(bass@cell_name)
cluster_vec <- unlist(bass@results$c, use.names = FALSE)

if (length(final_cell_names) != length(cluster_vec)) {
  stop(sprintf(
    "❌ BASS post-run cell count (%d) does not match cluster result count (%d).",
    length(final_cell_names),
    length(cluster_vec)
  ))
}

names(cluster_vec) <- final_cell_names
seu$bass_clusters <- NA
seu$bass_clusters[names(cluster_vec)] <- factor(cluster_vec)

## 5. 绘图 (同主脚本第 8 步) -------------------------------------------------
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

# 从可调参数获取 nct 和 nsd，以保持文件名一致性
nct <- bass@C
nsd <- bass@R
outfile <- sprintf("c1.nct%s.nsd%s.%s.pdf", nct, nsd, level)

pdf(outfile, 7, 7)
p <- SpatialDimPlot(
  seu, group.by = "bass_clusters",
  label = FALSE, pt.size.factor = pt_size, # 使用可调参数
  cols  = cluster_cols,
  image.alpha = 0,
  crop = FALSE
)
print(p)
dev.off()

message("✅ 快速绘图完成！PDF 已保存至：", outfile)
