import scanpy as sc
import sys
import os
import matplotlib.pyplot as plt
import pandas as pd
import squidpy as sq
# 读取输入参数
save_path = sys.argv[2]
mat = sys.argv[1]

# 定义保存图像的函数
def save_figure(adata, genes, save_path, filename, plot_type, dpi=300, Finally=False):
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    fig = None
    if plot_type == 'umap':
        fig = sc.pl.umap(
            adata,
            color=[genes, "leiden"],
            frameon=False,
            ncols=3,
            size=5, show=False, return_fig=True
        )
    elif plot_type == 'dotplot':
        # 绘制 dotplot 图形
        sc.pl.dotplot(adata, [genes] if isinstance(genes, str) else genes, groupby="leiden", show=False)
        fig = plt.gcf()  # 获取当前的 Matplotlib 图像
    elif plot_type == 'spatial':
        sq.pl.spatial_scatter(adata, color=genes, size=1, shape=None, edges_color="black", return_fig=True)
    else:
        return

    if fig is not None:
        # 保存图像并关闭
        full_save_path = os.path.join(save_path, f'{filename}_{plot_type}.pdf')
        fig.savefig(full_save_path, format='pdf', dpi=dpi)
        plt.close(fig)

    if Finally:
        print(save_path)



# 读取AnnData对象
sc_test = sc.read_h5ad(f"{mat}/sc_test.h5ad")

# 检查是否已经存在 "leiden" 列，如果不存在，执行聚类和降维
if "leiden" not in sc_test.obs:
    # 创建计数层
    sc_test.layers["counts"] = sc_test.X.copy()

    # 数据标准化和对数转换
    sc.pp.normalize_total(sc_test, inplace=True)
    sc.pp.log1p(sc_test)

    # 筛选高度变异基因
    sc.pp.highly_variable_genes(sc_test, flavor="seurat", n_top_genes=2000)

    # PCA降维
    sc.pp.pca(sc_test)

    # 计算邻居图
    sc.pp.neighbors(sc_test, n_neighbors=10, n_pcs=30)

    # 运行UMAP嵌入
    sc.tl.umap(sc_test)

    # 运行t-SNE嵌入，使用随机初始化避免PCA和稀疏矩阵冲突
    sc.tl.tsne(sc_test, use_rep='X_pca')

    # 运行Leiden聚类
    sc.tl.leiden(sc_test, resolution=0.5)

# 检查文件并生成相应图像
gene_id_file_path = f'{mat}/geneid.txt'

if os.path.exists(gene_id_file_path):
    # 如果基因ID文件存在，读取基因并生成多个基因的UMAP图
    gene_ids = []
    with open(gene_id_file_path, 'r') as file:
        next(file)  # 跳过标题行
        for line in file:
            gene_ids.append(line.strip())

    # 生成多个基因的UMAP图
    for idx, gene_id in enumerate(gene_ids):
        is_last = idx == len(gene_ids) - 1  # 检查是否是最后一个基因ID
        save_figure(sc_test, gene_id, save_path, f'gene_projection_{gene_id}', plot_type='umap')
        save_figure(sc_test, gene_id, save_path, f'dotplot_gene_projection_{gene_id}', plot_type='dotplot')
        save_figure(sc_test,gene_id, save_path, f'spatial_dotplot_gene_projection_{gene_id}',
                    plot_type='spatial',Finally=is_last)
