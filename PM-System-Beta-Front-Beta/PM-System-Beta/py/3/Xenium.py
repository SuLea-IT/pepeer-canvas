import scanpy as sc
import sys
import os
import matplotlib.pyplot as plt
import pandas as pd
import squidpy as sq
save_path = sys.argv[2]
mat = sys.argv[1]

# 定义一个保存图像的函数
def save_figure(adata, genes, save_path, filename, plot_type, dpi=300, Finally=False):
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    if plot_type == 'umap':
        fig_list = sc.pl.umap(
            adata,
            color=["my_score", "leiden"],
            frameon=False,
            ncols=3,
            size=5,
            show=False,
            return_fig=True  # Make sure it returns the figure
        )
        # Check if it's a list and loop through it to save each figure
        if isinstance(fig_list, list):
            for idx, fig in enumerate(fig_list):
                full_save_path = os.path.join(save_path, f'{filename}_{plot_type}_{idx}.pdf')
                fig.savefig(full_save_path, format='pdf', dpi=dpi)
                plt.close(fig)
        else:
            full_save_path = os.path.join(save_path, f'{filename}_{plot_type}.pdf')
            fig_list.savefig(full_save_path, format='pdf', dpi=dpi)
            plt.close(fig_list)

    elif plot_type == 'dotplot':
        # Dotplot only returns a single figure
        sc.pl.dotplot(adata, "my_score", groupby="leiden", standard_scale="var", show=False)
        fig = plt.gcf()  # Get the current figure
        full_save_path = os.path.join(save_path, f'{filename}_{plot_type}.pdf')
        fig.savefig(full_save_path, format='pdf', dpi=dpi)
        plt.close(fig)
    elif plot_type == 'spatial':
        sq.pl.spatial_scatter(adata, color=genes, size=1, shape=None, edges_color="black")
        fig =plt.gcf()
        full_save_path = os.path.join(save_path, f'{filename}_{plot_type}.pdf')
        fig.savefig(full_save_path, format='pdf', dpi=dpi)
        plt.close(fig)

    if Finally:
        print(save_path)


# 读取数据并进行预处理
gene_id_file_path = f'{mat}/geneid.txt'
xenium_test = sc.read_10x_h5(filename=f"{mat}/cell_feature_matrix.h5")
df = pd.read_csv(f"{mat}/cells.csv.gz")
df.set_index(xenium_test.obs_names, inplace=True)
xenium_test.obs = df.copy()
xenium_test.obsm["spatial"] = xenium_test.obs[["x_centroid", "y_centroid"]].copy().to_numpy()

# 创建计数层
xenium_test.layers["counts"] = xenium_test.X.copy()

# 数据标准化和对数转换
sc.pp.normalize_total(xenium_test, inplace=True)
sc.pp.log1p(xenium_test)

# 筛选高度变异基因
sc.pp.highly_variable_genes(xenium_test, flavor="seurat", n_top_genes=2000)

# PCA降维
sc.pp.pca(xenium_test)

# 计算邻居图
sc.pp.neighbors(xenium_test, n_neighbors=10, n_pcs=30)

# 运行UMAP嵌入
sc.tl.umap(xenium_test)
xenium_test.X = xenium_test.X.toarray()
# 运行t-SNE嵌入
sc.tl.tsne(xenium_test, use_rep='X_pca')


# 运行Leiden聚类
sc.tl.leiden(xenium_test, resolution=0.5)

gene_ids = []
# 检查基因ID文件并生成相应的图像
if os.path.exists(gene_id_file_path):
    with open(gene_id_file_path, 'r') as file:
        next(file)  # 跳过标题行
        for line in file:
            gene_ids.append(line.strip())
sc.tl.score_genes(xenium_test, gene_ids, ctrl_size=50, n_bins=25, score_name='my_score')

save_figure(xenium_test, gene_ids, save_path, f'spatial_gene_projection', plot_type='umap')
save_figure(xenium_test, gene_ids, save_path, f'spatial_gene_projection', plot_type='dotplot')
save_figure(xenium_test, gene_ids, save_path, f'spatial_gene_projection', plot_type='spatial', Finally=True)
