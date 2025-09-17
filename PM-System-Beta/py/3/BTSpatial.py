import sys
import scanpy as sc
import pandas as pd
import numpy as np
from pandas.api.types import CategoricalDtype
import os
import matplotlib.pyplot as plt
import squidpy as sq
save_path = sys.argv[2]
mat = sys.argv[1]

# 定义一个保存图像的函数
def save_figure(adata, genes, save_path, filename, plot_type, dpi=300, Finally=False):
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    fig = None

    if plot_type == 'umap':
        fig = sc.pl.umap(
            adata,
            color=["my_score", "leiden"],
            frameon=False,
            ncols=3,
            size=5,
            show=False,
            return_fig=True
        )
    elif plot_type == 'dotplot':
        sc.pl.dotplot(adata, "my_score", groupby="leiden", standard_scale="var",show=False)
        fig = plt.gcf()  # Get the current figure
    elif plot_type == 'spatial':
        sq.pl.spatial_scatter(adata, color=genes, size=1, shape=None, edges_color="black")
        fig = plt.gcf()  # Get the current figure
    else:
        return

    if fig is not None:
        # 保存图像并关闭
        full_save_path = os.path.join(save_path, f'{filename}_{plot_type}.pdf')
        fig.savefig(full_save_path, format='pdf', dpi=dpi)
        plt.close(fig)

    if Finally:
        print(save_path)




# 读取数据并进行预处理
gene_id_file_path = f'{mat}/geneid.txt'
# 读取百迈克空间转录组数据
coord = f"{mat}/barcodes_pos.tsv.gz"
adata = sc.read_10x_mtx(path=mat)
pos = pd.read_csv(coord, sep="\t", names=["cellID", "row", "col"])

# 处理坐标数据
pos2 = pos[pos["cellID"].isin(adata.to_df().index.to_list())]
cat_size_order = CategoricalDtype(adata.to_df().index.to_list(), ordered=False)
pos2.loc[:, 'cellID'] = pos2['cellID'].astype(cat_size_order)
pos2.set_index('cellID', inplace=True)

# 移除缺失细胞的数据
cells_in_adata = set(adata.obs_names)
cells_in_pos2 = set(pos2.index)
missing_cells = cells_in_adata - cells_in_pos2
adata = adata[~adata.obs_names.isin(missing_cells)]

# 将坐标数据添加到AnnData对象
adata.obsm['spatial'] = np.array(pos2[['row', 'col']])

# 过滤数据
min_genes = 100  # 每个细胞的最小基因数
min_cells = 3    # 每个基因的最小细胞数
sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)

# 数据预处理步骤
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
sc.pp.pca(adata)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata)
sc.tl.tsne(adata)
sc.tl.leiden(adata, resolution=0.5)

# 检查基因ID文件并生成相应的图像
if os.path.exists(gene_id_file_path):
    gene_ids = []
    with open(gene_id_file_path, 'r') as file:
        next(file)  # 跳过标题行
        for line in file:
            gene_ids.append(line.strip())

    # 为空间数据生成多个基因的UMAP、dotplot和空间散点图
sc.tl.score_genes(adata, gene_ids, ctrl_size=50, n_bins=25, score_name='my_score')
save_figure(adata, gene_ids, save_path, f'spatial_bin_gene_projection_gene', plot_type='umap')
save_figure(adata, gene_ids, save_path, f'spatial_bin_dotplot_gene_projection_gene', plot_type='dotplot')
save_figure(adata, gene_ids, save_path, f'spatial_bin_spatial_projection_gene', plot_type='spatial',Finally=True)
