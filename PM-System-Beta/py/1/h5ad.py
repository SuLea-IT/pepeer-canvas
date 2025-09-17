# 导入必要的库
# 1/h5ad.py
import sys
import pandas as pd
import numpy as np
import gzip
from pandas.api.types import CategoricalDtype
import scanpy as sc
import os
import matplotlib.pyplot as plt
import squidpy as sq
save_path = sys.argv[2]
def save_figure(adata, save_path, filename, plot_type='umap', key=None, dpi=300,Finally=False):
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    if plot_type == 'umap':
        sc.pl.umap(adata, color=key,show=False)
    elif plot_type == 'tsne':
        sc.pl.tsne(adata, color=key,show=False)
    elif plot_type == 'spatial':
        sq.pl.spatial_scatter(adata, size=1, shape=None, edges_color="black", color=key)
    else:
        return

    fig = plt.gcf()  # 获取当前图形
    full_save_path = os.path.join(save_path, f'{filename}_{plot_type}.pdf')
    fig.savefig(full_save_path, format='pdf', dpi=dpi)
    plt.close(fig)
    if Finally:
        print(save_path)
mat = sys.argv[1]
#anndata格式
adata = sc.read_h5ad(f"{mat}/sc_all_adata.h5ad")
#过滤数据
min_genes = 100  # 每个细胞的最小基因数
min_cells = 3 # 每个基因的最小细胞数
sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)
has_umap = 'X_umap' in adata.obsm
has_tsne = 'X_tsne' in adata.obsm
if "leiden" not in adata.obs:
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


elif has_umap and not has_tsne:
    sc.tl.tsne(adata)
elif not has_umap and has_tsne:
    sc.tl.umap(adata)
save_figure(adata, save_path, 'clustered_data_umap', plot_type='umap')  # 无leiden
save_figure(adata, save_path, 'clustered_data_tsne', plot_type='tsne')  # 无leiden
save_figure(adata, save_path, 'clustered_data_umap_leiden', plot_type='umap', key='leiden')
save_figure(adata, save_path, 'clustered_data_tsne_leiden', plot_type='tsne', key='leiden',Finally=True)
#
# save_figure(adata, save_path, 'clustered_data_spatial', plot_type='spatial', key='leiden',Finally=True)
