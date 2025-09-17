import sys
import pandas as pd
import numpy as np
import gzip
from pandas.api.types import CategoricalDtype
import scanpy as sc
import os
import matplotlib.pyplot as plt
import squidpy as sq
import subprocess
from PIL import Image

Image.MAX_IMAGE_PIXELS = None
save_path = sys.argv[2]


# 定义保存图片的函数
def save_figure(adata, save_path, filename, plot_type='umap', key=None, dpi=300, Finally=False):
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    fig = None
    if plot_type == 'umap':
        fig = sc.pl.umap(adata, color=key, show=False, return_fig=True)
    elif plot_type == 'tsne':
        sc.pl.tsne(adata, color=key, show=False, return_fig=True)
    elif plot_type == 'spatial':
        sq.pl.spatial_scatter(adata, size=1, shape=None, edges_color="black", color=key)
        fig = plt.gcf()
    elif plot_type == 'cluster':
        cluster_labels = adata.obs['leiden']
        csv_file_path = os.path.normpath(os.path.join(save_path, "cluster.csv"))
        save_dir = os.path.normpath(os.path.join(save_path, "result"))

        # Save cluster labels to CSV
        pd.DataFrame(cluster_labels).to_csv(csv_file_path, sep=",", index=True)
        command = [
            "python", "py/util/cluster_plot_black_lineR.py",
            "--save_dir", save_dir,
            "--cells_path", f"{mat}/S_37.npy",
            "--cluster_path", csv_file_path,
            "--redo"
        ]
        subprocess.run(command, capture_output=True, text=True, check=True)
        tif_file_path = os.path.join(save_dir, "cell_cluster_with_legend_img.tif")
        pdf_file_path = os.path.join(save_dir, "cell_cluster_with_legend_img.jpg")
        if os.path.exists(tif_file_path):
            with Image.open(tif_file_path) as img:
                img.save(pdf_file_path, "JPEG", resolution=300)
        for filename in os.listdir(save_dir):
            file_path = os.path.join(save_dir, filename)
            if filename != "cell_cluster_with_legend_img.jpg" and os.path.isfile(file_path):
                os.remove(file_path)
    else:
        return
    if fig is not None:
        # 保存图像并关闭
        full_save_path = os.path.join(save_path, f'{filename}_{plot_type}.pdf')
        fig.savefig(full_save_path, format='pdf', dpi=dpi)
        plt.close(fig)
    if Finally:
        print(save_path)


# 从命令行参数获取路径
mat = sys.argv[1]
result_dir = os.path.join(save_path, "result")
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
coord = f"{mat}/barcodes_pos.tsv.gz"
st_test = sc.read_10x_mtx(path=mat)
pos = pd.read_csv(coord, sep="\t", names=["cellID", "row", "col"])
pos2 = pos[pos["cellID"].isin(st_test.to_df().index.to_list())]
cat_size_order = CategoricalDtype(st_test.to_df().index.to_list(), ordered=False)
pos2['cellID'] = pos2['cellID'].astype(cat_size_order)
pos2.set_index('cellID', inplace=True)
cells_in_adata = set(st_test.obs_names)
cells_in_pos2 = set(pos2.index)
missing_cells = cells_in_adata - cells_in_pos2
st_test = st_test[~st_test.obs_names.isin(missing_cells)]
st_test.obsm['spatial'] = np.array(pos2)
# 过滤数据
min_genes = 100  # 每个细胞的最小基因数
min_cells = 3  # 每个基因的最小细胞数
sc.pp.filter_cells(st_test, min_genes=min_genes)
sc.pp.filter_genes(st_test, min_cells=min_cells)
st_test.layers["counts"] = st_test.X.copy()
sc.pp.normalize_total(st_test, inplace=True)
sc.pp.log1p(st_test)
sc.pp.highly_variable_genes(st_test, flavor="seurat", n_top_genes=2000)
sc.pp.pca(st_test)
sc.pp.neighbors(st_test, n_neighbors=10, n_pcs=30)
sc.tl.umap(st_test)
sc.tl.tsne(st_test)
sc.tl.leiden(st_test, resolution=0.5)
# 输出聚类图（五张）
save_figure(st_test, save_path, 'clustered_data_umap', plot_type='umap')  # 无leiden
save_figure(st_test, save_path, 'clustered_data_tsne', plot_type='tsne')  # 无leiden
save_figure(st_test, save_path, 'clustered_data_umap_leiden', plot_type='umap', key='leiden')
save_figure(st_test, save_path, 'clustered_data_tsne_leiden', plot_type='tsne', key='leiden')
save_figure(st_test, save_path, 'clustered_data_spatial', plot_type='spatial', key='leiden')
save_figure(st_test, save_path, 'clustered_data', plot_type='cluster', key='leiden', Finally=True)
