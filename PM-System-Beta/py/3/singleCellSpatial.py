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
import stat
import platform
import pwd
Image.MAX_IMAGE_PIXELS = None
save_path = sys.argv[2]
mat = sys.argv[1]

# 定义一个保存图像的函数
def save_figure(adata, genes, save_path, filename, plot_type, dpi=300, Finally=False):
    if not os.path.exists(save_path):
        os.makedirs(save_path)

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
    elif plot_type == 'cluster':
        genes_score=adata.obs['my_score'].reset_index()
        genes_score['index'] = genes_score['index'].apply(lambda x: 'cell_' + x.split('_')[-1])
        csv_file_path = os.path.normpath(os.path.join(save_path, f"genes_score.csv"))
        genes_score.to_csv(csv_file_path, sep=",", index=False)
        save_dir = os.path.normpath(os.path.join(save_path, "result"))

        script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'util', 'gene_plot_s.py'))

        cells_npy_path = os.path.normpath(os.path.join(mat, "S_37.npy"))

        command = [
            "python", script_path,
            "--csv", csv_file_path,
            "--outdir", save_dir,
            "--cmap", 'viridis',
            "--cells_npy", cells_npy_path,
            "--background_color", '#000000',
            "--line_color", '#000000',
            "--name", "gene_score"
        ]
        config_file_path = os.path.join(save_path, "run_config.txt")
        with open(config_file_path, "a") as config_file:
            config_file.write(f"Running gene: genes\n")
            config_file.write(f"Command: {' '.join(command)}\n")
            config_file.write(f"Save directory: {save_dir}\n\n")
        subprocess.run(command, capture_output=True, text=True, check=True)
        tif_file_path = os.path.join(save_dir, f"gene_score_with_legend.tif")
        jpg_file_path = os.path.join(save_dir, f"gene_score_with_legend.jpg")
        if os.path.exists(tif_file_path):
            with Image.open(tif_file_path) as img:
                img.save(jpg_file_path, "JPEG", resolution=300)
            os.remove(tif_file_path)
        # 清理结果目录中的其他文件
        for filename in os.listdir(save_dir):
            file_path = os.path.join(save_dir, filename)
            if not filename.lower().endswith('.jpg') and os.path.isfile(file_path):
                os.remove(file_path)

        fig = None  # No figure to return in this case
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
result_dir = os.path.join(save_path, "result")
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
    if platform.system() == "Linux" or platform.system() == "Darwin":  # 如果是类 Unix 系统（Linux 或 macOS）
        os.chmod(result_dir, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
        uid = pwd.getpwnam("xjh").pw_uid
        gid = pwd.getpwnam("songxiehai").pw_gid
        os.chown(result_dir, uid, gid)
sc_test = sc.read_10x_mtx(path=mat)
coord = f"{mat}/barcodes_pos.tsv.gz"
sc_test = sc.read_10x_mtx(path=mat)
pos = pd.read_csv(coord, sep="\t", names=["cellID", "row", "col"])
pos2 = pos[pos["cellID"].isin(sc_test.to_df().index.to_list())]
cat_size_order = CategoricalDtype(sc_test.to_df().index.to_list(), ordered=False)
pos2['cellID'] = pos2['cellID'].astype(cat_size_order)
pos2.set_index('cellID', inplace=True)
cells_in_adata = set(sc_test.obs_names)
cells_in_pos2 = set(pos2.index)
missing_cells = cells_in_adata - cells_in_pos2
sc_test = sc_test[~sc_test.obs_names.isin(missing_cells)]
sc_test.obsm['spatial'] = np.array(pos2)
# 数据过滤步骤
min_genes = 100  # 每个细胞的最小基因数
min_cells = 3    # 每个基因的最小细胞数
sc.pp.filter_cells(sc_test, min_genes=min_genes)
sc.pp.filter_genes(sc_test, min_cells=min_cells)

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

# 运行t-SNE嵌入
sc.tl.tsne(sc_test)

# 运行Leiden聚类
sc.tl.leiden(sc_test, resolution=0.5)


# 检查基因ID文件并生成相应的图像
if os.path.exists(gene_id_file_path):
    gene_ids = []
    with open(gene_id_file_path, 'r') as file:
        next(file)  # 跳过标题行
        for line in file:
            gene_ids.append(line.strip())
sc.tl.score_genes(sc_test, gene_ids, ctrl_size=50, n_bins=25, score_name='my_score')

save_figure(sc_test, gene_ids, save_path, f'spatial_gene_projection', plot_type='umap')
save_figure(sc_test, gene_ids, save_path, f'spatial_gene_projection', plot_type='dotplot')
save_figure(sc_test, gene_ids, save_path, f'spatial_gene_projection', plot_type='spatial')
save_figure(sc_test, gene_ids, save_path, filename=f'clustered_data_gene', plot_type='cluster', Finally=True)
