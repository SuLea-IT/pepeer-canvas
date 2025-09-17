from pandas.core.frame import DataFrame
from numba import NumbaDeprecationWarning
import os
os.getcwd()
os.chdir('/mnt/newdisk/xjh/RNA_V/pepper/output')
#import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import squidpy as sq
import dynamo as dyn
import warnings
warnings.filterwarnings("ignore", category=NumbaDeprecationWarning)

warnings.filterwarnings("ignore", category=FutureWarning)
import os
import glob

# 指定目录路径
directory_path = '/mnt/newdisk/RNA_V/pepper/output'
save_path = '/mnt/newdisk/xjh/RNA_data/results'
#
#RNA速率分析
#------------------------------------------------------------------

adata= sc.read_h5ad("/mnt/newdisk/xjh/RNA_data/All_dynamicModel_epidermis.h5ad")
adata_spl_sel = adata
adata_spl_sel
##作图代码调试
#自定义分类颜色字典
ident_color={"4":"#7495D3","6":"#6D1A9C","7":"#15821E"}

adata_spl_sel.obs.seurat_clusters = adata_spl_sel.obs.seurat_clusters.astype('str')
sc.pl.spatial(adata_spl_sel, color = 'seurat_clusters', size=5, spot_size=3)
sc.pl.umap(adata_spl_sel,color = 'seurat_clusters',size=20)
dyn.pl.umap(adata_spl_sel, color='seurat_clusters',pointsize=0.1,alpha=0.7,show_legend="on data")
adata_spl_sel.obsm['X_spatial'][:,1]=adata_spl_sel.obsm['X_spatial'][:,1]*(-1)
dyn.pl.space(adata_spl_sel, color='seurat_clusters',pointsize=0.45,alpha=0.7,show_legend="on data",
             figsize= (10,16),
             save_show_or_return='save',
             save_kwargs={"path": save_path, "prefix": 'dynamo', "dpi": None, "ext": 'pdf'})
             # 计算剪切与未剪切的 QC 指标
sc.pp.calculate_qc_metrics(adata_spl_sel, expr_type='unspliced', layer="unspliced", inplace=True)
sc.pp.calculate_qc_metrics(adata_spl_sel, expr_type='spliced', layer="spliced", inplace=True)

# 绘制小提琴图，使用 density_norm='width'
sc.pl.violin(adata_spl_sel, ["n_genes_by_spliced", "n_genes_by_unspliced"], density_norm='width')
sc.pl.violin(adata_spl_sel, ["log1p_total_spliced", "log1p_total_unspliced"], density_norm='width')
#标准化数据
dyn.pp.recipe_monocle(adata_spl_sel) #获得ntr,X_pca
adata_spl_sel
dyn.tl.dynamics(adata_spl_sel,model='stochastic', cores=3)
#UMAP降维
adata_spl_sel.obsm['X_umap_b'] = adata_spl_sel.obsm['X_umap']
dyn.tl.reduceDimension(adata_spl_sel, enforce = True)
adata_spl_sel.obsm['X_umap'] = adata_spl_sel.obsm['X_umap_b'] #替换为seurat聚类结果
adata_spl_sel.obs.seurat_clusters = adata_spl_sel.obs.seurat_clusters.astype('str')
dyn.pl.umap(adata_spl_sel, color='seurat_clusters',pointsize=0.1,alpha=0.7,show_legend="on data")
#dyn.tl.gene_wise_confidence(adata_spl_sel, group='seurat_clusters', lineage_dict={'6': ['4']})
dyn.pl.phase_portraits(adata_spl_sel, genes=adata.var_names[adata.var.use_for_dynamics][:4], figsize=(6, 4), color='seurat_clusters')
#速度投影
dyn.tl.cell_velocities(adata_spl_sel, method="pearson",enforce=True,  transition_genes = list(adata_spl_sel.var_names[adata_spl_sel.var.use_for_pca]))
#dyn.tl.cell_velocities(adata, method='pearson', other_kernels_dict={'transform': 'sqrt'})
dyn.tl.cell_wise_confidence(adata_spl_sel)
dyn.pl.cell_wise_vectors(adata_spl_sel, color=['seurat_clusters'], basis='umap', show_legend='on data', quiver_length=5, quiver_size=5, pointsize=1,alpha=1, show_arrowed_spines=False,
                         figsize= (10,8))
dyn.pl.cell_wise_vectors(adata_spl_sel, color=['seurat_clusters'], basis='spatial', show_legend='on data', quiver_length=2, quiver_size=2, pointsize=0.3, show_arrowed_spines=False,
                         figsize= (10,16))
dyn.pl.streamline_plot(adata_spl_sel, color=['seurat_clusters'], basis='umap', show_legend='on data', show_arrowed_spines=True,pointsize=1,alpha=1)
dyn.pl.streamline_plot(adata_spl_sel, color=['seurat_clusters'], basis='umap', show_legend='on data', show_arrowed_spines=True,pointsize=1,alpha=1,
                       figsize= (10,8),
                       save_show_or_return='save',
                       save_kwargs={"path": save_path, "prefix": 'streamline_umap', "dpi": 900, "ext": 'pdf'})
dyn.pl.streamline_plot(adata_spl_sel, color=['seurat_clusters'], basis='spatial', show_legend='on data', show_arrowed_spines=True, pointsize=0.5,
                       linewidth=2,
                       streamline_alpha=1,
                       density=1.5,
                       figsize= (10,16),
                       save_show_or_return='save',
                       background= 'white',
                       save_kwargs={"path": save_path, "prefix": 'streamline_spatial', "dpi": 900, "ext": 'png'})
