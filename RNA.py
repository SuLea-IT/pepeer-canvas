# 1. 导入必要的库
import os
import glob
import scanpy as sc
import pandas as pd
import numpy as np
import squidpy as sq
import dynamo as dyn
import warnings
from pandas.core.frame import DataFrame
from numba import NumbaDeprecationWarning

# 忽略警告信息
warnings.filterwarnings("ignore", category=NumbaDeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# 2. 数据加载和预处理
def load_loom_data(directory_path):
    """加载.loom文件并进行初始处理"""
    loom_files = glob.glob(os.path.join(directory_path, '*.loom'))
    if len(loom_files) == 0:
        raise FileNotFoundError("没有找到.loom文件")
    elif len(loom_files) == 1:
        loom_data = sc.read(loom_files[0], cache=False)
        # 处理barcode名称
        loom_data.obs = loom_data.obs.rename(index=lambda x: x.split(':', 1)[-1].replace('L1', 'L1'))
        return loom_data
    else:
        raise ValueError(f"找到多个.loom文件: {loom_files}")

# 3. 数据整合函数
def integrate_data(loom_data, meta_path):
    """整合loom数据与meta信息"""
    # 读取meta数据
    sample_obs = pd.read_csv(os.path.join(meta_path, "cell ID_obs.csv"))
    cell_umap = pd.read_csv(os.path.join(meta_path, "cell_embeddings.csv"), 
                           header=0, 
                           names=["Cell ID", "UMAP_1", "UMAP_2"])
    cell_spatial = pd.read_csv(os.path.join(meta_path, "cell_spatial.csv"), 
                              header=0, 
                              names=["Cell ID", "spatial_1", "spatial_2"])
    cell_clusters = pd.read_csv(os.path.join(meta_path, "cell_clusters.csv"), 
                               header=0, 
                               names=["Cell ID", "seurat_clusters"])
    
    # 数据整合
    sample_one = loom_data[np.isin(loom_data.obs.index, sample_obs)]
    
    # 确保所有数据使用相同的细胞ID索引
    common_cells = set(sample_one.obs.index).intersection(
        set(cell_umap["Cell ID"]),
        set(cell_spatial["Cell ID"]),
        set(cell_clusters["Cell ID"])
    )
    
    # 过滤数据，只保留共同的细胞
    sample_one = sample_one[sample_one.obs.index.isin(common_cells)]
    cell_umap = cell_umap[cell_umap["Cell ID"].isin(common_cells)]
    cell_spatial = cell_spatial[cell_spatial["Cell ID"].isin(common_cells)]
    cell_clusters = cell_clusters[cell_clusters["Cell ID"].isin(common_cells)]
    
    # 确保数据按照相同顺序排列
    cell_umap = cell_umap.set_index("Cell ID").loc[sample_one.obs.index]
    cell_spatial = cell_spatial.set_index("Cell ID").loc[sample_one.obs.index]
    cell_clusters = cell_clusters.set_index("Cell ID").loc[sample_one.obs.index]
    
    return sample_one, cell_umap, cell_spatial, cell_clusters

# 4. RNA速率分析函数
def perform_velocity_analysis(adata, save_path='result'):
    """执行RNA速率分析"""
    # 预处理
    dyn.pp.recipe_monocle(adata)
    dyn.tl.dynamics(adata, model='stochastic', cores=3)
    
    # 降维和可视化
    adata.obsm['X_umap_b'] = adata.obsm['X_umap']
    dyn.tl.reduceDimension(adata, enforce=True)
    adata.obsm['X_umap'] = adata.obsm['X_umap_b']

    
    # 计算细胞速度
    dyn.tl.cell_velocities(adata, 
                          method="pearson",
                          enforce=True,
                          transition_genes=list(adata.var_names[adata.var.use_for_pca]))
    
    # 移除 velocity_embedding 调用，改用 dyn.tl.cell_velocities 的结果
    # 速度信息已经包含在 adata.obsm['velocity'] 中
    
    dyn.tl.cell_wise_confidence(adata)
    
    return adata

# 5. 可视化函数
def create_visualizations(adata, save_path='result'):
    """创建各种可视化图表"""
    os.makedirs(save_path, exist_ok=True)
    sc.settings.figdir = save_path
    
    # 设置 dynamo 的图片保存路径
    dyn.configuration.set_figure_params('dynamo')
    dyn.configuration.save_fig_path = save_path

    # 基础空间和UMAP可视化
    sc.pl.spatial(adata, 
                 color='seurat_clusters', 
                 size=5, 
                 spot_size=3,
                 save='spatial_clusters.pdf')
    
    sc.pl.umap(adata, 
               color='seurat_clusters', 
               size=20,
               save='umap_clusters.pdf')
    
    # UMAP向量场可视化
    try:
        dyn.pl.streamline_plot(adata, 
                              color=['seurat_clusters'], 
                              basis='umap',
                              show_legend='on data',
                              show_arrowed_spines=True,
                              pointsize=1,
                              alpha=1,
                              figsize=(10,8),
                              background='white',
                              save_show_or_return='save',
                              save_kwargs={"path": os.path.join(save_path, 'results'), 
                                         "prefix": 'streamline_umap', 
                                         "dpi": 900, 
                                         "ext": 'pdf'})
    except Exception as e:
        print(f"生成UMAP向量场时出错: {str(e)}")

    # 空间向量场可视化
    try:
        # 直接使用 cell_wise_vectors 和 streamline_plot
        print(adata)
        dyn.pl.cell_wise_vectors(adata, 
                                color=['seurat_clusters'], 
                                basis='spatial',
                                show_legend='on data',
                                quiver_length=2,
                                quiver_size=2,
                                pointsize=0.3,
                                show_arrowed_spines=False,
                                figsize=(10,16),
                                save_show_or_return='save',
                                save_kwargs={"path": os.path.join(save_path, 'results'), 
                                           "prefix": 'spatial_vectors', 
                                           "dpi": 900, 
                                           "ext": 'pdf'})
        
        dyn.pl.streamline_plot(adata, 
                              color=['seurat_clusters'], 
                              basis='spatial',
                              show_legend='on data',
                              show_arrowed_spines=True,
                              pointsize=0.5,
                              linewidth=2,
                              streamline_alpha=1,
                              density=1.5,
                              figsize=(10,16),
                              background='white',
                              save_show_or_return='save',
                              save_kwargs={"path": os.path.join(save_path, 'results'), 
                                         "prefix": 'streamline_spatial', 
                                         "dpi": 900, 
                                         "ext": 'png'})
    except Exception as e:
        print(f"生成空间向量场时出错: {str(e)}")
        print("错误的详细信息：")
        import traceback
        traceback.print_exc()

def part_visualization(adata, save_path='result'):
    """部分可视化处理"""
    # 设置 dynamo 的图片保存路径
    dyn.configuration.set_figure_params('dynamo')
    dyn.configuration.save_fig_path = save_path
    
    adata.obs.seurat_clusters = adata.obs.seurat_clusters.astype('str')
    
    dyn.pl.umap(adata, color='seurat_clusters',pointsize=0.1,alpha=0.7,show_legend="on data"
                ,figsize= (10,16),save_show_or_return='save',
             save_kwargs={"path": save_path, "prefix": 'dynUMAP', "dpi": None, "ext": 'pdf'})
    # 将 'spatial' 重命名为 'X_spatial'
    adata.obsm['X_spatial'][:,1]=adata.obsm['X_spatial'][:,1]*(-1)
    
    dyn.pl.space(adata, color='seurat_clusters',pointsize=0.45,alpha=0.7,show_legend="on data",
             figsize= (10,16),
             save_show_or_return='save',
             save_kwargs={"path": os.path.join(save_path, 'results'), "prefix": 'dynamo', "dpi": None, "ext": 'pdf'})
    
    sc.pp.calculate_qc_metrics(adata, expr_type='unspliced', layer="unspliced", inplace=True)
    sc.pp.calculate_qc_metrics(adata, expr_type='spliced', layer="spliced", inplace=True)
    return adata

# 6. 主函数
def main():
    """主函数：执行完整的数据处理和分析流程"""
    # 设置数据路径
    base_dir = '/mnt/newdisk/xjh/RNA_data'
    meta_dir = os.path.join(base_dir, 'meta')
    result_dir = os.path.join(base_dir, 'results')
    
    # 确保结果目录存在
    os.makedirs(result_dir, exist_ok=True)
    
    try:
        # 1. 加载loom数据
        print("正在加载loom数据...")
        adata = load_loom_data(base_dir)
        print(f"成功加载数据，形状为: {adata.shape}")

        # 2. 整合meta数据
        print("正在整合meta数据...")
        adata, cell_umap, cell_spatial, cell_clusters = integrate_data(adata, meta_dir)
        
        
        # 将meta数据添加到adata对象中
        adata.obsm['X_umap'] = cell_umap[['UMAP_1', 'UMAP_2']].values
        adata.obsm['X_spatial'] = cell_spatial[['spatial_1', 'spatial_2']].values  # 直接存储为 X_spatial
        adata.obsm['spatial'] = adata.obsm['X_spatial'].copy()  # 保留一个副本作为 spatial
        adata.obs['seurat_clusters'] = cell_clusters['seurat_clusters']
        print("数据整合完成")
        adata = part_visualization(adata,save_path=result_dir)

        # # 3. 执行RNA速率分析
        print("正在进行RNA速率分析...")
        adata = perform_velocity_analysis(adata, save_path=result_dir)
        print("RNA速率分析完成")

        # 4. 生成可视化结果
        print("正在生成可视化结果...")
        create_visualizations(adata, save_path=result_dir)
        print(f"分析完成！结果已保存到: {result_dir}")

    except Exception as e:
        print(f"处理过程中出现错误: {str(e)}")

if __name__ == "__main__":
    main()
