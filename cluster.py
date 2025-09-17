import sys
import pandas as pd
import numpy as np
import gzip
from pandas.api.types import CategoricalDtype
import scanpy as sc
import os
import matplotlib.pyplot as plt
import squidpy as sq
import time
import leidenalg
import igraph as ig

class LeidenCluster:
    def __init__(self, adata):
        self.adata = adata
        # 缓存图结构
        self.graph = None
        self.partition_type = leidenalg.RBConfigurationVertexPartition
        
    def initialize_graph(self):
        # 只在第一次调用时初始化图
        if self.graph is None:
            neighbors = self.adata.obsp['connectivities']
            # 转换为 igraph 格式
            sources, targets = neighbors.nonzero()
            weights = neighbors[sources, targets]
            self.graph = ig.Graph(
                n=neighbors.shape[0],
                edges=list(zip(sources, targets)),
                edge_attrs={'weight': weights}
            )
    
    def cluster(self, resolution: float) -> np.ndarray:
        self.initialize_graph()
        partition = leidenalg.find_partition(
            self.graph,
            self.partition_type,
            resolution_parameter=resolution,
        )
        return np.array(partition.membership)

# 主要处理流程
def process_data(mat_path, save_path):
    # 读取数据
    adata = sc.read_10x_mtx(path=mat_path)
    
    # 过滤数据
    min_genes = 100
    min_cells = 3
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    # 检查是否已经计算了UMAP和t-SNE
    has_umap = 'X_umap' in adata.obsm
    has_tsne = 'X_tsne' in adata.obsm
    
    if not has_umap and not has_tsne:
        # 数据预处理步骤
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, inplace=True)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
        sc.pp.pca(adata)
        
    # 确保计算neighbors - 添加调试信息
    print("计算邻居关系...")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
    
    # 添加调试检查
    print("检查 adata.obsp 的键：", list(adata.obsp.keys()))
    
    if not has_umap and not has_tsne:
        start_time = time.time()
        sc.tl.umap(adata)
        print(f"UMAP 耗时: {time.time() - start_time:.2f} 秒")
        sc.tl.tsne(adata)
    
    # 再次检查
    print("聚类前检查 adata.obsp 的键：", list(adata.obsp.keys()))
    
    # 初始化 LeidenCluster
    leiden_cluster = LeidenCluster(adata)
    
    # 测试两个不同resolution的聚类时间
    print("开始Leiden聚类性能测试...")
    
    # 第一次聚类 (resolution=0.5) - 包含图的初始化时间
    start_time = time.time()
    clusters_0_5 = leiden_cluster.cluster(resolution=0.5)
    time_0_5 = time.time() - start_time
    print(f"Resolution=0.5 总耗时: {time_0_5:.2f} 秒")
    
    # 第二次聚类 (resolution=0.4) - 复用已有图结构
    start_time = time.time()
    clusters_0_4 = leiden_cluster.cluster(resolution=0.4)
    time_0_4 = time.time() - start_time
    print(f"Resolution=0.4 耗时: {time_0_4:.2f} 秒")
    
    # 使用resolution=0.5的结果作为最终聚类结果
    adata.obs['leiden'] = pd.Categorical(clusters_0_5.astype(str))
    
    # 创建 UMAP DataFrame
    umap_coordinates = adata.obsm['X_umap']
    umap_df = pd.DataFrame(umap_coordinates, columns=['UMAP1', 'UMAP2'])
    umap_df['cell'] = adata.obs.index
    umap_df['cluster'] = adata.obs['leiden'].values
    
    # 保存结果
    umap_df.to_csv(os.path.join(save_path, '10DPA.csv'), index=False)
    
    return adata, leiden_cluster

if __name__ == "__main__":
    mat = "/home/xjh/gene_xy/data/zyx"
    save_path = "/home/xjh/gene_xy/data/zyx"
    
    # 处理数据并获取聚类器实例
    adata, leiden_cluster = process_data(mat, save_path)
    
    # 示例：如何快速更改 resolution 并获取新的聚类结果
    def update_clustering(resolution):
        start_time = time.time()
        new_clusters = leiden_cluster.cluster(resolution)
        print(f"更新 resolution={resolution} 耗时: {time.time() - start_time:.2f} 秒")
        return new_clusters
