import sys
import pandas as pd
import struct

def main():
    # 获取命令行参数
    mat = sys.argv[1]         # 输入的CSV文件路径
    save_path = sys.argv[2]   # 输出文件保存的目录
    gene = sys.argv[3]        # 要保留的基因列名称
    color = "#460b5e"         # 颜色代码

    cols_to_keep = ['UMAP1', 'UMAP2', 'cell', gene]
    df = pd.read_csv(mat, usecols=cols_to_keep)

    # 准备二进制数据
    # 文件头
    magic = b'PMB1'
    version = 1
    mode = 2  # Gene mode
    num_groups = 1
    header = struct.pack('<4siii', magic, version, mode, num_groups)
    
    # 元数据
    group_name = gene.encode('utf-8').ljust(32, b'\\0')
    color_str = color.encode('utf-8').ljust(24, b'\\0')
    point_start_index = 0
    num_points = len(df)
    meta = struct.pack('<32s24sii', group_name, color_str, point_start_index, num_points)

    # 数据点
    points_data = bytearray()
    for index, row in df.iterrows():
        points_data.extend(struct.pack('<fff', float(row["UMAP1"]), float(row["UMAP2"]), float(row[gene])))

    # 将所有部分写入标准输出
    sys.stdout.buffer.write(header)
    sys.stdout.buffer.write(meta)
    sys.stdout.buffer.write(points_data)

if __name__ == "__main__":
    main()
