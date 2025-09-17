import os
import sys
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from typing import Optional
from tqdm import tqdm
import re

# ============================================================
# 配置
# ============================================================
# 默认的 CSV 输入目录
DEFAULT_CSV_ROOT = "E:/Downloads/object_unintegrated/clusters.csv"
# 默认的 BIN 输出目录
DEFAULT_BIN_ROOT = "E:/Downloads/object_unintegrated/"
# ============================================================

# ------------------------------------------------------------
# 单文件转换：一次性写入
# ------------------------------------------------------------
def convert_one(csv_path: str, bin_path: str) -> str:
    try:
        # 1️⃣ 读取CSV文件
        df = pd.read_csv(csv_path)

        # 2️⃣ 检查必需的列
        need_cols = ["x", "y", "cluster"]
        if not set(need_cols) <= set(df.columns):
            # 如果没有x,y,cluster列，检查是否存在barcode和cluster列（非空间数据）
            if set(["barcode", "cluster"]) <= set(df.columns):
                # 对于非空间数据，我们只关心barcode和cluster
                df_out = df[["barcode", "cluster"]]
                # 在这种情况下，我们可能需要不同的二进制格式或只是保存为不同的csv/json
                # 为了与现有流程兼容，这里只打印一条消息，实际应用中可能需要不同的处理
                # return f"非空间数据，跳过二进制转换: {csv_path}"
                # 或者，如果只需要cluster值，可以这样做：
                arr = df[["cluster"]].to_numpy(np.float32, copy=False)
            else:
                return f"缺列: {csv_path} (需要 'x', 'y', 'cluster' 或 'barcode', 'cluster')"


        # 3️⃣ 逐列数值化
        if "x" in df.columns and "y" in df.columns:
            for col in ("x", "y"):
                df[col] = pd.to_numeric(df[col], errors="coerce")

        # cluster 可能是 L7_0 这一类符号
        def _to_float(v):
            try:
                return float(v)
            except (ValueError, TypeError):
                m = re.search(r"-?\d+\.?\d*", str(v))
                return float(m.group()) if m else np.nan
        df["cluster"] = df["cluster"].apply(_to_float)

        # 4️⃣ 去掉任何包含 NaN 的行
        df = df.dropna(subset=need_cols)
        if df.empty:
            return f"{csv_path} 没有有效数据行"

        # 5️⃣ 写二进制
        arr = df[need_cols].to_numpy(np.float32, copy=False)
        os.makedirs(os.path.dirname(bin_path), exist_ok=True)
        arr.tofile(bin_path)
        return ""

    except Exception as e:
        return f"{csv_path} 出错: {e}"

# ------------------------------------------------------------
# 主流程
# ------------------------------------------------------------
def main(input_path: str, output_path: str, workers: Optional[int] = None):
    input_path, output_path = map(os.path.abspath, (input_path, output_path))

    # 收集任务
    tasks = []

    # 检查输入是文件还是目录
    if os.path.isfile(input_path) and input_path.endswith('.csv'):
        # --- 单文件模式 ---
        print("检测到单文件模式。")
        # 如果输出路径是目录，或没有扩展名，则视为目录
        if os.path.isdir(output_path) or not os.path.splitext(output_path)[1]:
             os.makedirs(output_path, exist_ok=True)
             bin_path = os.path.join(output_path, os.path.basename(input_path)[:-4] + ".bin")
        else: # 否则视为文件
            bin_path = output_path
            os.makedirs(os.path.dirname(bin_path), exist_ok=True)
        tasks.append((input_path, bin_path))

    elif os.path.isdir(input_path):
        # --- 目录模式 ---
        print("检测到目录模式。")
        os.makedirs(output_path, exist_ok=True)
        
        cluster_csv = os.path.join(input_path, "clusters.csv")
        if os.path.exists(cluster_csv):
            tasks.append((cluster_csv, os.path.join(output_path, "clusters.bin")))
        else:
            print(f"警告: 在 '{input_path}' 中未找到 'clusters.csv'。")

    else:
        print(f"错误: 输入路径 '{input_path}' 不是一个有效的 CSV 文件或目录。")
        return

    if not tasks:
        print("未找到可转换的 'clusters.csv' 文件。")
        return

    # 并行执行
    err_list = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(convert_one, src, dst): src for src, dst in tasks}
        for fut in tqdm(as_completed(futures), total=len(futures), desc="转换"):
            msg = fut.result()
            if msg:
                err_list.append(msg)

    # 汇总
    if err_list:
        print("\n=== 有问题的文件 ===")
        for m in err_list:
            print(m)
    print(f"\n全部完成：成功 {len(tasks)-len(err_list)} / {len(tasks)}")

# ------------------------------------------------------------
# 脚本执行
# ------------------------------------------------------------
if __name__ == "__main__":
    # 使用脚本顶部定义的默认路径
    input_path = DEFAULT_CSV_ROOT
    output_path = DEFAULT_BIN_ROOT
    workers_arg = None  # 使用默认的 CPU 核心数

    print("使用脚本中定义的固定路径:")
    print(f"  - 输入路径: {input_path}")
    print(f"  - 输出路径: {output_path}")
    
    # 检查默认路径是否存在
    if not os.path.exists(input_path):
        print(f"错误: 默认输入路径不存在: {input_path}")
        sys.exit(1)

    main(input_path, output_path, workers_arg)