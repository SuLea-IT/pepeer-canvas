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
DEFAULT_CSV_ROOT = "D:/xjh/laj/PM-System-Beta-Front-Beta/R/csv_exports/10DPA_L7"
# 默认的 BIN 输出目录
DEFAULT_BIN_ROOT = "D:/xjh/laj/PM-System-Beta-Front-Beta/PM-System-Beta/data/10DPA_L7"
# ============================================================

# ------------------------------------------------------------
# 单文件转换：一次性写入
# ------------------------------------------------------------
def convert_one(csv_path: str, bin_path: str) -> str:
    try:
        # 1️⃣ 先按默认类型读，保留字符串列
        df = pd.read_csv(csv_path)

        # 2️⃣ 确定 value 列
        value_col = "cluster" if "cluster" in df.columns else "expression"
        need_cols = ["x", "y", value_col]
        if not set(need_cols) <= set(df.columns):
            return f"缺列: {csv_path}"

        # 3️⃣ 逐列数值化
        for col in ("x", "y"):
            df[col] = pd.to_numeric(df[col], errors="coerce")

        if value_col == "cluster":
            # cluster 可能是 L7_0 这一类符号
            def _to_float(v):
                try:
                    return float(v)
                except ValueError:
                    m = re.search(r"-?\d+\.?\d*", str(v))
                    return float(m.group()) if m else np.nan
            df[value_col] = df[value_col].apply(_to_float)
        else:
            df[value_col] = pd.to_numeric(df[value_col], errors="coerce")

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

        genes_dir = os.path.join(input_path, "genes")
        genes_bin_dir = os.path.join(output_path, "genes")
        if os.path.isdir(genes_dir):
            for fn in os.listdir(genes_dir):
                if fn.endswith(".csv"):
                    tasks.append((
                        os.path.join(genes_dir, fn),
                        os.path.join(genes_bin_dir, fn[:-4] + ".bin")
                    ))
    else:
        print(f"错误: 输入路径 '{input_path}' 不是一个有效的 CSV 文件或目录。")
        return

    if not tasks:
        print("未找到可转换的 CSV 文件。")
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
# CLI
# ------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        print("用法: python convert_csv_to_bin.py <input_path> <output_path> [workers]")
        print("  <input_path>: 可以是单个 .csv 文件或包含 .csv 文件的目录。")
        print("  <output_path>: 可以是单个 .bin 文件名或输出目录。")
        sys.exit(1)

    # 如果提供了命令行参数，则使用它们；否则使用文件顶部的默认值
    if len(sys.argv) >= 3:
        input_path, output_path = sys.argv[1:3]
        workers_arg = int(sys.argv[3]) if len(sys.argv) == 4 else None
        print("使用命令行参数指定的路径。")
    else:
        input_path = DEFAULT_CSV_ROOT
        output_path = DEFAULT_BIN_ROOT
        workers_arg = None  # 使用默认的 CPU 核心数
        print("未提供命令行参数，使用脚本中定义的默认路径。")

    print(f"  - 输入路径: {input_path}")
    print(f"  - 输出路径: {output_path}")
    
    main(input_path, output_path, workers_arg)
