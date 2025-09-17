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
# 閰嶇疆
# ============================================================
# 榛樿鐨?CSV 杈撳叆鐩綍
DEFAULT_CSV_ROOT = "D:/xjh/laj/PM-System-Beta-Front-Beta/R/csv_exports/10DPA_L7"
# 榛樿鐨?BIN 杈撳嚭鐩綍
DEFAULT_BIN_ROOT = "D:/xjh/laj/PM-System-Beta/data/10DPA_L7"
# ============================================================

# ------------------------------------------------------------
# 鍗曟枃浠惰浆鎹細涓€娆℃€у啓鍏?
# ------------------------------------------------------------
def convert_one(csv_path: str, bin_path: str) -> str:
    try:
        # 1锔忊儯 鍏堟寜榛樿绫诲瀷璇伙紝淇濈暀瀛楃涓插垪
        df = pd.read_csv(csv_path)

        # 2锔忊儯 纭畾 value 鍒?
        value_col = "cluster" if "cluster" in df.columns else "expression"
        need_cols = ["x", "y", value_col]
        if not set(need_cols) <= set(df.columns):
            return f"缂哄垪: {csv_path}"

        # 3锔忊儯 閫愬垪鏁板€煎寲
        for col in ("x", "y"):
            df[col] = pd.to_numeric(df[col], errors="coerce")

        if value_col == "cluster":
            # cluster 鍙兘鏄?L7_0 杩欎竴绫荤鍙?
            def _to_float(v):
                try:
                    return float(v)
                except ValueError:
                    m = re.search(r"-?\d+\.?\d*", str(v))
                    return float(m.group()) if m else np.nan
            df[value_col] = df[value_col].apply(_to_float)
        else:
            df[value_col] = pd.to_numeric(df[value_col], errors="coerce")

        # 4锔忊儯 鍘绘帀浠讳綍鍖呭惈 NaN 鐨勮
        df = df.dropna(subset=need_cols)
        if df.empty:
            return f"{csv_path} 娌℃湁鏈夋晥鏁版嵁琛?

        # 5锔忊儯 鍐欎簩杩涘埗
        arr = df[need_cols].to_numpy(np.float32, copy=False)
        os.makedirs(os.path.dirname(bin_path), exist_ok=True)
        arr.tofile(bin_path)
        return ""

    except Exception as e:
        return f"{csv_path} 鍑洪敊: {e}"

# ------------------------------------------------------------
# 涓绘祦绋?
# ------------------------------------------------------------
def main(input_path: str, output_path: str, workers: Optional[int] = None):
    input_path, output_path = map(os.path.abspath, (input_path, output_path))

    # 鏀堕泦浠诲姟
    tasks = []

    # 妫€鏌ヨ緭鍏ユ槸鏂囦欢杩樻槸鐩綍
    if os.path.isfile(input_path) and input_path.endswith('.csv'):
        # --- 鍗曟枃浠舵ā寮?---
        print("妫€娴嬪埌鍗曟枃浠舵ā寮忋€?)
        # 濡傛灉杈撳嚭璺緞鏄洰褰曪紝鎴栨病鏈夋墿灞曞悕锛屽垯瑙嗕负鐩綍
        if os.path.isdir(output_path) or not os.path.splitext(output_path)[1]:
             os.makedirs(output_path, exist_ok=True)
             bin_path = os.path.join(output_path, os.path.basename(input_path)[:-4] + ".bin")
        else: # 鍚﹀垯瑙嗕负鏂囦欢
            bin_path = output_path
            os.makedirs(os.path.dirname(bin_path), exist_ok=True)
        tasks.append((input_path, bin_path))

    elif os.path.isdir(input_path):
        # --- 鐩綍妯″紡 ---
        print("妫€娴嬪埌鐩綍妯″紡銆?)
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
        print(f"閿欒: 杈撳叆璺緞 '{input_path}' 涓嶆槸涓€涓湁鏁堢殑 CSV 鏂囦欢鎴栫洰褰曘€?)
        return

    if not tasks:
        print("鏈壘鍒板彲杞崲鐨?CSV 鏂囦欢銆?)
        return

    # 骞惰鎵ц
    err_list = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(convert_one, src, dst): src for src, dst in tasks}
        for fut in tqdm(as_completed(futures), total=len(futures), desc="杞崲"):
            msg = fut.result()
            if msg:
                err_list.append(msg)

    # 姹囨€?
    if err_list:
        print("\n=== 鏈夐棶棰樼殑鏂囦欢 ===")
        for m in err_list:
            print(m)
    print(f"\n鍏ㄩ儴瀹屾垚锛氭垚鍔?{len(tasks)-len(err_list)} / {len(tasks)}")

# ------------------------------------------------------------
# CLI
# ------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        print("鐢ㄦ硶: python convert_csv_to_bin.py <input_path> <output_path> [workers]")
        print("  <input_path>: 鍙互鏄崟涓?.csv 鏂囦欢鎴栧寘鍚?.csv 鏂囦欢鐨勭洰褰曘€?)
        print("  <output_path>: 鍙互鏄崟涓?.bin 鏂囦欢鍚嶆垨杈撳嚭鐩綍銆?)
        sys.exit(1)

    # 濡傛灉鎻愪緵浜嗗懡浠よ鍙傛暟锛屽垯浣跨敤瀹冧滑锛涘惁鍒欎娇鐢ㄦ枃浠堕《閮ㄧ殑榛樿鍊?
    if len(sys.argv) >= 3:
        input_path, output_path = sys.argv[1:3]
        workers_arg = int(sys.argv[3]) if len(sys.argv) == 4 else None
        print("浣跨敤鍛戒护琛屽弬鏁版寚瀹氱殑璺緞銆?)
    else:
        input_path = DEFAULT_CSV_ROOT
        output_path = DEFAULT_BIN_ROOT
        workers_arg = None  # 浣跨敤榛樿鐨?CPU 鏍稿績鏁?
        print("鏈彁渚涘懡浠よ鍙傛暟锛屼娇鐢ㄨ剼鏈腑瀹氫箟鐨勯粯璁よ矾寰勩€?)

    print(f"  - 杈撳叆璺緞: {input_path}")
    print(f"  - 杈撳嚭璺緞: {output_path}")
    
    main(input_path, output_path, workers_arg)
