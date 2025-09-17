# scripts/convert_json_to_bin.py (v2)

import json
import struct
import sys
import os

def hex_to_rgba(hex_color):
    """Converts a hex color string to a 4-byte RGBA tuple."""
    hex_color = hex_color.lstrip('#')
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    return (r, g, b, 255) # Assuming full opacity

def convert_json_to_bin_v2(json_path, bin_path):
    """
    Converts a cluster JSON file to a custom binary format (v2)
    using coordinate quantization.
    """
    with open(json_path, 'r') as f:
        data = json.load(f)

    # --- First Pass: Find global min/max for coordinates ---
    min_x, max_x, min_y, max_y = float('inf'), float('-inf'), float('inf'), float('-inf')
    
    for cluster_id, cluster_info in data.items():
        for cell in cluster_info['cells']:
            min_x = min(min_x, cell["UMAP1"])
            max_x = max(max_x, cell["UMAP1"])
            min_y = min(min_y, cell["UMAP2"])
            max_y = max(max_y, cell["UMAP2"])

    # --- Second Pass: Write the binary file ---
    clusters = sorted(data.keys(), key=int)
    num_clusters = len(clusters)

    with open(bin_path, 'wb') as f:
        # --- Write Header (v2) ---
        magic = b'PMB2'
        version = 2
        mode = 1  # Cluster mode
        f.write(struct.pack('<4sii', magic, version, mode))
        f.write(struct.pack('<ffffi', min_x, max_x, min_y, max_y, num_clusters))

        # --- Write Metadata and Data ---
        all_points_data = bytearray()
        point_start_index = 0

        for cluster_id in clusters:
            cluster_info = data[cluster_id]
            
            # --- Metadata (v2) ---
            cluster_name_raw = cluster_id.encode('utf-8')
            padded_name = cluster_name_raw.ljust(32, b'\0')
            
            rgba = hex_to_rgba(cluster_info['color'])
            color_bytes = struct.pack('<BBBB', *rgba)

            num_points = len(cluster_info['cells'])
            
            # Pack metadata without color string's padding
            f.write(struct.pack('<32s4sI', padded_name, color_bytes, num_points))
            
            # --- Quantize and pack point data ---
            for cell in cluster_info['cells']:
                quant_x = int(((cell["UMAP1"] - min_x) / (max_x - min_x)) * 65535) if max_x > min_x else 0
                quant_y = int(((cell["UMAP2"] - min_y) / (max_y - min_y)) * 65535) if max_y > min_y else 0
                all_points_data.extend(struct.pack('<HH', quant_x, quant_y))
            
        f.write(all_points_data)

    print(f"Successfully converted {json_path} to {bin_path} (v2 format)")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_json_to_bin.py <input_json_path> <output_path>")
        sys.exit(1)
    
    json_file = sys.argv[1]
    output_path = sys.argv[2]
    
    # Check if the output path is a directory
    if os.path.isdir(output_path):
        # Create the full file path automatically
        base_name = os.path.basename(json_file)
        file_name_without_ext = os.path.splitext(base_name)[0]
        bin_file = os.path.join(output_path, f"{file_name_without_ext}.bin")
    else:
        # Assume it's a full file path
        bin_file = output_path

    convert_json_to_bin_v2(json_file, bin_file)
