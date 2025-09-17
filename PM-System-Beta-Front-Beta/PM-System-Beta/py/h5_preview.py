import sys
import h5py
import json
import numpy as np

def get_h5_preview(file_path):
    structure = {}
    preview = {}
    
    with h5py.File(file_path, 'r') as f:
        def get_structure(name, obj):
            if isinstance(obj, h5py.Dataset):
                structure[name] = {
                    "type": "dataset",
                    "shape": str(obj.shape),
                    "dtype": str(obj.dtype)
                }
                
                try:
                    if len(obj.shape) == 0:
                        preview[name] = obj[()].item() if isinstance(obj[()], (np.integer, np.floating)) else str(obj[()])
                    elif len(obj.shape) == 1:
                        preview[name] = obj[:10].tolist() if obj.shape[0] > 10 else obj[:].tolist()
                    elif len(obj.shape) == 2:
                        preview[name] = obj[:5, :5].tolist() if obj.shape[0] > 5 or obj.shape[1] > 5 else obj[:, :].tolist()
                    else:
                        preview[name] = obj.flatten()[:10].tolist()
                except Exception as e:
                    preview[name] = f"Preview failed: {str(e)}"
            else:
                structure[name] = {
                    "type": "group"
                }
        
        f.visititems(get_structure)
        
    return json.dumps({
        "structure": structure,
        "preview": preview
    })

if __name__ == '__main__':
    file_path = sys.argv[1]
    print(get_h5_preview(file_path))