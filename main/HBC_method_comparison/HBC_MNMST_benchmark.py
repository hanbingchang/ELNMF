import os
from pathlib import Path
import time
import random
import warnings

warnings.filterwarnings("ignore")
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd
import scanpy as sc
import torch
import mnmstpy as mnmst

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = str(REPO_ROOT / "data" / "Human_Breast_Cancer")
OUTPUT_DIR = os.environ.get("HBC_BENCH_OUTPUT", str(REPO_ROOT / "results" / "HBC_baseline_computational_benchmark" / "labels"))

# Input preparation and spot alignment are outside the timed region.
adata = sc.read_visium(
    path=DATA_DIR,
    count_file="filtered_feature_bc_matrix.h5",
)
metadata = pd.read_csv(os.path.join(DATA_DIR, "metadata.tsv"), sep="\t").set_index("ID")
positions = pd.read_csv(
    os.path.join(DATA_DIR, "spatial", "tissue_positions_list.csv"), header=None
)
positions.columns = ["barcode", "in_tissue", "array_row", "array_col", "x", "y"]
positions = positions.set_index("barcode")
valid_ids = metadata.index[metadata["ground_truth"].notna()]
tissue_ids = positions.index[positions["in_tissue"] == 1]
valid_set = set(valid_ids).intersection(tissue_ids)
final_spots = [barcode for barcode in adata.obs_names if barcode in valid_set]
if len(final_spots) != 3798:
    raise RuntimeError(f"Expected 3798 HBC spots, found {len(final_spots)}")
adata = adata[final_spots].copy()

random.seed(1234)
np.random.seed(1234)
torch.manual_seed(1234)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(1234)
    torch.cuda.reset_peak_memory_stats()

print("BENCHMARK_START|MNMST|spots=3798|target_clusters=20", flush=True)
timer_start = time.perf_counter()

# Preserve supplied MNMST preprocessing/model parameters.
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
hvg_filter = adata.var["highly_variable"]
sc.pp.normalize_total(adata, inplace=True)
adata = adata[:, hvg_filter].copy()

enhanced_adata, cell_spatial = mnmst.data_enhance(adata, k_nei=6)
sc.pp.pca(enhanced_adata)
low_dim_x = enhanced_adata.obsm["X_pca"]
init_graph = mnmst.construct_initial_graph(low_dim_x)
init_W_tensor = torch.from_numpy(init_graph).float().to(device)
cell_spatial_tensor = torch.from_numpy(cell_spatial.toarray()).float().to(device)
x_tensor = torch.from_numpy(low_dim_x.T).to(device)
C_tensor = mnmst.sparse_self_representation_torch(
    x_tensor, init_w=init_W_tensor, device=device
)
Z_tensor = mnmst.MNMST_representation_gpu(C_tensor, cell_spatial_tensor, device=device)
Z = Z_tensor.detach().cpu().numpy()

key_added = "representation"
adata.uns[key_added] = {
    "connectivities_key": key_added,
    "distances_key": key_added,
    "var_names_use": adata.var_names.to_numpy(),
    "params": {"method": "umap"},
}
adata.obsp[key_added] = Z

# The supplied HBC-family script uses resolution 1.15 and produced 12 groups.
# Use a nearby predeclared higher resolution to approach the requested 20 groups
# without performance-based parameter selection.
mnmst_resolution = 1.80
sc.tl.leiden(
    adata,
    neighbors_key=key_added,
    resolution=mnmst_resolution,
    key_added="MNMST_cluster",
    random_state=1234,
)
labels = adata.obs["MNMST_cluster"].astype(str)

runtime_sec = time.perf_counter() - timer_start
print(
    f"BENCHMARK_END|MNMST|runtime_sec={runtime_sec:.6f}|clusters={labels.nunique()}|resolution={mnmst_resolution:.2f}",
    flush=True,
)

os.makedirs(OUTPUT_DIR, exist_ok=True)
pd.DataFrame({"Barcode": adata.obs_names, "Cluster": labels.to_numpy()}).to_csv(
    os.path.join(OUTPUT_DIR, "HBC_MNMST_labels.csv"), index=False
)

