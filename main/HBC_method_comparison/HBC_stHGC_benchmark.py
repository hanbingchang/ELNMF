import os
from pathlib import Path
import sys
import time
import random
import warnings

warnings.filterwarnings("ignore")
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from scipy.spatial import cKDTree

# Set STHGC_SOURCE_DIR to the directory containing stHGC.py and utils.py,
# or install the same stHGC version used in the experiment into this Python environment.
SOURCE_DIR = os.environ.get("STHGC_SOURCE_DIR")
if SOURCE_DIR:
    sys.path.insert(0, SOURCE_DIR)
import stHGC
from utils import *

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = str(REPO_ROOT / "data" / "Human_Breast_Cancer")
OUTPUT_DIR = os.environ.get("HBC_BENCH_OUTPUT", str(REPO_ROOT / "results" / "HBC_baseline_computational_benchmark" / "labels"))

# Input preparation and deterministic alignment are outside the timed region.
adata = sc.read_10x_h5(os.path.join(DATA_DIR, "filtered_feature_bc_matrix.h5"))
adata.var_names_make_unique()
positions = pd.read_csv(
    os.path.join(DATA_DIR, "spatial", "tissue_positions_list.csv"), header=None
)
positions.columns = ["barcode", "in_tissue", "array_row", "array_col", "x", "y"]
positions = positions.set_index("barcode")
metadata = pd.read_csv(os.path.join(DATA_DIR, "metadata.tsv"), sep="\t").set_index("ID")
valid_ids = metadata.index[metadata["ground_truth"].notna()]
tissue_ids = positions.index[positions["in_tissue"] == 1]
valid_set = set(valid_ids).intersection(tissue_ids)
final_spots = [barcode for barcode in adata.obs_names if barcode in valid_set]
if len(final_spots) != 3798:
    raise RuntimeError(f"Expected 3798 HBC spots, found {len(final_spots)}")
adata = adata[final_spots].copy()
adata.obsm["spatial"] = positions.loc[final_spots, ["x", "y"]].to_numpy()

random.seed(2020)
np.random.seed(2020)
print("BENCHMARK_START|stHGC|spots=3798|target_clusters=20", flush=True)
timer_start = time.perf_counter()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(
    adata, flavor="seurat_v3", n_top_genes=5000, batch_key=None
)

coords = adata.obsm["spatial"]
tree = cKDTree(coords)
distances, _ = tree.query(coords, k=2)
mean_nn_dist = distances[:, 1].mean()
rad_cutoff = max(mean_nn_dist * 3, 30)
stHGC.cal_hybrid_graph(adata, rad_cutoff=rad_cutoff)
stHGC.stats_hybrid_graph(adata)

# Preserve supplied model defaults: [512, 30], 1000 epochs, lr=0.0009,
# weight_decay=0.0001, training seed=2020 and automatic CPU device selection.
adata = stHGC.stHGC_train(adata)
X_embed = adata.obsm["stHGC"].astype(np.float64)
n_pca = min(30, X_embed.shape[0] - 1, X_embed.shape[1])
X_pca = PCA(n_components=n_pca, random_state=0).fit_transform(X_embed)

# GMM is the supplied final-label generator; only n_components is set to the
# requested HBC spatial-domain count.
labels = GaussianMixture(
    n_components=20, random_state=0, reg_covar=1e-3
).fit_predict(X_pca)
adata.obs["stHGC_cluster"] = labels.astype(str)

runtime_sec = time.perf_counter() - timer_start
print(
    f"BENCHMARK_END|stHGC|runtime_sec={runtime_sec:.6f}|clusters={len(np.unique(labels))}",
    flush=True,
)

os.makedirs(OUTPUT_DIR, exist_ok=True)
pd.DataFrame({"Barcode": adata.obs_names, "Cluster": labels}).to_csv(
    os.path.join(OUTPUT_DIR, "HBC_stHGC_labels.csv"), index=False
)

