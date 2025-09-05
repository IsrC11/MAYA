# maya/similarity.py
import numpy as np
from rdkit import DataStructs
from joblib import Parallel, delayed
from typing import List

def compute_similarity_matrix(fps: List, n_jobs: int = -1) -> np.ndarray:
    """Compute Tanimoto similarity matrix in parallel."""
    def compute_sim(fp, fps_list):
        return [round(DataStructs.TanimotoSimilarity(fp, f), 3) for f in fps_list]

    sim_matrix = np.array(Parallel(n_jobs=n_jobs)(delayed(compute_sim)(fp, fps) for fp in fps))
    np.fill_diagonal(sim_matrix, 1.0)
    return sim_matrix
