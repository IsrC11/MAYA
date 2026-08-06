# maya/similarity.py
import numpy as np
from rdkit import DataStructs
from joblib import Parallel, delayed
from typing import List

def _row_similarities(i: int, fps: List):
    return DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:])

def compute_similarity_matrix(fps: List, n_jobs: int = -1) -> np.ndarray:
    """Compute Tanimoto similarity matrix in parallel."""
    n = len(fps)
    sim_matrix = np.ones((n, n), dtype=np.float64)

    if n > 1:
        rows = Parallel(n_jobs=n_jobs)(
            delayed(_row_similarities)(i, fps) for i in range(n - 1)
        )
        for i, sims in enumerate(rows):
            sim_matrix[i, i + 1:] = sims
            sim_matrix[i + 1:, i] = sims

    return sim_matrix
