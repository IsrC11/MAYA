# maya/metrics.py
import numpy as np
import pandas as pd
from sklearn.manifold import trustworthiness
from sklearn.metrics import pairwise_distances
from scipy.stats import pearsonr

def calculate_similarity_correlation(original_space: pd.DataFrame, reduced_space: pd.DataFrame, metric='euclidean'):
    """
    Calculates the correlation of pairwise similarities between the original and reduced spaces.
    Args:
        original_space: pd.DataFrame of shape (n_compounds, n_features)
        reduced_space: pd.DataFrame (n_compounds, n_components)
        metric (str): The distance metric to compute pairwise distances.
    Returns:
        correlation (float): Correlation between original and reduced pairwise distances.
    """
    original_distances = pairwise_distances(original_space, metric=metric)
    reduced_distances = pairwise_distances(reduced_space, metric=metric)

    # Flatten upper triangular parts to avoid redundancy
    original_flat = original_distances[np.triu_indices_from(original_distances, k=1)]
    reduced_flat = reduced_distances[np.triu_indices_from(reduced_distances, k=1)]

    correlation, _ = pearsonr(original_flat, reduced_flat)
    return correlation


def evaluate_reduction(original_space: pd.DataFrame, reduced_space: pd.DataFrame, metric: str = 'euclidean', n_neighbors: int = 5):
    """
    Evaluate dimensionality reduction with trustworthiness and similarity correlation.
    """
    trust = trustworthiness(original_space, reduced_space, n_neighbors=n_neighbors, metric=metric)
    corr = calculate_similarity_correlation(original_space, reduced_space, metric=metric)
    return {"trustworthiness": trust, "correlation": corr}
