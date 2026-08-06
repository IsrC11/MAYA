# maya/metrics.py
import numpy as np
import pandas as pd
from sklearn.manifold import trustworthiness
from sklearn.metrics import pairwise_distances
from scipy.stats import spearmanr

def calculate_similarity_correlation(original_space: pd.DataFrame, reduced_space: pd.DataFrame, metric='euclidean'):
    """
    Calculates the correlation of pairwise similarities between the original and reduced spaces.
    Args:
        original_space: pd.DataFrame of shape (n_compounds, n_features)
        reduced_space: pd.DataFrame (n_compounds, n_components)
        metric (str): distancia para el espacio original. Usa 'jaccard' cuando
            original_space son fingerprints binarios (espacio estructural) y
            'euclidean' cuando son descriptores fisicoquímicos ya escalados
            (espacio de propiedades). El espacio reducido (coordenadas 2D/3D
            continuas) siempre se evalúa con distancia euclidiana, que sí es
            apropiada ahí sin importar el espacio de origen.
    Returns:
        correlation (float): Correlación (Spearman) entre las distancias por
            pares del espacio original y del espacio reducido.
    """
    original_distances = pairwise_distances(original_space, metric=metric)
    reduced_distances = pairwise_distances(reduced_space, metric='euclidean')

    # Flatten upper triangular parts to avoid redundancy
    original_flat = original_distances[np.triu_indices_from(original_distances, k=1)]
    reduced_flat = reduced_distances[np.triu_indices_from(reduced_distances, k=1)]

    # CAMBIO: Pearson -> Spearman. Lo que interesa al evaluar una reducción
    # dimensional es si se preserva el ORDEN relativo de las similitudes
    # (¿los compuestos más parecidos en el espacio original siguen siendo los
    # más cercanos en 2D?), no si existe una relación estrictamente LINEAL
    # entre ambas distancias -- t-SNE y UMAP en particular no buscan preservar
    # distancias absolutas, solo vecindades/orden, así que Pearson las
    # penalizaría incorrectamente aunque la reducción sea buena.
    correlation, _ = spearmanr(original_flat, reduced_flat)
    return correlation


def evaluate_reduction(original_space: pd.DataFrame, reduced_space: pd.DataFrame, metric: str = 'euclidean', n_neighbors: int = 5):
    """
    Evaluate dimensionality reduction with trustworthiness and similarity correlation.

    Args:
        metric: métrica del ESPACIO ORIGINAL -- pasa 'jaccard' si evalúas una
            reducción hecha sobre fingerprints binarios (espacio estructural),
            o deja 'euclidean' (default) para descriptores escalados (espacio
            de propiedades). # CAMBIO: antes siempre se usaba 'euclidean' sin
            importar qué tipo de datos traía original_space, lo cual es
            estadísticamente incorrecto quando original_space son bits.
    """
    trust = trustworthiness(original_space, reduced_space, n_neighbors=n_neighbors, metric=metric)
    corr = calculate_similarity_correlation(original_space, reduced_space, metric=metric)
    return {"trustworthiness": trust, "correlation": corr}
