# maya/reduction.py
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
import umap


def apply_pca(fps, n_components: int = 2):
    """Apply PCA dimensionality reduction. Úsala solo sobre datos continuos/
    escalados (p.ej. descriptores fisicoquímicos), NO sobre fingerprints
    binarios crudos -- PCA asume distancia euclidiana, y la distancia natural
    entre fingerprints es 1 - Tanimoto. Para fingerprints usa apply_structure_pcoa().
    """
    pca = PCA(n_components=n_components)
    coords = pca.fit_transform(fps)
    return coords, pca.explained_variance_ratio_, pca.components_


def apply_structure_pcoa(sim_matrix: np.ndarray, n_components: int = 2):
    """Principal Coordinate Analysis (classical MDS) sobre una matriz de
    similitud de Tanimoto ya calculada.

    Returns:
        coords: array (n_samples, n_components)
        explained_variance_ratio: array (n_components,) -- análogo al de PCA,
            basado en los eigenvalores de la matriz de Gram doble-centrada.
    """
    sim_matrix = np.asarray(sim_matrix, dtype=np.float64)
    dist = 1.0 - sim_matrix
    n = dist.shape[0]

    
    dist_sq = dist ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ dist_sq @ J

    eigvals, eigvecs = np.linalg.eigh(B)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    eigvals_top = np.clip(eigvals[:n_components], a_min=0, a_max=None)
    coords = eigvecs[:, :n_components] * np.sqrt(eigvals_top)

    total_positive_variance = np.sum(np.clip(eigvals, a_min=0, a_max=None))
    if total_positive_variance > 0:
        explained_variance_ratio = eigvals_top / total_positive_variance
    else:
        explained_variance_ratio = np.zeros(n_components)

    return coords, explained_variance_ratio


def apply_tsne(fps, n_components: int = 2, perplexity: int = 30, metric: str = 'euclidean'):
    """Apply t-SNE dimensionality reduction.

    Args:
        metric: 'jaccard' para fingerprints binarios (estructura), 'euclidean'
            (default) para descriptores continuos ya escalados (propiedades).
            # CAMBIO: antes 'metric' no existía como parámetro y TSNE siempre
            # usaba 'euclidean' por default, incluso cuando 'fps' eran bits
            # binarios -- estadísticamente incorrecto para ese caso.
    """
    tsne = TSNE(n_components=n_components, perplexity=perplexity, metric=metric, random_state=42)
    return tsne.fit_transform(fps)


def apply_umap(fps, n_components: int = 2, metric: str = 'euclidean'):
    """Apply UMAP dimensionality reduction.

    Args:
        metric: 'jaccard' para fingerprints binarios (estructura), 'euclidean'
            (default) para descriptores continuos ya escalados (propiedades).
            # CAMBIO: mismo caso que en apply_tsne -- antes siempre euclidiana.
    """
    reducer = umap.UMAP(n_components=n_components, metric=metric, random_state=42)
    return reducer.fit_transform(fps)


def scale_descriptors(props_df) -> np.ndarray:
    """Estandariza (z-score) descriptores fisicoquímicos antes de PCA/t-SNE/UMAP.

    # CAMBIO: función nueva. MolWt (~150-600), TPSA (~0-150), LogP (~-3-6) y
    # HBA/HBD (enteros 0-10) tienen rangos muy distintos. Sin escalar, PCA
    # maximiza varianza en unidades crudas y MolWt domina los primeros
    # componentes solo por su rango numérico, no por ser más informativo
    # químicamente -- es el error estadístico clásico de "features sin escalar".
    """
    return StandardScaler().fit_transform(props_df)
