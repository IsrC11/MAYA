# maya/reduction.py
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap

def apply_pca(fps, n_components: int = 2) -> pd.DataFrame:
    """Apply PCA dimensionality reduction."""
    pca = PCA(n_components=n_components)
    coords = pca.fit_transform(fps)
    return pd.DataFrame(coords, columns=[f"PCA{i+1}" for i in range(n_components)])

def apply_tsne(fps, n_components: int = 2, perplexity: int = 30) -> pd.DataFrame:
    """Apply t-SNE dimensionality reduction."""
    tsne = TSNE(n_components=n_components, perplexity=perplexity, random_state=42)
    coords = tsne.fit_transform(fps)
    return pd.DataFrame(coords, columns=[f"Dimension{i+1}" for i in range(n_components)])

def apply_umap (fps, n_components=2):
    reducer = umap.UMAP(n_components=n_components, random_state=42)
    coords = reducer.fit_transforms(fps)
    return pd.DataFrame(coords, columns=[f'Dimension{i+1}' for i in range(n_components)])
