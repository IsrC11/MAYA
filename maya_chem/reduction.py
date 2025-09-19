# maya/reduction.py
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap

def apply_pca(fps, n_components: int = 2) -> pd.DataFrame:
    """Apply PCA dimensionality reduction."""
    pca = PCA(n_components=n_components)
    return pca.fit_transform(fps)

def apply_tsne(fps, n_components: int = 2, perplexity: int = 30) -> pd.DataFrame:
    """Apply t-SNE dimensionality reduction."""
    tsne = TSNE(n_components=n_components, perplexity=perplexity, random_state=42)
    return tsne.fit_transform(fps)

def apply_umap (fps, n_components=2):
    reducer = umap.UMAP(n_components=n_components, random_state=42)
    return reducer.fit_transform(fps)
