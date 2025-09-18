# maya/visualization.py
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def plot_similarity_heatmap(sim_matrix, labels, output_path: str | None = None, show: bool = True, title: str | None = None):
    """Plot heatmap of similarity matrix."""
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(sim_matrix, xticklabels=False, yticklabels=False, cmap="magma", ax=ax)
    ax.set_title(title if title else "Tanimoto Similarity Heatmap")
    if output_path:
        fig.savefig(output_path, dpi=600, bbox_inches="tight")
    if show:
        plt.show()
    return fig

def plot_scatter(df: pd.DataFrame, x: str, y: str, hue: str | None=None, output_path: str | None=None, show: bool = True, title: str | None = None):
    """Scatter plot for chemical space visualization."""
    fig, ax =plt.subplots(figsize=(8, 6))
    sns.scatterplot(data=df, x=x, y=y, hue=hue, alpha=0.7, ax=ax)
    ax.set_title(title)
    if output_path:
        fig.savefig(output_path, dpi=600, bbox_inches="tight")
    if show:
        plt.show()
    return fig
