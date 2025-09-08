# maya/visualization.py
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def plot_similarity_heatmap(sim_matrix, labels, output_path: str | None = None, show: bool = True):
    """Plot heatmap of similarity matrix."""
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(sim_matrix, xticklabels=False, yticklabels=False, cmap="viridis", ax=ax)
    ax.set_title("Tanimoto Similarity Heatmap")
    if output_path:
        fig.savefig(output_path, dpi=600, bbox_inches="tight")
    if show:
        plt.show()
    return fig

def plot_scatter(df: pd.DataFrame, x: str, y: str, hue: str | None=None, output_path: str | None=None, show: bool = True):
    """Scatter plot for chemical space visualization."""
    fig, ax =plt.subplots(figsize=(8, 6))
    sns.scatterplot(data=df, x=x, y=y, hue=hue, alpha=0.7, ax=ax)
    ax.set_title('Chemical Space')
    if output_path:
        fig.savefig(output_path, dpi=600, bbox_inches="tight")
    if show:
        plt.show()
    return fig
