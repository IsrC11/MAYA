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
    else:
        plt.close(fig)
    return fig

def plot_scatter(df: pd.DataFrame, x: str, y: str, hue: str | None=None, palette= None, output_path: str | None=None, show: bool = True, title: str | None = None):
    """Scatter plot for chemical space visualization."""
    fig, ax =plt.subplots(figsize=(8, 6))
    sns.scatterplot(data=df, x=x, y=y, hue=hue, palette=palette if hue else None, alpha=0.7, ax=ax)
    ax.set_title(title)
    if output_path:
        fig.savefig(output_path, dpi=600, bbox_inches="tight")
    plt.close(fig)
    return fig


def plot_pca_biplot(coords_df: pd.DataFrame, loadings, feature_names: list,
                     x_col: str, y_col: str, output_path: str | None = None,
                     show: bool = True, title: str | None = "PCA Biplot",
                     arrow_scale: float = 1.0):
    """
    Args:
        coords_df: DataFrame con las coordenadas reducidas (debe incluir x_col, y_col).
        loadings: array (n_components, n_features) -- analyzer.pca_loadings.
        feature_names: nombres de los descriptores en el mismo orden que loadings
            (analyzer.pca_feature_names).
        x_col, y_col: nombres de columna en coords_df a graficar (p.ej. 'PCA1', 'PCA2').
        arrow_scale: factor para alargar/acortar visualmente las flechas
            (no cambia la dirección/proporción relativa entre descriptores).
    """
    if x_col not in coords_df.columns or y_col not in coords_df.columns:
        raise ValueError(f"'{x_col}' y/o '{y_col}' no están en coords_df")

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(coords_df[x_col], coords_df[y_col], alpha=0.4, s=20, color="steelblue")

    max_coord = max(coords_df[x_col].abs().max(), coords_df[y_col].abs().max())
    if max_coord == 0 or pd.isna(max_coord):
        max_coord = 1.0

    x_idx = 0  # loadings[0, :] corresponde al componente graficado en x_col
    y_idx = 1  # loadings[1, :] corresponde al componente graficado en y_col

    for i, feat in enumerate(feature_names):
        vx = loadings[x_idx, i] * max_coord * arrow_scale
        vy = loadings[y_idx, i] * max_coord * arrow_scale
        ax.annotate(
            "", xy=(vx, vy), xytext=(0, 0),
            arrowprops=dict(arrowstyle="->", color="crimson", lw=1.5),
        )
        ax.text(vx * 1.15, vy * 1.15, feat, color="crimson", ha="center",
                 va="center", fontsize=10, fontweight="bold")

    ax.axhline(0, color="grey", lw=0.5, linestyle="--")
    ax.axvline(0, color="grey", lw=0.5, linestyle="--")
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title(title)

    if output_path:
        fig.savefig(output_path, dpi=600, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig
