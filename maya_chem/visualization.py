# maya/visualization.py
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull

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


# CAMBIO: función nueva -- biplot de PCA. Muestra los compuestos en el espacio
# reducido (igual que plot_scatter) MÁS un vector por cada descriptor original
# indicando qué tanto y en qué dirección contribuye a los ejes mostrados. Es la
# forma estándar de interpretar qué significa "moverse" en un PC: por ejemplo,
# si el vector de MolWt apunta fuerte hacia +PC1, los compuestos con PC1 alto
# tienden a ser de mayor peso molecular.
#
# Solo es válido para reducciones con space='properties' (PCA real sobre
# descriptores escalados) -- necesita 'loadings' (analyzer.pca_loadings), que
# no existen para space='structure' (PCoA parte de una matriz de distancias,
# no de features originales, así que no hay "carga" de un descriptor sobre un eje).
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


# CAMBIO: función nueva -- contornos de densidad (KDE) sobre el espacio químico
# reducido. Complementa al scatter: con muchos compuestos el scatter se satura
# (overplotting) y es difícil distinguir "zona densamente poblada" de "puñado
# de outliers sueltos" solo por la nube de puntos. El KDE responde eso
# directamente marcando isolíneas de densidad estimada.
def plot_density_contours(df: pd.DataFrame, x: str, y: str,
                           output_path: str | None = None, show: bool = True,
                           title: str | None = "Densidad del Espacio Químico",
                           overlay_points: bool = True):
    """
    Args:
        overlay_points: si True, dibuja los puntos individuales debajo de los
            contornos (recomendado -- el KDE solo es un resumen, no reemplaza
            ver los compuestos reales, incluidos los outliers que caen fuera
            de cualquier contorno).
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    if overlay_points:
        ax.scatter(df[x], df[y], s=10, alpha=0.25, color="grey", zorder=1)
    sns.kdeplot(data=df, x=x, y=y, fill=True, cmap="magma", alpha=0.6,
                thresh=0.05, levels=8, ax=ax, zorder=2)
    ax.set_title(title)
    if output_path:
        fig.savefig(output_path, dpi=600, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig


# CAMBIO: función nueva -- scatter coloreado por clúster (analyzer.cluster_compounds)
# con una envolvente convexa (convex hull) dibujada alrededor de cada clúster,
# para que la agrupación sea visualmente evidente incluso con muchos puntos
# superpuestos. Los clusters se calculan en el espacio original (estructura o
# propiedades, ver cluster_compounds); esta función solo los DIBUJA sobre las
# coordenadas 2D ya reducidas -- la envolvente es una ayuda visual, no una
# frontera estadística real (dos clusters pueden solaparse en 2D aunque estén
# bien separados en el espacio original, precisamente por lo que se explicó
# arriba sobre por qué no se clusteriza directo en 2D).
def plot_cluster_hulls(df: pd.DataFrame, x: str, y: str, cluster_col: str = "cluster",
                        output_path: str | None = None, show: bool = True,
                        title: str | None = "Clusters en el Espacio Químico",
                        palette: str = "tab10"):
    fig, ax = plt.subplots(figsize=(8, 6))
    clusters = sorted(df[cluster_col].unique())
    colors = sns.color_palette(palette, n_colors=len(clusters))

    for cluster_id, color in zip(clusters, colors):
        subset = df[df[cluster_col] == cluster_id]
        ax.scatter(subset[x], subset[y], s=25, alpha=0.7, color=color,
                    label=f"Cluster {cluster_id}" if cluster_id != -1 else "Ruido / outliers")

        # La envolvente convexa necesita al menos 3 puntos no colineales;
        # clusters muy pequeños se dejan sin hull en vez de fallar.
        if len(subset) >= 3:
            points = subset[[x, y]].to_numpy()
            try:
                hull = ConvexHull(points)
                hull_pts = points[hull.vertices]
                hull_pts = np.vstack([hull_pts, hull_pts[0]])  # cerrar el polígono
                ax.plot(hull_pts[:, 0], hull_pts[:, 1], color=color, lw=1.5, alpha=0.8)
                ax.fill(hull_pts[:, 0], hull_pts[:, 1], color=color, alpha=0.08)
            except Exception:
                pass  # puntos colineales u otro caso degenerado: se omite el hull

    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_title(title)
    ax.legend(loc="best", fontsize=8)

    if output_path:
        fig.savefig(output_path, dpi=600, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig
