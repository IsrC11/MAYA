# maya/interactive.py
import plotly.express as px
import plotly.graph_objects as go
import molplotly
import pandas as pd
import plotly.io as pio

def plot_interactive_scatter(df: pd.DataFrame, x: str, y: str,
                             smiles_col: str = "Canonical_Smiles",
                             color_col: str | None = None,
                             id_col: str | None = None,
                             title: str = "Interactive Chemical Space",
                             output_path: str | None = None):
    """
    Create interactive scatterplot with molecules shown on hover.
    """
    fig = px.scatter(
        df, x=x, y=y,
        color=color_col,
        hover_data=[smiles_col] + ([id_col] if id_col else []),
        title=title,
        opacity=0.7
    )

    # Añade moléculas interactivas con MolPlotly
    fig = molplotly.add_molecules(
        fig=fig,
        df=df,
        smiles_col=smiles_col,
        title_col=id_col if id_col else smiles_col,
        color_col=color_col,
        caption_cols=[smiles_col, x, y]
    )

    if output_path:
        pio.write_html(fig, output_path, auto_open=False)

    return fig


def plot_descriptor_radar(df: pd.DataFrame, id_col: str, compound_ids: list,
                           descriptor_cols: list | None = None,
                           output_path: str | None = None, show: bool = True,
                           title: str = "Perfil de Descriptores"):
    """
    Args:
        df: DataFrame con los descriptores ya calculados (self.data del analyzer).
        id_col: columna de identificador (p.ej. config.data['id_col']).
        compound_ids: lista de IDs a graficar (1 o varios, se sobreponen en el
            mismo radar para comparar).
        descriptor_cols: columnas a incluir en los ejes del radar. Default:
            MolWt, LogP, HBA, HBD, TPSA (los que ya calcula compute_descriptors).

    Nota de normalización: cada descriptor se reescala a [0,1] usando el
    mínimo/máximo de TODO el dataset (no solo de los compuestos seleccionados),
    para que la forma del radar sea comparable en el contexto de la librería
    completa -- un compuesto con LogP "alto" se ve alto respecto al resto del
    dataset, no respecto a sí mismo.
    """
    if descriptor_cols is None:
        descriptor_cols = ['MolWt', 'LogP', 'HBA', 'HBD', 'TPSA']

    missing = [c for c in descriptor_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Columnas de descriptor no encontradas en el dataset: {missing}")

    subset = df[df[id_col].isin(compound_ids)]
    if subset.empty:
        raise ValueError(f"Ningún compuesto con {id_col} en {compound_ids} fue encontrado")

    mins = df[descriptor_cols].min()
    maxs = df[descriptor_cols].max()
    ranges = (maxs - mins).replace(0, 1)  # evita división por cero si un descriptor es constante

    categories = descriptor_cols + [descriptor_cols[0]]  # cierra el polígono

    fig = go.Figure()
    for _, row in subset.iterrows():
        raw_values = [row[c] for c in descriptor_cols]
        norm_values = [(row[c] - mins[c]) / ranges[c] for c in descriptor_cols]
        norm_values += [norm_values[0]]
        hover_text = [f"{c}: {v:.2f}" for c, v in zip(descriptor_cols, raw_values)] + [f"{descriptor_cols[0]}: {raw_values[0]:.2f}"]

        fig.add_trace(go.Scatterpolar(
            r=norm_values,
            theta=categories,
            fill='toself',
            name=str(row[id_col]),
            text=hover_text,
            hoverinfo='text+name',
        ))

    fig.update_layout(
        polar=dict(radialaxis=dict(visible=True, range=[0, 1])),
        title=title,
        showlegend=True,
    )

    if output_path:
        pio.write_html(fig, output_path, auto_open=False)
    if show:
        fig.show()

    return fig
