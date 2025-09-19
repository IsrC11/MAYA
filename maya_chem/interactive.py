# maya/interactive.py
import plotly.express as px
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
