# maya/analyzer.py
import pandas as pd
import os
import numpy as np
from rdkit import DataStructs
from .config import MayaConfig
from .utils import load_data
from . import utils, descriptors, similarity, reduction, visualization, interactive

class MayaAnalyzer:
    """Main class to run MAYA pipeline."""

    def __init__(self, config: MayaConfig):
        self.config = config
        self.data: pd.DataFrame | None = None
        self.fps = None
        self.sim_matrix: np.ndarray | None = None

    def load_data(self):
        self.data = utils.load_data(self.config.data_path, id_col=self.config.data['id_col'], smiles_col=self.config.data['smiles_col'])
        return self.data

    def compute_descriptors(self, fp_type: str):
        smiles_col = self.config.data['smiles_col']

        props_series = self.data[smiles_col].apply(descriptors.compute_physicochemical_properties)
        props_df = pd.DataFrame(list(props_series))
        self.data = pd.concat([self.data, props_df], axis=1)
        
        self.fps = [descriptors.compute_fingerprint(smi, method=fp_type) for smi in self.data[smiles_col]]

    def compute_similarity(self, n_jobs: int = -1):
        self.sim_matrix = similarity.compute_similarity_matrix(self.fps, n_jobs)
        return self.sim_matrix

    def reduce_dimensions(self, method: str = 'pca', n_components: int = 2):
        
        x = np.array([np.frombuffer(fp.ToBitString().encode('utf-8'), dtype='S1') for fp in self.fps])
        x = (x.view(np.uint8) - ord('0')).reshape(len(self.fps), -1)
        method_lower = method.lower()
        if method_lower == 'pca':
            coords = reduction.apply_pca(x, n_components=n_components)
            prefix = 'PCA'
        elif method_lower == 'tsne':
            coords = reduction. apply_tsne(x, n_components=n_components) 
            prefix = 'Dim'
        elif method_lower == 'umap':
            coords = reduction.apply_umap(x, n_components=n_components)
            prefix = 'Dim'
        else:
            raise ValueError(f'Unknown dimentionallity reduction method: {method}')

        coords = pd.DataFrame(coords, index=self.data.index, columns=[f'{prefix}{i+1}' for i in range(coords.shape[1])])
        self.data=pd.concat([self.data, coords], axis=1)
        return coords
        
    def visualize(self, show: bool = True, save_prefix: str | None = None, title:str | None = None, heatmap_title: str = 'Tanimoto Heatmap', interactive_mode: bool = False, port: int = 8050):
        coords_cols = [col for col in self.data.columns if col.startswith('PCA') or col.startswith('Dim')]
        import plotly.express as px
        from molplotly import add_molecules
        from google.colab.output import serve_kernel_port_as_iframe
        from .visualization import plot_similarity_heatmap
        from dash import Dash, dcc, html, Input, Output
        
        if len(coords_cols) < 2:
            raise ValueError('At least two reduction columns (PC1/PC2 or Dim1/Dim2) were not found')

        x_col, y_col = coords_cols[:2]

        hover_cols = [self.config.data['smiles_col'], x_col, y_col]

        palette = self.config.viz.get('palette', None)
        color_col = self.config.viz.get('color_by',None)
        #if color_col not in self.data.columns:
            #color_col = None

        if save_prefix:
            from . import visualization
            color_by = self.config.viz.get('color_by', None)
            palette = self.config.viz.get('palette', None)
            visualization.plot_scatter(self.data, x=x_col, y=y_col, hue=color_by if color_by in self.data.columns else None, palette = palette, output_path=f'{save_prefix}_scatter.png', show=False, title=title)

        if interactive_mode:
            import molplotly
            from jupyter_dash import JupyterDash
            from plotly import graph_objects as go
 
            try:
                fig= px.scatter(self.data, x=x_col, y=y_col, color=color_col, title=title, width=900, height=700, color_continuous_scale=palette if color_col else None)
                fig = molplotly.add_molecules(fig=fig, df=self.data, smiles_col=self.config.data['smiles_col'], title_col=self.config.data['id_col'], color_col=color_col)
                serve_kernel_port_as_iframe('localhost')
                fig.run(port=port)
                return fig
                
            except Exception as e:
                raise RuntimeError(f'Error to generate interactive graph with Dash: {e}')
      
        return None

    
    def run(self):
        fingerprints = self.config.analysis['fingerprint']
        reductions = self.config.analysis['reduction_method']

        if isinstance(fingerprints, str):
            fingerprints = [fingerprints]
        if isinstance(reductions, str):
            reductions = [reductions]
        
        os.makedirs(self.config.viz['output_dir'], exist_ok=True)
        
        original_data = self.load_data()
        results = []
        port=8050

        for fp in fingerprints:
            self.data = original_data.copy()
            self.compute_descriptors(fp_type=fp)
            self.compute_similarity()
            heatmap_title = f'Tanimoto Heatmap - {fp.upper()}'
            heatmap_path = f"{self.config.viz['output_dir']}/{fp}_heatmap.png"
            heatmap_figure = visualization.plot_similarity_heatmap(self.sim_matrix, labels=False, output_path = heatmap_path, show=True, title=heatmap_title)
            results.append((fp, 'heatmap', heatmap_figure))
            
            for red in reductions:
                self.data = original_data.copy()
                self.compute_descriptors(fp_type=fp)
                reduced = self.reduce_dimensions(method=red)
                save_prefix = f'{fp}_{red}'
                heatmap_title = f'Tanimoto Heatmap - {fp.upper()}'
                scatter_title = f'{fp.upper()} + {red.upper()}'
                figs = self.visualize(save_prefix=save_prefix, show=False, title=scatter_title, heatmap_title=heatmap_title, interactive_mode=True, port=port)
                results.append((fp, red, figs))
                port+=3
        
        return results
