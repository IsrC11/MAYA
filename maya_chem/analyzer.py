# maya/analyzer.py
import pandas as pd
import os
import numpy as np
from rdkit import DataStructs
from .config import MayaConfig
from .utils import load_data
from . import utils, descriptors, similarity, reduction, visualization

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
        
        if method.lower() == 'pca':
            coords = reduction.apply_pca(x, n_components=n_components)
        elif method.lower() == 'tsne':
            coords = reduction. apply_tsne(x, n_components=n_components)    
        elif method.lower() == 'umap':
            coords = reduction.apply_umap(x, n_components=n_components)
        else:
            raise ValueError(f'Unknown dimentionallity reduction method: {method}')

        coords = pd.DataFrame(coords, index=self.data.index, columns=[f'{method.upper()}_{i+1}' for i in range(coords.shape[1])])
        self.data=pd.concat([self.data, coords], axis=1)
        return coords
        
    def visualize(self, show: bool = True, save_prefix: str | None = None, title:str = 'Chemical Space', heatmap_title: str = 'Tanimoto Heatmap'):
        coords_cold = [col for col in self.data.columns if col,startswitth('PCA') or col.startswith('Dim')]

        if len(coords_cols) < 2:
            raise ValueError('At least two reduction columns (PC1/PC2 or Dim1/Dim2) were not found')

        x_col, y_col = coord_cols[:2]
        
        fig1 =visualization.plot_similarity_heatmap(self.sim_matrix, labels=self.data[self.config.data['id_col']], output_path=f'{save_prefix}_heatmap.png' if save_prefix else None, show=show, title =heatmap_title)

        fig2 = visualization.plot_scatter(self.data, x='PC1', y='PC2', hue='MolWt' if 'MolWt' in self.data.columns else None, output_path=f'{save_prefix}_scatter.png' if save_prefix else None, show=show, title=title)
        return fig1, fig2

    
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

        for fp in fingerprints:
            self.data = original_data.copy()
            self.compute_descriptors(fp_type=fp)
            self.compute_similarity()
            heatmap_title = f'Tanimoto Heatmap - {fp.upper()}'
            heatmap_path = f'{self.config.viz['output_dir']}/{fp}_heatmap.png'
            heatmap_figure = visualization.plot_similarity_heatmap(self.sim_matrix, labels=False, output_path = heatmap_path, show=True, title=heatmap_title)
            results.append((fp, 'heatmap', heatmap_figure))
            
            for red in reductions:
                reduced = self.reduce_dimensions(method=red)
                save_prefix = f'{fp}_{red}'
                heatmap_title = f'Tanimoto Heatmap - {fp.upper()}'
                scatter_title = f'{fp.upper()} + {red.upper()}'
                figs = self.visualize(save_prefix=save_prefix, show=True, title=scatter_title, heatmap_title=heatmap_title)
                results.append((fp, red, figs))
        
        return results
