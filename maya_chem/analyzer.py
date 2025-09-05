# maya/analyzer.py
import pandas as pd
import numpy as np
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

    def compute_descriptors(self):
        smiles_col = self.config.data['smiles_col']
        self.data["Properties"] = self.data[smiles_col].apply(descriptors.compute_physicochemical_properties)
        desc_type = self.config.analysis.get('descriptor', 'morgan')
        self.fps = [descriptors.compute_fingerprint(smi, method=desc_type) for smi in self.data[smiles_col]]

    def compute_similarity(self, n_jobs: int = -1):
        self.sim_matrix = similarity.compute_similarity_matrix(self.fps, n_jobs)
        return self.sim_matrix

    def reduce_dimensions(self, method: str = None, n_components: int = 2):
        
        method = method or self.config.analysis.get('reduction_method', 'pca')
        
        x = np.array([np.frombuffer(fp.ToBitString().encode('utf-8'), dtype='S1') for fp in self.fps])
        x = (x.view(np.uint8) - ord('0')).reshape(len(self.fps), -1)
        
        if method.lower() == 'pca':
            coords = reduction.apply_pca(x, n_components=n_components)
        elif method.lower() == 'tsne':
            coords = reduction. apply_tsne(x, n_components=n_components)    
        elif method.lower() == 'umap':
            coords = reduction.apply_umap(x, n_components=n_componets)
        else:
            raise ValueError(f'Unknown dimentionallity reduction method: {method}')

        coords.columns = [f'PC{i+1}' for i in range(coords.shape[1])]
        self.data=pd.concat([self.data, coords], axis=1)
        return coords
        
    def visualize(self, show: bool = True, save_prefix: str | None = None):
        fig1 =visualization.plot_similarity_heatmap(self.sim_matrix, labels=self.data[self.config.data['id_col']], output_path=f'{save_prefix}_heatmap.png' if save_prefix else None, show=show)

        fig2 = visualization.plot_scatter(self.data, x='PC1', y='PC2', hue='MolWt' if 'MolWt' in self.data.columns else None, output_path=f'{save_prefix}_scatter.png' if save_prefix else None, show=show)
        return fig1, fig2

    
    def run(self):
        self.load_data()
        self.compute_descriptors()
        self.compute_similarity()
        self.reduce_dimensions(method=self.config.data.get('reduction', 'pca'))
        return self.visualize(show=True)
