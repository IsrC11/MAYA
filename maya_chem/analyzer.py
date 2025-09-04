# maya/analyzer.py
import pandas as pd
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
        self.fps = [descriptors.compute_morgan_fingerprint(smi) for smi in self.data[smiles_col]]
        return self.fps, self.data

    def compute_similarity(self, n_jobs: int = -1):
        self.sim_matrix = similarity.compute_similarity_matrix(self.fps, n_jobs)
        return self.sim_matrix

    def reduce_dimensions(self, method: str = 'pca', n_components: int = 2):
        if method.lower() == 'pca':
            coords = reduction.apply_pca(sim_matrix, n_components=n_components)
        elif method.lower() == 'tsne':
            coords = reduction. apply_tsne(sim_matrix, n_components=n_components)    
        elif method.lower() == 'umap':
            coords = reduction.apply_umap(sim_matrix, n_components=n_componets)
        else:
            raise ValueError(f'Unknown dimentionallity reduction method: {method}')
        self.data=pd.concat([self.data, coords], axis=1)
        return coords
        
    def visualize(self, show: bool = True, save_prefix: str | None = None):
        fig1 =visualization.plot_similarity_heatmap(self.sim_matrix, label=self.data[self.config.data[id_col]], output_path=f'{save_prefix}_heatmap.png' if save_prefix else None, show=how)

        fig2 = visualization.plot_scatter(self.data, x='Dimention1', y='Dimention2', hue='MolWt' if 'MolWt' in self.data.columns else None, output_path=f'{save_prefix}_scatter.png' if save_prefix else None, show=show)
        return self.visualize(show=True)

    
    def run(self):
        self.load_data()
        self.compute_descriptors()
        self.compute_similarity()
        self.reduce_dimentions(method=self.config.data.get('reduction', 'pca'))
        return self.visualize(show=True)
