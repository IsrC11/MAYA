# maya/analyzer.py
import pandas as pd
import os
import numpy as np
from rdkit import DataStructs
from .curation import MayaCuration, clean_dataset
from .config import MayaConfig
from .utils import load_data
from . import utils, descriptors, similarity, reduction, visualization, interactive, metrics
from .metrics import evaluate_reduction

class MayaAnalyzer:
    """Main class to run MAYA pipeline."""

    def __init__(self, config: MayaConfig):
        self.config = config
        self.data: pd.DataFrame | None = None
        self.fps = None
        self.sim_matrix: np.ndarray | None = None
        

    def load_data(self):
        self.data = utils.load_data(self.config.data_path, id_col=self.config.data['id_col'], smiles_col=self.config.data['smiles_col'], canonicalize=False)
        self.data = clean_dataset(self.data, smiles_col=self.config.data['smiles_col'], curation_config=self.config.curation)
        return self.data

    def compute_descriptors(self, fp_type: str):
        smiles_col = self.config.data['smiles_col']

        props_series = self.data[smiles_col].apply(descriptors.compute_physicochemical_properties)
        props_df = pd.DataFrame(list(props_series))
        self.data = pd.concat([self.data, props_df], axis=1)

        if fp_type.lower() == 'map4':
            smiles_list = self.data[smiles_col].tolist()
            self.fps = descriptors.compute_map4_fingerprint(smiles_list)
        else:
            self.fps = [descriptors.compute_fingerprint(smi, method=fp_type) for smi in self.data[smiles_col]]

    def compute_similarity(self, n_jobs: int = -1):
        self.sim_matrix = similarity.compute_similarity_matrix(self.fps, n_jobs)
        return self.sim_matrix

    def reduce_dimensions(self, method: str = 'pca', n_components: int = 2, space: str = 'structure'):
        """
        Args:
            space: 'structure' (default) reduce el espacio de fingerprints
                usando la métrica de Tanimoto/Jaccard, correcta para bits
                binarios. 'properties' reduce el espacio de descriptores
                fisicoquímicos (MolWt, LogP, HBA, HBD, TPSA...) escalados con
                z-score, correcto para PCA/t-SNE/UMAP euclidianos.
                # CAMBIO: antes solo existía un camino -- PCA/t-SNE/UMAP con
                # distancia euclidiana SIEMPRE sobre el array binario del
                # fingerprint. Esto mezclaba dos preguntas distintas
                # ("¿qué tan diversa es la librería estructuralmente?" vs.
                # "¿cómo se distribuyen las propiedades tipo Lipinski?") en un
                # solo método, y encima con una métrica (euclidiana) que no es
                # la apropiada para bits binarios.
        """
        method_lower = method.lower()
        explained_variance = None
        coords = None
        prefix = ''

        if space == 'structure':
            x = np.array([np.frombuffer(fp.ToBitString().encode('utf-8'), dtype='S1') for fp in self.fps])
            x = (x.view(np.uint8) - ord('0')).reshape(len(self.fps), -1)
            eval_space = x
            eval_metric = 'jaccard'

            if method_lower == 'pca':
                if self.sim_matrix is None:
                    self.compute_similarity()
                coords, explained_variance = reduction.apply_structure_pcoa(self.sim_matrix, n_components=n_components)
                self.explained_variance = explained_variance
                prefix = 'PCoA'
            elif method_lower == 'tsne':
                coords = reduction.apply_tsne(x, n_components=n_components, metric='jaccard')
                prefix = 'Dim'
                if hasattr(self, 'explained_variance'):
                    del self.explained_variance
            elif method_lower == 'umap':
                coords = reduction.apply_umap(x, n_components=n_components, metric='jaccard')
                prefix = 'Dim'
                if hasattr(self, 'explained_variance'):
                    del self.explained_variance
            else:
                raise ValueError(f'Unknown dimentionallity reduction method: {method}')

        elif space == 'properties':
            prop_cols = [c for c in self.config.analysis['properties'] if c in self.data.columns]
            if len(prop_cols) < 2:
                raise ValueError(
                    f"Se necesitan al menos 2 descriptores calculados para space='properties'. "
                    f"Encontrados: {prop_cols}. ¿Corriste compute_descriptors() antes?"
                )
            x_scaled = reduction.scale_descriptors(self.data[prop_cols])
            eval_space = x_scaled
            eval_metric = 'euclidean'

            if method_lower == 'pca':
                coords, explained_variance, components = reduction.apply_pca(x_scaled, n_components=n_components)
                self.explained_variance = explained_variance
                self.pca_loadings = components
                self.pca_feature_names = prop_cols
                prefix = 'PCA'
            elif method_lower == 'tsne':
                coords = reduction.apply_tsne(x_scaled, n_components=n_components, metric='euclidean')
                prefix = 'Dim'
                if hasattr(self, 'explained_variance'):
                    del self.explained_variance
            elif method_lower == 'umap':
                coords = reduction.apply_umap(x_scaled, n_components=n_components, metric='euclidean')
                prefix = 'Dim'
                if hasattr(self, 'explained_variance'):
                    del self.explained_variance
            else:
                raise ValueError(f'Unknown dimentionallity reduction method: {method}')
        else:
            raise ValueError(f"Unknown space: '{space}'. Usa 'structure' o 'properties'.")

        coords = pd.DataFrame(coords, index=self.data.index, columns=[f'{prefix}{i+1}' for i in range(coords.shape[1])])
        self.data = pd.concat([self.data, coords], axis=1)

        results_eval = evaluate_reduction(eval_space, coords, metric=eval_metric)
        trust = results_eval['trustworthiness']
        coor = results_eval['correlation']

        print(f'Trustworthiness ({method}, space={space}): {trust:.3f}')
        print(f'Correlation ({method}, space={space}): {coor:.3f}')

        return coords, results_eval, trust, coor, explained_variance

    def visualize(self, show: bool = True, save_prefix: str | None = None, title: str = 'Chemical Space', heatmap_title: str = 'Tanimoto Heatmap', interactive_mode: bool = False, port: int = 8060, color_by: str | None = None):
        
        coords_cols = [col for col in self.data.columns if col.startswith('PCA') or col.startswith('PCoA') or col.startswith('Dim')]
        import plotly.express as px
        from molplotly import add_molecules
        try:
            from google.colab.output import serve_kernel_port_as_iframe
            _in_colab = True
        except ImportError:
            serve_kernel_port_as_iframe = None
            _in_colab = False

        from .visualization import plot_similarity_heatmap
        from dash import Dash, dcc, html, Input, Output

        if len(coords_cols) < 2:
            raise ValueError('At least two reduction columns (PC1/PC2 or Dim1/Dim2) were not found')

        x_col, y_col = coords_cols[:2]

        hover_cols = [self.config.data['smiles_col'], x_col, y_col]

        palette = self.config.viz.get('palette', None)
        color_col = self.config.viz.get('color_by', None)

        color_col = color_by if color_by else self.config.viz.get('color_by', None)
        if color_col not in self.data.columns:
            print(f'Warning: "{color_col}" not exist in dataset')
            #color_col = None

        if save_prefix:
            from . import visualization
            visualization.plot_scatter(self.data, x=x_col, y=y_col, hue=color_col, palette=palette, output_path=f'{save_prefix}_scatter.png', show=False, title=title)

        if interactive_mode:
            import molplotly
            from plotly import graph_objects as go

            try:
                fig = px.scatter(self.data, x=x_col, y=y_col, color=color_col, title=title, width=900, height=700, color_continuous_scale=palette)

                if hasattr(self, 'explained_variance') and self.explained_variance is not None:
                    x_label = f'PC1({self.explained_variance[0]*100:.2f}%)'
                    y_label = f'PC2({self.explained_variance[1]*100:.2f}%)'
                    fig.update_layout(xaxis_title=x_label, yaxis_title=y_label)
                else:
                    x_label = 'Dim1'
                    y_label = 'Dim2'

                fig.update_layout(plot_bgcolor='white', paper_bgcolor='white', xaxis=dict(showgrid=False, zeroline=False, mirror=True), yaxis=dict(showgrid=False, zeroline=False, mirror=True))
                fig = molplotly.add_molecules(fig=fig, df=self.data, smiles_col=self.config.data['smiles_col'], title_col=self.config.data['id_col'], color_col=color_col, caption_cols=self.config.data['activities'])

                if _in_colab:
                    serve_kernel_port_as_iframe('localhost')
                fig.run(port=port)
                return fig

            except Exception as e:
                raise RuntimeError(f'Error to generate interactive graph with Dash: {e}')

        return None

    def run(self, color_by: str | None = None):
        fingerprints = self.config.analysis['fingerprint']
        reductions = self.config.analysis['reduction_method']

        if isinstance(fingerprints, str):
            fingerprints = [fingerprints]
        if isinstance(reductions, str):
            reductions = [reductions]

        os.makedirs(self.config.viz['output_dir'], exist_ok=True)

        original_data = self.load_data()
        results = []
        port = 8050

        for fp in fingerprints:
            self.data = original_data.copy()
            self.compute_descriptors(fp_type=fp)
            self.compute_similarity()
            heatmap_title = f'Tanimoto Heatmap - {fp.upper()}'
            heatmap_path = f"{self.config.viz['output_dir']}/{fp}_heatmap.png"
            heatmap_figure = visualization.plot_similarity_heatmap(self.sim_matrix, labels=False, output_path=heatmap_path, show=True, title=heatmap_title)
            results.append((fp, 'heatmap', heatmap_figure))

            
            data_with_descriptors = self.data.copy()
            fps_for_fp = self.fps

            for red in reductions:
                self.data = data_with_descriptors.copy()
                self.fps = fps_for_fp
                reduced = self.reduce_dimensions(method=red)
                save_prefix = f'{fp}_{red}'
                heatmap_title = f'Tanimoto Heatmap - {fp.upper()}'
                scatter_title = f'{fp.upper()} + {red.upper()}'
                figs = self.visualize(save_prefix=save_prefix, show=False, title=scatter_title, heatmap_title=heatmap_title, interactive_mode=True, port=port, color_by=color_by)
                coords, results_eval, trust, corr, explained_variance = reduced
                results.append((fp, red, figs, results_eval))
                port += 3

        return results
