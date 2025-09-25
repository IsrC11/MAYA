# maya/config.py
import os
from dataclasses import dataclass

@dataclass
class MayaConfig:
    def __init__(self, data_path: str = None, output_dir: str = './results', fingerprint: str | list = 'morgan', reduction_method: str | list = 'pca', properties: list = None, color_by: str = 'LogP', palette: str = 'viridis'):
        self.data_path = data_path
        self.data = {'id_col': None, 'smiles_col': None, 'activities': [], 'eval_metric': None, 'metric_value': None}
        self.analysis = {'fingerprint': fingerprint, 'reduction_method': reduction_method, 'properties': properties or ['MolWt', 'LogP', 'HBA', 'HBD', 'TPSA']}
        self.viz = {'output_dir': output_dir, 'color_by':color_by, 'palette':palette}
