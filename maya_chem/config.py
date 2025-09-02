# maya/config.py
import os
from dataclasses import dataclass

@dataclass
class MayaConfig:
    def __init__(self, data_path: str = None):
        
        self.data_path = data_path
        
        self.data = {'id_col': None, 'smiles_col': None, 'activities': [], 'eval_metric': None, 'metric_value': None}

        self.analysis = {'descriptors': [], 'reduction_method': 'pca'}

        self.viz = {'output_dir': './results'}
