# maya/config.py
import os
from dataclasses import dataclass

@dataclass
class MayaConfig:
    def __init__(self, data_path: str = None, output_dir: str = './results', fingerprint: str | list = 'morgan', reduction_method: str | list = 'pca', properties: list = None, color_by: str = 'LogP', palette: str = 'RdBu_r', n_jobs: int = -1, id_col: str = 'ID', smiles_col: str = 'SMILES'):
        self.data_path = data_path
        self.n_jobs = n_jobs
        self.data = {'id_col': id_col, 'smiles_col': smiles_col, 'activities': [], 'eval_metric': None, 'metric_value': None}
        self.curation = {"standardize": True, "largest_fragment": True, "neutralize": False, "canonical_tautomer": False}
        self.analysis = {'fingerprint': fingerprint, 'reduction_method': reduction_method, 'properties': properties or ['MolWt', 'LogP', 'HBA', 'HBD', 'TPSA']}
        self.viz = {'output_dir': output_dir, 'color_by':color_by, 'palette':palette}
