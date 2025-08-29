# maya/config.py
from dataclasses import dataclass, field
from typing import List, Dict, Optional

@dataclass
class MayaConfig:
    """Configuration for MAYA chemical multiverse analysis."""
    data: Dict = field(default_factory=lambda: {
        'id_col': 'ID', 'smiles_col': 'SMILES', 'activities': [],
        'eval_metric': lambda df, activities: -np.log10(df[activities]).mean(axis=1)
    })
    
    descriptors: List[str] = field(default_factory=lambda: ['MACCS', 'ECFP', 'druglikeness', 'MAP4'])
    reduction: Dict = field(default_factory=lambda: {
        'methods': ['PCA', 't-SNE'], 'perplexity': 30, 'n_iter': 1000
    })
    viz: Dict = field(default_factory=lambda: {
        'palette': 'RdBu_r', 'point_size': 12.0, 'size_repr': 'normal_deviation',
        'shape': 'circle', 'output_dir': 'maya_results'
    })
    parallel: Dict = field(default_factory=lambda: {'n_jobs': -1})
    ecfp_radius: int = 3
    signaturizer_codes: List[str] = field(default_factory=list)
