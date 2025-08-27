import pandas as pd
import os
from typing import List,Tuple
from .config import MayaConfig
from .utils import compute_activity_metrics

def load_dataset(file_path: str, config: MayaConfig) -> pd.DataFrame:
    """Load dataset from various formats and validate."""
    ext = os.path.splitext(file_path)[1].lower().lstrip('.')
    readers = {
        'csv': pd.read_csv, 'xlsx': pd.read_excel, 'tsv': lambda f: pd.read_csv(f, sep='\t'),
        'json': pd.read_json, 'xml': pd.read_xml
    }
    if ext not in readers:
        raise ValueError(f"Unsupported format: {ext}")
    df = readers[ext](file_path)
    required_cols = [config.data['id_col'], config.data['smiles_col']] + config.data['activities']
    missing = set(required_cols) - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}")
    if not df[config.data['id_col']].is_unique:
        raise ValueError("Duplicate IDs")
    for col in required_cols:
        if df[col].isnull().any():
            raise ValueError(f"Nulls in {col}")
        if col in config.data['activities'] and not pd.api.types.is_numeric_dtype(df[col]):
            raise TypeError(f"Non-numeric values in {col}")

    # Compute activity metrics
    df, metric_col = compute_activity_metrics(
        df=df,
        activities=config.data['activities'],
        eval_metric=config.data.get('eval_metric', 'mpIC50_value'),
        metric_name=config.data.get('metric_name', 'custom_metric'),
        use_custom_metric=True,
        create_new_column=True
    )

    return df
