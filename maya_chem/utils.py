import pandas as pd
import numpy as np
import re
from typing import List, Tuple

def compute_activity_metrics(
    df: pd.DataFrame,
    activities: List[str],
    eval_metric: str = 'mpIC50_value',
    metric_name: str = 'custom_metric',
    use_custom_metric: bool = True,
    create_new_column: bool = True
) -> Tuple[pd.DataFrame, str]:
    """Compute activity metrics for the dataset, including default and custom metrics.

    Args:
        df: Input DataFrame with activity columns.
        activities: List of column names containing activity values.
        eval_metric: Expression or column name for custom metric (default: 'mpIC50_value').
        metric_name: Name for the custom metric column (default: 'custom_metric').
        use_custom_metric: Whether to compute the custom metric (default: True).
        create_new_column: Whether to create a new column for the custom metric (default: True).

    Returns:
        Tuple containing the updated DataFrame and the name of the final metric column.
    """
    final_metric_col = 'mpIC50_value'  # Default metric column

    if activities:
        df['mpIC50_value'] = -np.log10(df[activities]).mean(axis=1)
        df['std'] = df[activities].std(axis=1)
        df['norma'] = (df['std'] - df['std'].min()) / (df['std'].max() - df['std'].min() + 1e-6)
        df['normal_deviation'] = df['norma'] * 2 + 9
        df['normal_deviation'] = df['normal_deviation'].fillna(10)
        final_metric_col = 'mpIC50_value'  # Update if calculated

    if use_custom_metric and eval_metric != 'mpIC50_value':
        cols = re.findall(r'\b\w+\b', eval_metric)
        missing = [c for c in cols if c not in df.columns]
        if missing:
            raise ValueError(f"Missing columns: {missing}")

        # Replace 'activities' in the expression
        expr = eval_metric.replace('activities', f"df[{activities}]")

        # Check if eval_metric is a single column name
        is_single_column = len(cols) == 1 and cols[0] in df.columns and expr == cols[0]

        if is_single_column and not create_new_column:
            # Use existing column directly
            final_metric_col = eval_metric
        else:
            # Evaluate and create new column
            df[metric_name] = pd.eval(expr, engine='python')
            final_metric_col = metric_name

    # If not using custom metric, keep default
    if not use_custom_metric:
        final_metric_col = 'mpIC50_value'

    return df, final_metric_col

# Other utility functions (if any) remain unchanged
