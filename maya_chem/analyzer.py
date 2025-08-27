from typing import Optional
from .config import MayaConfig
from .data import load_dataset
from .utils import compute_activity_metrics
import pandas as pd

class MayaAnalyzer:
    def __init__(self, config: MayaConfig):
        self.config = config

    def run(self, file_path: str) -> Optional[pd.DataFrame]:
        """Run the MAYA analysis pipeline."""
        # Load dataset
        df = load_dataset(file_path, self.config)

        # Compute activity metrics
        df, metric_col = compute_activity_metrics(
            df=df,
            activities=self.config.data['activities'],
            eval_metric=self.config.data.get('eval_metric', 'mpIC50_value'),
            metric_name=self.config.data.get('metric_name', 'custom_metric'),
            use_custom_metric=True,
            create_new_column=True
        )

        # Rest of the pipeline (e.g., descriptors, similarity, visualization)
        # Use metric_col for further processing
        # ...

        return df
