# maya/curation.py

from rdkit import Chem
from molvs import Standardizer, LargestFragmentChooser, Uncharger, TautomerCanonicalizer

class MayaCuration:
    """
    Configurable molecular curation pipeline for MAYA.
    """

    def __init__(self, curation_config: dict):
        self.config = curation_config

        # Initialize MolVS tools only if needed
        self.standardizer = Standardizer() if self.config.get("standardize", False) else None
        self.lfc = LargestFragmentChooser() if self.config.get("largest_fragment", False) else None
        self.uncharger = Uncharger() if self.config.get("neutralize", False) else None
        self.tautomer = TautomerCanonicalizer() if self.config.get("canonical_tautomer", False) else None

    def process_molecule(self, smiles: str) -> str | None:
        """
        Apply configured curation steps to a SMILES string.
        Returns canonical curated SMILES or None if invalid.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        try:
            # Standardization (normalize, reionize, disconnect metals)
            if self.standardizer:
                mol = self.standardizer.standardize(mol)

            # Keep largest fragment (remove salts/solvents)
            if self.lfc:
                mol = self.lfc.choose(mol)

            # Neutralize if requested
            if self.uncharger:
                mol = self.uncharger.uncharge(mol)

            # Canonical tautomer if requested
            if self.tautomer:
                mol = self.tautomer.canonicalize(mol)

            return Chem.MolToSmiles(mol, canonical=True)

        except Exception:
            return None


def clean_dataset(df, smiles_col: str, curation_config: dict):
    """
    Apply curation pipeline to entire dataframe.
    """
    curator = MayaCuration(curation_config)

    df["Canonical_Smiles"] = df[smiles_col].apply(curator.process_molecule)

    # Remove invalid molecules
    df = df.dropna(subset=["Canonical_Smiles"])

    # Remove duplicates after curation
    df = df.drop_duplicates(subset=["Canonical_Smiles"])

    return df
