# maya/curation.py

from rdkit import Chem
from molvs import Standardizer, LargestFragmentChooser, Uncharger, TautomerCanonicalizer

class MayaCuration:

    def __init__(self, config):
        self.config = config
        self.standardizer = Standardizer()
        self.lfc = LargestFragmentChooser()
        self.uncharger = Uncharger()
        self.tautomer = TautomerCanonicalizer()

    def process_molecule(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        if self.config["standardize"]:
            mol = self.standardizer.standardize(mol)

        if self.config["largest_fragment"]:
            mol = self.lfc.choose(mol)

        if self.config["neutralize"]:
            mol = self.uncharger.uncharge(mol)

        if self.config["canonical_tautomer"]:
            mol = self.tautomer.canonicalize(mol)

        return Chem.MolToSmiles(mol, canonical=True)
