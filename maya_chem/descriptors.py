# maya/descriptors.py
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from skfp.fingerprints import MAPFingerprint

def compute_morgan_fingerprint(smiles: str, radius: int = 2, nBits: int = 2048):
    """Compute Morgan fingerprint."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    gab = AllChem.GetMorganGenerator(radius=radius, fpSize=nBits)
    return gab.GetFingerprint(mol)

def coompute_maccs_fingerprin(smiles: str):
    """Compute MACCS kesys 166 bits"""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    fp_maccs = MACCSkeys.GenMACCSKeys(mol)
    return fp_maccs

def compute_map4_fingerprint(smiles: str):
    """compute MAP4 fingerprints"""
    map_4 =　MAPFingerprint()
    fg_map4 = map_4.transform(smiles)
    fp_map4_bitvect = [numpy_to_bitvect(fg_map4)]
    return fp_map4_bitvect

def compute_physicochemical_properties(smiles: str):
    """Compute basic molecular descriptors."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    return {
        "MolWt": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "TPSA": Descriptors.TPSA(mol),
    }
