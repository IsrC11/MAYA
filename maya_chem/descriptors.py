# maya/descriptors.py
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, MACCSkeys
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from skfp.fingerprints import MAPFingerprint
import numpy as np

def compute_morgan_fingerprint(smiles: str, radius: int = 2, nBits: int = 2048):
    """Compute Morgan fingerprint."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    gen = AllChem.GetMorganGenerator(radius=radius, fpSize=nBits)
    return gen.GetFingerprint(mol)

def compute_maccs_fingerprint(smiles: str):
    """Compute MACCS kesys 166 bits"""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    return MACCSkeys.GenMACCSKeys(mol)

map_4 = MAPFingerprint()

def numpy_to_bitvect(arr:np.ndarray) -> ExplicitBitVect:
    bv = ExplicitbitVecct(len(arr))
    for i, bit in enumerate(arr):
        if bit:
            bv.SetBit(i)
    return bv

def compute_map4_fingerprint(smiles: str):
    """compute MAP4 fingerprints"""
    
    if not smiles:
        return None

    try:
       arr = map_4.transform([smiles])[0]
       return numpy_to_bitvect(arr)
    except Exception:
       return None

def compute_physicochemical_properties(smiles: str, selected_props=None):
    """Compute basic molecular descriptors."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    all_props = {
        "MolWt": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "TPSA": Descriptors.TPSA(mol),
    }

    if selected_props:
        return {k: v for k, v in all_props.items() if k in selected_props}
    return all_props

def compute_fingerprint(smiles: str, method: str = 'morgan'):
    method = method.lower()
    if method == 'morgan' or method == 'ecfp':
        return compute_morgan_fingerprint(smiles)
    elif method =='maccs':
        return compute_maccs_fingerprint(smiles)
    elif method == 'map4':
        return compute_map4_fingerprint(smiles)
    else:
        raise ValueError(f'Unknown fingerprint method:{method}')

