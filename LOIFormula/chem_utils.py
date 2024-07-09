from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

def smi2mol(smiles: str, kekulize=False, sanitize=True):
    '''turn smiles to molecule'''
    mol = Chem.MolFromSmiles(smiles, sanitize=sanitize)
    if kekulize:
        Chem.Kekulize(mol, True)
    return mol

def mol2smi(mol, canonical=True):
    return Chem.MolToSmiles(mol, canonical=canonical)

def sma2mol(sma: str, kekulize=False):
    '''turn smarts to molecule'''
    mol = Chem.MolFromSmarts(sma)
    if kekulize:
        Chem.Kekulize(mol, True)
    return mol

def DelEndGrupIdx(mol, kekulize=True):
    if kekulize:
        Chem.Kekulize(mol)
    new_mol = AllChem.DeleteSubstructs(mol, Chem.MolFromSmarts('[Fr]'), onlyFrags=False)
    new_mol = AllChem.DeleteSubstructs(new_mol, Chem.MolFromSmarts('[Rb]'), onlyFrags=False)
    return new_mol

def InChI2mol(inchi: str, kekulize=False):
    '''turn inchi to molecule'''
    mollist = [Chem.MolFromInchi(inchi)]
    return mollist

def SDF2mol(sdf):
    '''turn sdf to molecule'''
    SDF_dict = {}
    supplier = Chem.SDMolSupplier(sdf)
    for index, mol in enumerate(supplier):
        SDF_dict[index] = mol
    return SDF_dict

def MDL2mol(mdl):
    '''turn mdl to molecule'''
    MDL_dict = {}
    supplier = Chem.MolFromMolBlock(mdl)
    for index, mol in enumerate(supplier):
        MDL_dict[index] = mol
    return MDL_dict

def csv2mol(csv, smiles_col):
    '''turn csv to molecule'''
    df = pd.read_csv(csv)
    CSV_dict = {}
    for index, smi in enumerate(df[smiles_col]):
        mol = smi2mol(smi)
        CSV_dict[index] = mol
    return CSV_dict

def txt2mol(txt):
    '''turn txt to molecule'''
    with open(txt, 'r') as fin:
        lines = fin.read().strip().split('\n')
    TXT_dict = {}
    for index, smi in enumerate(lines):
        mol = smi2mol(smi)
        TXT_dict[index] = mol
    return TXT_dict

def Get_group(mol, group_list):
    group_dict = {}
    if mol is None:
        return None
    else:
        if mol.HasSubstructMatch(Chem.MolFromSmarts('[#0]')):
            mol = DelEndGrupIdx(mol)
            Chem.SanitizeMol(mol)
        if mol is not None:
            for group in group_list:
                group_feat = len(mol.GetSubstructMatches(Chem.MolFromSmarts(group)))
                group_dict[group] = group_feat

        return group_dict

