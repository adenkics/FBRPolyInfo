from rdkit import Chem
from rdkit.Chem import AllChem

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

def FRstr(smi: str):
    FR_str = ['F', 'Cl', 'Br', 'I', 'S','P']
    Endstr = ['Fr', 'Rb']

    for end in Endstr:
        smi = smi.replace(end, '')
    for fr_elements in FR_str:
        if fr_elements in smi:
            return False
    return True

def Get_vocab(smi, vocab_list):
    mol = smi2mol(smi, kekulize=FRstr(smi))

    if mol.HasSubstructMatch(Chem.MolFromSmarts('[Fr]')):
        mol = DelEndGrupIdx(mol)
        Chem.SanitizeMol(mol)

    vocab_dict = {}
    if mol is not None:
        for vocab in vocab_list:
            vocab_feat = len(mol.GetSubstructMatches(Chem.MolFromSmarts(vocab)))
            vocab_dict[vocab] = vocab_feat

    return vocab_dict