from rdkit import Chem
from rdkit.Chem import Draw,rdFMCS, AllChem
from rdkit.Chem import DataStructs
from rdkit.Chem import rdFingerprintGenerator

MAX_VALENCE = {'B': 3, 'Br': 1, 'C': 4, 'Cl': 1, 'F': 1, 'I': 1, 'N': 3, 'O': 2, 'P': 5, 'S': 6}
Bond_List = [None, Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE, Chem.BondType.AROMATIC,
             Chem.BondType.IONIC]

def smi2mol(smiles: str, kekulize = True):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        if kekulize:
            Chem.Kekulize(mol)

    return mol

def mol2smi(mol):
    return Chem.MolToSmiles(mol)

def sma2mol(sma):
    return Chem.MolFromSmarts(sma)

def mol2sma(mol):
    return Chem.MolToSmarts(mol)

def rec(mol):
    return smi2mol(mol2smi(mol), kekulize=False)

def CountProdAtoms(mol_list):
    atom_cnt_dict = {}
    for idx, mol in enumerate(mol_list):
        AtomNum = mol.GetNumAtoms()
        atom_cnt_dict[idx] = AtomNum

    return atom_cnt_dict

def CountProdsBonds(mol_list):
    bond_cnt_dict = {}
    for idx, mol in enumerate(mol_list):
        BondNum = prod.GetNumBonds()
        bond_cnt_dict[idx] = BondNum

    return bond_cnt_dict


def MoleSimilr(mol1, mol2):

    AllChem.Compute2DCoords(mol1)
    AllChem.Compute2DCoords(mol2)
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius = 4, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius = 4, nBits=2048)

    if mol1.GetNumAtoms() > 5 and mol2.GetNumAtoms() > 5:
        if mol1.HasSubstructMatch(mol2) or mol2.HasSubstructMatch(mol1):
            return True
        elif rdFMCS.FindMCS((mol1, mol2)).numAtoms > min(mol1.GetNumAtoms(), mol2.GetNumAtoms()) / 3 * 2 and DataStructs.TanimotoSimilarity(fp1,fp2) > 0.45:
            return True

    return False

def cnt_atom(smi):
    cnt = 0
    for c in smi:
        if c in MAX_VALENCE:
            cnt += 1
    return cnt

def shortest_path_len(i, j, mol):
    queue = Queue()
    queue.put((mol.GetAtomWithIdx(i), 1))
    visited = {}
    visited[i] = True
    while not queue.empty():
        atom, dist = queue.get()
        neis = []
        for nei in atom.GetNeighbors():
            idx = nei.GetIdx()
            if idx == j:
                return dist + 1
            if idx not in visited:
                visited[idx] = True
                neis.append(idx)
                queue.put((mol.GetAtomWithIdx(idx), dist + 1))
    return None


def cycle_check(i, j, mol):
    cycle_len = shortest_path_len(i, j, mol)
    return cycle_len is None or (cycle_len > 4 and cycle_len < 7)






###########   Test   ############

# mol = Chem.MolFromSmiles('CCOc1ccccc1C(=O)CC')
# mol2 = Chem.MolFromSmiles('CCOc1ccccc1C(=O)OCCOc1ccccc1')
#
# fpgen = AllChem.GetMorganFingerprint
# fp1 = fpgen(mol, radius=4)
# fp2 = fpgen(mol2, radius=4)
# print(DataStructs.DiceSimilarity(fp1,fp2))
# result = mol2.HasSubstructMatch(mol)
# print(result)
#
# res = rdFMCS.FindMCS((mol, mol2))
# print(res.numAtoms)
# print(MoleSimilr(mol, mol2))
