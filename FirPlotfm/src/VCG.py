import copy
from rdkit import Chem
from rdkit.Chem import rdChemReactions as Reactions
from rdkit.Chem.FilterCatalog import FilterMatcher
import os
import pickle
from FirPlotfm.src.utils.chem_utils import MoleSimilr, smi2mol,mol2smi

pkl_file_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'reactDefs.pkl')


with open(pkl_file_path, 'rb') as file:
    loaded_data = pickle.load(file)

reactDefs = loaded_data

def is_valid_molecule(mol):
    return mol is not None and FilterMatcher.IsValid(mol)

def VCGene(mol):
    PyProds = set()
    for reaction in reactDefs:
        templ = Reactions.ReactionFromSmarts(reaction)
        templ.Initialize()
        if templ.IsMoleculeReactant(mol) is True:
            new_prods = templ.RunReactants((mol, ))

            for prods in new_prods:
                for prod in prods:
                    if Chem.MolToSmiles(prod) not in PyProds:
                        PyProds.add(Chem.MolToSmiles(prod))
                    # else:
                    #     final_prods.add(Chem.MolToSmiles(prod))
    if not PyProds:
        return None
    return PyProds

def DropSimiSmiles(prods_list):

    unique_smiles = set()
    to_keep_indices = []

    for idx, mol in enumerate(prods_list):
        smi = mol2smi(mol)

        if smi not in unique_smiles:
            unique_smiles.add(smi)
            to_keep_indices.append(idx)

    updated_prods_list = [prods_list[idx] for idx in to_keep_indices]

    return updated_prods_list

def VCGout(flat_new_prods: list):
    Prods_dict = {index: element for index, element in enumerate(flat_new_prods)}

    to_remove = []

    for idx1, mol1 in Prods_dict.items():
        for idx2, mol2 in Prods_dict.items():
            if idx1 != idx2:
                if MoleSimilr(mol1, mol2):
                    to_remove.append(idx2)

    for idx in reversed(sorted(to_remove)):
        del flat_new_prods[idx]

    return flat_new_prods


###########  test  ############
'''
mol = Chem.MolFromSmiles('COCCCCOC(=O)c1ccc(cc1)C(=O)OCCCCOC(=O)c1ccc(cc1)C(=O)C')
prods = VCGene(mol)
for smiles in prods:
    print(mol2smi(smiles))
    
>>> COC(=O)c1ccc(C(=O)OC)cc1
>>> CCCOC(=O)c1ccc(C(C)=O)cc1
>>> COC(=O)c1ccc(C(=O)OC)cc1
>>> CCCOC(=O)c1ccc(C(C)=O)cc1
>>> CCCOC(=O)c1ccc(C(C)=O)cc1
>>> COC(=O)c1ccc(C(=O)OC)cc1
>>> CCCOC(=O)c1ccc(C(C)=O)cc1
>>> COC(=O)c1ccc(C(=O)OC)cc1
'''

# VCG_data_path = os.path.join(os.path.dirname(__file__),'..', 'data', 'VCGdata.pkl')
#
# with open(VCG_data_path, 'rb') as file:
#     moldata = pickle.load(file)
#
# prods_list = []
# for smiles in moldata['SMILES']:
#     mol = Chem.MolFromSmiles(smiles)
#     prods = VCGene(mol)
#     prods_list.append(prods)
#
# with open(r'C:\Users\86135\Desktop\PolyFBR\FirPlotfm\data\VCGBurnInfo.pkl', 'wb') as output_file:
#     pickle.dump(prods_list, output_file)
#     print("All data outputed")
#
#
# mol1 = Chem.MolFromSmiles('O=C(c1ccc(F)cc1)c1cc(-c2ccc(P(=O)(c3ccc(F)cc3)c3ccc(-c4ccc(-c5ccc(P(=O)(c6ccc(F)cc6)c6ccc(-c7ccc(-c8ccc(P(=O)(c9ccc(F)cc9)c9ccc([Rb])cc9)cc8)cc7C(=O)c7ccc(F)cc7)cc6)cc5)cc4C(=O)c4ccc(F)cc4)cc3)cc2)ccc1[Fr]')
