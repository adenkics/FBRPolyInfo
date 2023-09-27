import os
import pickle
from rdkit import Chem
from FirPlotfm.src.VCG import VCGene
def Desc_list():

    pkl_file_path = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'RDKit_desc_list.pkl')

    with open(pkl_file_path, 'rb') as file:
        loaded_list = pickle.load(file)

    if not isinstance(loaded_list, list):
        loaded_list = [loaded_list]

    return loaded_list

class BED_data:
    def __init__(self):

        self.train_path = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'BED_train.pkl')
        self.test_path = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'BED_test.pkl')

        self.train_data = self.load_data(self.train_path)
        self.test_data = self.load_data(self.test_path)

    def load_data(self, file_path):

        with open(file_path, 'rb') as file:
            data = pickle.load(file)
        return data

    def TrainMol(self):
        train_mol = self.train_data['SMILES'].apply(smi2mol)
        return train_mol

    def TrainBurnProps(self):

        return self.train_data['LOI']

    def TestMol(self):
        test_mol = self.test_data['SMILES'].apply(smi2mol)
        return self.test_data['SMILES']

    def TestBurnProps(self):

        return self.test_data['LOI']

# bed_data = BED_data()
# train_smiles = bed_data.TrainMol()
# train_burn_props = bed_data.TrainBurnProps()
# test_smiles = bed_data.TestMol()
# test_burn_props = bed_data.TestBurnProps()

class VCGData:
    def __init__(self):
        self.VCG_data = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'VCGdata_test.pkl')

    def VCGLoader(self):
        with open(self.VCG_data, 'rb') as file:
            loaded_data = pickle.load(file)
        # VCGMole_list = loaded_data['SMILES'].apply(Chem.MolFromSmiles)
        VCGMole_list = [Chem.MolFromSmiles(smi) for smi in loaded_data['SMILES'] if Chem.MolFromSmiles(smi) is not None]
        smilist = list(loaded_data['SMILES'])
        return smilist, VCGMole_list

    def AplyVCG(self, loaded_data=None):
        if loaded_data is None:
            smi_list, mol_list = self.VCGLoader()
        else:
            smi_list = loaded_data['SMILES']
            mol_list = loaded_data['SMILES'].apply(Chem.MolFromSmiles)

        molecule_data = {}
        for smiles, mol in zip(smi_list, mol_list):
            flat_new_prods = VCGene(mol)

            molecule_data[smiles] = flat_new_prods

        return molecule_data



# class VCGData:
#     def __init__(self):
#         self.VCG_data = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'VCGdata_test.pkl')
#
#     def VCGLoader(self):
#         with open(self.VCG_data, 'rb') as file:
#             loaded_data = pickle.load(file)
#         VCGMole_list = loaded_data['SMILES'].apply(smi2mol)
#         return VCGMole_list
#
#     def AplyVCG(self, molecule_list=None):
#         if molecule_list is None:
#             molecule_list = self.VCGLoader()
#
#         molecule_data = {}
#
#         for mol, smiles in zip(molecule_list, loaded_data['SMILES']):
#             flat_new_prods = VCGene(mol)
#
#             molecule_data[smiles] = flat_new_prods
#
#         return molecule_data

    # def AplyVCG(self, molecule_list=None):
    #     if molecule_list is None:
    #         molecule_list = self.VCGLoader()
    #
    #     molecule_data = {}
    #
    #     for mol in molecule_list:
    #         flat_new_prods = VCGene(mol)
    #
    #         molecule_data[mol] = flat_new_prods
    #
    #     return molecule_data


# vcg_data = VCGData()
# VCGMol = vcg_data.AplyVCG()
# print(VCGMol)


