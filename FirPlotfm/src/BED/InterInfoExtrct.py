
from rdkit import Chem
from rdkit.Chem import AllChem
from FirPlotfm.src.modules.Desc import FP2npAry
from FirPlotfm.src.utils.chem_utils import sma2mol, smi2mol,mol2smi
from FirPlotfm.src.utils.data_utils import VCGData
from rdkit.Chem import rdFMCS, rdMolDescriptors
from rdkit import DataStructs

def BEDInput(mol_dict: dict):
    VCGBEDIn_dict = {}

    for key, value in mol_dict.items():
        VCG_dict = {}
        keymol = Chem.MolFromSmiles(key)
        keyfp = AllChem.GetMorganFingerprintAsBitVect(keymol, radius=4, nBits=2048)

        for VCGmol in value:
            VCGfp = AllChem.GetMorganFingerprintAsBitVect(VCGmol, radius=4, nBits=2048)
            VCGSimilarity = DataStructs.TanimotoSimilarity(keyfp, VCGfp)
            BitString = VCGfp.ToBitString()
            VCGarray = FP2npAry(BitString)
            VCG_dict[VCGSimilarity] = VCGarray

        VCGBEDIn_dict[key] = VCG_dict

    return VCGBEDIn_dict




def FPDiffer(mol1, mol2):
    if mol2 is None:
        return None
    else:
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius = 4)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius = 4)

        return DataStructs.TanimotoSimilarity(fp1, fp2)

class BurIntrEx:
    @staticmethod
    def SplitDifferMol(OriginMol, ComparMol):
        if None in (OriginMol, ComparMol):
            return None

        mcs = rdFMCS.FindMCS([OriginMol, ComparMol])
        mcs_mol = sma2mol(mcs.smartsString)
        UniFrag = Chem.DeleteSubstructs(ComparMol, mcs_mol)

        return UniFrag

    @staticmethod
    def GetDifferMol(vcg_data_dict):
        differ_mol_dict = {}

        for key, value in vcg_data_dict.items():
            differ_mols = [BurIntrEx.SplitDifferMol(key, BuIntMol) for BuIntMol in value]
            differ_mol_dict[key] = differ_mols

        return differ_mol_dict

# class BurIntrEx:
#     @staticmethod
#     def SplitDifferMol(OriginMol, ComparMol):
#         if None in (OriginMol, ComparMol):
#             return None
#
#         mcs = rdFMCS.FindMCS([OriginMol, ComparMol])
#         mcs_mol = smi2mol(mcs.smartsString)
#         UniFrag = Chem.DeleteSubstructs(ComparMol, mcs_mol)
#
#         return UniFrag
#
#     @staticmethod
#     def GetDifferMol(vcg_data_dict):
#         differ_mol_dict = {}
#
#         for key, value in vcg_data_dict.items():
#             differ_mols = [BurIntrEx.SplitDifferMol(key, BuIntMol) for BuIntMol in value]
#             differ_mol_dict[key] = differ_mols
#
#         return differ_mol_dict

# def SplitDifferMol(OriginMol, ComparMol):
#     mcs = rdFMCS.FindMCS([OriginMol, ComparMol])
#     mcs_mol =sma2mol(mcs.smartsString)
#     UniFrag = Chem.DeleteSubstructs(ComparMol, mcs_mol)
#
#     return UniFrag
#
# def GetDifferMol(vcg_data_dict):
#     differ_mol_dict = {}
#
#     for key, value in vcg_data_dict.items():
#         differ_mols = []
#         for BuIntMol in value:
#             uni_frag = SplitDifferMol(key, BuIntMol)
#             differ_mols.append(uni_frag)
#
#         differ_mol_dict[key] = differ_mols
#
#     return differ_mol_dict

# vcg_data = VCGData()
# VCGMol = vcg_data.AplyVCG()
# print(GetDifferMol(VCGMol))


###### testt ######

# mol1 = smi2mol('CC1=CC=C(C=C1)C(=O)O')
# mol2 = smi2mol('CC1=CC=C(C=C1)C(=O)OCCC')
# print(FPDiffer(mol1, mol2))
#
# from rdkit.Chem import rdFMCS
# mol1 = Chem.MolFromSmiles("O=C(NCc1cc(OC)c(O)cc1)CCCC/C=C/C(C)C")
# mol2 = Chem.MolFromSmiles("CC(C)CCCCCC(=O)NCC1=CC(=C(C=C1)O)OC")
# mol3 = Chem.MolFromSmiles("c1(C=O)cc(OC)c(O)cc1")
# mols = [mol1,mol2,mol3]
# res=rdFMCS.FindMCS(mols)




