import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

bitsize = [1024, 2048]


def Desc2D(mol, desc_list:list):
    Desc_list = []
    calc_desc = MoleculeDescriptors.MolecularDescriptorCalculator(desc_list)
    header = list(calc_desc.GetDescriptorNames())
    if mol is None:
        Desc_list.append(['NaN'])
    else:
        descriptors=calc_desc.CalcDescriptors(mol)
        Desc_list.append(descriptors)

    desc_array = np.array(Desc_list, dtype=float)
    desc_array = desc_array.reshape(1, -1)

    return desc_array

def MorganFP2bit(mol, radius, nBits):
    BitString = str()
    if mol:
        Morganfp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
        BitString = Morganfp.ToBitString()
    return BitString if BitString else ''


def RDKFP2bit(mol, maxPath, fpSize, nBitsPerHash):
    if mol:
        RDKitfp = Chem.RDKFingerprint(mol, maxPath=maxPath, fpSize=fpSize, nBitsPerHash=nBitsPerHash)
        bitlist = [int(bit) for bit in RDKitfp.ToBitString()]

    return bitlist

def PatFP2bit(mol, fpSize):
    if mol:
        Patfp = Chem.PatternFingerprint(mol, fpSize=fpSize)
        bitlist = [int(bit) for bit in Patfp.ToBitString()]

    return bitlist

def TopTri2bit(mol, nBits, targetSize):
    if mol:
        TTfp = Chem.rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=nBits, targetSize=targetSize)
        bitlist = [int(bit) for bit in TTfp.ToBitString()]

    return bitlist


def FP2npAry(bitstr: str):
    if not bitstr:
        return np.zeros(shape=(1, 2048), dtype=int)
    npfp = np.frombuffer(bitstr.encode(), 'u1') - ord('0')
    # if np.all(npfp == 0):
    #     raise ValueError("Converted NumPy array contains only zeros")
    return npfp
