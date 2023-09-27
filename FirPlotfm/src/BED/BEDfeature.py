import os
import pickle
import torch
from rdkit import Chem
from FirPlotfm.src.BED.BEDvae import VAE, LatntVectReSum
from FirPlotfm.src.BED.InterInfoExtrct import BEDInput
from FirPlotfm.src.modules.Desc import MorganFP2bit, FP2npAry
from FirPlotfm.src.utils.data_utils import VCGData


def BEDVAEFeat(pkl_file_path, vae_dim_hidden_list,vae_dim_latent):

    with open(pkl_file_path, 'rb') as file:
        BED_Data = pickle.load(file)

    vcg_data = VCGData()
    VCGinfo = vcg_data.AplyVCG(BED_Data)
    VCG_FP_tensor = {}
    VCGmol_dict = {}

    for key, value in VCGinfo.items():
        mollist = [Chem.MolFromSmiles(smiles) for smiles in value if Chem.MolFromSmiles(smiles) is not None]
        VCGmol_dict[key] = mollist
        VCGBitStr = [MorganFP2bit(mol, radius = 4, nBits = 2048) for mol in mollist]
        VCGFParray = [FP2npAry(Bitstr) for Bitstr in VCGBitStr]
        VCGFPtensor = [torch.from_numpy(arr) for arr in VCGFParray]
        VCG_FP_tensor[key] = VCGFPtensor

    BED_input = BEDInput(VCGmol_dict)
    vae = VAE(dim_hidden_list=vae_dim_hidden_list, dim_latent=vae_dim_latent)
    latent_dict = {}
    for smiles, weit_ary_dict in BED_input.items():
        latent_vectors_dict = {}
        for weit, feature in weit_ary_dict.items():
            feature = feature.reshape(1, -1)
            mu, _ = vae.encoder(torch.from_numpy(feature).float())
            latent_vectors_dict[weit] = mu
        latent_resum = LatntVectReSum(latent_vectors_dict, dim_latent=vae_dim_latent)
        latent_dict[smiles] = latent_resum

    BED_new_tensor = {}
    for smiles, new_latent_vector in latent_dict.items():
        new_vector = vae.decoder(new_latent_vector)
        BED_new_tensor[smiles] = new_vector

    return BED_new_tensor


#####################   test code  #####################
# script_dir = os.path.dirname(os.path.abspath(__file__))
# relative_path = os.path.join('..', '..', 'data', 'BED_test_code.pkl')
# pkl_file_path = os.path.abspath(os.path.join(script_dir, relative_path))
#
# vae_dim_hidden_list = [1024, 512, 512]
# dim_latent = 256
# BED_new_tensor = BEDVAEFeat(pkl_file_path,vae_dim_hidden_list, dim_latent)