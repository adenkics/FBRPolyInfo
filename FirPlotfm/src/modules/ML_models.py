import numpy as np
import torch
import os
import pickle
from FirPlotfm.src.modules.Desc import Desc2D, MorganFP2bit, FP2npAry
from rdkit import Chem
from FirPlotfm.src.BED.BEDfeature import BEDVAEFeat

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score


script_dir = os.path.dirname(os.path.abspath(__file__))
relative_path1 = os.path.join('..', '..', 'data', 'RDKit_desc_list.pkl')
pkl_file_path1 = os.path.abspath(os.path.join(script_dir, relative_path1))

relative_path2 = os.path.join('..', '..', 'data', 'BED_train.pkl')
pkl_file_path2 = os.path.abspath(os.path.join(script_dir, relative_path2))

relative_path3 = os.path.join('..', '..', 'data', 'BED_test.pkl')
pkl_file_path3 = os.path.abspath(os.path.join(script_dir, relative_path2))

with open(pkl_file_path1, 'rb') as file:
    Desc_list = pickle.load(file)

with open(pkl_file_path2, 'rb') as file:
    PolymerData = pickle.load(file)

with open(pkl_file_path2, 'rb') as file:
    Polymer_test_Data = pickle.load(file)

labels_dict = PolymerData.set_index('SMILES')['LOI'].to_dict()
labels_test_dict = Polymer_test_Data.set_index('SMILES')['LOI'].to_dict()

Desc_dict = {}
Morgan_dict = {}
for smiles in PolymerData['SMILES']:
    mol = Chem.MolFromSmiles(smiles)
    desc_data = Desc2D(mol, Desc_list['RDKit_Desc_list'])
    Desc_dict[smiles] = desc_data
    MorganFPbitStr = MorganFP2bit(mol, radius= 4, nBits=2048)
    MorganFPAry = FP2npAry(MorganFPbitStr)
    MorganFPAry = np.reshape(MorganFPAry, (1, 2048))
    Morgan_dict[smiles] = MorganFPAry

vae_dim_hidden_list = [1024, 1024, 512]
dim_latent = 256
BED_new_tensor = BEDVAEFeat(pkl_file_path2,vae_dim_hidden_list, dim_latent)
BEDary = {}
for smiles, value in BED_new_tensor.items():
    BED_tensor2ary = value.detach().numpy()
    BEDary[smiles] = BED_tensor2ary
concat_feat = {}
for key in Desc_dict.keys():
    concat_feat[key] = np.concatenate((Desc_dict[key], Morgan_dict[key], BEDary[key]),axis=1)


train_features = np.array(list(concat_feat.values()))
original_shape = train_features.shape
reshaped_features = train_features.reshape(original_shape[0], -1)

train_labels = np.array(list(labels_dict.values()))



rf_regressor = RandomForestRegressor()


param_grid = {
    'n_estimators': [100, 200, 300],  
    'max_depth': [None, 10, 20],     
    'min_samples_split': [2, 5, 10],  
    'min_samples_leaf': [1, 2, 4],    
    'max_features': ['auto', 'sqrt', 'log2'],
}

grid_search = GridSearchCV(estimator=rf_regressor, param_grid=param_grid,
                           scoring='neg_mean_squared_error', cv=5)

grid_search.fit(reshaped_features, train_labels)

best_params = grid_search.best_params_
best_score = np.sqrt(-grid_search.best_score_)
print(f"Best Parameters: {best_params}")
print(f"Best RMSE Score: {best_score:.4f}")


best_rf_regressor = RandomForestRegressor(**best_params)


kf = KFold(n_splits=5, shuffle=True, random_state=42) 


rmse_scores = np.sqrt(-cross_val_score(best_rf_regressor, reshaped_features, train_labels, scoring='neg_mean_squared_error', cv=kf))

r2_scores = cross_val_score(best_rf_regressor, reshaped_features, train_labels, scoring='r2', cv=kf)

for i, (rmse, r2) in enumerate(zip(rmse_scores, r2_scores), start=1):
    print(f"Fold {i}: RMSE = {rmse:.4f}, R^2 = {r2:.4f}")

print(f"Average RMSE = {np.mean(rmse_scores):.4f}")
print(f"Average R^2 = {np.mean(r2_scores):.4f}")


