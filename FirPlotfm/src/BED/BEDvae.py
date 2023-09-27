import numpy as np
import torch
import torch.nn as nn

class BurnIntrEncoder(nn.Module):
    def __init__(self, dim_hidden_list, dim_latent):
        super(BurnIntrEncoder, self).__init__()
        self.dim_latent = dim_latent

        self.hidden_layers = nn.ModuleList()
        self.hidden_layers.append(nn.Linear(2048, dim_hidden_list[0]))  # 修改这里
        self.hidden_layers.append(nn.ReLU())
        for i in range(len(dim_hidden_list) - 1):
            self.hidden_layers.append(nn.Linear(dim_hidden_list[i], dim_hidden_list[i + 1]))
            self.hidden_layers.append(nn.ReLU())

        self.fc_mu = nn.Linear(dim_hidden_list[-1], dim_latent)
        self.fc_var = nn.Linear(dim_hidden_list[-1], dim_latent)

    def forward(self, x):
        for layer in self.hidden_layers:
            x = layer(x)

        mu = self.fc_mu(x)
        log_var = self.fc_var(x)

        return mu, log_var

class BurnIntrDecoder(nn.Module):
    def __init__(self, dim_latent, dim_hidden_list):
        super(BurnIntrDecoder, self).__init__()

        # 构建隐藏层
        self.hidden_layers = nn.ModuleList()
        prev_dim = dim_latent
        for dim in dim_hidden_list:
            self.hidden_layers.append(nn.Linear(prev_dim, dim))
            self.hidden_layers.append(nn.ReLU())
            prev_dim = dim

        # 输出层
        self.fc_out = nn.Linear(prev_dim, 2048)

    def forward(self, x):
        for layer in self.hidden_layers:
            x = layer(x)

        return torch.sigmoid(self.fc_out(x))

class VAE(nn.Module):
    def __init__(self, dim_hidden_list, dim_latent):
        super(VAE, self).__init__()

        self.encoder = BurnIntrEncoder(dim_hidden_list, dim_latent)
        self.decoder = BurnIntrDecoder(dim_latent, dim_hidden_list)

    def forward(self, x):
        mu, log_var = self.encoder(x)

        return mu


def LatntVectReSum(BED_weit_array_dict, dim_latent):

    if not BED_weit_array_dict:
        return torch.zeros(1, dim_latent)

    First_tensor = next(iter(BED_weit_array_dict.values()))
    Resum_tensor = torch.zeros_like(First_tensor)
    for weight, BED_latent_tensor in BED_weit_array_dict.items():
        weight_tensor = weight * BED_latent_tensor
        Resum_tensor += weight_tensor

    return Resum_tensor



def BEDLoss(InptFeat_dict, NewVect_dict):
    mse_loss = nn.MSELoss()
    losses = []

    for smiles, new_vector in NewVect_dict.items():
        for InptFeat_list in InptFeat_dict[smiles]:
            for InptFeat_tensor in InptFeat_list:
                loss = mse_loss(InptFeat_tensor, new_vector)
                losses.append(loss)
    average_loss = sum(losses) / len(losses)

    return average_loss



