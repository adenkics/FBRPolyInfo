import torch
from rdkit import Chem
from torch_geometric.data import Data

MAX_VALENCE = {'B': 3, 'Br': 1, 'C': 4, 'Cl': 1, 'F': 1, 'I': 1, 'N': 3, 'O': 2, 'P': 5, 'S': 6, 'Si': 4, 'Fr': 1, 'Rb': 1}
Bond_List = [None, Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE, Chem.BondType.AROMATIC,
             Chem.BondType.IONIC]


def Mol2Graph(mol):
    if mol is None:
        return None

    num_atoms = mol.GetNumAtoms()
    edge_indices = []
    edge_types = []
    node_features = []

    for atom in mol.GetAtoms():
        atom_type = atom.GetSymbol()
        valence = MAX_VALENCE.get(atom_type, 4)  # Default valence is 4 for unknown atoms
        node_features.append(valence)

    for bond in mol.GetBonds():
        start_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        bond_type = bond.GetBondType()
        edge_indices.extend([(start_idx, end_idx), (end_idx, start_idx)])  # Add edges in both directions
        edge_types.extend([Bond_List.index(bond_type), Bond_List.index(bond_type)])

    edge_indices = torch.tensor(edge_indices, dtype=torch.long).t().contiguous()
    edge_types = torch.tensor(edge_types, dtype=torch.long)
    node_features = torch.tensor(node_features, dtype=torch.float)

    graph_data = Data(x=node_features, edge_index=edge_indices, edge_attr=edge_types)

    return graph_data

#############  Test  ################
# Example of Mol2Graph

# mol = Chem.MolFromSmiles('Oc1ccc(Cc2ccc3c(c2)CC(c2ccccc2)NO3)c([Rb])c1CN(Cc1c(Cc2ccc3c(c2)CC(c2ccccc2)NO3)ccc(O)c1CN(Cc1c(Cc2ccc3c(c2)CC(c2ccccc2)NO3)ccc(O)c1CN(C[Fr])c1ccccc1)c1ccccc1)c1ccccc1')  # 用您自己的 SMILES 字符串替换此处的示例
# graph_data = Mol2Graph(mol)
#
# if graph_data is not None:
#     print("Node Features (Max Valence):", graph_data.x)
#     print("Edge Indices:", graph_data.edge_index)
#     print("Edge Types:", graph_data.edge_attr)
# else:
#     print("Invalid SMILES string.")
