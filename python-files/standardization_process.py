
#importing the cheminformatics module needed
import molvs
import pandas as pd
import numpy as np
from molvs import Standardizer, normalize
from rdkit import Chem
from tqdm import tqdm
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.drawOptions.comicMode=True
import rdkit
print(rdkit.__version__)

#import the datasets
# replicase = pd.read_csv('replicase_data_preprocessed.csv')
# replicase.head()
# M_pro = pd.read_csv('3cl-pro_data_preprocessed.csv')
# replicase.head()

#Define a canonical function to canonicalize the smiles
def Canonical(smiles):
    try:
        if smiles is not None:
            mol= Chem.MolFromSmiles(smiles, sanitize=True)
            if mol is not None:
            #convert the sanitized molecule to canonical smiles
                canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
                return canonical_smiles
            else:
                #log problematic smiles
                print(f"Invalid SMILES: {smiles}")
    except Exception as e:
        # log exceptions and problematics smiles
        print(f"Exception: {str(e)}, SMILES: {smiles}")
    return None
# Generate the canonical smiles for M_pro
# M_pro['canonical_smiles'] = M_pro['canonical_smiles'].apply(Canonical)
# M_pro = M_pro.dropna(subset=['canonical_smiles'])
# M_pro

# # Generate the canonical smiles for replicase
# replicase['canonical_smiles'] = replicase['canonical_smiles'].apply(Canonical)
# replicase

# #filter out rows where canonical smiles is not none
# replicase = replicase.dropna(subset=['canonical_smiles'])
# replicase

# M_pro.reset_index(drop=True, inplace=True)
# #reset the index
# replicase.reset_index(drop=True, inplace=True)

# # Create a Standardizer object for the standardization
# standardizer = Standardizer()

# # Modify the normalization rules
# norms = list(normalize.NORMALIZATIONS)
# for i in range(len(norms) - 1, 0, -1):
#     if norms[i].name == "Sulfoxide to -S+(O-)":
#         del norms[i]
# norms.append(normalize.Normalization("[S+]-[O-] to S=O", "[S+:1]([O-:2])>>[S+0:1](=[O-0:2])"))

# # Set the modified normalization rules
# standardizer.normalizations = norms

# # Set the prefer_organic option to True
# standardizer.prefer_organic = True
# #standardize the canonical smiles and the rdkit mol object of 3cl_pro
# stand_mol1 = []
# stand_smi1 = []
# for smi in M_pro['canonical_smiles'].tolist():
#     mol = Chem.MolFromSmiles(smi)  # Convert SMILES to RDKit Mol object
#     standardized_mol1 = standardizer.standardize(mol)
#     stand_mol1.append(standardized_mol1)
#     standardized_smi1 = Chem.MolToSmiles(standardized_mol1)
#     stand_smi1.append(standardized_smi1)

# #standardize the canonical smiles and the rdkit mol object of replicase
# stand_mol2 = []
# stand_smi2 = []
# for smi in replicase['canonical_smiles'].tolist():
#     mol = Chem.MolFromSmiles(smi)  # Convert SMILES to RDKit Mol object
#     standardized_mol2 = standardizer.standardize(mol)
#     stand_mol2.append(standardized_mol2)
#     standardized_smi2 = Chem.MolToSmiles(standardized_mol2)
#     stand_smi2.append(standardized_smi2)

# #put the standardized smiles in the two dataframe
# M_pro['standardized_smiles'] = pd.DataFrame(stand_smi1)
# replicase['standardized_smiles'] = pd.DataFrame(stand_smi2)
# #reset the index
# M_pro.reset_index(drop=True, inplace=True)
# #reset the index
# replicase.reset_index(drop=True, inplace=True)
# M_pro.to_csv("3cl-pro_stand_smi_data1.csv", index=False)
# replicase.to_csv("Replicase_stand_smi_data1.csv", index=False)

# replicase

def standard_preprocess(df):
    df['canonical_smiles'] = df['canonical_smiles'].apply(Canonical)
    df = df.dropna(subset=['canonical_smiles'])
    df.reset_index(drop=True, inplace=True)
    #standardize the canonical smiles and the rdkit mol object
    stand_mol2 = []
    stand_smi2 = []

    norms = list(normalize.NORMALIZATIONS)
    for i in range(len(norms) - 1, 0, -1):
        if norms[i].name == "Sulfoxide to -S+(O-)":
            del norms[i]
    norms.append(normalize.Normalization("[S+]-[O-] to S=O", "[S+:1]([O-:2])>>[S+0:1](=[O-0:2])"))

    standardizer = Standardizer()
    standardizer.normalizations = norms
    standardizer.prefer_organic = True

    for smi in df['canonical_smiles'].tolist():
        mol = Chem.MolFromSmiles(smi)  # Convert SMILES to RDKit Mol object
        standardized_mol2 = standardizer.standardize(mol)
        stand_mol2.append(standardized_mol2)
        standardized_smi2 = Chem.MolToSmiles(standardized_mol2)
        stand_smi2.append(standardized_smi2)
    df['standardized_smiles'] = pd.DataFrame(stand_smi2)
    df.reset_index(drop=True, inplace=True)
    # df.to_csv(output_filename, index=False)

    return df
