
# Import necessary libraries
import pandas as pd
import numpy as np
from chembl_webresource_client.new_client import new_client

def main():
    # Target search for coronavirus
    target = new_client.target
    target_query = target.search('coronavirus')
    targets = pd.DataFrame.from_dict(target_query)
    targets

    #We will assign 6 (3cl-pro) to the *selected_target* variable
    selected_target0 = targets.target_chembl_id[6]
    print("selected_target0 :", selected_target0)

    #We will assign the 8,9 (Replicase polyprotein 1ab) to the *selected_target* variable
    selected_target1 = targets.target_chembl_id[8] #binding assays (CHEMBL5118)
    selected_target2 = targets.target_chembl_id[9] #functional assays (CHEMBL4523582)
    print("selected_target1 :", selected_target1)
    print("selected_target2 :", selected_target2)

    # Here, we will retrieve only bioactivity data for coronavirus 3cl-pro that are reported as IC
    # values in nM (nanomolar) unit.
    activity = new_client.activity
    res0 = activity.filter(target_chembl_id=selected_target0).filter(standard_type="IC50")

    # Here, we will retrieve only bioactivity data for coronavirus replicase polyprotein 1ab that are reported as IC
    # values in nM (nanomolar) unit.
    activity = new_client.activity
    res1 = activity.filter(target_chembl_id=selected_target1).filter(standard_type="IC50")
    res2 = activity.filter(target_chembl_id=selected_target2).filter(standard_type="IC50")

    #put them in a dataframe
    df = pd.DataFrame.from_dict(res0) #3cl-pro
    df1 = pd.DataFrame.from_dict(res1)# binding assays of replicase polyprotein
    df2 = pd.DataFrame.from_dict(res2)# functional assays of replicase polyprotein

    #Concatenate the dataframes on repliacase polyprotein vertically
    data = pd.concat([df1, df2], axis=0)
    # Reset the index of the resulting dataframe
    data = data.reset_index(drop=True)
    data.head(5)
    data.shape

    # #Finally we will save the resulting raw data of 3cl-pro and replicase polyprotein
    # df.to_csv('3cl-pro_data_raw.csv', index=False)
    #data.to_csv('replicase_data_raw.csv', index=False)

    #check to see if there are duplicates
    len(df.canonical_smiles.unique())
    #Drop the duplicates
    df = df.drop_duplicates(['canonical_smiles'])
    df.shape
    #Reset the index
    df.reset_index(drop=True, inplace=True)
    df.shape

    #check to see if there are duplicates
    len(data.canonical_smiles.unique())
    #Drop the duplicates
    data = data.drop_duplicates(['canonical_smiles'])
    data.shape
    #Reset the index
    data.reset_index(drop=True, inplace=True)
    data.shape

    #add 4.00 to the missing values in the pchembl_value column in the 3cl-pro datasets
    df['pchembl_value'] = df['pchembl_value'].fillna(4.00)

    #add 4.00 to the missing values in the pchembl_value column in the replicase datasets
    data['pchembl_value'] = data['pchembl_value'].fillna(4.00)

    # define the activity class using 4 as the pchembl value cut-off for 3cl-pro dataset
    bioactivity_class1 = []
    for i in df.pchembl_value:
        if float(i) <= 4:
            bioactivity_class1.append("no activity")
        elif float(i) > 4 and float(i) <= 5.99:
            bioactivity_class1.append("low activity")
        elif float(i) >= 6 and float(i) <= 7.99:
            bioactivity_class1.append("moderate activity")
        else:
            bioactivity_class1.append("high activity")

    len(bioactivity_class1)

    # define the activity class using 4 as the pchembl value cut-off for replicase dataset
    bioactivity_class2 = []
    for i in data.pchembl_value:
        if float(i) <= 4:
            bioactivity_class2.append("no activity")
        elif float(i) > 4 and float(i) <= 5.99:
            bioactivity_class2.append("low activity")
        elif float(i) >= 6 and float(i) <= 7.99:
            bioactivity_class2.append("moderate activity")
        else:
            bioactivity_class2.append("high activity")

    len(bioactivity_class2)

    #slice your data for 3cl-pro
    df = df[['molecule_chembl_id','canonical_smiles','pchembl_value']]
    df.shape
    #serialize the bioactivity list
    bioactivity_class1 = pd.Series(bioactivity_class1, name='bioactivity_class')
    # Reset the index of 'df' and 'bioactivity_class' to ensure consistent indices
    df.reset_index(drop=True, inplace=True)
    bioactivity_class1.reset_index(drop=True, inplace=True)

    # Concatenate 'mynewdata' and 'bioactivity_class' horizontally (axis=1)
    dataframe1 = pd.concat([df, bioactivity_class1], axis=1)
    print(dataframe1.columns)

    #slice your data for replicase
    data= data[['molecule_chembl_id','canonical_smiles','pchembl_value']]
    data.shape
    #serialize the bioactivity list
    bioactivity_class2 = pd.Series(bioactivity_class2, name='bioactivity_class')
    # Reset the index of 'data' and 'bioactivity_class' to ensure consistent indices
    data.reset_index(drop=True, inplace=True)
    bioactivity_class2.reset_index(drop=True, inplace=True)

    # Concatenate 'mynewdata' and 'bioactivity_class' horizontally (axis=1)
    dataframe2 = pd.concat([data, bioactivity_class2], axis=1)
    print(dataframe2.columns)

    # Get the unique canonical smiles of 3cl-pro
    unique_smiles = dataframe1['canonical_smiles'].unique()
    len(unique_smiles)
    #filter out the rows matching smiles in the 3cl-pro from the replicase datasets
    dataframe2 = dataframe2[~dataframe2['canonical_smiles'].isin(unique_smiles)]
    #reset the index of the dataframe
    dataframe2.reset_index(drop=True, inplace=True)
    print(dataframe2.shape)

    #save the two dataframes that are now the preprocessed data

    # dataframe1.to_csv('3cl-pro_data_preprocessed.csv', index=False)
    # dataframe2.to_csv('replicase_data_preprocessed.csv', index=False)


    return dataframe2, dataframe1

    # dataframe1.to_csv(_3clpro_file_name, index=False)
    # dataframe2.to_csv(replicase_file_name, index=False)