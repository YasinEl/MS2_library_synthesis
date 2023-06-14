import molbloom
import pandas as pd
from molbloom import buy
from rdkit import Chem
import rdkit

#Read table
df = pd.read_csv("database_accumulation.csv", low_memory=False)

#Iterate over each row and use the SMILES to determine if the compound is commercially available according to molbloom
for index,row in df.iterrows():
    smiles = row['smiles']
    try:
        df.at[index, 'commercially_available'] = buy(smiles, catalog='zinc-instock', canonicalize=True)
    except: 
        df.at[index, 'commercially_available'] = "SMILES_error"
        continue
    

#Uncomment if you want to save the intermediate step: the table that contains all the commercially available and not
df.to_csv("with_commercial_availability.csv", index=False)

#Exclude those that are not commercially available
#commercially_available_df = df[df['commercially_available'] == True]

#Save table of commercially available compounds
#commercially_available_df.to_csv("only_commercially_available.csv", index=False)