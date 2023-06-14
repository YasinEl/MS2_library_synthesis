import molbloom
import pandas as pd
from molbloom import buy

#Read table
df = pd.read_csv("combined_datasources_yasin_harmonizedSMILES_hasSubst_summarized.csv", low_memory=False)

#Iterate over each row and use the SMILES to determine if the compound is commercially available according to molbloom
for index,row in df.iterrows():
    smiles = row['smiles_harmonized']
    df.at[index, 'commercially_available'] = buy(smiles, catalog='zinc20')

Uncomment if you want to save the intermedi
#df.to_csv("with_commercial_availability.csv", index=False)

#Exclude those that are not commercially available
commercially_available_df = df[df['commercially_available'] == True]

commercially_available_df.to_csv("only_commercially_available.csv", index=False)