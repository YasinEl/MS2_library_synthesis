import pandas as pd
from chemparse import parse_formula
import pyteomics.mass
from itertools import product

df = pd.read_csv('output2.csv', low_memory=False)
delta_masses = pd.read_csv('delta_mass_counts_rationale.csv', low_memory=False)

#Exclude all rows with ACP, which correspond to enzymes or enzyme-metabolite complexes
filtered_df = df[~df['id'].str.contains('ACP')].copy()

for index,row in filtered_df.iterrows():
    charge = abs(int(row['charge']))
    formula_str = row['formula']
    formula_dict = parse_formula(formula_str)
    
    if charge == 0: #This is a charged molecule
        filtered_df.loc[index, 'neutral_formula'] = ''.join([f"{element}{round(count)}" for element, count in formula_dict.items()])
        monoisotopic_mass = pyteomics.mass.calculate_mass(formula_str)
        filtered_df.loc[index, 'monoisotopic_mass'] = round(monoisotopic_mass, 6)
    else:
        formula_dict['H']+= charge
        neutral_formula = ''.join([f"{element}{round(count)}" for element, count in formula_dict.items()])
        filtered_df.loc[index, 'neutral_formula'] = neutral_formula
        monoisotopic_mass = pyteomics.mass.calculate_mass(neutral_formula)
        filtered_df.loc[index, 'monoisotopic_mass'] = round(monoisotopic_mass, 6)

matched_df.to_csv("with_neutral_formulas.csv", index=False)
        
# Set absolute mass tolerance
mass_tolerance = 0.02

# Generate all possible combinations of rows between df1 and df2
combinations = list(product(filtered_df.iterrows(), delta_masses.iterrows()))

# Create an empty list to store matched rows
matched_rows = []

# Iterate over the combinations and compare the monoisotopic masses
for (index1, row1), (index2, row2) in combinations:
    mass1 = row1['monoisotopic_mass']
    mass2 = row2['mz_delta']
    mass_difference = abs(mass1 - mass2)

    if mass_difference <= mass_tolerance:
        matched_rows.append((index1, index2))

# Create a DataFrame from the matched rows
matched_df = pd.DataFrame(matched_rows, columns=['Index1', 'Index2'])
matched_df = pd.merge(filtered_df, matched_df, left_index=True, right_on='Index1')
matched_df = pd.merge(delta_masses, matched_df, left_index=True, right_on='Index2')

matched_df.to_csv("delta_masses_matched.csv", index=False)
