import os
import requests
import numpy as np
from rdkit.Chem import MolStandardize, rdMolDescriptors, MolFromSmiles
from rdkit import Chem
from tqdm import tqdm
import functools
from multiprocessing import Pool
import pandas as pd
import argparse
import warnings

def check_folder_access(paths):
    folder_paths = [os.path.dirname(path) for path in paths]  # Extract parent folder paths
    inaccessible_folders = []
    for path in folder_paths:
        if not os.path.isdir(path):
            inaccessible_folders.append(path)
        elif not os.access(path, os.W_OK):
            inaccessible_folders.append(path)

    if inaccessible_folders:
        error_message = "The following folders are inaccessible for writing:\n"
        for folder in inaccessible_folders:
            error_message += f"- {folder}\n"
        raise ValueError(error_message)

def check_csv_files(csv_paths):
    unreadable_files = []
    for path in csv_paths:
        if not os.path.exists(path):
            unreadable_files.append(path)
        else:
            try:
                pd.read_csv(path)
            except pd.errors.ParserError:
                unreadable_files.append(path)

    if unreadable_files:
        error_message = "The following CSV files could not be read:\n"
        for file in unreadable_files:
            error_message += f"- {file}\n"
        raise ValueError(error_message)


def count_substructure(smiles, smarts):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0

    # Define the SMARTS pattern
    pattern = Chem.MolFromSmarts(smarts)

    matches = mol.GetSubstructMatches(pattern)
    return len(matches)

def get_properties(smiles):
    mol = MolFromSmiles(smiles)
    if mol is None:
        return (None, None)
    mf = rdMolDescriptors.CalcMolFormula(mol)
    monomass = rdMolDescriptors.CalcExactMolWt(mol)
    return (mf, monomass)

#probably dont need to use this anymore. use harmonize_smiles_rdkit instead
def neutralize_atoms(smiles, remove_stereo = True):
    mol = Chem.MolFromSmiles(smiles)
    if remove_stereo == True:
        Chem.RemoveStereochemistry(mol)
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
            mol = Chem.rdmolops.RemoveHs(mol)
    return Chem.MolToSmiles(mol)

#probably dont need to use this anymore. use harmonize_smiles_rdkit instead
def harmonize_smiles(smiles_string):
    try:
        # take the largest covalently bound molecule
        smiles_string_new = MolStandardize.fragment.LargestFragmentChooser(smiles).prefer_organic
        # neutralize atoms
        new_smiles_string = neutralize_atoms(smiles_string_new)
        # remove implicit and explicit hydrogens and generate formulas
        mol = Chem.rdmolops.RemoveHs(Chem.MolFromSmiles(new_smiles_string))
        # make sure the SMILES is in rdkit format
        new_smiles_string = Chem.MolToSmiles(mol)

        return new_smiles_string
    except Exception as e:
        print(f"An error occurred with input {smiles_string}: {e}")
        return ""


def harmonize_smiles_rdkit(smiles, tautomer_limit = 900):
    try:
        # take the largest covalently bound molecule
        smiles_largest = MolStandardize.fragment.LargestFragmentChooser(smiles).prefer_organic
        mol = Chem.MolFromSmiles(smiles_largest)

        monomass = rdMolDescriptors.CalcExactMolWt(mol)
        # standardize tautomer
        if monomass < tautomer_limit:
            smiles_largest = MolStandardize.canonicalize_tautomer_smiles(smiles_largest)
            mol = Chem.MolFromSmiles(smiles_largest)

        # remove unnecessary charges
        uc = MolStandardize.charge.Uncharger()
        uncharged_mol = uc.uncharge(mol)

        # standardize the molecule
        lfc = MolStandardize.fragment.LargestFragmentChooser()
        standard_mol = lfc.choose(uncharged_mol)

        # remove stereochemistry
        Chem.RemoveStereochemistry(standard_mol)

        # get the standardized SMILES
        standard_smiles = Chem.MolToSmiles(standard_mol)
        return standard_smiles

    except Exception as e:
        print(f"An error occurred with input {smiles}: {e}")
        return ""

#not in use yet
def get_smiles_from_inchikey(inchikey):
    if pd.isnull(inchikey):
        return np.nan
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return data["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
    else:
        return np.nan


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Take a csv with SMILES and add columns with harmonized SMILES. There are options to also add columns with molecular formulas, exact masses and substructure counts.'
                    'Note that we will always select the longest organic part of every SMILES. Other parts of any SMILES will be ignored.')
    parser.add_argument("--input_csv", "-input", type=str,
                        help='Give a csv with at minimum one column named "smiles". Other columns are kept.',
                        required=True)
    parser.add_argument("--output_csv", "-output", type=str,
                        help='Give a path to a output csv.',
                        required=True)
    parser.add_argument("--intermedResults", "-interR", type=str,
                        help='If you want intermediate results to be exported after each step give a path where those csvs should be put. They are named automatically.',
                        default = '',
                        required=False)
    parser.add_argument("--mf_em", "-mf_em", type=bool,
                        help='Should columns with exact mass and molecular formula be added?',
                        default=True,
                        required=False)
    parser.add_argument("--substructure_csv", "-subSt", type=str,
                        help='Give a csv with substructure SMARTS to be counted. Needs to have a column "smarts_name" and a column "smarts".',
                        default='',
                        required=False)
    parser.add_argument("--TautomerLimit", "-tautLim", type=int,
                        help='For SMILES with a mass above this we wont attempt to harmonize tautomers. We dont recommend to attempt molecules > 1000 Da because this would take very long.'
                             'Note that changing this value might lead to slightly different harmonized SMILES!',
                        default=900,
                        required=False)
    parser.add_argument("--test", "-test", type=bool,
                        help='Is this just a test? If yes you can just process 1 percent of the input rows.',
                        default=False,
                        required=False)
    parser.add_argument("--warnings", "-warnings", type=bool,
                        help='Should warnings be displayed?',
                        default=False,
                        required=False)

    args = parser.parse_args()

    input_csv_path = args.input_csv
    output_csv_path = args.output_csv
    intermedResults_path = args.intermedResults
    add_mf_em = args.mf_em
    substructure_csv_path = args.substructure_csv
    tautomer_limit = args.TautomerLimit
    is_test = args.test
    display_warnings = args.warnings

    if display_warnings == False:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")

    #make sure all files can be accessed
    check_csv_files([input_csv_path, substructure_csv_path])
    check_folder_access([output_csv_path, intermedResults_path])

    #read input csv
    df = pd.read_csv(input_csv_path)

    if is_test:
        # Set the random seed for reproducibility
        np.random.seed(42)
        # take only 1% of dataframe rows
        df = df.sample(frac=0.005)


    #start process

    #
    #Harmonize SMILES. (normalize tautomers, uncharge and standardize)
    #

    # get unique values from the 'smiles' column
    unique_smiles = df.loc[
        (df['smiles'].notna()) & (df['smiles'] != ''), 'smiles'].unique().tolist()

    if len(unique_smiles) > 0:


        partial_func = functools.partial(harmonize_smiles_rdkit, tautomer_limit=tautomer_limit)
        # Create a pool of workers. The number of workers is usually set to the number of CPUs
        # but you can adjust this number based on your specific situation.
        pool = Pool()

        # Use the pool's map function to apply your_function to each unique smiles string.
        # We'll use tqdm to show a progress bar.
        results = list(tqdm(pool.imap(partial_func, unique_smiles), total=len(unique_smiles)))

        # create a dictionary mapping from the unique 'smiles' values to their function outputs
        mapping_dict = dict(zip(unique_smiles, results))

        # create a new column by mapping the 'smiles' column using the dictionary
        df['smiles_harmonized'] = df['smiles'].map(mapping_dict)

        # Make sure to close the pool when you're done to free up system resources.
        pool.close()
        pool.join()
        if len(intermedResults_path) > 1:
            df.to_csv(os.path.join(intermedResults_path, 'harmonized_smiles.csv'), index=False)

    #
    # Add molecular formula and exact mass
    #
    if add_mf_em:

        if 'mf' not in df.columns:
            df['mf'] = np.nan
            df['mf'] = df['mf'].astype(str)
        else:
            df['mf'] = df['mf'].astype(str)
        if 'monomass' not in df.columns:
            df['monomass'] = np.nan
        else:
            df['monomass'] = df['monomass'].astype(float)

        mask = (
                (df['mf'].isna() | df['monomass'].isna() | df['mf'] == '')
                & df['smiles_harmonized'].notna()
                & (df['smiles_harmonized'] != '')
        )

        unique_smiles = df.loc[mask, 'smiles_harmonized'].unique()
        print(len(unique_smiles))
        if len(unique_smiles) > 0:
            pool = Pool()

            results = list(tqdm(pool.imap(get_properties, unique_smiles), total=len(unique_smiles)))

            properties_dict = dict(zip(unique_smiles, results))

            df.loc[mask, 'mf'] = df.loc[mask, 'smiles_harmonized'].map(lambda x: properties_dict[x][0])
            df.loc[mask, 'monomass'] = df.loc[mask, 'smiles_harmonized'].map(lambda x: properties_dict[x][1])

            pool.close()
            pool.join()
            if len(intermedResults_path) > 1:
                df.to_csv(os.path.join(intermedResults_path, 'added_mf_and_monomass.csv'), index=False)


    #
    #Search for subgroups and add counts!
    #
    if len(substructure_csv_path) > 1:

        df_smarts = pd.read_csv(substructure_csv_path)

        if len(df_smarts) > 0:

            unique_smiles = df['smiles_harmonized'].dropna().unique()

            if len(unique_smiles) > 0:
                smarts_dict = df_smarts.set_index('smarts_name')['smarts'].to_dict()

                # Loop over the dictionary
                for name_query, smarts_query in smarts_dict.items():
                    print(name_query, smarts_query)


                    # Create a partial function with additional arguments pre-set
                    partial_func = functools.partial(count_substructure, smarts=smarts_query)

                    # Create a pool of workers. The number of workers is usually set to the number of CPUs
                    # but you can adjust this number based on your specific situation.
                    pool = Pool()

                    # Use the pool's map function to apply your_function to each unique smiles string.
                    # We'll use tqdm to show a progress bar.
                    results = list(tqdm(pool.imap(partial_func, unique_smiles), total=len(unique_smiles)))

                    # create a dictionary mapping from the unique 'smiles' values to their function outputs
                    mapping_dict = dict(zip(unique_smiles, results))

                    # create a new column by mapping the 'smiles' column using the dictionary
                    df['has_' + name_query] = df['smiles_harmonized'].map(mapping_dict)

                    # Make sure to close the pool when you're done to free up system resources.
                    pool.close()
                    pool.join()

    df.to_csv(output_csv_path, index=False)
