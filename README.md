# MS2 library synthesis
Some functions to summarize different molecular databases.

# Instructions to Run Functions from the Command Line

First, set the path to the project directory. For example:

```bash
cd C:\Users\Dorrestein Lab\MS2_library_synthesis
```

To combine all different databases into one CSV file, you have to have first download them (e.g. from the google drive). You must not change the column names in these csv files, and also not the csv file itself! You can also add more databases if you want as csv or tsv files. If you do you need to extend the file in data/database_nomenclature.csv. Just add a column with the appropriate column names as has been done for the others. The first row should be a substring of the csv file name (make sure it appears ONLY in the csvs you want this substring to match to; can be multiple). The second row should be the name of the database. I hope the rest is intuitive. Afterwards you can assemble all databases like this:
```bash
Rscript r_module\assemble_molecule_databases.R -i Z:\5M_molecules\datasources\prepared -o Z:\5M_molecules\database_accumulation.csv -n data\database_nomenclature.csv
```
This merges different IDs, SMILES, etc., into the same columns and removes unused columns.

To harmonize SMILES and search structural subgroups in the SMILES, as well as calculate molecular formula and masses from SMILES, run:


```bash
python python_module\main.py --input_csv Z:/5M_molecules/database_accumulation.csv --output_csv Z:/5M_molecules/combined_datasources_yasin_update_command.csv --intermedResults Z:/5M_molecules --substructure_csv Z:/5M_molecules/smarts_lib.csv --TautomerLimit 500
```
If you want to test the function on less than 1% of your data, you can set the --test flag to True:


```bash
python python_module\main.py --input_csv Z:/5M_molecules/database_accumulation.csv --output_csv Z:/5M_molecules/combined_datasources_yasin_update_command.csv --intermedResults Z:/5M_molecules --substructure_csv Z:/5M_molecules/smarts_lib.csv --TautomerLimit 500 --test True
```
If you do not want to repeat the SMILES harmonization and want to perform the substructure search (a smiles_harmonized column must be present in the input), you can set --OnlySubgroups to True:


```bash
python python_module\main.py --input_csv Z:/5M_molecules/database_accumulation.csv --output_csv Z:/5M_molecules/combined_datasources_yasin_update_command.csv --intermedResults Z:/5M_molecules --substructure_csv Z:/5M_molecules/smarts_lib.csv --OnlySubgroups True
```
Finally, to summarize everything along the harmonized SMILES (the final table will have one row per harmonized SMILES), run:


```bash
Rscript r_module\summarize_harmonized_smiles.R -i Z:/5M_molecules/combined_datasources_yasin_update_command.csv -o Z:/5M_molecules/combined_datasources_yasin_summary_update_command.csv

```
