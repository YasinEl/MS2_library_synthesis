library(optparse)
library(data.table)


option_list = list(
  make_option(c("-i", "--input_folder"), type="character", default=NULL,
              help="In which folder are csv or tsv files to be used?", metavar="character"),
  make_option(c("-o", "--output_csv"), type="character", default=NULL,
              help="Path to output csv [default= %default]", metavar="character")
)


parser = OptionParser(option_list=option_list,
                      description = "This allows to assembl molecules from different databases. The different databases allowed are currently hardcoded.
                      The same is true for the csv file names. Thiw will remove most columns and keep only a few.")
arguments = parse_args(parser)

# Extract the input and output file paths
input_folder <- arguments$input_folder
output_file <- arguments$output_csv



databases_csvs = list.files(input_folder, full.names = TRUE, pattern = '\\.(csv|tsv)$')


li_collection = list()
length(li_collection) = length(databases_csvs)

for(i in seq(length(databases_csvs))){
  databases_csv = databases_csvs[i]
  
  print(paste0(c(basename(databases_csv), ' (', i, '/', length(databases_csvs), ')'), collapse = ''))
  
  
  if(grepl('PubChem_compound_list', basename(databases_csv))){
    
    dt_tmp = fread(databases_csv)
    dt_tmp[, data_source_is := 'PubChem']
    dt_tmp[, csv_file_name := basename(databases_csv)]
    dt_tmp[, dbid_PubChem := cid]
    dt_tmp = dt_tmp[, c('data_source_is', 'dbid_PubChem', 'cmpdname', 'mf', 'exactmass', 'canonicalsmiles', 'inchikey', 'csv_file_name')]
    setnames(dt_tmp, 'cmpdname', 'name')
    setnames(dt_tmp, 'exactmass', 'monomass')
    setnames(dt_tmp, 'canonicalsmiles', 'smiles')
    
  } else if(grepl('HMDB', basename(databases_csv))){
    dt_tmp = fread(databases_csv)
    dt_tmp[, data_source_is := 'HMDB']
    dt_tmp[, csv_file_name := basename(databases_csv)]
    dt_tmp[, dbid_HMDB := HMDB_ID]
    dt_tmp = dt_tmp[, c('data_source_is', 'dbid_HMDB', 'NAME', 'CHEMICAL_FORMULA', 'MONO_MASS', 'SMILES', 'INCHIKEY', 'csv_file_name')]
    setnames(dt_tmp, 'NAME', 'name')
    setnames(dt_tmp, 'CHEMICAL_FORMULA', 'mf')
    setnames(dt_tmp, 'MONO_MASS', 'monomass')
    setnames(dt_tmp, 'SMILES', 'smiles')
    setnames(dt_tmp, 'INCHIKEY', 'inchikey')
    
  } else if(grepl('CyanoMetDB', basename(databases_csv))){
    dt_tmp = fread(databases_csv)
    dt_tmp[, data_source_is := 'CyanoMetDB']
    dt_tmp[, csv_file_name := basename(databases_csv)]
    dt_tmp[, dbid_CyanoMetDB := `Compound identifier`]
    dt_tmp = dt_tmp[, c('data_source_is', 'dbid_CyanoMetDB', 'Compound name', 'Molecular formula (neutral)', 'Monoisotopic mass (neutral)', 'SMILES (canonical or isomeric)', 'InChlKey', 'csv_file_name')]
    setnames(dt_tmp, 'Compound name', 'name')
    setnames(dt_tmp, 'Molecular formula (neutral)', 'mf')
    setnames(dt_tmp, 'Monoisotopic mass (neutral)', 'monomass')
    setnames(dt_tmp, 'SMILES (canonical or isomeric)', 'smiles')
    setnames(dt_tmp, 'InChlKey', 'inchikey')
    
  } else if(grepl('PubChemLite', basename(databases_csv))){
    dt_tmp = fread(databases_csv)
    dt_tmp[, data_source_is := 'PubChem Lite']
    dt_tmp[, csv_file_name := basename(databases_csv)]
    dt_tmp[, dbid_PubChemLite := Identifier]
    dt_tmp = dt_tmp[, c('data_source_is', 'dbid_PubChemLite', 'Synonym', 'MolecularFormula', 'MonoisotopicMass', 'SMILES', 'InChIKey', 'csv_file_name')]
    setnames(dt_tmp, 'Synonym', 'name')
    setnames(dt_tmp, 'MolecularFormula', 'mf')
    setnames(dt_tmp, 'MonoisotopicMass', 'monomass')
    setnames(dt_tmp, 'SMILES', 'smiles')
    setnames(dt_tmp, 'InChIKey', 'inchikey')
    
  } else if(grepl('susdat_', basename(databases_csv))){
    dt_tmp = fread(databases_csv)
    dt_tmp[, data_source_is := 'NORMAN']
    dt_tmp[, csv_file_name := basename(databases_csv)]
    dt_tmp[, dbid_NORMAN := Norman_SusDat_ID]
    dt_tmp = dt_tmp[, c('data_source_is', 'dbid_NORMAN', 'Name', 'Molecular_Formula', 'Monoiso_Mass', 'MS_Ready_SMILES', 'MS_Ready_StdInChIKey', 'csv_file_name')]
    setnames(dt_tmp, 'Name', 'name')
    setnames(dt_tmp, 'Molecular_Formula', 'mf')
    setnames(dt_tmp, 'Monoiso_Mass', 'monomass')
    setnames(dt_tmp, 'MS_Ready_SMILES', 'smiles')
    setnames(dt_tmp, 'MS_Ready_StdInChIKey', 'inchikey')
    
  }  else if(grepl('toxins.csv', basename(databases_csv), fixed = TRUE)){
    dt_tmp = fread(databases_csv)
    dt_tmp[, data_source_is := 'T3DB']
    dt_tmp[, csv_file_name := basename(databases_csv)]
    dt_tmp[, dbid_T3DB := `T3DB ID`]
    dt_tmp = dt_tmp[, c('data_source_is', 'dbid_T3DB', 'Name', 'Chemical Formula', 'Monoisotopic Mass', 'SMILES', 'InChI Key', 'csv_file_name')]
    setnames(dt_tmp, 'Name', 'name')
    setnames(dt_tmp, 'Chemical Formula', 'mf')
    setnames(dt_tmp, 'Monoisotopic Mass', 'monomass')
    setnames(dt_tmp, 'SMILES', 'smiles')
    setnames(dt_tmp, 'InChI Key', 'inchikey')
    
  }  else if(grepl('drugCentral', basename(databases_csv))){
    dt_tmp = fread(databases_csv)
    dt_tmp[, data_source_is := 'drugCentral']
    dt_tmp[, csv_file_name := basename(databases_csv)]
    dt_tmp[, dbid_drugCentral := id]
    dt_tmp = dt_tmp[, c('data_source_is', 'dbid_drugCentral', 'name', 'smiles', 'inchikey', 'csv_file_name')]
    setnames(dt_tmp, 'name', 'name')
    setnames(dt_tmp, 'smiles', 'smiles')
    setnames(dt_tmp, 'inchikey', 'inchikey')
    
  }  else if(grepl('Broad_Institute', basename(databases_csv))){
    dt_tmp = fread(databases_csv)
    dt_tmp[, data_source_is := 'Broad Institute Drug List']
    dt_tmp[, csv_file_name := basename(databases_csv)]
    dt_tmp[, dbid_Broad := Broad_ID]
    dt_tmp = dt_tmp[, c('data_source_is', 'dbid_Broad', 'pert_iname', 'InChIKey', 'csv_file_name')]
    setnames(dt_tmp, 'pert_iname', 'name')
    setnames(dt_tmp, 'InChIKey', 'inchikey')
    
  }  else if(grepl('chembl_final', basename(databases_csv))){
    dt_tmp = fread(databases_csv)
    dt_tmp[, data_source_is := 'chembl']
    dt_tmp[, csv_file_name := basename(databases_csv)]
    dt_tmp[, dbid_chembl := chembl_id]
    dt_tmp = dt_tmp[, c('data_source_is', 'dbid_chembl', 'compound_name', 'canonical_smiles', 'csv_file_name')]
    setnames(dt_tmp, 'compound_name', 'name')
    setnames(dt_tmp, 'canonical_smiles', 'smiles')
    
  }  else if(grepl('drugbank_final', basename(databases_csv))){
    dt_tmp = fread(databases_csv)
    dt_tmp[, data_source_is := 'drugbank']
    dt_tmp[, csv_file_name := basename(databases_csv)]
    dt_tmp[, dbid_drugbank := drugbank_id]
    dt_tmp = dt_tmp[, c('data_source_is', 'dbid_drugbank', 'name', 'smiles', 'inchi_key', 'csv_file_name')]
    setnames(dt_tmp, 'name', 'name')
    setnames(dt_tmp, 'smiles', 'smiles')
    setnames(dt_tmp, 'inchi_key', 'inchikey')
    
  }

  li_collection[[i]] = dt_tmp
  
}

dt_overall_table = rbindlist(li_collection, fill = TRUE, use.names = TRUE)


write.csv(dt_overall_table, output_file, fileEncoding = "UTF-8")