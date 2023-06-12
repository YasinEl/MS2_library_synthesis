library(optparse)
library(data.table)


option_list = list(
  make_option(c("-i", "--input_csv"), type="character", default=NULL,
              help="This should be the output csv from the previous python function.", metavar="character"),
  make_option(c("-o", "--output_csv"), type="character", default=NULL,
              help="Path to output csv [default= %default]", metavar="character")
)


parser = OptionParser(option_list=option_list,
                      description = "This summarizes the output by smiles_harmonized.")
arguments = parse_args(parser)

# Extract the input and output file paths
input_csv <- arguments$input_csv
output_csv <- arguments$output_csv





dt_tmp = fread(input_csv)

dt_tmp = dt_tmp[!is.na(smiles_harmonized) & smiles_harmonized != '', c("data_source_is", "dbid_Broad", 
                                                                       "name", "inchikey", "csv_file_name", "dbid_chembl", "smiles", "dbid_CyanoMetDB", 
                                                                       "mf", "monomass", "dbid_drugCentral", "dbid_HMDB", "dbid_PubChem", "dbid_PubChemLite", 
                                                                       "dbid_NORMAN", "dbid_T3DB", "smiles_harmonized", "has_primaryAmine", "has_secondaryAmine", 
                                                                       "has_carboxylicAcid", "has_alcohol", "has_phenol", "has_alkene"
)]

dt_tmp_2 = dt_tmp[, .(data_source_is = {cat("progress",.GRP/.NGRP*100,"%\n"); paste0(unique(data_source_is), collapse = ", ")},
                      #name = as.character(as.data.frame(table(name))$'Var1'[which.max(as.data.frame(table(name))$'Freq')]), 
                      synonyms = paste0(unique(name), collapse = ", "),
                      mf = unique(mf[!is.na(mf) & mf != ''])[1],
                      monomass = unique(monomass[!is.na(monomass) & monomass != ''])[1],
                      inchikey = unique(inchikey[!is.na(inchikey) & inchikey != ''])[1],
                      csv_file_name = paste0(unique(csv_file_name), collapse = ', '),
                      dbid_chembl = paste0(unique(dbid_chembl[!is.na(dbid_chembl) & dbid_chembl != '']), collapse = ', '),
                      dbid_Broad = paste0(unique(dbid_Broad[!is.na(dbid_Broad) & dbid_Broad != '']), collapse = ', '),
                      dbid_CyanoMetDB = paste0(unique(dbid_CyanoMetDB[!is.na(dbid_CyanoMetDB) & dbid_CyanoMetDB != '']), collapse = ', '),
                      dbid_NORMAN = paste0(unique(dbid_NORMAN[!is.na(dbid_NORMAN) & dbid_NORMAN != '']), collapse = ', '),
                      dbid_T3DB = paste0(unique(dbid_T3DB[!is.na(dbid_T3DB) & dbid_T3DB != '']), collapse = ', '),
                      dbid_drugCentral = paste0(unique(dbid_drugCentral[!is.na(dbid_drugCentral) & dbid_drugCentral != '']), collapse = ', '),
                      dbid_HMDB = paste0(unique(dbid_HMDB[!is.na(dbid_HMDB) & dbid_HMDB != '']), collapse = ', '),
                      dbid_PubChem = paste0(unique(dbid_PubChem[!is.na(dbid_PubChem) & dbid_PubChem != '']), collapse = ', '),
                      dbid_PubChemLite = paste0(unique(dbid_PubChemLite[!is.na(dbid_PubChemLite) & dbid_PubChemLite != '']), collapse = ', ')
), by = .(smiles_harmonized, has_primaryAmine, has_secondaryAmine, 
          has_carboxylicAcid, has_alcohol, has_phenol, has_alkene)]


fwrite(dt_tmp_2, output_csv)


