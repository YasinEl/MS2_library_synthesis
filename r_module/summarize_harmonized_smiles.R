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


#read data.table
dt_tmp = fread(input_csv)

#retain columns and remove empty harmonized_smiles
all_cols = colnames(dt_tmp)
general_cols = c('smiles_harmonized', 'data_source_is', 'name', 'inchikey', 'mf' ,'monomass', 'csv_file_name')
id_cols = all_cols[grepl('dbid_', all_cols)]
structure_cols = all_cols[grepl('has_', all_cols)]
use_cols = c(general_cols, structure_cols, id_cols)
dt_tmp = dt_tmp[!is.na(smiles_harmonized) & smiles_harmonized != '', ..use_cols]


# Filter column names that contain 'dbid_' or 'has_'
dbid_cols <- all_cols[grep('^dbid_', all_cols)]
by_cols <- c('smiles_harmonized', all_cols[grep('^has_', all_cols)])

# Initialize string to store expressions
expr_str <- "list(
  data_source_is = {
    cat('progress', .GRP/.NGRP*100, '%\n'); 
    paste0(unique(data_source_is), collapse = ', ')
  },
  synonyms = paste0(unique(name), collapse = ', '),
  mf = unique(mf[!is.na(mf) & mf != ''])[1],
  monomass = unique(monomass[!is.na(monomass) & monomass != ''])[1],
  inchikey = unique(inchikey[!is.na(inchikey) & inchikey != ''])[1],
  csv_file_name = paste0(unique(csv_file_name), collapse = ', '),
  zinc20 = any(zinc20 == 'TRUE'),
  `zinc-instock` = any(`zinc-instock` == 'TRUE'),
  surechembl = any(surechembl == 'TRUE'),
  commerical = any(c(zinc20 == 'TRUE', `zinc-instock` == 'TRUE', surechembl == 'TRUE'))
"

# Add dbid column expressions to the string to summarize all IDs
for(col in dbid_cols){
  expr_str <- paste0(expr_str, ", ", col, " = paste0(unique(", col, "[!is.na(", col, ") & ", col, " != '']), collapse = ', ')")
}

# Assemble code elements
expr_str <- paste0(expr_str, ")")
by_cols <- paste0("list(", paste(by_cols, collapse = ", "), ")")

#Summarize results table
dt_tmp_2 = dt_tmp[, eval(parse(text = expr_str)), by = eval(parse(text = by_cols))]

#write output
fwrite(dt_tmp_2, output_csv)


