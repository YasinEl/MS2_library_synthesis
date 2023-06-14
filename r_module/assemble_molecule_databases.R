library(optparse)
library(data.table)

option_list = list(
  make_option(c("-i", "--input_folder"), type="character", default=NULL,
              help="In which folder are csv or tsv files to be used?", metavar="character"),
  make_option(c("-o", "--output_csv"), type="character", default=NULL,
              help="Path to output csv [default= %default]", metavar="character"),
  make_option(c("-n", "--nomenclature_csv"), type="character", default=NULL,
              help="Path to the database nomenclature csv [default= %default]", metavar="character")
)


parser = OptionParser(option_list=option_list,
                      description = "This allows to assembl molecules from different databases. The different databases allowed are currently hardcoded.
                      The same is true for the csv file names. Thiw will remove most columns and keep only a few.")
arguments = parse_args(parser)


# Extract the input and output file paths
input_folder <- arguments$input_folder
output_file <- arguments$output_csv
nomenclature_file <- arguments$nomenclature_csv


databases_csvs = list.files(input_folder, full.names = TRUE, pattern = '\\.(csv|tsv)$')

dt_nomenclature = fread(nomenclature_file)

csv_triggers = colnames(dt_nomenclature)[-1]

li_collection = list()
length(li_collection) = length(databases_csvs)

for(i in seq(length(databases_csvs))){
  databases_csv = databases_csvs[i]
  
  print(paste0(c(basename(databases_csv), ' (', i, '/', length(databases_csvs), ')'), collapse = ''))
  
  trigger_hit = csv_triggers[sapply(csv_triggers , function(x) grepl(x, basename(databases_csv)))]
  
  if(length(trigger_hit) > 0){
    dt_tmp = fread(databases_csv)
    
    #put the database name
    dt_tmp[, data_source_is := suppressWarnings(unname(dt_nomenclature[Info_by_database == 'database', ..trigger_hit]))]
    
    #put the csv file name
    dt_tmp[, csv_file_name := basename(databases_csv)]
    
    #add the ID from the database
    db_id = paste0(c('dbid_', suppressWarnings(unname(dt_nomenclature[Info_by_database == 'database', ..trigger_hit]))), collapse = '')
    dbid_col = suppressWarnings(unlist(unname(dt_nomenclature[Info_by_database == 'database id', ..trigger_hit])))
    setnames(dt_tmp, dbid_col, db_id)
    
    #set the columns that should be kept
    db_cols =  suppressWarnings(unlist(unname(dt_nomenclature[!(Info_by_database %in% c('database','database id')), ..trigger_hit])))
    db_cols = db_cols[!is.na(db_cols) & db_cols != '']
    cols_to_keep = c('data_source_is', db_id, 'csv_file_name', db_cols)
    dt_tmp = dt_tmp[, ..cols_to_keep]
    
    #harmonize column names
    setnames(dt_tmp, 
             db_cols,
             dt_nomenclature[!(Info_by_database %in% c('database','database id')) & !is.na(get(trigger_hit)) & get(trigger_hit) != '']$Info_by_database)
    } else {
    warning('Ignore ', basename(databases_csv))
  }

  li_collection[[i]] = dt_tmp
  
}
#assemble everything into one table
dt_overall_table = rbindlist(li_collection, fill = TRUE, use.names = TRUE)


write.csv(dt_overall_table, output_file, fileEncoding = "UTF-8")