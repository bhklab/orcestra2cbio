options(timeout=1000)

library(MultiAssayExperiment)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(biomaRt)
library(stringr)
library(data.table)
library(jsonlite)

# Root director of the scripts
work_dir <- '~/Code/bhklab/orcestra2cbio'

# Directory to download the data objects
input_dir <- '~/Documents/bhklab/orcestra2cbio/data_objects'

# Directory to save the files to upload to cBioPortal
output_dir <- '~/Documents/bhklab/orcestra2cbio/results'

source(file.path(work_dir, 'common.R'))
source(file.path(work_dir, 'clinical_patient_data.R'))
source(file.path(work_dir, 'process_assay_data.R'))
source(file.path(work_dir, 'format_snv_data.R'))

# snv_files <- list.files('~/Documents/bhklab/orcestra2cbio/SNV')
# for(file in snv_files){
#   obj_name <- str_replace(file, '_SNV.csv', '')
#   
#   snv <- read.csv(file.path('~/Documents/bhklab/orcestra2cbio/SNV', file), sep=';')
#   genome_ref <- 'GRCh38'
#   if(obj_name %in% c('ICB_VanDenEnde', 'ICB_Miao1', 'ICB_Snyder', 'ICB_Van_Allen', 'ICB_Puch', 'ICB_Jerby_Arnon')){
#     genome_ref <- 'GRCh37'
#   }
#   
#   tsv_table <- format_snv_data(snv, genome_ref)
#   write.table(
#     tsv_table, 
#     file=file.path('~/Documents/bhklab/orcestra2cbio/SNV', paste0(obj_name, '_SNV.tsv')), 
#     sep='\t', 
#     quote=FALSE, 
#     row.names = FALSE
#   )
# }

# amino acid changes available: ICB_Braun, ICB_Mariathasan, ICB_Miao1, ICB_Nathanson, ICB_Riaz, ICB_Roh, ICB_Van_Allen, ICB_Jung
# missing amino acid changes: ICB_Liu 


dir.create(input_dir)
dir.create(output_dir)

# Download all available ICB data objects from ORCESTRA.
icb_objects <- fromJSON(txt='https://www.orcestra.ca/api/clinical_icb/canonical')
for(obj_name in icb_objects$name){
  link <- icb_objects$downloadLink[icb_objects$name == obj_name]
  filename <- paste0(obj_name, '.rds')
  if(!file.exists(file.path(input_dir, filename))){
    download.file(url=link, destfile=file.path(input_dir, filename))
  }
}

# Construct a table to indicate which data type is present in each study (expr, cna or snv)
obj_metadata <- data.frame(matrix(nrow=length(icb_objects$name), ncol=5))
colnames(obj_metadata) <- c('name', 'expr', 'cna', 'snv', 'description')  
obj_metadata$name <- icb_objects$name
obj_metadata$description <- paste('Description available at', paste0('https://www.orcestra.ca/', icb_objects$doi))
expr <- c()
cna <- c()
snv <- c()

for(object in icb_objects$name){
  obj <- readRDS(file.path(input_dir, paste0(object, '.rds')))
  assays <- names(experiments(obj))
  if('expr' %in% assays | 'expr_gene_tpm' %in% assays){
    expr <- c(expr, 1) 
  }else{
    expr <- c(expr, NA)
  }
  if('cna' %in% assays){
    cna <- c(cna, 1) 
  }else{
    cna <- c(cna, NA)
  }
  if('snv' %in% assays){
    snv <- c(snv, 1) 
  }else{
    snv <- c(snv, NA)
  }
}
obj_metadata$expr <- expr
obj_metadata$cna <- cna
obj_metadata$snv <- snv

# Iterate though the studies to generate the files.
for(obj_name in icb_objects$name){
  print(obj_name)
  # obj_name <- 'ICB_Riaz'
  
  dir.create(file.path(output_dir, obj_name))
  dir.create(file.path(output_dir, obj_name, 'case_lists'))
  
  # Process Clinical data
  print('processing clinical data')
  obj <- readRDS(file.path(input_dir, paste0(obj_name, '.rds')))
  clin <- data.frame(colData(obj))
  
  format_clinical_patient_data(output_dir, obj_name, clin)
  
  meta_clin <- list(
    cancer_study_identifier=str_to_lower(obj_name),
    genetic_alteration_type='CLINICAL',
    datatype='PATIENT_ATTRIBUTES',
    data_filename='data_clinical_patient.txt'
  )
  format_meta_data(file.path(output_dir, obj_name, 'meta_clinical_patient.txt'), meta_clin)
  
  cases_all <- list(
    cancer_study_identifier=str_to_lower(obj_name),
    stable_id=paste0(str_to_lower(obj_name), '_all'),
    case_list_name='All samples',
    case_list_description=paste0('All samples (', length(clin$patientid), ')'),
    case_list_category='all_cases_in_study',
    case_list_ids=paste(clin$patientid, collapse='\t')
  )
  format_meta_data(file.path(output_dir, obj_name, 'case_lists', 'cases_all.txt'), cases_all)
  
  obj_meta <- obj_metadata[obj_metadata$name == obj_name, ]
  
  # Process RNA-seq data
  if(!is.na(obj_meta$expr)){
    print('processing expr data')
    assay <- 'expr'
    assays <- names(experiments(obj))
    if('expr_gene_tpm' %in% assays){
      assay <- 'expr_gene_tpm'
    }
    
    expr <- data.frame(assays(obj)[[assay]])
    expr_gene_metadata <- data.frame(rowData(experiments(obj)[[assay]]))
    expr_gene_metadata <- get_gene_metadata(rownames(expr), expr_gene_metadata)
    
    format_assay_data(output_dir, obj_name, 'data_rna_seq_mrna.txt', expr, expr_gene_metadata)
    
    meta_expr <- list(
      cancer_study_identifier=str_to_lower(obj_name),
      stable_id='rna_seq_mrna',
      genetic_alteration_type='MRNA_EXPRESSION',
      datatype='CONTINUOUS',
      show_profile_in_analysis_tab='FALSE',
      profile_name='mRNA Expression',
      profile_description='mRNA Expression',
      data_filename='rna_seq_mrna.txt'
    )
    format_meta_data(file.path(output_dir, obj_name, 'meta_rna_seq_mrna.txt'), meta_expr) 
    
    cases_expr <- list(
      cancer_study_identifier=str_to_lower(obj_name),
      stable_id=paste0(str_to_lower(obj_name), '_rna_seq_mrna'),
      case_list_name='Samples with mRNA (RNA-seq) data',
      case_list_description=paste0('Samples with mRNA (RNA-seq) data (', length(colnames(expr)), ' samples)'),
      case_list_category='all_cases_with_mrna_rnaseq_data',
      case_list_ids=paste(colnames(expr), collapse='\t')
    )
    format_meta_data(file.path(output_dir, obj_name, 'case_lists', 'cases_rna_seq_mrna.txt'), cases_expr)
  }
  
  # Process CNA data
  if(!is.na(obj_meta$cna)){
    print('processing cna data')
    cna <- data.frame(assays(obj)[['cna']])
    cna_gene_metadata <- data.frame(rowData(experiments(obj)[['cna']]))
    cna_gene_metadata <- get_gene_metadata(rownames(cna), cna_gene_metadata)
    
    format_assay_data(output_dir, obj_name, 'data_cna.txt', cna, cna_gene_metadata)
    
    meta_cna <- list(
      cancer_study_identifier=str_to_lower(obj_name),
      stable_id='cna',
      genetic_alteration_type='COPY_NUMBER_ALTERATION',
      datatype='DISCRETE',
      show_profile_in_analysis_tab='TRUE',
      profile_name='Putative copy-number alterations from GISTIC',
      profile_description='Putative copy-number from GISTIC 2.0. Values: -2 = homozygous deletion; -1 = hemizygous deletion; 0 = neutral / no change; 1 = gain; 2 = high level amplification.',
      data_filename='data_cna.txt'
    )
    format_meta_data(file.path(output_dir, obj_name, 'meta_cna.txt'), meta_cna)  
    
    cases_cna <- list(
      cancer_study_identifier=str_to_lower(obj_name),
      stable_id=paste0(str_to_lower(obj_name), '_cna'),
      case_list_name='Samples with CNA data',
      case_list_description=paste0('Samples with CNA data (', length(colnames(cna)), ' samples)'),
      case_list_category='all_cases_with_cna_data',
      case_list_ids=paste(colnames(cna), collapse='\t')
    )
    format_meta_data(file.path(output_dir, obj_name, 'case_lists', 'cases_cna.txt'), cases_cna)
  }
  
  # Process SNV data
  if(!is.na(obj_meta$snv)){
    print('processing snv data')
    genome_ref <- 'GRCh38'
    if(obj_name %in% c('ICB_VanDenEnde', 'ICB_Miao1', 'ICB_Snyder', 'ICB_Van_Allen', 'ICB_Puch', 'ICB_Jerby_Arnon')){
      genome_ref <- 'GRCh37'
    }
    
    snv <- read.csv(file.path(work_dir, 'SNV', paste0(obj_name, '_SNV.csv')), sep=';')
    format_snv_data(output_dir, obj_name, 'data_mutations.txt', snv, genome_ref)
    
    meta_snv <- list(
      cancer_study_identifier=str_to_lower(obj_name),
      stable_id='mutations',
      genetic_alteration_type='MUTATION_EXTENDED',
      datatype='MAF',
      show_profile_in_analysis_tab='TRUE',
      profile_name='SNV',
      profile_description='SNV data.',
      data_filename='data_mutations.txt'
    )
    format_meta_data(file.path(output_dir, obj_name, 'meta_mutations.txt'), meta_snv)
    
    patients <- unique(snv$Sample)
    cases_sequenced <- list(
      cancer_study_identifier=str_to_lower(obj_name),
      stable_id=paste0(str_to_lower(obj_name), '_sequenced'),
      case_list_name='Samples with mutation data',
      case_list_description=paste0('Samples with mutation data (', length(patients), ' samples)'),
      case_list_category='all_cases_with_mutation_data',
      case_list_ids=paste(patients, collapse='\t')
    )
    format_meta_data(file.path(output_dir, obj_name, 'case_lists', 'cases_sequenced.txt'), cases_sequenced)
  }
  
  # Create study metadata
  meta_study <- list(
    type_of_cancer='mixed',
    cancer_study_identifier=str_to_lower(obj_name),
    name=obj_name,
    description=obj_meta$description,
    groups='PUBLIC;PANCAN'
  )
  format_meta_data(file.path(output_dir, obj_name, 'meta_study.txt'), meta_study)
}
