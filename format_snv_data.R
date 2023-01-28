library(stringr)
library(data.table)

format_snv_data <- function(work_dir, obj_name, filename, snv, genome_ref){
  colnames(snv) <- c(
    'Tumor_Sample_Barcode', 
    'Hugo_Symbol', 
    'Chromosome', 
    'Start_Position', 
    'Reference_Allele', 
    'Tumor_Seq_Allele2', 
    'Variant_Classification', 
    'Variant_Type',
    'HGVSp_Short'
  )
  
  snv$End_Position <- unlist(lapply(seq_along(rownames(snv)), function(index){
    len <- str_length(snv$Reference_Allele[index]) - 1
    return(snv$Start_Position[index] + len)
  }))
  
  snv$Chromosome <- str_replace(snv$Chromosome, 'chr', '')
  
  snv$NCBI_Build <- genome_ref
  
  snv <- snv[, c(
    'Tumor_Sample_Barcode', 
    'Hugo_Symbol', 
    'Chromosome', 
    'Start_Position',
    'End_Position',
    'Reference_Allele', 
    'Tumor_Seq_Allele2', 
    'Variant_Classification', 
    'Variant_Type',
    'HGVSp_Short',
    'NCBI_Build'
  )]
  
  write.table(
    snv,
    file=file.path(work_dir, obj_name, filename),
    quote=FALSE,
    row.names = FALSE,
    sep='\t',
    na=""
  )
}