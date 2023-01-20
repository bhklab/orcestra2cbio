library(data.table)

format_assay_data <- function(work_dir, obj_name, filename, assay, gene_metadata){
  assay <- assay[rownames(assay)[rownames(assay) %in% rownames(gene_metadata)], ]
  
  assay$Hugo_Symbol <- unlist(lapply(rownames(assay), function(gene){
    filtered <- gene_metadata[gene_metadata$gene_id == gene, ]
    return(filtered$hugo_symbol)
  }))
  
  assay$Entrez_Gene_Id <- unlist(lapply(rownames(assay), function(gene){
    filtered <- gene_metadata[rownames(gene_metadata) == gene, ]
    return(filtered$entrez_id)
  }))
  
  assay <- assay[, c('Hugo_Symbol', 'Entrez_Gene_Id', colnames(assay)[!colnames(assay) %in% c('Hugo_Symbol', 'Entrez_Gene_Id')])]
  
  write.table(
    assay,
    file=file.path(work_dir, obj_name, filename),
    quote=FALSE,
    row.names = FALSE,
    sep='\t'
  )
}