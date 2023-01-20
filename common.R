library(biomaRt)
library(stringr)

format_meta_data <- function(filepath, meta){
  contents <- unlist(lapply(names(meta), function(name){
    return(paste0(name, ': ', meta[[name]], '\n'))
  }))
  cat(contents, sep="", file=filepath)
}

get_gene_metadata <- function(genes, gene_metadata){
  
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  filters = listFilters(ensembl)
  
  genes_no_ver <- gene_metadata$gene_id_no_ver
  
  cols <- c(
    'ensembl_gene_id',
    'ensembl_gene_id_version',
    'hgnc_symbol',
    'entrezgene_id'
  )
  filtered <- select(
    ensembl,
    keys=genes_no_ver,
    columns=cols,
    keytype='ensembl_gene_id'
  )
  
  missing <- genes_no_ver[!genes_no_ver %in% filtered$ensembl_gene_id]
  filtered_gene_metadata <- gene_metadata[gene_metadata$gene_id_no_ver[!gene_metadata$gene_id_no_ver %in% missing], ]
  filtered_gene_metadata <- filtered_gene_metadata[rownames(filtered_gene_metadata)[str_detect(rownames(filtered_gene_metadata), 'ENSG')], ]
  
  filtered_gene_metadata$hugo_symbol <- unlist(lapply(filtered_gene_metadata$gene_id_no_ver, function(gene){
    found <- filtered[filtered$ensembl_gene_id == gene, c('hgnc_symbol')]
    return(found[1])
  }))
  
  filtered_gene_metadata$entrez_id <- unlist(lapply(filtered_gene_metadata$gene_id_no_ver, function(gene){
    found <- filtered[filtered$ensembl_gene_id == gene, c('entrezgene_id')]
    return(found[1])
  }))
  
  print(paste(length(missing), 'genes omitted'))
  print(paste(length(rownames(filtered_gene_metadata)), 'genes included'))
  
  return(filtered_gene_metadata)
}