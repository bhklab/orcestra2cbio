library(stringr)
library(data.table)
library(seqinr)

work_dir <- '~/Documents/bhklab/orcestra2cbio/SNV'

# ICB_Braun
braun <- readRDS(file.path(work_dir, 'data_objects/ICB_Braun.rds'))
braun_clin <- data.frame(colData(braun))
braun_snv_original <- data.frame( fread( file.path(work_dir, 'ICB_Braun_SNV.txt.gz') , sep="\t" , stringsAsFactors=FALSE  ))
braun_snv <- read.csv(file.path(work_dir, 'ICB_Braun_SNV.csv'), sep=';')

braun_clin <- braun_clin[, c('patientid', 'MAF_Tumor_ID')]
braun_snv$Chr <- str_replace(braun_snv$Chr, 'chr', '')
braun_snv$Sample <- str_replace_all(braun_snv$Sample, '-', '_')

unique_id <- unique(braun_clin$MAF_Tumor_ID)
unique_id <- unique_id[!is.na(unique_id)]
braun_snv_original <- braun_snv_original[braun_snv_original$Tumor_Sample_Barcode %in% unique_id, ]
patient_id <- unlist(lapply(braun_snv_original$Tumor_Sample_Barcode, function(id){
  filtered <- braun_clin[braun_clin$MAF_Tumor_ID %in% id, ]
  return(filtered$patientid)
}))
braun_snv_original$patientid <- patient_id

protein_changes <- c()
for(i in 1:nrow(braun_snv)){
  found <- braun_snv_original[
    braun_snv_original$patientid %in% braun_snv[i, 1] & braun_snv_original$Hugo_Symbol %in% braun_snv[i, 2] & braun_snv_original$Chromosome %in% braun_snv[i, 3] & braun_snv_original$Start_position %in% braun_snv[i, 4], 
    c('Protein_Change')]
  if(length(found) > 1){
    found <- found[1]
  }
  protein_changes <- append(protein_changes, found)
}
braun_snv$HGVSp_Short <- protein_changes
write.table( braun_snv , file=file.path(work_dir, "ICB_Braun_SNV.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )

# ICB_Jung
jung_snv_original <- data.frame( fread( file.path(work_dir, 'ICB_Jung_SNV.csv.gz') , sep=";" , stringsAsFactors=FALSE  ))
jung_snv <- read.csv(file.path(work_dir, 'ICB_Jung_SNV.csv'), sep=';')

protein_changes <- str_replace_all(jung_snv_original$Gene, '\\(', ':')
protein_changes <- str_replace_all(protein_changes, '\\)', ':')
protein_changes <- unlist(lapply(protein_changes, function(change){
  str <- str_split(change, ',')[[1]][1]
  return(str)
}))
protein_changes <- unlist(lapply(protein_changes, function(change){
  val <- NA
  if(str_detect(change, 'p\\.')){
    val <- str_extract(change, 'p\\..*')
  }
  return(val)
}))

jung_snv_original$HGVSp_Short <- protein_changes
jung_snv_original$Sample.ID <- paste0('P.', jung_snv_original$Sample.ID)

protein_changes <- c()
for(i in 1:nrow(jung_snv)){
  found <- jung_snv_original[
    jung_snv_original$Sample.ID %in% jung_snv[i, 1] & jung_snv_original$Chromosome %in% jung_snv[i, 3] & jung_snv_original$Start %in% jung_snv[i, 4], 
    c('HGVSp_Short')]
  if(length(found) > 1){
    found <- found[1]
  }
  protein_changes <- append(protein_changes, found)
}
jung_snv$HGVSp_Short <- protein_changes
write.table( jung_snv , file=file.path(work_dir, "ICB_Jung_SNV.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )

# ICB_Liu
liu_snv <- read.csv(file.path(work_dir, 'ICB_Liu_SNV.csv'), sep=';')
liu_snv$HGVSp_Short <- NA
write.table( liu_snv , file=file.path(work_dir, "ICB_Liu_SNV.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )

# ICB_Mariathasan
mariathasan_snv_original <- data.frame( fread( file.path(work_dir, 'ICB_Mariathasan_SNV.txt.gz') , sep=";" , stringsAsFactors=FALSE  ))
mariathasan_snv <- read.csv(file.path(work_dir, 'ICB_Mariathasan_SNV.csv'), sep=';')
mariathasan_snv_original <- mariathasan_snv_original[mariathasan_snv_original$patient %in% unique(mariathasan_snv$Sample), ]

protein_changes <- unlist(lapply(mariathasan_snv_original$mutation, function(change){
  changes <- unlist(str_split(change, '_'))
  protein <- changes[str_detect(changes, 'p\\.[A-Z]')]
  if(length(protein) > 0){
    return(protein)
  }
  return(NA)
}))
mariathasan_snv_original$HGVSp_Short <- protein_changes

protein_changes <- c()
for(i in 1:nrow(mariathasan_snv)){
  found <- mariathasan_snv_original[
    mariathasan_snv_original$patient %in% mariathasan_snv[i, 1] & mariathasan_snv_original$gene %in% mariathasan_snv[i, 2],
    c('HGVSp_Short')]
  if(length(found) > 1){
    found <- found[1]
  }
  protein_changes <- append(protein_changes, found)
}
mariathasan_snv$HGVSp_Short <- protein_changes
write.table( mariathasan_snv , file=file.path(work_dir, "ICB_Mariathasan_SNV.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )

# ICB_Miao1
miao1_snv_original <- data.frame( fread( file.path(work_dir, 'ICB_Miao1_SNV.txt.gz') , sep="\t" , stringsAsFactors=FALSE  ))
miao1_snv <- read.csv(file.path(work_dir, 'ICB_Miao1_SNV.csv'), sep=';')
miao1_snv$Chr <- str_replace(miao1_snv$Chr, 'chr', '')
miao1_snv_original$patientid <- str_replace(miao1_snv_original$Tumor_Sample_Barcode, '-.*', '')

protein_changes <- c()
for(i in 1:nrow(miao1_snv)){
  found <- miao1_snv_original[
    miao1_snv_original$patientid %in% miao1_snv[i, 1] & miao1_snv_original$Hugo_Symbol %in% miao1_snv[i, 2] & miao1_snv_original$Chromosome %in% miao1_snv[i, 3] & miao1_snv_original$Start_position %in% miao1_snv[i, 4], 
    c('Protein_Change')]
  if(length(found) > 1){
    found <- found[1]
  }
  protein_changes <- append(protein_changes, found)
}

miao1_snv$HGVSp_Short <- protein_changes 
write.table( miao1_snv , file=file.path(work_dir, "ICB_Miao1_SNV.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )

# ICB_Nathanson
nathanson_snv_original <- data.frame( fread( file.path(work_dir, 'ICB_Nathanson_SNV.txt.gz') , sep="\t" , stringsAsFactors=FALSE  ))
nathanson_snv <- read.csv(file.path(work_dir, 'ICB_Nathanson_SNV.csv'), sep=';')
nathanson_snv$Chr <- str_replace(nathanson_snv$Chr, 'chr', '')
nathanson_snv_original$sample <- paste0('P', nathanson_snv_original$sample)

protein_changes <- c()
for(i in 1:nrow(nathanson_snv)){
  found <- nathanson_snv_original[
    nathanson_snv_original$sample %in% nathanson_snv[i, 1] & nathanson_snv_original$chr %in% nathanson_snv[i, 3] & nathanson_snv_original$pos %in% nathanson_snv[i, 4], 
    c('effect')]
  if(length(found) > 1){
    found <- found[1]
  }
  if(found == ""){
    found <- NA
  }
  protein_changes <- append(protein_changes, found)
}

nathanson_snv$HGVSp_Short <- protein_changes
write.table( nathanson_snv , file=file.path(work_dir, "ICB_Nathanson_SNV.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )

# ICB_Riaz
riaz_snv_original <- data.frame( fread( file.path(work_dir, 'ICB_Riaz_SNV.txt.gz') , sep=";" , stringsAsFactors=FALSE  ))
riaz_snv <- read.csv(file.path(work_dir, 'ICB_Riaz_SNV.csv'), sep=';')

protein_changes <- c()
for(i in 1:nrow(riaz_snv)){
  found <- riaz_snv_original[
    riaz_snv_original$Patient %in% riaz_snv[i, 1] & riaz_snv_original$Hugo.Symbol %in% riaz_snv[i, 2] & riaz_snv_original$Start %in% riaz_snv[i, 4], 
    c('HGVS_p')]
  if(length(found) > 1){
    found <- found[1]
  }
  if(found == "."){
    found <- NA
  }
  protein_changes <- append(protein_changes, found)
}

protein_changes <- lapply(protein_changes, function(change){
  if(!is.na(change)){
    num <- str_extract(change, '\\d{1,4}')
    str <- unlist(str_split(change, '\\d{1,4}'))
    first <- str_replace(str[1], 'p\\.', '')
    second <- str[2]
    first <- a(first)
    if(is.na(a(second))){
      second <- '*'
    }else{
      second <- a(second)
    }
    change <- paste0('p.', first, num, second)
  }
  return(change)
})
protein_changes <- unlist(protein_changes)

riaz_snv$HGVSp_Short <- protein_changes
write.table( riaz_snv , file=file.path(work_dir, "ICB_Riaz_SNV.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )

# ICB_Roh
roh_snv_original <- data.frame( fread( file.path(work_dir, 'ICB_Roh_SNV.txt.gz') , sep="\t" , stringsAsFactors=FALSE  ))
roh_snv <- read.csv(file.path(work_dir, 'ICB_Roh_SNV.csv'), sep=';')
roh_snv$Chr <- str_replace(roh_snv$Chr, 'chr', '')
roh_snv_original$Sample <- paste0('P.', str_replace(roh_snv_original$Sample, 'A', ''))
roh_snv_original <- roh_snv_original[roh_snv_original$Sample %in% unique(roh_snv$Sample), ]
roh_snv_original$Hugo_Symbol <- unlist(lapply(roh_snv_original$Hugo_Symbol, function(gene){
  gene <- unlist(str_split(gene, ';'))
  return(gene[1])
}))

protein_changes <- c()
for(i in 1:nrow(roh_snv)){
  found <- roh_snv_original[
    roh_snv_original$Sample %in% roh_snv[i, 1] & roh_snv_original$Hugo_Symbol %in% roh_snv[i, 2] & roh_snv_original$start %in% roh_snv[i, 4], 
    c('aaannotation')]
  if(length(found) > 1){
    found <- found[1]
  }
  protein_changes <- append(protein_changes, found)
}

roh_snv$HGVSp_Short <- protein_changes
write.table( roh_snv , file=file.path(work_dir, "ICB_Roh_SNV.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )

# ICB_Van_Allen
van_allen_snv_original <- data.frame( fread( file.path(work_dir, 'ICB_Van_Allen_SNV.txt.gz') , sep="\t" , stringsAsFactors=FALSE  ))
van_allen_snv <- read.csv(file.path(work_dir, 'ICB_Van_Allen_SNV.csv'), sep=';')
van_allen_snv$Chr <- str_replace(van_allen_snv$Chr, 'chr', '')

protein_changes <- c()
for(i in 1:nrow(van_allen_snv)){
  found <- van_allen_snv_original[
    van_allen_snv_original$patient %in% van_allen_snv[i, 1] & van_allen_snv_original$Hugo_Symbol %in% van_allen_snv[i, 2] & van_allen_snv_original$Start_position %in% van_allen_snv[i, 4], 
    c('Protein_Change')]
  if(length(found) > 1){
    found <- found[1]
  }
  protein_changes <- append(protein_changes, found)
}

van_allen_snv$HGVSp_Short <- protein_changes
write.table( van_allen_snv , file=file.path(work_dir, "ICB_Van_Allen_SNV.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )
