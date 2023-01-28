library(stringr)
library(data.table)

clin_colnames <- c("patientid","sex","age","cancer_type","histo","tissueid","treatmentid","stage","response.other.info","recist","response","treatment","dna","rna","survival_time_pfs","event_occurred_pfs","survival_time_os","event_occurred_os","survival_unit","survival_type","TMB_raw","nsTMB_raw","indel_TMB_raw","indel_nsTMB_raw","TMB_perMb","nsTMB_perMb","indel_TMB_perMb","indel_nsTMB_perMb","CIN","CNA_tot","AMP","DEL")

replacement_map <- list(
  patientid="PATIENT_ID",
  tissueid="TISSUE_ID",
  treatmentid="TREATMENT_ID",
  event_occurred_pfs="PFS_STATUS",
  survival_time_pfs="PFS_MONTHS",
  event_occurred_os="OS_STATUS",
  survival_time_os="OS_MONTHS"
)

# attribute lines
attr_display_names <- c("patient identifier","sex","age","cancer type","histo","tissue id","treatment id","stage","response other info","recist","response","treatment","dna","rna","survival time pfs","event occurred pfs","survival time os","event occurred os","survival unit","survival type","TMB_raw","nsTMB_raw","indel_TMB_raw","indel_nsTMB_raw","TMB_perMb","nsTMB_perMb","indel_TMB_perMb","indel_nsTMB_perMb","CIN","CNA_tot","AMP","DEL")

attr_desc <- c(
  "patient identifier",
  "sex","age",
  "Type of cancer tissue	source",
  "Histological info such as subtype	source",
  "Cancer type standardized based on the lab's nomenclature from Oncotree.",
  "Treatment regimen of each patient.",
  "Cancer stage	from source",
  "response other info",
  "Annotated using RECIST. The most commonly used responses are CR,PR,SD, PD.",
  "Response status of the patients to the given treatment - Responders (R) and Non-responders (NR)",
  "Drug target or drug name",
  "DNA sequencing type",
  "Type of rna processed data",
  "The time starting from taking the treatment to the occurrence of the event of interest",
  "Binary measurement showing whether the event of interest occurred (1) or not (0)",
  "The time starting from taking the treatment to the occurrence of the event of interest",
  "Binary measurement showing whether the event of interest occurred (1) or not (0)",
  "The unit in which the survival time is measured",
  "PFS or OS or both",
  "Tumor Mutation Burden raw values",
  "nsTMB_raw","indel_TMB_raw","indel_nsTMB_raw","TMB_perMb","nsTMB_perMb","indel_TMB_perMb","indel_nsTMB_perMb",
  "Calculated from CNA values",
  "Sum of total CNA/coverage, calculated from CNA values",
  "Sum of total AMP/coverage, calculated from CNA values",
  "Sum of total DEL/coverage, calculated from CNA values"
)

attr_datatype <- c("STRING","STRING","STRING","STRING","STRING","STRING","STRING","STRING","STRING","STRING","STRING","STRING","STRING","STRING","NUMBER","event_occurred_pfs","NUMBER","event_occurred_os","STRING","STRING","NUMBER","NUMBER","NUMBER","NUMBER","NUMBER","NUMBER","NUMBER","NUMBER","NUMBER","NUMBER","NUMBER","NUMBER")

format_clinical_patient_data <- function(work_dir, obj_name, clin){
  clin <- clin[, clin_colnames]
  new_colnames <- clin_colnames
  
  for(name in names(replacement_map)){
    new_colnames[new_colnames == name] <- replacement_map[[name]]
  }
  
  new_colnames <- str_to_upper(new_colnames)
  new_colnames <- str_replace_all(new_colnames, '\\.', '_')
  
  colnames(clin) <- new_colnames
  
  clin$PFS_STATUS[clin$PFS_STATUS == 1] <- '1:DECEASED'
  clin$PFS_STATUS[clin$PFS_STATUS == 0] <- '0:LIVING'
  clin$OS_STATUS[clin$OS_STATUS == 1] <- '1:Recurred/Progressed'
  clin$OS_STATUS[clin$OS_STATUS == 0] <- '0:DiseaseFree'
  
  attr_priority <- rep(1, length(new_colnames))
  
  attr_display_names <- paste(attr_display_names, collapse='\t')
  attr_desc <- paste(attr_desc, collapse='\t')
  attr_datatype <- paste(attr_datatype, collapse='\t')
  attr_priority <- paste(attr_priority, collapse='\t')
  
  attr_display_names <- paste0('#', attr_display_names, '\t', '\n')
  attr_desc <- paste0('#', attr_desc, '\t', '\n')
  attr_datatype <- paste0('#', attr_datatype, '\t', '\n')
  attr_priority <- paste0('#', attr_priority, '\t', '\n')
  
  file <- file.path(work_dir, obj_name, 'data_clinical_patients.txt')
  file.create(file)
  cat(
    attr_display_names, 
    attr_desc, 
    attr_datatype,
    attr_priority, 
    sep = "",
    file=file
  )
  
  write.table(
    clin,
    file=file,
    append=TRUE,
    quote=FALSE,
    sep='\t',
    na=""
  )
}
