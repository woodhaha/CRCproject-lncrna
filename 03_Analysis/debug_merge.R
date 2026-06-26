suppressMessages(library(dplyr))
suppressMessages(library(data.table))

# Load clinical
clinical_raw <- read.csv("../02_Data/raw/tcga_coad_clinical_followup_data.csv", header=TRUE, stringsAsFactors=FALSE)
clinical <- data.frame(
  id = clinical_raw$tcga_participant_barcode,
  OS = ifelse(!is.na(clinical_raw$days_to_death), clinical_raw$days_to_death, clinical_raw$days_to_last_followup),
  status = ifelse(clinical_raw$vital_status == "dead", 1, 0),
  stage = tolower(gsub("stage ", "", clinical_raw$pathologic_stage)),
  age = clinical_raw$age_at_initial_pathologic_diagnosis,
  gender = clinical_raw$gender,
  stringsAsFactors = FALSE
)

# Load expression
Exprs <- data.table::fread("../02_Data/raw/COAD_HTSeq_FPKM.txt", sep="\t", header=TRUE, data.table=FALSE, stringsAsFactors=FALSE)
colnames(Exprs) <- gsub("-", ".", colnames(Exprs))
rownames(Exprs) <- Exprs$ID
Exprs <- Exprs[, -1, drop=FALSE]
# Remove .11 non-tumor
idx_11 <- grep("\\.11", colnames(Exprs))
if(length(idx_11) > 0) Exprs <- Exprs[, -idx_11, drop=FALSE]
cat("Samples after removing .11:", ncol(Exprs), "\n")

colnames(Exprs) <- substr(colnames(Exprs), 1, 12)
expr_ids <- unique(colnames(Exprs))
cat("Unique patient IDs:", length(expr_ids), "\n")

# Check match
clinicdata <- clinical %>% dplyr::select(id, stage, OS, status)
clinicdata <- na.omit(clinicdata)
cat("Clinicdata rows:", nrow(clinicdata), "\n")
cat("Matching IDs:", sum(clinicdata$id %in% expr_ids), "/", nrow(clinicdata), "\n")
cat("Clinical first 5:", head(clinicdata$id, 5), "\n")
cat("Expression first 5:", head(expr_ids, 5), "\n")
