# data preparation --------------------------------------------------------

# packages loading
library(SummarizedExperiment)
library(SingleCellExperiment)
library(CATALYST)
library(dplyr)


# 01-data preparation --------------------------------------------------------

# read cell phenotype data
sce <- readRDS("results/afterTSNE_UMAP.rds")

# celldata including sample information
celldata <- as.data.frame(colData(sce))

# coreid and phenotype
library(tidyverse)

pheno = data.frame(table(celldata[,c("sample_id", "phenotype")]))

# transform to wider format
pheno <- pheno %>% 
  spread(key = phenotype, value = Freq)

# core information
col_names <- c(
  "patient_id", "sample_id", "TMA_ID", "CoreID", 
  "NumDet", "Area", "CD8G3_byDFS",
  "No.Diameter", "No.AFP", "No.ALT", "No.Path.stage", "No.Child", "No.Tumors", "No.Tumor", "No.Vasc.invad",
  "No.BCLC", "No.Cirhosis", "No.Capsule", "No.ALBI", "No.HBsAgn", "HBsAb", "HBeAg", "HBeAb", "HbcAb",
  "No.HCVAb", "No.Sex", "No.Age", 
   "No.T.1", "No.N", "No.M", "No.TNM",
  "DFS.time.5y", "DFS.outcome.5y"
)

colnames(celldata)
coreInf <- unique(celldata[,col_names])

# merge data
pheno <- merge(coreInf,pheno, by = "sample_id")

# total immune cell counts
ct <- c("mono/mf", "Neutrophil", "DC", "mast cell", "CD8T", "CD3T", "NK", "B cell", "Treg", "other")

pheno$im <- rowSums(pheno[,ct])
pheno$im_perA <- (pheno$im)/(pheno$Area)
# proportion
pheno <- pheno %>%
  mutate(100* across(.cols=ct,.names = '{.col}_inIC')/rowSums(across(.cols=ct)))

# cell number per mm2
pheno <- pheno %>%
  mutate(across(.cols=ct,.names = '{.col}_perA')/(Area))

# !! subset myeloid cells
mct <- c("mono/mf", "Neutrophil", "DC", "mast cell")

pheno$My <- rowSums(pheno[,mct])

# myeloid cell per mm2
pheno$My_perA <- (pheno$My)/(pheno$Area)

# total myeloid cell proportion
pheno$My_prop <- pheno$My/pheno$im

# myeloid subsets proportion
pheno <- pheno %>%
  mutate(100* across(.cols=mct,.names = '{.col}_inMy')/rowSums(across(.cols=mct)))

# !! subset lymphocytes
lct <- c("CD8T", "CD3T", "B cell", "NK", "Treg")

pheno$Ly <- rowSums(pheno[,lct])

# lymphocytes per mm2
pheno$Ly_perA <- (pheno$Ly)/(pheno$Area)

# total lymphocytes proportion
pheno$Ly_prop <- pheno$Ly/pheno$im

# myeloid subsets proportion
pheno <- pheno %>%
  mutate(100* across(.cols=lct,.names = '{.col}_inLy')/rowSums(across(.cols=lct)))


# !! export
write.csv(pheno,file = 'results/coredata.csv')
saveRDS(pheno,file = 'results/coredata.rds')

# patientdata
colnames(pheno)
col_names <- c("patient_id",
  "mono/mf_inIC", "Neutrophil_inIC", "DC_inIC", "mast cell_inIC", "CD8T_inIC", "CD3T_inIC",
  "NK_inIC", "B cell_inIC", "Treg_inIC", "other_inIC", "mono/mf_perA", "Neutrophil_perA",
  "DC_perA", "mast cell_perA", "CD8T_perA", "CD3T_perA", "NK_perA", "B cell_perA",
  "Treg_perA", "other_perA", "im_perA", "My", "My_perA", "My_prop", "mono/mf_inMy",
  "Neutrophil_inMy", "DC_inMy", "mast cell_inMy", "Ly", "Ly_perA", "Ly_prop", "CD8T_inLy",
  "CD3T_inLy", "B cell_inLy", "NK_inLy", "Treg_inLy"
)
pheno_numeric <- pheno[,col_names]


library(data.table)
setDT(pheno_numeric)
patdata <- pheno_numeric[,lapply(.SD, mean), by= patient_id]

col_names <- c(
   "patient_id", "TMA_ID",
  "CD8G3_byDFS",  "No.Diameter", "No.AFP", "No.ALT", "No.Path.stage",
  "No.Child", "No.Tumors", "No.Tumor", "No.Vasc.invad", "No.BCLC", "No.Cirhosis", "No.Capsule",
  "No.ALBI", "No.HBsAgn", "HBsAb", "HBeAg", "HBeAb", "HbcAb", "No.HCVAb", "No.Sex", "No.Age",
    "No.T.1", "No.N", "No.M", "No.TNM", "DFS.time.5y", "DFS.outcome.5y"
)

patInf <- unique(coreInf[,col_names])

# merge data
patdata <- merge(patInf, patdata, by = "patient_id")

# !! export
write.csv(patdata,file = 'results/20240121-patientdata.csv')
saveRDS(patdata,file = 'results/202401219-patientdata.rds')

# pause save dat
rm(list = ls())









