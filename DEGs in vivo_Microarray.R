#############################
######### Load data #########
#############################

#### Install necessary R packages ###
## if you have not installed the necessary R packages, run the following codes ##
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# BiocManager::install(c("limma","biomaRt","ExpressionNormalizationWorkflow","sva","snm", "ggplot2","clusterProfiler","GSEABase","enrichplot","org.Hs.eg.db","GOSemSim","RColorBrewer","cowplot","magrittr"))
library(limma)

# set work path
setwd("...")

# array information
SDRF <- readRDS("dataforGit/Microarray data_sdrf.rds")

# loading patientdata
library(readxl)
Whole <- read_excel("dataforGit/pheno.xlsx")

## Read in raw data and other info ##
x <- read.maimages(SDRF[, "ArrayDataFile"],
    path = "data/Raw data",
    annotation = c("FeatureNum", "Row", "Col", "ProbeName"),
    source = "agilent", green.only = TRUE, other.columns = "gIsWellAboveBG"
)

## check if you properly have 62976 probes and 19 samples ##
dim(x)

#############################
###### Gene Anotation #######
#############################

## import the array annotation info from a file ##
annot.info <- read.csv("dataforGit/GPL20844-27159.csv",
    header = TRUE,
    comment.char = "#", na.strings = c("", "NA")
)

## check if you properly have 62976 probe annotation ##
length(unique(annot.info$ID))
# 62976

## merge expression data and annotation info ##
x$genes <- merge(x$genes, annot.info,
    by.x = "FeatureNum",
    by.y = "ID", sort = FALSE
)

#############################
### Background correction ###
#############################

y <- backgroundCorrect(x, method = "normexp")
y$E <- log2(y$E + 1)
boxplot(y$E, main = "expression value", las = 2)
# View(y$genes)

#############################
####### Gene filtering ######
#############################

## Gene filtering ##
isControl <- y$genes$CONTROL_TYPE != "FALSE"
yfilt <- y[!isControl, ]
dim(yfilt)

## clean up the SampleNames ##
colnames(yfilt$E) <- sub("_T@.*", "", colnames(yfilt$E))
colnames(yfilt$E) <- gsub("HCC_", "", colnames(yfilt$E))

## combine data measured by probes with same ProbeIDs ##
UniProbe <- aggregate(
    x = yfilt$E,
    by = list(yfilt$genes$ProbeName), FUN = mean
)
colnames(UniProbe)[1] <- "ProbeNames"
# View(UniProbe)
dim(UniProbe)
# [1] 58201    20
boxplot(UniProbe[, -1], las = 2, col = colnames(UniProbe[, -1]), main = "")


#############################
####### Normalization #######
#############################


### First, evaluate the batch effects ###
#########################################

library(ExpressionNormalizationWorkflow)
library(sva)
library(Biobase)

## Create an Expression Set object ##
# Re-arrang expression data
exprs <- UniProbe
row.names(exprs) <- exprs[, 1]
exprs[, 1] <- NULL
# View(exprs)

# import patient=(sample) information
covrts <-  Whole
rownames(covrts) <- covrts$SN

# re-order the expression matrix and the sample matrix
exprs <- exprs[,order(as.numeric(colnames(exprs)))]
covrts <- covrts[order(as.numeric(row.names(covrts))),]
rownames(covrts) <- covrts$SN
all(colnames(exprs)==rownames(covrts))

# construct ExpressionSet object
inpData <- expSetobj(exprs, covrts)
annotation(inpData) <- "org.Hs.eg.db"
dim(inpData)
colnames(covrts)

## Set the covariates whose effect size on the data needs to be calculated (optional)
cvrts_eff_var <- c(
    "Batch", "MRS_CD68_3", "MRS.group", "hs.CD68_img_X1", "RnaQC", "GoodPrognosis",
    "No.Sex", "PDL1.Mf", "PDL1.TC"
)
## Set a PVCA Threshold Value between 0 & 1
## PVCA Threshold Value is the percentile value of the minimum amount of the variabilities that the selected principal components need to explain, here requiring 75% of the expression variance to be captured by the PCs
pct_thrsh <- 0.75
## Perform the PVCA analysis
pvcAnaly(inpData, pct_thrsh, cvrts_eff_var)

### use "sva" method to removing batch   ###
### effects and other unwanted variation ###
############################################

### Setting up the data ###
pheno <- pData(inpData)
edata <- exprs(inpData)

## create the full model matrix - including both the adjustment variables and the variable of interest
mod <- model.matrix(~ as.factor(hs.CD68_img_X1) +
    as.factor(Batch) +
    as.factor(No.Sex) +
    as.factor(PDL1.Mf) +
    as.factor(PDL1.TC) +
    as.factor(MRS.group), data = pheno)
## the null model contains only the adjustment variables
mod0 <- model.matrix(~ as.factor(Batch) +
    as.factor(No.Sex) +
    as.factor(PDL1.Mf) +
    as.factor(PDL1.TC) +
    as.factor(MRS.group), data = pheno)

### Applying the sva function to  ########
### estimate batch and other artifacts ###
set.seed(1)
svobj <- sva(edata, mod, mod0)
colnames(svobj$sv) <- paste("sv", 1:svobj$n.sv, sep = "")
phenoSv <- cbind.data.frame(pheno, svobj$sv)
rownames(phenoSv) <- rownames(pheno)

### Computing the correlation between the  ###
### surrogate variables and the covariates ###
##############################################

glm.sv1 <- glm(phenoSv[, "sv1"] ~
    phenoSv[, "hs.CD68_img_X1"] + phenoSv[, "GoodPrognosis"] +
    phenoSv[, "Batch"] + phenoSv[, "RnaQC"] +
    phenoSv[, "No.Sex"] + phenoSv[, "PDL1.Mf"] +
    phenoSv[, "PDL1.TC"] + phenoSv[, "MRS.group"])
summary(glm.sv1)

glm.sv2 <- glm(phenoSv[, "sv2"] ~
    phenoSv[, "hs.CD68_img_X1"] + phenoSv[, "GoodPrognosis"] +
    phenoSv[, "Batch"] + phenoSv[, "RnaQC"] +
    phenoSv[, "No.Sex"] + phenoSv[, "PDL1.Mf"] +
    phenoSv[, "PDL1.TC"] + phenoSv[, "MRS.group"])
summary(glm.sv2)

glm.sv3 <- glm(phenoSv[, "sv3"] ~
    phenoSv[, "hs.CD68_img_X1"] + phenoSv[, "GoodPrognosis"] +
    phenoSv[, "Batch"] + phenoSv[, "RnaQC"] +
    phenoSv[, "No.Sex"] + phenoSv[, "PDL1.Mf"] +
    phenoSv[, "PDL1.TC"] + phenoSv[, "MRS.group"])
summary(glm.sv3)

glm.sv4 <- glm(phenoSv[, "sv4"] ~
    phenoSv[, "hs.CD68_img_X1"] + phenoSv[, "GoodPrognosis"] +
    phenoSv[, "Batch"] + phenoSv[, "RnaQC"] +
    phenoSv[, "No.Sex"] + phenoSv[, "PDL1.Mf"] +
    phenoSv[, "PDL1.TC"] + phenoSv[, "MRS.group"])
summary(glm.sv4)

### The results show that sv1 is correlated with Batch
### sv4  [, "PDL1.Mf"]  p.value = 0.0675

### Principal variance component analysis of the raw data to  ###
### estimate the contributions of the surrogate variables to  ###
### the overall expression variance                           ###
#################################################################

## First discretize the continuous surrogate variables
rownames(phenoSv)
colnames(edata)
inpData_sv <- expSetobj(edata, phenoSv)
var_names <- c("sv1", "sv2", "sv3", "sv4")
pData(inpData_sv) <- conTocat(pData(inpData_sv), var_names)

save.image(file = "dataforGit/before svn.RData")

## Include the surrogate variables as covariates in addition to BMI3, Rin3, CAD and Study (be sure to use categorical measures of BMI and Rin rather than the continuous measure)
cvrts_eff_var <- c(
    "Batch",
    "MRS.group", "hs.CD68_img_X1", "RnaQC", "GoodPrognosis",
    "No.Sex", "PDL1.Mf", "PDL1.TC",
    "sv1_cat", "sv2_cat", "sv3_cat", "sv4_cat"
)

## Again set the PVCA Threshold to explain 75% of the expression variation
pct_thrsh <- 0.75

## Perform PVCA
pvcAnaly(inpData_sv, pct_thrsh, cvrts_eff_var)

# rm(list = ls())

### Supervised normalization of Microarrays ###
###############################################

load("dataforGit/before svn.RData")

### Use independent SVs as adjustment variable
## Choose the biological variableof interest
bv <- c("hs.CD68_img_X1")
## Choose your adjustment variable of interest,
#  combine 'Batch' and independent SVs
av <- c("Batch", "sv2", "sv3", "sv4", "MRS.group")
## The intensity-dependent adjustment variables adjust for array effects
iv <- c("Array")
## Run SNM
sv_snmObj <- snmAnaly(edata, phenoSv, bv, av, iv)
## After comparison with other models, we chose this one for normalization
graphics.off()


## Create an expressionSet object of the normalized dataset(s)
sv_snmNorm_data <- sv_snmObj$norm.dat
# View(phenoSv)
dim(sv_snmNorm_data)
colnames(sv_snmNorm_data) <- rownames(phenoSv)
# View(sv_snmNorm_data)
sv_snm_data <- expSetobj(sv_snmNorm_data, phenoSv)
boxplot(sv_snmNorm_data, las = 2)

## Write this dateset to a table with rows as genes and columns as samples (with names the same as that from the input file)
## By doing this, you may pause the analaysis here
# !! export
write.csv(sv_snmNorm_data, file = "dataforGit/NormalizedData.csv")
save(sv_snmNorm_data, file = "dataforGit/NormalizedData.rdata")
write.csv(phenoSv, file = "dataforGit/Sample info for NormalizedData.csv")
save(phenoSv, file = "dataforGit/Sample info for NormalizedData.rdata")

### Pause point. You may run remove(list = ls()) to clean the runing environment. ###
# remove(list = ls())

##########################################
### Differential Expresssion Anlaysis ####
##########################################

## Load necessary R packages for this section
library(ExpressionNormalizationWorkflow)
library(Biobase)
library(limma)

## Set your working directory (or use Control + Shift + H to choose it)
setwd("...")

## Read phenotype data
load("dataforGit/Sample info for NormalizedData.rdata")
Pheno <- phenoSv

## Read normalized expression data
load("dataforGit/NormalizedData.rdata")

### Construct ExpresssionSet ###
sv_snm_data <- expSetobj(sv_snmNorm_data, Pheno)
annotation(sv_snm_data) <- "org.Hs.eg.db"

### Differential Expresssion Anlaysis ###
design1 <- model.matrix(~ 0 + as.factor(hs.CD68_img_X1), data = pData(sv_snm_data))
colnames(design1) <- c("sparse", "clustered")
design1
fit1 <- lmFit(exprs(sv_snm_data), design1)

cont.matrix1 <- makeContrasts("CLvsSP" = clustered - sparse, levels = design1)
cfit1 <- contrasts.fit(fit1, cont.matrix1)
efit1 <- eBayes(cfit1, trend = TRUE, robust = TRUE)
results1 <- decideTests(efit1)
summary(results1)
Output1 <- topTable(efit1, coef = 1, number = Inf)
dim(Output1)

#########################
### Annotate Results ####
#########################

Output1$NAME <- row.names(Output1)
Output1$Reg <- Output1$logFC > 0
Output1$Reg <- sub("TRUE", "Up", Output1$Reg)
Output1$Reg <- sub("FALSE", "Down", Output1$Reg)
Output1$FC <- 2^Output1$logFC
# View(Output1)

## Read probes annotation info
annot.info <- read.csv("data/GPL20844-27159.csv",
    header       = TRUE,
    comment.char = "#",  na.strings = c("", "NA")
)
colnames(annot.info)
UseAnnotInfo <- annot.info[, c(
    "NAME", "UNIGENE_ID", "GENE_SYMBOL",
    "GENE_NAME", "ENSEMBL_ID"
)]
dim(UseAnnotInfo)
UseAnnotInfo <- UseAnnotInfo[!duplicated(UseAnnotInfo), ]
dim(UseAnnotInfo)

## annotate differentially expressing probes
## (transfer probe names into unigene IDs)
PrAnnotResult1 <- merge(Output1, UseAnnotInfo,
    all.x = TRUE,
    sort = FALSE
)
dim(PrAnnotResult1)
colnames(PrAnnotResult1)[1] <- "ProbeID"

## exclude differentially expressing probes without Unigene IDs
## exclude differentially expressing probes without Unigene IDs
length(which(!is.na(PrAnnotResult1$ENSEMBL_ID)))
# 33390
length(which(!is.na(PrAnnotResult1$GENE_SYMBOL)))
# 48908
length(which(!is.na(PrAnnotResult1$GENE_NAME)))
# 47656
PrAnnotResult1 <- PrAnnotResult1[!is.na(PrAnnotResult1$GENE_SYMBOL), ]
length(PrAnnotResult1$GENE_SYMBOL)
length(unique(PrAnnotResult1$GENE_SYMBOL))

## annotate genes with EntrezGene IDs
library(AnnotationDbi)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
EnAnnotResult1 <- PrAnnotResult1

EnAnnotResult1$ENTREZID <- mapIds(org.Hs.eg.db,
    keys      = as.character(EnAnnotResult1$GENE_SYMBOL),
    column    = "ENTREZID",
    keytype   = "SYMBOL",
    multiVals = "first"
)
EnAnnotResult1$GeneName <- mapIds(org.Hs.eg.db,
    keys      = as.character(EnAnnotResult1$GENE_SYMBOL),
    column    = "GENENAME",
    keytype   = "SYMBOL",
    multiVals = "first"
)


## exclude differentially expressing probes without EntrezGene IDs
EnAnnotResult1 <- EnAnnotResult1[!is.na(EnAnnotResult1$ENTREZID), ]
# View(EnAnnotResult1)
dim(EnAnnotResult1)


## keep only one result data related to each EntrezID
EnAnnotResult1 <- EnAnnotResult1[order(EnAnnotResult1$adj.P.Val), ]
EnAnnotResult1 <- EnAnnotResult1[!duplicated(EnAnnotResult1$ENTREZID), ]
dim(EnAnnotResult1)

OutputGenes1 <- EnAnnotResult1

## re-order columns
OutputGenes1 <- cbind(
    EnAnnotResult1$ENTREZID, EnAnnotResult1$UNIGENE_ID,
    EnAnnotResult1$GENE_SYMBOL, EnAnnotResult1$GENE_NAME,
    EnAnnotResult1$ENSEMBL_ID, EnAnnotResult1[1:9]
)
colnames(OutputGenes1) <- gsub(pattern = "EnAnnotResult1\\$", replacement = "", colnames(OutputGenes1))
# View(OutputGenes1)

### Export the result of the differential expression analysis ###
write.csv(OutputGenes1,
    file = "dataforGit/Differential expression analysis result.csv",
    row.names = FALSE
)
save(OutputGenes1, file = "dataforGit/Differential expression analysis result.rdata")


### Pause point. You may run remove(list = ls()) to clean the runing environment. ###
# remove(list = ls())

##################
### load DEGs ####
##################

### load Expression Microarray result data
load("dataforGit/Differential expression analysis result.rdata")
Genes <- OutputGenes1
rm(OutputGenes1)
######################
### Volcanol plot ####
######################
library(ggplot2)
library(ggrepel)
Genes$Group <- NA
Genes$Group <- ifelse(Genes$adj.P.Val < 0.05 & abs(Genes$logFC) > 1, "Sig", "notSig")
Genes$Group[Genes$adj.P.Val < 0.01 & abs(Genes$logFC) > 1] <- ifelse(Genes$logFC[Genes$adj.P.Val < 0.01 & abs(Genes$logFC) > 1] > 0, "Upstr", "Dnstr")
Genes$Group[Genes$adj.P.Val < 0.05 & abs(Genes$logFC) > 1] <- ifelse(Genes$logFC[Genes$adj.P.Val < 0.05 & abs(Genes$logFC) > 1] > 0, "Up", "Dn")
Genes$Group <- factor(Genes$Group, levels = c("notSig", "Dn", "Up", "Dnstr", "Upstr"))
table(Genes$Group)
### Export the data of DEG volcano ###
write.csv(Genes,
    file = "dataforGit/DEG volcano.csv",
    row.names = FALSE
)
save(Genes, file = "dataforGit/DEG volcano.rdata")

## plot volcanol

# add label column
Genes$Label <- ""

select_G <- c(
    "ILDR2", "CDKN2C", "FOS", "HLA-DOB",
    "MAP3K15", "SPDYA", "ADRB1", "FLNC", "DCN",
    "WLS", "CD36", "ILDR2", "CCNE2", "KLRG1", "F3", "SYK",
    "CA12", "CXCL2", "MUC13", "HAS1", "HILPDA", "KLF4",
    "MMP2", "MRC2", "IL2RG", "CD40", "CD27", "KLRD1", "TNF", "NKG7", "SLAMF7"
)

# select_G <- c("CD36", "KLF4")

Genes$Label[match(select_G, Genes$GENE_SYMBOL)] <- select_G

ggplot(Genes, aes(x = logFC, y = -log10(adj.P.Val), color = Group)) +
    geom_point(alpha = 0.5, size = 1.5) +
    theme_bw(base_size = 12) +
    xlab("Log2(Fold change)") +
    ylab("-Log10(P.Val)") +
    theme(plot.title = element_text(size = 5, hjust = 0.5)) +
    scale_colour_manual(values = c("grey70", "#04a445", "#0090ef")) +
    geom_hline(yintercept = -log10(0.05), lty = 4) +
    geom_vline(xintercept = c(-1, 1), lty = 4) +
    # labs(title = this_tile)+
    scale_x_continuous(limits = c(-4.5, 5)) +
    scale_y_continuous(limits = c(0, 3.5)) +
    # geom_label_repel(data = Genes, aes(label = Label),
    #                   size          = 5,
    #                 box.padding   = unit(1.5, "lines"),
    #                   point.padding = unit(0.5, "lines"),
    #                 min.segment.length = 0.5,
    #                 # segment.color = "black",
    #                   # show.legend   = FALSE,
    #                   max.overlaps  = 900000
    #                   )
    geom_text_repel(
        data = Genes, aes(label = Label),
        box.padding = unit(1.5, "lines"),
        point.padding = unit(0.1, "lines"),
        max.overlaps = 900000
    )

ggsave("dataforGit/DEG volcano.pdf", width = 5.5, height = 5, units = "in")
dev.off()


# save all data -----------------------------------------------

save.image(file = "dataforGit/gene expression analysis.rdata")
