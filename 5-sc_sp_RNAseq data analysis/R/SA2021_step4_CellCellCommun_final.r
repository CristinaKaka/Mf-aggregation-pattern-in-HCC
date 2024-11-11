### Section 0: Set up environment --------------------------------------------------


# Clear the environment
rm(list = ls())
graphics.off()
gc()

# Source custom R profile for specific configurations/settings
source("~/.radian_profile")

# Set working directory
workdir <- "/path_to_data/HCC-sp-RNAseq"
setwd(workdir)

# Create directories for figures and results if they don't exist
# if (!file.exists(file.path(workdir, "figures"))) {
#     dir.create(file.path(workdir, "figures"))
# }
# if (!file.exists(file.path(workdir, "results"))) {
#     dir.create(file.path(workdir, "results"))
# }


### Section 1: Load scRNA-seq data (cytospace results) --------------------------------------------------
# Load packages
library(Seurat)
library(SeuratObject)
library(Matrix)
# options(Seurat.object.assay.version = "v5")
# library(future)
library(patchwork)
library(ggplot2)
library(dplyr)
library(readr)
set.seed(123)
# Parallel processing setup (uncomment if needed)
workers <- 16
# plan("multicore", workers = workers)
plan("sequential")

# ! Load the total scRNA-seq data
seu_all <- readRDS("SA2021_02_all_cytospace_result_seu")
seu_all
colnames(seu_all@meta.data)

### Section 2: Define cell types --------------------------------------------------
# combine mono-like macrphages and macrophages
seu_all$CellType_l1 <- seu_all$CellType
table(seu_all$CellType_l1)
seu_all$CellType_l1[seu_all$CellType %in% c("Mono-like", "Mph")]  <- "Mf"
table(seu_all$CellType_l1)
seu_all$Mf_total <- seu_all$"Mono-like" + seu_all$Mph
table(seu_all$Mf_total)

# rename cells
seu_all
# Replace "-" with "_" in cell names
seu_all <- RenameCells(seu_all, new.names = gsub("-", "_", Cells(seu_all)))
seu_all <- RenameCells(seu_all, new.names = gsub("-", "_", Cells(seu_all)))

# ### Section 3: Split samples and write gene expression in mtx format --------------------------------------------------
# Creat a unique ST_sample_spot_ID
seu_all$ST_spot_ID <- paste0(seu_all$ST_sample, "_", seu_all$SpotID)

# spots with aggregated macrophages
Agg_spot <- unique(seu_all$ST_spot_ID[seu_all$Mf_total >=4])
length(Agg_spot)
seu_agg <- seu_all[, seu_all$ST_spot_ID %in% Agg_spot]
seu_agg

# SCT
seu_agg <- PercentageFeatureSet(seu_agg, "^MT-", col.name = "percent_mito")
seu_agg <- PercentageFeatureSet(seu_agg, "^RP[SL]", col.name = "percent_ribo")
seu_agg <- NormalizeData(seu_agg)
seu_agg <- FindVariableFeatures(seu_agg, selection.method = "vst", nfeatures = 3000)
vars_to_regress <- c("percent_mito", "percent_ribo") 
seu_agg <- SCTransform(seu_agg,
  vst.flavor = "v2",
  verbose = T,
  return.only.var.genes = F,
  vars.to.regress = vars_to_regress
)
saveRDS(seu_agg,
    file = "SA2021_04_seu_agg.rds"
)


# spots with scattered macrophages
Sca_spot <- unique(seu_all$ST_spot_ID[seu_all$Mf_total == 1])
length(Sca_spot)
seu_sca <- seu_all[, seu_all$ST_spot_ID %in% Sca_spot]
seu_sca

# SCT
seu_sca <- PercentageFeatureSet(seu_sca, "^MT-", col.name = "percent_mito")
seu_sca <- PercentageFeatureSet(seu_sca, "^RP[SL]", col.name = "percent_ribo")
seu_sca <- NormalizeData(seu_sca)
seu_sca <- FindVariableFeatures(seu_sca, selection.method = "vst", nfeatures = 3000)
vars_to_regress <- c("percent_mito", "percent_ribo") 
seu_sca <- SCTransform(seu_sca,
  vst.flavor = "v2",
  verbose = T,
  return.only.var.genes = F,
  vars.to.regress = vars_to_regress
)
saveRDS(seu_sca,
    file = "SA2021_04_seu_sca.rds"
)


### Section 4: Cellchat for Agg_spots --------------------------------------------------
library(CellChat)
ptm = Sys.time()
cellchat_agg <- createCellChat(object = seu_agg,  group.by = "CellType_l1", assay = "SCT")
cellchat_agg <- setIdent(cellchat_agg, ident.use = "CellType_l1") 
levels(cellchat_agg@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_agg@idents)) # number of cells in each cell group

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat_agg@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat_agg <- subsetData(cellchat_agg) # This step is necessary even if using the whole database
future::plan("multisession", workers = 16) # do parallel
cellchat_agg <- identifyOverExpressedGenes(cellchat_agg)
cellchat_agg <- identifyOverExpressedInteractions(cellchat_agg)

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

# Compute the communication probability and infer cellular communication network
ptm <- Sys.time()
cellchat_agg <- computeCommunProb(cellchat_agg, type = "triMean")
cellchat_agg <- filterCommunication(cellchat_agg, min.cells = 10)

# Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat_agg)

# Infer the cell-cell communication at a signaling pathway level
cellchat_agg <- computeCommunProbPathway(cellchat_agg)

# Calculate the aggregated cell-cell communication network
cellchat_agg <- aggregateNet(cellchat_agg)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

# Compute the network centrality scores
cellchat_agg <- netAnalysis_computeCentrality(cellchat_agg, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways


# ! save the results
saveRDS(cellchat_agg, "SA2021_step4_cellchat_agg.rds")
# cellchat_agg <- readRDS("SA2021_step4_cellchat_agg.rds")


### Section 5: Cellchat for Sca_spots --------------------------------------------------
ptm = Sys.time()
cellchat_sca <- createCellChat(object = seu_sca,  group.by = "CellType_l1", assay = "SCT")
cellchat_sca <- setIdent(cellchat_sca, ident.use = "CellType_l1") 
levels(cellchat_sca@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_sca@idents)) # number of cells in each cell group

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat_sca@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat_sca <- subsetData(cellchat_sca) # This step is necessary even if using the whole database
future::plan("multisession", workers = 16) # do parallel
cellchat_sca <- identifyOverExpressedGenes(cellchat_sca)
cellchat_sca <- identifyOverExpressedInteractions(cellchat_sca)

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

# Compute the communication probability and infer cellular communication network
ptm <- Sys.time()
cellchat_sca <- computeCommunProb(cellchat_sca, type = "triMean")
cellchat_sca <- filterCommunication(cellchat_sca, min.cells = 10)

# Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat_sca)

# Infer the cell-cell communication at a signaling pathway level
cellchat_sca <- computeCommunProbPathway(cellchat_sca)

# Calculate the aggregated cell-cell communication network
cellchat_sca <- aggregateNet(cellchat_sca)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

# Compute the network centrality scores
cellchat_sca <- netAnalysis_computeCentrality(cellchat_sca, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

# ! save the results
saveRDS(cellchat_sca, "SA2021_step4_cellchat_sca.rds")
# cellchat_sca <- readRDS("SA2021_step4_cellchat_sca.rds")