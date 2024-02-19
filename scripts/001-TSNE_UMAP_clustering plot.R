# package preparation -----------------------------------------------------

# if(!require('devtools')) {install.packages('devtools')}
# library('devtools')
# install_github("immunedynamics/spectre")
#
# library("Spectre")
# ## Check if all required packages have been installed
# Spectre::package.check(type = 'spatial')
#
# ## Load all required packages
# Spectre::package.load(type = 'spatial')

# data subset and exploration --------------------------------------------------------------

## loading packages
# package preparation -----------------------------------------------------

# if(!require('devtools')) {install.packages('devtools')}
# library('devtools')
# install_github("immunedynamics/spectre")
#
# library("Spectre")
# ## Check if all required packages have been installed
# Spectre::package.check(type = 'spatial')
#
# ## Load all required packages
# Spectre::package.load(type = 'spatial')


# 01-data subset and exploration --------------------------------------------------------------

## loading packages
library(SummarizedExperiment)
library(SingleCellExperiment)
library(CATALYST)
library(diffcyt)
library(data.table)
# BiocManager::install('diffcyt')

# set work directory
setwd("E:/BaiduSyncdisk/WYL/1-ForMfdata/")

# loading sce data
sce <- readRDS(file = "data/sce_iMM_patInf.rds")
dim(sce)
# transform to singlecellexperiment object
sce <- as(sce, "SingleCellExperiment")
n_cells(sce)

unique(data.table(colData(sce)$CoreID))
# 57
unique(data.table(colData(sce)$patient_id))
# 32
table(sce$im_G)
#   iMM
# 94345
table(sce$phenotype)
dim(sce)


# 02-run t-SNE/UMAP ----------------------------------------------------------
set.seed(1)
sce_redu <- runDR(sce, "TSNE", cells = 500, features = c(
  "CD68", "CD14", "CD45",
  "CD15", "S100", "CD57", "FOXP3", "CD8", "Tryptase",
  "CD3", "CD20"
))
set.seed(1)
sce_redu <- runDR(sce_redu, "UMAP", cells = 1e3, features = c(
  "CD68", "CD14", "CD45",
  "CD15", "S100", "CD57", "FOXP3", "CD8", "Tryptase",
  "CD3", "CD20"
))


colData(sce_redu)$phenotype <- factor(colData(sce_redu)$phenotype, levels = c(
  "mono/mf", "Neutrophil", "DC",
  "mast cell", "CD8T", "CD3T", "NK", "B cell",
  "Treg", "other"
))

sce_redu@metadata[["cluster_codes"]][, "merging1"] <-
  factor(sce_redu@metadata[["cluster_codes"]][, "merging1"], levels = c(
    "mono/mf", "Neutrophil", "DC",
    "mast cell", "CD8T", "CD3T", "NK", "B cell", "Treg", "other"
  ))


# add subpopulation (Myeloid cells and Lymphocytes)
library(dplyr)
pheno_temp <- as.data.frame(sce_redu@metadata[["cluster_codes"]]$merging1)
colnames(pheno_temp) <- "phenotype"

pheno_temp <- pheno_temp %>%
  mutate(
    Mye_lym = {
      x <- rep("other", nrow(pheno_temp))
      x[phenotype %in% c("mono/mf", "Neutrophil", "DC", "mast cell")] <- "Myeloid"
      x[phenotype %in% c("CD8T", "CD3T", "NK", "B cell", "Treg")] <- "Lymphocyte"
      x
    }
  )
table(pheno_temp$Mye_lym)
sce_redu@metadata[["cluster_codes"]]$merging2 <- factor(pheno_temp$Mye_lym, levels = c("Myeloid", "Lymphocyte", "other"))

rm(pheno_temp)

# add phenotype data to colData
colData(sce_redu)$phenotypeML <- cluster_ids(sce_redu, "merging2")
# View(colData(sce_redu))

# save data after clustering
saveRDS(sce_redu, "results/001-TSNE_UMAP_clustering plot.rds")
load(file, envir = parent.frame(), verbose = FALSE)
# rm(list = ls())

# 03-Plots -------------------------------------------------------------------

# import data
sce_redu <- readRDS("results/001-TSNE_UMAP_clustering plot.rds")


# plot clusters
library(ggplot2)
clusters_tsne <- plotDR(sce_redu, "TSNE", color_by = "phenotype")

# color panel
clusters_tsne <- clusters_tsne + scale_color_manual(values = c(
  "#f21b3f", "#ffca3a", "#ff5129", "#ea6aa3",
  "#33a1fd", "#3dc6c9", "#29bf12", "#b7e4c7", "#9655cc", "#ccb9de"
))
clusters_tsne

clusters_umap <- plotDR(sce_redu, "UMAP", color_by = "phenotype")

# color panel
clusters_umap <- clusters_umap + scale_color_manual(values = c(
  "#f21b3f", "#ffca3a", "#ff5129", "#ea6aa3",
  "#33a1fd", "#3dc6c9", "#29bf12", "#b7e4c7", "#9655cc", "#ccb9de"
))
clusters_umap

# !! export
ggsave(filename = "figures/TSNE plot.pdf", clusters_tsne, device = "pdf", width = 4, height = 4, dpi = 300)



# clusters heatmap
cluster_heat <- plotExprHeatmap(sce_redu,
  features = c(
    "CD68", "CD14", "CD45",
    "CD15", "S100", "CD57", "FOXP3", "CD8", "Tryptase",
    "CD3", "CD20"
  ),
  by = "cluster_id", k = "merging1", k_pal = c(
    "#f21b3f", "#ffca3a", "#ff5129", "#ea6aa3",
    "#33a1fd", "#3dc6c9", "#29bf12", "#b7e4c7", "#9655cc", "#ccb9de"
  ),
  scale = "last", distance = "euclidean", bars = T, perc = T
)

pdf("figures/clusters heatmap.pdf", width = 7, height = 4)
cluster_heat
dev.off()

