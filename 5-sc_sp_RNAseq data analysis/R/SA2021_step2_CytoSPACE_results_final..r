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


### Section 1: Load spatial transcriptomic data (cytospace results) --------------------------------------------------
# Load packages
library(Seurat)
# options(Seurat.object.assay.version = "v5")
# library(future)
library(patchwork)
library(ggplot2)
library(dplyr)
library(readr)
set.seed(123)
# Parallel processing setup (uncomment if needed)
# workers <- 4
# plan("multisession", workers = workers)
# plan("sequential")

# ! Load the results from cytospaces
# Specify directory path
dir_path <- "/path_to_data/HCC-sp-RNAseq/cytospace/output"

# List all "_assigned_locations.csv" files
location_files <- list.files(path = dir_path, pattern = "_assigned_locations.csv", full.names = TRUE)

# Initialize a list to store results
result_list <- list()
# Loop through each sample
for (file in location_files) {
    # Extract sample name
  sample_name <- gsub("_assigned_locations.csv", "", basename(file))
  # Read "_assigned_locations.csv" file
  df_locations <- read_csv(file, show_col_types = T)
  # Path for the corresponding "_cell_type_assignments_by_spot.csv" file
  cell_type_file <- file.path(dir_path, paste0(sample_name, "_cell_type_assignments_by_spot.csv"))
  # Read "_cell_type_assignments_by_spot.csv" file
  df_cell_type <- read_csv(cell_type_file, show_col_types = T)
  # Perform left_join based on "SpotID" column
  merged_df <- df_locations %>%
    left_join(df_cell_type, by = "SpotID")
  # Store the merged dataframe in the list
  result_list[[sample_name]] <- merged_df
}

# Add a column to indicate the sample name
for (sample_name in names(result_list)) {
  result_list[[sample_name]]$ST_sample <- sample_name
  result_list[[sample_name]]$SC_sample <- str_extract(result_list[[sample_name]]$OriginalCID, "^[^_]+_[^_]+")
  result_list[[sample_name]]$OldUniqueCID <- result_list[[sample_name]]$UniqueCID
  result_list[[sample_name]]$UniqueCID <- paste0(result_list[[sample_name]]$ST_sample,"_",result_list[[sample_name]]$OldUniqueCID)
}


# At this point, result_list contains the merged results for each sample, with one entry per sample.

# Combine all tibbles in the list into one tibble
combined_tibble <- bind_rows(result_list) 

# Print the number of duplicated OriginalCID values
total_ids <- length(combined_tibble$UniqueCID)
unique_ids <- length(unique(combined_tibble$UniqueCID))
duplicated_ids_count <- total_ids - unique_ids

cat("Total OriginalCID:", total_ids, "\n")
cat("Unique OriginalCID:", unique_ids, "\n")
cat("Duplicated OriginalCID:", duplicated_ids_count, "\n")

# length(unique(combined_tibble$OriginalCID))


# ! save the combined results
saveRDS(combined_tibble,
    file = "02_SA2021_cytospace_combined_tibble.rds"
)
# combined_tibble <- readRDS("02_SA2021_cytospace_combined_tibble.rds")
# View(combined_tibble)


### Section 2: Load HCC single-cell (SC) transcriptomic data --------------------------------------------------
# Load packages
library(Seurat)
# options(Seurat.object.assay.version = "v5")
# library(future)
library(patchwork)
library(ggplot2)
library(dplyr)
library(readr)
set.seed(123)
# Parallel processing setup (uncomment if needed)
workers <- 8
plan("multicore", workers = workers)
# plan("sequential")

# ! Load the processed results from cytospaces
combined_tibble <- readRDS("02_SA2021_cytospace_combined_tibble.rds")
# View(combined_tibble)

combined_tibble$ST_sample <- factor(combined_tibble$ST_sample)
levels(combined_tibble$ST_sample)

# select pre-treatment samples
pre_tb <- combined_tibble[grepl("HCC-[0-9]T|HCC-5", combined_tibble$ST_sample),]
pre_tb$ST_sample <- factor(pre_tb$ST_sample)
levels(pre_tb$ST_sample)
dim(pre_tb)

# Define the base directory and sample names
base_dir <- "/path_to_data/HCC-sp-RNAseq/cytospace/output"
sample_names <- levels(pre_tb$ST_sample)

# Loop through each sample and read the 10X data
seurat_list <- list()
for (sample in sample_names) {
  expression_dir <- file.path(base_dir, paste0(sample, "_assigned_expression"))
  seurat_obj <- CreateSeuratObject(counts = Read10X(expression_dir))
  seurat_obj <- RenameCells(object = seurat_obj, add.cell.id = sample)
  cells_to_keep <- intersect(Cells(seurat_obj), pre_tb$UniqueCID)
  
  # Print the number of duplicated OriginalCID value
  total_ids <- length(Cells(seurat_obj))
  intersect_ids <- length(cells_to_keep)
  cat(sample, "\n")
  cat("Total cells:", total_ids, "\n")
  cat("Remained cells:", intersect_ids, "\n")

  if (length(cells_to_keep) > 0) {
    seurat_obj <- seurat_obj[, Cells(seurat_obj) %in% cells_to_keep]
  } else {
    seurat_obj <- NULL
  }
  
  seurat_list[[sample]] <- seurat_obj
}

# Merge Seurat objects in the list
seu <- scCustomize::Merge_Seurat_List(
  seurat_list,
  add.cell.ids = NULL,
  merge.data = FALSE,
  project = "SeuratProject"
)

# Remove the first underscore from the cell names
seu <- RenameCells(seu, new.names = sub("^_", "", Cells(seu)))

# Check the first 5 cell names
print(Cells(seu)[1:5])

# Check cell names format and the unqiueCIDs in the pre_tb
head(seu@meta.data)
table(Cells(seu) %in% pre_tb$UniqueCID)

# Add pre_tb to the meta.data of the combined Seurat object
seu$UniqueCID <- Cells(seu)
meta_data <- dplyr::left_join(seu@meta.data, pre_tb, by = c("UniqueCID" = "UniqueCID"))
# View(meta_data)
# View(seu@meta.data)
seu@meta.data <- cbind(seu@meta.data, meta_data[ ,(ncol(seu@meta.data)+1):ncol(meta_data)])
head(seu@meta.data)
seu
# Now, combined_seurat is your desired Seurat object with pre_tb in its meta.data


# calculate QC
seu <- PercentageFeatureSet(seu, "^MT-", col.name = "percent_mito")
seu <- PercentageFeatureSet(seu, "^RP[SL]", col.name = "percent_ribo")

# plot QC
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo")
# VlnPlot(seu, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 2) +
#     NoLegend()

# preprocess data after merging
seu <- NormalizeData(seu)

# filter seurat object to keep unduplicated cells only
table(!duplicated(seu$OriginalCID))
seu_uni <- seu[, !duplicated(seu$OriginalCID)]
seu_uni <- NormalizeData(seu_uni)
seu_uni <- FindVariableFeatures(seu_uni, selection.method = "vst", nfeatures = 3000)
all_genes <- rownames(seu_uni)

# Apply SCTransform to each Seurat object 
vars_to_regress <- c("percent_mito", "percent_ribo") 
seu_uni <- SCTransform(seu_uni,
  vst.flavor = "v2",
  verbose = T,
  return.only.var.genes = F,
  vars.to.regress = vars_to_regress
)

# Set the default assay
DefaultAssay(seu_uni) <- "SCT"

# ! save the combined results
saveRDS(seu_uni,
    file = "SA2021_02_uni_cytospace_result_seu.rds"
)
saveRDS(seu,
    file = "SA2021_02_all_cytospace_result_seu.rds"
)
# seu <- readRDS("SA2021_02_all_cytospace_result_seu.rds")


### Section 3: Total cell dimension reduction and clustering --------------------------------------------------
## Run the standard workflow for visualization and clustering
# remain unique cells only

# PCA
# Define PCA features
# stress_genes <- c("G0S2", "JUN", "JUNB", "JUND", "FOS", "DUSP1", "CDKN1A", "FOSB", "BTG2", "KLF6", "KLF4")
# cc_genes <- unique((unlist(cc.genes.updated.2019)))
stress_genes <- cc_genes <- c()
hist_genes <- grep("HIST", rownames(seu_uni@assays$RNA), v = T)
hb_genes <- c(grep("^HB[^(P)]", rownames(seu_uni@assays$RNA), v = T))
bad_features <- unique(c(
  hist_genes, cc_genes, stress_genes, hb_genes,
  grep("^MT-|^MTRNR2L|MTRNR2L|RP[SL]|^RP[SL]|^HSP|^DNAJ|^HSP|^DNAJ|RIK|AL|-RS|-PS|MIR|ATP|GM|UQC",
  rownames(seu_uni@assays$RNA),
        v = T
      )
))

PCA_features <-  setdiff(seu_uni@assays$SCT@var.features, bad_features)

# run PCA
seu_uni <- RunPCA(seu_uni,
    assay = "SCT",
    npcs = 50, verbose = T,
    features = PCA_features
)

# data integration using Harmony
library(harmony)
colnames(seu_uni@meta.data)
table(seu_uni$SC_sample)

seu_uni <- RunHarmony(
  seu_uni,
  c("SC_sample")
)

## UMAP and TSNE Clustering
library(future)
nworkers <- 16
# plan("sequential")
plan("multicore", workers = nworkers)
set.seed(123)

# seu_int <- RunTSNE(seu_int, reduction = "harmony", dims = 1:30)
seu_uni <- RunUMAP(seu_uni, reduction = "harmony", dims = 1:30)
seu_uni <- FindNeighbors(seu_uni, reduction = "harmony", dims = 1:30)

## estimate clustering
# try different resolutions
# seu_uni@meta.data[,grep("SCT_snn_res.",colnames(seu_uni@meta.data),value = T)] <- NULL
seu_uni <- FindClusters(
    seu_uni,
    resolution = c(seq(0.01, 0.1, 0.01), seq(0.2, 1, 0.1))
)
# relevel cluster names (we do not what cluster "0")
for (i in 1:length(grep("SCT_snn_res.", colnames(seu_uni@meta.data)))) {
    j <- grep("SCT_snn_res.", colnames(seu_uni@meta.data))[i]
    k <- seu_uni@meta.data[, j]
    levels(k) <- as.character(1:length(levels(k)))
    seu_uni@meta.data[, j] <- k
}

# estimate clustering in different resolustion
library(clustree)
clustree(seu_uni@meta.data, prefix = "SCT_snn_res.", return = "plot")

# ! resolution choice
res <- 0.03

# Plot 1: see group difference in clusters-UMAP
DimPlot(seu_uni,
    reduction = "umap", 
    # split.by = "Mf_agg",
    group.by = paste0("SCT_snn_res.", res), 
    label = TRUE
) +
    coord_fixed(ratio = 1)


# ! set resolution
Idents(seu_uni) <- seu_uni$seurat_clusters <- factor(seu_uni@meta.data[[paste0("SCT_snn_res.", res)]])
# seu_uni@meta.data[,grep("SCT_snn_res.",colnames(seu_uni@meta.data),value = T)] <- NULL

VlnPlot(seu_uni, features = c("nCount_RNA","nFeature_RNA", "percent_mito", "percent_ribo"), pt.size = 0)

colnames(seu_uni@meta.data)

# View cluster annotation
cluster_annot <- tibble(
    seu_uni$seurat_clusters,
    seu_uni$CellType
) %>%
    set_names("Cluster", "CellType") %>%
    group_by(Cluster, CellType) %>%
    summarise(no.cell = n()) %>%
    group_by(Cluster) %>%
    mutate(
        total.no = sum(no.cell),
        perc = 100 * no.cell / total.no
    ) %>%
    arrange(Cluster, dplyr::desc(perc)) %>%
    top_n(n = 5, wt = perc)

View(cluster_annot)

# Plot 1: see group difference in clusters-UMAP
DimPlot(seu_uni,
    reduction = "umap", 
    group.by = "CellType", 
    # group.by = paste0("SCT_snn_res.", res), 
    label = TRUE
) +
    coord_fixed(ratio = 1)


# ! save the combined results
saveRDS(seu_uni,
    file = "SA2021_02_uni_cytospace_result_seu.rds"
)
# seu_uni <- readRDS("SA2021_02_uni_cytospace_result_seu.rds")