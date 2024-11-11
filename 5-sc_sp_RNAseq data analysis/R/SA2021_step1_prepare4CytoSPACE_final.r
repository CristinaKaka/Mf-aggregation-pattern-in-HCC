### Section 0: Set up environment --------------------------------------------------

# clear the environment
rm(list = ls())
graphics.off()
gc()
source("~/.radian_profile")

# !set your work diretory
workdir <- "/path_to_data/HCC-sp-RNAseq"
setwd(workdir)

# creat diretories for figures and results
if (file.exists(file.path(workdir, "figures"))) {
} else {
    dir.create(file.path(workdir, "figures"))
}
if (file.exists(file.path(workdir, "results"))) {
} else {
    dir.create(file.path(workdir, "results"))
}


### Section 1: Data preprocessing for spatial (SP) transcriptomic data --------------------------------------------------

## load packages
library(Seurat)
library(patchwork)
library(ggplot2)
set.seed(123)

## Load data
root_path <- "/path_to_data/HCC-sp-RNAseq/STDATA_SA"
sub_dirs <- list.dirs(path = root_path, full.names = TRUE, recursive = FALSE)
sub_dirs <- sub_dirs[1:20]

# a editted version of Load10X_Spatial
my_Load10X_Spatial <- function(
    data.dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "slice1",
    filter.matrix = TRUE,
    to.upper = FALSE,
    image = NULL,
    ...) {
  data <- Read10X_h5(filename = file.path(data.dir, filename))
  if (to.upper) {
    rownames(x = data) <- toupper(x = rownames(x = data))
  }
  object <- CreateSeuratObject(counts = data, assay = assay)
  image <- Read10X_Image(
    image.dir = file.path(data.dir, "spatial"),
    filter.matrix = filter.matrix
  )
  image <- image[Cells(x = object)]
  DefaultAssay(object = image) <- assay
  object[[slice]] <- image
  return(object)
}

# Read each file and save it to a named list
spt_list <- list()
for (dir in sub_dirs) {
  sample_name <- basename(dir)
  seurat_obj <- my_Load10X_Spatial (
    data.dir = dir
  )
  spt_list[[sample_name]] <- seurat_obj
  rm(seurat_obj)
}


### Section 2: Preparing input files from these ST data for CytoSPACE--------------------------------------------------
# Load necessary libraries
library(Matrix)

# Source the functions from the provided script
source('/path_to_cytospace/cytospace/Prepare_input_files/generate_cytospace_from_seurat_object.R')

# Set the output directory for the CytoSPACE input files
output_dir <- "/path_to_data/HCC-sp-RNAseq/cytospace/input/ST"

# Iterate over the spt_list and generate CytoSPACE input files for each ST Seurat object
for (sample_name in names(spt_list)) {
    cat("Processing sample:", sample_name, "\n")
    
    # Determine the appropriate slice name for the current Seurat object
    slice_name <- names(spt_list[[sample_name]]@images)[1]
    
    # Generate CytoSPACE input files for the current ST Seurat object
    generate_cytospace_from_ST_seurat_object(
        st_seurat = spt_list[[sample_name]], 
        dir_out = output_dir, 
        fout_prefix = paste0(sample_name, "_"), 
        write_sparse = FALSE, 
        slice = slice_name
    )
}


### Section 3: Preparing input files from sc-RNAseq data from ZHANG Ning et al ------------------------------------------------

# ! load data
# Load the Seurat object containing scRNA-seq data
seu <- readRDS("/path_to_SCdata/seu_A124.counts_anno.rds")

# Load the UMAP embeddings for the data
allUMAP <- readRDS("/path_to_SCdata/umap_embedding/seu_A124_all.umap.rds")

# Check the metadata columns in the Seurat object
colnames(seu@meta.data)

# Ensure that the cell order in the Seurat object matches the order in the UMAP embeddings
all(colnames(seu) == rownames(allUMAP))

# Add the UMAP embeddings to the Seurat object
seu[["umap"]] <- allUMAP  # Add the DimReduc object to the `reductions` slot of the Seurat object

# Adjust the cluster annotations
seu$clusters <- factor(seu$clusters)  # Convert 'clusters' to a factor
levels(seu$clusters)  # Display the levels of the 'clusters' factor

# Extract the primary cell type labels by removing any sub-type information
seu$CellType_l1 <- gsub("_.+$", "", seu$clusters)
table(seu$CellType_l1)  # Display the levels of the primary cell type labels
table(seu$Cancer_type, seu$CellType_l1)

# Recode the levels of 'CellType_l1' to merge specific cell types
# seu$CellType_l1[seu$CellType_l1 %in% c("MonoDC", "DC")] <- "DC"
# seu$CellType_l1[seu$CellType_l1 %in% c("Mono-like", "Mph")] <- "Mph"
# seu$CellType_l1[seu$CellType_l1 %in% c("Mu", "Fb")] <- "MC"
seu$CellType_l1[grepl("FOXP3", seu$clusters)] <- "Treg"


# Subset the seu object to retain only samples containing "HCC" in the Sample column
table(grepl("HCC", unique(seu$Sample)))
table(grepl("HCC", unique(seu$Sample)))

length(unique(seu$Sample)) # see how many samples are included in the original dataset
table(seu$Sample, seu$Cancer_type) # see Sample and Cancer_type catogories
length(unique(seu$Sample[seu$Cancer_type %in% c("AL", "HCC")])) # see how many samples are from the adjacent liver and tumor site of HCCs 
table(seu$Cancer_type %in% c("AL", "HCC"))
table(seu$Cancer_type %in% c("AL", "HCC") & seu$CellType_l1 %in% c("Mo", "MonoDC")) # see how many cells from HCCs are annotated as "Mo" or "MonoDC"
table(seu$Cancer_type %in% c("AL", "HCC") & seu$CellType_l1 %in% c("Mast")) # see how many cells from HCCs are annotated as "Mo" or "MonoDC"


seu_filtered <- seu[, seu$Cancer_type %in% c("AL", "HCC") & !(seu$CellType_l1 %in% c("Mo", "Mast", "MonoDC"))]

# Check the updated counts for each cancer type and cell type in the filtered Seurat object
table(seu_filtered$Cancer_type, seu_filtered$CellType_l1)

# Randomly sample 100,000 cells from the filtered Seurat object
# set.seed(123)  # Setting a seed for reproducibility
# cells_to_keep <- sample(colnames(seu_filtered), 100000)

# Subset the Seurat object to only keep the randomly sampled cells
# seu_filtered <- subset(seu_filtered, cells = cells_to_keep)

# Update the factor levels after recoding
seu_filtered$CellType_l1 <- factor(seu_filtered$CellType_l1)

# Set the identity class for the Seurat object to 'CellType_l1'
Idents(seu_filtered) <- "CellType_l1"

# !  save the filtered scRNA-seq seurat object
saveRDS(seu_filtered,
  file = "01_sc_filtered.rds"
)
# seu_filtered <- readRDS("01_sc_filtered.rds")

# Set the output directory for the CytoSPACE input files
output_dir_scRNA <- "/path_to_data/HCC-sp-RNAseq/cytospace/input/SC"

# Load necessary libraries
library(Matrix)

# Source the functions from the provided script
source('/path_to_cytospace/cytospace/Prepare_input_files/generate_cytospace_from_seurat_object.R')

# Generate CytoSPACE input files for the scRNA-seq Seurat object
generate_cytospace_from_scRNA_seurat_object(
    seu_filtered,
    dir_out = output_dir_scRNA,
    fout_prefix = "scRNA_",
    write_sparse = FALSE
)

# find cluster markers
Idents(seu_filtered) <- "CellType_l1"
levels(Idents(seu_filtered))
cl_markers <- FindAllMarkers(
    seu_filtered,
    test.use = "roc",
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.25,
    verbose = TRUE,
    densify = TRUE
)


# View the cluster markers
# View(cl_markers)

# !  save the markers
saveRDS(cl_markers,
  file = "01_SC_cluster_markers.rds"
)

# for exporting results
library(rJava)
library(xlsx)
library(stringr)
jgc <- function()
{
  .jcall("java/lang/System", method ="gc")
  gc()
  
}

# export DEG results
write.xlsx(cl_markers %>% as.data.frame(), 
           file = "results/01_SC_cluster_markers.xlsx", 
           row.names = T,
           sheetName = "roc", 
           append = F)
jgc()