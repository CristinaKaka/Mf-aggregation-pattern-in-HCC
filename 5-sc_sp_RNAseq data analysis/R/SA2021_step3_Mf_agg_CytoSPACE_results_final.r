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
workers <- 16
# plan("multicore", workers = workers)
plan("sequential")

# ! Load the processed results from cytospaces
combined_tibble <- readRDS("02_SA2021_cytospace_combined_tibble.rds")
combined_tibble_Mf <- combined_tibble %>%
    dplyr::filter(CellType %in% c("Mono-like", "Mph"))
# View(combined_tibble_Mf)
table(combined_tibble_Mf$CellType)

combined_tibble_Mf$ST_sample <- factor(combined_tibble_Mf$ST_sample)
levels(combined_tibble_Mf$ST_sample)

# select pre-treatment samples
pre_tb <- combined_tibble_Mf[grepl("HCC-[0-9]T|HCC-5", combined_tibble_Mf$ST_sample),]
pre_tb$ST_sample <- factor(pre_tb$ST_sample)
levels(pre_tb$ST_sample)

### Section 2: Load HCC single-cell (SC) transcriptomic data --------------------------------------------------
# Load packages
library(Seurat)
# options(Seurat.object.assay.version = "v5")

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
VlnPlot(seu, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 2) +
    NoLegend()

# preprocess data after merging
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
all_genes <- rownames(seu)
seu <- ScaleData(seu, features = all_genes)

# Apply SCTransform to each Seurat object 
vars_to_regress <- c("percent_mito", "percent_ribo") 
seu <- SCTransform(seu,
  vst.flavor = "v2",
  verbose = T,
  return.only.var.genes = F,
  vars.to.regress = vars_to_regress
)

# Set the default assay
DefaultAssay(seu) <- "SCT"

# # ! save the combined results
# saveRDS(seu,
#     file = "SA2021_03_cytospace_result_Mf_seu.rds"
# )
# # seu <- readRDS("SA2021_03_cytospace_result_Mf_seu.rds")


### Section 3: spatial spot_jitter --------------------------------------------------
library(Seurat)
library(ggplot2)

# Define macrophage aggregation
# table((seu$Mph_total >= 4), (seu$adj_Mph_count >= 8))
seu$Mf_agg <- factor(ifelse((seu$Mph_total >= 4), "Agg", "Sca"))

tb <- as.data.frame(table(seu$ST_sample, seu$Mf_agg)) %>%
  tidyr::spread(key = Var2, value = Freq) %>%
  dplyr::mutate(Sum = Agg + Sca) %>%
  dplyr::mutate(Agg_perc = Agg / Sum)
tb

# Define a function to add jitter
rand_jitter <- function(data, scale_factor) {
    interval <- diff(range(data)) / length(unique(data)) * scale_factor
    data + runif(length(data), -interval/2, interval/2)
}

  
#   Function to generate and save plots for each ST_sample
generate_plot <- function(sample_data) {
  # Apply jitter
  sample_data$Y_jitter <- rand_jitter(1-sample_data$row, scale_factor = 1.5)*1.75
  sample_data$X_jitter <- rand_jitter(sample_data$col, scale_factor = 1.5)
  
  # Define color based on aggregation category
  sample_data$color <- ifelse(sample_data$Mf_agg == "Sca", "#71C9DD", "#ca2d15")
  
  # Plot using ggplot2
  p <- ggplot(sample_data, aes(x = X_jitter, y = Y_jitter, color = color)) + 
    geom_point(size = 0.1, alpha = 0.9) +
    scale_color_identity() +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    coord_fixed() +  # Keep aspect ratio fixed
    labs(color = "Mf_agg category")
  
  # Define output path and filenames
  output_prefix <- unique(sample_data$ST_sample)
  png_filename <- paste0("figures/", output_prefix, "_Mph_distrib_jitter.png")
  pdf_filename <- paste0("figures/", output_prefix, "_Mph_distrib_jitter.pdf")
  
  # Save as PNG image
  ggsave(filename = png_filename, plot = p, width = 300/72, height = 300/72, dpi = 300)
  
  # Save as PDF
  ggsave(filename = pdf_filename, plot = p, width = 300/72, height = 300/72)
}

# Ensure directory exists
if (!dir.exists("figures")) {
  dir.create("figures")
}

# Loop through each ST_sample and apply the function
samples <- unique(seu@meta.data$ST_sample)
for(sample in samples) {
  sample_data <- subset(seu@meta.data, ST_sample == sample)
  generate_plot(sample_data)
}


### Section 4: DEG analysis and gene set enrichment --------------------------------------------------
library(rJava)
library(xlsx)
library(stringr)
library(dplyr)
library(tidyr)
library(Seurat)
 plan("multicore", workers = workers)

# Function to perform garbage collection in Java and R
jgc <- function() {
  .jcall("java/lang/System", method = "gc")
  gc()
}

# stress_genes <- c("G0S2", "JUN", "JUNB", "JUND", "FOS", "DUSP1", "CDKN1A", "FOSB", "BTG2", "KLF6", "KLF4")
    # cc_genes <- unique((unlist(cc.genes.updated.2019)))
    stress_genes <- cc_genes <- c()
    hist_genes <- grep("HIST", rownames(seu@assays$RNA), v = T)
    hb_genes <- c(grep("^HB[^(P)]", rownames(seu@assays$RNA), v = T))
    bad_features <- unique(c(
        hist_genes, cc_genes, stress_genes, hb_genes,
        # "MALAT1", "NEAT1", "XIST", "ACTB",
        grep("^MT-|^MTRNR2L|MTRNR2L|RP[SL]|^RP[SL]|^HSP|^DNAJ|^HSP|^DNAJ|RIK|AL|-RS|-PS|MIR|ATP|GM|UQC",
            rownames(seu@assays$RNA),
            v = T
        )
    ))


## ! DEG
meta <- seu@meta.data %>%
  mutate(Mf_agg = case_when(
    Mph_total >= 4 ~ "Agg",
    Mph_total == 1 ~ "Sca",
    Mph_total > 1 & Mph_total < 4 ~ "Int",
    TRUE ~ NA_character_
  ))
table(meta$Mf_agg)
seu$Mf_agg <- factor(meta$Mf_agg)

Idents(seu) <- seu$Mf_agg
# Perform differential expression analysis
# seu <- PrepSCTFindMarkers(seu)
features <- setdiff(rownames(seu@assays$SCT), bad_features)

sc_results <- FindMarkers(
  object = seu,
  ident.1 = "Agg",
  ident.2 = "Sca",
  features = features,
  test.use = "wilcox",
  assay = "SCT",
  only.pos = FALSE,
  min.pct = 0.1,
  logfc.threshold = 0,
  verbose = TRUE,
  densify = TRUE
)

write.xlsx(sc_results %>% as.data.frame(),
  file = "results/SA2021_AggSpa_DEGs_scRNAseq.xlsx",
  row.names = T,
  sheetName = "DEGs",
  append = FALSE
)
jgc()


# ! enriched GO terms
library(clusterProfiler)
library(org.Hs.eg.db) # for human
# library(org.Mm.eg.db) # for mouse
de_results <- sc_results
View(de_results)

universe_genes <- na.omit(unique(rownames(seu@assays$SCT)))
markersUP <- de_results %>%
    dplyr::filter(p_val_adj <= 0.05) %>%
    dplyr::filter(avg_log2FC >= 0.1) %>%
    arrange(dplyr::desc(avg_log2FC)) %>%
    rownames() %>%
    na.omit() %>%
    unique()
goUP <- enrichGO(
    gene = markersUP,
    universe = universe_genes,
    OrgDb = org.Hs.eg.db, # for human
    # OrgDb = org.Mm.eg.db, # for mouse
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    minGSSize = 50,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    readable = TRUE
)
write.xlsx(goUP,
        file = "results/SA2021_AggSpa_DEGs_scRNAseq.xlsx",
        row.names = T,
        sheetName = "goUP",
        append = TRUE
)

markersDN <- de_results %>%
    dplyr::filter(p_val_adj <= 0.05) %>%
    dplyr::filter(avg_log2FC <= -0.1) %>%
    arrange(dplyr::desc(avg_log2FC)) %>%
    rownames() %>%
    na.omit() %>%
    unique()
goDN <- enrichGO(
    gene = markersDN,
    universe = universe_genes,
    OrgDb = org.Hs.eg.db, # for human
    # OrgDb = org.Mm.eg.db, # for mouse
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    minGSSize = 50,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    readable = TRUE
)
write.xlsx(goDN,
        file = "results/SA2021_AggSpa_DEGs_scRNAseq.xlsx",
        row.names = T,
        sheetName = "goDN",
        append = TRUE
)

# ! GSEA
# library loading
library(fgsea)
library(msigdbr)
set.seed(123)
# check species for selection
msigdbr_species()

# check data in The Molecular Signatures Database (MSigDB) 
m_df = msigdbr(species = "Homo sapiens")
print(m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat), n=23)


# Prepare gene list 
de_results$Symbol <- rownames(de_results)
colnames(de_results)


# all genes
logFC_M <- de_results %>%
  dplyr::filter(p_val_adj <= 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(Symbol, avg_log2FC)

logFC <- logFC_M$avg_log2FC; names(logFC) = logFC_M$Symbol
gsea_list <- sort(logFC, decreasing = T)
gsea_list <- gsea_list[!is.na(names(gsea_list))]
length(gsea_list)
# [1] 2103

# significant genes
# select genes for go term analysis
# log2FC_cutoff = log2(2)
padj_cutoff = 0.05

print(m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat), n=23)
# # A tibble: 23 × 2
# gs_cat gs_subcat        
# <chr>  <chr>            
#   1 C1     ""               
# 2 C2     "CGP"            
# 3 C2     "CP"             
# 4 C2     "CP:BIOCARTA"    
# 5 C2     "CP:KEGG"        
# 6 C2     "CP:PID"         
# 7 C2     "CP:REACTOME"    
# 8 C2     "CP:WIKIPATHWAYS"
# 9 C3     "MIR:MIRDB"      
# 10 C3     "MIR:MIR_Legacy" 
# 11 C3     "TFT:GTRD"       
# 12 C3     "TFT:TFT_Legacy" 
# 13 C4     "CGN"            
# 14 C4     "CM"             
# 15 C5     "GO:BP"          
# 16 C5     "GO:CC"          
# 17 C5     "GO:MF"          
# 18 C5     "HPO"            
# 19 C6     ""               
# 20 C7     "IMMUNESIGDB"    
# 21 C7     "VAX"            
# 22 C8     ""               
# 23 H      "" 

# C2 CP:BIOCARTA
pathwaysDF <- msigdbr("human", category="C2", subcategory = "CP:BIOCARTA")
pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

fgseaRes <- fgsea(pathways, gsea_list, minSize = 10, maxSize = 500, nproc = 1)
fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ", "))
CPres <- fgseaRes
head(CPres)

# export xlsx
library(xlsx)
write.xlsx(CPres,
           file = "results/SA2021_AggSpa_DEGs_scRNAseq.xlsx",
           row.names = T,
           sheetName = "GSEA_CPBIOCARTA",
           append = TRUE
)
gc()

# C2 KEGG CP
pathwaysDF <- msigdbr("human", category="C2", subcategory = "CP:KEGG")
pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

fgseaRes <- fgsea(pathways, gsea_list, minSize = 10, maxSize = 500, nproc = 1)
fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ", "))
KEGGres <- fgseaRes
head(KEGGres)

library(xlsx)  
write.xlsx(KEGGres,
           file = "results/SA2021_AggSpa_DEGs_scRNAseq.xlsx",
           row.names = T,
           sheetName = "GSEA_KEGGCP",
           append = TRUE
)
jgc()

# CP:REACTOME
pathwaysDF <- msigdbr("human", category="C2", subcategory = "CP:REACTOME")
pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

fgseaRes <- fgsea(pathways, gsea_list, minSize = 10, maxSize = 500, nproc = 1)
fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ", "))
Reactres <- fgseaRes
head(Reactres)

library(xlsx)
write.xlsx(Reactres,
           file = "results/SA2021_AggSpa_DEGs_scRNAseq.xlsx",
           row.names = T,
           sheetName = "GSEA_CPReactome",
           append = TRUE
)
jgc()



# GO:BP
pathwaysDF <- msigdbr("human", category="C5", subcategory = "GO:BP")
pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

fgseaRes <- fgsea(pathways, gsea_list, minSize = 10, maxSize = 500, nproc = 1)
fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ", "))
GOBPres <- fgseaRes
head(GOBPres)

write.xlsx(GOBPres,
           file = "results/SA2021_AggSpa_DEGs_scRNAseq.xlsx",
           row.names = T,
           sheetName = "GSEA_GO_BP",
           append = TRUE
)
jgc()

# GO:MF
pathwaysDF <- msigdbr("human", category="C5", subcategory = "GO:MF")
pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

fgseaRes <- fgsea(pathways, gsea_list, minSize = 10, maxSize = 500, nproc = 1)
fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ", "))
GOMFres <- fgseaRes
head(GOMFres)

write.xlsx(GOMFres,
           file = "results/SA2021_AggSpa_DEGs_scRNAseq.xlsx",
           row.names = T,
           sheetName = "GSEA_GO_MF",
           append = TRUE
)
jgc()

# IMMUNESIGDB
pathwaysDF <- msigdbr("human", category="C7", subcategory = "IMMUNESIGDB")
pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

fgseaRes <- fgsea(pathways, gsea_list, minSize = 10, maxSize = 500, nproc = 1)
fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ", "))
IMMUNESIGDBres <- fgseaRes
head(IMMUNESIGDBres)

write.xlsx(IMMUNESIGDBres,
           file = "results/SA2021_AggSpa_DEGs_scRNAseq.xlsx",
           row.names = T,
           sheetName = "IMMUNESIGDBres",
           append = TRUE
)
jgc()

# HALLMARK 
pathwaysDF <- msigdbr("human", category="H", subcategory = "")
pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

fgseaRes <- fgsea(pathways, gsea_list, minSize = 10, maxSize = 500, nproc = 1)
fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ", "))
HALLres <- fgseaRes
head(HALLres)

write.xlsx(HALLres,
           file = "results/SA2021_AggSpa_DEGs_scRNAseq.xlsx",
           row.names = T,
           sheetName = "HALLres",
           append = TRUE
)
jgc()


### Section 5: Prepare for SCENIC --------------------------------------------------
# load libraries
library(Seurat)
library(tidyverse)
library(magrittr)
library(SCENIC)
library(SingleCellExperiment)
library(SCopeLoomR)
set.seed(123)


# Filtering genes for calculating time
exprMat <- seu@assays$RNA@data
cellInfo <- data.frame(seu@meta.data)
dim(cellInfo)
colnames(cellInfo)

loci1 <- which(rowSums(exprMat) > 1 * .01 * ncol(exprMat))
table(rowSums(exprMat) > 1 * .01 * ncol(exprMat))
exprMat_filter <- exprMat[loci1, ]
dim(exprMat_filter)

# Tranfer into loom files
loom <- build_loom("R/SA2021_step3_Mf_seu.loom", dgem = exprMat_filter)
# loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

### Section 6: Load SCENIC results --------------------------------------------------
# Required packages:
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)

# !reset your work diretory
SCENICdir <- "/path_to_data/HCC-sp-RNAseq/SCENIC/results"
scenicLoomPath <- file.path(SCENICdir, 'SA2021_step3_Mf_seu.scenic.loom')
motifEnrichmentFile <- file.path(SCENICdir, 'SA2021_step3_Mf_seu.motifs.csv')
regulonAucFile <- file.path(SCENICdir, 'SA2021_step3_Mf_seu.auc.csv')
BinarymatFile <- file.path(SCENICdir, 'SA2021_step3_Mf_seu.bin.csv')
regulonAucThresholdsFile <- file.path(SCENICdir, 'SA2021_step3_Mf_seu.thresholds.csv')

all(file.exists(scenicLoomPath), file.exists(motifEnrichmentFile), file.exists(regulonAucFile),file.exists(BinarymatFile), file.exists(regulonAucThresholdsFile))

# Retrieve AUC scores per cell
loom <- open_loom(scenicLoomPath)
  # Read information from loom file:
  regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
  regulons <- regulonsToGeneLists(regulons_incidMat)
  regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
close_loom(loom)

motifEnrichment <- data.table::fread(motifEnrichmentFile, header=T, skip=1)[-1,]
colnames(motifEnrichment)[1:2] <- c("TF", "MotifID")
motifEnrichment[1:5,1:5]
length(unique(motifEnrichment$TF))
# View(motifEnrichment)

regulonAUC2  <- data.table::fread(regulonAucFile, header=T, skip=0) %>%
column_to_rownames("Cell") %>%
t()
rownames(regulonAUC2) <- gsub("[(+)]", "", rownames(regulonAUC2))
regulonAUC2[1:5,1:5]
ncol(regulonAUC2)

regulonBin  <- data.table::fread(BinarymatFile, header=T, skip=0) %>%
column_to_rownames("Cell") %>%
t()
rownames(regulonBin) <- gsub("[(+)]", "", rownames(regulonBin))
regulonBin[1:5,1:5]

# get AUC matrix
AUCmat <- AUCell::getAUC(regulonAUC)
rownames(AUCmat) <- gsub("[(+)]", "", rownames(AUCmat))
AUCmat[1:5,1:5]

save(regulons_incidMat, regulons, regulonAUC, motifEnrichment, regulonAUC2, regulonBin,AUCmat,
  file = "SA2021_03_Mf_pySCENIC_results.RData"
)
# load(file = "SA2021_03_Mf_pySCENIC_results.RData")
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

# Incorporate data into seurat object 
# seu <- readRDS("SA2021_03_cytospace_result_Mf_seu.rds")

# merge AUC matrix into seurat object
seu[['AUC']] <- CreateAssayObject(data = regulonAUC2)

seu <- ScaleData(seu, assay = 'AUC', features = rownames(regulonAUC2))

# merge AUC-bin matrix into seurat object
seu[['Bin']] <- CreateAssayObject(data = regulonBin)
seu <- ScaleData(seu, assay = 'Bin', features = rownames(regulonBin))

DefaultAssay(seu) <- 'SCT'
saveRDS(seu, file = "SA2021_03_cytospace_result_Mf_seu.rds")

### Section 7: Regulon DEG analysis --------------------------------------------------
library(rJava)
library(xlsx)
library(stringr)
library(dplyr)
library(tidyr)
library(Seurat)

Idents(seu) <- seu$Mf_agg
# Function to perform garbage collection in Java and R
jgc <- function() {
  .jcall("java/lang/System", method = "gc")
  gc()
}

# regulon
seu_regulon <- FindMarkers(
  object = seu,
  ident.1 = "Agg",
  ident.2 = "Sca",
  features = NULL,
  test.use = "wilcox",
  assay = "AUC",
  logfc.threshold = 0.005,
  densify = TRUE
)
write.xlsx(seu_regulon %>% as.data.frame(),
  file = "results/SA2021_AggSpa_DEGs_scRNAseq.xlsx",
  row.names = T,
  sheetName = "regulon",
  append = TRUE
)
jgc()

# regulon_Bin
seu_regulon_bin <- FindMarkers(
  object = seu,
  ident.1 = "Agg",
  ident.2 = "Sca",
  features = NULL,
  test.use = "wilcox",
  assay = "Bin",
  logfc.threshold = 0.005,
  densify = TRUE
)
write.xlsx(seu_regulon_bin %>% as.data.frame(),
  file = "results/SA2021_AggSpa_DEGs_scRNAseq.xlsx",
  row.names = T,
  sheetName = "regulon_bin",
  append = TRUE
)
jgc()


# ! save the combined results
saveRDS(seu,
    file = "SA2021_03_cytospace_result_Mf_seu.rds"
)
# seu <- readRDS("SA2021_03_cytospace_result_Mf_seu.rds")