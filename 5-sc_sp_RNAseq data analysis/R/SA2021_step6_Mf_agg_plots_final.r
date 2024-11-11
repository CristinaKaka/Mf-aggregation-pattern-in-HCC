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


### Section 1: Load spatial transcriptomic data  --------------------------------------------------
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

# ! load the combined results
seu <- readRDS("SA2021_03_cytospace_result_Mf_seu.rds")

# Define macrophage aggregation
table(seu$Mph_total)

# table((seu$Mph_total >= 4), (seu$adj_Mph_count >= 8))
meta <- seu@meta.data %>%
  mutate(Mf_agg = case_when(
    Mph_total >= 4 ~ "Agg",
    Mph_total == 1 ~ "Sca",
    Mph_total > 1 & Mph_total < 4 ~ "Int",
    TRUE ~ NA_character_
  ))
table(meta$Mf_agg)

seu$Mf_agg <- factor(meta$Mf_agg)
rm(meta)

Idents(seu) <- seu$Mf_agg


### Section 2: spatial spot_jitter --------------------------------------------------
library(Seurat)
library(ggplot2)

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
  sample_data$color <- ifelse(sample_data$Mf_agg == "Sca", "#71C9DD",
      ifelse(sample_data$Mf_agg == "Agg", "#ca2d15", "#ac15ca")
  )
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
  png_filename <- paste0("figures/", output_prefix, "_Mph_distrib_jitter_bin3.png")
  pdf_filename <- paste0("figures/", output_prefix, "_Mph_distrib_jitter_bin3.pdf")
  
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


### Section 3: Feature expression dot plots --------------------------------------------------
Idents(seu) <- seu$Mf_agg
feature_genes <- c(
    # M2 like or “remodeler-like”
    "CD163", "MSR1", "MRC1",
    # lipid metabolism and LA-TAM (from Ma et al, "Macrophage diversity in cancer revisited in the era of single-cell omics", 2022)
    "APOC1","APOE", "FABP5","GPNMB","ACP5","PLA2G7","TREM2", 
    # complement components
    "C1QA", "C1QB", "C1QC", 
    # cathepsin proteases
    "CTSB", "CTSD", "CTSL",
    # # Mφ reprogramming and an immunosuppressive state
    # "SLC40A1","FOLR2","HES1",
    # Mf-c1-THBS1 from Landscape and Dynamics of Single Immune Cells in Hepatocellular Carcinoma, 2019
    "SERPINB2","VCAN", "FCN1", "S100A12", "S100A8", "S100A9","S100A6"
)
p1 <- DotPlot(seu,
    assay = "SCT",
    c("white", "#f02364"),
    features = feature_genes
) +
#   coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_blank(),
        legend.direction = "vertical",
        legend.position = "bottom")
p1 #+ scale_color_gradient2(low = "#5bacdb", mid = "white", high = "#ed7fa2", midpoint = 0)

png_filename <- paste0("figures/", "Mph_feature_dot_plot.png")
pdf_filename <- paste0("figures/", "Mph_feature_dot_plot.pdf")
ggsave(file = png_filename, plot = p1, width = 200, height = 100, units = "mm", dpi = 300, device = "png")
ggsave(file = pdf_filename, plot = p1, width = 200, height = 100, units = "mm", device = "pdf", bg = 'transparent')

# Trascription factors / regulons
tf_expression <- FetchData(seu, assay = "Bin", vars = c("MAF","MAFB", "TCF4"))
tf_expression$Mf_agg <- seu$Mf_agg
head(tf_expression)
# save as a csv file
write.csv(tf_expression, file = "results/TF_expression.csv", row.names = TRUE)


### Section 4: Ligand–receptor pair expression analysis --------------------------------------------------

# Load necessary libraries for data manipulation and analysis.
library(scater)
library(SingleCellExperiment)
library(rJava)
library(xlsx)
library(stringr)


# Extract narmalized counts and metadata to create SingleCellExperiment object
Idents(seu) <- seu$Mf_agg
counts <- LayerData(seu, assay = "SCT", layer = "counts")
metadata <- seu@meta.data

avg_expression <- as.data.frame(AverageExpression(seu, verbose = T, assays = "SCT", slot = "scale.data"))
head(avg_expression)


# Section 5: Comparison analysis of Agg_spot vs Sca_spt datasets using CellChat --------------------------------------------------
# Load the required CellChat objects and merge them together
library(CellChat)
library(rJava)
library(xlsx)
library(stringr)

# Function to perform garbage collection in Java and R
jgc <- function() {
  .jcall("java/lang/System", method = "gc")
  gc()
}

# read cell chat data
cellchat_sca <- readRDS("SA2021_step4_cellchat_sca_20231229.rds")
cellchat_agg <- readRDS("SA2021_step4_cellchat_agg_20231229.rds")
object.list <- list(Sca = cellchat_sca, Agg = cellchat_agg)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
output_dir <- "/ifs1/User/yangbinVM/01.Project/wuchong/data/WYL/HCC-sp-RNAseq/figures/cellchat"

## 1. Differential number of interactions or interaction strength among different cell populations 
# Visualize and save interaction strength using heatmap
par(mfrow = c(1,1), xpd=TRUE)
pdf(file.path(output_dir, "Fig.4a_Differential_Interaction_Heatmap.pdf"))
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

## 2. plot the signaling sent from or received by macrophages 
wt_mat_sca <- cellchat_sca@net$weight
wt_mat_agg <- cellchat_agg@net$weight
groupSize_sca <- as.numeric(table(cellchat_sca@idents))
groupSize_agg <- as.numeric(table(cellchat_agg@idents))

# extract macrophage-related communications
# scattered
mat_mf_sca <- matrix(0, nrow = nrow(wt_mat_sca), ncol = ncol(wt_mat_sca), dimnames = dimnames(wt_mat_sca))
mat_mf_sca["Mf",] <- wt_mat_sca["Mf",]
mat_mf_sca[,"Mf"] <- wt_mat_sca[,"Mf"]
# aggregated
mat_mf_agg <- matrix(0, nrow = nrow(wt_mat_agg), ncol = ncol(wt_mat_agg), dimnames = dimnames(wt_mat_agg))
mat_mf_agg["Mf",] <- wt_mat_agg["Mf",]
mat_mf_agg[,"Mf"] <- wt_mat_agg[,"Mf"]

# control the parameter edge.weight.max so that we can compare edge weights between groups
max_wt <- max(mat_mf_sca, mat_mf_agg)

# Compare the interaction strength
par(mfrow = c(1,2), xpd=TRUE)
pdf(file.path(output_dir, "Fig.4b_Mf_interaction_strength.pdf"))
netVisual_circle(mat_mf_sca,
  vertex.weight = groupSize_sca,
  weight.scale = T,
  edge.weight.max = max_wt,
  title.name = "Sca-Mf spots"
)
netVisual_circle(mat_mf_agg,
  vertex.weight = groupSize_agg,
  weight.scale = T,
  edge.weight.max = max_wt,
  title.name = "Agg-Mf spots"
)
dev.off()

## 3. export cell-cell communication data
net_sca <- subsetCommunication(cellchat_sca,
  slot.name = "net",
  # sources.use = source_cells, targets.use = target_cells,
  thresh = 0.05
)
write.xlsx(cellchat_sca@net$count %>% as.data.frame(),
  file = "results/SA2021_AggSpa_cellchat_net.xlsx",
  row.names = TRUE,
  sheetName = "net_sca_count",
  append = FALSE
)
write.xlsx(cellchat_sca@net$weight %>% as.data.frame(),
  file = "results/SA2021_AggSpa_cellchat_net.xlsx",
  row.names = TRUE,
  sheetName = "net_sca_weight",
  append = TRUE
)
write.xlsx(net_sca %>% as.data.frame(),
  file = "results/SA2021_AggSpa_cellchat_net.xlsx",
  row.names = TRUE,
  sheetName = "net_sca",
  append = TRUE
)

names(cellchat_sca@net)
net_agg <- subsetCommunication(cellchat_agg,
  slot.name = "net",
  # sources.use = source_cells, targets.use = target_cells,
  thresh = 0.05
)
write.xlsx(cellchat_agg@net$count %>% as.data.frame(),
  file = "results/SA2021_AggSpa_cellchat_net.xlsx",
  row.names = TRUE,
  sheetName = "net_agg_count",
  append = TRUE
)
write.xlsx(cellchat_agg@net$weight %>% as.data.frame(),
  file = "results/SA2021_AggSpa_cellchat_net.xlsx",
  row.names = TRUE,
  sheetName = "net_agg_weight",
  append = TRUE
)
write.xlsx(net_agg %>% as.data.frame(),
  file = "results/SA2021_AggSpa_cellchat_net.xlsx",
  row.names = TRUE,
  sheetName = "net_agg",
  append = TRUE
)