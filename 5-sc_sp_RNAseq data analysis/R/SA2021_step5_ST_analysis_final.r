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

### Section 1: Data preprocessing for spatial (SP) transcriptomic data --------------------------------------------------

## load packages
library(httpgd)
library(Seurat)
# options(Seurat.object.assay.version = "v5")
# library(future)
library(patchwork)
library(ggplot2)
set.seed(123)
# workers <- 4
# plan("multisession", workers = workers)
# plan("sequential")


## Load data
root_path <- "/path_to_data/HCC-sp-RNAseq/STDATA_SA"
sub_dirs <- list.dirs(path = root_path, full.names = TRUE, recursive = FALSE)
sub_dirs <- sub_dirs[grepl("HCC-[0-9]T$|HCC-5",sub_dirs)]

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


### Section 2: ST data normalization --------------------------------------------------

# Normalize each Seurat object in spt_list using SCTransform
for (name in names(spt_list)) {
  cat("Processing:", name, "\n")
  tryCatch({
    spt_list[[name]] <- UpdateSeuratObject(spt_list[[name]])
    spt_list[[name]] <- SCTransform(spt_list[[name]], assay = "Spatial", verbose = F)
  }, error = function(e) {
    cat("Error processing", name, ":", e$message, "\n")
  })
}

# ! Save normalized data
saveRDS(spt_list,
    file = "05_normalized_spt_list.rds"
)
# spt_list <- readRDS("05_normalized_spt_list.rds")


### Section 3: Load in cytospace results --------------------------------------------------
library(readr)

# cytospace outcome dir
output_dir <- "/path_to_data/HCC-sp-RNAseq/cytospace/output"

# Function: Add meta.data to a given Seurat object
add_metadata_from_csv <- function(seurat_obj, seurat_name, output_dir) {
  # Construct the full path to the CSV file
  csv_file_path <- file.path(output_dir, paste0(seurat_name, "_cell_type_assignments_by_spot.csv"))

  # Check if the file exists
  if (!file.exists(csv_file_path)) {
    stop("File does not exist:", csv_file_path)
  }

  # Read the CSV file
  csv_data <- read_csv(csv_file_path, col_types = cols())
  csv_data$Mf_total <- csv_data$"Mph" + csv_data$"Mono-like"

  # Convert SpotID to factor for easier merging
  csv_data$SpotID <- as.factor(csv_data$SpotID)

  # Confirm if there's a column in the Seurat object's meta.data corresponding to SpotID
  # This example assumes Seurat object's cell names are the index of the meta.data
  if (!"cells" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj@meta.data$cells <- rownames(seurat_obj@meta.data)
  }

  # Merge the CSV data with the meta.data of the Seurat object
  # Use left_join to ensure all cells are retained in the Seurat object
  new_meta_data <- left_join(
    x = seurat_obj@meta.data,
    y = csv_data,
    by = c("cells" = "SpotID")
  ) 

  # Update the meta.data of the Seurat object
  seurat_obj@meta.data <- new_meta_data
  # restore meta.data rownames
  rownames(seurat_obj@meta.data) <- Cells(seurat_obj)

  # Return the updated Seurat object
  return(seurat_obj)
}

# Loop through spt_list to add meta.data to each Seurat object
for (seurat_name in names(spt_list)) {
  cat("Processing:", seurat_name, "\n")

  # Attempt to add meta.data, print error message if an error occurs
  tryCatch(
    {
      spt_list[[seurat_name]] <- add_metadata_from_csv(
        seurat_obj = spt_list[[seurat_name]],
        seurat_name = seurat_name,
        output_dir = output_dir
      )
    },
    error = function(e) {
      cat("Error processing", seurat_name, ":", e$message, "\n")
    }
  )
}



### Section 4: Gene set enrichment analysis --------------------------------------------------
# View(spt_list[["HCC-1T"]]@meta.data)

# define gene sets 
# exhausted T cell gene set from Zheng, C., et al. (2017). "Landscape of Infiltrating T Cells in Liver Cancer Revealed by Single-Cell Sequencing." Cell 169(7): 1342-1356 e1316.
Tex <- c("RAB27A", "HLA-DRA", "UBE2F-SCLY", "IGFLR1", "CXCL13", "IFNG", "UBE2F", "ITM2A", "ID3", "CD2BP2", "CHST12", "CTSD", "STAT3", "BST2", "CXCR6", "RALGDS", "VCAM1", "TRAFD1", "SYNGR2", "VAPA", "IFI35", "CD63", "NAB1", "PARK7", "MIR4632", "YARS", "PRKAR1A", "IL2RB", "DFNB31", "HMGN3", "MTHFD2", "NDFIP2", "PRF1", "PRDX5", "LAG3", "MS4A6A", "LINC00299")
Tc <- c("GZMA", "GZMB", "PRF1", "GZMK", "IL2", "IFNG", "FASLG", "CCL5", "NKG7", "CD28")

# Calculate T cell exhaustion and cytotoxicity score
# Loop over each Seurat object in the list
for (name in names(spt_list)) {
  cat("Processing:", name, "\n")
  tryCatch({
    # Calculate module score for Tex gene set
    spt_list[[name]] <- AddModuleScore(
      object = spt_list[[name]],
      features = list(Tex = Tex),  # Features should be a named list
      name = 'Tex_Score'  # Provide a new name for the module score column
    )

    # Calculate module score for Tc gene set
    spt_list[[name]] <- AddModuleScore(
      object = spt_list[[name]],
      features = list(Tc = Tc),  # Features should be a named list
      name = 'Tc_Score'  # Provide a new name for the module score column
    )
    
    spt_list[[name]]$Tex_Tc_ratio <-  spt_list[[name]]$Tex_Score1 / spt_list[[name]]$Tc_Score1

  }, error = function(e) {
    cat("Error processing", name, ":", e$message, "\n")
  })
}


# ! Save data
saveRDS(spt_list,
    file = "05_spt_list_scores.rds"
)
# spt_list <- readRDS("05_spt_list_scores.rds")


### Section 5: Plot  --------------------------------------------------
library(ggplot2)
library(cowplot)

colnames(spt_list[[1]]@meta.data)


# Directory to save the plots
output_dir <- "/path_to_data/HCC-sp-RNAseq/figures/ST_plots"

# calculate maximum and minimun 
min_mf <- min(sapply(spt_list, function(x) min(x$Mf_total, na.rm = TRUE)))
max_mf <- max(sapply(spt_list, function(x) max(x$Mf_total, na.rm = TRUE)))

min_tex <- min(sapply(spt_list, function(x) min(x$Tex_Score1, na.rm = TRUE)))
max_tex <- max(sapply(spt_list, function(x) max(x$Tex_Score1, na.rm = TRUE)))

# Function to plot spatial feature plots for given Seurat object with adjusted sizes
plot_spatial_features_adjusted <- function(seurat_obj) {

  # Generate SpatialFeaturePlot for each feature with adjusted sizes
  seurat_obj$Mf_total <- ifelse(is.na(seurat_obj$Mf_total), 0, seurat_obj$Mf_total)
  plot1 <- SpatialFeaturePlot(seurat_obj, features = "Mf_total", crop = TRUE, pt.size.factor = 1.6) +
    scale_fill_viridis_c(option = "D")

  plot2 <- SpatialFeaturePlot(seurat_obj, features = "Tex_Score1", crop = TRUE, pt.size.factor = 1.6, min.cutoff = "q50") +
    scale_fill_viridis_c(option = "D")

  plot3 <- SpatialFeaturePlot(seurat_obj, features = "Tc_Score1", crop = TRUE, pt.size.factor = 1.6, min.cutoff = "q1") +
    scale_fill_viridis_c(option = "D")

  plot4 <- SpatialFeaturePlot(seurat_obj, features = "Tex_Tc_ratio", crop = TRUE, pt.size.factor = 1.6, min.cutoff = "q10") +
    scale_fill_viridis_c(option = "D")
  plots <- list(plot1, plot2, plot3, plot4)
  # Combine plots into a single plot with adjusted sizes and reduced gaps
  combined_plot <- cowplot::plot_grid(
    plotlist = plots, ncol = 2,
    rel_widths = rep(1.2, length(plots)), rel_heights = rep(1.2, length(plots)), hjust = 0, vjust = 0
  )


  return(combined_plot)
}

# Loop through each Seurat object in the spt_list and generate/save the adjusted plots
for (sample_name in names(spt_list)) {
  cat("Plotting:", sample_name, "\n")
  tryCatch({
  # Generate the combined adjusted plot
  combined_plot <- plot_spatial_features_adjusted(spt_list[[sample_name]])
  
  # Save the plot as PDF and PNG
  pdf_file <- file.path(output_dir, paste0(sample_name, "_SpatialFeaturePlots.pdf"))
  png_file <- file.path(output_dir, paste0(sample_name, "_SpatialFeaturePlots.png"))
  
  ggsave(pdf_file, combined_plot, device = "pdf", width = 16, height = 12)
  ggsave(png_file, combined_plot, device = "png", width = 16, height = 12, dpi = 300)
  }, error = function(e) {
    cat("Error plotting", sample_name, ":", e$message, "\n")
  })
}



### Section 5: Correlation between Mf aggregation and T cell exhaustion --------------------------------------------------
library(dplyr)
library(ggplot2)
library(Seurat)

# 1. Add ST_sample column to each object's meta.data and combine them
all_metadata <- lapply(names(spt_list), function(sample_name) {
  # Add the ST_sample column
  spt_list[[sample_name]]@meta.data$ST_sample <- sample_name

  # Return the modified metadata
  return(spt_list[[sample_name]]@meta.data)
}) %>% bind_rows()

# 2. Plotting the scatter plot with Mf_total on x-axis and Tex_Score1 on y-axis
p <- ggplot(all_metadata, aes(x = as.factor(Mf_total), y = Tex_Score1)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +  # 添加盒形图，不显示异常点
  geom_point(aes(color = as.factor(case_when(
    Mf_total == 0 ~ "Black",
    Mf_total == 1 ~ "Blue",
    Mf_total > 1 & Mf_total < 4 ~ "Purple",
    Mf_total >= 4 ~ "Red"
  ))), position = position_jitter(width = 0.2), alpha = 0.6) +  # 添加散点图
  scale_color_manual(
    values = c("Black" = "black", "Blue" = "blue", "Purple" = "purple", "Red" = "red"),
    name = "Mf_total Category"
  ) +
  theme_bw() +
  labs(
    x = "Macrophage number per spot", 
    y = "T cell exhaustion score", 
    title = "Association between Macrophage aggregation and T cell exhaustion"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )
p

write.csv(all_metadata,
file = "results/SA2021_05_ST_spots_Association_Mf_Tex.csv")

# Save the plot
pdf_file <- file.path(output_dir, paste0("Association_Mf_Tex_dotplot.pdf"))
png_file <- file.path(output_dir, paste0("Association_Mf_Tex_dotplot.pdf"))

ggsave(pdf_file, p, device = "pdf", width = 16, height = 12)
ggsave(png_file, p, device = "png", width = 16, height = 12, dpi = 300)