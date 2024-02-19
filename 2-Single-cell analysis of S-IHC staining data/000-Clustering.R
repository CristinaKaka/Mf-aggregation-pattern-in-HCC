# package preparation -----------------------------------------------------

#if(!require('devtools')) {install.packages('devtools')}
#library('devtools')
#install_github("immunedynamics/spectre")
# 
# library("Spectre")
# ## Check if all required packages have been installed
# Spectre::package.check(type = 'spatial')
# 
# ## Load all required packages
# Spectre::package.load(type = 'spatial')
library(data.table)
library(HDCytoData)
library(data.table)
library(dplyr)
library(flowCore)

# data preparation --------------------------------------------------------
setwd('...')

# stratify intensity  ------------------------------------------------------

## load cell data: cell.dat
load(file = "..")

head(cell.dat)


# stratify stain intensity into 4 categories(<=0.05,0.05< x <=0.1,0.1< x <0.3,>=0.3)
cell.dat.cat <- cell.dat[,c("CD68","CD14","CD45""CD15","S100","CD57","FOXP3","CD8","Tryptase", "CD3", "CD20")] 
                          
cell.dat.cat<- lapply(cell.dat.cat,cut,breaks = c(-Inf,0.05,0.1,0.3,Inf),right = T,labels = c(0,1,2,3))
cell.dat.cat <- as.data.frame(cell.dat.cat)
View(head(cell.dat.cat))
names(cell.dat.cat)

names(cell.dat.cat) <- paste(names(cell.dat.cat),"_G",sep = "")
cell.dat <- cbind.data.frame(cell.dat,cell.dat.cat)

# adding group subset immune cell
cell.dat <- cell.dat %>%
    mutate(im_G = {
        x <- rep("Non-immune", nrow(cell.dat))
        x[
            LAMP3_G %in% c(2, 3) |
                CD204_G %in% c(2, 3) |
                CD68_G %in% c(2, 3) |
                CD163_G %in% c(2, 3) |
                CD169_G %in% c(2, 3) |
                CD45_G %in% c(2, 3) |
                CD15_G %in% c(2, 3) |
                S100_G %in% c(2, 3) |
                CD57_G %in% c(2, 3) |
                CD11c_G %in% c(2, 3) |
                FOXP3_G %in% c(1, 2, 3) |
                CD8_G %in% c(2, 3) |
                Tryptase_G %in% c(2, 3) |
                CD3_G %in% c(2, 3) |
                CD20_G %in% c(2, 3)
        ] <- "iMM"
        x
    })

expr <- data.frame(t(cell.dat[,c("CD68","CD14","CD45","CD15","S100","CD57","FOXP3","CD8","Tryptase","CD3","CD20",
                                 "CD68_G","CD14_G","CD45_G","CD15_G","S100_G","CD57_G","FOXP3_G","CD8_G","Tryptase_G",
                                 "CD3_G","CD20_G")]))

View(expr[,1:10])
expr <- as.data.frame(lapply(expr, as.numeric))

# cellID to colnames
colnames(expr) <- cell.dat$CellID

row.names(expr)
antigen <- c(
    "CD68", "CD14", "CD45", "CD15", "S100", "CD57", "FOXP3", "CD8", "Tryptase", "CD3", "CD20",
    "CD68_G", "CD14_G", "CD45_G", "CD15_G", "S100_G", "CD57_G", "FOXP3_G", "CD8_G", "Tryptase_G",
    "CD3_G", "CD20_G"
)

rownames(expr) <- antigen
ihc_colname <- antigen
marker_class <- rep("type", length(antigen))

length(ihc_colname)
length(antigen)
length(marker_class)
panel <- cbind.data.frame(ihc_colname,antigen,marker_class)
rownames(panel) <- panel$ihc_colname


# loading metadata 
names(cell.dat)
md <- cell.dat[,c("sample_id","CellID","TMA_ID","x_x","y_y","CoreID","tissue_type","patient_id",
                  "im_G")]
View(md[1:10,])
dim(md)

# specify levels for conditions & sample IDs to assure desired ordering
md$tissue_type <- factor(md$tissue_type)
md$sample_id <- factor(md$sample_id)
rownames(md) <- md$CellID


# spot check that all panel columns are in the expr object
all(panel$ihc_colname %in% rownames(expr))
all(rownames(md) %in% colnames(expr))
dim(panel)
dim(md)
dim(expr)

# data subset and exploration --------------------------------------------------------------

## loading packages
library(SummarizedExperiment)
library(SingleCellExperiment)
library(CATALYST)
library(diffcyt)
#BiocManager::install('diffcyt')

# set work directory
setwd('...')

# Create sumarizedexperiment object
sce <- SummarizedExperiment(assays = list(exprs=expr),colData = md,rowData = panel)

# transform to singlecellexperiment
sce <- as(sce, "SingleCellExperiment")
n_cells(sce)


# select immune cell data
sce_IM <- sce[,sce$im_G == "iMM"]
# rowData(sce1)


## Such plots show similarities between samples measured in an unsupervised way 
## and give a sense of how much differential expression can be detected before conducting any formal tests

plotCounts(sce_IM, group_by = "sample_id")
pbMDS(sce_IM, label_by = "sample_id")

## plot each patient with each marker expression levels
plotExprHeatmap(sce_IM,
    features = c(
        "CD68", "CD14", "CD45", "CD15", "S100", "CD57", "FOXP3", "CD8", "Tryptase", "CD3", "CD20",
        "CD68_G", "CD14_G", "CD45_G", "CD15_G", "S100_G", "CD57_G", "FOXP3_G", "CD8_G", "Tryptase_G",
        "CD3_G", "CD20_G"
    ),
    scale = "last",
    hm_pal = rev(hcl.colors(10, "YlGnBu"))
)


# Marker ranking based on the non-redundancy score(seem to be a socre of expression levels)
# Markers with higher score explain a larger portion of variability present in a given sample
plotNRS(sce_IM, features = c(
    "CD68", "CD14", "CD45", "CD15", "S100", "CD57", "FOXP3", "CD8", "Tryptase", "CD3", "CD20",
    "CD68_G", "CD14_G", "CD45_G", "CD15_G", "S100_G", "CD57_G", "FOXP3_G", "CD8_G", "Tryptase_G",
    "CD3_G", "CD20_G"
), color_by = "condition")



# Clustering --------------------------------------------------------------

setwd("...")

set.seed(4399)
sce_IM <- cluster(sce_IM, features = c("CD68_G","CD14_G","CD45_G",
                                       "CD15_G","S100_G","CD57_G","FOXP3_G","CD8_G","Tryptase_G",
                                       "CD3_G","CD20_G"),
                  xdim = 10, ydim = 10, maxK = 60, seed = 4399)

plotExprHeatmap(sce_IM, features = c("CD68","CD14","CD45",
                                     "CD15","S100","CD57","FOXP3","CD8","Tryptase",
                                     "CD3","CD20"), 
                by = "cluster_id", k = "meta60", 
                bars = TRUE, perc = TRUE,
                scale="never")

plotExprHeatmap(sce_IM, features = c("CD68","CD14","CD45",
                                     "CD15","S100","CD57","FOXP3","CD8","Tryptase",
                                     "CD3","CD20"), 
                by = "cluster_id", k = "som100", 
                bars = TRUE, perc = TRUE,
                scale="never")

plotExprHeatmap(sce_IM, features = c("CD68_G","CD14_G","CD45_G",
                                     "CD15_G","S100_G","CD57_G","FOXP3_G","CD8_G","Tryptase_G",
                                     "CD3_G","CD20_G"), 
                by = "cluster_id", k = "meta20", 
                row_clust = T,col_clust = T,
                hm_pal = rev(hcl.colors(10, "YlGnBu")),
                bars = TRUE, perc = TRUE,
                scale="never")

plotExprHeatmap(sce_IM, features = c("CD68","CD14","CD45",
                                     "CD15","S100","CD57","FOXP3","CD8","Tryptase",
                                     "CD3","CD20"), 
                by = "cluster_id", k = "meta20", 
                bars = TRUE, perc = TRUE,
                scale="never")



# cell subsets verification -----------------------------------------------

# cell subset
library(dplyr)

# meta60--cluster_id
a <- as.data.frame(cluster_codes(sce_IM))
a <- a[,c("som100","meta60")]

a <- a$som100[a$meta60=="28"]

cells1 <- as.data.frame(colData(sce_IM[,sce_IM$cluster_id %in% a]))

View(cells1)

# merge cluster -----------------------------------------------------------
library(readxl)

merging_table1 <- "20230116-cluster60.xlsx"
merging_table1 <- read_excel(merging_table1)
head(data.frame(merging_table1))
colnames(merging_table1)[1:2] <- c("original_cluster","new_cluster")

# convert to factor with merged clusters in desired order
unique(merging_table1$new_cluster)

merging_table1$new_cluster <- factor(merging_table1$new_cluster, 
                                     levels = c("CD3T","B cell","CD8T","DC","Neutrophil",
                                                "other","mono/mf","mast cell","Treg","NK"))
# apply manual merging
sce_IMme <- mergeClusters(sce_IM, k = "meta60", 
                      table = merging_table1[,1:2], id = "merging1")


plotExprHeatmap(sce_IMme, features = c("LAMP3","CD68","CD14","CD45",
                                   "CD15","S100","CD57","FOXP3","CD8","Tryptase",
                                   "CD3","CD20"), 
                by = "cluster_id", k = "meta60",m="merging1",
                scale = 'last',distance = "euclidean")



# export data after merging clusters
saveRDS(sce_IMme,file = '...')


