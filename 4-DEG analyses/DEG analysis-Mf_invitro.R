
# Section 1: data preparation ---------------------------------------------

# set work directory
setwd('D:/Clustered Mf/WYL/Spatial Analysis/Qupath/18-MRS&macrophage cluster/differentially expressed genes v1/Mf in vitro')

# install packages
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#  options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# BiocManager::install(version = '3.16',lib = '...')
# BiocManager::install(c("rnaseqGene","glmpca","sva",
#                        "DEGreport","rnaseqDTU","msigdbr","PoiClaClu"),lib = '...')
# 
# BiocManager::install(c("tximeta"),lib = '...')
# BiocManager::install(c("airway"),lib = '...')
# BiocManager::install(c("DESeq2"),lib = '...')

# packages loading
library("tximeta")
library("airway")

# loading data
load("20221108_Mf_RNAseq_Part_1.RData")

dim(se)
## [1] 245261     12
head(rownames(se))

# setTximetaBFC('E:/software/R-4.2.2/BiocFileCache')
gse <- summarizeToGene(se)
dim(gse)

## [1] 61125    12 
head(rownames(gse))

# data(gse)

# data frame viewing
assayNames(gse)
rowRanges(gse)
seqinfo(rowRanges(gse))

# annotate groups
groupdata <- as.data.frame(colData(gse))
group <- c('clustered','sparse')
groupdata <- cbind.data.frame(groupdata,group)

gse$condition <- as.factor(groupdata$group)
gse$donor <- as.factor(substr(groupdata$names,0,1))
colData(gse)

# check the millions of fragments that could be mapped by Salmon 
# to the genes (the second argument of round tells how many decimal points to keep)
round( colSums(assay(gse)) / 1e6, 1 )


# construct a DESeqDataSet object
library("DESeq2")

dds <- DESeqDataSet(gse, design = ~ condition+donor)

# Section 2: Export gene expression data -------------------------

# export gene expression level (TPM) 
TPM <- gse@assays@data$abundance
library("AnnotationDbi")
library("org.Hs.eg.db")
ens_TPM <- substr(rownames(TPM), 1, 15) # note that the rownames should be shorten before mapping
sym_TPM <- mapIds(org.Hs.eg.db,
                  keys=ens_TPM,
                  column="SYMBOL",
                  keytype="ENSEMBL",
                  org.Hs.eg.db="first")


library(tidyverse)
library(tibble)
TPM <- TPM %>%
  as_tibble() %>%
  mutate(GeneSymbol = sym_TPM) %>%
  mutate(Ensembl = ens_TPM) %>%
  filter(!is.na(GeneSymbol)) %>%
  dplyr::select(Ensembl, GeneSymbol, everything()) %>%
  arrange(Ensembl)
write.csv(TPM, file = "20221109-Mf_in_vitro_TPM.csv")

# export gene expression level (normalized counts) 
counts <- counts(dds)

library("AnnotationDbi")
library("org.Hs.eg.db")
## note that the rownames should be shorten before mapping
ens_counts <- substr(rownames(counts), 1, 15)
sym_counts <- mapIds(org.Hs.eg.db,
                     keys=ens_counts,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     org.Hs.eg.db="first")
library(tidyverse)
library(tibble)
counts <- counts %>%
  as_tibble() %>%
  mutate(GeneSymbol = sym_counts) %>%
  mutate(Ensembl = ens_counts) %>%
  filter(!is.na(GeneSymbol)) %>%
  dplyr::select(Ensembl, GeneSymbol, everything()) %>%
  arrange(Ensembl)
# View(counts)
## export data
write.csv(counts, file = "20221109-Mf_in_vitro_counts.csv")


# Section 3: Exploratory analysis and visualization ----------------------------------

# pre-filtering the dataset
nrow(dds)
## [1] 61125

keep <- rowSums(counts(dds)) > 10
dds <- dds[keep,]
nrow(dds)

# The variance stabilizing transformation and the rlog
# DESeq2 offers two transformations for count data that stabilize the variance across the mean: the variance stabilizing transformation (VST) for negative binomial data with a dispersion-mean trend
# The rlog tends to work well on small datasets (n < 30)
rld <- rlog(dds, blind = FALSE)##blind = T, which means that differences the variables in the design will contribute to the expected variance-mean trend of the experiment
head(assay(rld), 3)

## VST is recommended for medium-to-large datasets (n > 30)
vsd <- vst(dds, blind = T) 
head(assay(vsd), 3)

# show the effect of the transformation
library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:3]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

# Scatterplot of transformed counts from two samples.
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation) 

# Sample distances
sampleDists <- dist(t(assay(rld)))
sampleDists

# visualize the distances in a heatmap
library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$condition, rld$donor, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

## Another option for calculating sample distances is to use the Poisson Distance
## This measure of dissimilarity between counts also takes the inherent variance structure of counts into consideration when calculating the distances between samples
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
## plot the heatmap of Poisson Distance
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- dds$names
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

# visualize sample-to-sample distances by principal components analysis (PCA)
plotPCA(rld, intgroup = c("condition", "donor"))

pcaData <- plotPCA(rld, intgroup = c("condition", "donor"), returnData = TRUE)
pcaData

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape = donor)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")

# PCA plot using Generalized PCA
library("glmpca")
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$condition <- dds$condition
gpca.dat$donor <- dds$donor

ggplot(gpca.dat, aes(x = dim1, y = dim2, color = condition, shape = donor)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")


# MDS plot
mds <- as.data.frame(colData(rld))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = condition, shape = donor)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with rlog data")

# MDS plot using the Poisson Distance.
mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = condition, shape = donor)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")



# Section 4: Formal Differential expression analysis -----------------------------
#using SVA with DESeq2 to Removing hidden batch effects

# construct a DESeqDataSet object
library("DESeq2")
library('sva')

dds <- DESeqDataSet(gse, design = ~ condition+donor)

# pre-filtering the dataset
nrow(dds)
## [1] 61125

keep <- rowSums(counts(dds)) > 10
dds <- dds[keep,]
nrow(dds)

dds <- estimateSizeFactors(dds)
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ condition+donor, colData(dds))
mod0 <- model.matrix(~   donor, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv = 1)

# Surrogate variables 1  plotted over donor.
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:1) {
  stripchart(svseq$sv[, i] ~ dds$donor, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
}

# redesign formula
ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
design(ddssva) <- ~ SV1 + donor + condition


# Running the differential expression pipeline
ddssva <- DESeq(ddssva)

# Building the results table
res <- results(ddssva, contrast=c("condition","clustered","sparse"),
               alpha = 0.05)
table(res$padj < 0.05)
summary(res)

## annotate results
library("AnnotationDbi")
library("org.Hs.eg.db")
## note that the rownames should be shorten before mapping
ens.res <- substr(rownames(res), 1, 15)
res_output <- as.data.frame(res)
res_output$symbol <- mapIds(org.Hs.eg.db,
                         keys=ens.res,
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         org.Hs.eg.db="first")
res_output$entrez <- mapIds(org.Hs.eg.db,
                         keys=ens.res,
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         org.Hs.eg.db="first")
dim(res_output)
res_output <- res_output %>%
  filter(!is.na(padj))
dim(res_output)
res_output$Group <- NA
res_output$Group <- ifelse(res_output$padj<0.05 & abs(res_output$log2FoldChange)>1,"Sig","notSig")
res_output$Group[res_output$padj<0.05 & abs(res_output$log2FoldChange)>1] <- ifelse(res_output$log2FoldChange[res_output$padj<0.05 & abs(res_output$log2FoldChange)>1]>0,"Up","Dn")
res_output$Group[res_output$padj<0.01 & abs(res_output$log2FoldChange)>1] <- ifelse(res_output$log2FoldChange[res_output$padj<0.01 & abs(res_output$log2FoldChange)>1]>0,"Upstr","Dnstr")
res_output$Group <- factor(res_output$Group, levels=c("notSig","Dn","Up","Dnstr","Upstr"))
res_output <- res_output[order(res_output$padj),]
## export results
write.csv(res_output, 
          file = "20221109-Mf_in_vitro_DEG.csv")

# save environment
save.image(file = "20221109-Mf_invitro.Rdata")

rm(list = ls())

# Section 5: Volcanol plot ---------------------------------------------
# options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# BiocManager::install("EnhancedVolcano")
library('EnhancedVolcano')
library('dplyr')

res <- read.csv('20221109-Mf_in_vitro_DEG.csv',row.names = 1)
colnames(res)
dim(res)
res <- res %>%
  filter(!is.na(symbol)) 

res <- res%>%
  aggregate(by = list(res$symbol),FUN = max) %>%
  arrange(padj)

rownames(res) <- res$Group.1
res$Group.1 <- NULL
res$FC <- 2^(res$log2FoldChange)

# select genes to label-20221109
# selectgene <- c("IL7R","CD163","MAF","CCL13","CA12",
#                 "FCGR2A","FCGR2B","","SLC40A1","VCAN","TREM2",
#                 "TNF","IL6","CCL3","CCL4","CCL13")

# 20240102
selectgene <- c(
  # lipid metabolism
  "ABCG1", "APOC1", "TREM2","CD163", "MRC1",
  # tissue remodeling
  "CCL18", "VEGFA", "MMP9", "MMP19", "CD163",
  # fetal-like
  "FOLR2", "HES1", "C1QA",
  # pro-tumoral TAM phenotype
  "SLC40A1", "CA12",
  # overlap with LA-TAM
  # "CD163", "FOLR2", "MRC1", "MMP9", "C1QA","SLC40A1", "TREM2", "APOC", "RNASE1",
  # high expression and in vivo upregulation
  "CCL13", "CA12", "IL7R",
  # TF of LA-TAM 
  "MAF",
  # Mf-THBS1
  # "VCAN", "S100A8", "S100A9", "S10012",
  # "ID3", "RUNX2",
  # inflammatory cytokine
  "TNF", "IL6",
  "CCL3","CCL4"
)

#
p1 <- EnhancedVolcano(res,
    lab              = rownames(res),
    x                = "log2FoldChange",
    y                = "padj",
    xlim             = c(-5,6),
    selectLab        = selectgene,
  # lab              = NULL,
    xlab             = bquote(~ Log[2] ~ "fold change"),
    pCutoff          = 0.05,
    FCcutoff         = 1,
    pointSize        = 2,
    labSize          = 5,
  # labCol           = 'black',
  # labFace          = 'bold',
    boxedLabels      = T,
    colAlpha         = 3 / 5,
    legendPosition   = "right",
    legendLabSize    = 6,
    legendIconSize   = 4.0,
    drawConnectors   = TRUE,
    arrowheads       = F,
    max.overlaps     = 10000000,
    widthConnectors  = 1,
    lengthConnectors = unit(1, "npc"),
  col = c("black", "black", "black", "red3"),
  # colConnectors = 'black'
)
# + coord_flip()

p1
help(EnhancedVolcano)
ggsave(file = 'Plots/20240102-DEG volcano plot.png', plot = p1, width = 10, height = 8)
ggsave(file = 'Plots/20240102-DEG volcano plot.pdf', plot = p1, width = 10, height = 8)


# Section 6: go term analysis ---------------------------------------------
options(stringsAsFactors = F) 
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(org.Hs.eg.db)
library(DOSE)
# BiocManager::install("pathview",lib = 'E:/software/R-4.2.2/library',ask = F,update = F)
library(pathview)

# select genes for go term analysis
log2FC_cutoff = log2(2)
pvalue_cutoff = 0.01
padj_cutoff = 0.01

# assessing upregulated or downregulated genes
colnames(res)
need_DEG <- res[,c("log2FoldChange",'pvalue','padj')]  

gene_up=rownames(need_DEG[with(need_DEG,log2FoldChange>log2FC_cutoff & pvalue<pvalue_cutoff & padj<padj_cutoff),])
gene_down=rownames(need_DEG[with(need_DEG,log2FoldChange < -log2FC_cutoff & pvalue<pvalue_cutoff & padj<padj_cutoff),])

# transfer symbol to entrezID
gene_up_entrez <- as.character(na.omit(bitr(gene_up,
                                            fromType="SYMBOL",
                                            toType="ENTREZID", 
                                            OrgDb="org.Hs.eg.db")[,2])) #"org.Hs.eg.db" "org.Mm.eg.db"
gene_down_entrez <- as.character(na.omit(bitr(gene_down,
                                              fromType="SYMBOL", 
                                              toType="ENTREZID", #
                                              OrgDb="org.Hs.eg.db")[,2])) #"org.Hs.eg.db" "org.Mm.eg.db"

gene_diff_entrez <- unique(c(gene_up_entrez ,gene_down_entrez))




# Section 7: KEGG GO enrichment ------------------------------------------------------


kegg_enrich_results <- enrichKEGG(gene  = gene_up_entrez,
                                  organism  = "hsa", #human hsa mous mmu
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.2)

kegg_enrich_results <- DOSE::setReadable(kegg_enrich_results, 
                                         OrgDb="org.Hs.eg.db", 
                                         keyType='ENTREZID')#ENTREZID to gene Symbol


write.csv(kegg_enrich_results@result,'20221110-KEGG_gene_up_enrichresults.csv') 
# save(kegg_enrich_results, file ='KEGG_gene_up_enrichresults.Rdata')

# GO ontology analysis bp
go_enrich_results <- enrichGO(gene = gene_up_entrez,
                              OrgDb = "org.Hs.eg.db",
                              ont   = "BP"  ,     #One of "BP", "MF"  "CC"  "ALL" 
                              pvalueCutoff  = 0.05,
                              qvalueCutoff  = 0.2,
                              readable      = TRUE)
write.csv(go_enrich_results@result, '20221110-GO_gene_up_BP_enrichresults.csv') 
# save(go_enrich_results, file ='GO_gene_up_enrichresults.Rdata')

go_down_results <- enrichGO(gene = gene_down_entrez,
                              OrgDb = "org.Hs.eg.db",
                              ont   = "BP"  ,     #One of "BP", "MF"  "CC"  "ALL" 
                              pvalueCutoff  = 0.05,
                              qvalueCutoff  = 0.2,
                              readable      = TRUE)
write.csv(go_down_results@result, '20221110-GO_gene_down_BP_enrichresults.csv')



# GO ontology analysis bp
go_enrich_results <- enrichGO(gene = gene_up_entrez,
                              OrgDb = "org.Hs.eg.db",
                              ont   = "MF"  ,     #One of "BP", "MF"  "CC"  "ALL" 
                              pvalueCutoff  = 0.05,
                              qvalueCutoff  = 0.2,
                              readable      = TRUE)
write.csv(go_enrich_results@result, '20221110-GO_gene_up_MF_enrichresults.csv') 
# save(go_enrich_results, file ='GO_gene_up_enrichresults.Rdata')

go_down_results <- enrichGO(gene = gene_down_entrez,
                            OrgDb = "org.Hs.eg.db",
                            ont   = "MF"  ,     #One of "BP", "MF"  "CC"  "ALL" 
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.2,
                            readable      = TRUE)
write.csv(go_down_results@result, '20221110-GO_gene_down_MF_enrichresults.csv') 


# Section 8: GO analysis plot ---------------------------------------------
library(ggplot2) 
library(ggrepel)

go_enrich_30 <- data.frame(go_enrich_results) %>%
  filter(Count>=30) 

go_enrich_30 <- go_enrich_30[1:10,]
go_enrich_30$number <- factor(go_enrich_30$Description,levels = rev(go_enrich_30$Description))

gotermplot <- ggplot(data = go_enrich_30, 
       aes(x = Count,y = reorder(number,Count)))+ 
  geom_point(aes(size = Count,color = p.adjust))+ 
  theme_bw()+ # 去除背景色
  scale_colour_gradient(low = "red",high = "blue")+ 
  labs(x = "Gene Number", y = "",title = "Gene ontologies", 
       color = expression(p.adjust),size = "Count")+ 
  scale_y_discrete(labels = function(x) str_wrap(x,width = 43))+
  theme(axis.title = element_text(size = 8),axis.text = element_text(size = 8),
        plot.title = element_text(size = 8,hjust = 0.5,face = "bold"),
        legend.title = element_text(size = 8),legend.text = element_text(size = 8),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))

ggsave(gotermplot,filename = "Plots/20221109-gene ontologies-dotplot.pdf",height = 3.3, width =4.9, units = "in")

