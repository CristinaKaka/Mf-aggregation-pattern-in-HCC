# packages loading
library(dplyr)

# set work directory
setwd("...")

#!! import data 
patdata <- read.csv("results/patientdata.csv",row.names = 1)

# 02-patient proportion data -------------------------------------------------
# 02-change to longer table
library(reshape2)
library(ggplot2)

names(patdata)
col_names <- c("patient_id","CD8G3_byDFS",
  "mono.mf_inIC", "Neutrophil_inIC", "DC_inIC", "mast.cell_inIC", "CD8T_inIC",
  "CD3T_inIC", "NK_inIC", "B.cell_inIC", "Treg_inIC", "other_inIC"
)

propdata <- patdata %>% 
  select(col_names) %>% 
  melt( id.vars = c("patient_id","CD8G3_byDFS"),variable.name = "celltype", value.name = "Freq") %>% 
  mutate(celltype = case_when(
    celltype == "mono.mf_inIC" ~ "Mono/mf",
    celltype == "Neutrophil_inIC" ~ "Neutrophil",
    celltype == "DC_inIC" ~ "DC",
    celltype == "mast.cell_inIC" ~ "Mast cell",
    celltype == "CD8T_inIC" ~ "CTL",
    celltype == "CD3T_inIC" ~ "Th",
    celltype == "NK_inIC" ~ "NK",
    celltype == "B.cell_inIC" ~ "B cell",
    celltype == "Treg_inIC" ~ "Treg",
    celltype == "other_inIC" ~ "other",
    TRUE ~ celltype
  ))


# as.factor
propdata$celltype <- factor(propdata$celltype,levels = c("Mono/mf", "Neutrophil", "DC", "Mast cell", "CTL", "Th", "NK", "B cell", "Treg", "other"))
propdata$CD8G3_byDFS <- factor(propdata$CD8G3_byDFS,levels = c("0","1","2"))

# order the plot
patdata <- patdata[order(patdata$patient_id), ]
patdata <- patdata[order(patdata$CD8G3_byDFS), ]
propdata$patient_id <- factor(propdata$patient_id, levels = patdata$patient_id)


k_pal <-  c("#f21b3f","#ffca3a","#ff5129","#ea6aa3",
            "#33a1fd","#3dc6c9","#29bf12","#b7e4c7","#9655cc","#ccb9de")

# plot by patients
proplot <- ggplot(propdata, aes_string(y = "Freq")) + 
  labs(x = NULL,  y = "Proportion [%]") + theme_bw() + 
  theme(panel.grid = element_blank(),  strip.text = element_text(face = "bold"), strip.background = element_rect(fill = NA,  color = NA), axis.text = element_text(color = "black"), 
                                                                                          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
                                                                                          legend.key.height = unit(0.8, "lines"))
proplot <- proplot+#facet_wrap(~CD8G3_byDFS, scales = "free_x",) + 
  geom_bar(aes_string(x = "patient_id",  fill = "celltype"), position = "fill", stat = "identity") + 
         scale_fill_manual("celltype", values = k_pal) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) +
  theme(panel.border = element_blank(), panel.spacing.x = unit(1, "lines"))

proplot

ggsave(filename = "figures/celltype_patients.pdf",proplot,width = 8,height = 3)
