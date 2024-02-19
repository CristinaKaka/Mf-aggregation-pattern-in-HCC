
# Section 01: data preparation --------------------------------------------------  
# packages loading
library(SingleCellExperiment)

# # loading single cell data after clustering
# sce <- readRDS("results/20231203-afterTSNE_UMAP.rds")
# 
# # cell data for used
# celldata <- as.data.frame(colData(sce))
# colnames(celldata)

# loading cell coordinate
sce <- readRDS("data/sce_iMM_patInf.rds")

loc <- as.data.frame(colData(sce))
# head(loc)

celldata <- loc
# colnames(celldata)
celldata$phenotype <- factor(celldata$phenotype,levels = c("mono/mf", "Neutrophil", "DC", "mast cell", "CD8T", "CD3T", "NK", "B cell", "Treg", "other","Non-immune cell"))


### Section 02: plot cell type in situ --------------------------------------------------

#### 02-plot cell type in situ-1310406 -----------------------------------------------

# packages loading
library(ggplot2)
library(scales)
library(RColorBrewer)
library(dplyr)

# select patient
dataplot1 <- celldata %>%
  dplyr::filter(sample_id == '1310406') %>%
  dplyr::filter(phenotype !='Non-immune cell')
# View(dataplot1)

# dataplot1$phenotype
# display.brewer.all(colorblindFriendly = TRUE)
# display.brewer.pal(5,"Set4")

fig1 <- dataplot1 %>%
  ggplot(aes(x=x_x,y=y_y,colour=phenotype,fill=phenotype))+
  geom_point(alpha=0.9,size=1)+
  # scale_color_brewer(palette = "Set2")+
  # scale_fill_brewer(palette = 'Set2') +
  ggtitle("1310406")+
  scale_color_manual(values =c("#f21b3f","#ffca3a","#ff5129","#ea6aa3",
                               "#33a1fd","#3dc6c9","#29bf12","#b7e4c7","#9655cc","#ccb9de"),
                     breaks = levels(dataplot1$phenotype))+
  # scale_x_continuous(limits = c(3500, 4800))+
  # scale_y_continuous(limits = c(11100, 15000))+
  # scale_y_reverse(limits = c(12300, 11100)) +
  theme(panel.background = element_rect(fill = "black"))+
  # theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),panel.border = element_blank(),panel.grid = element_blank())+
  coord_fixed()
# scale_x_continuous(breaks = seq(1,7,0.5))+
# scale_y_continuous(breaks = seq(4,8,0.5))

fig1

#### 02-plot cell type in situ-1311008 -----------------------------------------------------------------

# select patient
dataplot4 <- celldata %>%
  dplyr::filter(sample_id == '1311008') %>% 
  dplyr::filter(phenotype !='Non-immune cell')


# display.brewer.all(colorblindFriendly = TRUE)
# display.brewer.pal(5,"Set4")

fig4 <- dataplot4 %>%
  ggplot(aes(x=x_x,y=y_y,colour=phenotype,fill=phenotype))+
  geom_point(alpha=0.9,size=1)+
  # scale_color_brewer(palette = "Set2")+
  # scale_fill_brewer(palette = 'Set2') +
  ggtitle("1311008")+
  scale_color_manual(values =c("#f21b3f","#ffca3a","#ff5129","#ea6aa3",
                               "#33a1fd","#3dc6c9","#29bf12","#b7e4c7","#9655cc","#ccb9de"),
                     breaks = levels(dataplot4$phenotype))+
  # scale_x_continuous(limits = c(3500, 4800))+
  # scale_y_continuous(limits = c(11100, 15000))+
  # scale_y_reverse(limits = c(12300, 11100)) +
  theme(panel.background = element_rect(fill = "black"))+
  # theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),panel.border = element_blank(),panel.grid = element_blank())+
  coord_fixed()
# scale_x_continuous(breaks = seq(1,7,0.5))+
# scale_y_continuous(breaks = seq(4,8,0.5))

fig4


### Section 03: contour plot mf --------------------------------------------------

# packages loading
library(viridis)
library(rayshader)
library(patchwork)
source('source/render_snapshot.R')

# devtools::install_github("tylermorganwall/rayshader")

densityplot <- celldata %>%
  filter(sample_id %in% c( "1310406")) %>%
  filter(phenotype == "mono/mf")

aggregate <- ggplot(densityplot, aes(x = x_x, y = y_y)) +
  # raster plot
  layer(
    stat     = "density_2d",
    geom     = "raster",
    mapping  = aes(fill = after_stat(count) * 200),
    params   = list(contour = FALSE, interpolate = T, h = 100, alpha = 0.9),
    position = "identity"
  ) +
  # group by phenotype
  # facet_wrap(phenotype ~ .) +
    # scale color
    # scale_fill_viridis(
    #     option   = "A",
    #     # limits   = c(0.25, 1),
    #     na.value = "transparent"
    #   # ,         oob         = scales::squish
    # ) +
    scale_fill_gradient(
    low    = "pink",
    high   = "#f21b3f",
    limits = c(0.2, 1.2)
    # ,      oob        = scales::squish
    ,       na.value   = "transparent"
    ) +
    theme_bw()+
    theme(
      plot.title = element_blank(),
      axis.title=element_blank(),
      panel.grid      = element_blank(),
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      panel.border    = element_blank(),
      legend.position = "none",

    ) +
    coord_fixed()



plot_gg(aggregate,multicore=TRUE,width=5,height=5,scale=300,shadow_intensity = 1
, windowsize = c(800, 800)
, soliddepth =0
# ,           alpha              = 0.5
# ,           save_shadow_matrix = T
 )

# # capture the 3d plot
render_camera(zoom = 0.7, theta = 10, phi = 30)
Sys.sleep(0.2)
render_snapshot(clear = F)





densityplot <- celldata %>%
  filter(sample_id %in% c( "1311008")) %>%
  filter(phenotype == "mono/mf")

scatter <- ggplot(densityplot, aes(x = x_x, y = y_y)) +
  # raster plot
  layer(
    stat     = "density_2d",
    geom     = "raster",
    mapping  = aes(fill = after_stat(count) * 200),
    params   = list(contour = FALSE, interpolate = T, h = 100, alpha = 0.9),
    position = "identity"
  ) +
  # group by phenotype
  # facet_wrap(phenotype ~ .) +
    # scale color
    # scale_fill_viridis(
    #     option   = "A",
    #     # limits   = c(0.25, 1),
    #     na.value = "transparent"
    #   # ,         oob         = scales::squish
    # ) +
    scale_fill_gradient(
      low       = "pink",
      high      = "#f21b3f",
      limits    = c(0.2, 1.2)
    #,      oob = scales::squish
    , na.value  = "transparent"
    ) +
    theme_bw()+
    theme(
      plot.title = element_blank(),
      axis.title=element_blank(),
      panel.grid      = element_blank(),
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      panel.border    = element_blank(),
      legend.position = "none",

    ) +
    coord_fixed()


plot_gg(scatter,multicore=TRUE,width=5,height=5,scale=300,shadow_intensity = 1
, windowsize = c(800, 800)
, soliddepth =0
# ,           alpha              = 0.5
# ,           save_shadow_matrix = T
 )

# # capture the 3d plot
render_camera(zoom = 0.7, theta = 10, phi = 30)
Sys.sleep(0.2)
render_snapshot(clear = F)


### Section 04: density plot of all cell type --------------------------------------------------
library(viridis)
library(ggplot2)

densityplot <- celldata %>%
  filter(sample_id %in% c( "1310406")) %>%
  filter(phenotype != "Non-immune cell")

ggplot(densityplot, aes(x = x_x, y = y_y)) +
  # raster plot
  layer(
    stat     = "density_2d",
    geom     = "raster",
    mapping  = aes(fill = after_stat(count) * 200),
    params   = list(contour = FALSE, interpolate = T, h = 100, alpha = 0.9),
    position = "identity"
  ) +
  # group by phenotype
  facet_wrap(phenotype ~ .) +
    # scale color
    scale_fill_viridis(
        option   = "H",
        limits   = c(0, 1.2),
        na.value = "transparent"
      ,oob     = scales::squish
    ) +
    # scale_fill_gradient(
    # low    = "pink",
    # high   = "#ff0a7c",
    # # limits = c(0.2, 1.2)
    # #,      oob        = scales::squish
    # ,       na.value   = "transparent"
    # ) +
    theme_bw()+
    theme(
      plot.title = element_blank(),
      axis.title=element_blank(),
      panel.grid      = element_blank(),
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      panel.border    = element_blank(),
      # legend.position = "none",

    ) +
    coord_fixed()



densityplot <- celldata %>%
  filter(sample_id %in% c( "1311008")) %>%
  filter(phenotype != "Non-immune cell")

ggplot(densityplot, aes(x = x_x, y = y_y)) +
  # raster plot
  layer(
    stat     = "density_2d",
    geom     = "raster",
    mapping  = aes(fill = after_stat(count) * 200),
    params   = list(contour = FALSE, interpolate = T, h = 100, alpha = 0.9),
    position = "identity"
  ) +
  # group by phenotype
  facet_wrap(phenotype ~ .) +
    # scale color
    scale_fill_viridis(
        option   = "H",
        limits   = c(0, 1.2),
        na.value = "transparent"
      ,         oob         = scales::squish
    ) +
    # scale_fill_gradient(
    # low    = "pink",
    # high   = "#ff0a7c",
    # # limits = c(0.2, 1.2)
    # #,      oob        = scales::squish
    # ,       na.value   = "transparent"
    # ) +
    theme_bw()+
    theme(
      plot.title = element_blank(),
      axis.title=element_blank(),
      panel.grid      = element_blank(),
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      panel.border    = element_blank(),
      # legend.position = "none",

    ) +
    coord_fixed()
