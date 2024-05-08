# Section 1: distance calculation -----------------------------------------

# loading packages
library(tidyverse)
library(phenoptr)
library(rtree)
# library(devtools)
# devtools::install_github('akoyabio/rtree')
# devtools::install_github("akoyabio/phenoptr")


# read single cell data including 'CellID','Cell X Position','Cell Y Position','Phenotype'
csd <- readRDS('')

# quphath ID (sample ID)
qupathID <- unique(csd$QupathID)

csdD_all <- data.frame()

for (i in qupathID) {
  
  # subset sample
  csd1 <- csd[csd$QupathID==i,]
  
  # find nearest distance to each phenotypes
  distances <- find_nearest_distance(csd1)
  
  # create a combined data frame
  csd_with_distance <- bind_cols(csd1, distances)
  
  # count cell number within given value
  count_within1 <- function(csd, from, to, radius, category = NA, dst = NULL) {
    stopifnot(length(radius) > 0, all(radius > 0))
    if (!is.na(category)) {
      category_cells <- csd$`Tissue Category` == category
      csd <- csd[category_cells, ]
      if (!is.null(dst)) {
        dst <- dst[category_cells, category_cells, drop = FALSE]
      }
    }
    if (is.null(dst)) {
      dst <- distance_matrix(csd)
    }
    dst <- subset_distance_matrix(csd, dst, from, to)
    if (prod(dim(dst)) > 0) {
      purrr::map_df(radius, function(rad) {
        within <- apply(dst, 1, function(r) {
          sum(r > 0 & r <=
            rad)
        })
      })
    } else {
      NA
    }
  }
  
    
  within <- count_within1(csd1, from=unique(csd1$Phenotype), to='Mf', radius=50)
  within <- t(within)
  colnames(within) <- '50um_Mf'
  
  csd_with_distance <- bind_cols(csd_with_distance, within)
}


# renew the nearest distance to CD8 CTL
min_col_idx <- apply(csdD_all[, c("Distance to CD8Tc", "Distance to PD1_CD8Tc", "Distance to LAG3_PD1_CD8Tc")], 1, function(x) {
  min_col_name <- names(x)[which.min(x)]
  which(names(csdD_all) == min_col_name)
})

csdD_all$`Cell ID CD8Tc` <- sapply(1:nrow(csdD_all), function(i) {
  min_col <- min_col_idx[[i]]
  csdD_all[i, min_col + 1]
})

csdD_all$`Distance to CD8Tc` <- sapply(1:nrow(csdD_all), function(i) {
  min_col <- min_col_idx[[i]]
  csdD_all[i, min_col]
})



# Section 2: distance normalization ---------------------------------------------

# subset macrophages 
csd_Mf <- csdD_all %>% 
  filter(Phenotype=='Mf') 


csd_Mf  <- csd_Mf  %>%
  group_by(QupathID) %>%
  mutate(normalized_Mf_Dis = (`Distance to Mf`/max(`Distance to Mf`))*100)

csd_Mf  <- csd_Mf  %>%
  group_by(QupathID) %>%
  mutate(normalized_CD8Tc = (`Distance to CD8Tc`/max(`Distance to CD8Tc`))*100)

csd_Mf  <- csd_Mf  %>%
  group_by(QupathID) %>%
  mutate(normalized_LAG3PD1Tc = (`Distance to LAG3_PD1_CD8Tc`/max(`Distance to LAG3_PD1_CD8Tc`))*100)

csdD_Mf <- csdD_all %>%
  `50um_LAG3_PD1_CD8Tc_R` = ((`50um_LAG3_PD1_CD8Tc`)/(`50um_CD8Tc`))*100
  )
