# set work directory
setwd("...")

# load data after GScore calculation
load("results/GScore/n1/CD68.rdata")

# hotspot cutoff

library(dplyr)

# function
hotspot <- function(x) {
  data <- x
  sq.count <- nrow(data)
  if (sq.count >= 1 & sq.count < 50) {
    g <- 1.645
  }
  if (sq.count >= 50 & sq.count < 100) {
    g <- 3.083
  }
  if (sq.count >= 100 & sq.count < 1000) {
    g <- 3.289
  }
  if (sq.count >= 1000) {
    g <- 3.886
  }
  temp <- data %>%
    mutate(hs.CD68 = ifelse(G.CD68 >= g, 1, 0))
  temp
}

temp1 <- split(data.wg, data.wg$TMAloc)
data <- do.call(rbind, lapply(temp1, function(z) hotspot(z)))
row.names(data) <- (1:nrow(data))

# export
if (file.exists("results/hotspot") == FALSE) {
  dir.create("results/hotspot")
} else {
  print("Already existed")
}

save(data, file = "results/hotspot/CD68_hs.rdata")
rm(list = ls())
