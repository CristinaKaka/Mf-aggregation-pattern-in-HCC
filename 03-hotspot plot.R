# packages loading
library(stringr)
library(data.table)
library(dplyr)

# data loading
load("results/hotspot/CD68_hs.rdata")
unique(data$TMAloc)


# Replace with actual core IDs
coreids <- c(
  "TMA131-f9", "TMA135-g1",
  "TMA135-g3", "TMA135-g4",
  "TMA147-d2", "TMA147-e2",
  "TMA147-h1", "TMA119-j9", "TMA139-a4", "TMA143-b4"
)
# Replace with actual cutoff values
hotspot_cuts <- c(3.886, 2.5) 

for(coreid in coreids) {
  for(hotspot_cut in hotspot_cuts) {
    coreplot <- paste(coreid, ".rdata", sep = "")
    filename <- paste("figures/", coreid, "-", hotspot_cut, ".png", sep = "")
    # Load data if it's stored in an .rdata file
    # load(coreplot)

    data.naomit <- subset(data, data$TMAloc == coreplot)
    l.zm.l <- matrix(nrow = 55, ncol = 55)

    for (i in 1:nrow(data.naomit)) {
      l.zm.l[data.naomit$x[i], data.naomit$y[i]] <- data.naomit$G.CD68[i]
    }
    
    png(file = filename, width = 5000, height = 5000)
    x.i <- 1 * (1:(nrow(l.zm.l)))
    y.i <- 1 * (1:(ncol(l.zm.l)))
    par(mar = c(0, 0, 0, 0))
    
    image(x.i, y.i, l.zm.l, xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
          col = c("mistyrose", "white", "mistyrose", "navy"), 
          breaks = c(min(data.naomit$G.CD68) - 10, -0.0000000000000000001, 
                     0.0000000000000000001, hotspot_cut, 
                     max(data.naomit$G.CD68) + 10))
    
    dev.off()
  }
}
