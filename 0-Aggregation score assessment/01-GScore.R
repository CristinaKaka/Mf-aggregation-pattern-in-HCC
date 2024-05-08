### Section 01: preparation --------------------------------------------------

# set work directory
wdir <- "..."
setwd(wdir)

# library preparation
library(data.table)
library(dplyr)
library(spdep)
library(sp)

# data list (data for analysis)
filelist <- list.files(path = "data/", pattern = ".rdata")

### Section 02: calculate Gscore of each tile --------------------------------------------------
for (i in filelist) {
  # load data
  load(paste("data/", i, sep = "", collapse = ""))

  # function
  get_Gscore <- function(nb.size) {
    # Path to save
    grid_file <- paste("results/Gscore/n", as.character(nb.size), sep = "", collapse = "")

    # create directory
    if (file.exists(paste("results/Gscore/n", as.character(nb.size), sep = "", collapse = "")) == FALSE) {
      dir.create(paste("results/Gscore/n", as.character(nb.size), sep = "", collapse = ""))
    } else {
      print("Already existed")
    }

    # neighborhood identification
    data <- Pergrid
    colnames(data)[1] <- "TMAloc"
    min.d <- 0
    max.d <- (nb.size * sqrt(2))

    # Gscore
    GScore <- function(y) {
      data.copy <- y
      coordinates(y) <- c("x", "y")
      nlist <- dnearneigh(y, d1 = min.d, d2 = max.d)
      nbs <- nb2listw(include.self(nlist), style = "B", zero.policy = TRUE)
      GScore <- localG(y[[3]], nbs)
      class(GScore) <- "numeric"
      assign(paste("G.", names(y)[3], sep = "", collapse = ""), GScore)
      temp <- cbind(data.copy, get(paste("G.", names(y)[3], sep = "", collapse = "")))
      names(temp)[6] <- paste("G.", names(y)[3], sep = "", collapse = "")
      temp
    }
    temp1 <- split(data, data$TMAloc)
    data.wg <- do.call(rbind, lapply(temp1, function(z) GScore(z)))
    row.names(data.wg) <- (1:nrow(data.wg))
    data.wg
    # save(data.wg,file = paste(grid_file,"/",x))
  }

  # !! change neighborhood size
  nb.size <- 1

  data.wg <- get_Gscore(nb.size) # edit

  save(data.wg, file = paste("results/Gscore/n", nb.size, "/", i, sep = "", collapse = "")) # edit
}
