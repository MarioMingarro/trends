packages.to.use <- c("readxl",
                     "tidyverse",
                     "sp",
                     "raster",
                     "strucchange",
                     "tictoc",
                     "doParallel",
                     "foreach")

packages.to.use <- unique(packages.to.use)

for(package in packages.to.use) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages(package ) }
  if( ! package %in% rownames(installed.packages()) ) { stop("Error on package instalation") }
  suppressWarnings( library(package, character.only = TRUE) )
}
