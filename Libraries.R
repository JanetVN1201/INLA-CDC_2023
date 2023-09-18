##Libraries necessary for the course
#Install INLA
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
library("INLA")

#Other libraries
libs = c("DiagrammeR","car","ggpubr","spdep","RColorBrewer",
         "spatstat","rgdal","sp","maptools","latticeExtra",
         "gridExtra","gstat","raster","ggplot2","ggfortify",
         "survival","joineR","BayesSurvival","icenReg","nloptr",
         "faraway","lme4","boot","sf","coda","spBayesSurv",
         "BayesX", "R2BayesX", "fields")
ix <- which(!sapply(libs, require, char = TRUE))
if (length(ix) > 0) {install.packages(pkgs[ix], repos = "https://cloud.r-project.org/")
  sapply(pkgs[ix], require, char = TRUE)}

