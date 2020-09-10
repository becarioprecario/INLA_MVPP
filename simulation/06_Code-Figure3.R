#########################################################################
#################### Figure 3 code of the manuscript ####################
#########################################################################

# In this script we are going to show the code for ploting the Figure 3 of
# the manuscript which belons to the simulation study

library(sp)
library(INLA)
library(maptools)
library(spatstat)
library(spatstat.utils)
library(rgeos)
library(RColorBrewer)
library(ggmap)
library(rgdal)
library(OpenStreetMap)
library(tmaptools)
library(mapproj)
library(grDevices)
library(lattice)
library(gridExtra)


####################################################################
###### Figure 3: Estimate of the spatial effect of the sim. study
####################################################################

# Model 2
mi = 2

#Creating the pdf archive for all the plots
pdf(file = "M2_SPDE1D_SimStu_phi3.pdf", width = 7, height = 7)
par(mfrow=c(3,2))

# Adjusting the model for the different number of cases
for(m in 1:6){
  
  #Different number of cases
  n.cases <- c(50, 100, 300, 500, 1000, 2000)
  
  #Loading results for models M2 and M3 (for phi = 3)
  cargar <- paste0("M2_M3_SPDE1D_cases", n.cases[m] , "phi3.RData" )
  load(cargar)
  
  l.model <- list(1, pp.res.M2.SPDE1D, pp.res.M3.SPDE1D)
  
  #Choosing between model 2 and model 3 
  pp.res <- l.model[[mi]]
  
  ######## Lung cancer 
  N <- length(pp.res$summary.random[[mi]][, 2])
  
  # Plotting the IC 95 with the effect, its median and its mode. 
  plot(mesh_cov$mid , pp.res$summary.random[[mi]][, 2], 
    main=paste0("SPDE 1D (Model ", mi, " with ", n.cases[m], " cases)"),
    ylab="Posterior estimate of the covariate effect", xlab="Distance (Km)",
    ylim =c( min(pp.res$summary.random[[mi]][, 4]),
      max(pp.res$summary.random[[mi]][, 6]) ) 
  )
  abline(h=0, col="grey")
  lines(mesh_cov$mid, 
    pp.res$summary.random[[mi]][, 4], col="red1")
  lines(mesh_cov$mid, 
    pp.res$summary.random[[mi]][, 6], col="red1")
  lines(mesh_cov$mid, 
    pp.res$summary.random[[mi]][, 5], col="deepskyblue")# median
  lines(mesh_cov$mid, 
    pp.res$summary.random[[mi]][, 7], col="darkmagenta")# mode
  legend("topright", legend = c(c("Effect", "95% CI"), c("Median", "Mode")),
    cex = 0.75, ncol = 2, 
    col=c("black", "red1",  "deepskyblue", "darkmagenta"), 
    pch =c(1, 1, NA, NA), lty = c(NA, NA, 1, 1))
  
  
  rm(list=setdiff(ls(), "mi"))
  
}

dev.off()



rm(list = ls())
####################################################################
