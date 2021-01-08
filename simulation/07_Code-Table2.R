############################################################################
#################### Code for Table 2 of the manuscript ####################
############################################################################

# In this script we are going to show the code for reproducing Table 2 of 
# the manuscript. 

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

############################################################################
# First we will define a function to select the DIC value


pickDIC <- function(res){
  v <- res$dic$dic
  return(v)
}

load.and.pickM0M1 <- function(phi){
  
  resul <-paste0("M0_M1_cases", j, "phi", phi,".Rdata" )
  load(resul)

  x <- pickDIC(pp.res.M0)
  y <- pickDIC(pp.res.M1)
  v <- c(x, y); return(v)
}

load.and.pickM2M3Rw1 <- function(phi){
  
  resul <-paste0("M2_M3_Rw1_cases", j, "phi", phi,".Rdata" )
  load(resul)
  
  x <- pickDIC(pp.res.M2.Rw1)
  y <- pickDIC(pp.res.M3.Rw1)
  v <- c(x, y); return(v)
}

load.and.pickM2M3SPDE1D <- function(phi){
  
  resul <-paste0("M2_M3_SPDE1D_cases", j, "phi", phi,".Rdata" )
  load(resul)
  
  x <- pickDIC(pp.res.M2.SPDE1D)
  y <- pickDIC(pp.res.M3.SPDE1D)
  v <- c(x, y); return(v)
}

col1 <- matrix(0, nrow = 2, ncol = 2)
col2 <- matrix(0, nrow = 2, ncol = 2)
col3 <- matrix(0, nrow = 2, ncol = 2)


for(j in c(50,100,300, 500, 1000, 2000)){
  
  v1 <- t(sapply(c(1, 3, 6), load.and.pickM0M1))
  col1 <- rbind(col1, v1)
  
  v2 <- t(sapply(c(1, 3, 6), load.and.pickM2M3Rw1))
  col2 <- rbind(col2, v2)
  
  v3 <- t(sapply(c(1, 3, 6), load.and.pickM2M3SPDE1D))
  col3 <- rbind(col3, v3)

}

#Leaving 2 decimals
col1 <- round(col1, 2)
col2 <- round(col2, 2)
col3 <- round(col3, 2)

#Headings of the table
col1[1, ] <- c(" ", " ")
col1[2, ] <- c("Model 0", "Model 1")

col2[1, ] <- c("RW1", " ")
col2[2, ] <- c("Model 2", "Model 3")

col3[1, ] <- c("SPDE1D", " ")
col3[2, ] <- c("Model 2", "Model 3")

col0 <- cbind(rep(c(50, 100, 300, 500, 1000, 2000), rep(3, 6)), 
  rep(c(1, 3, 6), 6))
col0 <- rbind(c("Settings", " "), c("#Cases", "Phi"), col0)

Table2 <- as.data.frame(cbind(col0, col1, col2, col3)); Table2

#Table2[3:20, ] <- as.numeric(Table2[3:20, ])

write.table(Table2, file = "Table2.csv", row.names =FALSE, sep = ";", 
          col.names = FALSE)

rm(list=ls())