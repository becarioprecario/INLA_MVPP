############################################################################
################## Code for the tables of the manuscript. ##################
############################################################################

# In this script we are going to show the code for the tables 3 and 4 of 
# of the manuscript. 

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
#install.packages("remotes")
#remotes::install_github("tshmak/Tmisc")
library(Tmisc)

############################################################################
########## Table 3 
############################################################################

DicWaicMlik <- function(res){
  
  cargar <- paste0("Results/Res-", res, ".Rdata")
  load(cargar)
  
  v <- numeric(3)
  v[1] <- pp.res$dic$dic
  v[2] <- pp.res$waic$waic
  v[3] <- pp.res$mlik[2]
  return(v)
}


##################### Top Part

head1 <- c(" ", " ", "Model", "Choice", "Criteria", rep(" ", 3))  
head2 <- c(" ", "Model 0",  rep(" ", 3), "Model 1",  rep(" ", 2))
head3 <- c(rep(c( "DIC", "WAIC", "Marg. lik.", " "), 2) )   
col1 <- rep(" ", 5)
col2 <- c(rep(" ", 2), "Conf. factors", "No", "Yes")

#Model 0 and Model 1 WITHOUT confounding factors
row1 <- c(round(DicWaicMlik("M0-NoConFact"), 2), " ",
  round(DicWaicMlik("M1-NoConFact"), 2), " ")

#Model 0 and Model 1 WITH confounding factors
row2 <- c(round(DicWaicMlik("M0-ConFact"), 2), " ",
  round(DicWaicMlik("M1-ConFact"), 2), " " )

top <- cbind(col1, col2, rbind(head1, head2, head3, row1, row2))

rm(list=setdiff(ls(), "top"))

################## Middle Part

pickConfact <- function(res){
  cargar <- paste0("Results/Res-", res, ".Rdata")
  load(cargar)
  
  mat <- pp.res$summary.fixed[ -(1:4), c(1, 2, 3, 5)]
  return(mat)
}

col1 <- c(" ", "Variable", " ", rep(c("UNEMP2059", "AVGSOC", "AVGEDU", "PCTSCH"), 
                               rep(3, 4) ) )
col2 <- c(" ", "Cancer", " ", rep(c("lung", "stomach", "kidney"), 4))
head1 <- c(rep(" ", 1), "Effects", "of", "Confounding", "factors",
  rep(" ", 3))
head2 <- c( "Model 0", "with", "conf.", "factors", "Model 1", "with", 
  "conf.", "factors")
head3 <- c(rep(c("Mean", "St. dev.", "95% C.I.", " "), 2))

numb <- round(cbind(pickConfact("M0-ConFact"), pickConfact("M1-ConFact")), 2)

middle <- as.matrix(cbind(col1, col2, rbind(head1, head2, head3, numb)))
rm(list= ls()[!(ls() %in% c('top','middle'))])

################## Bottom Part

pickConfact <- function(res){
  cargar <- paste0("Results/Res-", res, ".Rdata")
  load(cargar)
  
  mat <- pp.res$summary.hyperpar[ , c(1, 2)]
  
  return(mat)
}

col1 <- c(" ", " ", "Pattern", " ", 
          rep(c("Controls", "Lung", "Stomach", "Kidney"), 2 ) )
col2 <- c(" ", " ", "Conf. factors", " ", rep(c("No", "Si"), c(4, 4)))
head1 <- c(rep(" ", 2), "Spatial", "Effects", rep(" ", 4))
head2 <- c( " ", "Model 0", " ", " ", " ", "Model 1", " ", " ")
head3 <- rep(c("Nominal", "range", "Nominal", "St. Dev."), 2)
head4 <- rep(c("Mean", "St. dev."), 4)

res1 <- as.vector(t(pickConfact("M0-NoConFact")))
row1 <- c(round(res1, 2), 
  round(as.vector(t(pickConfact("M1-NoConFact")[1:2, ])), 2))
row2 <- c(rep(" ", 4), 
  round(as.vector(t(pickConfact("M1-NoConFact")[3:4, ])), 2))
row3 <- c(rep(" ", 4), 
  round(as.vector(t(pickConfact("M1-NoConFact")[5:6, ])), 2))
row4 <- c(rep(" ", 4), 
  round(as.vector(t(pickConfact("M1-NoConFact")[7:8, ])), 2))
res2 <- as.vector(t(pickConfact("M0-ConFact")))
row5 <- c(round(res2, 2), 
  round(as.vector(t(pickConfact("M1-ConFact")[1:2, ])), 2))
row6 <- c(rep(" ", 4), 
  round(as.vector(t(pickConfact("M1-ConFact")[3:4, ])), 2))
row7 <- c(rep(" ", 4), 
  round(as.vector(t(pickConfact("M1-ConFact")[5:6, ])), 2))
row8 <- c(rep(" ", 4), 
   round(as.vector(t(pickConfact("M1-ConFact")[7:8, ])), 2))

mat <- rbind(row1, row2, row3, row4, row5, row6, row7, row8)

bottom <- cbind(col1, col2, rbind(head1, head2, head3, head4, mat))
rm(list= ls()[!(ls() %in% c('top','middle','bottom'))])
#Joining the three parts
  
Table3 <- as.data.frame(rbind(top, middle, bottom)); colnames(Table3)<- NULL
rownames(Table3)<- NULL; Table3 


write.table(Table3, file = "Table3.csv", row.names =FALSE, sep = ";", 
            col.names = FALSE)


############################################################################

rm(list=ls())
############################################################################
########### Table 4
############################################################################

dicwaic <- function(res){
  cargar <- paste0("Results/", res, ".Rdata")
  load(cargar)
  
  v <- numeric(2)
  v[1] <- pp.res$dic$dic
  v[2] <- pp.res$waic$waic
  v <- round(v, 2)
  return(v)
}

##################### Top Part (Without confounding factors)

col1 <- c(" ", " ", "Source", rep(c("Industry 1", "Industry 2", "Industry 5",
  "Industry 6", "Industry 7", "Industry 8", "Industry 9"), 3))
col2 <- c(" ", " ", "Effect", rep(c("Fixed", "RW1", "SPDE1D"), rep(7, 3)) )

head1 <- c(rep(" ", 1), "WITHOUT", "CONFOUNDING", "FACTORS", rep(" ", 2))
head2 <- c(" ", "Model 2", " ", " ", "Model 3", " ")
head3 <- rep(c("DIC", "WAIC", "Cancer"), 2)

mat<- rbind(
  c(dicwaic("M2-Var1_LinEffect_NoConFact"), "L,S", 
    dicwaic("M3-Var1_LinEffect_NoConFact"), "L,S"),
  c(dicwaic("M2-Var2_LinEffect_NoConFact"), "L,S", 
    dicwaic("M3-Var2_LinEffect_NoConFact"), "L,S"),
  c(dicwaic("M2-Var5_LinEffect_NoConFact"), "L,S", 
    dicwaic("M3-Var5_LinEffect_NoConFact"), "-"),
  c(dicwaic("M2-Var6_LinEffect_NoConFact"), "L,S, K", 
    dicwaic("M3-Var6_LinEffect_NoConFact"), "L"),
  c(dicwaic("M2-Var7_LinEffect_NoConFact"), "-", 
    dicwaic("M3-Var7_LinEffect_NoConFact"), "L"),
  c(dicwaic("M2-Var8_LinEffect_NoConFact"), "L,K", 
    dicwaic("M3-Var8_LinEffect_NoConFact"), "-"),
  c(dicwaic("M2-Var9_LinEffect_NoConFact"), "L,K", 
    dicwaic("M3-Var9_LinEffect_NoConFact"), "-"),
  c(dicwaic("M2_Var1_RW1_NoConFact"), "L", 
    dicwaic("M3_Var1_RW1_NoConFact"), "L,S"),
  c(dicwaic("M2_Var2_RW1_NoConFact"), "L", 
    dicwaic("M3_Var2_RW1_NoConFact"), "L,S"),
  c(dicwaic("M2_Var5_RW1_NoConFact"), "L,S", 
    dicwaic("M3_Var5_RW1_NoConFact"), "L,S"),
  c(dicwaic("M2_Var6_RW1_NoConFact"), "L", 
    dicwaic("M3_Var6_RW1_NoConFact"), "L,S"),
  c(dicwaic("M2_Var7_RW1_NoConFact"), "-", 
    dicwaic("M3_Var7_RW1_NoConFact"), "-"),
  c(dicwaic("M2_Var8_RW1_NoConFact"), "-", 
    dicwaic("M3_Var8_RW1_NoConFact"), "-"),
  c(dicwaic("M2_Var9_RW1_NoConFact"), "-", 
    dicwaic("M3_Var9_RW1_NoConFact"), "-"),
  c(dicwaic("M2_Var1_SPDE_NoConFact"), "L,S", 
    dicwaic("M3_Var1_SPDE_NoConFact"), "L,S"),
  c(dicwaic("M2_Var2_SPDE_NoConFact"), "L,S", 
    dicwaic("M3_Var2_SPDE_NoConFact"), "L,S"),
  c(dicwaic("M2_Var5_SPDE_NoConFact"), "L,S", 
    dicwaic("M3_Var5_SPDE_NoConFact"), "L,S"),
  c(dicwaic("M2_Var6_SPDE_NoConFact"), "L", 
    dicwaic("M3_Var6_SPDE_NoConFact"), "-"),
  c(dicwaic("M2_Var7_SPDE_NoConFact"), "-", 
    dicwaic("M3_Var7_SPDE_NoConFact"), "-"),
  c(dicwaic("M2_Var8_SPDE_NoConFact"), "-", 
    dicwaic("M3_Var8_SPDE_NoConFact"), "-"),
  c(dicwaic("M2_Var9_SPDE_NoConFact"), "-", 
    dicwaic("M3_Var9_SPDE_NoConFact"), "-")
)

top <- cbind(col1, col2, rbind(head1, head2, head3, mat))

rm(list= ls()[!(ls() %in% c('top','dicwaic'))])
##################### Bottom Part (With confounding factors)

col1 <- c(" ", " ", "Source", rep(c("Industry 1", "Industry 2", "Industry 5",
  "Industry 6", "Industry 7", "Industry 8", "Industry 9"), 3))
col2 <- c(" ", " ", "Effect", rep(c("Fixed", "RW1", "SPDE1D"), rep(7, 3)) )

head1<- c(rep(" ", 1), "WITH", "CONFOUNDING", "FACTORS", rep(" ", 2))
head2 <- c(" ", "Model 2", " ", " ", " ", " ")
head3 <- c("DIC", "WAIC", "Cancer", rep(" ", 3)) 

mat<- rbind(
  c(dicwaic("M2-Var1_LinEffect"), "L,S,K", rep(" ", 3)),
  c(dicwaic("M2-Var2_LinEffect"), "L,S,K", rep(" ", 3)),
  c(dicwaic("M2-Var5_LinEffect"), "L,S,K", rep(" ", 3)),
  c(dicwaic("M2-Var6_LinEffect"), "L,S,K", rep(" ", 3)),
  c(dicwaic("M2-Var7_LinEffect"), "K", rep(" ", 3)),
  c(dicwaic("M2-Var8_LinEffect"), "L,K", rep(" ", 3)),
  c(dicwaic("M2-Var9_LinEffect"), "L,K", rep(" ", 3)),
  c(dicwaic("M2_Var1_RW1"), "L,S", rep(" ", 3)),
  c(dicwaic("M2_Var2_RW1"), "L,K", rep(" ", 3)),
  c(dicwaic("M2_Var5_RW1"), "L,S", rep(" ", 3)),
  c(dicwaic("M2_Var6_RW1"), "L,S", rep(" ", 3)),
  c(dicwaic("M2_Var7_RW1"), "-", rep(" ", 3)),
  c(dicwaic("M2_Var8_RW1"), "L", rep(" ", 3)),
  c(dicwaic("M2_Var9_RW1"), "L", rep(" ", 3)),
  c(dicwaic("M2_Var1_SPDE"), "L,S", rep(" ", 3)),
  c(dicwaic("M2_Var2_SPDE"), "L,S", rep(" ", 3)),
  c(dicwaic("M2_Var5_SPDE"), "L,S", rep(" ", 3)),
  c(dicwaic("M2_Var6_SPDE"), "-", rep(" ", 3)),
  c(dicwaic("M2_Var7_SPDE"), "-", rep(" ", 3)),
  c(dicwaic("M2_Var8_SPDE"), "-", rep(" ", 3)),
  c(dicwaic("M2_Var9_SPDE"), "-", rep(" ", 3))
)

bottom <- cbind(col1, col2, rbind(head1, head2, head3, mat))
rm(list= ls()[!(ls() %in% c('top','bottom'))])
#Joining the three parts

Table4 <- as.data.frame(rbind(top, bottom)); colnames(Table4)<- NULL
rownames(Table4)<- NULL; Table4 


write.table(Table4, file = "Table4.csv", row.names =FALSE, sep = ";", 
            col.names = FALSE)

############################################################################
rm(list=ls())