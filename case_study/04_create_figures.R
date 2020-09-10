############################################################################
################### Code for the plots of the manuscript ###################
############################################################################

# In this script we are going to show the code for the different plots of 
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



####################################################################
###### Figure 1: Plotting the data
####################################################################

#Loading the mesh
load("Data/Stack-SPDE1-Models.RData")
rm(list=setdiff(ls(), "mesh"))

#Load data with variability
load("Data/Simulated-CaseControldata-github.RData")
rm(air.ind, metals.ind)

#Changing to SpatialPoints object
cont.vir <- as.SpatialPoints.ppp(cls.sim)
pulmon <- as.SpatialPoints.ppp(lung.sim)
estom<- as.SpatialPoints.ppp(stom.sim)
rinon <- as.SpatialPoints.ppp(kid.sim)

#Changing from km to meters (note that this breaks the CRS information)
cont.vir@coords <- cont.vir@coords * 1000
pulmon@coords <- pulmon@coords * 1000
estom@coords <- estom@coords * 1000
rinon@coords <- rinon@coords * 1000


### Covariates data
metals.ind <- 
  read.table(file = "Data/focos_metales_pesados.txt", header = TRUE,
    sep = "", dec = ".", quote = "\"", row.names = NULL)

air.ind <-
  read.table(file = "Data/focos_aire_contam.txt", header = TRUE,
    sep = "", dec = ".", quote = "\"", row.names = NULL, encoding = "UTF8")

coordinates(air.ind) <- ~ xed50utm + yed50utm
coordinates(metals.ind) <- ~ xed50utm + yed50utm


proj4string(air.ind) <-
  CRS("+proj=utm +zone=30 +init=epsg:23030 +towgs84=-87,-96,-120,0,0,0,0")
air.ind.proj <- spTransform(air.ind, CRS("+proj=longlat +datum=WGS84"))
air <- as.data.frame(coordinates(air.ind.proj))
names(air) <- c("lon", "lat")

proj4string(metals.ind) <-
  CRS("+proj=utm +zone=30 +init=epsg:23030 +towgs84=-87,-96,-120,0,0,0,0")
metals.ind.proj <- spTransform(metals.ind, CRS("+proj=longlat +datum=WGS84"))
metals <- as.data.frame(coordinates(metals.ind.proj))
names(metals) <- c("lon", "lat")


## Getting the border of the mesh
owin2Polygons <- function(x, id="1") {
  stopifnot(is.owin(x))
  x <- as.polygonal(x)
  closering <- function(df) { df[c(seq(nrow(df)), 1), ] }
  pieces <- lapply(x$bdry,
                   function(p) {
                     Polygon(coords=closering(cbind(p$x,p$y)),
                             hole=is.hole.xypolygon(p))  })
  z <- Polygons(pieces, id)
  return(z)
}

owin2SP <- function(x) {
  stopifnot(is.owin(x))
  y <- owin2Polygons(x)
  z <- SpatialPolygons(list(y))
  return(z)
}

# Changing now the coords of the border of the mesh to meters
overlay <- convexhull.xy(mesh$loc[,1]*1000, mesh$loc[,2]*1000);
overlay <- owin2SP(overlay)
proj4string(overlay) <- 
  CRS("+proj=utm +zone=30 +init=epsg:23030 +towgs84=-87,-96,-120,0,0,0,0")
overlay <- spTransform(overlay, CRS("+proj=longlat +datum=WGS84"))
overlay <- fortify(overlay, region="COMMUNITY")

###### Controls 
#Assign CRS
proj4string(cont.vir) <- 
  CRS("+proj=utm +zone=30 +init=epsg:23030 +towgs84=-87,-96,-120,0,0,0,0")

#Changing CRS of the data
cont.vir.proj <- spTransform(cont.vir, CRS("+proj=longlat +datum=WGS84"))
pt <- as.data.frame(coordinates(cont.vir.proj))
names(pt) <- c("lon", "lat")

#Finally, plotting controls
png(filename = "Results/Controls_data.png", width = 720, height = 720)
ggmap(get_stamenmap(
  rbind(as.numeric(paste(geocode_OSM("Alcala de Henares")$bbox))), 
    zoom = 13), legend="bottomright") +
  geom_point(data=pt, colour="black", size=0.5) + 
  geom_point(data=air, colour="red",  pch=17, size=5) +
  geom_point(data=metals, colour="chartreuse", pch=16, size=3) + 
  geom_polygon(data=overlay, aes(x=long, y=lat, group=group), 
    color="black", size= 1.5, alpha=0) + coord_map() + theme_bw()
dev.off()


#####  Lung cancer

#Assign CRS
proj4string(pulmon) <- 
  CRS("+proj=utm +zone=30 +init=epsg:23030 +towgs84=-87,-96,-120,0,0,0,0")

#Changing CRS of the data
lung.proj <- spTransform(pulmon, CRS("+proj=longlat +datum=WGS84"))
pt <- as.data.frame(coordinates(lung.proj))
names(pt) <- c("lon", "lat")

png(filename = "Results/Lung_data.png", width = 720, height = 720)
ggmap(get_stamenmap(rbind(
    as.numeric(paste(geocode_OSM("Alcala de Henares")$bbox))), zoom = 13)) + 
  geom_point(data=pt, colour="blue", size=1)+ 
  geom_point(data=air, colour="red",  pch=17, size=5 ) +
  geom_point(data=metals, colour="chartreuse", pch=16, size=3) + 
  geom_polygon(data=overlay, aes(x=long, y=lat, group=group), 
    color="black", size= 1.5, alpha=0) + coord_map() + theme_bw()
dev.off()

###### Stomach cancer
#Assign CRS
proj4string(estom) <- 
  CRS("+proj=utm +zone=30 +init=epsg:23030 +towgs84=-87,-96,-120,0,0,0,0")

#Changing CRS of the data
stomach.proj <- spTransform(estom, CRS("+proj=longlat +datum=WGS84"))
pt <- as.data.frame(coordinates(stomach.proj))
names(pt) <- c("lon", "lat")

png(filename = "Results/Stomach_data.png", width = 720, height = 720)
ggmap(get_stamenmap(rbind(
    as.numeric(paste(geocode_OSM("Alcala de Henares")$bbox))), zoom = 13)) + 
  geom_point(data=pt, colour="purple", size=1) + 
  geom_point(data=air, colour="red",  pch=17, size=5 ) +
  geom_point(data=metals, colour="chartreuse", pch=16, size=3)  + 
  geom_polygon(data=overlay, aes(x=long, y=lat, group=group), 
    color="black", size= 1.5, alpha=0) + coord_map() + theme_bw()
dev.off()


###### Kidney cancer
#Assign CRS
proj4string(rinon) <- 
  CRS("+proj=utm +zone=30 +init=epsg:23030 +towgs84=-87,-96,-120,0,0,0,0")

#Changin CRS of the data
kidney.proj <- spTransform(rinon, CRS("+proj=longlat +datum=WGS84"))
pt <- as.data.frame(coordinates(kidney.proj))
names(pt) <- c("lon", "lat")

png(filename = "Results/Kidney_data.png", width = 720, height = 720)
ggmap(get_stamenmap(rbind(
    as.numeric(paste(geocode_OSM("Alcala de Henares")$bbox))), zoom = 13)) + 
  geom_point(data=pt, colour="brown", size=1) + 
  geom_point(data=air, colour="red",  pch=17, size=5 ) +
  geom_point(data=metals, colour="chartreuse", pch=16, size=3)  + 
  geom_polygon(data=overlay, aes(x=long, y=lat, group=group), 
    color="black", size= 1.5, alpha=0) + coord_map() + theme_bw()
dev.off()

rm(list = ls())
####################################################################

####################################################################
###### Figure 2: Estimating intensity and RR with spatstat
####################################################################

#Loading the mesh
load("Data/Stack-SPDE1-Models.RData")
rm(list=setdiff(ls(), "mesh"))

#Load simulated data
load("Data/Simulated-CaseControldata-github.RData")

#Aproppriate names
controls_ppp <- cls.sim
lung_ppp<- lung.sim
sto_ppp <- stom.sim
kid_ppp <- kid.sim
pts_ppp <- superimpose(controls_ppp,lung_ppp, sto_ppp, kid_ppp)

#Define the border
barrier <- convexhull.xy(mesh$loc[,1], mesh$loc[,2]) 
pts_ppp$window <- barrier
controls_ppp$window <- barrier
lung_ppp$window <- barrier
sto_ppp$window <- barrier
kid_ppp$window <- barrier


# Using density to estimate the density
dens <- density(pts_ppp, dimyx = c(300, 300))
pts.dens <- 
  as(density(controls_ppp, xy = dens), "SpatialGridDataFrame")
pts.dens <- 
  cbind(pts.dens,as(density(lung_ppp, xy = dens), "SpatialGridDataFrame"))
pts.dens <- 
  cbind(pts.dens,as(density(sto_ppp, xy = dens), "SpatialGridDataFrame"))
pts.dens <- 
  cbind(pts.dens,as(density(kid_ppp, xy = dens), "SpatialGridDataFrame"))

names(pts.dens) <- 
  c("dens.controls", "dens.lung", "dens.stomach", "dens.kidney")

pts.dens <- as(pts.dens, "SpatialPixelsDataFrame")


## We are going to plot the intensity of the cases and the controls 

#intensity.ppp compute the estimator of the intensity
# assuming stationarity 
controls.int <- intensity.ppp(controls_ppp); controls.int


#Now we calculate and plot the density of all cancers 
s.bw <- bw.diggle(controls_ppp); s5 <- 0.3; 
png(filename = "Results/Intensity_Controlss03.png", 
    width = 480, height = 480)
plot(density(controls_ppp, sigma=s5), col=gray.colors(30), main="Controls")
dev.off()

png(filename = "Results/Intensity_Lungs03.png", width = 480, height = 480)
plot(density(lung_ppp, sigma=s5), col=gray.colors(30), main="Lung cancer")
dev.off()

png(filename = "Results/Intensity_Stomachs03.png", width = 480, height = 480)
plot(density(sto_ppp, sigma=s5), col=gray.colors(30), main="Stomach cancer")
dev.off()

png(filename = "Results/Intensity_Kidneys03.png", width = 480, height = 480)
plot(density(kid_ppp, sigma=s5), col=gray.colors(30), main="Kidney cancer")
dev.off()

#Changing to SpatialPolygons
bdy.SP <- as(owin(poly = bdy), "SpatialPolygons")

m <- 
  c(rep(0, controls_ppp$n), rep(1, lung_ppp$n),
    rep(2, sto_ppp$n), rep(3, kid_ppp$n))
m <- factor(m, levels=0:3)
levels(m) <- 
  c("Controls", "Lung cancer", "Stomach cancer", "Kidney cancer")
mpp <- ppp(pts_ppp$x,pts_ppp$y, window =pts_ppp$window , marks=m)
mpp$window <- barrier

## We are going to estimate the relative risk 

rr <- relrisk(mpp, relative=TRUE, control="Controls", sigma = s5)

png(filename = "Results/RR_Lung_s03.png", width = 480, height = 480)
plot(rr$`Lung cancer`, col=gray.colors(30), main="Lung Cancer")
dev.off()

png(filename = "Results/RR_Stomach_s03.png", width = 480, height = 480)
plot(rr$`Stomach cancer`, col=gray.colors(30), main="Stomach Cancer")
dev.off()

png(filename = "Results/RR_Kidney_s03.png", width = 480, height = 480)
plot(rr$`Kidney cancer`, col=gray.colors(30), main="Kidney Cancer")
dev.off()

rm(list = ls())
####################################################################

####################################################################
###### Figure 4: Spatial specific effect estimates
####################################################################

#Load results model 1
load("Results/Res-M1-ConFact.RData")


################ Controls 
hp.res <-inla.spde2.result(inla=pp.res, name="spatial.field1", 
  spde=spde, do.transform = T)

# Projecting the mesh and the surfaces into the study region
gproj <- 
  inla.mesh.projector(mesh, xlim=c(464, 472.5), ylim =c(4478.5, 4487.5))
g.mean <- inla.mesh.project(gproj, hp.res$summary.values$mean)
g.sd <- inla.mesh.project(gproj, hp.res$summary.values$sd)

# Palette
gre <- brewer.pal(9, "Greys"); 
gre <- colorRampPalette(gre)(18); gre <- gre[-(1:2)]

#Plotting posterior spatial effects 
png(filename = "Results/Si_controls_M1.png", width = 720, height = 360)
grid.arrange(
  levelplot(g.mean, scales=list(draw=F), xlab='', ylab='',  
    main=expression("Posterior mean of S "[i]), col.regions=gre),
  levelplot(g.sd, scal=list(draw=F), xlab='', ylab='', 
    main=expression("Posterior sd of S "[i]), col.regions=gre), 
  nrow=1)
dev.off()

# Plotting the ordered estimation of the posterior mean of the 
# sptaial effect 
lowIC <- hp.res$summary.values$`0.025quant`
m     <- hp.res$summary.values$mean
UpIC  <- hp.res$summary.values$`0.975quant`


mt <- cbind(lowIC, m, UpIC); 
mt.or <- mt[order(mt[,2]),]

png(filename = "Results/Si_ordered_controls_M1.png", 
    width = 720, height = 480)
plot(mt.or[, 2], type="l", ylim=c(-25.0, 25.0), main = "Controls",
  ylab=expression("Ordered estimations of S "[i]))
abline(h=0, col="grey")
points(mt.or[, 1], col="darkgrey", type="l")
points(mt.or[, 3], col="darkgrey", type="l")
legend("topleft", legend = c(expression("S "[i]), "95% Credible Interval"),
  col=c("black", "darkgrey"), lty = c(1, 1))
dev.off()

rm(hp.res, g.mean, g.sd, gproj, lowIC, UpIC, m, mt, mt.or)



############### Lung cancer
hp.res <-inla.spde2.result(inla=pp.res, name="spatial.field2", spde=spde,
  do.transform = T)

# Projecting the mesh and the surfaces into the study region
gproj <- 
  inla.mesh.projector(mesh, xlim=c(464, 472.5), ylim =c(4478.5, 4487.5))
g.mean <- inla.mesh.project(gproj, hp.res$summary.values$mean)
g.sd <- inla.mesh.project(gproj, hp.res$summary.values$sd)


blu <- brewer.pal(9, "Blues"); 
blu <- colorRampPalette(blu)(18); blu <- blu[-(1:2)]

png(filename = "Results/Si_lung_M1.png", width = 720, height = 360)
grid.arrange(
  levelplot(g.mean, scales=list(draw=F), xlab='', ylab='', 
    main=expression("Posterior mean of S "[i]), col.regions=blu),
  levelplot(g.sd, scal=list(draw=F), xlab='', ylab='', 
    main=expression("Posterior sd of S "[i]), col.regions=blu), 
  nrow=1)
dev.off()


# Plotting the ordered estimation of the posterior mean of the 
# sptaial effect 
lowIC <- hp.res$summary.values$`0.025quant`
m     <- hp.res$summary.values$mean
UpIC  <- hp.res$summary.values$`0.975quant`

mt <- cbind(lowIC, m, UpIC); 
mt.or <- mt[order(mt[,2]),]

png(filename = "Results/Si_ordered_lung_M1.png", width = 720, height = 480)
plot(mt.or[, 2], type="l", ylim=c(-0.55, 0.55), col="dodgerblue4", 
  main = "Lung cancer", ylab=expression("O_rdered estimations of S "[i]))
abline(h=0,  col="grey") 
points(mt.or[, 1], col="blue", type="l") #col="gray64"
points(mt.or[, 3], col="blue", type="l") #col="gray64"
legend("topleft", legend = c(expression("S "[i]), "95% Credible Interval"),
  col=c("dodgerblue4", "blue"), lty = c(1, 1))
dev.off()

rm(hp.res, g.mean, g.sd, gproj, lowIC, UpIC, m, mt, mt.or)



############### Stomach cancer
hp.res <-inla.spde2.result(inla=pp.res, name="spatial.field3", spde=spde,
  do.transform = T)

# Projecting the mesh and the surfaces into the study region
gproj <- 
  inla.mesh.projector(mesh, xlim=c(464, 472.5), ylim =c(4478.5, 4487.5))
g.mean <- 
  inla.mesh.project(gproj, hp.res$summary.values$mean)
g.sd <- 
  inla.mesh.project(gproj, hp.res$summary.values$sd)

#Palette
pur <- brewer.pal(9, "Purples");
pur <- colorRampPalette(pur)(18); pur <- pur[-(1:2)] 

png(filename = "Results/Si_stomach_M1.png", width = 720, height = 360)
grid.arrange(
  levelplot(g.mean, scales=list(draw=F), xlab='', ylab='', 
    main=expression("Posterior mean of S "[i]), col.regions=pur),
  levelplot(g.sd, scal=list(draw=F), xlab='', ylab='', 
    main=expression("Posterior sd of S "[i]), col.regions=pur), 
  nrow=1)
dev.off()

# Plotting the ordered estimation of the posterior mean of the 
# sptaial effect 
lowIC <- hp.res$summary.values$`0.025quant`
m     <- hp.res$summary.values$mean
UpIC  <- hp.res$summary.values$`0.975quant`

mt <- cbind(lowIC, m, UpIC); 
mt.or <- mt[order(mt[,2]),]

png(filename = "Results/Si_ordered_stomach_M1.png", 
  width = 720, height = 480)
plot(mt.or[, 2], type="l", ylim=c(-0.55, 0.55), col="darkorchid1",
  main = "Stomach cancer", ylab=expression("Ordered estimations of S "[i]))
abline(h=0, col="grey")
points(mt.or[, 1], col="purple", type="l") #col="gray64"
points(mt.or[, 3], col="purple", type="l") #col="gray64"
legend("topleft", legend = c(expression("S "[i]), "95% Credible Interval"),
  col=c("darkorchid1", "purple"), lty = c(1, 1)) #col=c("black", "gray64") 
dev.off()


rm(hp.res, g.mean, g.sd, gproj, lowIC, UpIC, m, mt, mt.or)



############### Kidney cancer
hp.res <-inla.spde2.result(inla=pp.res, name="spatial.field4", spde=spde,
  do.transform = T)

# Projecting the mesh and the surfaces into the study region
gproj <- 
  inla.mesh.projector(mesh, xlim=c(464, 472.5), ylim =c(4478.5, 4487.5))
g.mean <- 
  inla.mesh.project(gproj, hp.res$summary.values$mean)
g.sd <- 
  inla.mesh.project(gproj, hp.res$summary.values$sd)

#Palette
bro <- brewer.pal(9, "YlOrBr"); 
bro <- colorRampPalette(bro)(18); bro <- bro[-(1:2)]

png(filename = "Results/Si_kidney_M1.png", width = 720, height = 360)
grid.arrange(
  levelplot(g.mean, scales=list(draw=F), xlab='', ylab='', 
    main=expression("Posterior mean of S "[i]), col.regions=bro),
  levelplot(g.sd, scal=list(draw=F), xlab='', ylab='', 
    main=expression("Posterior sd of S "[i]), col.regions=bro), 
  nrow=1)
dev.off()

# Plotting the ordered estimation of the posterior mean of the 
# sptaial effect 
lowIC <- hp.res$summary.values$`0.025quant`
m     <- hp.res$summary.values$mean
UpIC  <- hp.res$summary.values$`0.975quant`
mt <- cbind(lowIC, m, UpIC); 
mt.or <- mt[order(mt[,2]),]

png(filename = "Results/Si_ordered_kidney_M1.png", 
  width = 720, height = 480)
plot(mt.or[, 2], type="l", ylim=c(-0.55, 0.55), col="chocolate", 
  main = "Kidney cancer", ylab=expression("Ordered estimations of S "[i]))
abline(h=0, col="grey")
points(mt.or[, 1], col="darkorange1", type="l") #col="gray64"
points(mt.or[, 3], col="darkorange1", type="l") #col="gray64"
legend("topleft", legend = c(expression("S "[i]), "95% Credible Interval"),
  col=c("chocolate", "darkorange1"), lty = c(1, 1)) 
dev.off()

rm(hp.res, g.mean, g.sd, gproj, lowIC, UpIC, m, mt, mt.or)

rm(list = ls())

####################################################################

####################################################################
###### Figure 5: Estimated exposure effect of the Industry 1
####################################################################

## Note: This code has been design in order to plot all the 
## combinations (not only the effect of Industry I with Model 2)

#Labels for cancer
lab.cancer <- c("lung", "stomach", "kidney")

#Plot results from RW1 model
#obj: INLA model
#cancer: Type of cancer (1:3)
#m: model
#p: Industry
# mesh_cov: mesh_cov when defining 1D SPDE
plot0 <- function(obj, cancer, m = "", p = "", mesh_cov = NULL) {
  
  if(substr(m, 1, 1) == 2) {
    tab <- obj$summary.random[[1 + cancer]]
  } else { #M3
    tab <- obj$summary.random[[4 + cancer]]
  }
  # Plotting the IC 95 with the effect, its median and its mode. 
  if(length(grep("RW1", m)) == 1) {#RW1
    xx <- tab[, 1]
  }  else {
    xx <- mesh_cov$mid
  }
  plot(xx, tab[, 2], 
    main = paste0(lab.cancer[cancer], " cancer (M", m, ", Ind. ", p, ")"),
    type = "l", #xlim = c(0, 6),
    ylab = "Effect estimate", xlab = "Distance (km)",
    ylim = c(min(tab[, 4]), max(tab[, 6]) )
  )
  abline(h = 0, col = "grey")
  lines(xx, tab[, 4], col = "red1")
  lines(xx, tab[, 6], col = "red1")
  legend("topright", legend = c("Mean", "95% CI"), cex = 0.75, ncol = 2,
    col = c("black", "red1"), lty = 1, bty = "n")
  
}

# Create several plots
# filename: Function to return the name of the RW1 and SPDE1 models given p
summary0 <- function(filename, m) {
  
  effect <- c("RW1", "SPDE1")
  
  tbl <- sapply(1:13, function(p, filename) {
    
    print(paste0("INDUSTRY ", p))
    
    par(mfcol = c(3, 2))
    for(s in c(1:2)[file.exists(filename(p))]) {
      #Load results of the model for the variable p
      load(filename(p)[s])
      for(cancer in 1:3) {
        if(s == 2) { #SPDE
          plot0(pp.res, cancer, paste0(m, "-", effect[s]), p, mesh_cov = mesh_cov) 
        } else {
          plot0(pp.res, cancer, paste0(m, "-", effect[s]), p) 
        }
      }
    }
    
  }, filename = filename)
}

filename <- function(p) {
  if(p == 10) { p <- 9} #We do not have data for source 10
  c(
    paste0("Results/M2_Var", p, "_RW1.RData"),
    paste0("Results/M2_Var", p, "_SPDE.RData")
  )
}

summary0(filename, 2)

####################################################################


