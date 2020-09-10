##########################################################################
## Creating the data stack for all the models with confounding covariates ##
##########################################################################

# In this script we are going to build the data stack to run all the models  
# since model 0 without covariates and specific effects to model 3 with 
# covariates and fixed effects. We will consider penalized complexity priors 
# and adding 4 different confounding factors.

#Loading necessary packages
library(sp)
library(INLA)
library(maptools)
library(spatstat)
library(rgeos)
library(RColorBrewer)
library(ggmap)
library(rgdal)



###################################################################
# Stack for all models except the one with SPDE 1D effects
###################################################################

#Load data 
load("Data/Simulated-CaseControldata-github.RData")
load("Data/census_points_mesh.RData")


######## A bit of data treatment

#Asigning convinient names 
controls_ppp <- cls.sim
lung_ppp <- lung.sim
sto_ppp <- stom.sim
kid_ppp <- kid.sim


#Taking only the coordinates of each point pattern
pts.controls <- as.matrix(coords(controls_ppp))
pts.lung <- as.matrix(coords(lung_ppp))
pts.sto  <- as.matrix(coords(sto_ppp))
pts.kid  <- as.matrix(coords(kid_ppp))

n.controls <- dim(pts.controls)[1]
n.lung <- dim(pts.lung)[1]
n.sto  <- dim(pts.sto)[1]
n.kid  <- dim(pts.kid)[1]

# Putting all point patterns in the same object
pts <- rbind(pts.controls, pts.lung, pts.sto, pts.kid)

#Number of points in pattern
n <- nrow(pts)
n.pts <- dim(pts)[1]

# Number of air polluting industries(the heavy metal industries are a subset)
n.ind <- nrow(air.ind) 

cov.dist <- matrix(NA, nrow = n.pts, ncol = n.ind)
for(i in 1:n.ind){
  p1 <- coordinates(air.ind)[i,];
  for(j in 1:n.pts){
    p2 <- t(as.matrix(pts[j,]));
    cov.dist[j, i] <- spDistsN1(p2, p1);
    rm(p2)
  }
  rm(p1)
}
cov.dist <- as.data.frame(cov.dist)
colnames(cov.dist) <- c("i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8", 
                        "i9", "i10", "i11", "i12", "i13")
head(cov.dist)


# An approximation of the density using spatstat package

#Creating a window to estimate the density with spatstat
pts_ppp <- ppp(pts[, 1], pts[, 2], window = cls.sim$window)
dens <- density(pts_ppp, dimyx = c(30, 30))

# Estimation of the density of both cases and controls with 
# density function of spatstat
pts.dens <- 
  as(density(controls_ppp, xy = dens), "SpatialGridDataFrame")
pts.dens <- 
  cbind(pts.dens,as(density(lung_ppp, xy = dens), "SpatialGridDataFrame"))
pts.dens <- 
  cbind(pts.dens,as(density(sto_ppp, xy = dens), "SpatialGridDataFrame"))
pts.dens <- 
  cbind(pts.dens,as(density(kid_ppp, xy = dens), "SpatialGridDataFrame"))

names(pts.dens) <- c("dens.controls", "dens.lung", 
  "dens.stomach", "dens.kidney")

pts.dens <- as(pts.dens, "SpatialPixelsDataFrame")

# Create mesh
# mesh in the data file
# IMPORTANT: The mesh is created in file 01_...R when census data is
#    matched to the data.
#mesh <- inla.mesh.2d(loc.domain = bdy, 
#  max.edge = c(0.4, 0.4), offset = c(0.2, 0.2))

#number of initial nodes
nv <- mesh$n

# Create spde basis functions using PC-priors
spde <- inla.spde2.pcmatern(
  # Mesh and smoothness parameter
  mesh = mesh, alpha = 2,
  # P(practic.range < 5) = 0.95
  prior.range = c(5, 0.95),
  # P(sigma > 10) = 0.01
  prior.sigma = c(10, 0.01)) 

#Index of sptial field
s.index1 <- 
  inla.spde.make.index(name = "spatial.field1", n.spde = spde$n.spde)
s.index2 <- 
  inla.spde.make.index(name = "spatial.field2", n.spde = spde$n.spde)
s.index3 <- 
  inla.spde.make.index(name = "spatial.field3", n.spde = spde$n.spde)
s.index4 <- 
  inla.spde.make.index(name = "spatial.field4", n.spde = spde$n.spde)

#Voronoi triangulation
require(deldir)
dd = deldir(mesh$loc[,1],mesh$loc[,2])
#create a list of tiles in a tessellation
mytiles = tile.list(dd)

if (!require(gpclib)) install.packages("gpclib", type="source")
#require(gpclib)
pl.study = as(bdy, "gpc.poly") #class for polygons
area.poly(pl.study) #compute the area of the polygon

#compute weight as area of the polygon
#given as interaction between voronoi tiles and domain polygon
w = unlist(lapply(mytiles,
  function(p) area.poly(intersect(as(cbind(p$x,p$y),"gpc.poly"), pl.study))
  )
)
sum(w)


#Calculating distance between the industries and the mesh points
grid.pts <- mesh$loc[, 1:2]
n.grid.pts <- dim(grid.pts)[1]


# Number of air polluting industries(the heavy metal industries are a subset)
n.ind <- nrow(air.ind) 

grid.dist <- matrix(NA, nrow = nv, ncol = n.ind)
for(i in 1:n.ind){
  p1 <- coordinates(air.ind)[i,];
  for(j in 1:n.grid.pts){
    p2 <- t(as.matrix(grid.pts[j,]));
    grid.dist[j, i] <- spDistsN1(p2, p1);
    rm(p2)
  }
  rm(p1)
}
grid.dist <- as.data.frame(grid.dist)
colnames(grid.dist) <- c("i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8", 
  "i9", "i10", "i11", "i12", "i13")
head(grid.dist)

#Computing distances between the industries and the mesh points
n.pred.pts <- dim(pts.dens)[1]

pred.dist <- matrix(NA, nrow = n.pred.pts, ncol = n.ind)
for(i in 1:n.ind){
  p1 <- coordinates(air.ind)[i,];
  for(j in 1:n.pred.pts){
    p2 <- t(as.matrix(coordinates(pts.dens)[j,]));
    pred.dist[j, i] <- spDistsN1(p2, p1);
    rm(p2)
  }
  rm(p1)
}
pred.dist <- as.data.frame(pred.dist)
colnames(pred.dist) <- c("i1", "i2", "i3", "i4", "i5", "i6", "i7",
  "i8", "i9", "i10", "i11", "i12", "i13")
head(pred.dist)

###################Function to stack a point pattern

#pts: Points in the point pattern
#mesh: Mesh to use 
#s.index: s.index for mesh
#w: weights for expected cases
#n.pp: Number of points patterns
#idx.pts: Index of this point pattern
#tag: Tag to use
#pts.pred: Points for prediction
#bl: argument indicating if the pp is the baseline pp (TRUE) or not (FALSE)
#c.dist: distance between the observed points and the covariables
#g.dist: distance between grid points and covariables
#p.dist: distance between prediction points and covariables
#cf.cov: matrix with confounding factor values for each point of the pp
#cf.mesh: matrix with confounding factor values for the mesh vertices
#cf.pred: matrix with confounding factor values the prediction grid points
stack.pp.cfcov <- function(pts, mesh, s.index, w, n.pp, idx.pts, tag, 
  pts.pred, bl, c.dist, g.dist, p.dist, cf.cov, cf.mesh){#, cf.pred) {
  
  #Number of points in pp
  n <- nrow(pts)
  
  #Number of points in mesh
  nv <- mesh$n
  
  y.pp = rep(0:1, c(nv, n))
  e.pp = c(w, rep(0, n))
  
  lmat <- inla.spde.make.A(mesh, pts)
  
  imat <- Diagonal(nv, rep(1, nv))
  A.pp <- rBind(imat, lmat)
  
  y <- matrix(NA, nrow = nv + n, ncol = n.pp)
  e <- y
  
  y[, idx.pts] <- y.pp
  #e[, idx.pts] <- e.pp
  
  A.pred <- inla.spde.make.A (mesh = mesh, loc = pts.pred)
  
  if(bl ==TRUE){
    Intercept <- list(Intercept = rep(1, n+nv))
    names(Intercept) <- paste0("Intercept.", idx.pts)
    
    Intercept.pred <- list(Intercept = rep(1, length(pts.dens)))
    names(Intercept.pred) <- paste0("Intercept.", idx.pts)
    
    stk <- inla.stack(
      data = list(y = y, e = e.pp),
      A = list(A.pp, 1),
      effects = list(s.index, Intercept=Intercept),
      tag = tag)
    
    stk.pred <- inla.stack(
      data = list(y = matrix(NA, nrow = nrow(pts.pred), ncol = n.pp)),
      A = list(A.pred, 1),
      effects = list(s.index, Intercept= Intercept.pred),
      tag = paste0(tag, ".pred"))
  }
  else{
    
    mat <- rbind(g.dist, c.dist); mat <- cbind(rep(1, nv+n), mat)
    names(mat)[1] <- paste0("Intercept.", idx.pts)
    for(i in 2:ncol(mat)){names(mat)[i] <- paste0("var", i-1, ".", idx.pts)}
    cen <- rbind(cf.mesh, cf.cov); mat <- cbind(mat, cen)
    for(i in (ncol(mat)-ncol(cen)+1):ncol(mat))
      {names(mat)[i] <- paste0(names(mat)[i], ".", idx.pts)}
    
    mat.pred <- cbind(Intercept= rep(1, nrow(pts.dens)), p.dist)
    names(mat.pred)[1] <- paste0("Intercept.", idx.pts)
    for(i in 2:ncol(mat.pred)){
      names(mat.pred)[i] <- paste0("var", i-1, ".", idx.pts)}

    
    l <- list(s.index, 1:nv, mat)
    names(l)[2] <- paste0("cont.copy", idx.pts)
    
    stk <- inla.stack(
      data = list(y = y, e = e.pp),
      A = list(A.pp, A.pp, 1),
      effects = l,
      tag = tag)
    
    l <- list(s.index, 1:nv, mat.pred) 
    names(l)[2] <- paste0("cont.copy", idx.pts)
    
    stk.pred <- inla.stack(
      data = list(y = matrix(NA, nrow = nrow(pts.pred), ncol = n.pp)),
      A = list(A.pred, A.pred, 1),
      effects = l,
      tag = paste0(tag, ".pred"))
  }
  
  return(inla.stack(stk, stk.pred))
}



stk.controls <- stack.pp.cfcov(pts.controls, mesh, s.index1, w, 4, 1, 
  "pp.controls", coordinates(pts.dens), TRUE, cov.dist[1:n.controls, ], 
  grid.dist, pred.dist)

cov.lung <- cov.dist[(n.controls+1):(n.controls+n.lung), ]
stk.lung <-  stack.pp.cfcov(pts.lung, mesh, s.index2, w, 4, 2, "pp.lung",
  coordinates(pts.dens), FALSE, cov.lung, grid.dist, pred.dist, 
  census.lung[, -1], census.mesh[, -1])#, censo.pred[, -1])


cov.sto <- cov.dist[(n.controls+n.lung+1):(n.controls+n.lung+n.sto),]
stk.sto <-  stack.pp.cfcov(pts.sto, mesh, s.index3, w, 4, 3, "pp.stomach",
  coordinates(pts.dens), FALSE, cov.sto, grid.dist, pred.dist,
  census.stom[, -1], census.mesh[, -1])#, censo.pred[, -1])

cov.kid <- 
  cov.dist[(n.controls+n.lung+n.sto+1):(n.controls+n.lung+n.sto+n.kid),]
stk.kid <-  stack.pp.cfcov(pts.kid , mesh, s.index4, w, 4, 4, "pp.kidney",
  coordinates(pts.dens), FALSE, cov.kid, grid.dist, pred.dist,
  census.kidney[, -1], census.mesh[, -1])#, censo.pred[, -1])


#Merging all the info in one stack
join.stack <- inla.stack(stk.controls, stk.lung, stk.sto, stk.kid)

# remove non-relevant objects
rm(bdy, pts_ppp, dens, w, stack.pp.cfcov, pts, 
   n.controls, dd, mytiles, n, nv, pl.study, stk.lung, stk.controls,
   stk.kid, stk.sto, n.lung, n.kid, n.sto, i, j)


save.image("Data/Stack-LinRW-models.RData")

rm(list = ls())

###################################################################


###################################################################
# Stack for models with SPDE 1D effects
###################################################################

#Load data 
load("Data/Simulated-CaseControldata-github.RData")
load("Data/census_points_mesh.RData")


######### A bit of data treatment

#Asigning convinient names 
controls_ppp <- cls.sim
lung_ppp <- lung.sim
sto_ppp <- stom.sim
kid_ppp <- kid.sim


#Taking only the coordinates of each point pattern
pts.controls <- as.matrix(coords(controls_ppp))
n.controls <- dim(pts.controls)[1]
pts.lung <- as.matrix(coords(lung_ppp))
n.lung <- dim(pts.lung)[1]
pts.sto  <- as.matrix(coords(sto_ppp))
n.sto  <- dim(pts.sto)[1]
pts.kid  <- as.matrix(coords(kid_ppp))
n.kid  <- dim(pts.kid)[1]

# Putting all point patterns in the same object
pts <- rbind(pts.controls, pts.lung, pts.sto, pts.kid)

#Number of points in pattern
n.pts <- dim(pts)[1]

# Number of air polluting industries(the heavy metal industries are a subset)
n.ind <- nrow(air.ind)

cov.dist <- matrix(NA, nrow = n.pts, ncol = n.ind)
for(i in 1:n.ind){
  p1 <- coordinates(air.ind)[i,];
  for(j in 1:n.pts){
    p2 <- t(as.matrix(pts[j,]));
    cov.dist[j, i] <- spDistsN1(p2, p1);
    rm(p2)
  }
  rm(p1)
}
cov.dist <- as.data.frame(cov.dist)
colnames(cov.dist) <- c("i1", "i2", "i3", "i4", "i5", "i6", "i7",
  "i8", "i9", "i10", "i11", "i12", "i13")
head(cov.dist)
sum(cov.dist>0)

#Creating a window to estimate the density with spatstat
pts_ppp <- ppp(pts[, 1], pts[, 2], window = cls.sim$window)
dens <- density(pts_ppp, dimyx = c(30, 30))

pts.dens <- 
  as(density(controls_ppp, xy = dens), "SpatialGridDataFrame")
pts.dens <- 
  cbind(pts.dens,as(density(lung_ppp, xy = dens), "SpatialGridDataFrame"))
pts.dens <- 
  cbind(pts.dens,as(density(sto_ppp, xy = dens), "SpatialGridDataFrame"))
pts.dens <- 
  cbind(pts.dens,as(density(kid_ppp, xy = dens), "SpatialGridDataFrame"))

names(pts.dens) <- c("dens.controls", "dens.lung", 
  "dens.stomach", "dens.kidney")

pts.dens <- as(pts.dens, "SpatialPixelsDataFrame")

#Create mesh
mesh <- inla.mesh.2d(loc.domain = bdy, 
  max.edge = c(0.4, 0.4), offset = c(0.2, 0.2))

#number of initial nodes
nv <- mesh$n

# Create spde basis functions using PC-priors
spde <- inla.spde2.pcmatern(
  # Mesh and smoothness parameter
  mesh = mesh, alpha = 2,
  # P(practic.range < 5) = 0.95
  prior.range = c(5, 0.95),
  # P(sigma > 10) = 0.01
  prior.sigma = c(10, 0.01)) 


#Index of spatial field
s.index1 <- 
  inla.spde.make.index(name = "spatial.field1", n.spde = spde$n.spde)
s.index2 <- 
  inla.spde.make.index(name = "spatial.field2", n.spde = spde$n.spde)
s.index3 <- 
  inla.spde.make.index(name = "spatial.field3", n.spde = spde$n.spde)
s.index4 <- 
  inla.spde.make.index(name = "spatial.field4", n.spde = spde$n.spde)

## Here we define the SPDE for the covariate effect using 
# a SPDE aproximation of one dimension
mesh_cov <- inla.mesh.1d(seq(0, round(max(cov.dist[-(1:3000), ])),
  length=20), degree=2, boundary="free")

#Here we create the spatial index for the covariate SPDE effect 
spde_cov <- inla.spde2.matern(mesh = mesh_cov, alpha=2, constr = TRUE)


#Voronoi triangulation
require(deldir)
dd = deldir(mesh$loc[,1],mesh$loc[,2])
#create a list of tiles in a tessellation
mytiles = tile.list(dd)

if (!require(gpclib)) install.packages("gpclib", type="source")
#require(gpclib)
pl.study = as(bdy, "gpc.poly") #class for polygons
area.poly(pl.study) #compute the area of the polygon

## Compute weight as area of the polygon given as interaction 
# between voronoi tiles and domain polygon
w = unlist(lapply(mytiles,
  function(p) area.poly(intersect(as(cbind(p$x,p$y),"gpc.poly"), pl.study)
                  ))
)
sum(w)

#Number of points in pattern
n <- nrow(pts)


#Calculating distance between the industries and mesh points
grid.pts <- mesh$loc[, 1:2]
n.grid.pts <- dim(grid.pts)[1]

# Number of air polluting industries(the heavy metal industries are a subset)
n.ind <- nrow(air.ind)

grid.dist <- matrix(NA, nrow = nv, ncol = n.ind)
for(i in 1:n.ind){
  p1 <- coordinates(air.ind)[i,];
  for(j in 1:n.grid.pts){
    p2 <- t(as.matrix(grid.pts[j,]));
    grid.dist[j, i] <- spDistsN1(p2, p1);
    rm(p2)
  }
  rm(p1)
}
grid.dist <- as.data.frame(grid.dist)
colnames(grid.dist) <- c("i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8",
  "i9", "i10", "i11", "i12", "i13")
grid.dist <- grid.dist

#Calculating distance between the industries and the mesh points
n.pred.pts <- dim(pts.dens)[1]

pred.dist <- matrix(NA, nrow = n.pred.pts, ncol = n.ind)
for(i in 1:n.ind){
  p1 <- coordinates(air.ind)[i,];
  for(j in 1:n.pred.pts){
    p2 <- t(as.matrix(coordinates(pts.dens)[j,]));
    pred.dist[j, i] <- spDistsN1(p2, p1);
    rm(p2)
  }
  rm(p1)
}
pred.dist <- as.data.frame(pred.dist)
colnames(pred.dist) <- c("i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8",
  "i9", "i10", "i11", "i12", "i13")
pred.dist <- pred.dist
rm(i, j)


## Here we design a function which create the suitable stacks

#Function to stack a point pattern
#pts: Points in the point pattern
#mesh: Mesh to use in the spatial 2D approximation
#mesh_cov: Mesh for the covariate effect SPDE in 1D
#s.index: s.index for mesh
#w: weights for expected cases
#n.pp: Number of points patterns
#idx.pts: Index of this point pattern
#tag: Tag to use
#pts.pred: Points for prediction
#bl: argument indicating if the pp is the baseline pp (TRUE) or not (FALSE)
#c.dist: distance between the observed points and the covariables
#g.dist: distance between grid points and covariables
#p.dist: distance between prediction points and covariables
#cf.cov: matrix with confounding factor values for each point of the pp
#cf.mesh: matrix with confounding factor values for the mesh vertices
#cf.pred: matrix with confounding factor values for the grid points
stack.pp.cfcov.SPDE1 <- function(pts, mesh, mesh_cov, spde_cov, s.index,
  w, n.pp, idx.pts, tag, pts.pred, bl, c.dist, g.dist, p.dist, 
  cf.cov, cf.mesh) {
  
  #Number of points in pp
  n <- nrow(pts)
  
  #Number of points in mesh
  nv <- mesh$n
  
  y.pp = rep(0:1, c(nv, n))
  e.pp = c(w, rep(0, n))
  
  y <- matrix(NA, nrow = nv + n, ncol = n.pp)
  e <- y
  
  y[, idx.pts] <- y.pp
  #e[, idx.pts] <- e.pp
  
  # Building A for the 2D mesh
  lmat <- inla.spde.make.A(mesh, pts)
  
  imat <- Diagonal(nv, rep(1, nv))
  A.pp <- rBind(imat, lmat)
  
  # Building A for the prediction
  A.pred <- inla.spde.make.A (mesh = mesh, loc = pts.pred)
  
  
  if(bl ==TRUE){
    Intercept <- list(Intercept = rep(1, n+nv))
    names(Intercept) <- paste0("Intercept.", idx.pts)
    
    Intercept.pred <- list(Intercept = rep(1, length(pts.dens)))
    names(Intercept.pred) <- paste0("Intercept.", idx.pts)
    
    stk <- inla.stack(data = list(y = y, e = e.pp), 
      A = list(A.pp, 1),
      effects = list(s.index, Intercept=Intercept),
      tag = tag)
    
    stk.pred <- inla.stack(
      data = list(y = matrix(NA, nrow = nrow(pts.pred), ncol = n.pp)),
      A = list(A.pred, 1),
      effects = list(s.index, Intercept= Intercept.pred),
      tag = paste0(tag, ".pred"))
  }
  else{
    
    ## Here we create the list of the projector matrices 
    # for the diferent spde effects
    mat <- rbind(g.dist, c.dist); 
    A_list <-list(A.pp, A.pp)
    for(p in 1:ncol(mat)){A_cov <- inla.spde.make.A(mesh_cov, mat[, p])
    A_list <- c(A_list, A_cov); rm(A_cov)}
    A_list <- c(A_list, 1)#rep(1, 5))
    
    #Here we create the list of the spde effects 
    l1 <- list(s.index, 1:nv)
    for(p in 1:13){l1 <- c(l1, 
      list(inla.spde.make.index(name = paste0("ind_var", p, ".", idx.pts),
        n.spde = spde_cov$n.spde))) }
    
    cen <- rbind(cf.mesh, cf.cov); cen <- cbind(rep(1, nv+n), cen)
    names(cen)[1] <- paste0("Intercept.", idx.pts)
    for(i in 2:ncol(cen)){
      names(cen)[i] <- paste0(names(cen)[i], ".", idx.pts)}
    
    l <- c(l1, list(cen)); 
    names(l)[-length(l)] <- NA; names(l)[2] <- paste0("cont.copy", idx.pts)
    
    stk <- inla.stack(data = list(y = y, e = e.pp),
      A = A_list,
      effects = l,
      tag = tag)
    
    ## Here we create the list of the projector matrices for 
    # the diferent spde effects
    mat <- p.dist; 
    A.pred_list <-list(A.pred, A.pred)
    for(p in 1:13){A_cov <- inla.spde.make.A(mesh_cov, mat[, p])
    A.pred_list <- c(A.pred_list, A_cov); rm(A_cov)}
    A.pred_list <- c(A.pred_list, 1)
    
    l <- c(l1, list(rep(1, nrow(pts.dens)))) 
    names(l)[length(l)] <- paste0("Intercept.", idx.pts)
    names(l)[2] <- paste0("cont.copy", idx.pts)
    
    
    stk.pred <- inla.stack(
      data = list(y = matrix(NA, nrow = nrow(pts.pred), ncol = n.pp)),
      A = A.pred_list,
      effects = l,
      tag = paste0(tag, ".pred"))
  }
  
  return(inla.stack(stk, stk.pred))
}

stk.controls <- stack.pp.cfcov.SPDE1(pts.controls, mesh, mesh_cov, 
  spde_cov, s.index1, w, 4, 1, "pp.controls", coordinates(pts.dens),
  TRUE, cov.dist[1:n.controls, ], grid.dist, pred.dist)

cov.lung <- cov.dist[(n.controls+1):(n.controls+n.lung), ]
stk.lung <-  stack.pp.cfcov.SPDE1(pts.lung, mesh, mesh_cov, spde_cov, 
  s.index2, w, 4, 2, "pp.lung", coordinates(pts.dens), FALSE, cov.lung, 
  grid.dist, pred.dist, census.lung[, -1], census.mesh[, -1])

cov.sto <- cov.dist[(n.controls+n.lung+1):(n.controls+n.lung+n.sto),]
stk.sto <-  stack.pp.cfcov.SPDE1(pts.sto, mesh, mesh_cov, spde_cov, 
  s.index3, w, 4, 3, "pp.stomach", coordinates(pts.dens), FALSE, cov.sto,
  grid.dist, pred.dist, census.stom[, -1], census.mesh[, -1])

cov.kid <- 
  cov.dist[(n.controls+n.lung+n.sto+1):(n.controls+n.lung+n.sto+n.kid),]
stk.kid <-  stack.pp.cfcov.SPDE1(pts.kid, mesh, mesh_cov, spde_cov, 
  s.index4, w, 4, 4, "pp.kidney", coordinates(pts.dens), FALSE, cov.kid,
  grid.dist, pred.dist, census.kidney[, -1], census.mesh[, -1])

#Merging the data in one stack
join.stack <- inla.stack(stk.controls, stk.lung, stk.sto, stk.kid)

#Remove non-relevant objects
rm(bdy, pts_ppp, dens, w, stack.pp.cfcov.SPDE1, pts, n.controls, dd, 
   mytiles, n, nv, pl.study, stk.lung, stk.controls,
   stk.kid, stk.sto, n.lung, n.kid, n.sto)


save.image("Data/Stack-SPDE1-Models.RData")

rm(list = ls())

###################################################################
