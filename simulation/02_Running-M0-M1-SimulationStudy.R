############################################################################
#### Simulation study IIa: Running models with lineal covariate effect. ####
############################################################################

# In this script we run model 0 and  model 1, i.e. the ones which not 
# consider the exposure effect. This models will be computed using the 
# different groups of cases in order to study the ability of each model
# to highlight the effect depending on the available sample size. 

#Loading necessary libraries
library(sp)
library(INLA)
library(maptools)
library(spatstat)
library(rgeos)
library(RColorBrewer)
library(ggmap)
library(rgdal)
library(lattice)
library(gridExtra)
library(xtable)

for(m in 1:6){
  #Load data with variability (for phi = 3)
  load("SimulatedDataphi3.RData")
  
  #Setting number of cases
  n.cases <- c(50, 100, 300, 500, 1000, 2000)
  
  ########################################################################
  # Establishing number of cases and controls
  controls_ppp <- cls.sim
  cases_ppp <- cases[1:n.cases[m], ]
  
  #Taking the coordinates of the whole grup of points
  pts.controls <- as.matrix(coords(controls_ppp))
  n.controls <- dim(pts.controls)[1]
  pts.cases <- as.matrix(coords(cases_ppp))
  pts <- rbind(pts.controls, pts.cases)
  n.pts <- dim(pts)[1]
  
  pts.ind <- as.matrix(coordinates(air.ind)[, ])
  n.ind <- dim(pts.ind)[1]
  
  
  cov.dist <- matrix(NA, nrow = n.pts, ncol = n.ind)
  for(i in 1:n.ind){
    p1 <- pts.ind[i,];
    for(j in 1:n.pts){
      p2 <- t(pts[j,]);
      cov.dist[j, i] <- spDistsN1(p2, p1);
      rm(p2)
    }
    rm(p1)
  }
  
  cov.dist <- as.data.frame(cov.dist)
  colnames(cov.dist) <- c("i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8",
    "i9", "i10", "i11", "i12", "i13")[1:ncol(cov.dist)]
  head(cov.dist)
  sum(cov.dist<0)
  
  # Estimate the density for both cases and controls 
  pts_ppp <- ppp(pts[, 1], pts[, 2], window = cls.sim$window)
  dens <- density(pts_ppp, dimyx = c(30, 30))
  
  pts.dens <- 
    as(density(controls_ppp, xy = dens), "SpatialGridDataFrame")
  pts.dens <- 
    cbind(pts.dens,as(density(cases_ppp, xy = dens), "SpatialGridDataFrame"))
  
  names(pts.dens) <- c("dens.controls", "dens.cases")
  
  pts.dens <- as(pts.dens, "SpatialPixelsDataFrame")
  
  #Create mesh
  mesh <- inla.mesh.2d(loc.domain = bdy, 
    max.edge = c(0.4, 0.4), offset = c(0.2, 0.2))

  #number of initial nodes
  nv <- mesh$n
  
  #Create spde
  spde <- inla.spde2.matern(mesh = mesh, alpha=2, constr = TRUE)
  
  #Index of sptial field
  s.index1 <- 
    inla.spde.make.index(name = "spatial.field1", n.spde = spde$n.spde)
  s.index2 <- 
    inla.spde.make.index(name = "spatial.field2", n.spde = spde$n.spde)
  
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
    function(p) area.poly(intersect(as(cbind(p$x,p$y),"gpc.poly"), pl.study)
                    ))
  )
  sum(w)
  
  #Number of points in pattern
  n <- nrow(pts)
  
  
  #Calculating distance between the industries and the mesh points
  grid.pts <- mesh$loc[, 1:2]
  n.grid.pts <- dim(grid.pts)[1]
  
  grid.dist <- matrix(NA, nrow = nv, ncol = n.ind)
  for(i in 1:n.ind){
    p1 <- pts.ind[i,];
    for(j in 1:n.grid.pts){
      p2 <- t(as.matrix(grid.pts[j,]));
      grid.dist[j, i] <- spDistsN1(p2, p1);
      rm(p2)
    }
    rm(p1)
  }
  grid.dist <- as.data.frame(grid.dist)
  colnames(grid.dist) <-  c("i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8",
    "i9", "i10", "i11", "i12", "i13")[1:ncol(grid.dist)]
  head(grid.dist)
  sum(grid.dist<0)
  
  #Calculating distance between the industries and the mesh points
  n.pred.pts <- dim(pts.dens)[1]
  
  pred.dist <- matrix(NA, nrow = n.pred.pts, ncol = n.ind)
  for(i in 1:n.ind){
    p1 <- pts.ind[i,];
    for(j in 1:n.pred.pts){
      p2 <- t(as.matrix(coordinates(pts.dens)[j,]));
      pred.dist[j, i] <- spDistsN1(p2, p1);
      rm(p2)
    }
    rm(p1)
  }
  pred.dist <- as.data.frame(pred.dist)
  colnames(pred.dist) <- c("i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8",
    "i9", "i10", "i11", "i12", "i13")[1:ncol(pred.dist)]
  head(pred.dist)
  sum(pred.dist<0)
  
  
  ##############################################################
  ####### Creating a stack with covariate effect as Rw1 effect.
  ##############################################################
  
  #Function to stack a point pattern
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
  stack.pp.M3 <- function(pts, mesh, s.index, w, n.pp, idx.pts, tag, 
    pts.pred, bl, c.dist, g.dist, p.dist) {
    
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
      for(i in 2:ncol(mat)){
        names(mat)[i] <- paste0("var", i-1, ".", idx.pts)}
      
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
  
  stk.controls <- stack.pp.M3(pts.controls, mesh, s.index1, w, 2, 1,
    "pp.controls", coordinates(pts.dens), TRUE, 
    cov.dist[1:n.controls, ], grid.dist, pred.dist)
  
  cov.cases <- cov.dist[(n.controls+1):(n.controls+n.cases[m]), ]
  stk.cases <-  stack.pp.M3(pts.cases, mesh, s.index2, w, 2, 2, 
    "pp.cases", coordinates(pts.dens), FALSE, cov.cases, 
    grid.dist, pred.dist)
  
  #Merge all data together in one stack
  join.stack <- inla.stack(stk.controls, stk.cases)
  
  ### Now we are going to run the models 
  
  ############################################################
  #Model 0: No exposure effect and no specific spatial effect
  ############################################################
  
  t1 <-Sys.time()
  E <- inla.stack.data(join.stack)$e
  
  #formula
  form <- y ~ 0 + Intercept.1 + Intercept.2 +
     f(spatial.field1, model = spde) +
     f(cont.copy2, copy = "spatial.field1", fixed = TRUE)  
  
    #Estimation
    Sys.time()
    pp.res <- inla(formula = form, 
      family = rep('poisson', 2), data = inla.stack.data(join.stack),
      control.predictor = list(A = inla.stack.A(join.stack), link = 1,
        compute = TRUE),
      E = inla.stack.data(join.stack)$e, 
      control.compute = list(dic=TRUE, waic=TRUE, cpo=TRUE, 
        mlik=TRUE, po=TRUE))
    Sys.time()
    
    summary(pp.res)
    t2 <-Sys.time()
    
    ti <- t2-t1; 
    print(paste0("Model 0 with ", n.cases[m],  " cases has used ", ti,
      " cpu time. And have finished at ", Sys.time()))
    
    pp.res.M0<- pp.res
    rm(t1, t2, E, form, pp.res)
  ############################################################

  ############################################################
  #Model 1: Not exposure effect, but specific spatial effects
  ############################################################
  
  t1 <-Sys.time()
  E <- inla.stack.data(join.stack)$e
  
  #formula
  form <- y ~ 0 + Intercept.1 + Intercept.2 +
  f(spatial.field1, model = spde) +
  f(spatial.field2, model = spde) + 
  f(cont.copy2, copy = "spatial.field1", fixed = TRUE)  
    
  #Estimation
  Sys.time()
  pp.res <- inla(formula = form, verbose = FALSE,
    family = rep('poisson', 2), data = inla.stack.data(join.stack),
    control.predictor = list(A = inla.stack.A(join.stack), link = 1,
      compute = TRUE),
    E = inla.stack.data(join.stack)$e, 
    control.compute = list(dic=TRUE, waic=TRUE, cpo=TRUE, 
      mlik=TRUE, po=TRUE))
  Sys.time()
  
  summary(pp.res)
  t2 <-Sys.time()
  
  ti <- t2-t1
  print(paste0("Model 1 with ", n.cases[m],  " cases has used ", ti,
    " cpu time. And have finished at ", Sys.time()))
  
  pp.res.M1 <- pp.res
  
  rm(t1, t2, E, form, pp.res)
  join.stack.SPDE1 <- join.stack
  rm(bdy, cov.cases, grid.pts, pts, pts.cases, pts.controls, pts.ind,
     cases.1, cases.2, cases.3, dd, dens, f.coord, i, j, join.stack, 
     mytiles, n, n.grid.pts, n.ind, n.pred.pts, n.pts, nv, phi, pl.study,
     s.index1, s.index2, spde, stk.cases, stk.controls, w, stack.pp.M3)
  
  guardar <- paste0("M0_M1_cases", n.cases[m] , "phi3.RData" )
  
  save.image(guardar)
  
  
  rm(cases, cases_ppp, cls.sim, controls_ppp, guardar, mesh, ti, 
     n.cases, n.controls, pp.res.M0, pp.res.M1, pts_ppp, pts.dens, 
     cov.dist, grid.dist, metals.ind, pred.dist, air.ind)
  
}

rm(list = ls())

