############################################################################
#### Simulation study I: Simulating data form the controls of the study ####
############################################################################

# In this script, the data of the simulation study is simulated. First of
# all, we simulate the data. For this, we estimate the intensity of the 
# real controls. Then with this intensity we define the intensity of the 
# cases adding a covariate effect. Finally, we simulate cases and controls 
# of each intensities.

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

# Phi is a scale factor
#  Change this value to create other simulated datasets
#  In this code (and the other R file) we only consider phi = 3
#  but it is straightforward to simaulate for other values of phi
#  and fit the models.
phi <- 3

#Load data with variability
load("Data-with-VariabilityAdded-1Controls-Less.RData")


#Load industries
metals.ind <- read.table(file = "focos_metales_pesados.txt", header = TRUE,
  sep = "", dec = ".", quote = "\"", row.names = NULL)

air.ind <- read.table(file = "focos_aire_contam.txt", header = TRUE,
  sep = "", dec = ".", quote = "\"", row.names = NULL, encoding = "UTF8")

#Here we are going to compute the distance between covariates and points
# We are changing from meters to km
air.ind$xed50utm <- air.ind$xed50utm/1000
air.ind$yed50utm <- air.ind$yed50utm/1000
coordinates(air.ind) <- ~ xed50utm + yed50utm
n.ind <- nrow(air.ind)
pts.ind <- coordinates(air.ind)

n.pts <- dim(pts)[1]

# REAL CONTROLS
controls <- controls_ppp

# Estimate density
d0 <- density(controls); plot(d0)

# Sample SIMULATED controls
n.controls <- 3000
set.seed(1882011)
cls.sim <- rpoispp( 5 * d0)[1:n.controls, ]; 

# Distance to the pollution source (Industry 1) 
f.coord <- coordinates(air.ind)[1, ]
p.source <- ppp(f.coord[1], f.coord[2], window = controls$window)
ps.dist <- distmap(p.source); plot(ps.dist)

# Compute inverse distance
ps.dist2 <- 1 / (ps.dist + 0.5 * max(ps.dist))
plot(ps.dist2)

#Re-scale ps.dist to have a maximum at 1 
ps.dist2 <- 1 * ps.dist2 / max(ps.dist2); plot(ps.dist2)

# Effect on distance
ef.dist <- exp(-ps.dist / phi ); plot(ef.dist)


# Intensity of cases
d1 <- d0 * ef.dist
plot(d1)

# Sample cases
n.cases <- 2000
set.seed(1882011)
cases <- rpoispp(10 * d1)[1:n.cases, ]

cases.1 <- cases[1:500, ]; 
cases.2 <- cases[1:1000, ]; 
cases.3 <- cases[1:2000, ];

png(file = "SimulatedDataphi3.png", width = 720, height = 480)
par(mfrow=c(2,4))
plot(cls.sim)
plot(cases.1)
plot(cases.2)
plot(cases.3)
plot(density(cls.sim, adjust=0.5))
plot(density(cases.1, adjust=0.5))
plot(density(cases.2, adjust=0.5))
plot(density(cases.3, adjust=0.5))
dev.off()

plot(density(cls.sim))
plot(density(cls.sim, adjust=0.5))

rm(bdy.orig, pts, pts.controls, pts.ind, pts.kid, pts.lung, pts.sto, 
   bdy.owin, bdy.SP, cont.vir, controls, controls_ppp, d0, d1, ef.dist,
   estom, kid_ppp, lung_ppp,n.cases, n.controls, n.ind, n.kid, n.lung, 
   n.pts, n.sto, p.source, ps.dist, ps.dist2, pts_ppp, pulmon, 
   rinon, sto_ppp)

save.image("SimulatedDataphi3.RData")

rm(list = ls())


