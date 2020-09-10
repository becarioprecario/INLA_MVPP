############################################################################
######################### Adding the census data.  #########################
############################################################################

library(sp)
library(rgdal)
library(INLA)

#Census data
load("Data/mapa_alcala_UTM_rescaled.RData")
#Controls and cases
load("Data/Simulated-CaseControldata-github.RData")


# Reproject using km
mapa.alcala.UTM.rescaled <- spTransform(mapa.alcala.UTM.rescaled, 
  CRS("+proj=utm +zone=30 +init=epsg:23030 +towgs84=-87,-96,-120,0,0,0,0 +ellps=intl +units=km +no_defs"))

# Test map projection
plot(stom.sim, pch = ".", col = "red")
plot(mapa.alcala.UTM.rescaled, add = TRUE)

# Variables of interest; COD is not used in the models
vars <- c("COD", "PARO2059", "SOCMEDIA", "ESTMEDIO", "ESTPRE")

# Get census data from map
library(maptools)
get_census_data <- function(pts) {
  pts.sp <- as(pts, "SpatialPoints")
  proj4string(pts.sp) <- CRS(proj4string(mapa.alcala.UTM.rescaled))
#SpatialPoints(pts[, c("COORD_X", "COORD_Y")], 
#    proj4string = CRS(proj4string(mapa.alcala.UTM)))
  datos <- over(pts.sp, mapa.alcala.UTM.rescaled)

  datos[, vars]
}

# Data
# Controls
census.cnts <- get_census_data(cls.sim)
# Stomach cancer
census.stom <- get_census_data(stom.sim)
# Lung cancer
census.lung <- get_census_data(lung.sim)
# Kidney cancer
census.kidney <- get_census_data(kid.sim)

#Check values
summary(census.cnts)
summary(census.stom)
summary(census.lung)
summary(census.kidney)


# Get census data for mesh points
mesh <- inla.mesh.2d(loc.domain = bdy, 
  max.edge = c(0.4, 0.4), offset = c(0.2, 0.2))

pts.malla <- mesh$loc[, 1:2] 

malla.sp <- SpatialPoints(pts.malla,
    proj4string = CRS(proj4string(mapa.alcala.UTM.rescaled)))
datos.malla <- over(malla.sp, mapa.alcala.UTM.rescaled)

# NA's are replaced by the value at the neares point
notNA <- which(!is.na(datos.malla[, 1]))
for(i in which(is.na(datos.malla[, 1]))) {
  # aux takes values between 1 and length (notNA)
  aux <- which.min( (pts.malla[i, 1] - pts.malla[notNA, 1])^2 +
    (pts.malla[i, 2] - pts.malla[notNA, 2])^2 )
  datos.malla[i, ] <- datos.malla[notNA[aux], ] 
}

census.mesh <- datos.malla [, vars]


save(file = "Data/census_points_mesh.RData",
  list = c("census.cnts", "census.stom", "census.lung", "census.kidney",
    "mesh", "census.mesh"))


