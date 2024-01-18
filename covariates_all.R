library(randomForest)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(exactextractr)
library(raster)
library(snowfall)
library(viridisLite)
library(gbm)
library(spatialRF)
#library(MODIStsp)
library(sf)
library(geodata)
library(ecospat)
library(fields)
library(ncdf4)

setwd("D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/Raster") # window
source("fixNAs_adj.R")
source("D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/plotRaster.R")
options(scipen = 999)
chik_foi<- read.csv("chik_foi.csv")

raster.scale <- function(ras, logT = TRUE){
  vec = as.vector(ras)
  if(logT){
    zinf <- min(vec[vec > 0], na.rm = T)
    print(paste("zinf:", zinf))
    vec = log(vec + 0.5 * zinf)
    print(head(vec))
  }
  vec = scale(vec)
  values(ras) = vec
  return(ras)
}

# 1. covariate1: mean temperature ----------------------------------------------
## 2010-2020 
temp <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/Tmean_TerraClim_2010_2020_005dg_masked_.tif"
# inspect data distribution to determine if logT should be used 
temp.hist <- raster(temp)
hist(values(temp.hist)) # left skewed
temp <- raster.scale(raster(temp), logT = FALSE)
#temp <- raster(temp)
plot(temp)

points(p_covs$Longitude, p_covs$Latitude, pch = 20, col = 'red')
temp_colors <- colorRampPalette(c("blue", "green", "yellow", "red"))(length(unique(p_covs$Temp)))
point_colors <- temp_colors[as.factor(p_covs$Temp)]
points(p_covs$Longitude, p_covs$Latitude, pch = 20, col = point_colors)
image.plot(legend.only = TRUE, zlim = c(min(p_covs$Temp), max(p_covs$Temp)), col = temp_colors, 
           horizontal = FALSE, legend.lab = "Temperature")

# 2. covariate2: mean precipitation --------------------------------------------
## 2010-2020 
precip <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/PRCP_TerraClim_2010_2020_005dg_masked_.tif"
# inspect data distribution to determine if logT should be used 
precip.hist <- raster(precip)
hist(values(precip.hist)) # right skewed
precip <- raster.scale(raster(precip)) # scaling 
#precip <- raster(precip) # not scaling
plot(precip)

prcp_colors <- colorRampPalette(c("blue", "green", "yellow", "red"))(length(unique(p_covs$PRCP)))
point_colors <- temp_colors[as.factor(p_covs$PRCP)]
points(p_covs$Longitude, p_covs$Latitude, pch = 20, col = point_colors)
image.plot(legend.only = TRUE, zlim = c(min(p_covs$PRCP), max(p_covs$PRCP)), col = temp_colors, 
           horizontal = FALSE, legend.lab = "PRCP")

# 3. covariate3: elevation -----------------------------------------------------
## 1970-2000 average WorldClim
elev <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/wc2.1_5m_elev.tif"
elev <- raster.scale(raster(elev), logT = FALSE)
plot(elev)

elev_colors <- colorRampPalette(c("blue", "green", "yellow", "red"))(length(unique(p_covs$Elev)))
point_colors <- temp_colors[as.factor(p_covs$Elev)]
points(p_covs$Longitude, p_covs$Latitude, pch = 20, col = point_colors)
image.plot(legend.only = TRUE, zlim = c(min(p_covs$Elev), max(p_covs$Elev)), col = temp_colors, 
           horizontal = FALSE, legend.lab = "Elev")

# 4. covariate4: population density --------------------------------------------
## 2000, 2005, 2010, 2015, 2020 (all years combined average; SEDAC)
pop_dens <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/pop_dens.nc"
pop_dens <- raster.scale(raster(pop_dens))
#pop_dens <- raster(pop_dens)
plot(pop_dens)

pop_colors <- colorRampPalette(c("blue", "green", "yellow", "red"))(length(unique(p_covs$Pop_dens)))
point_colors <- temp_colors[as.factor(p_covs$Pop_dens)]
points(p_covs$Longitude, p_covs$Latitude, pch = 20, col = point_colors)
image.plot(legend.only = TRUE, zlim = c(min(p_covs$Pop_dens), max(p_covs$Pop_dens)), col = temp_colors, 
           horizontal = FALSE, legend.lab = "pop")

# 5. covariate5: GDP national --------------------------------------------------
## 2009-2019
gdp_nat <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/GDP_2009_2019_National_masked.tif"
gdp_nat <- raster.scale(raster(gdp_nat))
#gdp_nat <- raster(gdp_nat)
plot(gdp_nat)


# 6. covariate6: Albopictus ----------------------------------------------------
## 2020
albo <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/Albopictus_mean_2020_rcp60_spreadXsuit_masked_.tif"
albo <- raster.scale(raster(albo), logT = FALSE)
#albo <- raster(albo)
plot(albo)

# 7. covariate7: Aegypti -------------------------------------------------------
## 2020
aegyp <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/Aegypti_mean_2020_rcp60_spreadXsuit_masked_.tif"
aegyp <- raster.scale(raster(aegyp), logT = FALSE)
#aegyp <- raster(aegyp)
plot(aegyp)

# 8. covariate8: NDVI ----------------------------------------------------------
## 2010-2020
ndvi <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/NDVI_2010_2020_005dg_masked_.tif"
ndvi <- raster.scale(raster(ndvi), logT = FALSE)
#ndvi <- raster(ndvi)
plot(ndvi)

# 9. aligning covariates -------------------------------------------------------
# extent, ncells, resolutions
temp   <- projectRaster(temp, elev, method="bilinear")
precip <- projectRaster(precip, elev, method="bilinear")
gdp_nat   <- projectRaster(gdp_nat, elev, method="bilinear")
albo   <- projectRaster(albo, elev, method="bilinear")
aegyp   <- projectRaster(aegyp, elev, method="bilinear")
ndvi   <- projectRaster(ndvi, elev, method="bilinear")
pop_dens <- projectRaster(pop_dens, elev, method="bilinear")

#10. stack covariates ----------------------------------------------------------
covlist <- stack(temp, precip, pop_dens, gdp_nat, albo, aegyp, ndvi, elev)
names(covlist) <- c("Temp", "PRCP", "Pop_dens", "GDP", "Albo", "Aegyp", "NDVI", "Elev")

#11. pcov data -----------------------------------------------------------------
# extract covariate values
p_covs <- extract(covlist, data.frame(chik_foi$long, chik_foi$lat), df= T, na.rm=TRUE)
p_covs$FOI <- chik_foi$mfoi
p_covs$Longitude <- chik_foi$long
p_covs$Latitude <- chik_foi$lat
point_idx <- which(apply(p_covs[, names(covlist)], 1, function(row) any(is.na(row))))
p_covs <- fixNAs(p_covs, covlist)

#11. save data------------------------------------------------------------------
save("temp", "precip", "elev", "pop_dens", "gdp_nat", "albo", "aegyp", "ndvi", "covlist", "p_covs", file = "covariates.RData")

