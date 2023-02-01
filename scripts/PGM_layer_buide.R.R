
#' @title Cost-benefit of conservation actions in Amazon
#' @description script to build exploratory variables from 
#' land use - land cover, secondary forest, edge, 
#' degradation (loggin and fire), temperature, precipitation, 
#' elevation and distances to road and water body
#' in Paragominas municipality - PA;
#' this set of exploratory variables is used in fit 
#' regional species distribution models 

#### setting working directory ####
setwd("~/projetos/lancaster/conserv_opportunities_jamesthomson/Rscripts")
memory.limit(1000000)

##### loading required packages ####
library(tidyverse)
library(raster)
library(rgdal)
library(rgeos)
library(sf)
library(spatialEco)
library(scales)


#### importing input rasters ####

# shapefile paragominas
pgm.shp <- readOGR(dsn = "shapes", layer = "Paragominas_Mask_R3")

# land use land cover from mapbiomas collection 7 [2010 and 2020]
# [PGM] paragominas
pgm.lulc <- stack(c("rasters/PGM/input/pgm-2010-lulc-mapbiomas-brazil-collection-70.tif",
                    "rasters/PGM/input/pgm-2020-lulc-mapbiomas-brazil-collection-70.tif"))
names(pgm.lulc) <- c("pgm.lulc.2010real", "pgm.lulc.2020real")
#checking
#pgm.lulc
#plot(pgm.lulc)
#sort(unique(values(pgm.lulc[["PGM.LULC.2010"]])))

values(pgm.lulc)[values(pgm.lulc) <= 0] = NA
# land ues land cover pixel values and codes
# 0  == NA
# 3  == Forest Formation      == Forest
# 4  == Savanna Formation     == Forest
# 9  == Forest Plantation     == Farming
# 11 == Wetland               == Non Forest Natural Formation
# 12 == Grassland             == Non Forest Natural Formation
# 15 == Pasture               == Farming
# 24 == Urban Area            == Non vegetated area
# 30 == Mining                == Non vegetated area
# 33 == River, Lake and Ocean == Water
# 39 == Soybean               == Farming
# 41 == Other Temporary Crops == Farming
# 48 == Other Perennial Crops == Farming

# isolating forest class pixels
pgm.lulc.2010.forest.class <- pgm.lulc[["pgm.lulc.2010real"]]
pgm.lulc.2010.forest.class[pgm.lulc.2010.forest.class==3] <- 1
pgm.lulc.2010.forest.class[pgm.lulc.2010.forest.class>1] <- NA

pgm.lulc.2020.forest.class <- pgm.lulc[["pgm.lulc.2020real"]]
pgm.lulc.2020.forest.class[pgm.lulc.2020.forest.class==3] <- 1
pgm.lulc.2020.forest.class[pgm.lulc.2020.forest.class>1] <- NA

#
#



#candidate areas for restoration scenarios

#isolating deforestation class pixels (crops, pasture)

deforestatio.class.list <- c(15,39,41,48)

candidate.areas.total <- pgm.2010[[1]]

values(candidate.areas.total)[values(candidate.areas.total) %in% deforestatio.class.list] = 1
values(candidate.areas.total)[values(candidate.areas.total) > 1] = 0
names(candidate.areas.total) <- "restoration.candidate.areas"
#plot(candidate.areas.total)

#select pixels based on proximity to water (<500m), slope (>25°) and proximity to forest (<1000m)

dist.river <- raster("rasters/PGM/2010_real/dist_river_pgm.tif")
values(dist.river)[values(dist.river) <= 500] = 1
values(dist.river)[values(dist.river) > 500] = NA
dist.river <- projectRaster(dist.river, crs = "+proj=longlat +datum=WGS84 +no_defs")
dist.river <- resample(dist.river, candidate.areas.water, method='ngb')
#plot(dist.river)

candidate.areas.water <- candidate.areas.total
candidate.areas.water <- mask(candidate.areas.water, dist.river)
values(candidate.areas.water)[is.na(values(candidate.areas.water))] = 0
candidate.areas.water <- mask(candidate.areas.water, pgm.shp)
#plot(candidate.areas.water)



elevation <- raster("rasters/PGM/2010_real/elevation_pgm.tif")
#plot(elevation)
slope <- terrain(elevation, opt = 'slope', unit = 'degrees', neighbors=8)
values(slope)[values(slope) < 25] = NA
values(slope)[values(slope) >= 25] = 1
slope <- projectRaster(slope, crs = "+proj=longlat +datum=WGS84 +no_defs")
slope <- resample(slope, candidate.areas.slope, method='ngb')
#plot(slope)

candidate.areas.slope <- candidate.areas.total
candidate.areas.slope <- mask(candidate.areas.slope, slope)
values(candidate.areas.slope)[is.na(values(candidate.areas.slope))] = 0
candidate.areas.slope <- mask(candidate.areas.slope, pgm.shp)
#plot(candidate.areas.slope)

#
#



# time since degradation 2010 data from RAS 
# quantitative comparison of manual inspection of satellite images
# and field observations done by two observers (TG and SN)
# see RAS environmental explanatory variable guideline document for details
# 2020 data from DETER
# 
pgm.degrad.2010 <- raster("rasters/PGM/input/pgm-2010-deg_tsince0_150m.grd")

names(pgm.degrad.2010) <- c("pgm.degrad.2010real")

#checking
#pgm.degrad.2010
#plot(pgm.degrad.2010)
#range(values(pgm.degrad.2010), na.rm=T)

# Conversion of rasters into same extent
pgm.degrad.2010 <- projectRaster(pgm.degrad.2010, crs = "+proj=longlat +datum=WGS84 +no_defs")
pgm.degrad.2010 <- resample(pgm.degrad.2010, pgm.lulc, method='ngb')


# calculating time since degradation for 2020
pgm.degrad.temp <- pgm.degrad.2010


# deter data between 2011 and 2015
library(datazoom.amazonia)

deter.2011.15 <- load_degrad(dataset = "degrad", raw_data = T, time_period = 2011:2015)

for (year in 1:5) {   #1=2011; 5=2015
  
  deter.yearx <- deter.2011.15[[year]] 
  deter.yearx <- sf:::as_Spatial(deter.yearx$geometry)
  pgm.deter.yearx <- crop(deter.yearx, extent(pgm.shp))
  
  pgm.deter.yearx <- rasterize(pgm.deter.yearx, pgm.lulc[[1]], field=999)
  pgm.deter.yearx[is.na(pgm.deter.yearx)]<-0
  pgm.deter.yearx <- mask(pgm.deter.yearx, pgm.shp)
  
  pgm.degrad.temp <- pgm.degrad.temp+1
  pgm.degrad.temp[get("pgm.deter.yearx")[] == 999] <- 0
  
  cat("\n> year", year, "done! <\n")
}


#plot(pgm.degrad.temp)
#plot(pgm.deter.yearx, add=T)


# deter data between 2016 and 2020
deter.2016.20 <- readOGR(dsn = "rasters/PGM/input", layer = "deter_public")
pgm.deter.2016.20 <- crop(deter.2016.20, extent(pgm.shp))

rm(deter.2016.20)

for (year in 2016:2020) {
  
  pgm.deter.yearx <- pgm.deter.2016.20[grep(year, pgm.deter.2016.20$VIEW_DATE),]
  pgm.deter.yearx <- rasterize(pgm.deter.yearx, pgm.lulc[[1]], field=999)
  pgm.deter.yearx[is.na(pgm.deter.yearx)]<-0
  pgm.deter.yearx <- mask(pgm.deter.yearx, pgm.shp)
  
  pgm.degrad.temp <- pgm.degrad.temp+1
  pgm.degrad.temp[get("pgm.deter.yearx")[] == 999] <- 0
  
  cat("\n> year", year, "done! <\n")
}



pgm.degrad.temp <- mask(pgm.degrad.temp, pgm.lulc.2020.forest.class)

pgm.degrad.2020 <- pgm.degrad.temp

names(pgm.degrad.2020) <- "pgm.degrad.2020real"

pgm.degrad <- stack(pgm.degrad.2010, pgm.degrad.2020)

#checking
#pgm.degrad
#plot(pgm.degrad)
#range(values(pgm.degrad[["PGM.Degrad.2010"]]), na.rm=T)

#removing undegraded forest
pgm.degrad[["pgm.degrad.2010real"]][pgm.degrad[["pgm.degrad.2010real"]]>23] <- NA
pgm.degrad[["pgm.degrad.2020real"]][pgm.degrad[["pgm.degrad.2020real"]]>33] <- NA

# isolating degraded forest class pixels
pgm.degrad.2010.forest.class <- pgm.degrad[["pgm.degrad.2010real"]]
pgm.degrad.2010.forest.class[!is.na(pgm.degrad.2010.forest.class)] <- 1
pgm.degrad.2010.forest.class<-sum(pgm.lulc.2010.forest.class, pgm.degrad.2010.forest.class, na.rm=T)
pgm.degrad.2010.forest.class[pgm.degrad.2010.forest.class<2]<-0
pgm.degrad.2010.forest.class[pgm.degrad.2010.forest.class==2]<-1

pgm.degrad.2020.forest.class <- pgm.degrad[["pgm.degrad.2020real"]]
pgm.degrad.2020.forest.class[!is.na(pgm.degrad.2020.forest.class)] <- 1
pgm.degrad.2020.forest.class<-sum(pgm.lulc.2020.forest.class, pgm.degrad.2020.forest.class, na.rm=T)
pgm.degrad.2020.forest.class[pgm.degrad.2020.forest.class<2]<-0
pgm.degrad.2020.forest.class[pgm.degrad.2020.forest.class==2]<-1


#
#



# secondary forest age from Silva Jr. et al 2020  [2010 and 2020]
# [DOI: 10.1038/s41597-020-00600-4]
# [PGM] paragominas
pgm.sfage <- stack(c("rasters/PGM/input/pgm-2010-sfage-mapbiomas-brazil-collection-60.tif",
                      "rasters/PGM/input/pgm-2020-sfage-mapbiomas-brazil-collection-60.tif"))
names(pgm.sfage) <- c("pgm.sfage.2010real", "pgm.sfage.2020real")

#checking
#pgm.sfage
#plot(pgm.sfage)
#range(values(pgm.sfage[["pgm.sfage.2010real"]]), na.rm = T)

# Conversion of rasters into same extent
pgm.sfage <- resample(pgm.sfage, pgm.lulc, method='ngb')

# isolating all secondary forest class pixels
pgm.sfage.2010.all.class <- pgm.sfage[["pgm.sfage.2010real"]]
pgm.sfage.2010.all.class[pgm.sfage.2010.all.class > 0] <- 1

# isolating secondary forest with more than two years old
pgm.sfage.2010.t2.class <- pgm.sfage[["pgm.sfage.2010real"]]
pgm.sfage.2010.t2.class[pgm.sfage.2010.all.class < 2] <- 0
pgm.sfage.2010.t2.class[pgm.sfage.2010.all.class >= 2] <- 1

# isolating secondary forest with more than ten years old
pgm.sfage.2010.t10.class <- pgm.sfage[["pgm.sfage.2010real"]]
pgm.sfage.2010.t10.class[pgm.sfage.2010.all.class < 10] <- 0
pgm.sfage.2010.t10.class[pgm.sfage.2010.all.class >= 10] <- 1

#
pgm.sfage.2020.all.class <- pgm.sfage[["pgm.sfage.2020real"]]
pgm.sfage.2020.all.class[pgm.sfage.2020.all.class > 0] <- 1

pgm.sfage.2020.t2.class <- pgm.sfage[["pgm.sfage.2020real"]]
pgm.sfage.2020.t2.class[pgm.sfage.2020.all.class < 2] <- 0
pgm.sfage.2020.t2.class[pgm.sfage.2020.all.class >= 2] <- 1

pgm.sfage.2020.t10.class <- pgm.sfage[["pgm.sfage.2020real"]]
pgm.sfage.2020.t10.class[pgm.sfage.2020.all.class < 10] <- 0
pgm.sfage.2020.t10.class[pgm.sfage.2020.all.class >= 10] <- 1

#
#



# function to make sure all raster are fully stacked
#intersect_mask <- function(x){
#  values_x <- getValues(x)
#  inter_x <- values_x %*% rep(1,nlayers(x))
#  mask <- setValues(subset(x,1),values = (inter_x>0))
#  return(mask)
#}

# stacking by year
pgm.2010 <- stack(pgm.lulc[[1]], pgm.sfage[[1]], pgm.degrad[[1]])
#pgm.2010 <- stack(mask(pgm.2010, intersect_mask(pgm.2010)))

pgm.2020 <- stack(pgm.lulc[[2]], pgm.sfage[[2]], pgm.degrad[[2]])
#pgm.2020 <- stack(mask(pgm.2020, intersect_mask(pgm.2020)))


rm(list=ls()[!ls() %in% c("pgm.shp", "pgm.2010", "pgm.2020",
                          "pgm.lulc.2010.forest.class", "pgm.lulc.2020.forest.class",
                          "pgm.sfage.2010.all.class", "pgm.sfage.2020.all.class",
                          "pgm.sfage.2010.t2.class", "pgm.sfage.2020.t2.class",
                          "pgm.sfage.2010.t10.class", "pgm.sfage.2020.t10.class",
                          "pgm.degrad.2010.forest.class", "pgm.degrad.2020.forest.class")]) #keeping only raster stack
gc()
#



#######################################################################################################################

#### setting explanatory variables ####
# [UPF] undisturbed primary forest -- pixel: 5x5 (150m); and landscape: 35x35 (1050m)
# this variable includes all forest pixels in LULC raster (value == 3)
# excluding those with age < 25 in 2010 SF raster or age < 35 in 2020 SF raster
# excluding pixels degraded

dir.create("rasters/PGM/2010_real", recursive = T)
dir.create("rasters/PGM/2020_real", recursive = T)
dir.create("rasters/PGM/2020_avoiddeforest", recursive = T)
dir.create("rasters/PGM/2020_avoiddegrad", recursive = T)
dir.create("rasters/PGM/2020_avoidboth", recursive = T)
dir.create("rasters/PGM/2020_restor_wo_avoid", recursive = T)
dir.create("rasters/PGM/2020_restor_n_avoid", recursive = T)



# scenario 2010
UPF2010<-sum(pgm.lulc.2010.forest.class, pgm.sfage.2010.all.class, pgm.degrad.2010.forest.class, na.rm = T)
UPF2010[UPF2010>1]<-0
##cheking
#unique(UPF2010[])
#plot(UPF2010)

# mean upf cover in pixel scale (150m)
UPF2010.px <- focal(UPF2010, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#UPF2010.px
#anyNA(UPF2010.px[])
#plot(UPF2010.px)

names(UPF2010.px)<-"UPFpx"
UPF2010.px[is.nan(UPF2010.px)] <- 0
UPF2010.px <- mask(UPF2010.px, pgm.shp)

#saving
writeRaster(UPF2010.px, "rasters/PGM/2010_real/UPFpx.tif", format="GTiff", overwrite=T)
writeRaster(UPF2010.px, "rasters/PGM/2020_avoidboth/UPFpx.tif", format="GTiff", overwrite=T)
writeRaster(UPF2010.px, "rasters/PGM/2020_restor_n_avoid/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (1050m)
UPF2010.ls <- focal(UPF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#UPF2010.ls
#anyNA(UPF2010.ls[])
#plot(UPF2010.ls)

names(UPF2010.ls)<-"UPFls"
UPF2010.ls[is.nan(UPF2010.ls)] <- 0
UPF2010.ls <- mask(UPF2010.ls, pgm.shp)

#saving
writeRaster(UPF2010.ls, "rasters/PGM/2010_real/UPFls.tif", format="GTiff", overwrite=T)
writeRaster(UPF2010.ls, "rasters/PGM/2020_avoidboth/UPFls.tif", format="GTiff", overwrite=T)
writeRaster(UPF2010.ls, "rasters/PGM/2020_restor_n_avoid/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[!ls() %in% c("pgm.shp", "pgm.2010", "pgm.2020",
                          "pgm.lulc.2010.forest.class", "pgm.lulc.2020.forest.class",
                          "pgm.sfage.2010.all.class", "pgm.sfage.2020.all.class",
                          "pgm.sfage.2010.t2.class", "pgm.sfage.2020.t2.class",
                          "pgm.sfage.2010.t10.class", "pgm.sfage.2020.t10.class",
                          "pgm.degrad.2010.forest.class", "pgm.degrad.2020.forest.class",
                          "UPF2010")]) #keeping only raster stack
gc()
#

# scenario avoid degradation
UPF.avoiddegrad<-sum(pgm.lulc.2020.forest.class, pgm.sfage.2020.all.class, pgm.degrad.2010.forest.class, na.rm = T)
UPF.avoiddegrad[UPF.avoiddegrad>1]<-0
##cheking
#unique(UPF.avoiddegrad[])
#plot(UPF.avoiddegrad)

# mean upf cover in pixel scale (150m)
UPF.avoiddegrad.px <- focal(UPF.avoiddegrad, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#UPF.avoiddegrad.px
#anyNA(UPF.avoiddegrad.px[])
#plot(UPF.avoiddegrad.px)

names(UPF.avoiddegrad.px)<-"UPFpx"
UPF.avoiddegrad.px[is.nan(UPF.avoiddegrad.px)] <- 0
UPF.avoiddegrad.px <- mask(UPF.avoiddegrad.px, pgm.shp)

#saving
writeRaster(UPF.avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (1050m)
UPF.avoiddegrad.ls <- focal(UPF.avoiddegrad, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#UPF.avoiddegrad.ls
#anyNA(UPF.avoiddegrad.ls[])
#plot(UPF.avoiddegrad.ls)

names(UPF.avoiddegrad.ls)<-"UPFls"
UPF.avoiddegrad.ls[is.nan(UPF.avoiddegrad.ls)] <- 0
UPF.avoiddegrad.ls <- mask(UPF.avoiddegrad.ls, pgm.shp)

#saving
writeRaster(UPF.avoiddegrad.ls, "rasters/PGM/2020_avoiddegrad/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[!ls() %in% c("pgm.shp", "pgm.2010", "pgm.2020",
                          "pgm.lulc.2010.forest.class", "pgm.lulc.2020.forest.class",
                          "pgm.sfage.2010.all.class", "pgm.sfage.2020.all.class",
                          "pgm.sfage.2010.t2.class", "pgm.sfage.2020.t2.class",
                          "pgm.sfage.2010.t10.class", "pgm.sfage.2020.t10.class",
                          "pgm.degrad.2010.forest.class", "pgm.degrad.2020.forest.class",
                          "UPF2010", "UPF.avoiddegrad")]) #keeping only raster stack
gc()
#

# scenario avoid deforestation
UPF.avoiddefor<-sum(pgm.lulc.2010.forest.class, pgm.sfage.2010.all.class, pgm.degrad.2020.forest.class, na.rm = T)
UPF.avoiddefor[UPF.avoiddefor>1]<-0
##cheking
#unique(UPF.avoiddegrad[])
#plot(UPF.avoiddegrad)

# mean upf cover in pixel scale (150m)
UPF.avoiddefor.px <- focal(UPF.avoiddefor, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#UPF.avoiddefor.px
#anyNA(UPF.avoiddefor.px[])
#plot(UPF.avoiddefor.px)

names(UPF.avoiddefor.px)<-"UPFpx"
UPF.avoiddefor.px[is.nan(UPF.avoiddefor.px)] <- 0
UPF.avoiddefor.px <- mask(UPF.avoiddefor.px, pgm.shp)

#saving
writeRaster(UPF.avoiddefor.px, "rasters/PGM/2020_avoiddeforest/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (1050m)
UPF.avoiddefor.ls <- focal(UPF.avoiddefor, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#UPF.avoiddefor.ls
#anyNA(UPF.avoiddefor.ls[])
#plot(UPF.avoiddefor.ls)

names(UPF.avoiddefor.ls)<-"UPFls"
UPF.avoiddefor.ls[is.nan(UPF.avoiddefor.ls)] <- 0
UPF.avoiddefor.ls <- mask(UPF.avoiddefor.ls, pgm.shp)

#saving
writeRaster(UPF.avoiddefor.ls, "rasters/PGM/2020_avoiddeforest/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[!ls() %in% c("pgm.shp", "pgm.2010", "pgm.2020",
                          "pgm.lulc.2010.forest.class", "pgm.lulc.2020.forest.class",
                          "pgm.sfage.2010.all.class", "pgm.sfage.2020.all.class",
                          "pgm.sfage.2010.t2.class", "pgm.sfage.2020.t2.class",
                          "pgm.sfage.2010.t10.class", "pgm.sfage.2020.t10.class",
                          "pgm.degrad.2010.forest.class", "pgm.degrad.2020.forest.class",
                          "UPF2010", "UPF.avoiddegrad", "UPF.avoiddefor")]) #keeping only raster stack
gc()
#

# scenario 2020
UPF2020<-sum(pgm.lulc.2020.forest.class, pgm.sfage.2020.all.class, pgm.degrad.2020.forest.class, na.rm = T)
UPF2020[UPF2020>1]<-0
##cheking
#UPF2020
#plot(UPF2020)

# mean upf cover in pixel scale (150m)
UPF2020.px <- focal(UPF2020, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#UPF2020.px
#anyNA(UPF2020.px[])
#plot(UPF2020.px)

names(UPF2020.px)<-"UPFpx"
UPF2020.px[is.nan(UPF2020.px)] <- 0
UPF2020.px <- mask(UPF2020.px, pgm.shp)

#saving
writeRaster(UPF2020.px, "rasters/PGM/2020_real/UPFpx.tif", format="GTiff", overwrite=T)
writeRaster(UPF2020.px, "rasters/PGM/2020_restor_wo_avoid/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (1050m)
UPF2020.ls <- focal(UPF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#UPF2020.ls
#anyNA(UPF2020.ls[])
#plot(UPF2020.ls)

names(UPF2020.ls)<-"UPFls"
UPF2020.ls[is.nan(UPF2020.ls)] <- 0
UPF2020.ls <- mask(UPF2020.ls, pgm.shp)

#saving
writeRaster(UPF2020.ls, "rasters/PGM/2020_real/UPFls.tif", format="GTiff", overwrite=T)
writeRaster(UPF2020.ls, "rasters/PGM/2020_restor_wo_avoid/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[!ls() %in% c("pgm.shp", "pgm.2010", "pgm.2020",
                          "pgm.lulc.2010.forest.class", "pgm.lulc.2020.forest.class",
                          "pgm.sfage.2010.all.class", "pgm.sfage.2020.all.class",
                          "pgm.sfage.2010.t2.class", "pgm.sfage.2020.t2.class",
                          "pgm.sfage.2010.t10.class", "pgm.sfage.2020.t10.class",
                          "pgm.degrad.2010.forest.class", "pgm.degrad.2020.forest.class",
                          "UPF2010", "UPF.avoiddegrad", "UPF.avoiddefor", "UPF2020")]) #keeping only raster stack
gc()
#



#######################################################################################################################

# [DPF] degraded primary forest -- pixel: 5x5 (150m); and landscape: 35x35 (1050m)
# this variable includes forest pixels in LULC raster (value == 3)
# which overlaps with pixels with fire (burned) and/or pixels degraded (burned and logged / logged)

# scenario 2010
DPF2010 <- pgm.degrad.2010.forest.class
#plot(DPF2010)

# mean dpf cover in pixel scale (150m)
DPF2010.px <- focal(DPF2010, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#DPF2010.px
#anyNA(DPF2010.px[])
#plot(DPF2010.px)

names(DPF2010.px)<-"DPFpx"
DPF2010.px[is.nan(DPF2010.px)] <- 0
DPF2010.px <- mask(DPF2010.px, pgm.shp)

#saving
writeRaster(DPF2010.px, "rasters/PGM/2010_real/DPFpx.tif", format="GTiff", overwrite=T)
writeRaster(DPF2010.px, "rasters/PGM/2020_avoiddegrad/DPFpx.tif", format="GTiff", overwrite=T)
writeRaster(DPF2010.px, "rasters/PGM/2020_avoidboth/DPFpx.tif", format="GTiff", overwrite=T)
writeRaster(DPF2010.px, "rasters/PGM/2020_restor_n_avoid/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (1050m)
DPF2010.ls <- focal(DPF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#DPF2010.ls
#anyNA(DPF2010.ls[])
#plot(DPF2010.ls)

names(DPF2010.ls)<-"DPFls"
DPF2010.ls[is.nan(DPF2010.ls)] <- 0
DPF2010.ls <- mask(DPF2010.ls, pgm.shp)

#saving
writeRaster(DPF2010.ls, "rasters/PGM/2010_real/DPFls.tif", format="GTiff", overwrite=T)
writeRaster(DPF2010.ls, "rasters/PGM/2020_avoiddegrad/DPFls.tif", format="GTiff", overwrite=T)
writeRaster(DPF2010.ls, "rasters/PGM/2020_avoidboth/DPFls.tif", format="GTiff", overwrite=T)
writeRaster(DPF2010.ls, "rasters/PGM/2020_restor_n_avoid/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[!ls() %in% c("pgm.shp", "pgm.2010", "pgm.2020",
                          "pgm.lulc.2010.forest.class", "pgm.lulc.2020.forest.class",
                          "pgm.sfage.2010.all.class", "pgm.sfage.2020.all.class",
                          "pgm.sfage.2010.t2.class", "pgm.sfage.2020.t2.class",
                          "pgm.sfage.2010.t10.class", "pgm.sfage.2020.t10.class",
                          "DPF2010", "pgm.degrad.2020.forest.class",
                          "UPF2010", "UPF.avoiddegrad", "UPF.avoiddefor", "UPF2020")])

gc()
#

# scenario 2020
DPF2020 <- pgm.degrad.2020.forest.class
#plot(DPF2020)

# mean dpf cover in pixel scale (150m)
DPF2020.px <- focal(DPF2020, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#DPF2020.px
#anyNA(DPF2020.px[])
#plot(DPF2020.px)

names(DPF2020.px)<-"DPFpx"
DPF2020.px[is.nan(DPF2020.px)] <- 0
DPF2020.px <- mask(DPF2020.px, pgm.shp)

#saving
writeRaster(DPF2020.px, "rasters/PGM/2020_real/DPFpx.tif", format="GTiff", overwrite=T)
writeRaster(DPF2020.px, "rasters/PGM/2020_avoiddeforest/DPFpx.tif", format="GTiff", overwrite=T)
writeRaster(DPF2020.px, "rasters/PGM/2020_restor_wo_avoid/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (1050m)
DPF2020.ls <- focal(DPF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#DPF2020.ls
#anyNA(DPF2020.ls[])
#plot(DPF2020.ls)

names(DPF2020.ls)<-"DPFls"
DPF2020.ls[is.nan(DPF2020.ls)] <- 0
DPF2020.ls <- mask(DPF2020.ls, pgm.shp)

#saving
writeRaster(DPF2020.ls, "rasters/PGM/2020_real/DPFls.tif", format="GTiff", overwrite=T)
writeRaster(DPF2020.ls, "rasters/PGM/2020_avoiddeforest/DPFls.tif", format="GTiff", overwrite=T)
writeRaster(DPF2020.ls, "rasters/PGM/2020_restor_wo_avoid/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[!ls() %in% c("pgm.shp", "pgm.2010", "pgm.2020",
                          "pgm.lulc.2010.forest.class", "pgm.lulc.2020.forest.class",
                          "pgm.sfage.2010.all.class", "pgm.sfage.2020.all.class",
                          "pgm.sfage.2010.t2.class", "pgm.sfage.2020.t2.class",
                          "pgm.sfage.2010.t10.class", "pgm.sfage.2020.t10.class",
                          "DPF2010", "DPF2020", "UPF2010", "UPF.avoiddegrad", "UPF.avoiddefor", "UPF2020")])

gc()
#



#######################################################################################################################

# [TSD] time since degradation -- pixel: 5x5 (150m); and landscape: 35x35 (1050m)
# this variable is the mean time since a degradation event

# scenario 2010
TSD2010 <- pgm.2010[["pgm.degrad.2010real"]]
#plot(TSD2010)

# mean tsd cover in pixel scale (150m)
TSD2010.px <- focal(TSD2010, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#TSD2010.px
#anyNA(TSD2010.px[])
#plot(TSD2010.px)

names(TSD2010.px)<-"TSDpx"
TSD2010.px[is.nan(TSD2010.px)] <- 0
TSD2010.px <- mask(TSD2010.px, pgm.shp)

#saving
writeRaster(TSD2010.px, "rasters/PGM/2010_real/TSDpx.tif", format="GTiff", overwrite=T)

# mean tsd cover in landscape scale (1050m)
TSD2010.ls <- focal(TSD2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#TSD2010.ls
#anyNA(TSD2010.ls[])
#plot(TSD2010.ls)

names(TSD2010.ls)<-"TSDls"
TSD2010.ls[is.nan(TSD2010.ls)] <- 0
TSD2010.ls <- mask(TSD2010.ls, pgm.shp)

#saving
writeRaster(TSD2010.ls, "rasters/PGM/2010_real/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[!ls() %in% c("pgm.shp", "pgm.2010", "pgm.2020",
                          "pgm.lulc.2010.forest.class", "pgm.lulc.2020.forest.class",
                          "pgm.sfage.2010.all.class", "pgm.sfage.2020.all.class",
                          "pgm.sfage.2010.t2.class", "pgm.sfage.2020.t2.class",
                          "pgm.sfage.2010.t10.class", "pgm.sfage.2020.t10.class",
                          "DPF2010", "DPF2020", "UPF2010", "UPF.avoiddegrad", "UPF.avoiddefor", "UPF2020",
                          "TSD2010")])

gc()
#

# scenario 2020 -- avoid degradation
TSD2010.recovery10 <- pgm.2010[["pgm.degrad.2010real"]]
TSD2010.recovery10 <- calc(TSD2010.recovery10, fun=function(x){ifelse(is.na(x), x, x+10)})
#plot(TSD2010.recovery10)

# mean dpf cover in pixel scale (150m)
TSD2010.recovery10.px <- focal(TSD2010.recovery10, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#TSD2010.recovery10.px
#anyNA(TSD2010.recovery10.px[])
#plot(TSD2010.recovery10.px)

names(TSD2010.recovery10.px)<-"TSDpx"
TSD2010.recovery10.px[is.nan(TSD2010.recovery10.px)] <- 0
TSD2010.recovery10.px <- mask(TSD2010.recovery10.px, pgm.shp)

#saving
writeRaster(TSD2010.recovery10.px, "rasters/PGM/2020_avoiddegrad/TSDpx.tif", format="GTiff", overwrite=T)
writeRaster(TSD2010.recovery10.px, "rasters/PGM/2020_avoidboth/TSDpx.tif", format="GTiff", overwrite=T)
writeRaster(TSD2010.recovery10.px, "rasters/PGM/2020_restor_n_avoid/TSDpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (1050m)
TSD2010.recovery10.ls <- focal(TSD2010.recovery10, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#TSD2010.recovery10.ls
#anyNA(TSD2010.recovery10.ls[])
#plot(TSD2010.recovery10.ls)

names(TSD2010.recovery10.ls)<-"TSDls"
TSD2010.recovery10.ls[is.nan(TSD2010.recovery10.ls)] <- 0
TSD2010.recovery10.ls <- mask(TSD2010.recovery10.ls, pgm.shp)

#saving
writeRaster(TSD2010.recovery10.ls, "rasters/PGM/2020_avoiddegrad/TSDls.tif", format="GTiff", overwrite=T)
writeRaster(TSD2010.recovery10.ls, "rasters/PGM/2020_avoidboth/TSDls.tif", format="GTiff", overwrite=T)
writeRaster(TSD2010.recovery10.ls, "rasters/PGM/2020_restor_n_avoid/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[!ls() %in% c("pgm.shp", "pgm.2010", "pgm.2020",
                          "pgm.lulc.2010.forest.class", "pgm.lulc.2020.forest.class",
                          "pgm.sfage.2010.all.class", "pgm.sfage.2020.all.class",
                          "pgm.sfage.2010.t2.class", "pgm.sfage.2020.t2.class",
                          "pgm.sfage.2010.t10.class", "pgm.sfage.2020.t10.class",
                          "DPF2010", "DPF2020", "UPF2010", "UPF.avoiddegrad", "UPF.avoiddefor", "UPF2020",
                          "TSD2010", "TSD2010.recovery10")])

gc()
#

# scenario 2020
TSD2020 <- pgm.2020[["pgm.degrad.2020real"]]
#plot(TSD2020)

# mean tsd cover in pixel scale (150m)
TSD2020.px <- focal(TSD2020, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#TSD2020.px
#anyNA(TSD2020.px[])
#plot(TSD2020.px)

names(TSD2020.px)<-"TSDpx"
TSD2020.px[is.nan(TSD2020.px)] <- 0
TSD2020.px <- mask(TSD2020.px, pgm.shp)

#saving
writeRaster(TSD2020.px, "rasters/PGM/2020_real/TSDpx.tif", format="GTiff", overwrite=T)
writeRaster(TSD2020.px, "rasters/PGM/2020_avoiddeforest/TSDpx.tif", format="GTiff", overwrite=T)
writeRaster(TSD2020.px, "rasters/PGM/2020_restor_wo_avoid/TSDpx.tif", format="GTiff", overwrite=T)

# mean tsd cover in landscape scale (1050m)
TSD2020.ls <- focal(TSD2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#TSD2020.ls
#anyNA(TSD2020.ls[])
#plot(TSD2020.ls)

names(TSD2020.ls)<-"TSDls"
TSD2020.ls[is.nan(TSD2020.ls)] <- 0
TSD2020.ls <- mask(TSD2020.ls, pgm.shp)

#saving
writeRaster(TSD2020.ls, "rasters/PGM/2020_real/TSDls.tif", format="GTiff", overwrite=T)
writeRaster(TSD2020.ls, "rasters/PGM/2020_avoiddeforest/TSDls.tif", format="GTiff", overwrite=T)
writeRaster(TSD2020.ls, "rasters/PGM/2020_restor_wo_avoid/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[!ls() %in% c("pgm.shp", "pgm.2010", "pgm.2020",
                          "pgm.lulc.2010.forest.class", "pgm.lulc.2020.forest.class",
                          "pgm.sfage.2010.all.class", "pgm.sfage.2020.all.class",
                          "pgm.sfage.2010.t2.class", "pgm.sfage.2020.t2.class",
                          "pgm.sfage.2010.t10.class", "pgm.sfage.2020.t10.class",
                          "DPF2010", "DPF2020", "UPF2010", "UPF.avoiddegrad", "UPF.avoiddefor", "UPF2020",
                          "TSD2010", "TSD2020", "TSD2010.recovery10")])

gc()
#



#######################################################################################################################

# [SF] secondary forest -- pixel: 5x5 (150m); and landscape: 35x35 (1050m)
# this variable includes forest pixels in SFage raster
# which has less than 25 years for 2010 or less than 35 for 2020

# scenario 2010
SF2010 <- pgm.sfage.2010.all.class
#plot(SF2010)

# mean sf cover in pixel scale (150m)
SF2010.px <- focal(SF2010, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#SF2010.px
#anyNA(SF2010.px[])
#plot(SF2010.px)

names(SF2010.px)<-"SFpx"
SF2010.px[is.nan(SF2010.px)] <- 0
SF2010.px <- mask(SF2010.px, pgm.shp)

#saving
writeRaster(SF2010.px, "rasters/PGM/2010_real/SFpx.tif", format="GTiff", overwrite=T)
writeRaster(SF2010.px, "rasters/PGM/2020_avoiddeforest/SFpx.tif", format="GTiff", overwrite=T)
writeRaster(SF2010.px, "rasters/PGM/2020_avoidboth/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (1050m)
SF2010.ls <- focal(SF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SF2010.ls
#anyNA(SF2010.ls[])
#plot(SF2010.ls)

names(SF2010.ls)<-"SFls"
SF2010.ls[is.nan(SF2010.ls)] <- 0
SF2010.ls <- mask(SF2010.ls, pgm.shp)

#saving
writeRaster(SF2010.ls, "rasters/PGM/2010_real/SFls.tif", format="GTiff", overwrite=T)
writeRaster(SF2010.ls, "rasters/PGM/2020_avoiddeforest/SFls.tif", format="GTiff", overwrite=T)
writeRaster(SF2010.ls, "rasters/PGM/2020_avoidboth/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[!ls() %in% c("pgm.shp", "pgm.2010", "pgm.2020",
                          "pgm.lulc.2010.forest.class", "pgm.lulc.2020.forest.class",
                          "SF2010", "pgm.sfage.2020.all.class",
                          "pgm.sfage.2010.t2.class", "pgm.sfage.2020.t2.class",
                          "pgm.sfage.2010.t10.class", "pgm.sfage.2020.t10.class",
                          "DPF2010", "DPF2020", "UPF2010", "UPF.avoiddegrad", "UPF.avoiddefor", "UPF2020",
                          "TSD2010", "TSD2020", "TSD2010.recovery10")])

gc()
#

# scenario 2020
SF2020 <- pgm.sfage.2020.all.class
#plot(SF2020)

# mean sf cover in pixel scale (150m)
SF2020.px <- focal(SF2020, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#SF2020.px
#anyNA(SF2020.px[])
#plot(SF2020.px)

names(SF2020.px)<-"SFpx"
SF2020.px[is.nan(SF2020.px)] <- 0
SF2020.px <- mask(SF2020.px, pgm.shp)

#saving
writeRaster(SF2020.px, "rasters/PGM/2020_real/SFpx.tif", format="GTiff", overwrite=T)
writeRaster(SF2020.px, "rasters/PGM/2020_avoiddegrad/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (1050m)
SF2020.ls <- focal(SF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SF2020.ls
#anyNA(SF2020.ls[])
#plot(SF2020.ls)

names(SF2020.ls)<-"SFls"
SF2020.ls[is.nan(SF2020.ls)] <- 0
SF2020.ls <- mask(SF2020.ls, pgm.shp)

#saving
writeRaster(SF2020.ls, "rasters/PGM/2020_real/SFls.tif", format="GTiff", overwrite=T)
writeRaster(SF2020.ls, "rasters/PGM/2020_avoiddegrad/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[!ls() %in% c("pgm.shp", "pgm.2010", "pgm.2020",
                          "pgm.lulc.2010.forest.class", "pgm.lulc.2020.forest.class",
                          "SF2010", "SF2020",
                          "pgm.sfage.2010.t2.class", "pgm.sfage.2020.t2.class",
                          "pgm.sfage.2010.t10.class", "pgm.sfage.2020.t10.class",
                          "DPF2010", "DPF2020", "UPF2010", "UPF.avoiddegrad", "UPF.avoiddefor", "UPF2020",
                          "TSD2010", "TSD2020", "TSD2010.recovery10")])

gc()
#

# scenario -- restoration













#



#######################################################################################################################

# [SFage] secondary forest age -- pixel: 5x5 (150m); and landscape: 35x35 (1050m)
# this variable is the mean age of secondary forest

# scenario 2010
SFage2010 <- pgm.2010[["pgm.sfage.2010real"]]
#plot(SFage2010)

# mean sf age in pixel scale (150m)
SFage2010.px <- focal(SFage2010, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#SFage2010.px
#anyNA(SFage2010.px[])
#plot(SFage2010.px)

names(SFage2010.px)<-"SFagepx"
SFage2010.px[is.nan(SFage2010.px)] <- 0
SFage2010.px <- mask(SFage2010.px, pgm.shp)

#saving
writeRaster(SFage2010.px, "rasters/PGM/2010_real/SFagepx.tif", format="GTiff", overwrite=T)

# mean sf age in landscape scale (1050m)
SFage2010.ls <- focal(SFage2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SFage2010.ls
#anyNA(SFage2010.ls[])
#plot(SFage2010.ls)

names(SFage2010.ls)<-"SFagels"
SFage2010.ls[is.nan(SFage2010.ls)] <- 0
SFage2010.ls <- mask(SFage2010.ls, pgm.shp)

#saving
writeRaster(SFage2010.ls, "rasters/PGM/2010_real/SFagels.tif", format="GTiff", overwrite=T)

rm(list=ls()[!ls() %in% c("pgm.shp", "pgm.2010", "pgm.2020",
                          "pgm.lulc.2010.forest.class", "pgm.lulc.2020.forest.class",
                          "SF2010", "SF2020", "SFage2010",
                          "pgm.sfage.2010.t2.class", "pgm.sfage.2020.t2.class",
                          "pgm.sfage.2010.t10.class", "pgm.sfage.2020.t10.class",
                          "DPF2010", "DPF2020", "UPF2010", "UPF.avoiddegrad", "UPF.avoiddefor", "UPF2020",
                          "TSD2010", "TSD2020", "TSD2010.recovery10")])

gc()
#

# scenario 2020 -- avoid deforestation
SFage2010.recovery10 <- pgm.2010[["pgm.sfage.2010real"]]
SFage2010.recovery10 <- calc(SFage2010.recovery10, fun=function(x){ifelse(is.na(x), x, x+10)})
#plot(SFage2010.recovery10)

# mean dpf cover in pixel scale (150m)
SFage2010.recovery10.px <- focal(SFage2010.recovery10, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#SFage2010.recovery10.px
#anyNA(SFage2010.recovery10.px[])
#plot(SFage2010.recovery10.px)

names(SFage2010.recovery10.px)<-"SFagepx"
SFage2010.recovery10.px[is.nan(SFage2010.recovery10.px)] <- 0
SFage2010.recovery10.px <- mask(SFage2010.recovery10.px, pgm.shp)

#saving
writeRaster(SFage2010.recovery10.px, "rasters/PGM/2020_avoiddeforest/SFagepx.tif", format="GTiff", overwrite=T)
writeRaster(SFage2010.recovery10.px, "rasters/PGM/2020_avoidboth/SFagepx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (1050m)
SFage2010.recovery10.ls <- focal(SFage2010.recovery10, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SFage2010.recovery10.ls
#anyNA(SFage2010.recovery10.ls[])
#plot(SFage2010.recovery10.ls)

names(SFage2010.recovery10.ls)<-"SFagels"
SFage2010.recovery10.ls[is.nan(SFage2010.recovery10.ls)] <- 0
SFage2010.recovery10.ls <- mask(SFage2010.recovery10.ls, pgm.shp)

#saving
writeRaster(SFage2010.recovery10.ls, "rasters/PGM/2020_avoiddeforest/SFagels.tif", format="GTiff", overwrite=T)
writeRaster(SFage2010.recovery10.ls, "rasters/PGM/2020_avoidboth/SFagels.tif", format="GTiff", overwrite=T)

rm(list=ls()[!ls() %in% c("pgm.shp", "pgm.2010", "pgm.2020",
                          "pgm.lulc.2010.forest.class", "pgm.lulc.2020.forest.class",
                          "SF2010", "SF2020", "SFage2010", "SFage2010.recovery10",
                          "pgm.sfage.2010.t2.class", "pgm.sfage.2020.t2.class",
                          "pgm.sfage.2010.t10.class", "pgm.sfage.2020.t10.class",
                          "DPF2010", "DPF2020", "UPF2010", "UPF.avoiddegrad", "UPF.avoiddefor", "UPF2020",
                          "TSD2010", "TSD2020", "TSD2010.recovery10")])

gc()
#

# scenario 2020
SFage2020 <- pgm.2020[["pgm.sfage.2020real"]]
#plot(SFage2010)

# mean sf age in pixel scale (150m)
SFage2020.px <- focal(SFage2020, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#SFage2020.px
#anyNA(SFage2020.px[])
#plot(SFage2020.px)

names(SFage2020.px)<-"SFagepx"
SFage2020.px[is.nan(SFage2020.px)] <- 0
SFage2020.px <- mask(SFage2020.px, pgm.shp)

#saving
writeRaster(SFage2020.px, "rasters/PGM/2020_real/SFagepx.tif", format="GTiff")
writeRaster(SFage2020.px, "rasters/PGM/2020_avoiddegrad/SFagepx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (1050m)
SFage2020.ls <- focal(SFage2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SFage2020.ls
#anyNA(SFage2020.ls[])
#plot(SFage2020.ls)

names(SFage2020.ls)<-"SFagels"
SFage2020.ls[is.nan(SFage2020.ls)] <- 0
SFage2020.ls <- mask(SFage2020.ls, pgm.shp)

#saving
writeRaster(SFage2020.ls, "rasters/PGM/2020_real/SFagels.tif", format="GTiff", overwrite=T)
writeRaster(SFage2020.ls, "rasters/PGM/2020_avoiddegrad/SFagels.tif", format="GTiff", overwrite=T)

rm(list=ls()[!ls() %in% c("pgm.shp", "pgm.2010", "pgm.2020",
                          "pgm.lulc.2010.forest.class", "pgm.lulc.2020.forest.class",
                          "SF2010", "SF2020", "SFage2010", "SFage2020", "SFage2010.recovery10",
                          "pgm.sfage.2010.t2.class", "pgm.sfage.2020.t2.class",
                          "pgm.sfage.2010.t10.class", "pgm.sfage.2020.t10.class",
                          "DPF2010", "DPF2020", "UPF2010", "UPF.avoiddegrad", "UPF.avoiddefor", "UPF2020",
                          "TSD2010", "TSD2020", "TSD2010.recovery10")])

gc()
#

# scenario -- restoration













#



#######################################################################################################################

# [F3] Forest type 3 or Total forest -- pixel: 5x5 (150m)
# this variable includes forest pixels in LULC raster
# including degraded forest and secondary forest older than 2 years

# paragominas 2010
#names(pgm.2010)

# isolating forest class pixels
pgm.2010.forest.class <- pgm.2010[["PGM.LULC.2010"]]
pgm.2010.forest.class[pgm.2010.forest.class==3] <- 1
pgm.2010.forest.class[pgm.2010.forest.class>1] <- 0
##cheking
#pgm.2010.forest.class
#plot(pgm.2010.forest.class)

# isolating SF <= 2 years old
pgm.2010.SFyoung <- pgm.2010[["PGM.SFage.2010"]]
pgm.2010.SFyoung[pgm.2010.SFyoung<=2] <- 1
pgm.2010.SFyoung[pgm.2010.SFyoung>2] <- 0
##cheking
#pgm.2010.SFyoung
#plot(pgm.2010.SFyoung)

# excluding areas with SF <= 2
TF2010<-sum(pgm.2010.forest.class, pgm.2010.SFyoung, na.rm = T)
TF2010[TF2010>1]<-0
##cheking
#TF2010
#plot(TF2010)

## including areas with fire and degraded
#TF2010<-sum(TF2010, pgm.2010[["PGM.Fire.2010"]], pgm.2010[["PGM.Degrad.2010"]], na.rm = T)
#TF2010[TF2010>=1]<-1
###cheking
##TF2010
##plot(TF2010)

# mean upf cover in pixel scale (150m)
TF2010.px <- focal(TF2010, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#TF2010.px
#anyNA(TF2010.px[])
#plot(TF2010.px)

names(TF2010.px)<-"TFpx"
TF2010.px[is.nan(TF2010.px)] <- 0

#saving
writeRaster(TF2010.px, "rasters/PGM/2010/TFpx.tif", format="GTiff")

rm(list=ls()[!ls() %in% c("pgm.2010", "pgm.2020")]) #keeping only raster stack
gc()
#

# paragominas 2020
#names(pgm.2020)

# isolating forest class pixels
pgm.2020.forest.class <- pgm.2020[["PGM.LULC.2020"]]
pgm.2020.forest.class[pgm.2020.forest.class==3] <- 1
pgm.2020.forest.class[pgm.2020.forest.class>1] <- 0
##cheking
#pgm.2020.forest.class
#plot(pgm.2020.forest.class)


# isolating SF <= 2 years old
pgm.2020.SFyoung <- pgm.2020[["PGM.SFage.2020"]]
pgm.2020.SFyoung[pgm.2020.SFyoung<=2] <- 1
pgm.2020.SFyoung[pgm.2020.SFyoung>2] <- 0
##cheking
#pgm.2020.SFyoung
#plot(pgm.2020.SFyoung)

# excluding areas with SF <= 2
TF2020<-sum(pgm.2020.forest.class, pgm.2020.SFyoung, na.rm = T)
TF2020[TF2020>1]<-0
##cheking
#TF2020
#plot(TF2020)

## including areas with fire and degraded
#TF2020<-sum(TF2020, pgm.2020[["PGM.Fire.2020"]], pgm.2020[["PGM.Degrad.2020"]], na.rm = T)
#TF2020[TF2020>=1]<-1
###cheking
##TF2020
##plot(TF2020)

# mean upf cover in pixel scale (150m)
TF2020.px <- focal(TF2020, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#TF2020.px
#anyNA(TF2020.px[])
#plot(TF2020.px)

names(TF2020.px)<-"TFpx"
TF2020.px[is.nan(TF2020.px)] <- 0

#saving
writeRaster(TF2020.px, "rasters/PGM/2020/TFpx.tif", format="GTiff")

rm(list=ls()[!ls() %in% c("pgm.2010", "pgm.2020")]) #keeping only raster stack
gc()
#



#######################################################################################################################

# [F1] Forest type 1 or Mature forest -- pixel: 5x5 (150m)
# this variable includes forest pixels in LULC raster
# including degraded forest and secondary forest older than 10 years

# paragominas 2010
#names(pgm.2010)

# isolating forest class pixels
pgm.2010.forest.class <- pgm.2010[["PGM.LULC.2010"]]
pgm.2010.forest.class[pgm.2010.forest.class==3] <- 1
pgm.2010.forest.class[pgm.2010.forest.class>1] <- 0
##cheking
#pgm.2010.forest.class
#plot(pgm.2010.forest.class)

# isolating SF <= 10 years old
pgm.2010.SFold <- pgm.2010[["PGM.SFage.2010"]]
pgm.2010.SFold[pgm.2010.SFold<=10] <- 1
pgm.2010.SFold[pgm.2010.SFold>10] <- 0
##cheking
#pgm.2010.SFold
#plot(pgm.2010.SFold)

# excluding areas with SF <= 10
MF2010<-sum(pgm.2010.forest.class, pgm.2010.SFold, na.rm = T)
MF2010[MF2010>1]<-0
##cheking
#MF2010
#plot(MF2010)

## including areas with fire and degraded
#MF2010<-sum(MF2010, pgm.2010[["PGM.Fire.2010"]], pgm.2010[["PGM.Degrad.2010"]], na.rm = T)
#MF2010[MF2010>=1]<-1
###cheking
##MF2010
##plot(MF2010)

# mean mature forest cover in pixel scale (150m)
MF2010.px <- focal(MF2010, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#MF2010.px
#anyNA(MF2010.px[])
#plot(MF2010.px)

names(MF2010.px)<-"MFpx"
MF2010.px[is.nan(MF2010.px)] <- 0

#saving
writeRaster(MF2010.px, "rasters/PGM/2010/MFpx.tif", format="GTiff")

rm(list=ls()[!ls() %in% c("pgm.2010", "pgm.2020")]) #keeping only raster stack
gc()
#

# paragominas 2020
#names(pgm.2020)

# isolating forest class pixels
pgm.2020.forest.class <- pgm.2020[["PGM.LULC.2020"]]
pgm.2020.forest.class[pgm.2020.forest.class==3] <- 1
pgm.2020.forest.class[pgm.2020.forest.class>1] <- 0
##cheking
#pgm.2020.forest.class
#plot(pgm.2020.forest.class)

# isolating SF <= 10 years old
pgm.2020.SFold <- pgm.2020[["PGM.SFage.2020"]]
pgm.2020.SFold[pgm.2020.SFold<=10] <- 1
pgm.2020.SFold[pgm.2020.SFold>10] <- 0
##cheking
#pgm.2020.SFold
#plot(pgm.2020.SFold)

# excluding areas with SF <= 10
MF2020<-sum(pgm.2020.forest.class, pgm.2020.SFold, na.rm = T)
MF2020[MF2020>1]<-0
##cheking
#MF2020
#plot(MF2020)

## including areas with fire and degraded
#MF2020<-sum(MF2020, pgm.2020[["PGM.Fire.2020"]], pgm.2020[["PGM.Degrad.2020"]], na.rm = T)
#MF2020[MF2020>=1]<-1
###cheking
##MF2020
##plot(MF2020)

# mean mature forest cover in pixel scale (150m)
MF2020.px <- focal(MF2020, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#MF2020.px
#anyNA(MF2020.px[])
#plot(MF2020.px)

names(MF2020.px)<-"MFpx"
MF2020.px[is.nan(MF2020.px)] <- 0

#saving
writeRaster(MF2020.px, "rasters/PGM/2020/MFpx.tif", format="GTiff")

rm(list=ls()[!ls() %in% c("pgm.2010", "pgm.2020")]) #keeping only raster stack
gc()
#



##############################################################################################################################################################################################################################################

# [edgedist] distance to forest edge
# this variable is the mean distance of F1 to the edge

# paragominas 2010
#names(pgm.2010)

# isolating forest class pixels
pgm.2010.forest.class <- pgm.2010[["PGM.LULC.2010"]]
pgm.2010.forest.class[pgm.2010.forest.class==3] <- 1
pgm.2010.forest.class[pgm.2010.forest.class>1] <- 0
##cheking
#pgm.2010.forest.class
#plot(pgm.2010.forest.class)

# isolating SF <= 10 years old
pgm.2010.SFold <- pgm.2010[["PGM.SFage.2010"]]
pgm.2010.SFold[pgm.2010.SFold<=10] <- 1
pgm.2010.SFold[pgm.2010.SFold>10] <- 0
##cheking
#pgm.2010.SFold
#plot(pgm.2010.SFold)

# excluding areas with SF <= 10
MF2010<-sum(pgm.2010.forest.class, pgm.2010.SFold, na.rm = T)
MF2010[MF2010>1]<-0
##cheking
#MF2010
#plot(MF2010)

## including areas with fire and degraded
#MF2010<-sum(MF2010, pgm.2010[["PGM.Fire.2010"]], pgm.2010[["PGM.Degrad.2010"]], na.rm = T)
#MF2010[MF2010>=1]<-1
##MF2010[MF2010==0]<-NA
###cheking
##MF2010
##plot(MF2010)


# calculating distance, for all cells that are mature forest
# to the nearest cell that is not mf
pts2010 <- rasterToPoints(MF2010, spatial=TRUE)
core2010 <- pts2010[pts2010$layer == "0",]
#cheking
#core2010
#plot(core2010, add=T, col="red")

rm(pgm.2010.forest.class); rm(pgm.2010.SFold); gc() #keeping only raster in use

edge.dist.2010 <- rasterDistance(pts2010, core2010, reference = MF2010, scale=TRUE)
##cheking
#edge.dist.2010
#anyNA(edge.dist.2010[])
#plot(edge.dist.2010)

names(edge.dist.2010)<-"edgedist"

#saving
writeRaster(edge.dist.2010, "rasters/PGM/2010/edgedist.tif", format="GTiff")

rm(pts2010); rm(core2010); rm(edge.dist.2010); gc() #keeping only raster in use

# paragominas 2020
#names(pgm.2020)

# isolating forest class pixels
pgm.2020.forest.class <- pgm.2020[["PGM.LULC.2020"]]
pgm.2020.forest.class[pgm.2020.forest.class==3] <- 1
pgm.2020.forest.class[pgm.2020.forest.class>1] <- 0
##cheking
#pgm.2020.forest.class
#plot(pgm.2020.forest.class)

# isolating SF <= 10 years old
pgm.2020.SFold <- pgm.2020[["PGM.SFage.2020"]]
pgm.2020.SFold[pgm.2020.SFold<=10] <- 1
pgm.2020.SFold[pgm.2020.SFold>10] <- 0
##cheking
#pgm.2020.SFold
#plot(pgm.2020.SFold)

# excluding areas with SF <= 10
MF2020<-sum(pgm.2020.forest.class, pgm.2020.SFold, na.rm = T)
MF2020[MF2020>1]<-0
##cheking
#MF2020
#plot(MF2020)

## including areas with fire and degraded
#MF2020<-sum(MF2020, pgm.2020[["PGM.Fire.2020"]], pgm.2020[["PGM.Degrad.2020"]], na.rm = T)
#MF2020[MF2020>=1]<-1
##MF2020[MF2020==0]<-NA
###cheking
##MF2020
##plot(MF2020)


# calculating distance, for all cells that are mature forest
# to the nearest cell that is not mf
pts2020 <- rasterToPoints(MF2020, spatial=TRUE)
core2020 <- pts2020[pts2020$layer == "0",]
#cheking
#core2020
#plot(core2020, add=T, col="red")

rm(pgm.2020.forest.class); rm(pgm.2020.SFold); gc() #keeping only raster in use

edge.dist.2020 <- rasterDistance(pts2020, core2020, reference = MF2020, scale=TRUE)
##cheking
#edge.dist.2020
#anyNA(edge.dist.2020[])
#plot(edge.dist.2020)

names(edge.dist.2020)<-"edgedist"

#saving
writeRaster(edge.dist.2020, "rasters/PGM/2020/edgedist.tif", format="GTiff")

rm(pts2020); rm(core2020); rm(edge.dist.2020); gc() #keeping only raster in use

#

#

# [edge] forest edge -- pixel: 5x5 (150m); and landscape: 35x35 (1050m)
# this variable is the mean of edge area based on mature forest

# paragominas 2010

# marking edge and core areas
edge2010 <- focal(MF2010, matrix(1,ncol=3,nrow=3), fun = mean, na.rm = T, pad = T)
edge2010 <- edge2010 * MF2010
edge2010[edge2010 == 1] <- 0                  # core area
edge2010[edge2010 > 0 & edge2010 < 2] <- 1    # edge
#cheking
#edge2010
#unique(edge2010[])
#plot(edge2010)

# mean edge in pixel scale (150m)
edge2010.px <- focal(edge2010, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#edge2010.px
#anyNA(edge2010.px[])
#plot(edge2010.px)

names(edge2010.px)<-"edgepx"
edge2010.px[is.nan(edge2010.px)] <- 0

#saving
writeRaster(edge2010.px, "rasters/PGM/2010/edgepx.tif", format="GTiff")
#

# mean sf cover in landscape scale (1050m)
edge2010.ls <- focal(edge2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#edge2010.ls
#anyNA(edge2010.ls[])
#plot(edge2010.ls)

names(edge2010.ls)<-"edgels"
edge2010.ls[is.nan(edge2010.ls)] <- 0

#saving
writeRaster(edge2010.ls, "rasters/PGM/2010/edgels.tif", format="GTiff")
#

# paragominas 2020

# marking edge and core areas
edge2020 <- focal(MF2020, matrix(1,ncol=3,nrow=3), fun = mean, na.rm = T, pad = T)
edge2020 <- edge2020 * MF2020
edge2020[edge2020 == 1] <- 0                  # core area
edge2020[edge2020 > 0 & edge2020 < 2] <- 1    # edge
#cheking
#edge2020
#unique(edge2020[])
#plot(edge2020)

# mean edge in pixel scale (150m)
edge2020.px <- focal(edge2020, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#edge2020.px
#anyNA(edge2020.px[])
#plot(edge2020.px)

names(edge2020.px)<-"edgepx"
edge2020.px[is.nan(edge2020.px)] <- 0

#saving
writeRaster(edge2020.px, "rasters/PGM/2020/edgepx.tif", format="GTiff")
#

# mean sf cover in landscape scale (1050m)
edge2020.ls <- focal(edge2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#edge2020.ls
#anyNA(edge2020.ls[])
#plot(edge2020.ls)

names(edge2020.ls)<-"edgels"
edge2020.ls[is.nan(edge2020.ls)] <- 0

#saving
writeRaster(edge2020.ls, "rasters/PGM/2020/edgels.tif", format="GTiff")

rm(list=ls()[!ls() %in% c("pgm.2010", "pgm.2020")]) #keeping only raster stack
gc()
#

#

# [waterdist] distance to water
# this variable is the euclidean distance between pixels to water body

pgm.water.shp <- readOGR(dsn = "rasters/PGM/input", layer = "pgm-water") #reading shapefile
pgm.water.pts <- as(pgm.water.shp, "SpatialPointsDataFrame")
#cheking
#pgm.water.pts
#plot(pgm.water.pts)


pts2010 <- rasterToPoints(pgm.2020[["PGM.LULC.2020"]], spatial=TRUE)

water.dist <- rasterDistance(pts2010, pgm.water.pts, reference = pgm.2020[["PGM.LULC.2020"]], scale=TRUE)
##cheking
#water.dist
#anyNA(water.dist[])
#plot(water.dist)

names(water.dist)<-"waterdist"

#saving
writeRaster(water.dist, "rasters/PGM/2010/waterdist.tif", format="GTiff")
writeRaster(water.dist, "rasters/PGM/2020/waterdist.tif", format="GTiff")

rm(list=ls()[!ls() %in% c("pgm.2010", "pgm.2020")]) #keeping only raster stack
gc()
#

#

## [roaddist] distance to road
## this variable is the euclidean distance between pixels to road
#
#pgm.road.shp <- readOGR(dsn = "rasters/PGM/input", layer = "pgm-road") #reading shapefile
#pgm.road.pts <- as(pgm.road.shp, "SpatialPointsDataFrame")
##cheking
##pgm.road.pts
##plot(pgm.road.pts)
#
#
#pts2010 <- rasterToPoints(pgm.2020[["PGM.LULC.2020"]], spatial=TRUE)
#
#road.dist <- rasterDistance(pts2010, pgm.road.pts, reference = pgm.2020[["PGM.LULC.2020"]], scale=TRUE)
###cheking
##road.dist
##anyNA(road.dist[])
##plot(road.dist)
#
#names(road.dist)<-"roaddist"
#
##saving
#writeRaster(road.dist, "rasters/PGM/2010/roaddist.tif", format="GTiff")
#writeRaster(road.dist, "rasters/PGM/2020/roaddist.tif", format="GTiff")
#
#rm(list=ls()[!ls() %in% c("pgm.2010", "pgm.2020")]) #keeping only raster stack
#gc()
##
#
##

# [temp] historical mean annual temperature from worldclim -- 1970:2000
mean.temp <- raster("C:/Users/miral/Dropbox/GIS/wc2.1_30s_bio/wc2.1_30s_bio_1.tif")

# Conversion of rasters into same extent
pgm.mean.temp <- resample(mean.temp, pgm.2010[["PGM.LULC.2010"]], method='bilinear')
##cheking
#pgm.mean.temp
#anyNA(pgm.mean.temp[])
#plot(pgm.mean.temp)

names(pgm.mean.temp)<-"meantemp"

values(pgm.mean.temp) <- rescale(values(pgm.mean.temp)) # scale values from min=0 to max=1
pgm.mean.temp.clipp <- mask(pgm.mean.temp, pgm.2010[["PGM.LULC.2010"]]) # clipping to keep only municipilaty 

#saving
writeRaster(pgm.mean.temp.clipp, "rasters/PGM/2010/meantemp.tif", format="GTiff")
writeRaster(pgm.mean.temp.clipp, "rasters/PGM/2020/meantemp.tif", format="GTiff")
#

#

# [precip] historical mean annual precipitation from worldclim -- 1970:2000
mean.precip <- raster("C:/Users/miral/Dropbox/GIS/wc2.1_30s_bio/wc2.1_30s_bio_12.tif")

# Conversion of rasters into same extent
pgm.mean.precip <- resample(mean.precip, pgm.2010[["PGM.LULC.2010"]], method='bilinear')
##cheking
#pgm.mean.precip
#anyNA(pgm.mean.precip[])
#plot(pgm.mean.precip)

names(pgm.mean.precip)<-"meanprecip"

values(pgm.mean.precip) <- rescale(values(pgm.mean.precip)) # scale values from min=0 to max=1
pgm.mean.precip.clipp <- mask(pgm.mean.precip, pgm.2010[["PGM.LULC.2010"]]) # clipping to keep only municipilaty 

#saving
writeRaster(pgm.mean.precip.clipp, "rasters/PGM/2010/meanprecip.tif", format="GTiff")
writeRaster(pgm.mean.precip.clipp, "rasters/PGM/2020/meanprecip.tif", format="GTiff")
#

#

# [elev] elevation data derived from the SRTM
elev <- raster("C:/Users/miral/Dropbox/GIS/wc2.1_30s_bio/wc2.1_30s_elev.tif")

# Conversion of rasters into same extent
pgm.elev <- resample(elev, pgm.2010[["PGM.LULC.2010"]], method='bilinear')
##cheking
#pgm.elev
#anyNA(pgm.elev[])
#plot(pgm.elev)

names(pgm.elev)<-"elev"

values(pgm.elev) <- rescale(values(pgm.elev)) # scale values from min=0 to max=1
pgm.elev.clipp <- mask(pgm.elev, pgm.2010[["PGM.LULC.2010"]]) # clipping to keep only municipilaty 

#saving
writeRaster(pgm.elev.clipp, "rasters/PGM/2010/elev.tif", format="GTiff")
writeRaster(pgm.elev.clipp, "rasters/PGM/2020/elev.tif", format="GTiff")
#

#


rm(list=ls()) #end
gc()
#

#


#
