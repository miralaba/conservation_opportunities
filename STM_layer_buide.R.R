
#' @title Cost-effectiveness of conservation actions in Amazon
#' @description scrip to build exploratory variables from 
#' land use - land cover, secondary forest, edge, 
#' degradation (logguin and fire), temperature, precipitation, 
#' elevation and distances to road and water body
#' #' in Santarém - Belterra - Mojuí dos Campos region - PA
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


#### importing input rasters ####
#################################
# skip this if the explanatory  #
# variables to sdm have already #
# been created (see line xxx)   #
#################################

# land use land cover from mapbiomas collection 7 [2010 and 2020]
## [PGM] paragominas
#pgm.lulc <- stack(c("rasters/PGM/input/pgm-2010-lulc-mapbiomas-brazil-collection-70.tif",
#                    "rasters/PGM/input/pgm-2020-lulc-mapbiomas-brazil-collection-70.tif"))
#names(pgm.lulc) <- c("PGM.LULC.2010", "PGM.LULC.2020")
##checking
##pgm.lulc
##plot(pgm.lulc)
##sort(unique(values(pgm.lulc[["PGM.LULC.2010"]])))
#
#values(pgm.lulc)[values(pgm.lulc) <= 0] = NA
## land ues land cover pixel values and codes
## 0  == NA
## 3  == Forest Formation      == Forest
## 4  == Savanna Formation     == Forest
## 9  == Forest Plantation     == Farming
## 11 == Wetland               == Non Forest Natural Formation
## 12 == Grassland             == Non Forest Natural Formation
## 15 == Pasture               == Farming
## 24 == Urban Area            == Non vegetated area
## 30 == Mining                == Non vegetated area
## 33 == River, Lake and Ocean == Water
## 39 == Soybean               == Farming
## 41 == Other Temporary Crops == Farming
## 48 == Other Perennial Crops == Farming

# [STM] santarem
stm.lulc <- stack(c("rasters/STM/input/stm-2010-lulc-mapbiomas-brazil-collection-70.tif",
                    "rasters/STM/input/stm-2020-lulc-mapbiomas-brazil-collection-70.tif"))
names(stm.lulc) <- c("STM.LULC.2010", "STM.LULC.2020")
#checking
#stm.lulc
#plot(stm.lulc)
#sort(unique(values(stm.lulc[["STM.LULC.2010"]])))

values(stm.lulc)[values(stm.lulc) <= 0] = NA
# land ues land cover pixel values and codes
# 0  == NA
# 3  == Forest Formation      == Forest
# 4  == Savanna Formation     == Forest
# 11 == Wetland               == Non Forest Natural Formation
# 12 == Grassland             == Non Forest Natural Formation
# 15 == Pasture               == Farming
# 24 == Urban Area            == Non vegetated area
# 33 == River, Lake and Ocean == Water
# 39 == Soybean               == Farming
# 41 == Other Temporary Crops == Farming

# secondary forest age from Silva Jr. et al 2020  [2010 and 2020]
# [DOI: 10.1038/s41597-020-00600-4]
## [PGM] paragominas
#pgm.sfage <- stack(c("rasters/PGM/input/pgm-2010-sfage-mapbiomas-brazil-collection-60.tif",
#                      "rasters/PGM/input/pgm-2020-sfage-mapbiomas-brazil-collection-60.tif"))
#names(pgm.sfage) <- c("PGM.SFage.2010", "PGM.SFage.2020")
##checking
##pgm.sfage
##plot(pgm.sfage)
##range(values(pgm.sfage[["PGM.SFage.2010"]]))
#
#values(pgm.sfage)[values(pgm.sfage) <= 0] = NA

# Conversion of rasters into same extent
pgm.sfage_resampled <- resample(pgm.sfage, pgm.lulc, method='ngb')

# [STM] santarem
stm.sfage <- stack(c("rasters/STM/input/stm-2010-sfage-mapbiomas-brazil-collection-60.tif",
                     "rasters/STM/input/stm-2020-sfage-mapbiomas-brazil-collection-60.tif"))
names(stm.sfage) <- c("STM.SFage.2010", "STM.SFage.2020")
#checking
#stm.sfage
#plot(stm.sfage)
#range(values(stm.sfage[["STM.SFage.2010"]]))

values(stm.sfage)[values(stm.sfage) <= 0] = NA

# Conversion of rasters into same extent
stm.sfage_resampled <- resample(stm.sfage, stm.lulc, method='ngb')

# fire from mapbiomas fogo collection 1 [2010 and 2020]
## [PGM] paragominas
#pgm.fire <- stack(c("rasters/PGM/input/pgm-2010-fire-mapbiomas-brazil-collection-10.tif",
#                    "rasters/PGM/input/pgm-2020-fire-mapbiomas-brazil-collection-10.tif"))
#names(pgm.fire) <- c("PGM.Fire.2010", "PGM.Fire.2020")
##checking
##pgm.fire
##plot(pgm.fire)
##range(values(pgm.fire[["PGM.Fire.2010"]]))
#
#values(pgm.fire)[values(pgm.fire) <= 0] = NA
#
## Conversion of rasters into same extent
#pgm.fire_resampled <- resample(pgm.fire, pgm.lulc, method='bilinear')

# [STM] santarem
stm.fire <- stack(c("rasters/STM/input/stm-2010-fire-mapbiomas-brazil-collection-10.tif",
                    "rasters/STM/input/stm-2020-fire-mapbiomas-brazil-collection-10.tif"))
names(stm.fire) <- c("STM.Fire.2010", "STM.Fire.2020")
#checking
#stm.fire
#plot(stm.fire)
#range(values(stm.fire[["STM.Fire.2010"]]))

values(stm.fire)[values(stm.fire) <= 0] = NA

# Conversion of rasters into same extent
stm.fire_resampled <- resample(stm.fire, stm.lulc, method='ngb')

# degradation data from INPE 
# [load_degrad() function in datazoom.amazonia for 2010 and 
# deter for 2020]
## [PGM] paragominas
#pgm.degrad.2010.shp <- readOGR(dsn = "rasters/PGM/input", layer = "pgm-2010-degrad-inpe") #reading shapefile
#pgm.degrad.2010 <- rasterize(pgm.degrad.2010.shp, pgm.lulc[[1]], field=1) #rasterizing
#
#pgm.degrad.2020.shp <- readOGR(dsn = "rasters/PGM/input", layer = "pgm-2020-degrad-deter-inpe") #reading shapefile
#pgm.degrad.2020 <- rasterize(pgm.degrad.2020.shp, pgm.lulc[[1]], field=1) #rasterizing
#
#pgm.degrad <- stack(pgm.degrad.2010,pgm.degrad.2020) # stack
#names(pgm.degrad) <- c("PGM.Degrad.2010", "PGM.Degrad.2020")
#
#rm(list=ls()[ls() %in% c("pgm.degrad.2010.shp", "pgm.degrad.2010", "pgm.degrad.2020.shp", "pgm.degrad.2020")]) #keeping only raster stack
#
##checking
##pgm.degrad
##plot(pgm.degrad)
##unique(values(pgm.degrad[["PGM.Degrad.2010"]]))
#
## Conversion of rasters into same extent
#pgm.degrad_resampled <- resample(pgm.degrad, pgm.lulc, method='bilinear')

# [STM] santarem
stm.degrad.2010.shp <- readOGR(dsn = "rasters/STM/input", layer = "stm-2010-degrad-inpe") #reading shapefile
stm.degrad.2010 <- rasterize(stm.degrad.2010.shp, stm.lulc[[1]], field=1) #rasterizing

stm.degrad.2020.shp <- readOGR(dsn = "rasters/STM/input", layer = "stm-2020-degrad-deter-inpe") #reading shapefile
stm.degrad.2020 <- rasterize(stm.degrad.2020.shp, stm.lulc[[1]], field=1) #rasterizing

stm.degrad <- stack(stm.degrad.2010,stm.degrad.2020) # stack
names(stm.degrad) <- c("STM.Degrad.2010", "STM.Degrad.2020")

rm(list=ls()[ls() %in% c("stm.degrad.2010.shp", "stm.degrad.2010", "stm.degrad.2020.shp", "stm.degrad.2020")]) #keeping only raster stack

#checking
#stm.degrad
#plot(stm.degrad)
#unique(values(stm.degrad[["STM.Degrad.2010"]]))

# Conversion of rasters into same extent
stm.degrad_resampled <- resample(stm.degrad, stm.lulc, method='ngb')

## roads [???]
### [PGM] paragominas
##pgm.road <- readOGR()
##
###checking
##pgm.road
##plot(pgm.road)
#
## [STM] santarem
#stm.road <- readOGR()
#
##checking
#stm.road
#plot(stm.road)


# function to make sure all raster are fully stacked
#intersect_mask <- function(x){
#  values_x <- getValues(x)
#  inter_x <- values_x %*% rep(1,nlayers(x))
#  mask <- setValues(subset(x,1),values = (inter_x>0))
#  return(mask)
#}

# stacking by location by year
#pgm.2010 <- stack(pgm.lulc[[1]], pgm.sfage_resampled[[1]], pgm.fire_resampled[[1]], pgm.degrad_resampled[[1]])
##pgm.2010 <- stack(mask(pgm.2010, intersect_mask(pgm.2010)))
#
#pgm.2020 <- stack(pgm.lulc[[2]], pgm.sfage_resampled[[2]], pgm.fire_resampled[[2]], pgm.degrad_resampled[[2]])
##pgm.2020 <- stack(mask(pgm.2020, intersect_mask(pgm.2020)))

stm.2010 <- stack(stm.lulc[[1]], stm.sfage_resampled[[1]], stm.fire_resampled[[1]], stm.degrad_resampled[[1]])
#stm.2010 <- stack(mask(stm.2010, intersect_mask(stm.2010)))

stm.2020 <- stack(stm.lulc[[2]], stm.sfage_resampled[[2]], stm.fire_resampled[[2]], stm.degrad_resampled[[2]])
#stm.2020 <- stack(mask(stm.2020, intersect_mask(stm.2020)))

rm(list=ls()[!ls() %in% c("stm.2010", "stm.2020")]) #keeping only raster stack
gc()
#

#### setting explanatory variables ####
# [UPF] undisturbed primary forest -- pixel: 5x5 (150m); and landscape: 35x35 (1050m)
# this variable includes all forest pixels in LULC raster (value == 3)
# excluding those with age < 25 in 2010 SF raster or age < 35 in 2020 SF raster
# excluding pixels with fire
# excluding pixels degraded

# santarem 2010
#names(stm.2010)

# isolating forest class pixels
stm.2010.forest.class <- stm.2010[["STM.LULC.2010"]]
stm.2010.forest.class[stm.2010.forest.class==3] <- 1
stm.2010.forest.class[stm.2010.forest.class>1] <- 0
##cheking
#stm.2010.forest.class
#plot(stm.2010.forest.class)

# isolating SF < 25 years old
stm.2010.SFless25 <- stm.2010[["STM.SFage.2010"]]
stm.2010.SFless25[stm.2010.SFless25<25] <- 1
stm.2010.SFless25[stm.2010.SFless25==25] <- 0
##cheking
#stm.2010.SFless25
#plot(stm.2010.SFless25)

# excluding areas with SF > 25, with fire and degraded
UPF2010<-sum(stm.2010.forest.class, stm.2010.SFless25, stm.2010[["STM.Fire.2010"]], stm.2010[["STM.Degrad.2010"]], na.rm = T)
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

#saving
writeRaster(UPF2010.px, "rasters/STM/2010/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (1050m)
UPF2010.ls <- focal(UPF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#UPF2010.ls
#anyNA(UPF2010.ls[])
#plot(UPF2010.ls)

names(UPF2010.ls)<-"UPFls"
UPF2010.ls[is.nan(UPF2010.ls)] <- 0

#saving
writeRaster(UPF2010.ls, "rasters/STM/2010/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[!ls() %in% c("stm.2010", "stm.2020")]) #keeping only raster stack
gc()
#

# santarem 2020
#names(stm.2020)

# isolating forest class pixels
stm.2020.forest.class <- stm.2020[["STM.LULC.2020"]]
stm.2020.forest.class[stm.2020.forest.class==3] <- 1
stm.2020.forest.class[stm.2020.forest.class>1] <- 0
##cheking
#stm.2020.forest.class
#plot(stm.2020.forest.class)

# isolating SF < 35 years old
stm.2020.SFless25 <- stm.2020[["STM.SFage.2020"]]
stm.2020.SFless25[stm.2020.SFless25<35] <- 1
stm.2020.SFless25[stm.2020.SFless25>=35] <- 0
##cheking
#stm.2020.SFless25
#plot(stm.2020.SFless25)

# excluding areas with SF > 35, with fire and degraded
UPF2020<-sum(stm.2020.forest.class, stm.2020.SFless25, stm.2020[["STM.Fire.2020"]], stm.2020[["STM.Degrad.2020"]], na.rm = T)
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

#saving
writeRaster(UPF2020.px, "rasters/STM/2020/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (1050m)
UPF2020.ls <- focal(UPF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#UPF2020.ls
#anyNA(UPF2020.ls[])
#plot(UPF2020.ls)

names(UPF2020.ls)<-"UPFls"
UPF2020.ls[is.nan(UPF2020.ls)] <- 0

#saving
writeRaster(UPF2020.ls, "rasters/STM/2020/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[!ls() %in% c("stm.2010", "stm.2020")]) #keeping only raster stack
gc()
#

#

# [DPF] degraded primary forest -- pixel: 5x5 (150m); and landscape: 35x35 (1050m)
# this variable includes forest pixels in LULC raster (value == 3)
# which overlaps with pixels with fire (burned) and/or pixels degraded (burned and logged / logged)

# santarem 2010
#names(stm.2010)

# isolating forest class pixels
stm.2010.forest.class <- stm.2010[["STM.LULC.2010"]]
stm.2010.forest.class[stm.2010.forest.class==3] <- 1
stm.2010.forest.class[stm.2010.forest.class>1] <- 0
##cheking
#stm.2010.forest.class
#plot(stm.2010.forest.class)

# forested areas with fire and/or degraded
DPF2010<-sum(stm.2010.forest.class, stm.2010[["STM.Fire.2010"]], stm.2010[["STM.Degrad.2010"]], na.rm = T)
DPF2010[DPF2010<=1]<-0
DPF2010[DPF2010>1]<-1
##cheking
#DPF2010
#plot(DPF2010)

# mean dpf cover in pixel scale (150m)
DPF2010.px <- focal(DPF2010, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#DPF2010.px
#anyNA(DPF2010.px[])
#plot(DPF2010.px)

names(DPF2010.px)<-"DPFpx"
DPF2010.px[is.nan(DPF2010.px)] <- 0

#saving
writeRaster(DPF2010.px, "rasters/STM/2010/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (1050m)
DPF2010.ls <- focal(DPF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#DPF2010.ls
#anyNA(DPF2010.ls[])
#plot(DPF2010.ls)

names(DPF2010.ls)<-"DPFls"
DPF2010.ls[is.nan(DPF2010.ls)] <- 0

#saving
writeRaster(DPF2010.ls, "rasters/STM/2010/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[!ls() %in% c("stm.2010", "stm.2020")]) #keeping only raster stack
gc()
#

# santarem 2020
#names(stm.2020)

# isolating forest class pixels
stm.2020.forest.class <- stm.2020[["STM.LULC.2020"]]
stm.2020.forest.class[stm.2020.forest.class==3] <- 1
stm.2020.forest.class[stm.2020.forest.class>1] <- 0
##cheking
#stm.2020.forest.class
#plot(stm.2020.forest.class)

# forested areas with fire and/or degraded
DPF2020<-sum(stm.2020.forest.class, stm.2020[["STM.Fire.2020"]], stm.2020[["STM.Degrad.2020"]], na.rm = T)
DPF2020[DPF2020<=1]<-0
DPF2020[DPF2020>1]<-1
##cheking
#DPF2020
#plot(DPF2020)

# mean dpf cover in pixel scale (150m)
DPF2020.px <- focal(DPF2020, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#DPF2020.px
#anyNA(DPF2020.px[])
#plot(DPF2020.px)

names(DPF2020.px)<-"DPFpx"
DPF2020.px[is.nan(DPF2020.px)] <- 0

#saving
writeRaster(DPF2020.px, "rasters/STM/2020/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (1050m)
DPF2020.ls <- focal(DPF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#DPF2020.ls
#anyNA(DPF2020.ls[])
#plot(DPF2020.ls)

names(DPF2020.ls)<-"DPFls"
DPF2020.ls[is.nan(DPF2020.ls)] <- 0

#saving
writeRaster(DPF2020.ls, "rasters/STM/2020/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[!ls() %in% c("stm.2010", "stm.2020")]) #keeping only raster stack
gc()
#

#

# [SF] secondary forest -- pixel: 5x5 (150m); and landscape: 35x35 (1050m)
# this variable includes forest pixels in SFage raster
# which has less than 25 years for 2010 or less than 35 for 2020

# santarem 2010
#names(stm.2010)

# isolating SF < 25 years old
SF2010 <- stm.2010[["STM.SFage.2010"]]
SF2010[SF2010<25] <- 1
SF2010[SF2010>=25] <- 0
##cheking
#SF2010
#plot(SF2010)

# mean sf cover in pixel scale (150m)
SF2010.px <- focal(SF2010, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#SF2010.px
#anyNA(SF2010.px[])
#plot(SF2010.px)

names(SF2010.px)<-"SFpx"
SF2010.px[is.nan(SF2010.px)] <- 0

#saving
writeRaster(SF2010.px, "rasters/STM/2010/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (1050m)
SF2010.ls <- focal(SF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SF2010.ls
#anyNA(SF2010.ls[])
#plot(SF2010.ls)

names(SF2010.ls)<-"SFls"
SF2010.ls[is.nan(SF2010.ls)] <- 0

#saving
writeRaster(SF2010.ls, "rasters/STM/2010/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[!ls() %in% c("stm.2010", "stm.2020")]) #keeping only raster stack
gc()
#

# santarem 2020
#names(stm.2020)

# isolating SF < 35 years old
SF2020 <- stm.2020[["STM.SFage.2020"]]
SF2020[SF2020<35] <- 1
SF2020[SF2020>=35] <- 0
##cheking
#SF2020
#plot(SF2020)

# mean sf cover in pixel scale (150m)
SF2020.px <- focal(SF2020, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#SF2020.px
#anyNA(SF2020.px[])
#plot(SF2020.px)

names(SF2020.px)<-"SFpx"
SF2020.px[is.nan(SF2020.px)] <- 0

#saving
writeRaster(SF2020.px, "rasters/STM/2020/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (1050m)
SF2020.ls <- focal(SF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SF2020.ls
#anyNA(SF2020.ls[])
#plot(SF2020.ls)

names(SF2020.ls)<-"SFls"
SF2020.ls[is.nan(SF2020.ls)] <- 0

#saving
writeRaster(SF2020.ls, "rasters/STM/2020/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[!ls() %in% c("stm.2010", "stm.2020")]) #keeping only raster stack
gc()
#

#

# [F1] Forest type 1 or Mature forest -- pixel: 5x5 (150m)
# this variable includes forest pixels in LULC raster
# including degraded forest and secondary forest older than 10 years

# santarem 2010
#names(stm.2010)

# isolating forest class pixels
stm.2010.forest.class <- stm.2010[["STM.LULC.2010"]]
stm.2010.forest.class[stm.2010.forest.class==3] <- 1
stm.2010.forest.class[stm.2010.forest.class>1] <- 0
##cheking
#stm.2010.forest.class
#plot(stm.2010.forest.class)

# isolating SF <= 10 years old
stm.2010.SFold <- stm.2010[["STM.SFage.2010"]]
stm.2010.SFold[stm.2010.SFold<=10] <- 1
stm.2010.SFold[stm.2010.SFold>10] <- 0
##cheking
#stm.2010.SFold
#plot(stm.2010.SFold)

# excluding areas with SF <= 10
MF2010<-sum(stm.2010.forest.class, stm.2010.SFold, na.rm = T)
MF2010[MF2010>1]<-0
##cheking
#MF2010
#plot(MF2010)

## including areas with fire and degraded
#MF2010<-sum(MF2010, stm.2010[["STM.Fire.2010"]], stm.2010[["STM.Degrad.2010"]], na.rm = T)
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
writeRaster(MF2010.px, "rasters/STM/2010/MFpx.tif", format="GTiff")

rm(list=ls()[!ls() %in% c("stm.2010", "stm.2020")]) #keeping only raster stack
gc()
#

# santarem 2020
#names(stm.2020)

# isolating forest class pixels
stm.2020.forest.class <- stm.2020[["STM.LULC.2020"]]
stm.2020.forest.class[stm.2020.forest.class==3] <- 1
stm.2020.forest.class[stm.2020.forest.class>1] <- 0
##cheking
#stm.2020.forest.class
#plot(stm.2020.forest.class)

# isolating SF <= 10 years old
stm.2020.SFold <- stm.2020[["STM.SFage.2020"]]
stm.2020.SFold[stm.2020.SFold<=10] <- 1
stm.2020.SFold[stm.2020.SFold>10] <- 0
##cheking
#stm.2020.SFold
#plot(stm.2020.SFold)

# excluding areas with SF <= 10
MF2020<-sum(stm.2020.forest.class, stm.2020.SFold, na.rm = T)
MF2020[MF2020>1]<-0
##cheking
#MF2020
#plot(MF2020)

## including areas with fire and degraded
#MF2020<-sum(MF2020, stm.2020[["STM.Fire.2020"]], stm.2020[["STM.Degrad.2020"]], na.rm = T)
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
writeRaster(MF2020.px, "rasters/STM/2020/MFpx.tif", format="GTiff")

rm(list=ls()[!ls() %in% c("stm.2010", "stm.2020")]) #keeping only raster stack
gc()
#

#

# [F3] Forest type 3 or Total forest -- pixel: 5x5 (150m)
# this variable includes forest pixels in LULC raster
# including degraded forest and secondary forest older than 2 years

# santarem 2010
#names(stm.2010)

# isolating forest class pixels
stm.2010.forest.class <- stm.2010[["STM.LULC.2010"]]
stm.2010.forest.class[stm.2010.forest.class==3] <- 1
stm.2010.forest.class[stm.2010.forest.class>1] <- 0
##cheking
#stm.2010.forest.class
#plot(stm.2010.forest.class)

# isolating SF <= 2 years old
stm.2010.SFyoung <- stm.2010[["STM.SFage.2010"]]
stm.2010.SFyoung[stm.2010.SFyoung<=2] <- 1
stm.2010.SFyoung[stm.2010.SFyoung>2] <- 0
##cheking
#stm.2010.SFyoung
#plot(stm.2010.SFyoung)

# excluding areas with SF <= 2
TF2010<-sum(stm.2010.forest.class, stm.2010.SFyoung, na.rm = T)
TF2010[TF2010>1]<-0
##cheking
#TF2010
#plot(TF2010)

## including areas with fire and degraded
#TF2010<-sum(TF2010, stm.2010[["STM.Fire.2010"]], stm.2010[["STM.Degrad.2010"]], na.rm = T)
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
writeRaster(TF2010.px, "rasters/STM/2010/TFpx.tif", format="GTiff")

rm(list=ls()[!ls() %in% c("stm.2010", "stm.2020")]) #keeping only raster stack
gc()
#

# santarem 2020
#names(stm.2020)

# isolating forest class pixels
stm.2020.forest.class <- stm.2020[["STM.LULC.2020"]]
stm.2020.forest.class[stm.2020.forest.class==3] <- 1
stm.2020.forest.class[stm.2020.forest.class>1] <- 0
##cheking
#stm.2020.forest.class
#plot(stm.2020.forest.class)


# isolating SF <= 2 years old
stm.2020.SFyoung <- stm.2020[["STM.SFage.2020"]]
stm.2020.SFyoung[stm.2020.SFyoung<=2] <- 1
stm.2020.SFyoung[stm.2020.SFyoung>2] <- 0
##cheking
#stm.2020.SFyoung
#plot(stm.2020.SFyoung)

# excluding areas with SF <= 2
TF2020<-sum(stm.2020.forest.class, stm.2020.SFyoung, na.rm = T)
TF2020[TF2020>1]<-0
##cheking
#TF2020
#plot(TF2020)

## including areas with fire and degraded
#TF2020<-sum(TF2020, stm.2020[["STM.Fire.2020"]], stm.2020[["STM.Degrad.2020"]], na.rm = T)
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
writeRaster(TF2020.px, "rasters/STM/2020/TFpx.tif", format="GTiff")

rm(list=ls()[!ls() %in% c("stm.2010", "stm.2020")]) #keeping only raster stack
gc()
#

#

# [SFage] secondary forest age -- pixel: 5x5 (150m); and landscape: 35x35 (1050m)
# this variable is the mean age of secondary forest

# santarem 2010
#names(stm.2010)

# isolating SF < 25 years old
SF2010 <- stm.2010[["STM.SFage.2010"]]
SF2010[SF2010>24] <- NA
##cheking
#SF2010
#plot(SF2010)

# mean sf age in pixel scale (150m)
SFage2010.px <- focal(SF2010, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#SFage2010.px
#anyNA(SFage2010.px[])
#plot(SFage2010.px)

names(SFage2010.px)<-"SFagepx"
SFage2010.px[is.nan(SFage2010.px)] <- 0

#saving
writeRaster(SFage2010.px, "rasters/STM/2010/SFagepx.tif", format="GTiff")

# mean sf age in landscape scale (1050m)
SFage2010.ls <- focal(SF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SFage2010.ls
#anyNA(SFage2010.ls[])
#plot(SFage2010.ls)

names(SFage2010.ls)<-"SFagels"
SFage2010.ls[is.nan(SFage2010.ls)] <- 0

#saving
writeRaster(SFage2010.ls, "rasters/STM/2010/SFagels.tif", format="GTiff")

rm(list=ls()[!ls() %in% c("stm.2010", "stm.2020")]) #keeping only raster stack
gc()
#

# santarem 2020
#names(stm.2020)

# isolating SF < 35 years old
SF2020 <- stm.2020[["STM.SFage.2020"]]
SF2020[SF2020>34] <- NA
##cheking
#SF2020
#plot(SF2020)

# mean sf age in pixel scale (150m)
SFage2020.px <- focal(SF2020, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#SFage2020.px
#anyNA(SFage2020.px[])
#plot(SFage2020.px)

names(SFage2020.px)<-"SFagepx"
SFage2020.px[is.nan(SFage2020.px)] <- 0

#saving
writeRaster(SFage2020.px, "rasters/STM/2020/SFagepx.tif", format="GTiff")
#

# mean sf cover in landscape scale (1050m)
SFage2020.ls <- focal(SF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SFage2020.ls
#anyNA(SFage2020.ls[])
#plot(SFage2020.ls)

names(SFage2020.ls)<-"SFagels"
SFage2020.ls[is.nan(SFage2020.ls)] <- 0

#saving
writeRaster(SFage2020.ls, "rasters/STM/2020/SFagels.tif", format="GTiff")

rm(list=ls()[!ls() %in% c("stm.2010", "stm.2020")]) #keeping only raster stack
gc()
#

#

# [Edge.dist] distance to forest edge
# this variable is the mean distance of F1 to the edge

# santarem 2010
#names(stm.2010)

# isolating forest class pixels
stm.2010.forest.class <- stm.2010[["STM.LULC.2010"]]
stm.2010.forest.class[stm.2010.forest.class==3] <- 1
stm.2010.forest.class[stm.2010.forest.class>1] <- 0
##cheking
#stm.2010.forest.class
#plot(stm.2010.forest.class)

# isolating SF <= 10 years old
stm.2010.SFold <- stm.2010[["STM.SFage.2010"]]
stm.2010.SFold[stm.2010.SFold<=10] <- 1
stm.2010.SFold[stm.2010.SFold>10] <- 0
##cheking
#stm.2010.SFold
#plot(stm.2010.SFold)

# excluding areas with SF <= 10
MF2010<-sum(stm.2010.forest.class, stm.2010.SFold, na.rm = T)
MF2010[MF2010>1]<-0
##cheking
#MF2010
#plot(MF2010)

## including areas with fire and degraded
#MF2010<-sum(MF2010, stm.2010[["STM.Fire.2010"]], stm.2010[["STM.Degrad.2010"]], na.rm = T)
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

rm(stm.2010.forest.class); rm(stm.2010.SFold); gc() #keeping only raster in use

edge.dist.2010 <- rasterDistance(pts2010, core2010, reference = MF2010, scale=TRUE)
##cheking
#edge.dist.2010
#anyNA(edge.dist.2010[])
#plot(edge.dist.2010)

names(edge.dist.2010)<-"edgedist"

#saving
writeRaster(edge.dist.2010, "rasters/STM/2010/edgedist.tif", format="GTiff")

rm(pts2010); rm(core2010); rm(edge.dist.2010); gc() #keeping only raster in use

# santarem 2020
#names(stm.2020)

# isolating forest class pixels
stm.2020.forest.class <- stm.2020[["STM.LULC.2020"]]
stm.2020.forest.class[stm.2020.forest.class==3] <- 1
stm.2020.forest.class[stm.2020.forest.class>1] <- 0
##cheking
#stm.2020.forest.class
#plot(stm.2020.forest.class)

# isolating SF <= 10 years old
stm.2020.SFold <- stm.2020[["STM.SFage.2020"]]
stm.2020.SFold[stm.2020.SFold<=10] <- 1
stm.2020.SFold[stm.2020.SFold>10] <- 0
##cheking
#stm.2020.SFold
#plot(stm.2020.SFold)

# excluding areas with SF <= 10
MF2020<-sum(stm.2020.forest.class, stm.2020.SFold, na.rm = T)
MF2020[MF2020>1]<-0
##cheking
#MF2020
#plot(MF2020)

## including areas with fire and degraded
#MF2020<-sum(MF2020, stm.2020[["STM.Fire.2020"]], stm.2020[["STM.Degrad.2020"]], na.rm = T)
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

rm(stm.2020.forest.class); rm(stm.2020.SFold); gc() #keeping only raster in use

edge.dist.2020 <- rasterDistance(pts2020, core2020, reference = MF2020, scale=TRUE)
##cheking
#edge.dist.2020
#anyNA(edge.dist.2020[])
#plot(edge.dist.2020)

names(edge.dist.2020)<-"edgedist"

#saving
writeRaster(edge.dist.2020, "rasters/STM/2020/edgedist.tif", format="GTiff")

rm(pts2020); rm(core2020); rm(edge.dist.2020); gc() #keeping only raster in use

#

#

# [Edge] forest edge -- pixel: 5x5 (150m); and landscape: 35x35 (1050m)
# this variable is the mean of edge area based on mature forest

# santarem 2010

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
writeRaster(edge2010.px, "rasters/STM/2010/edgepx.tif", format="GTiff")
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
writeRaster(edge2010.ls, "rasters/STM/2010/edgels.tif", format="GTiff")
#

# santarem 2020

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
writeRaster(edge2020.px, "rasters/STM/2020/edgepx.tif", format="GTiff")
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
writeRaster(edge2020.ls, "rasters/STM/2020/edgels.tif", format="GTiff")

rm(list=ls()[!ls() %in% c("stm.2010", "stm.2020")]) #keeping only raster stack
gc()
#

#

# [waterdist] distance to water
# this variable is the euclidean distance between pixels to water body

stm.water.shp <- readOGR(dsn = "rasters/STM/input", layer = "stm-water") #reading shapefile
stm.water.pts <- as(stm.water.shp, "SpatialPointsDataFrame")
#cheking
#stm.water.pts
#plot(stm.water.pts)


pts2010 <- rasterToPoints(stm.2020[["STM.LULC.2020"]], spatial=TRUE)

water.dist <- rasterDistance(pts2010, stm.water.pts, reference = stm.2020[["STM.LULC.2020"]], scale=TRUE)
##cheking
#water.dist
#anyNA(water.dist[])
#plot(water.dist)

names(water.dist)<-"waterdist"

#saving
writeRaster(water.dist, "rasters/STM/2010/waterdist.tif", format="GTiff")
writeRaster(water.dist, "rasters/STM/2020/waterdist.tif", format="GTiff")

rm(list=ls()[!ls() %in% c("stm.2010", "stm.2020")]) #keeping only raster stack
gc()
#

#

## [roaddist] distance to road
## this variable is the euclidean distance between pixels to road
#
#stm.road.shp <- readOGR(dsn = "rasters/STM/input", layer = "stm-road") #reading shapefile
#stm.road.pts <- as(stm.road.shp, "SpatialPointsDataFrame")
##cheking
##stm.road.pts
##plot(stm.road.pts)
#
#
#pts2010 <- rasterToPoints(stm.2020[["STM.LULC.2020"]], spatial=TRUE)
#
#road.dist <- rasterDistance(pts2010, stm.road.pts, reference = stm.2020[["STM.LULC.2020"]], scale=TRUE)
###cheking
##road.dist
##anyNA(road.dist[])
##plot(road.dist)
#
#names(road.dist)<-"roaddist"
#
##saving
#writeRaster(road.dist, "rasters/STM/2010/roaddist.tif", format="GTiff")
#writeRaster(road.dist, "rasters/STM/2020/roaddist.tif", format="GTiff")
#
#rm(list=ls()[!ls() %in% c("stm.2010", "stm.2020")]) #keeping only raster stack
#gc()
##
#
##

# [temp] historical mean annual temperature from worldclim -- 1970:2000
mean.temp <- raster("C:/Users/miral/Dropbox/GIS/wc2.1_30s_bio/wc2.1_30s_bio_1.tif")

# Conversion of rasters into same extent
stm.mean.temp <- resample(mean.temp, stm.2010[["STM.LULC.2010"]], method='bilinear')
##cheking
#stm.mean.temp
#anyNA(stm.mean.temp[])
#plot(stm.mean.temp)

names(stm.mean.temp)<-"meantemp"

values(stm.mean.temp) <- rescale(values(stm.mean.temp)) # scale values from min=0 to max=1
stm.mean.temp.clipp <- mask(stm.mean.temp, stm.2010[["STM.LULC.2010"]]) # clipping to keep only municipilaty 

#saving
writeRaster(stm.mean.temp.clipp, "rasters/STM/2010/meantemp.tif", format="GTiff")
writeRaster(stm.mean.temp.clipp, "rasters/STM/2020/meantemp.tif", format="GTiff")
#

#

# [precip] historical mean annual precipitation from worldclim -- 1970:2000
mean.precip <- raster("C:/Users/miral/Dropbox/GIS/wc2.1_30s_bio/wc2.1_30s_bio_12.tif")

# Conversion of rasters into same extent
stm.mean.precip <- resample(mean.precip, stm.2010[["STM.LULC.2010"]], method='bilinear')
##cheking
#stm.mean.precip
#anyNA(stm.mean.precip[])
#plot(stm.mean.precip)

names(stm.mean.precip)<-"meanprecip"

values(stm.mean.precip) <- rescale(values(stm.mean.precip)) # scale values from min=0 to max=1
stm.mean.precip.clipp <- mask(stm.mean.precip, stm.2010[["STM.LULC.2010"]]) # clipping to keep only municipilaty 

#saving
writeRaster(stm.mean.precip.clipp, "rasters/STM/2010/meanprecip.tif", format="GTiff")
writeRaster(stm.mean.precip.clipp, "rasters/STM/2020/meanprecip.tif", format="GTiff")
#

#

# [elev] elevation data derived from the SRTM
elev <- raster("C:/Users/miral/Dropbox/GIS/wc2.1_30s_bio/wc2.1_30s_elev.tif")

# Conversion of rasters into same extent
stm.elev <- resample(elev, stm.2010[["STM.LULC.2010"]], method='bilinear')
##cheking
#stm.elev
#anyNA(stm.elev[])
#plot(stm.elev)

names(stm.elev)<-"elev"

values(stm.elev) <- rescale(values(stm.elev)) # scale values from min=0 to max=1
stm.elev.clipp <- mask(stm.elev, stm.2010[["STM.LULC.2010"]]) # clipping to keep only municipilaty 

#saving
writeRaster(stm.elev.clipp, "rasters/STM/2010/elev.tif", format="GTiff")
writeRaster(stm.elev.clipp, "rasters/STM/2020/elev.tif", format="GTiff")
#

#


rm(list=ls()) #end
gc()
#

#


#
