
#' @title Cost-benefit of conservation actions in Amazon
#' @description script to build exploratory variables from 
#' land use - land cover, secondary forest, edge, 
#' degradation (loggin and fire), temperature, precipitation, 
#' elevation and distances to road and water body
#' in Paragominas municipality - PA;
#' this set of exploratory variables is used in fit 
#' regional species distribution models 

#### pre-setting ####
memory.limit(1000000)

#### loading required packages ####
library(tidyverse)
library(raster)
library(rgdal)
library(rgeos)
library(sf)
library(spatialEco)
library(scales)
library(virtualspecies)
library(usdm)


#### creating directories ####
dir.create("rasters/PGM/2010_real", recursive = T)
dir.create("rasters/PGM/2020_real", recursive = T)
dir.create("rasters/PGM/2020_avoiddeforest", recursive = T)
dir.create("rasters/PGM/2020_avoiddegrad", recursive = T)
dir.create("rasters/PGM/2020_avoidboth", recursive = T)
dir.create("rasters/PGM/2020_restor_wo_avoid", recursive = T)
dir.create("rasters/PGM/2020_restor_n_avoid", recursive = T)



#### importing input rasters ####

# shapefile paragominas
pgm.shp <- readOGR(dsn = "shapes", layer = "Paragominas_Mask_R3")

#
#



# land use land cover from mapbiomas collection 7 [2010 and 2020]
# [PGM] paragominas
pgm.lulc <- stack(c("rasters/PGM/input/pgm-2010-lulc-mapbiomas-brazil-collection-70.tif",
                    "rasters/PGM/input/pgm-2020-lulc-mapbiomas-brazil-collection-70.tif"))
names(pgm.lulc) <- c("pgm.lulc.2010real", "pgm.lulc.2020real")
#checking
#pgm.lulc
#plot(pgm.lulc)
#sort(unique(values(pgm.lulc[["pgm.lulc.2010real"]])))

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
pgm.lulc.2010.forest.class[pgm.lulc.2010.forest.class>1] <- 0

pgm.lulc.2020.forest.class <- pgm.lulc[["pgm.lulc.2020real"]]
pgm.lulc.2020.forest.class[pgm.lulc.2020.forest.class==3] <- 1
pgm.lulc.2020.forest.class[pgm.lulc.2020.forest.class>1] <- 0

#
#



# distances to road and rivers, and elevation

dist.road <- raster("rasters/PGM/input/dist_road_pgm.tif")
dist.road <- projectRaster(dist.road, crs = "+proj=longlat +datum=WGS84 +no_defs")
dist.road <- resample(dist.road, pgm.lulc.2010.forest.class, method='bilinear')
values(dist.road)[is.na(values(dist.road))] = 0
dist.road <- mask(dist.road, pgm.shp)

#saving
writeRaster(dist.road, "rasters/PGM/2010_real/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/PGM/2020_real/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/PGM/2020_avoiddeforest/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/PGM/2020_avoiddegrad/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/PGM/2020_avoidboth/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/PGM/2020_restor_wo_avoid/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/PGM/2020_restor_n_avoid/distroad.tif", format="GTiff", overwrite=T)

#

dist.river <- raster("rasters/PGM/input/dist_river_pgm.tif")
dist.river <- projectRaster(dist.river, crs = "+proj=longlat +datum=WGS84 +no_defs")
dist.river <- resample(dist.river, pgm.lulc.2010.forest.class, method='bilinear')
values(dist.river)[is.na(values(dist.river))] = 0
dist.river <- mask(dist.river, pgm.shp)

#saving
writeRaster(dist.river, "rasters/PGM/2010_real/distriver.tif", format="GTiff", overwrite=T)
writeRaster(dist.river, "rasters/PGM/2020_real/distriver.tif", format="GTiff", overwrite=T)
writeRaster(dist.river, "rasters/PGM/2020_avoiddeforest/distriver.tif", format="GTiff", overwrite=T)
writeRaster(dist.river, "rasters/PGM/2020_avoiddegrad/distriver.tif", format="GTiff", overwrite=T)
writeRaster(dist.river, "rasters/PGM/2020_avoidboth/distriver.tif", format="GTiff", overwrite=T)
writeRaster(dist.river, "rasters/PGM/2020_restor_wo_avoid/distriver.tif", format="GTiff", overwrite=T)
writeRaster(dist.river, "rasters/PGM/2020_restor_n_avoid/distriver.tif", format="GTiff", overwrite=T)

#

elevation <- raster("rasters/PGM/input/elevation_pgm.tif")
elevation <- projectRaster(elevation, crs = "+proj=longlat +datum=WGS84 +no_defs")
elevation <- resample(elevation, pgm.lulc.2010.forest.class, method='bilinear')
values(elevation)[is.na(values(elevation))] = 0
elevation <- mask(elevation, pgm.shp)

#saving
writeRaster(elevation, "rasters/PGM/2010_real/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/PGM/2020_real/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/PGM/2020_avoiddeforest/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/PGM/2020_avoiddegrad/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/PGM/2020_avoidboth/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/PGM/2020_restor_wo_avoid/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/PGM/2020_restor_n_avoid/elevation.tif", format="GTiff", overwrite=T)

#


#### candidate areas for restoration scenarios ####

#isolating deforestation class pixels (crops, pasture)

deforestation.class.list <- c(15,39,41,48)

candidate.areas.total <- pgm.lulc[["pgm.lulc.2010real"]]

values(candidate.areas.total)[values(candidate.areas.total) %in% deforestation.class.list] = 1
values(candidate.areas.total)[values(candidate.areas.total) > 1] = 0
names(candidate.areas.total) <- "restoration.candidate.areas"
#plot(candidate.areas.total)

#select pixels based on proximity to water (<500m), slope (>25°) and proximity to forest (<1000m)
dist.river.all <- dist.river
values(dist.river)[values(dist.river) <= 500] = 1
values(dist.river)[values(dist.river) > 500] = NA
#plot(dist.river)

candidate.areas.water <- candidate.areas.total
candidate.areas.water <- mask(candidate.areas.water, dist.river)
values(candidate.areas.water)[is.na(values(candidate.areas.water))] = 0
candidate.areas.water <- mask(candidate.areas.water, pgm.shp)
#plot(candidate.areas.water)



slope <- terrain(elevation, opt = 'slope', unit = 'degrees', neighbors=8)
values(slope)[values(slope) < 45] = NA
values(slope)[values(slope) >= 45] = 1
#plot(slope)

candidate.areas.slope <- candidate.areas.total
candidate.areas.slope <- mask(candidate.areas.slope, slope)
values(candidate.areas.slope)[is.na(values(candidate.areas.slope))] = 0
candidate.areas.slope <- mask(candidate.areas.slope, pgm.shp)
#plot(candidate.areas.slope)


forest.class <- pgm.lulc[["pgm.lulc.2010real"]]
forest.class[forest.class==0] <- NA
forest.class[forest.class==3] <- 1
forest.class[forest.class>1] <- 0

deforest.pts <- rasterToPoints(forest.class, spatial=TRUE)
deforest.core <- deforest.pts[deforest.pts$pgm.lulc.2010real == "1",]

deforest.dist <- rasterDistance(deforest.pts, deforest.core, reference = candidate.areas.total, scale=TRUE)
values(deforest.dist)[values(deforest.dist) == 0] = NA
values(deforest.dist)[values(deforest.dist) > 0.15] = NA
values(deforest.dist)[values(deforest.dist) <= 0.15] = 1
#plot(deforest.dist)

candidate.areas.forest <- candidate.areas.total
candidate.areas.forest <- mask(candidate.areas.forest, deforest.dist)
values(candidate.areas.forest)[is.na(values(candidate.areas.forest))] = 0
candidate.areas.forest <- mask(candidate.areas.forest, pgm.shp)
#plot(candidate.areas.forest)


candidate.areas.final <- sum(candidate.areas.water,candidate.areas.slope,candidate.areas.forest)
values(candidate.areas.final)[values(candidate.areas.final) >= 1] = 1
#plot(candidate.areas.final)

#select rural properties with less than 50% of forest cover in 2010

#import rural properties shapefiles and data from SISCAR
#https://www.car.gov.br/publico/municipios/downloads
pgm.car <- readOGR(dsn = "rasters/PGM/input/SHAPE_1505502_CAR_Paragominas", layer = "AREA_IMOVEL")
#head(pgm.car@data)
#plot(pgm.car, add=T)


#calculate forest cover in each property and select properties with <50%

pgm.car@data$FOREST_COVER <- NA
j=nrow(pgm.car@data)
for (i in pgm.car$COD_IMOVEL) {
  
  rural.property <- pgm.car[pgm.car$COD_IMOVEL==i,]
  forest.cover <- crop(forest.class, extent(rural.property))
  forest.cover <- mask(forest.cover, rural.property)
  pgm.car[pgm.car$COD_IMOVEL==i,"FOREST_COVER"] <- tapply(area(forest.cover), forest.cover[], sum)[2]*100
  j=j-1
  cat("\n>", j, "out of", nrow(pgm.car@data), "properties left<\n")
  
}


pgm.car@data$FOREST_COVER_PP <- ceiling((pgm.car@data$FOREST_COVER/pgm.car@data$NUM_AREA)*100)

#select properties
pgm.car.restoration.candidates <- pgm.car[!is.na(pgm.car$FOREST_COVER_PP),]
pgm.car.restoration.candidates <- pgm.car.restoration.candidates[pgm.car.restoration.candidates$FOREST_COVER_PP <= 50,]
#head(pgm.car.restoration.candidates@data)
#nrow(pgm.car.restoration.candidates@data)


#filter candidate areas for restoration in properties with less than 50% forest cover
candidate.areas.final.copy <- candidate.areas.final
candidate.areas.final <- mask(candidate.areas.final, pgm.car.restoration.candidates)
values(candidate.areas.final)[is.na(values(candidate.areas.final))] = 0
candidate.areas.final <- mask(candidate.areas.final, pgm.shp)
#plot(candidate.areas.final)
#plot(pgm.car.restoration.candidates, add=T)


#forest cover increment -- adding the candidate areas to forest cover
pgm.car.restoration.candidates@data$FOREST_COVER_INCREMENT <- NA
j=nrow(pgm.car.restoration.candidates@data)
#i="PA-1505502-39CCE4418D2D487F9AC0FD3045A374CF"
for (i in pgm.car.restoration.candidates$COD_IMOVEL) {
  
  rural.property <- pgm.car.restoration.candidates[pgm.car.restoration.candidates$COD_IMOVEL==i,]
  
  forest.cover <- crop(forest.class, extent(rural.property))
  forest.cover <- mask(forest.cover, rural.property)
  
  restored.cover <- crop(candidate.areas.final, extent(rural.property))
  restored.cover <- mask(restored.cover, rural.property)
  
  forest.cover.increment <- sum(forest.cover, restored.cover)
  
  pgm.car.restoration.candidates[pgm.car.restoration.candidates$COD_IMOVEL==i,"FOREST_COVER_INCREMENT"] <- tapply(area(forest.cover.increment), forest.cover.increment[], sum)[2]*100
  j=j-1
  cat("\n>", j, "out of", nrow(pgm.car.restoration.candidates@data), "properties left<\n")
  
}

pgm.car.restoration.candidates@data$FOREST_COVER_INCREMENT_PP <- ceiling((pgm.car.restoration.candidates@data$FOREST_COVER_INCREMENT/pgm.car.restoration.candidates@data$NUM_AREA)*100)


#select properties with more than 80% forest cover
pgm.car.restoration.candidates <- pgm.car.restoration.candidates[!is.na(pgm.car.restoration.candidates$FOREST_COVER_INCREMENT_PP),]
pgm.car.restoration.candidates.m80 <- pgm.car.restoration.candidates[pgm.car.restoration.candidates$FOREST_COVER_INCREMENT_PP > 80,]
#head(pgm.car.restoration.candidates.m80@data)
#nrow(pgm.car.restoration.candidates.m80@data)



candidate.areas.final.copy <- candidate.areas.final
j=nrow(pgm.car.restoration.candidates.m80@data)
#i="PA-1505502-75B3625903EB4F6A9C36FF2B97B67282"
for (i in pgm.car.restoration.candidates.m80$COD_IMOVEL) {
  
  rural.property <- pgm.car.restoration.candidates.m80[pgm.car.restoration.candidates.m80$COD_IMOVEL==i,]
  
  forest.cover <- crop(forest.class, extent(rural.property))
  forest.cover <- mask(forest.cover, rural.property)
  
  restored.cover <- crop(candidate.areas.final, extent(rural.property))
  restored.cover <- mask(restored.cover, rural.property)
  
  
  while (pgm.car.restoration.candidates.m80@data[pgm.car.restoration.candidates.m80@data$COD_IMOVEL==i,"FOREST_COVER_INCREMENT_PP"] > 80) {
    
    restored.cover[restored.cover[]==1] <- sample(c(1,0), size = length(restored.cover[restored.cover[]==1]), replace = T, prob = c(0.9,0.1))
    forest.cover.increment <- sum(forest.cover, restored.cover)
    pgm.car.restoration.candidates.m80@data[pgm.car.restoration.candidates.m80@data$COD_IMOVEL==i,"FOREST_COVER_INCREMENT"] <- tapply(area(forest.cover.increment), forest.cover.increment[], sum)[2]*100
    pgm.car.restoration.candidates.m80@data[pgm.car.restoration.candidates.m80@data$COD_IMOVEL==i,"FOREST_COVER_INCREMENT_PP"] <- ceiling((pgm.car.restoration.candidates.m80@data[pgm.car.restoration.candidates.m80@data$COD_IMOVEL==i,"FOREST_COVER_INCREMENT"]/pgm.car.restoration.candidates.m80@data[pgm.car.restoration.candidates.m80@data$COD_IMOVEL==i,"NUM_AREA"])*100)
    
  }
  
  candidate.areas.final[restored.cover][candidate.areas.final[restored.cover]==0]<-0
  j=j-1
  cat("\n>", j, "out of", nrow(pgm.car.restoration.candidates.m80@data), "properties left<\n")
  
}

#length(candidate.areas.final.copy[candidate.areas.final.copy[]==1])
#length(candidate.areas.final[candidate.areas.final[]==1])


#select properties with still less than 80% forest cover

pgm.car.restoration.candidates.l80 <- pgm.car.restoration.candidates[pgm.car.restoration.candidates$FOREST_COVER_INCREMENT_PP < 80,]

#exclding properties with geometry problems and/or at municipality border
exclude <- c("PA-1505502-E52590A9031A4B67AD3C4F4AFF462E73", "PA-1505502-011CACBDA2D34B10AD9083F6089BF648",
             "PA-1505502-4A18BDFF0E3D407F9478AB1601EFF71A", "PA-1505502-461582E49CA240FB82B4180545842E33",
             "PA-1505502-DFEAB2B557664C6D9F30A3BE6A0ED860", "PA-1505502-96F0C6BCAED44ADBBCD438EE67747333",
             "PA-1505502-28DECEFD42BF46ECB91D202D82189AF9", "PA-1505502-1E8540456B3B4A9DA93434730C61DEA8",
             "PA-1505502-B97842AD547B44CD8FC8AE4A9B320EFB", "PA-1505502-6403C1FC323D44BE80C405867921410F",
             "PA-1505502-B3C59867F0EB48CA8BDBDA4C3BC286DF", "PA-1505502-87EAEADFA1A94AC49823A7B59811B385",
             "PA-1505502-5C865CF92FE9434A9E8B716147B3753A", "PA-1505502-94651AAF903F487DBC28C0D9783A0255",
             "PA-1505502-8B47713EF21441BB89CE6A2FAA3AE944", "PA-1505502-F533313E45BE45AF9B13EDF30A2F3479",
             "PA-1505502-FB027251AC494AECAECD646464AEF521")

pgm.car.restoration.candidates.l80 <- pgm.car.restoration.candidates.l80[-which(pgm.car.restoration.candidates.l80$COD_IMOVEL %in% exclude),]
#head(pgm.car.restoration.candidates.l80@data)
#nrow(pgm.car.restoration.candidates.l80@data)

#candidate.areas.final <- candidate.areas.final.copy

candidate.areas.final.copy <- candidate.areas.final
j=nrow(pgm.car.restoration.candidates.l80@data)
#i="PA-1505502-8B1A6744B90E42C6BC856664188651AF"
for (i in pgm.car.restoration.candidates.l80$COD_IMOVEL) {
  
  rural.property <- pgm.car.restoration.candidates.l80[pgm.car.restoration.candidates.l80$COD_IMOVEL==i,]
  
  forest.cover <- crop(forest.class, extent(rural.property))
  forest.cover <- mask(forest.cover, rural.property)
  
  restored.cover <- crop(candidate.areas.final, extent(rural.property))
  restored.cover <- mask(restored.cover, rural.property)
  
  forest.cover.increment <- sum(forest.cover, restored.cover)
  
  while (pgm.car.restoration.candidates.l80@data[pgm.car.restoration.candidates.l80@data$COD_IMOVEL==i,"FOREST_COVER_INCREMENT_PP"] < 80) {
    
    forest.cover.increment[forest.cover.increment[]==0] <- sample(c(0,1), size = length(forest.cover.increment[forest.cover.increment[]==0]), replace = T, prob = c(0.9,0.1))
    
    pgm.car.restoration.candidates.l80@data[pgm.car.restoration.candidates.l80@data$COD_IMOVEL==i,"FOREST_COVER_INCREMENT"] <- tapply(area(forest.cover.increment), forest.cover.increment[], sum)[2]*100
    pgm.car.restoration.candidates.l80@data[pgm.car.restoration.candidates.l80@data$COD_IMOVEL==i,"FOREST_COVER_INCREMENT_PP"] <- ceiling((pgm.car.restoration.candidates.l80@data[pgm.car.restoration.candidates.l80@data$COD_IMOVEL==i,"FOREST_COVER_INCREMENT"]/pgm.car.restoration.candidates.l80@data[pgm.car.restoration.candidates.l80@data$COD_IMOVEL==i,"NUM_AREA"])*100)
    
  }
  
  restored.cover.update <- forest.cover.increment-forest.cover
  candidate.areas.final[restored.cover.update][candidate.areas.final[restored.cover.update]==0]<-1
  j=j-1
  cat("\n>", j, "out of", nrow(pgm.car.restoration.candidates.l80@data), "properties left<\n")
  
}

#length(candidate.areas.final.copy[candidate.areas.final.copy[]==1])
#length(candidate.areas.final[candidate.areas.final[]==1])

candidate.areas.final <- mask(candidate.areas.final, pgm.shp)
#plot(candidate.areas.final)
#plot(pgm.car.restoration.candidates, add=T)
#candidate.areas.final <- raster("rasters/PGM/input/restoration_candidate_areas.tif")

rm(list=ls()[!ls() %in% c("pgm.shp", "pgm.lulc","pgm.lulc.2010.forest.class", "pgm.lulc.2020.forest.class", "dist.road",
                          "dist.river.all", "elevation", "candidate.areas.final")]) #keeping only raster stack
gc()



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
#plot(pgm.deter.yearx)


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

#non-degraded sites will be considered with 300 years following (BIB)
pgm.degrad[["pgm.degrad.2010real"]][pgm.degrad[["pgm.degrad.2010real"]]>23] <- 300
pgm.degrad[["pgm.degrad.2020real"]][pgm.degrad[["pgm.degrad.2020real"]]>33] <- 300

# isolating degraded forest class pixels
pgm.degrad.2010.forest.class <- pgm.degrad[["pgm.degrad.2010real"]]
pgm.degrad.2010.forest.class[pgm.degrad.2010.forest.class>23]<-NA
pgm.degrad.2010.forest.class[!is.na(pgm.degrad.2010.forest.class)] <- 1
pgm.degrad.2010.forest.class<-sum(pgm.lulc.2010.forest.class, pgm.degrad.2010.forest.class, na.rm=T)
pgm.degrad.2010.forest.class[pgm.degrad.2010.forest.class<2]<-0
pgm.degrad.2010.forest.class[pgm.degrad.2010.forest.class==2]<-1



pgm.degrad.2020.forest.class <- pgm.degrad[["pgm.degrad.2020real"]]
pgm.degrad.2020.forest.class[pgm.degrad.2020.forest.class>33]<-NA
pgm.degrad.2020.forest.class[!is.na(pgm.degrad.2020.forest.class)] <- 1
pgm.degrad.2020.forest.class<-sum(pgm.lulc.2020.forest.class, pgm.degrad.2020.forest.class, na.rm=T)
pgm.degrad.2020.forest.class[pgm.degrad.2020.forest.class<2]<-0
pgm.degrad.2020.forest.class[pgm.degrad.2020.forest.class==2]<-1

rm(list=ls()[!ls() %in% c("pgm.shp", "pgm.lulc","pgm.lulc.2010.forest.class", "pgm.lulc.2020.forest.class", "dist.road",
                          "dist.river.all", "elevation", "candidate.areas.final", "pgm.degrad", "pgm.degrad.2010.forest.class",
                          "pgm.degrad.2020.forest.class")]) #keeping only raster stack
gc()



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

# isolating secondary forest class pixels
pgm.sfage.2010.all.class <- pgm.sfage[["pgm.sfage.2010real"]]
pgm.sfage.2010.all.class[pgm.sfage.2010.all.class>0] <- 1
pgm.sfage.2010.all.class[pgm.sfage.2010.all.class<1] <- 0

pgm.sfage.2020.all.class <- pgm.sfage[["pgm.sfage.2020real"]]
pgm.sfage.2020.all.class[pgm.sfage.2020.all.class>0] <- 1
pgm.sfage.2020.all.class[pgm.sfage.2020.all.class<1] <- 0


rm(list=ls()[!ls() %in% c("pgm.shp", "pgm.lulc","pgm.lulc.2010.forest.class", "pgm.lulc.2020.forest.class", "dist.road",
                          "dist.river.all", "elevation", "candidate.areas.final", "pgm.degrad", "pgm.degrad.2010.forest.class",
                          "pgm.degrad.2020.forest.class", "pgm.sfage", "pgm.sfage.2010.all.class", "pgm.sfage.2020.all.class")]) #keeping only raster stack
gc()



#
#



# function to make sure all raster are fully stacked
#intersect_mask <- function(x){
#  values_x <- getValues(x)
#  inter_x <- values_x %*% rep(1,nlayers(x))
#  mask <- setValues(subset(x,1),values = (inter_x>0))
#  return(mask)
#}

#

#######################################################################################################################

#### setting explanatory variables ####
# [UPF] undisturbed primary forest -- pixel: 5x5 (150m); and landscape: 35x35 (1050m)
# this variable includes all forest pixels in LULC raster (value == 3)
# excluding those with age < 25 in 2010 SF raster or age < 35 in 2020 SF raster
# excluding pixels degraded

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

rm(list=ls()[ls() %in% c("UPF2010.px", "UPF2010.ls")]) #keeping only raster stack
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

rm(list=ls()[ls() %in% c("UPF.avoiddegrad.px", "UPF.avoiddegrad.ls")]) #keeping only raster stack
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

rm(list=ls()[ls() %in% c("UPF.avoiddefor.px", "UPF.avoiddefor.ls")]) #keeping only raster stack
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

rm(list=ls()[ls() %in% c("UPF2020.px", "UPF2020.ls")]) #keeping only raster stack
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

rm(list=ls()[ls() %in% c("pgm.degrad.2010.forest.class", "DPF2010.px", "DPF2010.ls")]) #keeping only raster stack
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

rm(list=ls()[ls() %in% c("pgm.degrad.2020.forest.class", "DPF2020.px", "DPF2020.ls")]) #keeping only raster stack
gc()



#


#######################################################################################################################

# [TSD] time since degradation -- pixel: 5x5 (150m); and landscape: 35x35 (1050m)
# this variable is the mean time since a degradation event

# scenario 2010
TSD2010 <- pgm.degrad[["pgm.degrad.2010real"]]
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

rm(list=ls()[ls() %in% c("TSD2010.px", "TSD2010.ls")]) #keeping only raster stack
gc()



#


# scenario 2020 -- avoid degradation
TSD2010.recovery10 <- calc(TSD2010, fun=function(x){ifelse(is.na(x), x, ifelse(x==0, x, ifelse(x==300, x, x+10)))})
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

rm(list=ls()[ls() %in% c("TSD2010.recovery10.px", "TSD2010.recovery10.ls")]) #keeping only raster stack
gc()



#


# scenario 2020
TSD2020 <- pgm.degrad[["pgm.degrad.2020real"]]
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

rm(list=ls()[ls() %in% c("TSD2020.px", "TSD2020.ls")]) #keeping only raster stack
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

rm(list=ls()[ls() %in% c("SF2010.px", "SF2010.ls")]) #keeping only raster stack
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

rm(list=ls()[ls() %in% c("SF2020.px", "SF2020.ls")]) #keeping only raster stack
gc()



#


# scenario -- restoration and avoid
SF2010.restore10 <- sum(SF2010, candidate.areas.final, na.rm = T)
SF2010.restore10[SF2010.restore10>1]<-1
#plot(SF2010.restore10)

# mean sf cover in pixel scale (150m)
SF2010.restore10.px <- focal(SF2010.restore10, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#SF2010.restore10.px
#anyNA(SF2010.restore10.px[])
#plot(SF2010.restore10.px)

names(SF2010.restore10.px)<-"SFpx"
SF2010.restore10.px[is.nan(SF2010.restore10.px)] <- 0
SF2010.restore10.px <- mask(SF2010.restore10.px, pgm.shp)

#saving
writeRaster(SF2010.restore10.px, "rasters/PGM/2020_restor_n_avoid/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (1050m)
SF2010.restore10.ls <- focal(SF2010.restore10, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SF2010.restore10.ls
#anyNA(SF2010.restore10.ls[])
#plot(SF2010.restore10.ls)

names(SF2010.restore10.ls)<-"SFls"
SF2010.restore10.ls[is.nan(SF2010.restore10.ls)] <- 0
SF2010.restore10.ls <- mask(SF2010.restore10.ls, pgm.shp)

#saving
writeRaster(SF2010.restore10.ls, "rasters/PGM/2020_restor_n_avoid/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2010.restore10.px", "SF2010.restore10.ls")]) #keeping only raster stack
gc()



#


# scenario -- restoration without avoid
SF2020.restore10 <- sum(SF2020, candidate.areas.final, na.rm = T)
SF2020.restore10[SF2020.restore10>1]<-1
#plot(SF2020.restore10)

# mean sf cover in pixel scale (150m)
SF2020.restore10.px <- focal(SF2020.restore10, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#SF2020.restore10.px
#anyNA(SF2020.restore10.px[])
#plot(SF2020.restore10.px)

names(SF2020.restore10.px)<-"SFpx"
SF2020.restore10.px[is.nan(SF2020.restore10.px)] <- 0
SF2020.restore10.px <- mask(SF2020.restore10.px, pgm.shp)

#saving
writeRaster(SF2020.restore10.px, "rasters/PGM/2020_restor_wo_avoid/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (1050m)
SF2020.restore10.ls <- focal(SF2020.restore10, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SF2020.restore10.ls
#anyNA(SF2020.restore10.ls[])
#plot(SF2020.restore10.ls)

names(SF2020.restore10.ls)<-"SFls"
SF2020.restore10.ls[is.nan(SF2020.restore10.ls)] <- 0
SF2020.restore10.ls <- mask(SF2020.restore10.ls, pgm.shp)

#saving
writeRaster(SF2020.restore10.ls, "rasters/PGM/2020_restor_wo_avoid/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020.restore10.px", "SF2020.restore10.ls")]) #keeping only raster stack
gc()



#


#######################################################################################################################

# [SFage] secondary forest age -- pixel: 5x5 (150m); and landscape: 35x35 (1050m)
# this variable is the mean age of secondary forest

# scenario 2010
SFage2010 <- pgm.sfage[["pgm.sfage.2010real"]]
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

rm(list=ls()[ls() %in% c("SFage2010.px", "SFage2010.ls")]) #keeping only raster stack
gc()



#


# scenario 2020 -- avoid deforestation
SFage2010.recovery10 <- calc(SFage2010, fun=function(x){ifelse(is.na(x), x, ifelse(x==0, x, x+10))})
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

rm(list=ls()[ls() %in% c("SFage2010.recovery10.px", "SFage2010.recovery10.ls")]) #keeping only raster stack
gc()



#


# scenario 2020
SFage2020 <- pgm.sfage[["pgm.sfage.2020real"]]
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
writeRaster(SFage2020.px, "rasters/PGM/2020_real/SFagepx.tif", format="GTiff", overwrite=T)
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

rm(list=ls()[ls() %in% c("SFage2020.px", "SFage2020.ls")]) #keeping only raster stack
gc()



#


# scenario -- restoration and avoid
candidate.areas.final.age <- calc(candidate.areas.final, fun=function(x){ifelse(x==1, x+9, x)})
#plot(candidate.areas.final.age)

SFAge2010.restore10 <- sum(SFage2010.recovery10, candidate.areas.final.age, na.rm = T)
values(SFAge2010.restore10)[values(SFAge2010.restore10)>=35]<-35
#plot(SFAge2010.restore10)

# mean sf cover in pixel scale (150m)
SFAge2010.restore10.px <- focal(SFAge2010.restore10, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#SFAge2010.restore10.px
#anyNA(SFAge2010.restore10.px[])
#plot(SFAge2010.restore10.px)

names(SFAge2010.restore10.px)<-"SFagepx"
SFAge2010.restore10.px[is.nan(SFAge2010.restore10.px)] <- 0
SFAge2010.restore10.px <- mask(SFAge2010.restore10.px, pgm.shp)

#saving
writeRaster(SFAge2010.restore10.px, "rasters/PGM/2020_restor_n_avoid/SFagepx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (1050m)
SFAge2010.restore10.ls <- focal(SFAge2010.restore10, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SFAge2010.restore10.ls
#anyNA(SFAge2010.restore10.ls[])
#plot(SFAge2010.restore10.ls)

names(SFAge2010.restore10.ls)<-"SFagels"
SFAge2010.restore10.ls[is.nan(SFAge2010.restore10.ls)] <- 0
SFAge2010.restore10.ls <- mask(SFAge2010.restore10.ls, pgm.shp)

#saving
writeRaster(SFAge2010.restore10.ls, "rasters/PGM/2020_restor_n_avoid/SFagels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2010.restore10.px", "SFAge2010.restore10.ls")]) #keeping only raster stack
gc()



#


# scenario -- restoration without avoid
SFAge2020.restore10 <- sum(SFage2020, candidate.areas.final.age, na.rm = T)
values(SFAge2010.restore10)[values(SFAge2010.restore10)>=35]<-35
#plot(SFAge2020.restore10)

# mean sf cover in pixel scale (150m)
SFAge2020.restore10.px <- focal(SFAge2020.restore10, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#SFAge2020.restore10.px
#anyNA(SFAge2020.restore10.px[])
#plot(SFAge2020.restore10.px)

names(SFAge2020.restore10.px)<-"SFagepx"
SFAge2020.restore10.px[is.nan(SFAge2020.restore10.px)] <- 0
SFAge2020.restore10.px <- mask(SFAge2020.restore10.px, pgm.shp)

#saving
writeRaster(SFAge2020.restore10.px, "rasters/PGM/2020_restor_wo_avoid/SFagepx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (1050m)
SFAge2020.restore10.ls <- focal(SFAge2020.restore10, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SFAge2020.restore10.ls
#anyNA(SFAge2020.restore10.ls[])
#plot(SFAge2020.restore10.ls)

names(SFAge2020.restore10.ls)<-"SFagels"
SFAge2020.restore10.ls[is.nan(SFAge2020.restore10.ls)] <- 0
SFAge2020.restore10.ls <- mask(SFAge2020.restore10.ls, pgm.shp)

#saving
writeRaster(SFAge2020.restore10.ls, "rasters/PGM/2020_restor_wo_avoid/SFagels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020.restore10.px", "SFAge2020.restore10.px")]) #keeping only raster stack
gc()



#


#######################################################################################################################

# [F3] Forest type 3 or Total forest -- pixel: 5x5 (150m)
# this variable includes forest pixels in LULC raster
# including degraded forest and secondary forest older than 2 years

# scenario 2010
SF2010.young <- SFage2010
SF2010.young[SF2010.young <= 2] <- 0
SF2010.young[SF2010.young > 2] <- 1

TF2010 <- sum(UPF2010, DPF2010, SF2010.young, na.rm = T)
TF2010[TF2010>1] <- 1
##cheking
#TF2010
#plot(TF2010)

# mean upf cover in pixel scale (150m)
TF2010.px <- focal(TF2010, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#TF2010.px
#anyNA(TF2010.px[])
#plot(TF2010.px)

names(TF2010.px)<-"TFpx"
TF2010.px[is.nan(TF2010.px)] <- 0
TF2010.px <- mask(TF2010.px, pgm.shp)

#saving
writeRaster(TF2010.px, "rasters/PGM/2010_real/TFpx.tif", format="GTiff", overwrite=T)
writeRaster(TF2010.px, "rasters/PGM/2020_avoidboth/TFpx.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TF2010.px")]) #keeping only raster stack
gc()



#


# scenario 2020
SF2020.young <- SFage2020
SF2020.young[SF2020.young <= 2] <- 0
SF2020.young[SF2020.young > 2] <- 1

TF2020 <- sum(UPF2020, DPF2020, SF2020.young, na.rm = T)
TF2020[TF2020>1] <- 1
##cheking
#TF2020
#plot(TF2020)

# mean upf cover in pixel scale (150m)
TF2020.px <- focal(TF2020, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#TF2020.px
#anyNA(TF2020.px[])
#plot(TF2020.px)

names(TF2020.px)<-"TFpx"
TF2020.px[is.nan(TF2020.px)] <- 0
TF2020.px <- mask(TF2020.px, pgm.shp)

#saving
writeRaster(TF2020.px, "rasters/PGM/2020_real/TFpx.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TF2020.px")]) #keeping only raster stack
gc()



#


# scenario avoid degradation
TF.avoiddegrad <- sum(UPF.avoiddegrad, DPF2010, SF2020.young, na.rm = T)
TF.avoiddegrad[TF.avoiddegrad>1] <- 1
##cheking
#TF.avoiddegrad
#plot(TF.avoiddegrad)

# mean upf cover in pixel scale (150m)
TF.avoiddegrad.px <- focal(TF.avoiddegrad, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#TF.avoiddegrad.px
#anyNA(TF.avoiddegrad.px[])
#plot(TF.avoiddegrad.px)

names(TF.avoiddegrad.px)<-"TFpx"
TF.avoiddegrad.px[is.nan(TF.avoiddegrad.px)] <- 0
TF.avoiddegrad.px <- mask(TF.avoiddegrad.px, pgm.shp)

#saving
writeRaster(TF.avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/TFpx.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TF.avoiddegrad.px")]) #keeping only raster stack
gc()



#


# scenario avoid deforestation
TF.avoiddefor <- sum(UPF.avoiddefor, DPF2020, SF2010.young, na.rm = T)
TF.avoiddefor[TF.avoiddefor>1] <- 1
##cheking
#TF.avoiddefor
#plot(TF.avoiddefor)

# mean upf cover in pixel scale (150m)
TF.avoiddefor.px <- focal(TF.avoiddefor, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#TF.avoiddefor.px
#anyNA(TF.avoiddefor.px[])
#plot(TF.avoiddefor.px)

names(TF.avoiddefor.px)<-"TFpx"
TF.avoiddefor.px[is.nan(TF.avoiddefor.px)] <- 0
TF.avoiddefor.px <- mask(TF.avoiddefor.px, pgm.shp)

#saving
writeRaster(TF.avoiddefor.px, "rasters/PGM/2020_avoiddeforest/TFpx.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TF.avoiddefor.px")]) #keeping only raster stack
gc()



#


# scenario restoration without avoid
SFAge2020.restore10.young <- SFAge2020.restore10
SFAge2020.restore10.young[SFAge2020.restore10.young <= 2] <- 0
SFAge2020.restore10.young[SFAge2020.restore10.young > 2] <- 1

TF.restore10.a <- sum(UPF2020, DPF2020, SFAge2020.restore10.young, na.rm = T)
TF.restore10.a[TF.restore10.a>1] <- 1
##cheking
#TF.restore10.a
#plot(TF.restore10.a)

# mean upf cover in pixel scale (150m)
TF.restore10.a.px <- focal(TF.restore10.a, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#TF.restore10.a.px
#anyNA(TF.restore10.a.px[])
#plot(TF.restore10.a.px)

names(TF.restore10.a.px)<-"TFpx"
TF.restore10.a.px[is.nan(TF.restore10.a.px)] <- 0
TF.restore10.a.px <- mask(TF.restore10.a.px, pgm.shp)

#saving
writeRaster(TF.restore10.a.px, "rasters/PGM/2020_restor_wo_avoid/TFpx.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TF.restore10.a.px")]) #keeping only raster stack
gc()



#


# scenario restoration and avoid
SFAge2010.restore10.young <- SFAge2010.restore10
SFAge2010.restore10.young[SFAge2010.restore10.young <= 2] <- 0
SFAge2010.restore10.young[SFAge2010.restore10.young > 2] <- 1

TF.restore10.b <- sum(UPF2020, DPF2020, SFAge2010.restore10.young, na.rm = T)
TF.restore10.b[TF.restore10.b>1] <- 1
##cheking
#TF.restore10.b
#plot(TF.restore10.b)

# mean upf cover in pixel scale (150m)
TF.restore10.b.px <- focal(TF.restore10.b, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#TF.restore10.b.px
#anyNA(TF.restore10.b.px[])
#plot(TF.restore10.b.px)

names(TF.restore10.b.px)<-"TFpx"
TF.restore10.b.px[is.nan(TF.restore10.b.px)] <- 0
TF.restore10.b.px <- mask(TF.restore10.b.px, pgm.shp)

#saving
writeRaster(TF.restore10.b.px, "rasters/PGM/2020_restor_n_avoid/TFpx.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TF.restore10.b.pxv")]) #keeping only raster stack
gc()



#


#######################################################################################################################

# [F1] Forest type 1 or Mature forest -- pixel: 5x5 (150m)
# this variable includes forest pixels in LULC raster
# including degraded forest and secondary forest older than 10 years

# scenario 2010
SF2010.mature <- SFage2010
SF2010.mature[SF2010.mature <= 10] <- 0
SF2010.mature[SF2010.mature > 10] <- 1

MF2010 <- sum(UPF2010, DPF2010, SF2010.mature, na.rm = T)
MF2010[MF2010>1] <- 1
##cheking
#MF2010
#plot(MF2010)

# mean upf cover in pixel scale (150m)
MF2010.px <- focal(MF2010, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#MF2010.px
#anyNA(MF2010.px[])
#plot(MF2010.px)

names(MF2010.px)<-"MFpx"
MF2010.px[is.nan(MF2010.px)] <- 0
MF2010.px <- mask(MF2010.px, pgm.shp)

#saving
writeRaster(MF2010.px, "rasters/PGM/2010_real/MFpx.tif", format="GTiff", overwrite=T)
writeRaster(MF2010.px, "rasters/PGM/2020_avoidboth/MFpx.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("MF2010.px")]) #keeping only raster stack
gc()



#


# scenario 2020
SF2020.mature <- SFage2020
SF2020.mature[SF2020.mature <= 10] <- 0
SF2020.mature[SF2020.mature > 10] <- 1

MF2020 <- sum(UPF2020, DPF2020, SF2020.mature, na.rm = T)
MF2020[MF2020>1] <- 1
##cheking
#MF2020
#plot(MF2020)

# mean upf cover in pixel scale (150m)
MF2020.px <- focal(MF2020, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#MF2020.px
#anyNA(MF2020.px[])
#plot(MF2020.px)

names(MF2020.px)<-"MFpx"
MF2020.px[is.nan(MF2020.px)] <- 0
MF2020.px <- mask(MF2020.px, pgm.shp)

#saving
writeRaster(MF2020.px, "rasters/PGM/2020_real/MFpx.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("MF2020.px")]) #keeping only raster stack
gc()



#


# scenario avoid degradation
MF.avoiddegrad <- sum(UPF.avoiddegrad, DPF2010, SF2020.mature, na.rm = T)
MF.avoiddegrad[MF.avoiddegrad>1] <- 1
##cheking
#MF.avoiddegrad
#plot(MF.avoiddegrad)

# mean upf cover in pixel scale (150m)
MF.avoiddegrad.px <- focal(MF.avoiddegrad, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#MF.avoiddegrad.px
#anyNA(MF.avoiddegrad.px[])
#plot(MF.avoiddegrad.px)

names(MF.avoiddegrad.px)<-"MFpx"
MF.avoiddegrad.px[is.nan(MF.avoiddegrad.px)] <- 0
MF.avoiddegrad.px <- mask(MF.avoiddegrad.px, pgm.shp)

#saving
writeRaster(MF.avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/MFpx.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("MF.avoiddegrad.px")]) #keeping only raster stack
gc()



#


# scenario avoid deforestation
MF.avoiddefor <- sum(UPF.avoiddefor, DPF2020, SF2010.mature, na.rm = T)
MF.avoiddefor[MF.avoiddefor>1] <- 1
##cheking
#MF.avoiddefor
#plot(MF.avoiddefor)

# mean upf cover in pixel scale (150m)
MF.avoiddefor.px <- focal(MF.avoiddefor, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#MF.avoiddefor.px
#anyNA(MF.avoiddefor.px[])
#plot(MF.avoiddefor.px)

names(MF.avoiddefor.px)<-"MFpx"
MF.avoiddefor.px[is.nan(MF.avoiddefor.px)] <- 0
MF.avoiddefor.px <- mask(MF.avoiddefor.px, pgm.shp)

#saving
writeRaster(MF.avoiddefor.px, "rasters/PGM/2020_avoiddeforest/MFpx.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("MF.avoiddefor.px")]) #keeping only raster stack
gc()



#


# scenario restoration without avoid
SFAge2020.restore10.mature <- SFAge2020.restore10
SFAge2020.restore10.mature[SFAge2020.restore10.mature <= 10] <- 0
SFAge2020.restore10.mature[SFAge2020.restore10.mature > 10] <- 1

MF.restore10.a <- sum(UPF2020, DPF2020, SFAge2020.restore10.mature, na.rm = T)
MF.restore10.a[MF.restore10.a>1] <- 1
##cheking
#MF.restore10.a
#plot(MF.restore10.a)

# mean upf cover in pixel scale (150m)
MF.restore10.a.px <- focal(MF.restore10.a, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#MF.restore10.a.px
#anyNA(MF.restore10.a.px[])
#plot(MF.restore10.a.px)

names(MF.restore10.a.px)<-"MFpx"
MF.restore10.a.px[is.nan(MF.restore10.a.px)] <- 0
MF.restore10.a.px <- mask(MF.restore10.a.px, pgm.shp)

#saving
writeRaster(MF.restore10.a.px, "rasters/PGM/2020_restor_wo_avoid/MFpx.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("MF.restore10.a.px")]) #keeping only raster stack
gc()



#


# scenario restoration and avoid
SFAge2010.restore10.mature <- SFAge2010.restore10
SFAge2010.restore10.mature[SFAge2010.restore10.mature <= 10] <- 0
SFAge2010.restore10.mature[SFAge2010.restore10.mature > 10] <- 1

MF.restore10.b <- sum(UPF2020, DPF2020, SFAge2010.restore10.mature, na.rm = T)
MF.restore10.b[MF.restore10.b>1] <- 1
##cheking
#MF.restore10.b
#plot(MF.restore10.b)

# mean upf cover in pixel scale (150m)
MF.restore10.b.px <- focal(MF.restore10.b, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#MF.restore10.b.px
#anyNA(MF.restore10.b.px[])
#plot(MF.restore10.b.px)

names(MF.restore10.b.px)<-"MFpx"
MF.restore10.b.px[is.nan(MF.restore10.b.px)] <- 0
MF.restore10.b.px <- mask(MF.restore10.b.px, pgm.shp)

#saving
writeRaster(MF.restore10.b.px, "rasters/PGM/2020_restor_n_avoid/MFpx.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("MF.restore10.b.px")]) #keeping only raster stack
gc()



#


##############################################################################################################################################################################################################################################

# [edgedist] distance to forest edge
# this variable is the euclidean distance between mature forest 
# to the nearest cell that is not MF

# scenario 2010
pts <- rasterToPoints(MF2010, spatial=TRUE)
core <- pts[pts$layer == "0",]
#cheking
#core
#plot(core, add=T, col="red")

edge.dist.2010 <- rasterDistance(pts, core, reference = MF2010, scale=TRUE)
##cheking
#edge.dist.2010
#anyNA(edge.dist.2010[])
#plot(edge.dist.2010)

names(edge.dist.2010)<-"edgedist"
edge.dist.2010[is.nan(edge.dist.2010)] <- 0
edge.dist.2010 <- mask(edge.dist.2010, pgm.shp)

#saving
writeRaster(edge.dist.2010, "rasters/PGM/2010_real/edgedist.tif", format="GTiff", overwrite=T)
writeRaster(edge.dist.2010, "rasters/PGM/2020_avoidboth/edgedist.tif", format="GTiff", overwrite=T)

gc()



#


# scenario 2020
pts <- rasterToPoints(MF2020, spatial=TRUE)
core <- pts[pts$layer == "0",]
#cheking
#core
#plot(core, add=T, col="red")

edge.dist.2020 <- rasterDistance(pts, core, reference = MF2020, scale=TRUE)
##cheking
#edge.dist.2020
#anyNA(edge.dist.2020[])
#plot(edge.dist.2020)

names(edge.dist.2020)<-"edgedist"
edge.dist.2020[is.nan(edge.dist.2020)] <- 0
edge.dist.2020 <- mask(edge.dist.2020, pgm.shp)

#saving
writeRaster(edge.dist.2020, "rasters/PGM/2020_real/edgedist.tif", format="GTiff", overwrite=T)

gc()



#


# scenario avoid degradation
pts <- rasterToPoints(MF.avoiddegrad, spatial=TRUE)
core <- pts[pts$layer == "0",]
#cheking
#core
#plot(core, add=T, col="red")

edge.dist.avoiddegrad <- rasterDistance(pts, core, reference = MF.avoiddegrad, scale=TRUE)
##cheking
#edge.dist.avoiddegrad
#anyNA(edge.dist.avoiddegrad[])
#plot(edge.dist.avoiddegrad)

names(edge.dist.avoiddegrad)<-"edgedist"
edge.dist.avoiddegrad[is.nan(edge.dist.avoiddegrad)] <- 0
edge.dist.avoiddegrad <- mask(edge.dist.avoiddegrad, pgm.shp)

#saving
writeRaster(edge.dist.avoiddegrad, "rasters/PGM/2020_avoiddegrad/edgedist.tif", format="GTiff", overwrite=T)

gc()



#


# scenario avoid deforestation
pts <- rasterToPoints(MF.avoiddefor, spatial=TRUE)
core <- pts[pts$layer == "0",]
#cheking
#core
#plot(core, add=T, col="red")

edge.dist.avoiddefor <- rasterDistance(pts, core, reference = MF.avoiddefor, scale=TRUE)
##cheking
#edge.dist.avoiddefor
#anyNA(edge.dist.avoiddefor[])
#plot(edge.dist.avoiddefor)

names(edge.dist.avoiddefor)<-"edgedist"
edge.dist.avoiddefor[is.nan(edge.dist.avoiddefor)] <- 0
edge.dist.avoiddefor <- mask(edge.dist.avoiddefor, pgm.shp)

#saving
writeRaster(edge.dist.avoiddefor, "rasters/PGM/2020_avoiddeforest/edgedist.tif", format="GTiff", overwrite=T)

gc()



#


# scenario restoration without avoid
pts <- rasterToPoints(MF.restore10.a, spatial=TRUE)
core <- pts[pts$layer == "0",]
#cheking
#core
#plot(core, add=T, col="red")

edge.dist.restore10.a <- rasterDistance(pts, core, reference = MF.restore10.a, scale=TRUE)
##cheking
#edge.dist.restore10.a
#anyNA(edge.dist.restore10.a[])
#plot(edge.dist.restore10.a)

names(edge.dist.restore10.a)<-"edgedist"
edge.dist.restore10.a[is.nan(edge.dist.restore10.a)] <- 0
edge.dist.restore10.a <- mask(edge.dist.restore10.a, pgm.shp)

#saving
writeRaster(edge.dist.restore10.a, "rasters/PGM/2020_restor_wo_avoid/edgedist.tif", format="GTiff", overwrite=T)

gc()



#


# scenario restoration and avoid
pts <- rasterToPoints(MF.restore10.b, spatial=TRUE)
core <- pts[pts$layer == "0",]
#cheking
#core
#plot(core, add=T, col="red")

edge.dist.restore10.b <- rasterDistance(pts, core, reference = MF.restore10.b, scale=TRUE)
##cheking
#edge.dist.restore10.b
#anyNA(edge.dist.restore10.b[])
#plot(edge.dist.restore10.b)

names(edge.dist.restore10.b)<-"edgedist"
edge.dist.restore10.b[is.nan(edge.dist.restore10.b)] <- 0
edge.dist.restore10.b <- mask(edge.dist.restore10.b, pgm.shp)

#saving
writeRaster(edge.dist.restore10.b, "rasters/PGM/2020_restor_n_avoid/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("pts", "core")]) #keeping only raster stack
gc()



#


##############################################################################################################################################################################################################################################

# [edge] forest edge -- pixel: 5x5 (150m); and landscape: 35x35 (1050m)
# this variable is the mean of edge area based on mature forest

# scenario 2010
# marking edge and core areas
edge2010 <- focal(MF2010, matrix(1,ncol=3,nrow=3), fun = mean, na.rm = T, pad = T)
edge2010 <- edge2010 + MF2010
edge2010[edge2010 == 2] <- 0                  # core area
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
edge2010.px <- mask(edge2010.px, pgm.shp)

#saving
writeRaster(edge2010.px, "rasters/PGM/2010_real/edgepx.tif", format="GTiff", overwrite=T)
writeRaster(edge2010.px, "rasters/PGM/2020_avoidboth/edgepx.tif", format="GTiff", overwrite=T)
#

# mean sf cover in landscape scale (1050m)
edge2010.ls <- focal(edge2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#edge2010.ls
#anyNA(edge2010.ls[])
#plot(edge2010.ls)

names(edge2010.ls)<-"edgels"
edge2010.ls[is.nan(edge2010.ls)] <- 0
edge2010.ls <- mask(edge2010.ls, pgm.shp)

#saving
writeRaster(edge2010.px, "rasters/PGM/2010_real/edgels.tif", format="GTiff", overwrite=T)
writeRaster(edge2010.px, "rasters/PGM/2020_avoidboth/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2010.px", "edge2010.ls")]) #keeping only raster stack
gc()



#


# scenario 2020
# marking edge and core areas
edge2020 <- focal(MF2020, matrix(1,ncol=3,nrow=3), fun = mean, na.rm = T, pad = T)
edge2020 <- edge2020 + MF2020
edge2020[edge2020 == 2] <- 0                  # core area
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
edge2020.px <- mask(edge2020.px, pgm.shp)

#saving
writeRaster(edge2020.px, "rasters/PGM/2020_real/edgepx.tif", format="GTiff", overwrite=T)
#

# mean sf cover in landscape scale (1050m)
edge2020.ls <- focal(edge2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#edge2020.ls
#anyNA(edge2020.ls[])
#plot(edge2020.ls)

names(edge2020.ls)<-"edgels"
edge2020.ls[is.nan(edge2020.ls)] <- 0
edge2020.ls <- mask(edge2020.ls, pgm.shp)

#saving
writeRaster(edge2020.ls, "rasters/PGM/2020_real/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020.px", "edge2020.ls")]) #keeping only raster stack
gc()



#


# scenario avoid degradation
# marking edge and core areas
edge.avoiddegrad <- focal(MF.avoiddegrad, matrix(1,ncol=3,nrow=3), fun = mean, na.rm = T, pad = T)
edge.avoiddegrad <- edge.avoiddegrad + MF.avoiddegrad
edge.avoiddegrad[edge.avoiddegrad == 2] <- 0                  # core area
edge.avoiddegrad[edge.avoiddegrad > 0 & edge.avoiddegrad < 2] <- 1    # edge
#cheking
#edge.avoiddegrad
#unique(edge.avoiddegrad[])
#plot(edge.avoiddegrad)

# mean edge in pixel scale (150m)
edge.avoiddegrad.px <- focal(edge.avoiddegrad, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#edge.avoiddegrad.px
#anyNA(edge.avoiddegrad.px[])
#plot(edge.avoiddegrad.px)

names(edge.avoiddegrad.px)<-"edgepx"
edge.avoiddegrad.px[is.nan(edge.avoiddegrad.px)] <- 0
edge.avoiddegrad.px <- mask(edge.avoiddegrad.px, pgm.shp)

#saving
writeRaster(edge.avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/edgepx.tif", format="GTiff", overwrite=T)
#

# mean sf cover in landscape scale (1050m)
edge.avoiddegrad.ls <- focal(edge.avoiddegrad, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#edge.avoiddegrad.ls
#anyNA(edge.avoiddegrad.ls[])
#plot(edge.avoiddegrad.ls)

names(edge.avoiddegrad.ls)<-"edgels"
edge.avoiddegrad.ls[is.nan(edge.avoiddegrad.ls)] <- 0
edge.avoiddegrad.ls <- mask(edge.avoiddegrad.ls, pgm.shp)

#saving
writeRaster(edge.avoiddegrad.ls, "rasters/PGM/2020_avoiddegrad/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge.avoiddegrad.px", "edge.avoiddegrad.ls")]) #keeping only raster stack
gc()



#


# scenario avoid deforestation
# marking edge and core areas
edge.avoiddefor <- focal(MF.avoiddefor, matrix(1,ncol=3,nrow=3), fun = mean, na.rm = T, pad = T)
edge.avoiddefor <- edge.avoiddefor + MF.avoiddefor
edge.avoiddefor[edge.avoiddefor == 2] <- 0                  # core area
edge.avoiddefor[edge.avoiddefor > 0 & edge.avoiddefor < 2] <- 1    # edge
#cheking
#edge.avoiddefor
#unique(edge.avoiddefor[])
#plot(edge.avoiddefor)

# mean edge in pixel scale (150m)
edge.avoiddefor.px <- focal(edge.avoiddefor, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#edge.avoiddefor.px
#anyNA(edge.avoiddefor.px[])
#plot(edge.avoiddefor.px)

names(edge.avoiddefor.px)<-"edgepx"
edge.avoiddefor.px[is.nan(edge.avoiddefor.px)] <- 0
edge.avoiddefor.px <- mask(edge.avoiddefor.px, pgm.shp)

#saving
writeRaster(edge.avoiddefor.px, "rasters/PGM/2020_avoiddeforest/edgepx.tif", format="GTiff", overwrite=T)
#

# mean sf cover in landscape scale (1050m)
edge.avoiddefor.ls <- focal(edge.avoiddefor, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#edge.avoiddefor.ls
#anyNA(edge.avoiddefor.ls[])
#plot(edge.avoiddefor.ls)

names(edge.avoiddefor.ls)<-"edgels"
edge.avoiddefor.ls[is.nan(edge.avoiddefor.ls)] <- 0
edge.avoiddefor.ls <- mask(edge.avoiddefor.ls, pgm.shp)

#saving
writeRaster(edge.avoiddefor.ls, "rasters/PGM/2020_avoiddeforest/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge.avoiddefor.px", "edge.avoiddefor.ls")]) #keeping only raster stack
gc()



#


# scenario restoration without avoid
# marking edge and core areas
edge.restore10.a <- focal(MF.restore10.a, matrix(1,ncol=3,nrow=3), fun = mean, na.rm = T, pad = T)
edge.restore10.a <- edge.restore10.a + MF.restore10.a
edge.restore10.a[edge.restore10.a == 2] <- 0                  # core area
edge.restore10.a[edge.restore10.a > 0 & edge.restore10.a < 2] <- 1    # edge
#cheking
#edge.restore10.a
#unique(edge.restore10.a[])
#plot(edge.restore10.a)

# mean edge in pixel scale (150m)
edge.restore10.a.px <- focal(edge.restore10.a, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#edge.restore10.a.px
#anyNA(edge.restore10.a.px[])
#plot(edge.restore10.a.px)

names(edge.restore10.a.px)<-"edgepx"
edge.restore10.a.px[is.nan(edge.restore10.a.px)] <- 0
edge.restore10.a.px <- mask(edge.restore10.a.px, pgm.shp)

#saving
writeRaster(edge.restore10.a.px, "rasters/PGM/2020_restor_wo_avoid/edgepx.tif", format="GTiff", overwrite=T)
#

# mean sf cover in landscape scale (1050m)
edge.restore10.a.ls <- focal(edge.restore10.a, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#edge.restore10.a.ls
#anyNA(edge.restore10.a.ls[])
#plot(edge.restore10.a.ls)

names(edge.restore10.a.ls)<-"edgels"
edge.restore10.a.ls[is.nan(edge.restore10.a.ls)] <- 0
edge.restore10.a.ls <- mask(edge.restore10.a.ls, pgm.shp)

#saving
writeRaster(edge.restore10.a.ls, "rasters/PGM/2020_restor_wo_avoid/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge.restore10.a.px", "edge.restore10.a.ls")]) #keeping only raster stack
gc()



#


# scenario restoration and avoid
# marking edge and core areas
edge.restore10.b <- focal(MF.restore10.b, matrix(1,ncol=3,nrow=3), fun = mean, na.rm = T, pad = T)
edge.restore10.b <- edge.restore10.b + MF.restore10.b
edge.restore10.b[edge.restore10.b == 2] <- 0                  # core area
edge.restore10.b[edge.restore10.b > 0 & edge.restore10.b < 2] <- 1    # edge
#cheking
#edge.restore10.b
#unique(edge.restore10.b[])
#plot(edge.restore10.b)

# mean edge in pixel scale (150m)
edge.restore10.b.px <- focal(edge.restore10.b, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#edge.restore10.b.px
#anyNA(edge.restore10.b.px[])
#plot(edge.restore10.b.px)

names(edge.restore10.b.px)<-"edgepx"
edge.restore10.b.px[is.nan(edge.restore10.b.px)] <- 0
edge.restore10.b.px <- mask(edge.restore10.b.px, pgm.shp)

#saving
writeRaster(edge.restore10.b.px, "rasters/PGM/2020_restor_n_avoid/edgepx.tif", format="GTiff", overwrite=T)
#

# mean sf cover in landscape scale (1050m)
edge.restore10.b.ls <- focal(edge.restore10.b, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#edge.restore10.b.ls
#anyNA(edge.restore10.b.ls[])
#plot(edge.restore10.b.ls)

names(edge.restore10.b.ls)<-"edgels"
edge.restore10.b.ls[is.nan(edge.restore10.b.ls)] <- 0
edge.restore10.b.ls <- mask(edge.restore10.b.ls, pgm.shp)

#saving
writeRaster(edge.restore10.b.ls, "rasters/PGM/2020_restor_n_avoid/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge.restore10.b.px", "edge.restore10.b.ls")]) #keeping only raster stack
gc()



#


##############################################################################################################################################################################################################################################

# [meantemp] annual average temperature from nasa earth observation

# download and save global data
#urls <- c("https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755469&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-01"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755470&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-02"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755471&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-03"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755472&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-04"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755473&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-05"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755474&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-06"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755475&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-07"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755476&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-08"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755477&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-09"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755478&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-10"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755479&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-11"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755480&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-12"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1784090&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-01"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1785058&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-02"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1785890&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-03"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1786979&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-04"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1794500&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-05"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1795357&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-06"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1796358&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-07"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1797155&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-08"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1799174&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-09"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1799941&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-10"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1800680&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-11"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1801650&cs=rgb&format=TIFF&width=3600&height=1800") #"2020-12"
#
#dir.create("rasters/PGM/input/climate")


# scenario 2010
temp.list <- list.files("rasters/PGM/input/climate", "LSTD", full.names = T, recursive = T)

temp2010.list <- grep("2010", temp.list, value = T)
temp2010 <- stack(temp2010.list)
#plot(temp2010)

meantemp2010 <- mean(temp2010, na.rm=T)
pgm.meantemp2010 <- crop(meantemp2010, extent(pgm.lulc.2010.forest.class))
#plot(pgm.meantemp2010)
#plot(pgm.shp, add=T)

pgm.meantemp2010 <- resample(pgm.meantemp2010, pgm.lulc.2010.forest.class, method='bilinear')
pgm.meantemp2010 <- mask(pgm.meantemp2010, pgm.shp)
#plot(pgm.meantemp2010)

#saving
writeRaster(pgm.meantemp2010, "rasters/PGM/2010_real/meantemps.tif", format="GTiff", overwrite=T)
#

# scenario 2020
temp2020.list <- grep("2020", temp.list, value = T)
temp2020 <- stack(temp2020.list)
#plot(temp2020)

meantemp2020 <- mean(temp2020, na.rm=T)
pgm.meantemp2020 <- crop(meantemp2020, extent(pgm.lulc.2020.forest.class))
#plot(pgm.meantemp2020)
#plot(pgm.shp, add=T)

pgm.meantemp2020 <- resample(pgm.meantemp2020, pgm.lulc.2020.forest.class, method='bilinear')
pgm.meantemp2020 <- mask(pgm.meantemp2020, pgm.shp)
#plot(pgm.meantemp2020)

#saving
writeRaster(pgm.meantemp2020, "rasters/PGM/2020_real/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meantemp2020, "rasters/PGM/2020_avoiddeforest/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meantemp2020, "rasters/PGM/2020_avoiddegrad/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meantemp2020, "rasters/PGM/2020_avoidboth/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meantemp2020, "rasters/PGM/2020_restor_wo_avoid/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meantemp2020, "rasters/PGM/2020_restor_n_avoid/meantemps.tif", format="GTiff", overwrite=T)
#

rm(list=ls()[ls() %in% c("temp.list", "temp2010.list", "meantemp2010", "temp2020.list", "meantemp2020")]) #keeping only raster stack
gc()



#


##############################################################################################################################################################################################################################################

# [meanprecip] annual average precipitation from nasa earth observation

# download and save global data
#urls <- c("https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843747&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-01"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843749&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-02"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843751&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-03"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843755&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-04"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843735&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-05"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843745&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-06"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843753&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-07"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843759&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-08"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843761&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-09"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843763&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-10"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843765&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-11"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843769&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-12"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843983&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-01"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843985&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-02"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843987&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-03"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843989&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-04"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843991&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-05"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843993&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-06"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843995&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-07"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843997&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-08"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843999&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-09"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1844001&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-10"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1844003&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-11"
#          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1844005&cs=rgb&format=TIFF&width=3600&height=1800") #"2020-12"
#


# scenario 2010
precip.list <- list.files("rasters/PGM/input/climate", "GPM", full.names = T, recursive = T)

precip2010.list <- grep("2010", precip.list, value = T)
precip2010 <- stack(precip2010.list)
#plot(precip2010)
rm(precip2010.list)

meanprecip2010 <- mean(precip2010, na.rm=T)
pgm.meanprecip2010 <- crop(meanprecip2010, extent(pgm.lulc.2010.forest.class))
#plot(pgm.meanprecip2010)
#plot(pgm.shp, add=T)

pgm.meanprecip2010 <- resample(pgm.meanprecip2010, pgm.lulc.2010.forest.class, method='bilinear')
pgm.meanprecip2010 <- mask(pgm.meanprecip2010, pgm.shp)
#plot(pgm.meanprecip2010)

#saving
writeRaster(pgm.meanprecip2010, "rasters/PGM/2010_real/meanprecips.tif", format="GTiff", overwrite=T)
#

# scenario 2020
precip2020.list <- grep("2020", precip.list, value = T)
precip2020 <- stack(precip2020.list)
#plot(precip2020)


meanprecip2020 <- mean(precip2020, na.rm=T)
pgm.meanprecip2020 <- crop(meanprecip2020, extent(pgm.lulc.2020.forest.class))
#plot(pgm.meanprecip2020)
#plot(pgm.shp, add=T)

pgm.meanprecip2020 <- resample(pgm.meanprecip2020, pgm.lulc.2020.forest.class, method='bilinear')
pgm.meanprecip2020 <- mask(pgm.meanprecip2020, pgm.shp)
#plot(pgm.meanprecip2020)

#saving
writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_real/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_avoiddeforest/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_avoiddegrad/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_avoidboth/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_restor_wo_avoid/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_restor_n_avoid/meanprecips.tif", format="GTiff", overwrite=T)
#

#optional
#save.image("~/conserv_opportunities_jamesthomson/github_repo/pgm_environment.RData")
rm(list=ls())
gc()



#


##############################################################################################################################################################################################################################################

#### detecting multicollinearity between exploratory variables ####
env.explanatory.var.list <- list.files("rasters/PGM/2010_real", pattern = ".tif", full.names = T, recursive = T)

env.explanatory.var <- stack(env.explanatory.var.list)
names(env.explanatory.var) <- unlist(strsplit(env.explanatory.var.list, "/|.tif"))[seq(4,80,4)]
##cheking
#env.explanatory.var
#plot(env.explanatory.var[[1:10]], nc=2)
#plot(env.explanatory.var[[11:20]], nc=2)

# visual inspection of aggregation using removeCollinearity() function from virtualspecies package
correlated.var <- removeCollinearity(env.explanatory.var, multicollinearity.cutoff = 0.7, sample.points = T, nb.points = 999999, method = "pearson", plot = T)
##cheking
#correlated.var

# variation inflation factor
inflated.var <- vifcor(env.explanatory.var, th = 0.7, maxobservations = 999999)
##cheking
#inflated.var@results
#inflated.var@excluded

sel.var.df <- data.frame(rbind(cbind(VAR=inflated.var@results$Variables, VIF=inflated.var@results$VIF),
                               cbind(VAR=inflated.var@excluded, VIF=NA)))

write.csv(sel.var.df, paste0("rasters/PGM/selected_environmental_explanatory_variables_byVIF.csv", sep=""), row.names = F)

##cheking
#env.explanatory.var <- env.explanatory.var[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
#plot(env.explanatory.var)


#

rm(list=ls()) #keeping only raster stack
gc()



#


