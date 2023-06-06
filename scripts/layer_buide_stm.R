
#' @title Cost-benefit of conservation actions in Amazon
#' @description script to build exploratory variables from 
#' land use - land cover, secondary forest, edge, 
#' degradation (loggin and fire), temperature, precipitation, 
#' elevation and distances to road and water body
#' in Santarem, Belterra and Mojui dos Campos municipalities - PA;
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
dir.create("rasters/STM/2010_real", recursive = T)
dir.create("rasters/STM/2020_real", recursive = T)
dir.create("rasters/STM/2020_avoiddeforest", recursive = T)
dir.create("rasters/STM/2020_avoiddegrad", recursive = T)
dir.create("rasters/STM/2020_avoidboth", recursive = T)
dir.create("rasters/STM/2020_restor_wo_avoid", recursive = T)
dir.create("rasters/STM/2020_restor_n_avoid_deforest", recursive = T)
dir.create("rasters/STM/2020_restor_n_avoid_both", recursive = T)



#### importing input rasters ####

# shapefile paragominas
stm.shp <- readOGR(dsn = "shapes", layer = "Santarem_Belterra_Mask_R3")
stm.shp <- spTransform(stm.shp, CRS("+proj=longlat +datum=WGS84 +no_defs"))
#
#



# land use land cover from mapbiomas collection 7 [2010 and 2020]
# [STM] paragominas
stm.lulc <- stack(c("rasters/STM/input/stm-2010-lulc-mapbiomas-brazil-collection-70.tif",
                    "rasters/STM/input/stm-2020-lulc-mapbiomas-brazil-collection-70.tif"))
names(stm.lulc) <- c("stm.lulc.2010real", "stm.lulc.2020real")
#checking
#stm.lulc
#plot(stm.lulc)
#sort(unique(values(stm.lulc[["stm.lulc.2010real"]])))

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

#stm.lulc.2010.df <- as.data.frame(stm.lulc[["stm.lulc.2010real"]], xy = TRUE)
#breakpoints <- sort(unique(stm.lulc.2010.df$stm.lulc.2010real))
#labels.legend <- c("Non Observed", "Forest Formation", "Savanna Formation", "Wetland", "Grassland",
#                   "Pasture", "Urban Area", "Water", "Soybean", "Other Temporary Crops")
#mapbiomas.legend <- c("#ffffff", "#006400", "#32CD32", "#45C2A5", "#B8AF4F", "#B8AF4F", "#af2a2a",
#                      "#0000FF", "#c59ff4", "#e787f8")
#
#ggplot() +
#  geom_raster(data = stm.lulc.2010.df , aes(x = x, y = y, fill = factor(stm.lulc.2010real))) + 
#  scale_fill_manual(breaks = breakpoints, values = mapbiomas.legend, labels = labels.legend, name = "LULC Classes") +
#  theme_void()


# isolating forest class pixels
stm.lulc.2010.forest.class <- stm.lulc[["stm.lulc.2010real"]]
stm.lulc.2010.forest.class[stm.lulc.2010.forest.class==3] <- 1
stm.lulc.2010.forest.class[stm.lulc.2010.forest.class>1] <- 0

stm.lulc.2020.forest.class <- stm.lulc[["stm.lulc.2020real"]]
stm.lulc.2020.forest.class[stm.lulc.2020.forest.class==3] <- 1
stm.lulc.2020.forest.class[stm.lulc.2020.forest.class>1] <- 0

#
#



# distances to road and rivers, and elevation

dist.road <- raster("rasters/STM/input/dist_road_stm.tif")
dist.road <- projectRaster(dist.road, crs = "+proj=longlat +datum=WGS84 +no_defs")
dist.road <- resample(dist.road, stm.lulc.2010.forest.class, method='bilinear')
values(dist.road)[is.na(values(dist.road))] = 0
dist.road <- mask(dist.road, stm.shp)

#saving
writeRaster(dist.road, "rasters/STM/2010_real/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/STM/2020_real/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/STM/2020_avoiddeforest/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/STM/2020_avoiddegrad/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/STM/2020_avoidboth/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/STM/2020_restor_wo_avoid/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/STM/2020_restor_n_avoid_both/distroad.tif", format="GTiff", overwrite=T)

#

dist.river <- raster("rasters/STM/input/dist_river_stm.tif")
dist.river <- projectRaster(dist.river, crs = "+proj=longlat +datum=WGS84 +no_defs")
dist.river <- resample(dist.river, stm.lulc.2010.forest.class, method='bilinear')
values(dist.river)[is.na(values(dist.river))] = 0
dist.river <- mask(dist.river, stm.shp)

#saving
writeRaster(dist.river, "rasters/STM/2010_real/distriver.tif", format="GTiff", overwrite=T)
writeRaster(dist.river, "rasters/STM/2020_real/distriver.tif", format="GTiff", overwrite=T)
writeRaster(dist.river, "rasters/STM/2020_avoiddeforest/distriver.tif", format="GTiff", overwrite=T)
writeRaster(dist.river, "rasters/STM/2020_avoiddegrad/distriver.tif", format="GTiff", overwrite=T)
writeRaster(dist.river, "rasters/STM/2020_avoidboth/distriver.tif", format="GTiff", overwrite=T)
writeRaster(dist.river, "rasters/STM/2020_restor_wo_avoid/distriver.tif", format="GTiff", overwrite=T)
writeRaster(dist.river, "rasters/STM/2020_restor_n_avoid_both/distriver.tif", format="GTiff", overwrite=T)

#

elevation <- raster("rasters/STM/input/elevation_stm.tif")
elevation <- projectRaster(elevation, crs = "+proj=longlat +datum=WGS84 +no_defs")
elevation <- resample(elevation, stm.lulc.2010.forest.class, method='bilinear')
values(elevation)[is.na(values(elevation))] = 0
elevation <- mask(elevation, stm.shp)

#saving
writeRaster(elevation, "rasters/STM/2010_real/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/STM/2020_real/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/STM/2020_avoiddeforest/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/STM/2020_avoiddegrad/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/STM/2020_avoidboth/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/STM/2020_restor_wo_avoid/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/STM/2020_restor_n_avoid_both/elevation.tif", format="GTiff", overwrite=T)

#

#import rural properties shapefiles and data from SISCAR
#https://www.car.gov.br/publico/municipios/downloads
stm.car <- readOGR(dsn = "rasters/STM/input/SHAPE_1506807_Santarem", layer = "AREA_IMOVEL")
#head(stm.car@data)
#plot(stm.car)

btr.car <- readOGR(dsn = "rasters/STM/input/SHAPE_1501451_Belterra", layer = "AREA_IMOVEL")
#head(btr.car@data)
#plot(btr.car, add=T)

mjc.car <- readOGR(dsn = "rasters/STM/input/SHAPE_1504752_Mojui_dos_Campos", layer = "AREA_IMOVEL")
#head(mjc.car@data)
#plot(mjc.car, add=T)

stm.car <- rbind(stm.car,btr.car,mjc.car)
stm.car <- crop(stm.car, extent(stm.shp))

#excluding protected areas
stm.car <- stm.car[-which(stm.car@data$TIPO_IMOVE=="PCT"),]
#head(stm.car@data)
#plot(stm.car, add=T)
stm.car.raster <- rasterize(stm.car, stm.lulc[[1]], field = "NUM_AREA", fun = "mean")

#1st edition -- replacing 0 and water to NA
stm.lulc.mask <- stm.lulc[[1]]
stm.lulc.mask[stm.lulc.mask == 0] <- NA
stm.lulc.mask[stm.lulc.mask == 33] <- NA

stm.car.raster <- mask(stm.car.raster, stm.lulc.mask)

#2nd edition -- excluding areas on islands
pts <- st_as_sf(gBuffer(SpatialPoints(coords = cbind(c(-54.7,-54.47),c(-2.22,-2.43))), width = .2))

mask <- st_bbox(stm.car.raster) %>% # take extent of your raster... 
  st_as_sfc() %>% # make it a sf object
  st_set_crs(st_crs(pts)) %>% # in CRS of your polygon 
  st_difference(pts) %>% # intersect with the polygon object
  st_as_sf() # interpret as sf (and not sfc) object

stm.car.raster <- mask(stm.car.raster, mask)


#saving
writeRaster(stm.car.raster, "rasters/STM/2010_real/property.tif", format="GTiff", overwrite=T)
writeRaster(stm.car.raster, "rasters/STM/2020_real/property.tif", format="GTiff", overwrite=T)
writeRaster(stm.car.raster, "rasters/STM/2020_avoiddeforest/property.tif", format="GTiff", overwrite=T)
writeRaster(stm.car.raster, "rasters/STM/2020_avoiddegrad/property.tif", format="GTiff", overwrite=T)
writeRaster(stm.car.raster, "rasters/STM/2020_avoidboth/property.tif", format="GTiff", overwrite=T)
writeRaster(stm.car.raster, "rasters/STM/2020_restor_wo_avoid/property.tif", format="GTiff", overwrite=T)
writeRaster(stm.car.raster, "rasters/STM/2020_restor_n_avoid_both/property.tif", format="GTiff", overwrite=T)

#




#### candidate areas for restoration scenarios ####

#isolating deforestation class pixels (crops, pasture)

deforestation.class.list <- c(15,39,41,48)

candidate.areas.total <- stm.lulc[["stm.lulc.2010real"]]

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
#candidate.areas.water <- mask(candidate.areas.water, stm.shp)
#plot(candidate.areas.water)



slope <- terrain(elevation, opt = 'slope', unit = 'degrees', neighbors=8)
values(slope)[values(slope) < 45] = NA
values(slope)[values(slope) >= 45] = 1
#plot(slope)

candidate.areas.slope <- candidate.areas.total
candidate.areas.slope <- mask(candidate.areas.slope, slope)
values(candidate.areas.slope)[is.na(values(candidate.areas.slope))] = 0
#candidate.areas.slope <- mask(candidate.areas.slope, stm.shp)
#plot(candidate.areas.slope)


forest.class <- stm.lulc[["stm.lulc.2010real"]]
forest.class[forest.class==0] <- NA
forest.class[forest.class==3] <- 1
forest.class[forest.class>1] <- 0

deforest.pts <- rasterToPoints(forest.class, spatial=TRUE)
deforest.core <- deforest.pts[deforest.pts$stm.lulc.2010real == "1",]

deforest.dist <- rasterDistance(deforest.pts, deforest.core, reference = candidate.areas.total, scale=TRUE)
values(deforest.dist)[values(deforest.dist) == 0] = NA
values(deforest.dist)[values(deforest.dist) > 0.15] = NA
values(deforest.dist)[values(deforest.dist) <= 0.15] = 1
#plot(deforest.dist)

candidate.areas.forest <- candidate.areas.total
candidate.areas.forest <- mask(candidate.areas.forest, deforest.dist)
values(candidate.areas.forest)[is.na(values(candidate.areas.forest))] = 0
#candidate.areas.forest <- mask(candidate.areas.forest, stm.shp)
#plot(candidate.areas.forest)


candidate.areas.final <- sum(candidate.areas.water,candidate.areas.slope,candidate.areas.forest)
values(candidate.areas.final)[values(candidate.areas.final) >= 1] = 1
#plot(candidate.areas.final)

#select rural properties with less than 50% of forest cover in 2010



#calculate forest cover in each property and select properties with <50%

stm.car@data$FOREST_COVER <- NA
j=nrow(stm.car@data)
for (i in stm.car$COD_IMOVEL) {
  
  rural.property <- stm.car[stm.car$COD_IMOVEL==i,]
  forest.cover <- crop(forest.class, extent(rural.property))
  forest.cover <- mask(forest.cover, rural.property)
  stm.car[stm.car$COD_IMOVEL==i,"FOREST_COVER"] <- tapply(area(forest.cover), forest.cover[], sum)[2]*100
  j=j-1
  cat("\n>", j, "out of", nrow(stm.car@data), "properties left<\n")
  
}


stm.car@data$FOREST_COVER_PP <- ceiling((stm.car@data$FOREST_COVER/stm.car@data$NUM_AREA)*100)

#select properties
stm.car.restoration.candidates <- stm.car[!is.na(stm.car$FOREST_COVER_PP),]
stm.car.restoration.candidates <- stm.car.restoration.candidates[stm.car.restoration.candidates$FOREST_COVER_PP <= 50,]
#head(stm.car.restoration.candidates@data)
#nrow(stm.car.restoration.candidates@data)


#filter candidate areas for restoration in properties with less than 50% forest cover
candidate.areas.final.copy <- candidate.areas.final
candidate.areas.final <- mask(candidate.areas.final, stm.car.restoration.candidates)
values(candidate.areas.final)[is.na(values(candidate.areas.final))] = 0
#candidate.areas.final <- mask(candidate.areas.final, stm.shp)
#plot(candidate.areas.final)
#plot(stm.car.restoration.candidates, add=T)


#forest cover increment -- adding the candidate areas to forest cover
stm.car.restoration.candidates@data$FOREST_COVER_INCREMENT <- NA
j=nrow(stm.car.restoration.candidates@data)
#i="PA-1505502-39CCE4418D2D487F9AC0FD3045A374CF"
for (i in stm.car.restoration.candidates$COD_IMOVEL) {
  
  rural.property <- stm.car.restoration.candidates[stm.car.restoration.candidates$COD_IMOVEL==i,]
  
  forest.cover <- crop(forest.class, extent(rural.property))
  forest.cover <- mask(forest.cover, rural.property)
  
  restored.cover <- crop(candidate.areas.final, extent(rural.property))
  restored.cover <- mask(restored.cover, rural.property)
  
  forest.cover.increment <- sum(forest.cover, restored.cover)
  
  stm.car.restoration.candidates[stm.car.restoration.candidates$COD_IMOVEL==i,"FOREST_COVER_INCREMENT"] <- tapply(area(forest.cover.increment), forest.cover.increment[], sum)[2]*100
  j=j-1
  cat("\n>", j, "out of", nrow(stm.car.restoration.candidates@data), "properties left<\n")
  
}

stm.car.restoration.candidates@data$FOREST_COVER_INCREMENT_PP <- ceiling((stm.car.restoration.candidates@data$FOREST_COVER_INCREMENT/stm.car.restoration.candidates@data$NUM_AREA)*100)


#select properties with more than 80% forest cover
stm.car.restoration.candidates <- stm.car.restoration.candidates[!is.na(stm.car.restoration.candidates$FOREST_COVER_INCREMENT_PP),]
stm.car.restoration.candidates.m80 <- stm.car.restoration.candidates[stm.car.restoration.candidates$FOREST_COVER_INCREMENT_PP > 80,]
#head(stm.car.restoration.candidates.m80@data)
#nrow(stm.car.restoration.candidates.m80@data)



candidate.areas.final.copy <- candidate.areas.final
j=nrow(stm.car.restoration.candidates.m80@data)
#i="PA-1505502-75B3625903EB4F6A9C36FF2B97B67282"
for (i in stm.car.restoration.candidates.m80$COD_IMOVEL) {
  
  rural.property <- stm.car.restoration.candidates.m80[stm.car.restoration.candidates.m80$COD_IMOVEL==i,]
  
  forest.cover <- crop(forest.class, extent(rural.property))
  forest.cover <- mask(forest.cover, rural.property)
  
  restored.cover <- crop(candidate.areas.final, extent(rural.property))
  restored.cover <- mask(restored.cover, rural.property)
  
  
  while (stm.car.restoration.candidates.m80@data[stm.car.restoration.candidates.m80@data$COD_IMOVEL==i,"FOREST_COVER_INCREMENT_PP"] > 80) {
    
    restored.cover[restored.cover[]==1] <- sample(c(1,0), size = length(restored.cover[restored.cover[]==1]), replace = T, prob = c(0.9,0.1))
    forest.cover.increment <- sum(forest.cover, restored.cover)
    stm.car.restoration.candidates.m80@data[stm.car.restoration.candidates.m80@data$COD_IMOVEL==i,"FOREST_COVER_INCREMENT"] <- tapply(area(forest.cover.increment), forest.cover.increment[], sum)[2]*100
    stm.car.restoration.candidates.m80@data[stm.car.restoration.candidates.m80@data$COD_IMOVEL==i,"FOREST_COVER_INCREMENT_PP"] <- ceiling((stm.car.restoration.candidates.m80@data[stm.car.restoration.candidates.m80@data$COD_IMOVEL==i,"FOREST_COVER_INCREMENT"]/stm.car.restoration.candidates.m80@data[stm.car.restoration.candidates.m80@data$COD_IMOVEL==i,"NUM_AREA"])*100)
    
  }
  
  candidate.areas.final[restored.cover][candidate.areas.final[restored.cover]==0]<-0
  j=j-1
  cat("\n>", j, "out of", nrow(stm.car.restoration.candidates.m80@data), "properties left<\n")
  
}

#length(candidate.areas.final.copy[candidate.areas.final.copy[]==1])
#length(candidate.areas.final[candidate.areas.final[]==1])


#select properties with still less than 80% forest cover

stm.car.restoration.candidates.l80 <- stm.car.restoration.candidates[stm.car.restoration.candidates$FOREST_COVER_INCREMENT_PP < 80,]

#exclding properties with geometry problems and/or at municipality border
exclude <- c("PA-1506807-FC422C55C19040B49873CD97D3ECCE2E", "PA-1506807-4356813EBAF549FBA918440953340295",
             "PA-1506807-44F27EDFADEC4E6FBBD67FC125F0CBA3", "PA-1506807-61419C70F3B84F2AA71B9531877E8EBC",
             "PA-1506807-0453989B639F49D19404CFBA83D777B0", "PA-1506807-4C4351BBFE7442F5887B98E8CC58E813",
             "PA-1506807-27D72F26BD6E4C208A654AC8115184A8", "PA-1506807-F7C55B21072E43578025F3A18198D1A6",
             "PA-1506807-18D6A81C7E27418C8935F50B7DE9C3B3", "PA-1506807-4B368AE2B276467A891EFB3397DFACD3",
             "PA-1506807-B812EB6FA13E4260947F9EBC93CC75D4", "PA-1506807-22B2C3AB789E417F935806F952F2496E",
             "PA-1506807-9B9B5EDCD5D74104B74AA861BDC21940", "PA-1506807-326E42CBC0E14E26AC9BBC709E15E0DE",
             "PA-1506807-F443F65E87104DA6B87EAF79C06523FB", "PA-1506807-D784B87346DC4B4AB4F7E9E9AC023521",
             "PA-1506807-929073F5F99A4349B27C76CF953E1454", "PA-1506807-791791E47AFE47698A653D5CDE6E58F1",
             "PA-1504752-64AD473767E74B778DB2004E40CFEA8A", "PA-1504752-F6B5CB69F58D4B33983A9A109695624F",
             "PA-1504752-DCEA75322C124C4AA0D1ACF1899F7C31", "PA-1504752-D0B7FE7337FD42C498C56BDC7D913C7F",
             "PA-1504752-512AAC3322C340F896BA9759220840C9", "PA-1504752-46305A62295848BCB6E870E6BC3E1B75",
             "PA-1501451-997D6B1AC088472F9ACB9E9A6486EA98", "PA-1506807-DE35C72AFF2343DCB0BCAAB10272DB95")

stm.car.restoration.candidates.l80 <- stm.car.restoration.candidates.l80[-which(stm.car.restoration.candidates.l80$COD_IMOVEL %in% exclude),]
#head(stm.car.restoration.candidates.l80@data)
#nrow(stm.car.restoration.candidates.l80@data)

#candidate.areas.final <- candidate.areas.final.copy

candidate.areas.final.copy <- candidate.areas.final
j=nrow(stm.car.restoration.candidates.l80@data)
#i="PA-1505502-8B1A6744B90E42C6BC856664188651AF"
for (i in stm.car.restoration.candidates.l80$COD_IMOVEL) {
  
  rural.property <- stm.car.restoration.candidates.l80[stm.car.restoration.candidates.l80$COD_IMOVEL==i,]
  
  forest.cover <- crop(forest.class, extent(rural.property))
  forest.cover <- mask(forest.cover, rural.property)
  
  restored.cover <- crop(candidate.areas.final, extent(rural.property))
  restored.cover <- mask(restored.cover, rural.property)
  
  forest.cover.increment <- sum(forest.cover, restored.cover)
  
  while (stm.car.restoration.candidates.l80@data[stm.car.restoration.candidates.l80@data$COD_IMOVEL==i,"FOREST_COVER_INCREMENT_PP"] < 80) {
    
    forest.cover.increment[forest.cover.increment[]==0] <- sample(c(0,1), size = length(forest.cover.increment[forest.cover.increment[]==0]), replace = T, prob = c(0.9,0.1))
    
    stm.car.restoration.candidates.l80@data[stm.car.restoration.candidates.l80@data$COD_IMOVEL==i,"FOREST_COVER_INCREMENT"] <- tapply(area(forest.cover.increment), forest.cover.increment[], sum)[2]*100
    stm.car.restoration.candidates.l80@data[stm.car.restoration.candidates.l80@data$COD_IMOVEL==i,"FOREST_COVER_INCREMENT_PP"] <- ceiling((stm.car.restoration.candidates.l80@data[stm.car.restoration.candidates.l80@data$COD_IMOVEL==i,"FOREST_COVER_INCREMENT"]/stm.car.restoration.candidates.l80@data[stm.car.restoration.candidates.l80@data$COD_IMOVEL==i,"NUM_AREA"])*100)
    
  }
  
  restored.cover.update <- forest.cover.increment-forest.cover
  candidate.areas.final[restored.cover.update][candidate.areas.final[restored.cover.update]==0]<-1
  j=j-1
  cat("\n>", j, "out of", nrow(stm.car.restoration.candidates.l80@data), "properties left<\n")
  
}

#length(candidate.areas.final.copy[candidate.areas.final.copy[]==1])
#length(candidate.areas.final[candidate.areas.final[]==1])

candidate.areas.final <- mask(candidate.areas.final, stm.shp)
#plot(candidate.areas.final)
#plot(stm.car.restoration.candidates, add=T)
#candidate.areas.final <- raster("rasters/STM/input/restoration_candidate_areas.tif")

rm(list=ls()[!ls() %in% c("stm.shp", "stm.lulc","stm.lulc.2010.forest.class", "stm.lulc.2020.forest.class", "dist.road",
                          "dist.river.all", "elevation", "stm.car.raster", "candidate.areas.final")]) #keeping only raster stack
gc()



#
#



# time since degradation 2010 data from RAS 
# quantitative comparison of manual inspection of satellite images
# and field observations done by two observers (TG and SN)
# see RAS environmental explanatory variable guideline document for details
# 2020 data from DETER
# 
stm.degrad.2010 <- raster("rasters/STM/input/stm-2010-deg_tsince0_150m.grd")

names(stm.degrad.2010) <- c("stm.degrad.2010real")

#checking
#stm.degrad.2010
#plot(stm.degrad.2010)
#range(values(stm.degrad.2010), na.rm=T)

# Conversion of rasters into same extent
stm.degrad.2010 <- projectRaster(stm.degrad.2010, crs = "+proj=longlat +datum=WGS84 +no_defs")
stm.degrad.2010 <- resample(stm.degrad.2010, stm.lulc, method='ngb')


# calculating time since degradation for 2020
stm.degrad.temp <- stm.degrad.2010


# deter data between 2011 and 2015
library(datazoom.amazonia)

deter.2011.15 <- load_degrad(dataset = "degrad", raw_data = T, time_period = 2011:2015)

for (year in 1:5) {   #1=2011; 5=2015
  
  deter.yearx <- deter.2011.15[[year]] 
  deter.yearx <- sf:::as_Spatial(deter.yearx$geometry)
  stm.deter.yearx <- crop(deter.yearx, extent(stm.shp))
  
  stm.deter.yearx <- rasterize(stm.deter.yearx, stm.lulc[[1]], field=999)
  stm.deter.yearx[is.na(stm.deter.yearx)]<-0
  stm.deter.yearx <- mask(stm.deter.yearx, stm.shp)
  
  stm.degrad.temp <- stm.degrad.temp+1
  stm.degrad.temp[get("stm.deter.yearx")[] == 999] <- 0
  
  cat("\n> year", year, "done! <\n")
}


#plot(stm.degrad.temp)
#plot(stm.deter.yearx)
#length(stm.deter.yearx[stm.deter.yearx[]==999])

# deter data between 2016 and 2020
deter.2016.20 <- readOGR(dsn = "rasters/STM/input", layer = "deter_public")
stm.deter.2016.20 <- crop(deter.2016.20, extent(stm.shp))

rm(deter.2016.20)

for (year in 2016:2020) {
  
  stm.deter.yearx <- stm.deter.2016.20[grep(year, stm.deter.2016.20$VIEW_DATE),]
  stm.deter.yearx <- rasterize(stm.deter.yearx, stm.lulc[[1]], field=999)
  stm.deter.yearx[is.na(stm.deter.yearx)]<-0
  stm.deter.yearx <- mask(stm.deter.yearx, stm.shp)
  
  stm.degrad.temp <- stm.degrad.temp+1
  stm.degrad.temp[get("stm.deter.yearx")[] == 999] <- 0
  
  cat("\n> year", year, "done! <\n")
}



stm.degrad.temp <- mask(stm.degrad.temp, stm.lulc.2020.forest.class)

stm.degrad.2020 <- stm.degrad.temp
#stm.degrad.2020 <- raster("rasters/STM/input/stm-2020-deg_tsince0.tif")

names(stm.degrad.2020) <- "stm.degrad.2020real"

stm.degrad <- stack(stm.degrad.2010, stm.degrad.2020)

#checking
#stm.degrad
#plot(stm.degrad)
#range(values(stm.degrad[["STM.Degrad.2010"]]), na.rm=T)

#non-degraded sites will be considered with 300 years following (BIB)
stm.degrad[["stm.degrad.2010real"]][stm.degrad[["stm.degrad.2010real"]]>23] <- 300
stm.degrad[["stm.degrad.2020real"]][stm.degrad[["stm.degrad.2020real"]]>33] <- 300

# isolating degraded forest class pixels
stm.degrad.2010.forest.class <- stm.degrad[["stm.degrad.2010real"]]
stm.degrad.2010.forest.class[stm.degrad.2010.forest.class>23]<-NA
stm.degrad.2010.forest.class[!is.na(stm.degrad.2010.forest.class)] <- 1
stm.degrad.2010.forest.class<-sum(stm.lulc.2010.forest.class, stm.degrad.2010.forest.class, na.rm=T)
stm.degrad.2010.forest.class[stm.degrad.2010.forest.class<2]<-0
stm.degrad.2010.forest.class[stm.degrad.2010.forest.class==2]<-1



stm.degrad.2020.forest.class <- stm.degrad[["stm.degrad.2020real"]]
stm.degrad.2020.forest.class[stm.degrad.2020.forest.class>33]<-NA
stm.degrad.2020.forest.class[!is.na(stm.degrad.2020.forest.class)] <- 1
stm.degrad.2020.forest.class<-sum(stm.lulc.2020.forest.class, stm.degrad.2020.forest.class, na.rm=T)
stm.degrad.2020.forest.class[stm.degrad.2020.forest.class<2]<-0
stm.degrad.2020.forest.class[stm.degrad.2020.forest.class==2]<-1

rm(list=ls()[!ls() %in% c("stm.shp", "stm.lulc","stm.lulc.2010.forest.class", "stm.lulc.2020.forest.class", "dist.road",
                          "dist.river.all", "elevation", "stm.car.raster", "candidate.areas.final", "stm.degrad", 
                          "stm.degrad.2010.forest.class", "stm.degrad.2020.forest.class")]) #keeping only raster stack
gc()



#
#



# secondary forest age from Silva Jr. et al 2020  [2010 and 2020]
# [DOI: 10.1038/s41597-020-00600-4]
# [STM] paragominas
stm.sfage <- stack(c("rasters/STM/input/stm-2010-sfage-mapbiomas-brazil-collection-60.tif",
                     "rasters/STM/input/stm-2020-sfage-mapbiomas-brazil-collection-60.tif"))
names(stm.sfage) <- c("stm.sfage.2010real", "stm.sfage.2020real")

#checking
#stm.sfage
#plot(stm.sfage)
#range(values(stm.sfage[["stm.sfage.2010real"]]), na.rm = T)

# Conversion of rasters into same extent
stm.sfage <- resample(stm.sfage, stm.lulc, method='ngb')

# isolating secondary forest class pixels
stm.sfage.2010.all.class <- stm.sfage[["stm.sfage.2010real"]]
stm.sfage.2010.all.class[stm.sfage.2010.all.class>0] <- 1
stm.sfage.2010.all.class[stm.sfage.2010.all.class<1] <- 0

stm.sfage.2020.all.class <- stm.sfage[["stm.sfage.2020real"]]
stm.sfage.2020.all.class[stm.sfage.2020.all.class>0] <- 1
stm.sfage.2020.all.class[stm.sfage.2020.all.class<1] <- 0


rm(list=ls()[!ls() %in% c("stm.shp", "stm.lulc","stm.lulc.2010.forest.class", "stm.lulc.2020.forest.class", "dist.road",
                          "dist.river.all", "elevation", "stm.car.raster", "candidate.areas.final", "stm.degrad", "stm.degrad.2010.forest.class",
                          "stm.degrad.2020.forest.class", "stm.sfage", "stm.sfage.2010.all.class", "stm.sfage.2020.all.class")]) #keeping only raster stack
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
UPF2010<-sum(stm.lulc.2010.forest.class, stm.sfage.2010.all.class, stm.degrad.2010.forest.class, na.rm = T)
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
UPF2010.px <- mask(UPF2010.px, stm.shp)

#saving
writeRaster(UPF2010.px, "rasters/STM/2010_real/UPFpx.tif", format="GTiff", overwrite=T)
writeRaster(UPF2010.px, "rasters/STM/2020_avoidboth/UPFpx.tif", format="GTiff", overwrite=T)
writeRaster(UPF2010.px, "rasters/STM/2020_restor_n_avoid_both/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (1050m)
UPF2010.ls <- focal(UPF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#UPF2010.ls
#anyNA(UPF2010.ls[])
#plot(UPF2010.ls)

names(UPF2010.ls)<-"UPFls"
UPF2010.ls[is.nan(UPF2010.ls)] <- 0
UPF2010.ls <- mask(UPF2010.ls, stm.shp)

#saving
writeRaster(UPF2010.ls, "rasters/STM/2010_real/UPFls.tif", format="GTiff", overwrite=T)
writeRaster(UPF2010.ls, "rasters/STM/2020_avoidboth/UPFls.tif", format="GTiff", overwrite=T)
writeRaster(UPF2010.ls, "rasters/STM/2020_restor_n_avoid_both/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2010.px", "UPF2010.ls")]) #keeping only raster stack
gc()



#


# scenario avoid degradation
UPF.avoiddegrad<-sum(stm.lulc.2020.forest.class, stm.sfage.2020.all.class, stm.degrad.2010.forest.class, na.rm = T)
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
UPF.avoiddegrad.px <- mask(UPF.avoiddegrad.px, stm.shp)

#saving
writeRaster(UPF.avoiddegrad.px, "rasters/STM/2020_avoiddegrad/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (1050m)
UPF.avoiddegrad.ls <- focal(UPF.avoiddegrad, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#UPF.avoiddegrad.ls
#anyNA(UPF.avoiddegrad.ls[])
#plot(UPF.avoiddegrad.ls)

names(UPF.avoiddegrad.ls)<-"UPFls"
UPF.avoiddegrad.ls[is.nan(UPF.avoiddegrad.ls)] <- 0
UPF.avoiddegrad.ls <- mask(UPF.avoiddegrad.ls, stm.shp)

#saving
writeRaster(UPF.avoiddegrad.ls, "rasters/STM/2020_avoiddegrad/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF.avoiddegrad.px", "UPF.avoiddegrad.ls")]) #keeping only raster stack
gc()



#


# scenario avoid deforestation
UPF.avoiddefor<-sum(stm.lulc.2020.forest.class, stm.sfage.2010.all.class, stm.degrad.2020.forest.class, na.rm = T)
UPF.avoiddefor[UPF.avoiddefor>1]<-0
##cheking
#unique(UPF.avoiddefor[])
#plot(UPF.avoiddefor)

# mean upf cover in pixel scale (150m)
UPF.avoiddefor.px <- focal(UPF.avoiddefor, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#UPF.avoiddefor.px
#anyNA(UPF.avoiddefor.px[])
#plot(UPF.avoiddefor.px)

names(UPF.avoiddefor.px)<-"UPFpx"
UPF.avoiddefor.px[is.nan(UPF.avoiddefor.px)] <- 0
UPF.avoiddefor.px <- mask(UPF.avoiddefor.px, stm.shp)

#saving
writeRaster(UPF.avoiddefor.px, "rasters/STM/2020_avoiddeforest/UPFpx.tif", format="GTiff", overwrite=T)
writeRaster(UPF.avoiddefor.px, "rasters/STM/2020_restor_n_avoid_deforest/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (1050m)
UPF.avoiddefor.ls <- focal(UPF.avoiddefor, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#UPF.avoiddefor.ls
#anyNA(UPF.avoiddefor.ls[])
#plot(UPF.avoiddefor.ls)

names(UPF.avoiddefor.ls)<-"UPFls"
UPF.avoiddefor.ls[is.nan(UPF.avoiddefor.ls)] <- 0
UPF.avoiddefor.ls <- mask(UPF.avoiddefor.ls, stm.shp)

#saving
writeRaster(UPF.avoiddefor.ls, "rasters/STM/2020_avoiddeforest/UPFls.tif", format="GTiff", overwrite=T)
writeRaster(UPF.avoiddefor.ls, "rasters/STM/2020_restor_n_avoid_deforest/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF.avoiddefor.px", "UPF.avoiddefor.ls")]) #keeping only raster stack
gc()



#


# scenario 2020
UPF2020<-sum(stm.lulc.2020.forest.class, stm.sfage.2020.all.class, stm.degrad.2020.forest.class, na.rm = T)
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
UPF2020.px <- mask(UPF2020.px, stm.shp)

#saving
writeRaster(UPF2020.px, "rasters/STM/2020_real/UPFpx.tif", format="GTiff", overwrite=T)
writeRaster(UPF2020.px, "rasters/STM/2020_restor_wo_avoid/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (1050m)
UPF2020.ls <- focal(UPF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#UPF2020.ls
#anyNA(UPF2020.ls[])
#plot(UPF2020.ls)

names(UPF2020.ls)<-"UPFls"
UPF2020.ls[is.nan(UPF2020.ls)] <- 0
UPF2020.ls <- mask(UPF2020.ls, stm.shp)

#saving
writeRaster(UPF2020.ls, "rasters/STM/2020_real/UPFls.tif", format="GTiff", overwrite=T)
writeRaster(UPF2020.ls, "rasters/STM/2020_restor_wo_avoid/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020.px", "UPF2020.ls")]) #keeping only raster stack
gc()



#


#######################################################################################################################

# [DPF] degraded primary forest -- pixel: 5x5 (150m); and landscape: 35x35 (1050m)
# this variable includes forest pixels in LULC raster (value == 3)
# which overlaps with pixels with fire (burned) and/or pixels degraded (burned and logged / logged)

# scenario 2010
DPF2010 <- stm.degrad.2010.forest.class
#plot(DPF2010)

# mean dpf cover in pixel scale (150m)
DPF2010.px <- focal(DPF2010, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#DPF2010.px
#anyNA(DPF2010.px[])
#plot(DPF2010.px)

names(DPF2010.px)<-"DPFpx"
DPF2010.px[is.nan(DPF2010.px)] <- 0
DPF2010.px <- mask(DPF2010.px, stm.shp)

#saving
writeRaster(DPF2010.px, "rasters/STM/2010_real/DPFpx.tif", format="GTiff", overwrite=T)
writeRaster(DPF2010.px, "rasters/STM/2020_avoiddegrad/DPFpx.tif", format="GTiff", overwrite=T)
writeRaster(DPF2010.px, "rasters/STM/2020_avoidboth/DPFpx.tif", format="GTiff", overwrite=T)
writeRaster(DPF2010.px, "rasters/STM/2020_restor_n_avoid_both/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (1050m)
DPF2010.ls <- focal(DPF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#DPF2010.ls
#anyNA(DPF2010.ls[])
#plot(DPF2010.ls)

names(DPF2010.ls)<-"DPFls"
DPF2010.ls[is.nan(DPF2010.ls)] <- 0
DPF2010.ls <- mask(DPF2010.ls, stm.shp)

#saving
writeRaster(DPF2010.ls, "rasters/STM/2010_real/DPFls.tif", format="GTiff", overwrite=T)
writeRaster(DPF2010.ls, "rasters/STM/2020_avoiddegrad/DPFls.tif", format="GTiff", overwrite=T)
writeRaster(DPF2010.ls, "rasters/STM/2020_avoidboth/DPFls.tif", format="GTiff", overwrite=T)
writeRaster(DPF2010.ls, "rasters/STM/2020_restor_n_avoid_both/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("stm.degrad.2010.forest.class", "DPF2010.px", "DPF2010.ls")])  #keeping only raster stack
gc()



#


# scenario 2020
DPF2020 <- stm.degrad.2020.forest.class
#plot(DPF2020)

# mean dpf cover in pixel scale (150m)
DPF2020.px <- focal(DPF2020, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#DPF2020.px
#anyNA(DPF2020.px[])
#plot(DPF2020.px)

names(DPF2020.px)<-"DPFpx"
DPF2020.px[is.nan(DPF2020.px)] <- 0
DPF2020.px <- mask(DPF2020.px, stm.shp)

#saving
writeRaster(DPF2020.px, "rasters/STM/2020_real/DPFpx.tif", format="GTiff", overwrite=T)
writeRaster(DPF2020.px, "rasters/STM/2020_avoiddeforest/DPFpx.tif", format="GTiff", overwrite=T)
writeRaster(DPF2020.px, "rasters/STM/2020_restor_wo_avoid/DPFpx.tif", format="GTiff", overwrite=T)
writeRaster(DPF2020.px, "rasters/STM/2020_restor_n_avoid_deforest/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (1050m)
DPF2020.ls <- focal(DPF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#DPF2020.ls
#anyNA(DPF2020.ls[])
#plot(DPF2020.ls)

names(DPF2020.ls)<-"DPFls"
DPF2020.ls[is.nan(DPF2020.ls)] <- 0
DPF2020.ls <- mask(DPF2020.ls, stm.shp)

#saving
writeRaster(DPF2020.ls, "rasters/STM/2020_real/DPFls.tif", format="GTiff", overwrite=T)
writeRaster(DPF2020.ls, "rasters/STM/2020_avoiddeforest/DPFls.tif", format="GTiff", overwrite=T)
writeRaster(DPF2020.ls, "rasters/STM/2020_restor_wo_avoid/DPFls.tif", format="GTiff", overwrite=T)
writeRaster(DPF2020.ls, "rasters/STM/2020_restor_n_avoid_deforest/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("stm.degrad.2020.forest.class", "DPF2020.px", "DPF2020.ls")]) #keeping only raster stack
gc()



#


#######################################################################################################################

# [TSD] time since degradation -- pixel: 5x5 (150m); and landscape: 35x35 (1050m)
# this variable is the mean time since a degradation event

# scenario 2010
TSD2010 <- stm.degrad[["stm.degrad.2010real"]]
#plot(TSD2010)

# mean tsd cover in pixel scale (150m)
TSD2010.px <- focal(TSD2010, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#TSD2010.px
#anyNA(TSD2010.px[])
#plot(TSD2010.px)

names(TSD2010.px)<-"TSDpx"
TSD2010.px[is.nan(TSD2010.px)] <- 0
TSD2010.px <- mask(TSD2010.px, stm.shp)

#saving
writeRaster(TSD2010.px, "rasters/STM/2010_real/TSDpx.tif", format="GTiff", overwrite=T)

# mean tsd cover in landscape scale (1050m)
TSD2010.ls <- focal(TSD2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#TSD2010.ls
#anyNA(TSD2010.ls[])
#plot(TSD2010.ls)

names(TSD2010.ls)<-"TSDls"
TSD2010.ls[is.nan(TSD2010.ls)] <- 0
TSD2010.ls <- mask(TSD2010.ls, stm.shp)

#saving
writeRaster(TSD2010.ls, "rasters/STM/2010_real/TSDls.tif", format="GTiff", overwrite=T)

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
TSD2010.recovery10.px <- mask(TSD2010.recovery10.px, stm.shp)

#saving
writeRaster(TSD2010.recovery10.px, "rasters/STM/2020_avoiddegrad/TSDpx.tif", format="GTiff", overwrite=T)
writeRaster(TSD2010.recovery10.px, "rasters/STM/2020_avoidboth/TSDpx.tif", format="GTiff", overwrite=T)
writeRaster(TSD2010.recovery10.px, "rasters/STM/2020_restor_n_avoid_both/TSDpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (1050m)
TSD2010.recovery10.ls <- focal(TSD2010.recovery10, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#TSD2010.recovery10.ls
#anyNA(TSD2010.recovery10.ls[])
#plot(TSD2010.recovery10.ls)

names(TSD2010.recovery10.ls)<-"TSDls"
TSD2010.recovery10.ls[is.nan(TSD2010.recovery10.ls)] <- 0
TSD2010.recovery10.ls <- mask(TSD2010.recovery10.ls, stm.shp)

#saving
writeRaster(TSD2010.recovery10.ls, "rasters/STM/2020_avoiddegrad/TSDls.tif", format="GTiff", overwrite=T)
writeRaster(TSD2010.recovery10.ls, "rasters/STM/2020_avoidboth/TSDls.tif", format="GTiff", overwrite=T)
writeRaster(TSD2010.recovery10.ls, "rasters/STM/2020_restor_n_avoid_both/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2010.recovery10.px", "TSD2010.recovery10.ls")]) #keeping only raster stack
gc()



#


# scenario 2020
TSD2020 <- stm.degrad[["stm.degrad.2020real"]]
#plot(TSD2020)

# mean tsd cover in pixel scale (150m)
TSD2020.px <- focal(TSD2020, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#TSD2020.px
#anyNA(TSD2020.px[])
#plot(TSD2020.px)

names(TSD2020.px)<-"TSDpx"
TSD2020.px[is.nan(TSD2020.px)] <- 0
TSD2020.px <- mask(TSD2020.px, stm.shp)

#saving
writeRaster(TSD2020.px, "rasters/STM/2020_real/TSDpx.tif", format="GTiff", overwrite=T)
writeRaster(TSD2020.px, "rasters/STM/2020_avoiddeforest/TSDpx.tif", format="GTiff", overwrite=T)
writeRaster(TSD2020.px, "rasters/STM/2020_restor_wo_avoid/TSDpx.tif", format="GTiff", overwrite=T)
writeRaster(TSD2020.px, "rasters/STM/2020_restor_n_avoid_deforest/TSDpx.tif", format="GTiff", overwrite=T)

# mean tsd cover in landscape scale (1050m)
TSD2020.ls <- focal(TSD2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#TSD2020.ls
#anyNA(TSD2020.ls[])
#plot(TSD2020.ls)

names(TSD2020.ls)<-"TSDls"
TSD2020.ls[is.nan(TSD2020.ls)] <- 0
TSD2020.ls <- mask(TSD2020.ls, stm.shp)

#saving
writeRaster(TSD2020.ls, "rasters/STM/2020_real/TSDls.tif", format="GTiff", overwrite=T)
writeRaster(TSD2020.ls, "rasters/STM/2020_avoiddeforest/TSDls.tif", format="GTiff", overwrite=T)
writeRaster(TSD2020.ls, "rasters/STM/2020_restor_wo_avoid/TSDls.tif", format="GTiff", overwrite=T)
writeRaster(TSD2020.ls, "rasters/STM/2020_restor_n_avoid_deforest/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020.px", "TSD2020.ls")]) #keeping only raster stack
gc()



#


#######################################################################################################################

# [SF] secondary forest -- pixel: 5x5 (150m); and landscape: 35x35 (1050m)
# this variable includes forest pixels in SFage raster
# which has less than 25 years for 2010 or less than 35 for 2020

# scenario 2010
SF2010 <- stm.sfage.2010.all.class
#plot(SF2010)

# mean sf cover in pixel scale (150m)
SF2010.px <- focal(SF2010, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#SF2010.px
#anyNA(SF2010.px[])
#plot(SF2010.px)

names(SF2010.px)<-"SFpx"
SF2010.px[is.nan(SF2010.px)] <- 0
SF2010.px <- mask(SF2010.px, stm.shp)

#saving
writeRaster(SF2010.px, "rasters/STM/2010_real/SFpx.tif", format="GTiff", overwrite=T)
writeRaster(SF2010.px, "rasters/STM/2020_avoiddeforest/SFpx.tif", format="GTiff", overwrite=T)
writeRaster(SF2010.px, "rasters/STM/2020_avoidboth/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (1050m)
SF2010.ls <- focal(SF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SF2010.ls
#anyNA(SF2010.ls[])
#plot(SF2010.ls)

names(SF2010.ls)<-"SFls"
SF2010.ls[is.nan(SF2010.ls)] <- 0
SF2010.ls <- mask(SF2010.ls, stm.shp)

#saving
writeRaster(SF2010.ls, "rasters/STM/2010_real/SFls.tif", format="GTiff", overwrite=T)
writeRaster(SF2010.ls, "rasters/STM/2020_avoiddeforest/SFls.tif", format="GTiff", overwrite=T)
writeRaster(SF2010.ls, "rasters/STM/2020_avoidboth/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2010.px", "SF2010.ls")]) #keeping only raster stack
gc()



#


# scenario 2020
SF2020 <- stm.sfage.2020.all.class
#plot(SF2020)

# mean sf cover in pixel scale (150m)
SF2020.px <- focal(SF2020, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#SF2020.px
#anyNA(SF2020.px[])
#plot(SF2020.px)

names(SF2020.px)<-"SFpx"
SF2020.px[is.nan(SF2020.px)] <- 0
SF2020.px <- mask(SF2020.px, stm.shp)

#saving
writeRaster(SF2020.px, "rasters/STM/2020_real/SFpx.tif", format="GTiff", overwrite=T)
writeRaster(SF2020.px, "rasters/STM/2020_avoiddegrad/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (1050m)
SF2020.ls <- focal(SF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SF2020.ls
#anyNA(SF2020.ls[])
#plot(SF2020.ls)

names(SF2020.ls)<-"SFls"
SF2020.ls[is.nan(SF2020.ls)] <- 0
SF2020.ls <- mask(SF2020.ls, stm.shp)

#saving
writeRaster(SF2020.ls, "rasters/STM/2020_real/SFls.tif", format="GTiff", overwrite=T)
writeRaster(SF2020.ls, "rasters/STM/2020_avoiddegrad/SFls.tif", format="GTiff", overwrite=T)

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
SF2010.restore10.px <- mask(SF2010.restore10.px, stm.shp)

#saving
writeRaster(SF2010.restore10.px, "rasters/STM/2020_restor_n_avoid_deforest/SFpx.tif", format="GTiff", overwrite=T)
writeRaster(SF2010.restore10.px, "rasters/STM/2020_restor_n_avoid_both/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (1050m)
SF2010.restore10.ls <- focal(SF2010.restore10, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SF2010.restore10.ls
#anyNA(SF2010.restore10.ls[])
#plot(SF2010.restore10.ls)

names(SF2010.restore10.ls)<-"SFls"
SF2010.restore10.ls[is.nan(SF2010.restore10.ls)] <- 0
SF2010.restore10.ls <- mask(SF2010.restore10.ls, stm.shp)

#saving
writeRaster(SF2010.restore10.ls, "rasters/STM/2020_restor_n_avoid_deforest/SFls.tif", format="GTiff", overwrite=T)
writeRaster(SF2010.restore10.ls, "rasters/STM/2020_restor_n_avoid_both/SFls.tif", format="GTiff", overwrite=T)

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
SF2020.restore10.px <- mask(SF2020.restore10.px, stm.shp)

#saving
writeRaster(SF2020.restore10.px, "rasters/STM/2020_restor_wo_avoid/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (1050m)
SF2020.restore10.ls <- focal(SF2020.restore10, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SF2020.restore10.ls
#anyNA(SF2020.restore10.ls[])
#plot(SF2020.restore10.ls)

names(SF2020.restore10.ls)<-"SFls"
SF2020.restore10.ls[is.nan(SF2020.restore10.ls)] <- 0
SF2020.restore10.ls <- mask(SF2020.restore10.ls, stm.shp)

#saving
writeRaster(SF2020.restore10.ls, "rasters/STM/2020_restor_wo_avoid/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020.restore10.px", "SF2020.restore10.ls")]) #keeping only raster stack
gc()



#


#######################################################################################################################

# [SFage] secondary forest age -- pixel: 5x5 (150m); and landscape: 35x35 (1050m)
# this variable is the mean age of secondary forest

# scenario 2010
SFage2010 <- stm.sfage[["stm.sfage.2010real"]]
#plot(SFage2010)

# mean sf age in pixel scale (150m)
SFage2010.px <- focal(SFage2010, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#SFage2010.px
#anyNA(SFage2010.px[])
#plot(SFage2010.px)

names(SFage2010.px)<-"SFagepx"
SFage2010.px[is.nan(SFage2010.px)] <- 0
SFage2010.px <- mask(SFage2010.px, stm.shp)

#saving
writeRaster(SFage2010.px, "rasters/STM/2010_real/SFagepx.tif", format="GTiff", overwrite=T)

# mean sf age in landscape scale (1050m)
SFage2010.ls <- focal(SFage2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SFage2010.ls
#anyNA(SFage2010.ls[])
#plot(SFage2010.ls)

names(SFage2010.ls)<-"SFagels"
SFage2010.ls[is.nan(SFage2010.ls)] <- 0
SFage2010.ls <- mask(SFage2010.ls, stm.shp)

#saving
writeRaster(SFage2010.ls, "rasters/STM/2010_real/SFagels.tif", format="GTiff", overwrite=T)

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
SFage2010.recovery10.px <- mask(SFage2010.recovery10.px, stm.shp)

#saving
writeRaster(SFage2010.recovery10.px, "rasters/STM/2020_avoiddeforest/SFagepx.tif", format="GTiff", overwrite=T)
writeRaster(SFage2010.recovery10.px, "rasters/STM/2020_avoidboth/SFagepx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (1050m)
SFage2010.recovery10.ls <- focal(SFage2010.recovery10, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SFage2010.recovery10.ls
#anyNA(SFage2010.recovery10.ls[])
#plot(SFage2010.recovery10.ls)

names(SFage2010.recovery10.ls)<-"SFagels"
SFage2010.recovery10.ls[is.nan(SFage2010.recovery10.ls)] <- 0
SFage2010.recovery10.ls <- mask(SFage2010.recovery10.ls, stm.shp)

#saving
writeRaster(SFage2010.recovery10.ls, "rasters/STM/2020_avoiddeforest/SFagels.tif", format="GTiff", overwrite=T)
writeRaster(SFage2010.recovery10.ls, "rasters/STM/2020_avoidboth/SFagels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFage2010.recovery10.px", "SFage2010.recovery10.ls")]) #keeping only raster stack
gc()



#


# scenario 2020
SFage2020 <- stm.sfage[["stm.sfage.2020real"]]
#plot(SFage2010)

# mean sf age in pixel scale (150m)
SFage2020.px <- focal(SFage2020, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#SFage2020.px
#anyNA(SFage2020.px[])
#plot(SFage2020.px)

names(SFage2020.px)<-"SFagepx"
SFage2020.px[is.nan(SFage2020.px)] <- 0
SFage2020.px <- mask(SFage2020.px, stm.shp)

#saving
writeRaster(SFage2020.px, "rasters/STM/2020_real/SFagepx.tif", format="GTiff", overwrite=T)
writeRaster(SFage2020.px, "rasters/STM/2020_avoiddegrad/SFagepx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (1050m)
SFage2020.ls <- focal(SFage2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SFage2020.ls
#anyNA(SFage2020.ls[])
#plot(SFage2020.ls)

names(SFage2020.ls)<-"SFagels"
SFage2020.ls[is.nan(SFage2020.ls)] <- 0
SFage2020.ls <- mask(SFage2020.ls, stm.shp)

#saving
writeRaster(SFage2020.ls, "rasters/STM/2020_real/SFagels.tif", format="GTiff", overwrite=T)
writeRaster(SFage2020.ls, "rasters/STM/2020_avoiddegrad/SFagels.tif", format="GTiff", overwrite=T)

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
SFAge2010.restore10.px <- mask(SFAge2010.restore10.px, stm.shp)

#saving
writeRaster(SFAge2010.restore10.px, "rasters/STM/2020_restor_n_avoid_deforest/SFagepx.tif", format="GTiff", overwrite=T)
writeRaster(SFAge2010.restore10.px, "rasters/STM/2020_restor_n_avoid_both/SFagepx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (1050m)
SFAge2010.restore10.ls <- focal(SFAge2010.restore10, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SFAge2010.restore10.ls
#anyNA(SFAge2010.restore10.ls[])
#plot(SFAge2010.restore10.ls)

names(SFAge2010.restore10.ls)<-"SFagels"
SFAge2010.restore10.ls[is.nan(SFAge2010.restore10.ls)] <- 0
SFAge2010.restore10.ls <- mask(SFAge2010.restore10.ls, stm.shp)

#saving
writeRaster(SFAge2010.restore10.ls, "rasters/STM/2020_restor_n_avoid_deforest/SFagels.tif", format="GTiff", overwrite=T)
writeRaster(SFAge2010.restore10.ls, "rasters/STM/2020_restor_n_avoid_both/SFagels.tif", format="GTiff", overwrite=T)

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
SFAge2020.restore10.px <- mask(SFAge2020.restore10.px, stm.shp)

#saving
writeRaster(SFAge2020.restore10.px, "rasters/STM/2020_restor_wo_avoid/SFagepx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (1050m)
SFAge2020.restore10.ls <- focal(SFAge2020.restore10, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#SFAge2020.restore10.ls
#anyNA(SFAge2020.restore10.ls[])
#plot(SFAge2020.restore10.ls)

names(SFAge2020.restore10.ls)<-"SFagels"
SFAge2020.restore10.ls[is.nan(SFAge2020.restore10.ls)] <- 0
SFAge2020.restore10.ls <- mask(SFAge2020.restore10.ls, stm.shp)

#saving
writeRaster(SFAge2020.restore10.ls, "rasters/STM/2020_restor_wo_avoid/SFagels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020.restore10.px", "SFAge2020.restore10.ls")]) #keeping only raster stack
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
TF2010.px <- mask(TF2010.px, stm.shp)

#saving
writeRaster(TF2010.px, "rasters/STM/2010_real/TFpx.tif", format="GTiff", overwrite=T)
writeRaster(TF2010.px, "rasters/STM/2020_avoidboth/TFpx.tif", format="GTiff", overwrite=T)

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
TF2020.px <- mask(TF2020.px, stm.shp)

#saving
writeRaster(TF2020.px, "rasters/STM/2020_real/TFpx.tif", format="GTiff", overwrite=T)

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
TF.avoiddegrad.px <- mask(TF.avoiddegrad.px, stm.shp)

#saving
writeRaster(TF.avoiddegrad.px, "rasters/STM/2020_avoiddegrad/TFpx.tif", format="GTiff", overwrite=T)

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
TF.avoiddefor.px <- mask(TF.avoiddefor.px, stm.shp)

#saving
writeRaster(TF.avoiddefor.px, "rasters/STM/2020_avoiddeforest/TFpx.tif", format="GTiff", overwrite=T)

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
TF.restore10.a.px <- mask(TF.restore10.a.px, stm.shp)

#saving
writeRaster(TF.restore10.a.px, "rasters/STM/2020_restor_wo_avoid/TFpx.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TF.restore10.a.px")]) #keeping only raster stack
gc()



#


# scenario restoration and avoid deforestation
SFAge2010.restore10.young <- SFAge2010.restore10
SFAge2010.restore10.young[SFAge2010.restore10.young <= 2] <- 0
SFAge2010.restore10.young[SFAge2010.restore10.young > 2] <- 1

TF.restore10.b <- sum(UPF.avoiddefor, DPF2020, SFAge2010.restore10.young, na.rm = T)
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
TF.restore10.b.px <- mask(TF.restore10.b.px, stm.shp)

#saving
writeRaster(TF.restore10.b.px, "rasters/STM/2020_restor_n_avoid_deforest/TFpx.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TF.restore10.b.px")]) #keeping only raster stack
gc()



#


# scenario restoration and avoid both
TF.restore10.c <- sum(UPF2010, DPF2010, SFAge2010.restore10.young, na.rm = T)
TF.restore10.c[TF.restore10.c>1] <- 1
##cheking
#TF.restore10.c
#plot(TF.restore10.c)

# mean upf cover in pixel scale (150m)
TF.restore10.c.px <- focal(TF.restore10.c, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#TF.restore10.c.px
#anyNA(TF.restore10.c.px[])
#plot(TF.restore10.c.px)

names(TF.restore10.c.px)<-"TFpx"
TF.restore10.c.px[is.nan(TF.restore10.c.px)] <- 0
TF.restore10.c.px <- mask(TF.restore10.c.px, stm.shp)

#saving
writeRaster(TF.restore10.c.px, "rasters/STM/2020_restor_n_avoid_both/TFpx.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TF.restore10.c.px")]) #keeping only raster stack
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
MF2010.px <- mask(MF2010.px, stm.shp)

#saving
writeRaster(MF2010.px, "rasters/STM/2010_real/MFpx.tif", format="GTiff", overwrite=T)
writeRaster(MF2010.px, "rasters/STM/2020_avoidboth/MFpx.tif", format="GTiff", overwrite=T)

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
MF2020.px <- mask(MF2020.px, stm.shp)

#saving
writeRaster(MF2020.px, "rasters/STM/2020_real/MFpx.tif", format="GTiff", overwrite=T)

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
MF.avoiddegrad.px <- mask(MF.avoiddegrad.px, stm.shp)

#saving
writeRaster(MF.avoiddegrad.px, "rasters/STM/2020_avoiddegrad/MFpx.tif", format="GTiff", overwrite=T)

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
MF.avoiddefor.px <- mask(MF.avoiddefor.px, stm.shp)

#saving
writeRaster(MF.avoiddefor.px, "rasters/STM/2020_avoiddeforest/MFpx.tif", format="GTiff", overwrite=T)

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
MF.restore10.a.px <- mask(MF.restore10.a.px, stm.shp)

#saving
writeRaster(MF.restore10.a.px, "rasters/STM/2020_restor_wo_avoid/MFpx.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("MF.restore10.a.px")]) #keeping only raster stack
gc()



#


# scenario restoration and avoid deforestation
SFAge2010.restore10.mature <- SFAge2010.restore10
SFAge2010.restore10.mature[SFAge2010.restore10.mature <= 10] <- 0
SFAge2010.restore10.mature[SFAge2010.restore10.mature > 10] <- 1

MF.restore10.b <- sum(UPF.avoiddefor, DPF2020, SFAge2010.restore10.mature, na.rm = T)
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
MF.restore10.b.px <- mask(MF.restore10.b.px, stm.shp)

#saving
writeRaster(MF.restore10.b.px, "rasters/STM/2020_restor_n_avoid_deforest/MFpx.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("MF.restore10.b.px")]) #keeping only raster stack
gc()



#


# scenario restoration and avoid both
MF.restore10.c <- sum(UPF2010, DPF2010, SFAge2010.restore10.mature, na.rm = T)
MF.restore10.c[MF.restore10.c>1] <- 1
##cheking
#MF.restore10.c
#plot(MF.restore10.c)

# mean upf cover in pixel scale (150m)
MF.restore10.c.px <- focal(MF.restore10.c, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#MF.restore10.c.px
#anyNA(MF.restore10.c.px[])
#plot(MF.restore10.c.px)

names(MF.restore10.c.px)<-"MFpx"
MF.restore10.c.px[is.nan(MF.restore10.c.px)] <- 0
MF.restore10.c.px <- mask(MF.restore10.c.px, stm.shp)

#saving
writeRaster(MF.restore10.c.px, "rasters/STM/2020_restor_n_avoid_both/MFpx.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("MF.restore10.c.px")]) #keeping only raster stack
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
edge.dist.2010 <- mask(edge.dist.2010, stm.shp)

#saving
writeRaster(edge.dist.2010, "rasters/STM/2010_real/edgedist.tif", format="GTiff", overwrite=T)
writeRaster(edge.dist.2010, "rasters/STM/2020_avoidboth/edgedist.tif", format="GTiff", overwrite=T)

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
edge.dist.2020 <- mask(edge.dist.2020, stm.shp)

#saving
writeRaster(edge.dist.2020, "rasters/STM/2020_real/edgedist.tif", format="GTiff", overwrite=T)

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
edge.dist.avoiddegrad <- mask(edge.dist.avoiddegrad, stm.shp)

#saving
writeRaster(edge.dist.avoiddegrad, "rasters/STM/2020_avoiddegrad/edgedist.tif", format="GTiff", overwrite=T)

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
edge.dist.avoiddefor <- mask(edge.dist.avoiddefor, stm.shp)

#saving
writeRaster(edge.dist.avoiddefor, "rasters/STM/2020_avoiddeforest/edgedist.tif", format="GTiff", overwrite=T)

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
edge.dist.restore10.a <- mask(edge.dist.restore10.a, stm.shp)

#saving
writeRaster(edge.dist.restore10.a, "rasters/STM/2020_restor_wo_avoid/edgedist.tif", format="GTiff", overwrite=T)

gc()



#


# scenario restoration and avoid deforestation
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
edge.dist.restore10.b <- mask(edge.dist.restore10.b, stm.shp)

#saving
writeRaster(edge.dist.restore10.b, "rasters/STM/2020_restor_n_avoid_deforest/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("pts", "core")]) #keeping only raster stack
gc()



#


# scenario restoration and avoid both
pts <- rasterToPoints(MF.restore10.b, spatial=TRUE)
core <- pts[pts$layer == "0",]
#cheking
#core
#plot(core, add=T, col="red")

edge.dist.restore10.c <- rasterDistance(pts, core, reference = MF.restore10.b, scale=TRUE)
##cheking
#edge.dist.restore10.c
#anyNA(edge.dist.restore10.c[])
#plot(edge.dist.restore10.c)

names(edge.dist.restore10.c)<-"edgedist"
edge.dist.restore10.c[is.nan(edge.dist.restore10.c)] <- 0
edge.dist.restore10.c <- mask(edge.dist.restore10.c, stm.shp)

#saving
writeRaster(edge.dist.restore10.c, "rasters/STM/2020_restor_n_avoid_both/edgedist.tif", format="GTiff", overwrite=T)

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
edge2010.px <- mask(edge2010.px, stm.shp)

#saving
writeRaster(edge2010.px, "rasters/STM/2010_real/edgepx.tif", format="GTiff", overwrite=T)
writeRaster(edge2010.px, "rasters/STM/2020_avoidboth/edgepx.tif", format="GTiff", overwrite=T)
#

# mean sf cover in landscape scale (1050m)
edge2010.ls <- focal(edge2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#edge2010.ls
#anyNA(edge2010.ls[])
#plot(edge2010.ls)

names(edge2010.ls)<-"edgels"
edge2010.ls[is.nan(edge2010.ls)] <- 0
edge2010.ls <- mask(edge2010.ls, stm.shp)

#saving
writeRaster(edge2010.px, "rasters/STM/2010_real/edgels.tif", format="GTiff", overwrite=T)
writeRaster(edge2010.px, "rasters/STM/2020_avoidboth/edgels.tif", format="GTiff", overwrite=T)

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
edge2020.px <- mask(edge2020.px, stm.shp)

#saving
writeRaster(edge2020.px, "rasters/STM/2020_real/edgepx.tif", format="GTiff", overwrite=T)
#

# mean sf cover in landscape scale (1050m)
edge2020.ls <- focal(edge2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#edge2020.ls
#anyNA(edge2020.ls[])
#plot(edge2020.ls)

names(edge2020.ls)<-"edgels"
edge2020.ls[is.nan(edge2020.ls)] <- 0
edge2020.ls <- mask(edge2020.ls, stm.shp)

#saving
writeRaster(edge2020.ls, "rasters/STM/2020_real/edgels.tif", format="GTiff", overwrite=T)

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
edge.avoiddegrad.px <- mask(edge.avoiddegrad.px, stm.shp)

#saving
writeRaster(edge.avoiddegrad.px, "rasters/STM/2020_avoiddegrad/edgepx.tif", format="GTiff", overwrite=T)
#

# mean sf cover in landscape scale (1050m)
edge.avoiddegrad.ls <- focal(edge.avoiddegrad, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#edge.avoiddegrad.ls
#anyNA(edge.avoiddegrad.ls[])
#plot(edge.avoiddegrad.ls)

names(edge.avoiddegrad.ls)<-"edgels"
edge.avoiddegrad.ls[is.nan(edge.avoiddegrad.ls)] <- 0
edge.avoiddegrad.ls <- mask(edge.avoiddegrad.ls, stm.shp)

#saving
writeRaster(edge.avoiddegrad.ls, "rasters/STM/2020_avoiddegrad/edgels.tif", format="GTiff", overwrite=T)

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
edge.avoiddefor.px <- mask(edge.avoiddefor.px, stm.shp)

#saving
writeRaster(edge.avoiddefor.px, "rasters/STM/2020_avoiddeforest/edgepx.tif", format="GTiff", overwrite=T)
#

# mean sf cover in landscape scale (1050m)
edge.avoiddefor.ls <- focal(edge.avoiddefor, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#edge.avoiddefor.ls
#anyNA(edge.avoiddefor.ls[])
#plot(edge.avoiddefor.ls)

names(edge.avoiddefor.ls)<-"edgels"
edge.avoiddefor.ls[is.nan(edge.avoiddefor.ls)] <- 0
edge.avoiddefor.ls <- mask(edge.avoiddefor.ls, stm.shp)

#saving
writeRaster(edge.avoiddefor.ls, "rasters/STM/2020_avoiddeforest/edgels.tif", format="GTiff", overwrite=T)

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
edge.restore10.a.px <- mask(edge.restore10.a.px, stm.shp)

#saving
writeRaster(edge.restore10.a.px, "rasters/STM/2020_restor_wo_avoid/edgepx.tif", format="GTiff", overwrite=T)
#

# mean sf cover in landscape scale (1050m)
edge.restore10.a.ls <- focal(edge.restore10.a, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#edge.restore10.a.ls
#anyNA(edge.restore10.a.ls[])
#plot(edge.restore10.a.ls)

names(edge.restore10.a.ls)<-"edgels"
edge.restore10.a.ls[is.nan(edge.restore10.a.ls)] <- 0
edge.restore10.a.ls <- mask(edge.restore10.a.ls, stm.shp)

#saving
writeRaster(edge.restore10.a.ls, "rasters/STM/2020_restor_wo_avoid/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge.restore10.a.px", "edge.restore10.a.ls")]) #keeping only raster stack
gc()



#


# scenario restoration and avoid deforestation
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
edge.restore10.b.px <- mask(edge.restore10.b.px, stm.shp)

#saving
writeRaster(edge.restore10.b.px, "rasters/STM/2020_restor_n_avoid_deforest/edgepx.tif", format="GTiff", overwrite=T)
#

# mean sf cover in landscape scale (1050m)
edge.restore10.b.ls <- focal(edge.restore10.b, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#edge.restore10.b.ls
#anyNA(edge.restore10.b.ls[])
#plot(edge.restore10.b.ls)

names(edge.restore10.b.ls)<-"edgels"
edge.restore10.b.ls[is.nan(edge.restore10.b.ls)] <- 0
edge.restore10.b.ls <- mask(edge.restore10.b.ls, stm.shp)

#saving
writeRaster(edge.restore10.b.ls, "rasters/STM/2020_restor_n_avoid_deforest/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge.restore10.b.px", "edge.restore10.b.ls")]) #keeping only raster stack
gc()



#


# scenario restoration and avoid both
# marking edge and core areas
edge.restore10.c <- focal(MF.restore10.c, matrix(1,ncol=3,nrow=3), fun = mean, na.rm = T, pad = T)
edge.restore10.c <- edge.restore10.c + MF.restore10.c
edge.restore10.c[edge.restore10.c == 2] <- 0                  # core area
edge.restore10.c[edge.restore10.c > 0 & edge.restore10.c < 2] <- 1    # edge
#cheking
#edge.restore10.c
#unique(edge.restore10.c[])
#plot(edge.restore10.c)

# mean edge in pixel scale (150m)
edge.restore10.c.px <- focal(edge.restore10.c, matrix(1,ncol=5,nrow=5), fun=mean, na.rm=T)
##cheking
#edge.restore10.c.px
#anyNA(edge.restore10.c.px[])
#plot(edge.restore10.c.px)

names(edge.restore10.c.px)<-"edgepx"
edge.restore10.c.px[is.nan(edge.restore10.c.px)] <- 0
edge.restore10.c.px <- mask(edge.restore10.c.px, stm.shp)

#saving
writeRaster(edge.restore10.c.px, "rasters/STM/2020_restor_n_avoid_both/edgepx.tif", format="GTiff", overwrite=T)
#

# mean sf cover in landscape scale (1050m)
edge.restore10.c.ls <- focal(edge.restore10.c, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
##cheking
#edge.restore10.c.ls
#anyNA(edge.restore10.c.ls[])
#plot(edge.restore10.c.ls)

names(edge.restore10.c.ls)<-"edgels"
edge.restore10.c.ls[is.nan(edge.restore10.c.ls)] <- 0
edge.restore10.c.ls <- mask(edge.restore10.c.ls, stm.shp)

#saving
writeRaster(edge.restore10.c.ls, "rasters/STM/2020_restor_n_avoid_both/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge.restore10.c.px", "edge.restore10.c.ls")]) #keeping only raster stack
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
#dir.create("rasters/STM/input/climate")


# scenario 2010
temp.list <- list.files("rasters/STM/input/climate", "LSTD", full.names = T, recursive = T)

temp2010.list <- grep("2010", temp.list, value = T)
temp2010 <- stack(temp2010.list)
#plot(temp2010)

meantemp2010 <- mean(temp2010, na.rm=T)
stm.meantemp2010 <- crop(meantemp2010, extent(stm.lulc.2010.forest.class))
#plot(stm.meantemp2010)
#plot(stm.shp, add=T)

stm.meantemp2010 <- resample(stm.meantemp2010, stm.lulc.2010.forest.class, method='bilinear')
stm.meantemp2010 <- mask(stm.meantemp2010, stm.shp)
#plot(stm.meantemp2010)

#saving
writeRaster(stm.meantemp2010, "rasters/STM/2010_real/meantemps.tif", format="GTiff", overwrite=T)
#

# scenario 2020
temp2020.list <- grep("2020", temp.list, value = T)
temp2020 <- stack(temp2020.list)
#plot(temp2020)
rm(temp2020.list)

meantemp2020 <- mean(temp2020, na.rm=T)
stm.meantemp2020 <- crop(meantemp2020, extent(stm.lulc.2020.forest.class))
#plot(stm.meantemp2020)
#plot(stm.shp, add=T)

stm.meantemp2020 <- resample(stm.meantemp2020, stm.lulc.2020.forest.class, method='bilinear')
stm.meantemp2020 <- mask(stm.meantemp2020, stm.shp)
#plot(stm.meantemp2020)

#saving
writeRaster(stm.meantemp2020, "rasters/STM/2020_real/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(stm.meantemp2020, "rasters/STM/2020_avoiddeforest/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(stm.meantemp2020, "rasters/STM/2020_avoiddegrad/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(stm.meantemp2020, "rasters/STM/2020_avoidboth/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(stm.meantemp2020, "rasters/STM/2020_restor_wo_avoid/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(stm.meantemp2020, "rasters/STM/2020_restor_n_avoid_deforest/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(stm.meantemp2020, "rasters/STM/2020_restor_n_avoid_both/meantemps.tif", format="GTiff", overwrite=T)
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
precip.list <- list.files("rasters/STM/input/climate", "GPM", full.names = T, recursive = T)

precip2010.list <- grep("2010", precip.list, value = T)
precip2010 <- stack(precip2010.list)
#plot(precip2010)

meanprecip2010 <- mean(precip2010, na.rm=T)
stm.meanprecip2010 <- crop(meanprecip2010, extent(stm.lulc.2010.forest.class))
#plot(stm.meanprecip2010)
#plot(stm.shp, add=T)

stm.meanprecip2010 <- resample(stm.meanprecip2010, stm.lulc.2010.forest.class, method='bilinear')
stm.meanprecip2010 <- mask(stm.meanprecip2010, stm.shp)
#plot(stm.meanprecip2010)

#saving
writeRaster(stm.meanprecip2010, "rasters/STM/2010_real/meanprecips.tif", format="GTiff", overwrite=T)
#

# scenario 2020
precip2020.list <- grep("2020", precip.list, value = T)
precip2020 <- stack(precip2020.list)
#plot(precip2020)
rm(precip2020.list)

meanprecip2020 <- mean(precip2020, na.rm=T)
stm.meanprecip2020 <- crop(meanprecip2020, extent(stm.lulc.2020.forest.class))
#plot(stm.meanprecip2020)
#plot(stm.shp, add=T)

stm.meanprecip2020 <- resample(stm.meanprecip2020, stm.lulc.2020.forest.class, method='bilinear')
stm.meanprecip2020 <- mask(stm.meanprecip2020, stm.shp)
#plot(stm.meanprecip2020)

#saving
writeRaster(stm.meanprecip2020, "rasters/STM/2020_real/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(stm.meanprecip2020, "rasters/STM/2020_avoiddeforest/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(stm.meanprecip2020, "rasters/STM/2020_avoiddegrad/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(stm.meanprecip2020, "rasters/STM/2020_avoidboth/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(stm.meanprecip2020, "rasters/STM/2020_restor_wo_avoid/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(stm.meanprecip2020, "rasters/STM/2020_restor_n_avoid_deforest/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(stm.meanprecip2020, "rasters/STM/2020_restor_n_avoid_both/meanprecips.tif", format="GTiff", overwrite=T)
#

#optional
#save.image("~/conserv_opportunities_jamesthomson/github_repo/stm_environment.RData")
rm(list=ls())
gc()



#


##############################################################################################################################################################################################################################################

# [distmarket] distance to municipality nucleus

stm.crs<-"+proj=longlat +datum=WGS84 +no_defs"

stm.munic.nucleus <- data.frame(ID = "stm", long = -54.7009, lat = -2.45063)

stm.munic.nucleus.coord <- SpatialPointsDataFrame(coords = stm.munic.nucleus[,c("long","lat")], 
                                                  data = stm.munic.nucleus, 
                                                  proj4string = CRS(stm.crs))

stm.munic.nucleus.coord <- spTransform(x = stm.munic.nucleus.coord, CRSobj = stm.crs)

distmarket <- distanceFromPoints(object = stm.lulc, xy = stm.munic.nucleus.coord)

distmarket <- mask(distmarket, stm.shp)

#saving
writeRaster(distmarket, "rasters/STM/2010_real/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/STM/2020_real/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/STM/2020_avoiddeforest/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/STM/2020_avoiddegrad/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/STM/2020_avoidboth/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/STM/2020_restor_wo_avoid/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/STM/2020_restor_n_avoid_deforest/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/STM/2020_restor_n_avoid_both/distmarket.tif", format="GTiff", overwrite=T)



#


#######################################################################################################################

# Creating forest mask for conservation action costs
# each land cover category is assigned a value

# scenario 2010
UPF2010.mask <- UPF2010

DPF2010.mask <- DPF2010

SF2010.mask <- SF2010

TF2010.mask <- sum(UPF2010.mask, DPF2010.mask, SF2010.mask, na.rm = T)
##cheking
#sort(unique(values(TF2010.mask)))
##handle overlaps
TF2010.mask[TF2010.mask>1] <- 1
#plot(TF2010.mask)
writeRaster(TF2010.mask, "rasters/STM/all_forest_mask/STM_2010_real.tif", format = "GTiff", overwrite = T)


#


# scenario 2020
UPF2020.mask <- UPF2020

DPF2020.mask <- DPF2020

SF2020.mask <- SF2020

TF2020.mask <- sum(UPF2020.mask, DPF2020.mask, SF2020.mask, na.rm = T)
##cheking
#sort(unique(values(TF2020.mask)))
##handle overlaps
TF2020.mask[TF2020.mask>1] <- 1 #sf over upf
#plot(TF2020.mask)
writeRaster(TF2020.mask, "rasters/STM/all_forest_mask/STM_2020_real.tif", format = "GTiff", overwrite = T)


#


# scenario avoid degradation
UPF.avoiddegrad.mask <- UPF.avoiddegrad

TF.avoiddegrad.mask <- sum(UPF.avoiddegrad.mask, DPF2010.mask, SF2020.mask, na.rm = T)
##cheking
#sort(unique(values(TF.avoiddegrad.mask)))
##handle overlaps
TF.avoiddegrad.mask[TF.avoiddegrad.mask>1] <- 1
#plot(TF.avoiddegrad.mask)
writeRaster(TF.avoiddegrad.mask, "rasters/STM/all_forest_mask/STM_2020_avoiddegrad.tif", format = "GTiff", overwrite = T)


#


# scenario avoid deforestation
UPF.avoiddefor.mask <- UPF.avoiddefor

TF.avoiddefor.mask <- sum(UPF.avoiddefor, DPF2020.mask, SF2010.mask, na.rm = T)
##cheking
#sort(unique(values(TF.avoiddefor.mask)))
##handle overlaps
TF.avoiddefor.mask[TF.avoiddefor.mask>1] <- 1
#plot(TF.avoiddefor.mask)
writeRaster(TF.avoiddefor.mask, "rasters/STM/all_forest_mask/STM_2020_avoiddeforest.tif", format = "GTiff", overwrite = T)


#


# scenario avoid both
UPF.avoidboth.mask <- UPF.avoidboth

TF.avoidboth.mask <- sum(UPF.avoidboth.mask, DPF2010.mask, SF2010.mask, na.rm = T)
##cheking
#sort(unique(values(TF.avoidboth.mask)))
##handle overlaps
TF.avoidboth.mask[TF.avoidboth.mask>1] <- 1
#plot(TF.avoidboth.mask)
writeRaster(TF.avoidboth.mask, "rasters/STM/all_forest_mask/STM_2020_avoidboth.tif", format = "GTiff", overwrite = T)


#


# scenario restoration without avoid
SF2020.restore10.mask <- SF2020.restore10

TF.restore10.a <- sum(UPF2020.mask, DPF2020.mask, SF2020.restore10.mask, na.rm = T)
##cheking
#sort(unique(values(TF.restore10.a)))
TF.restore10.a[TF.restore10.a>1] <- 1
#plot(TF.restore10.a)
writeRaster(TF.restore10.a, "rasters/STM/all_forest_mask/STM_2020_restor_wo_avoid.tif", format = "GTiff", overwrite = T)

#


# scenario restoration and avoid deforestation
SF2010.restore10.mask <- SF2010.restore10

TF.restore10.b <- sum(UPF.avoiddefor.mask, DPF2020.mask, SF2010.restore10.mask, na.rm = T)
##cheking
#sort(unique(values(TF.restore10.b)))
TF.restore10.b[TF.restore10.b>1] <- 1
#plot(TF.restore10.b)
writeRaster(TF.restore10.b, "rasters/STM/all_forest_mask/STM_2020_restor_n_avoid_deforest.tif", format = "GTiff", overwrite = T)

#


# scenario restoration and avoid both
TF.restore10.c <- sum(UPF.avoidboth.mask, DPF2010.mask, SF2010.restore10.mask, na.rm = T)
##cheking
#sort(unique(values(TF.restore10.c)))
TF.restore10.c[TF.restore10.c>1] <- 1
#plot(TF.restore10.c)
writeRaster(TF.restore10.c, "rasters/STM/all_forest_mask/STM_2020_restor_n_avoid_both.tif", format = "GTiff", overwrite = T)

#


##############################################################################################################################################################################################################################################

## detecting multicollinearity between exploratory variables
#env.explanatory.var.list <- list.files("rasters/STM/2010_real", pattern = ".tif", full.names = T, recursive = T)
#
#env.explanatory.var <- stack(env.explanatory.var.list)
#names(env.explanatory.var) <- unlist(strsplit(env.explanatory.var.list, "/|.tif"))[seq(4,80,4)]
###cheking
##env.explanatory.var
##plot(env.explanatory.var[[1:10]], nc=2)
##plot(env.explanatory.var[[11:20]], nc=2)
#
## visual inspection of aggregation using removeCollinearity() function from virtualspecies package
#correlated.var <- removeCollinearity(env.explanatory.var, multicollinearity.cutoff = 0.7, sample.points = T, nb.points = 999999, method = "pearson", plot = T)
###cheking
##correlated.var
#
## variation inflation factor
#inflated.var <- vifcor(env.explanatory.var, th = 0.7, maxobservations = 999999)
###cheking
##inflated.var@results
##inflated.var@excluded
#
#sel.var.df <- data.frame(rbind(cbind(VAR=inflated.var@results$Variables, VIF=inflated.var@results$VIF),
#                               cbind(VAR=inflated.var@excluded, VIF=NA)))
#
#write.csv(sel.var.df, paste0("rasters/STM/selected_environmental_explanatory_variables_byVIF.csv", sep=""), row.names = F)
#
###cheking
##env.explanatory.var <- env.explanatory.var[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
##plot(env.explanatory.var)
#
#
##
#
#rm(list=ls()) #keeping only raster stack
#gc()
#
#
#
##


