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
library(datazoom.amazonia)
library(stars)
library(lwgeom)

# importing raw rasters ===================================
## standard projection
std.proj <- "+proj=longlat +datum=WGS84 +units=m +no_defs"

## shapefile paragominas
pgm.shp <- readOGR(dsn = "shapes", layer = "Paragominas_Mask_R3")
proj4string(pgm.shp) <- CRS(std.proj)
pgm.shp <- spTransform(pgm.shp, crs(std.proj))
pgm.shp2 <- gBuffer(pgm.shp, width = 0.02700003)

### 100m resolution raster
#pgm.lulc.100 <- raster("rasters/PGM/raw/pgm-lulc-mapbiomas-brazil-collection-80-2022-100res.tif")
#
#



## land use land cover from mapbiomas collection 7, 30m res [2010 and 2020]
# [PGM] paragominas
pgm.lulc <- stack(c("rasters/PGM/raw/pgm-lulc-mapbiomas-brazil-collection-80-2010.tif",
                    "rasters/PGM/raw/pgm-lulc-mapbiomas-brazil-collection-80-2020.tif"))
names(pgm.lulc) <- c("pgm.lulc.2010real", "pgm.lulc.2020real")


# isolating forest class pixels
pgm.lulc.2010.forest.class <- pgm.lulc[["pgm.lulc.2010real"]]
pgm.lulc.2010.forest.class[pgm.lulc.2010.forest.class==3] <- 1
pgm.lulc.2010.forest.class[pgm.lulc.2010.forest.class>1] <- 0


pgm.lulc.2010.forest.mask <- pgm.lulc.2010.forest.class
pgm.lulc.2010.forest.mask[pgm.lulc.2010.forest.mask==0] <- NA


pgm.lulc.2020.forest.class <- pgm.lulc[["pgm.lulc.2020real"]]
pgm.lulc.2020.forest.class[pgm.lulc.2020.forest.class==3] <- 1
pgm.lulc.2020.forest.class[pgm.lulc.2020.forest.class>1] <- 0


pgm.lulc.2020.forest.mask <- pgm.lulc.2020.forest.class
pgm.lulc.2020.forest.mask[pgm.lulc.2020.forest.mask==0] <- NA
#
#





# secondary forest age from Silva Jr. et al 2020  [2010 and 2020]
# [DOI: 10.1038/s41597-020-00600-4]
# [PGM] paragominas
pgm.sfage <- stack(c("rasters/PGM/raw/pgm-2010-sfage-mapbiomas-brazil-collection-60.tif",
                     "rasters/PGM/raw/pgm-2020-sfage-mapbiomas-brazil-collection-60.tif"))
names(pgm.sfage) <- c("pgm.sfage.2010real", "pgm.sfage.2020real")

#checking
#st_crs(pgm.sfage)==st_crs(pgm.shp)
#plot(pgm.sfage)
#range(values(pgm.sfage[["pgm.sfage.2010real"]]), na.rm = T)

# Conversion of rasters into same extent
#pgm.sfage <- resample(pgm.sfage, pgm.lulc.100, method='ngb')

#excluding non-forest areas
pgm.sfage[["pgm.sfage.2010real"]] <- mask(pgm.sfage[["pgm.sfage.2010real"]], pgm.lulc.2010.forest.mask)
pgm.sfage[["pgm.sfage.2010real"]][is.na(pgm.sfage[["pgm.sfage.2010real"]])] <- 0


pgm.sfage[["pgm.sfage.2020real"]] <- mask(pgm.sfage[["pgm.sfage.2020real"]], pgm.lulc.2020.forest.mask)
pgm.sfage[["pgm.sfage.2020real"]][is.na(pgm.sfage[["pgm.sfage.2020real"]])] <- 0



# isolating secondary forest class pixels
pgm.sfage.2010.all.class <- pgm.sfage[["pgm.sfage.2010real"]]
pgm.sfage.2010.all.class[pgm.sfage.2010.all.class>0] <- 1
pgm.sfage.2010.all.class[pgm.sfage.2010.all.class<1] <- 0

pgm.sfage.2010.mask <- pgm.sfage.2010.all.class
pgm.sfage.2010.mask[pgm.sfage.2010.mask==0] <- NA

pgm.sfage.2020.all.class <- pgm.sfage[["pgm.sfage.2020real"]]
pgm.sfage.2020.all.class[pgm.sfage.2020.all.class>0] <- 1
pgm.sfage.2020.all.class[pgm.sfage.2020.all.class<1] <- 0

pgm.sfage.2020.mask <- pgm.sfage.2020.all.class
pgm.sfage.2020.mask[pgm.sfage.2020.mask==0] <- NA


#
#


# gathering time series degradation based on mapbiomas fire, degrad and deter ====================
#' @description mapbiomas fire is a time-series data from 1985 to present
#' degrad is the first monitoring system to detect degradation in brazil
#' it was functional from 2007 to 2016 but does not differentiate between classes of degradation.
#' deter is the current monitoring system, with data from 2016

## import, crop and save degrad data to study area

for (i in 2007:2016) {
  
  degrad.yearx <- readOGR(dsn = "C:/Users/miral/Dropbox/GIS/degradacao", layer = paste0("Degrad", i, "_Final_pol"))
  pgm.degrad.yearx <- crop(degrad.yearx, extent(pgm.shp2))
  
  pgm.degrad.yearx.raster <- rasterize(pgm.degrad.yearx, pgm.lulc.2010.forest.mask, field=1)
  pgm.degrad.yearx.raster <- mask(pgm.degrad.yearx.raster, pgm.shp2)
  pgm.degrad.yearx.raster[is.na(pgm.degrad.yearx.raster)]<-0
  
  writeRaster(pgm.degrad.yearx.raster, 
              paste0("C:/Users/miral/Dropbox/MAPBIOMAS-EXPORT/fogo-time-series/PGM/degrad-brazil-inpe-pgm-", i, ".tif"), format="GTiff", overwrite=T)
  
  cat("\n> year", i, "done! <\n")
}



## import, crop and save deter data to study area
### selecting only fire scars and selective logging classes

deter <- readOGR(dsn = "C:/Users/miral/Dropbox/GIS/degradacao", layer = "deter-amz-deter-public")
pgm.deter <- crop(deter, extent(pgm.shp2))
#unique(pgm.deter@data$CLASSNAME)
degradation_cat <- c('CICATRIZ_DE_QUEIMADA', 'CS_DESORDENADO', 'CS_GEOMETRICO') #'DEGRADACAO'
pgm.deter <- pgm.deter[pgm.deter$CLASSNAME %in% degradation_cat,]

for (i in 2016:2020) {
  
  pgm.deter.yearx <- pgm.deter[grep(i, pgm.deter$VIEW_DATE),]
  
  pgm.deter.yearx.raster <- rasterize(pgm.deter.yearx, pgm.lulc.2010.forest.mask, field=1)
  pgm.deter.yearx.raster <- mask(pgm.deter.yearx.raster, pgm.shp2)
  pgm.deter.yearx.raster[is.na(pgm.deter.yearx.raster)]<-0
  
  writeRaster(pgm.deter.yearx.raster, 
              paste0("C:/Users/miral/Dropbox/MAPBIOMAS-EXPORT/fogo-time-series/PGM/deter-brazil-inpe-pgm-", i, ".tif"), format="GTiff", overwrite=T)
  
  cat("\n> year", i, "done! <\n")
}

rm(list=ls()[ls() %in% c("degrad.yearx", "pgm.degrad.yearx", "pgm.degrad.yearx.raster",
                         "deter", "pgm.deter", "degradation_cat", "pgm.deter.yearx",
                         "pgm.deter.yearx.raster", "i")])
gc()
#
#





## import fire time series from mapbiomas
pgm.raw.fire.list <- list.files("C:/Users/miral/Dropbox/MAPBIOMAS-EXPORT/fogo-time-series/PGM",
                                pattern = "mapbiomas", full.names = T, recursive = T)


pgm.raw.fire <- raster::stack(pgm.raw.fire.list)

## import degradation time series from degrad [inpe 2007-2016]
pgm.degrad.list <- list.files("C:/Users/miral/Dropbox/MAPBIOMAS-EXPORT/fogo-time-series/PGM",
                                pattern = "degrad", full.names = T, recursive = T)


pgm.degrad <- raster::stack(pgm.degrad.list)

## import degradation time series from deter-b [inpe 201-2020]
pgm.deter.list <- list.files("C:/Users/miral/Dropbox/MAPBIOMAS-EXPORT/fogo-time-series/PGM",
                              pattern = "deter", full.names = T, recursive = T)


pgm.deter <- raster::stack(pgm.deter.list)






# building time since degradation =========================================
##1985-2006 (only fire)
#count for multiple degradation
pgm.freq.degrad <- pgm.raw.fire[[1]]

#assign time since degradation
pgm.time.since.degrad <- pgm.raw.fire[[1]]
pgm.time.since.degrad[pgm.time.since.degrad!=1] <- 0
pgm.time.since.degrad[pgm.time.since.degrad==1] <- 0

for (i in 1:22) {
  
  #selecting year-by-year
  degrad.yearx <- pgm.raw.fire[[i]]
  
  #count for multiple fires
  pgm.freq.degrad <- sum(pgm.freq.degrad, degrad.yearx)
  
  #assign time since fire
  degrad.yearx[degrad.yearx!=1] <- NA
  degrad.yearx[degrad.yearx==1] <- 1000
  degrad.yearx[is.na(degrad.yearx[])] <- 0
  
  #counting time since degradation
  pgm.time.since.degrad <- pgm.time.since.degrad+1
  pgm.time.since.degrad <- sum(pgm.time.since.degrad, degrad.yearx, na.rm = T)
  pgm.time.since.degrad[pgm.time.since.degrad[]>1000] <- 0  
  
  cat("\n> year", 1984+i, "done! <\n")
}


##2007-2010 
#combining fire and degradation for primary forests
#only fire for secondary forests
for (i in 23:26) {
  
  #selecting year-by-year
    #fire
    fire.yearx <- pgm.raw.fire[[i]]
    
    #degradation
    inpe.yearx <- pgm.degrad[[i-22]]
    
    
    degrad.yearx <- sum(inpe.yearx, fire.yearx)
    degrad.yearx[degrad.yearx>1] <- 1
    
  #count for multiple degradation events
  pgm.freq.fire <- sum(pgm.freq.degrad, fire.yearx)
  pgm.freq.degrad <- sum(pgm.freq.degrad, degrad.yearx)
  
  #assign time since degradation
  fire.yearx[fire.yearx!=1] <- NA
  fire.yearx[fire.yearx==1] <- 1000
  fire.yearx[is.na(fire.yearx[])] <- 0
  
  degrad.yearx[degrad.yearx!=1] <- NA
  degrad.yearx[degrad.yearx==1] <- 1000
  degrad.yearx[is.na(degrad.yearx[])] <- 0
  
  #counting time since degradation
  pgm.time.since.degrad.sf <- pgm.time.since.degrad+1
  pgm.time.since.degrad.sf <- sum(pgm.time.since.degrad.sf, fire.yearx, na.rm = T)
  pgm.time.since.degrad.sf[pgm.time.since.degrad.sf[]>1000] <- 0  
  
  pgm.time.since.degrad <- pgm.time.since.degrad+1
  pgm.time.since.degrad <- sum(pgm.time.since.degrad, degrad.yearx, na.rm = T)
  pgm.time.since.degrad[pgm.time.since.degrad[]>1000] <- 0  
  
 cat("\n> year", 1984+i, "done! <\n")
}

pgm.2010.time.since.degrad.sf <- pgm.time.since.degrad.sf
pgm.2010.time.since.degrad.sf <-  mask(pgm.2010.time.since.degrad.sf, pgm.sfage.2010.mask)
writeRaster(pgm.2010.time.since.degrad.sf, "rasters/PGM/raw/pgm-2010-dsf-tsince0.tif", format = "GTiff", overwrite = T)

pgm.2010.time.since.degrad <- pgm.time.since.degrad
pgm.2010.time.since.degrad <-  mask(pgm.2010.time.since.degrad, pgm.lulc.2010.forest.mask)
pgm.2010.time.since.degrad <-  mask(pgm.2010.time.since.degrad, pgm.sfage.2010.mask, inverse = T)
writeRaster(pgm.2010.time.since.degrad, "rasters/PGM/raw/pgm-2010-dpf-tsince0.tif", format = "GTiff", overwrite = T)


pgm.1985.2010.freq.degrad <- pgm.freq.degrad
pgm.1985.2010.freq.degrad <-  mask(pgm.1985.2010.freq.degrad, pgm.lulc.2010.forest.mask)
writeRaster(pgm.1985.2010.freq.degrad, "rasters/PGM/raw/pgm-degfreq-1985_2010.tif", format = "GTiff", overwrite = T)

##2011-2015
for (i in 27:31) {
  
  #selecting year-by-year
  #fire
  fire.yearx <- pgm.raw.fire[[i]]
  
  #degradation
  inpe.yearx1 <- pgm.degrad[[i-22]]
  
  
  degrad.yearx <- sum(inpe.yearx, fire.yearx)
  degrad.yearx[degrad.yearx>1] <- 1
  
  #count for multiple degradation events
  pgm.freq.fire <- sum(pgm.freq.fire, fire.yearx)
  pgm.freq.degrad <- sum(pgm.freq.degrad, degrad.yearx)
  
  #assign time since degradation
  fire.yearx[fire.yearx!=1] <- NA
  fire.yearx[fire.yearx==1] <- 1000
  fire.yearx[is.na(fire.yearx[])] <- 0
  
  degrad.yearx[degrad.yearx!=1] <- NA
  degrad.yearx[degrad.yearx==1] <- 1000
  degrad.yearx[is.na(degrad.yearx[])] <- 0
  
  #counting time since degradation
  pgm.time.since.degrad.sf <- pgm.time.since.degrad.sf+1
  pgm.time.since.degrad.sf <- sum(pgm.time.since.degrad.sf, fire.yearx, na.rm = T)
  pgm.time.since.degrad.sf[pgm.time.since.degrad.sf[]>1000] <- 0  
  
  pgm.time.since.degrad <- pgm.time.since.degrad+1
  pgm.time.since.degrad <- sum(pgm.time.since.degrad, degrad.yearx, na.rm = T)
  pgm.time.since.degrad[pgm.time.since.degrad[]>1000] <- 0  
  
  cat("\n> year", 1984+i, "done! <\n")
}

##2016-2020
for (i in 32:36) {
  
  #selecting year-by-year
  #fire
  fire.yearx <- pgm.raw.fire[[i]]
  
  #degradation
  inpe.yearx1 <- pgm.deter[[i-31]]
  
  
  degrad.yearx <- sum(inpe.yearx, fire.yearx)
  degrad.yearx[degrad.yearx>1] <- 1
  
  #count for multiple degradation events
  pgm.freq.fire <- sum(pgm.freq.fire, fire.yearx)
  pgm.freq.degrad <- sum(pgm.freq.degrad, degrad.yearx)
  
  #assign time since degradation
  fire.yearx[fire.yearx!=1] <- NA
  fire.yearx[fire.yearx==1] <- 1000
  fire.yearx[is.na(fire.yearx[])] <- 0
  
  degrad.yearx[degrad.yearx!=1] <- NA
  degrad.yearx[degrad.yearx==1] <- 1000
  degrad.yearx[is.na(degrad.yearx[])] <- 0
  
  #counting time since degradation
  pgm.time.since.degrad.sf <- pgm.time.since.degrad.sf+1
  pgm.time.since.degrad.sf <- sum(pgm.time.since.degrad.sf, fire.yearx, na.rm = T)
  pgm.time.since.degrad.sf[pgm.time.since.degrad.sf[]>1000] <- 0  
  
  pgm.time.since.degrad <- pgm.time.since.degrad+1
  pgm.time.since.degrad <- sum(pgm.time.since.degrad, degrad.yearx, na.rm = T)
  pgm.time.since.degrad[pgm.time.since.degrad[]>1000] <- 0  
  
  cat("\n> year", 1984+i, "done! <\n")
}

pgm.2020.time.since.degrad.sf <- pgm.time.since.degrad.sf
pgm.2020.time.since.degrad.sf <-  mask(pgm.2020.time.since.degrad.sf, pgm.sfage.2020.mask)
writeRaster(pgm.2020.time.since.degrad.sf, "rasters/PGM/raw/pgm-2020-dsf-tsince0.tif", format = "GTiff", overwrite = T)

pgm.2020.time.since.degrad <- pgm.time.since.degrad
pgm.2020.time.since.degrad <-  mask(pgm.2020.time.since.degrad, pgm.lulc.2020.forest.mask)
pgm.2020.time.since.degrad <-  mask(pgm.2020.time.since.degrad, pgm.sfage.2020.mask, inverse = T)
writeRaster(pgm.2020.time.since.degrad, "rasters/PGM/raw/pgm-2020-dpf-tsince0.tif", format = "GTiff", overwrite = T)


pgm.1985.2020.freq.degrad <- pgm.freq.degrad
pgm.1985.2020.freq.degrad <-  mask(pgm.1985.2020.freq.degrad, pgm.lulc.2020.forest.mask)
writeRaster(pgm.1985.2020.freq.degrad, "rasters/PGM/raw/pgm-degfreq-1985_2020.tif", format = "GTiff", overwrite = T)

rm(list=ls()[ls() %in% c("pgm.raw.fire.list", "pgm.raw.fire", "pgm.degrad.list", "pgm.degrad",
                         "pgm.deter.list", "pgm.deter", "degrad.yearx", "fire.yearx", "inpe.yearx", "i")])
gc()
#
#





# candidate areas for restoration scenarios ==============|
#' @description 
#' 
# land use land cover from mapbiomas collection 8 [2007, 2010 and 2022]
# 30m resolution raster
pgm.lulc <- stack(c("rasters/PGM/raw/pgm-lulc-mapbiomas-brazil-collection-80-2007.tif",
                    "rasters/PGM/raw/pgm-lulc-mapbiomas-brazil-collection-80-2010.tif"))
names(pgm.lulc) <- c("pgm.lulc.2007real", "pgm.lulc.2010real", "pgm.lulc.2022real")



# isolating forest class pixels
pgm.lulc.2007.forest.class <- pgm.lulc[["pgm.lulc.2007real"]]
pgm.lulc.2007.forest.class[pgm.lulc.2007.forest.class==3] <- 1
pgm.lulc.2007.forest.class[pgm.lulc.2007.forest.class>1] <- 0


pgm.lulc.2007.forest.mask <- pgm.lulc.2007.forest.class
pgm.lulc.2007.forest.mask[pgm.lulc.2007.forest.mask==0] <- NA


pgm.lulc.2010.forest.class <- pgm.lulc[["pgm.lulc.2010real"]]
pgm.lulc.2010.forest.class[pgm.lulc.2010.forest.class==3] <- 1
pgm.lulc.2010.forest.class[pgm.lulc.2010.forest.class>1] <- 0


pgm.lulc.2010.forest.mask <- pgm.lulc.2010.forest.class
pgm.lulc.2010.forest.mask[pgm.lulc.2010.forest.mask==0] <- NA
#
#



#deforestation class pixels (15 == pasture & 19 == soybean)
deforestation.class.list <- c(15,35,39,41)

converted.area.2007 <- pgm.lulc[["pgm.lulc.2007real"]]

converted.area.2007[converted.area.2007[] %in% deforestation.class.list] = 1
converted.area.2007[converted.area.2007[] > 1] = 0
names(converted.area.2007) <- "converted.area.2007"
#plot(converted.area.2007, col=c("#ffffff","brown"), legend = F)
#length(converted.area.2007[converted.area.2007[]==1])


converted.area.2010 <- pgm.lulc[["pgm.lulc.2010real"]]

converted.area.2010[converted.area.2010[] %in% deforestation.class.list] = 1
converted.area.2010[converted.area.2010[] > 1] = 0
names(converted.area.2010) <- "converted.area.2010"
#plot(converted.area.2010, col=c("#ffffff","brown"), legend = F)
#length(converted.area.2010[converted.area.2010[]==1])
converted.area.2010.mask <- converted.area.2010
converted.area.2010.mask[converted.area.2010.mask[] != 1] <- NA
#
#



#identifying consolidated areas (areas deforested before 2008)
consolidated.areas <- sum(converted.area.2007, converted.area.2010, na.rm = T)
consolidated.areas[consolidated.areas[] < 2] <- 0
consolidated.areas[consolidated.areas[] == 2] <- 1
names(consolidated.areas) <- "consolidated.areas"
#plot(consolidated.areas, col=c("#ffffff","brown"), legend = F)
#length(consolidated.areas[consolidated.areas[]==1])
consolidated.areas.mask <- consolidated.areas
consolidated.areas.mask[consolidated.areas.mask[] != 1] <- NA
#
#



# Rivers
# source1: ANA -  https://metadados.snirh.gov.br/geonetwork/srv/api/records/0f57c8a0-6a0f-4283-8ce3-114ba904b9fe
# source2: HydroSHEDS - https://www.hydrosheds.org/hydroatlas
# Permanent Protection Areas [APPs]
# According to Brazilian Forest Code (Law n. 12.651/2012)
# In general[*]:
# river width <10m, APP == 30
# river width 10-50m, APP == 50m
# river width 50-200m, APP == 100m
# river width 200-600m, APP == 200m
# river width +600m, APP == 500m
#
# [*] there are some exceptions, according to property size, 
# and if the area were converted before 2008 (see the law for details)
#

pgm.river <- readOGR(dsn = "rasters/PGM/raw", layer = "pgm_RiverATLAS_v10")
#head(pgm.river@data)
pgm.river@data <- pgm.river@data %>% dplyr::select(HYRIV_ID:LENGTH_KM, ria_ha_csu, ria_ha_usu) %>% 
  mutate(
    ril_m = LENGTH_KM * 1000, #converting to meters
    ria_m2 = ria_ha_csu * 10000, #converting to square meters
    riw_m = (ria_m2 / ril_m), #estimating the mean river segment width
    riapp = ifelse(riw_m<10, 30, #applying forest code rules
                   ifelse(riw_m>=10 & riw_m<50, 50,
                          ifelse(riw_m >=50 & riw_m<200, 100,
                                 ifelse(riw_m>=200 & riw_m<600, 200, 500)))), 
    buffer = (riw_m + (2*riapp))/111111 #buffering APP area on both river margins and converting to decimal
  )


pgm.appList <- vector("list", length(pgm.river))
for (i in 1:length(pgm.river)) {
  a <- gBuffer(pgm.river[i,], width = pgm.river$buffer[i])
  a$id = pgm.river$HYRIV_ID[i]
  pgm.appList[[i]] <- a
}

pgm.app <- do.call("rbind", pgm.appList)

#checking
#st_crs(pgm.app)==st_crs(pgm.app)
#plot(pgm.lulc[["pgm.lulc.2010real"]])
#plot(pgm.river, add=T)
#plot(pgm.app, add=T)

#rm(list=ls()[ls() %in% c("pgm.appList", "a", "i")])
#gc()



#
#



#import rural properties shapefiles and data from SISCAR
#https://www.car.gov.br/publico/municipios/downloads
pgm.car <- readOGR(dsn = "rasters/PGM/raw/SHAPE_1505502_CAR_Paragominas", layer = "AREA_IMOVEL")
pgm.car <- spTransform(pgm.car, crs(std.proj))

#checking
#st_crs(pgm.car)==st_crs(pgm.shp)


#creating a copy
pgm.car.copy <- pgm.car
#pgm.car <- pgm.car.copy

#reducing overlap
#properties have the same CAR code
pgm.car <- pgm.car[!duplicated(pgm.car@data$COD_IMOVEL),]
#properties occupy exactly the same area
coord_poly <- lapply(pgm.car@polygons, function(x){lapply(x@Polygons, function(x){coordinates(x)})}) 
pgm.car <- pgm.car[!duplicated(coord_poly),]


#
pgm.car@data$keep <- NA

combos <- combn(nrow(pgm.car),2)


suppressMessages(
  for(k in seq_along(combos[1,])){
    i <- combos[1,k]
    j <- combos[2,k]
    #print(paste("intersecting",i,j))
    
    si <- pgm.car[i,c(1,2,6)]
    sj <- pgm.car[j,c(1,2,6)]
    names(sj) <- c("COD_IMOVEL.1","NUM_AREA.1","TIPO_IMOVE.1")
    
    #if (nrow(st_intersection(st_as_sf(pgm.car[i,]), st_as_sf(pgm.car[j,]))) == 0) next
    
    teste <- try(intersect(si, sj), TRUE)
    
    if (class(teste)=="try-error") next
    
    teste$intersect_area <- area(teste)*0.0001
    teste$coverage <- teste$intersect_area/teste$NUM_AREA
    teste$coverage1 <- teste$intersect_area/teste$NUM_AREA.1
    
    
    #teste <- st_intersection(st_as_sf(pgm.car[i,]), st_as_sf(pgm.car[j,])) %>% # shape of overlap area
    #              mutate(intersect_area = st_area(.)*0.0001, # create new column with shape area
    #                     coverage = as.numeric(intersect_area/NUM_AREA), # calculate coverage
    #                     coverage1 = as.numeric(intersect_area/NUM_AREA.1)) %>%
    #              st_drop_geometry()
    
    #properties with an overlap larger than 30% with protected area or agrarian reform settlements were excluded;
    if(teste$TIPO_IMOVE=="IRU" & teste$TIPO_IMOVE.1!="IRU" & teste$coverage > .3) {
      print(paste0("excluding ", pgm.car@data[pgm.car@data$COD_IMOVEL==teste$COD_IMOVEL,"COD_IMOVEL"], " due to overlap with areas not suitable for registry: ", pgm.car@data[pgm.car@data$COD_IMOVEL==teste$COD_IMOVEL.1,"COD_IMOVEL"]))
      pgm.car@data[pgm.car@data$COD_IMOVEL==teste$COD_IMOVEL,"keep"] <- F
    }
    
    
    #where overlap was greater than 80%, the smallest property is excluded;
    if(teste$TIPO_IMOVE=="IRU" & teste$TIPO_IMOVE.1=="IRU" & teste$coverage > .8 & teste$NUM_AREA<teste$NUM_AREA.1){ 
      print(paste0("excluding ", pgm.car@data[pgm.car@data$COD_IMOVEL==teste$COD_IMOVEL,"COD_IMOVEL"], " due to overlap area is greater than 80%: ", pgm.car@data[pgm.car@data$COD_IMOVEL==teste$COD_IMOVEL.1,"COD_IMOVEL"]))
      pgm.car@data[pgm.car@data$COD_IMOVEL==teste$COD_IMOVEL,"keep"] <- F
    }
    
    if(teste$TIPO_IMOVE=="IRU" & teste$TIPO_IMOVE.1=="IRU" & teste$coverage1 > .8 & teste$NUM_AREA.1<teste$NUM_AREA){ 
      print(paste0("excluding ", pgm.car@data[pgm.car@data$COD_IMOVEL==teste$COD_IMOVEL,"COD_IMOVEL"], " due to overlap area is greater than 80%: ", pgm.car@data[pgm.car@data$COD_IMOVEL==teste$COD_IMOVEL,"COD_IMOVEL"]))
      pgm.car@data[pgm.car@data$COD_IMOVEL==teste$COD_IMOVEL.1,"keep"] <- F
    }
    
    
  }
)


#writeOGR(pgm.car, dsn="rasters/PGM", layer = "pgm.car.ed2", driver = "ESRI Shapefile")
#pgm.car <- readOGR(dsn = "rasters/PGM", layer = "pgm.car.ed2")

View(pgm.car@data)
nrow(pgm.car@data[is.na(pgm.car@data$keep),])

pgm.car <- pgm.car[is.na(pgm.car@data$keep),]

#writeOGR(pgm.car, dsn="rasters/PGM", layer = "pgm.car.ed3", driver = "ESRI Shapefile")
#pgm.car <- readOGR(dsn = "rasters/PGM", layer = "pgm.car.ed3")


#recalculating the size and number of fiscal modules and
#comparing property size between the value declared by the owner on attribute table 
#and the value calculated based on the polygon -- 0 == there is no significant difference; 1 == different by more than 3ha
pgm.car@data <- pgm.car@data %>% mutate(num_area_ha = raster::area(pgm.car)*0.0001,
                                        num_modulo_new = num_area_ha/55,
                                        num_area_flag = ifelse(abs(NUM_AREA-num_area_ha) < 3, 0, 1)) %>% 
                           dplyr::select(-keep)

#View(pgm.car@data)
#table(pgm.car@data$num_area_flag)

#checking how many properties have less than 30x30m
nrow(pgm.car@data[pgm.car@data$num_area_ha<.1,])
#excluding
pgm.car <- pgm.car[which(pgm.car@data$num_area_ha>=.1),]


#checking
#st_crs(pgm.car)==st_crs(pgm.shp)
#plot(pgm.lulc[["pgm.lulc.2010real"]])
#plot(pgm.car, add=T)


#calculating forest cover (ha) in each property
#adding variable for forest cover
pgm.car@data$FOREST_COVER_2007 <- NA
#creating layer with forest class
forest.class.2007 <- pgm.lulc.2007.forest.mask

j=nrow(pgm.car@data)
n=0
for (i in pgm.car$COD_IMOVEL) {
  
  rural.property <- pgm.car[pgm.car$COD_IMOVEL==i,]
  
  forest.cover.2007 <- crop(forest.class.2007, extent(rural.property))
  forest.cover.2007 <- mask(forest.cover.2007, rural.property)
  
  if(all(is.na(forest.cover.2007[])))
  {
    
    pgm.car[pgm.car$COD_IMOVEL==i,"FOREST_COVER_2007"] <- 0
    j=j-1
    n=n+1
    cat("\n>there are no forests in property", i, "in 2007 <\n")
    cat("\n>", j, "out of", nrow(pgm.car@data), "properties left<\n")
    
  } 
  else
  {
    
    pgm.car[pgm.car$COD_IMOVEL==i,"FOREST_COVER_2007"] <- tapply(raster::area(forest.cover.2007), forest.cover.2007[], sum, na.rm=T)*100
    j=j-1
    cat("\n>", j, "out of", nrow(pgm.car@data), "properties left<\n")
    
  }
  
}

#calculating the percentage of forest cover in relation to property size
pgm.car@data$FOREST_COVER_2007_PP <- pgm.car@data$FOREST_COVER_2007/pgm.car@data$num_area_ha

#View(pgm.car@data)


#calculating forest cover (ha) in each property
#adding variable for forest cover
pgm.car@data$FOREST_COVER_2010 <- NA
#creating layer with forest class
forest.class.2010 <- pgm.lulc.2010.forest.mask

j=nrow(pgm.car@data)
n=0
for (i in pgm.car$COD_IMOVEL) {
  
  rural.property <- pgm.car[pgm.car$COD_IMOVEL==i,]
  
  forest.cover.2010 <- crop(forest.class.2010, extent(rural.property))
  forest.cover.2010 <- mask(forest.cover.2010, rural.property)
  
  if(all(is.na(forest.cover.2010[])))
  {
    
    pgm.car[pgm.car$COD_IMOVEL==i,"FOREST_COVER_2010"] <- 0
    j=j-1
    n=n+1
    cat("\n>there are no forests in property", i, "in 2010 <\n")
    cat("\n>", j, "out of", nrow(pgm.car@data), "properties left<\n")
    
  } 
  else
  {
    
    pgm.car[pgm.car$COD_IMOVEL==i,"FOREST_COVER_2010"] <- tapply(raster::area(forest.cover.2010), forest.cover.2010[], sum, na.rm=T)*100
    j=j-1
    cat("\n>", j, "out of", nrow(pgm.car@data), "properties left<\n")
    
  }
  
}

#calculating the percentage of forest cover in relation to property size
pgm.car@data$FOREST_COVER_2010_PP <- pgm.car@data$FOREST_COVER_2010/pgm.car@data$num_area_ha

#View(pgm.car@data)

#checking
#anyNA(pgm.car@data$FOREST_COVER_2010_PP)
#length(which(is.na(pgm.car$FOREST_COVER_2010_PP)))
#pgm.car.restoration.candidates <- pgm.car[!is.na(pgm.car$FOREST_COVER_2010_PP),]


#creating a copy
pgm.car.restoration.candidates <- pgm.car


#select properties based on threshold
#if the area is big/medium properties with less than 50%
pgm.car.restoration.candidates@data$NEED_INCREMENT <- ifelse(pgm.car.restoration.candidates@data$num_modulo_new > 4 & 
                                                               pgm.car.restoration.candidates@data$FOREST_COVER_2010_PP < .5, 1, 0)

#if the area is small properties with less the 2007 forest cover
pgm.car.restoration.candidates@data$NEED_INCREMENT <- ifelse(pgm.car.restoration.candidates$num_modulo_new <= 4 &
                                                               pgm.car.restoration.candidates$FOREST_COVER_2007_PP >= .5 &
                                                               pgm.car.restoration.candidates$FOREST_COVER_2010_PP < .5, 1, pgm.car.restoration.candidates@data$NEED_INCREMENT)

pgm.car.restoration.candidates@data$NEED_INCREMENT <- ifelse(pgm.car.restoration.candidates$num_modulo_new <= 4 &
                                                               pgm.car.restoration.candidates$FOREST_COVER_2007_PP < .5 &
                                                               pgm.car.restoration.candidates$FOREST_COVER_2010_PP < pgm.car.restoration.candidates$FOREST_COVER_2007_PP, 1, pgm.car.restoration.candidates@data$NEED_INCREMENT)

pgm.car.restoration.candidates@data$NEED_INCREMENT <- ifelse(pgm.car.restoration.candidates$num_modulo_new <= 4 &
                                                               pgm.car.restoration.candidates$FOREST_COVER_2007_PP < .5 &
                                                               pgm.car.restoration.candidates$FOREST_COVER_2010_PP >= pgm.car.restoration.candidates$FOREST_COVER_2007_PP, 0, pgm.car.restoration.candidates@data$NEED_INCREMENT)

#if the area is an agrarian reform settlements, 
#it was treated as small property
pgm.car.restoration.candidates@data$NEED_INCREMENT <- ifelse(pgm.car.restoration.candidates$TIPO_IMOVE=="AST" &
                                                               pgm.car.restoration.candidates$FOREST_COVER_2007_PP >= .5 &
                                                               pgm.car.restoration.candidates$FOREST_COVER_2010_PP < .5, 1, pgm.car.restoration.candidates@data$NEED_INCREMENT)

#View(pgm.car.restoration.candidates@data)
#table(pgm.car.restoration.candidates$NEED_INCREMENT)



#
#




#### candidate areas for restoration considering app ####

#select pixels based on APPs 
candidate.areas.water <- converted.area.2010
candidate.areas.water <- mask(candidate.areas.water, pgm.app)
#excluding consolidated areas
candidate.areas.water <- mask(candidate.areas.water, consolidated.areas.mask, inverse=T)
candidate.areas.water[candidate.areas.water[]!=1] <- NA
#plot(candidate.areas.water, col=c("#ffffff","brown"), legend = F)
#plot(pgm.shp, add=T)
#length(candidate.areas.water[candidate.areas.water[]==1])
#tapply(raster::area(candidate.areas.water), candidate.areas.water[], sum, na.rm=T)*100


#select remaining pixels based on converted areas
candidate.areas.rl <- converted.area.2010
#excluding consolidated areas
candidate.areas.rl <- mask(candidate.areas.rl, consolidated.areas.mask, inverse=T)
#excluding app areas
candidate.areas.rl <- mask(candidate.areas.rl, candidate.areas.water, inverse=T)
candidate.areas.rl[candidate.areas.rl[]!=1] <- NA
#plot(candidate.areas.rl, col=c("#ffffff","brown"), legend = F)
#plot(pgm.shp, add=T)
#length(candidate.areas.rl[candidate.areas.rl[]==1])
#tapply(raster::area(candidate.areas.rl), candidate.areas.rl[], sum, na.rm=T)*100
#
#




################################################################################
##select pixels based on slope                                                 #
##obs: there are no deforested areas in PGM with slope >=45o                   #
##slope <- terrain(elevation, opt = 'slope', unit = 'degrees', neighbors=8)    #
##values(slope)[values(slope) < 45] = NA                                       #
##values(slope)[values(slope) >= 45] = 1                                       #
###plot(slope)                                                                 #
##                                                                             #
##candidate.areas.slope <- candidate.areas.total                               #
##candidate.areas.slope <- mask(candidate.areas.slope, slope)                  #
##values(candidate.areas.slope)[is.na(values(candidate.areas.slope))] = 0      #
###plot(candidate.areas.slope, col="#ffffff")                                  #
################################################################################




#forest cover increment -- adding the app candidate areas to forest cover
pgm.car.restoration.candidates@data$APP_FOREST_COVER_INCREMENT <- NA
j=nrow(pgm.car.restoration.candidates@data[pgm.car.restoration.candidates@data$NEED_INCREMENT==1,])
n=0
#i="PA-1505502-39CCE4418D2D487F9AC0FD3045A374CF"
for (i in pgm.car.restoration.candidates$COD_IMOVEL) {
  
  if(pgm.car.restoration.candidates@data[pgm.car.restoration.candidates@data$COD_IMOVEL==i,"NEED_INCREMENT"] != 1) next
  
  rural.property <- pgm.car.restoration.candidates[pgm.car.restoration.candidates$COD_IMOVEL==i,]
  
  forest.cover.2010 <- crop(forest.class.2010, extent(rural.property))
  forest.cover.2010 <- mask(forest.cover.2010, rural.property)
  
  restored.cover <- crop(candidate.areas.water, extent(rural.property))
  restored.cover <- mask(restored.cover, rural.property)
  
  forest.cover.increment <- sum(forest.cover.2010, restored.cover, na.rm = T)
  forest.cover.increment[forest.cover.increment[]>1]<-1
  forest.cover.increment[forest.cover.increment[]!=1]<-NA
  
  if(all(is.na(forest.cover.increment[])))
  {
    
    pgm.car.restoration.candidates@data[pgm.car.restoration.candidates@data$COD_IMOVEL==i,"APP_FOREST_COVER_INCREMENT"] <- pgm.car.restoration.candidates@data[pgm.car.restoration.candidates@data$COD_IMOVEL==i,"FOREST_COVER_2010"]
    j=j-1
    n=n+1
    cat("\n>there are no forests in property", i, "<\n")
    cat("\n>", j, "out of", nrow(pgm.car.restoration.candidates@data[pgm.car.restoration.candidates@data$NEED_INCREMENT==1,]), "properties left<\n")
    
  } 
  else
  {
    
    pgm.car.restoration.candidates@data[pgm.car.restoration.candidates@data$COD_IMOVEL==i,"APP_FOREST_COVER_INCREMENT"] <- tapply(raster::area(forest.cover.increment), forest.cover.increment[], sum, na.rm=T)*100
    j=j-1
    cat("\n>", j, "out of", nrow(pgm.car.restoration.candidates@data[pgm.car.restoration.candidates@data$NEED_INCREMENT==1,]), "properties left<\n")
    
  }
  
}

#checking
#View(pgm.car.restoration.candidates@data)
#anyNA(pgm.car.restoration.candidates@data$APP_FOREST_COVER_INCREMENT)
#length(which(is.na(pgm.car.restoration.candidates@data$APP_FOREST_COVER_INCREMENT)))

pgm.car.restoration.candidates@data$APP_FOREST_COVER_INCREMENT <- ifelse(is.na(pgm.car.restoration.candidates@data$APP_FOREST_COVER_INCREMENT),
                                                                         pgm.car.restoration.candidates@data$FOREST_COVER_2010, pgm.car.restoration.candidates@data$APP_FOREST_COVER_INCREMENT)

#calculating the percentage of forest cover in relation to property size
pgm.car.restoration.candidates@data$APP_FOREST_COVER_INCREMENT_PP <- pgm.car.restoration.candidates@data$APP_FOREST_COVER_INCREMENT/pgm.car.restoration.candidates@data$num_area_ha




#select properties based on threshold
#if the area is big/medium properties with less than 50%
pgm.car.restoration.candidates@data$NEED_INCREMENT_AFTER_APP <- ifelse(pgm.car.restoration.candidates@data$num_modulo_new > 4 &
                                                                         pgm.car.restoration.candidates@data$APP_FOREST_COVER_INCREMENT_PP < .5, 1, 0)

#if the area is small properties with less the 2007 forest cover
pgm.car.restoration.candidates@data$NEED_INCREMENT_AFTER_APP <- ifelse(pgm.car.restoration.candidates$num_modulo_new <= 4 &
                                                                         pgm.car.restoration.candidates$FOREST_COVER_2007_PP >= .5 &
                                                                         pgm.car.restoration.candidates$APP_FOREST_COVER_INCREMENT_PP < .5, 1, pgm.car.restoration.candidates@data$NEED_INCREMENT_AFTER_APP)

pgm.car.restoration.candidates@data$NEED_INCREMENT_AFTER_APP <- ifelse(pgm.car.restoration.candidates$num_modulo_new <= 4 &
                                                                         pgm.car.restoration.candidates$FOREST_COVER_2007_PP < .5 &
                                                                         pgm.car.restoration.candidates$APP_FOREST_COVER_INCREMENT_PP < pgm.car.restoration.candidates$FOREST_COVER_2007_PP, 1, pgm.car.restoration.candidates@data$NEED_INCREMENT_AFTER_APP)

pgm.car.restoration.candidates@data$NEED_INCREMENT_AFTER_APP <- ifelse(pgm.car.restoration.candidates$num_modulo_new <= 4 &
                                                                         pgm.car.restoration.candidates$FOREST_COVER_2007_PP < .5 &
                                                                         pgm.car.restoration.candidates$APP_FOREST_COVER_INCREMENT_PP >= pgm.car.restoration.candidates$FOREST_COVER_2007_PP, 0, pgm.car.restoration.candidates@data$NEED_INCREMENT_AFTER_APP)

#if the area is an agrarian reform settlements, 
#it was treated as small property
pgm.car.restoration.candidates@data$NEED_INCREMENT_AFTER_APP <- ifelse(pgm.car.restoration.candidates$TIPO_IMOVE=="AST" &
                                                                         pgm.car.restoration.candidates$FOREST_COVER_2007_PP >= .5 &
                                                                         pgm.car.restoration.candidates$APP_FOREST_COVER_INCREMENT_PP < .5, 1, pgm.car.restoration.candidates@data$NEED_INCREMENT_AFTER_APP)

#View(pgm.car.restoration.candidates@data)
#table(pgm.car.restoration.candidates$NEED_INCREMENT_AFTER_APP)




#creating a copy
pgm.car.restoration.candidates.app <- pgm.car.restoration.candidates
#pgm.car.restoration.candidates <- pgm.car.restoration.candidates.app


#writeOGR(pgm.car.restoration.candidates, dsn="rasters/PGM", layer = "pgm.car.ed4", driver = "ESRI Shapefile", morphToESRI=T)
#pgm.car.restoration.candidates <- readOGR(dsn = "rasters/PGM", layer = "pgm.car.ed4")
#names(pgm.car.restoration.candidates@data) <- c("COD_IMOVEL","NUM_AREA","COD_ESTADO","NOM_MUNICI","NUM_MODULO","TIPO_IMOVE","SITUACAO","CONDICAO_I",
#                                                "num_area_ha","num_modulo_new","num_area_flag","FOREST_COVER_2007","FOREST_COVER_2007_PP","FOREST_COVER_2010","FOREST_COVER_2010_PP",
#                                                "NEED_INCREMENT","APP_FOREST_COVER_INCREMENT","APP_FOREST_COVER_INCREMENT_PP","NEED_INCREMENT_AFTER_APP")


#forest cover increment -- adding the legal reserve candidate areas to forest cover
forest.class.202X <- forest.class.2010
pgm.car.restoration.candidates@data$ARL_FOREST_COVER_INCREMENT <- NA
j=nrow(pgm.car.restoration.candidates@data[pgm.car.restoration.candidates@data$NEED_INCREMENT_AFTER_APP==1,])
n=0
#i="PA-1505502-10BE5221CA6E414BAC1F28DC4739CBA4"
for (i in pgm.car.restoration.candidates$COD_IMOVEL) {
  
  if(pgm.car.restoration.candidates@data[pgm.car.restoration.candidates@data$COD_IMOVEL==i,"NEED_INCREMENT_AFTER_APP"] != 1) next
  
  rural.property <- pgm.car.restoration.candidates[pgm.car.restoration.candidates$COD_IMOVEL==i,]
  
  forest.cover.2010 <- crop(forest.class.2010, extent(rural.property))
  forest.cover.2010 <- mask(forest.cover.2010, rural.property)
  
  restored.cover.app <- crop(candidate.areas.water, extent(rural.property))
  restored.cover.app <- mask(restored.cover.app, rural.property)
  
  forest.cover.increment.app <- sum(forest.cover.2010, restored.cover.app, na.rm = T)
  forest.cover.increment.app[forest.cover.increment.app[]>1]<-1
  forest.cover.increment.app[forest.cover.increment.app[]!=1]<-NA
  
  restored.cover.arl <- crop(candidate.areas.rl, extent(rural.property))
  restored.cover.arl <- mask(restored.cover.arl, rural.property)
  
  forest.cover.increment.arl <- sum(forest.cover.increment.app, restored.cover.arl, na.rm = T)
  forest.cover.increment.arl[forest.cover.increment.arl[]>1]<-1
  forest.cover.increment.arl[forest.cover.increment.arl[]!=1]<-NA
  
  if(all(is.na(forest.cover.increment.arl[])))
  {
    
    pgm.car.restoration.candidates@data[pgm.car.restoration.candidates@data$COD_IMOVEL==i,"ARL_FOREST_COVER_INCREMENT"] <- pgm.car.restoration.candidates@data[pgm.car.restoration.candidates@data$COD_IMOVEL==i,"APP_FOREST_COVER_INCREMENT"]
    j=j-1
    n=n+1
    cat("\n>there are no forests in property", i, "<\n")
    cat("\n>", j, "out of", nrow(pgm.car.restoration.candidates@data[pgm.car.restoration.candidates@data$NEED_INCREMENT_AFTER_APP==1,]), "properties left<\n")
    
  } 
  else
  {
    
    rural.property.5kmbuffer <- gBuffer(rural.property, width = 0.04500005)
    
    landscape.forest.cover.2010 <- crop(forest.class.2010, extent(rural.property.5kmbuffer))
    #landscape.forest.cover.2010 <- mask(landscape.forest.cover.2010, rural.property)
    
    landscape.restored.cover.app <- crop(candidate.areas.water, extent(rural.property.5kmbuffer))
    #landscape.restored.cover.app <- mask(landscape.restored.cover.app, rural.property)
    
    landscape.forest.cover.increment.app <- sum(landscape.forest.cover.2010, landscape.restored.cover.app, na.rm = T)
    landscape.forest.cover.increment.app[landscape.forest.cover.increment.app[]>1]<-1
    landscape.forest.cover.increment.app[landscape.forest.cover.increment.app[]!=1]<-NA
    
    deforest.dist <- raster::distance(landscape.forest.cover.increment.app, doEdge=T)
    
    forest.increment.mask <- crop(deforest.dist, extent(rural.property))
    forest.increment.mask <- mask(forest.increment.mask, rural.property)
    #forest.increment.mask <- mask(forest.increment.mask, forest.cover.increment.app, inverse=T)
    
    
    if(rural.property@data$num_modulo_new > 4 & rural.property@data$APP_FOREST_COVER_INCREMENT_PP < .5)
    {
      
      cat("\n>", i, "is a big property with less than 50% forest cover after checking for app<\n")
      
      target_forest_cells <- 0.5 * ncell(forest.increment.mask[!is.na(forest.increment.mask[])])
      
      #sorted_cells <- order(forest.increment.mask[forest.increment.mask[]!=0])
      
      current_forest_cells <- sum(values(forest.cover.increment.app) == 1, na.rm = T)
      current_forest_raster <- forest.cover.increment.app
      
      distance_threshold <- 30
      
      while (current_forest_cells < target_forest_cells) {
        # Identify the cells within the distance threshold
        cells_to_add <- which(!is.na(values(forest.increment.mask)) & values(forest.increment.mask) != 0 & values(forest.increment.mask) <= distance_threshold)
        
        # Increment the forest cover for the selected cells
        current_forest_raster[cells_to_add] <- 1
        
        # Update the count of forested cells
        current_forest_cells <- sum(values(current_forest_raster) == 1, na.rm = T)
        
        # Update the distance threshold for the next iteration
        distance_threshold <- distance_threshold + 30
      }
      
    }
    
    
    if(rural.property@data$num_modulo_new <= 4 & rural.property@data$FOREST_COVER_2007_PP >= .5 & rural.property@data$APP_FOREST_COVER_INCREMENT_PP < .5)
    {
      
      cat("\n>", i, "is a small property with more than 50% forest cover in 2008 but less than 50% forest cover today after checking for app<\n")
      
      target_forest_cells <- 0.5 * ncell(forest.increment.mask[!is.na(forest.increment.mask[])])
      
      #sorted_cells <- order(forest.increment.mask[forest.increment.mask[]!=0])
      
      current_forest_cells <- sum(values(forest.cover.increment.app) == 1, na.rm = T)
      current_forest_raster <- forest.cover.increment.app
      
      distance_threshold <- 30
      
      while (current_forest_cells < target_forest_cells) {
        # Identify the cells within the distance threshold
        cells_to_add <- which(!is.na(values(forest.increment.mask)) & values(forest.increment.mask) != 0 & values(forest.increment.mask) <= distance_threshold)
        
        # Increment the forest cover for the selected cells
        current_forest_raster[cells_to_add] <- 1
        
        # Update the count of forested cells
        current_forest_cells <- sum(values(current_forest_raster) == 1, na.rm = T)
        
        # Update the distance threshold for the next iteration
        distance_threshold <- distance_threshold + 30
      }
      
    }
    
    
    if(rural.property@data$num_modulo_new <= 4 & rural.property@data$FOREST_COVER_2007_PP < .5 & rural.property@data$APP_FOREST_COVER_INCREMENT_PP < rural.property@data$FOREST_COVER_2007_PP)
    {
      
      cat("\n>", i, "is a small property with less than 50% forest cover in 2008 and less forest cover than 2008 after checking for app<\n")
      
      target_forest_cells <- rural.property@data$FOREST_COVER_2007_PP * ncell(forest.increment.mask[!is.na(forest.increment.mask[])])
      
      #sorted_cells <- order(forest.increment.mask[forest.increment.mask[]!=0])
      
      current_forest_cells <- sum(values(forest.cover.increment.app) == 1, na.rm = T)
      current_forest_raster <- forest.cover.increment.app
      
      distance_threshold <- 30
      
      while (current_forest_cells < target_forest_cells) {
        # Identify the cells within the distance threshold
        cells_to_add <- which(!is.na(values(forest.increment.mask)) & values(forest.increment.mask) != 0 & values(forest.increment.mask) <= distance_threshold)
        
        # Increment the forest cover for the selected cells
        current_forest_raster[cells_to_add] <- 1
        
        # Update the count of forested cells
        current_forest_cells <- sum(values(current_forest_raster) == 1, na.rm = T)
        
        # Update the distance threshold for the next iteration
        distance_threshold <- distance_threshold + 30
      }
      
    }
    
    
    
    
    
    
    
    pgm.car.restoration.candidates@data[pgm.car.restoration.candidates@data$COD_IMOVEL==i,"ARL_FOREST_COVER_INCREMENT"] <- tapply(raster::area(current_forest_raster), current_forest_raster[], sum, na.rm=T)*100
    j=j-1
    cat("\n>", j, "out of", nrow(pgm.car.restoration.candidates@data[pgm.car.restoration.candidates@data$NEED_INCREMENT_AFTER_APP==1,]), "properties left<\n")
    
    
  }
  
  
  forest.class.202X[!is.na(current_forest_raster) & current_forest_raster == 1] <- 1
  
  
}



#checking
#length(forest.class.2010[forest.class.2010[]==1])
#tapply(raster::area(forest.class.2010), forest.class.2010[], sum, na.rm=T)*100
#length(forest.class.202X[forest.class.202X[]==1])
#tapply(raster::area(forest.class.202X), forest.class.202X[], sum, na.rm=T)*100
# Plot the original forest and the updated forest
#par(mfrow = c(1, 2))
#plot(forest.class.2010, main = "Original Forested Area", col = c("white", "green"), legend=F)
#plot(forest.class.202X, main = "Updated Forested Area", col = c("white", "green"), legend=F)

writeRaster(forest.class.202X, "results/pgm_forest_cover_after_restoration.tif", fomrat = "GTiff")


#View(pgm.car.restoration.candidates@data)
#anyNA(pgm.car.restoration.candidates@data$ARL_FOREST_COVER_INCREMENT)
#length(which(is.na(pgm.car.restoration.candidates@data$ARL_FOREST_COVER_INCREMENT)))

pgm.car.restoration.candidates@data$ARL_FOREST_COVER_INCREMENT <- ifelse(is.na(pgm.car.restoration.candidates@data$ARL_FOREST_COVER_INCREMENT),
                                                                         pgm.car.restoration.candidates@data$APP_FOREST_COVER_INCREMENT, pgm.car.restoration.candidates@data$ARL_FOREST_COVER_INCREMENT)

#calculating the percentage of forest cover in relation to property size
pgm.car.restoration.candidates@data$ARL_FOREST_COVER_INCREMENT_PP <- pgm.car.restoration.candidates@data$ARL_FOREST_COVER_INCREMENT/pgm.car.restoration.candidates@data$num_area_ha




#select properties based on threshold
#if the area is big/medium properties with less than 50%
pgm.car.restoration.candidates@data$NEED_INCREMENT_AFTER_ARL <- ifelse(pgm.car.restoration.candidates@data$num_modulo_new > 4 &
                                                                         pgm.car.restoration.candidates@data$ARL_FOREST_COVER_INCREMENT_PP < .5, 1, 0)

#if the area is small properties with less the 2007 forest cover
pgm.car.restoration.candidates@data$NEED_INCREMENT_AFTER_ARL <- ifelse(pgm.car.restoration.candidates$num_modulo_new <= 4 &
                                                                         pgm.car.restoration.candidates$FOREST_COVER_2007_PP >= .5 &
                                                                         pgm.car.restoration.candidates$ARL_FOREST_COVER_INCREMENT_PP < .5, 1, pgm.car.restoration.candidates@data$NEED_INCREMENT_AFTER_ARL)

pgm.car.restoration.candidates@data$NEED_INCREMENT_AFTER_ARL <- ifelse(pgm.car.restoration.candidates$num_modulo_new <= 4 &
                                                                         pgm.car.restoration.candidates$FOREST_COVER_2007_PP < .5 &
                                                                         pgm.car.restoration.candidates$ARL_FOREST_COVER_INCREMENT_PP < pgm.car.restoration.candidates$FOREST_COVER_2007_PP, 1, pgm.car.restoration.candidates@data$NEED_INCREMENT_AFTER_ARL)

pgm.car.restoration.candidates@data$NEED_INCREMENT_AFTER_ARL <- ifelse(pgm.car.restoration.candidates$num_modulo_new <= 4 &
                                                                         pgm.car.restoration.candidates$FOREST_COVER_2007_PP < .5 &
                                                                         pgm.car.restoration.candidates$ARL_FOREST_COVER_INCREMENT_PP >= pgm.car.restoration.candidates$FOREST_COVER_2007_PP, 0, pgm.car.restoration.candidates@data$NEED_INCREMENT_AFTER_ARL)

#if the area is an agrarian reform settlements, 
#it was treated as small property
pgm.car.restoration.candidates@data$NEED_INCREMENT_AFTER_ARL <- ifelse(pgm.car.restoration.candidates$TIPO_IMOVE=="AST" &
                                                                         pgm.car.restoration.candidates$FOREST_COVER_2007_PP >= .5 &
                                                                         pgm.car.restoration.candidates$ARL_FOREST_COVER_INCREMENT_PP < .5, 1, pgm.car.restoration.candidates@data$NEED_INCREMENT_AFTER_ARL)

#View(pgm.car.restoration.candidates@data)
#table(pgm.car.restoration.candidates$NEED_INCREMENT_AFTER_ARL)







writeOGR(pgm.car.restoration.candidates, dsn="results", layer = "pgm_car_after_restoration", driver = "ESRI Shapefile", morphToESRI=T)
#pgm.car.restoration.candidates <- readOGR(dsn = "results/", layer = "pgm_car_after_restoration")
#names(pgm.car.restoration.candidates@data) <- c("COD_IMOVEL","NUM_AREA","COD_ESTADO","NOM_MUNICI","NUM_MODULO","TIPO_IMOVE","SITUACAO","CONDICAO_I",
#                                                "num_area_ha","num_modulo_new","num_area_flag","FOREST_COVER_2007","FOREST_COVER_2007_PP","FOREST_COVER_2010","FOREST_COVER_2010_PP",
#                                                "NEED_INCREMENT","APP_FOREST_COVER_INCREMENT","APP_FOREST_COVER_INCREMENT_PP","NEED_INCREMENT_AFTER_APP",
#                                                "ARL_FOREST_COVER_INCREMENT","ARL_FOREST_COVER_INCREMENT_PP","NEED_INCREMENT_AFTER_ARL")


rm(list=ls()[ls() %in% c("[...]")])
gc()


#

#Fig S1 and S2
#use this part of the script only after load objects in "layer_buide_[...].R" ==
## standard projection
std.proj <- "+proj=longlat +datum=WGS84 +units=m +no_defs"

## shapefile paragominas
pgm.shp <- readOGR(dsn = "shapes", layer = "Paragominas_Mask_R3")
proj4string(pgm.shp) <- CRS(std.proj)
pgm.shp <- spTransform(pgm.shp, crs(std.proj))

### shapefile santarem
#stm.shp <- readOGR(dsn = "shapes", layer = "Santarem")
#stm.shp <- spTransform(stm.shp, crs(std.proj))

## supplementary material figures 
###undegraded primary forest == 1
#UPF2010_real.sk <- UPF2010_real
UPF2010_real.sk <- raster("rasters/PGM/input/LULC/UPF2010_real.tif")
#UPF2010_real.sk[UPF2010_real.sk==1]<-1
UPF2010_real.sk <- mask(UPF2010_real.sk, pgm.shp)
UPF2010_real.sk[UPF2010_real.sk[]==0] <- 333
#plot(UPF2010_real.sk, main="undegraded primary forest", legend=F)

###degraded primary forest == 10
#uDPF2010_real.sk <- uDPF2010_real
uDPF2010_real.sk <- raster("rasters/PGM/input/LULC/uDPF2010_real.tif")
uDPF2010_real.sk[uDPF2010_real.sk==1]<-10
uDPF2010_real.sk <- mask(uDPF2010_real.sk, pgm.shp)
uDPF2010_real.sk[uDPF2010_real.sk[]==0] <- 333
#plot(uDPF2010_real.sk, main="degraded primary forest", legend=F)

###repeated degraded primary forest == 25
#RDPF2010_real.sk <- RDPF2010_real
RDPF2010_real.sk <- raster("rasters/PGM/input/LULC/RDPF2010_real.tif")
RDPF2010_real.sk[RDPF2010_real.sk==1]<-25
RDPF2010_real.sk <- mask(RDPF2010_real.sk, pgm.shp)
RDPF2010_real.sk[RDPF2010_real.sk[]==0] <- 333
#plot(RDPF2010_real.sk, main="degraded primary forest", legend=F)

###secondary forest == 100
#uSF2010_real.sk <- uSF2010_real
uSF2010_real.sk <- raster("rasters/PGM/input/LULC/uSF2010_real.tif")
uSF2010_real.sk[uSF2010_real.sk==1]<-100
uSF2010_real.sk <- mask(uSF2010_real.sk, pgm.shp)
uSF2010_real.sk[uSF2010_real.sk[]==0] <- 333
#plot(uSF2010_real.sk, main="secondary forest", legend=F)

###degraded secondary forest == 125
#DSF2010_real.sk <- DSF2010_real
DSF2010_real.sk <- raster("rasters/PGM/input/LULC/DSF2010_real.tif")
DSF2010_real.sk[DSF2010_real.sk==1]<-125
DSF2010_real.sk <- mask(DSF2010_real.sk, pgm.shp)
DSF2010_real.sk[DSF2010_real.sk[]==0] <- 333
#plot(DSF2010_real.sk, main="secondary forest", legend=F)

LULC2010_real <- sum(UPF2010_real.sk, uSF2010_real.sk, na.rm = T)
LULC2010_real <- sum(LULC2010_real, DSF2010_real.sk, na.rm = T)
LULC2010_real <- sum(LULC2010_real, uDPF2010_real.sk, na.rm = T)
LULC2010_real <- sum(LULC2010_real, RDPF2010_real.sk, na.rm = T)
LULC2010_real <- mask(LULC2010_real, pgm.shp)
sort(unique(LULC2010_real[]))


LULC2010_real[LULC2010_real==1333]<-1
LULC2010_real[LULC2010_real==1342]<-10
LULC2010_real[LULC2010_real==1357]<-25
LULC2010_real[LULC2010_real==1432]<-100
LULC2010_real[LULC2010_real==1457]<-125
LULC2010_real[LULC2010_real==1665 ]<-0
#plot(LULC2010_real, main="forest cover 2010", legend=F)

writeRaster(LULC2010_real, "rasters/PGM/input/LULC/_lulc_2010_real.tif", format="GTiff", overwrite=T)



### Create a data frame with the transition data
data_df2 <- data.frame(
  Period1 = factor(LULC2010_real[], levels = c(1, 10, 25, 125, 100, 0), labels = c("UPF", "DPF", "RDPF", "DSF", "SF", "D")),
  Period2 = factor(LULC2020_real[], levels = c(1, 10, 25, 125, 100, 0), labels = c("UPF", "DPF", "RDPF", "DSF", "SF", "D"))
  #Period2 = factor(LULC2020_avoiddegrad[], levels = c(1, 10, 25, 125, 100, 0), labels = c("UPF", "DPF", "RDPF", "DSF", "SF", "D"))
  #Period2 = factor(LULC2020_avoiddegrad2[], levels = c(1, 10, 25, 125, 100, 0), labels = c("UPF", "DPF", "RDPF", "DSF", "SF", "D"))
  #Period2 = factor(LULC2020_avoiddeforest[], levels = c(1, 10, 25, 125, 100, 0), labels = c("UPF", "DPF", "RDPF", "DSF", "SF", "D"))
  #Period2 = factor(LULC2020_avoiddeforest2[], levels = c(1, 10, 25, 125, 100, 0), labels = c("UPF", "DPF", "RDPF", "DSF", "SF", "D"))
  #Period2 = factor(LULC2020_restor_wo_avoid[], levels = c(1, 10, 25, 125, 100, 0), labels = c("UPF", "DPF", "RDPF", "DSF", "SF", "D"))
  #Period2 = factor(LULC2020_avoidboth[], levels = c(1, 10, 25, 125, 100, 0), labels = c("UPF", "DPF", "RDPF", "DSF", "SF", "D"))
  #Period2 = factor(LULC2020_avoidboth2[], levels = c(1, 10, 25, 125, 100, 0), labels = c("UPF", "DPF", "RDPF", "DSF", "SF", "D"))
  #Period2 = factor(LULC2020_restor_n_avoiddeforest[], levels = c(1, 10, 25, 125, 100, 0), labels = c("UPF", "DPF", "RDPF", "DSF", "SF", "D"))
  #Period2 = factor(LULC2020_restor_n_avoiddeforest2[], levels = c(1, 10, 25, 125, 100, 0), labels = c("UPF", "DPF", "RDPF", "DSF", "SF", "D"))
  #Period2 = factor(LULC2020_restor_n_avoidboth[], levels = c(1, 10, 25, 125, 100, 0), labels = c("UPF", "DPF", "RDPF", "DSF", "SF", "D"))
  #Period2 = factor(LULC2020_restor_n_avoidboth2[], levels = c(1, 10, 25, 125, 100, 0), labels = c("UPF", "DPF", "RDPF", "DSF", "SF", "D"))
)


###change the axis2 to check lulc transition in each scenario
set.seed(1237)
data_df2 %>% drop_na() %>% sample_n(size = 100000, replace = T) %>%
  ggplot(aes(axis1 = Period1, axis2 = Period2)) +
  geom_flow(aes(fill = Period2), width = .15, curve_type = "quintic") +
  geom_stratum(width = .15) +
  scale_x_discrete(limits = c("Period1", "Period2"), 
                   breaks=c("Period1", "Period2"), 
                   labels=addline_format(c("2010 Real", "2020 restoration and avoid both (PF_only)")), # (PF_only)
                   expand = c(.05, .05)) +
  scale_fill_manual(values = c("#294B29", "#50623A", "#76453B", 
                               "#B19470", "#789461", "#F97B22")) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal()+
  theme(axis.text.y= element_blank(), legend.position = "none")

#


# Fig S5
# importing species list and methods and evaluation metrics data frame
forestdep.spplist <- read.csv("data/species_summary.csv")


#checking
#individual.maps <- c("Pouteriamacrophylla", "Courataristellata", "Rinoreaguianensis", 
#                     "Handroanthusserratifolius", "Anibamegaphylla", "Manilkaraparaensis",
#                     "Hymenolobiumexcelsum", "Cedrelaodorata", "Automolusparaensis",
#                     "Buccotamatia", "Hylopezusmacularius", "Phlegopsisnigromaculata",
#                     "Geotrygonmontana", "Thamnomanescaesius", "Aburriacujubi",
#                     "Corythopistorquatus")

pgm.individual.maps.list <- list.files("models.output/biodiversity.maps/PGM/2010_real/", pattern = ".tif", full.names = T, recursive = T)
#pgm.individual.maps <- grep(paste(individual.maps, collapse = "|"), pgm.individual.maps.list, value = T)
#pgm.individual.maps <- stack(pgm.individual.maps)
pgm.individual.map1 <- raster(grep("Anibamegaphylla", pgm.individual.maps.list, value = T))
pgm.individual.map2 <- raster(grep("Phlegopsisnigromaculata", pgm.individual.maps.list, value = T))
pgm.ext <- extent(-47.35, -47.28, -3.35, -3.27)


stm.individual.maps.list <- list.files("models.output/biodiversity.maps/STM/2010_real/", pattern = ".tif", full.names = T, recursive = T)
#stm.individual.maps <- grep(paste(individual.maps, collapse = "|"), stm.individual.maps.list, value = T)
#stm.individual.maps <- stack(stm.individual.maps)
stm.individual.map1 <- raster(grep("Anibamegaphylla", stm.individual.maps.list, value = T))
stm.individual.map2 <- raster(grep("Phlegopsisnigromaculata", stm.individual.maps.list, value = T))
stm.ext <- extent(-54.83, -54.78, -2.98, -2.92)


## res = 1673 x 881
par(mfrow = c(4,2))
par(mar = c(2, 2.5, 4, .2))
plot(pgm.individual.map1, 
     col = c("#FBFCF8", "#FCF596", "#FBD28B", "#FF9C73", "#FF4545"), 
     breaks= seq(0, 1, by = .2), legend =F)

plot(pgm.ext, add=T)

par(mar = c(2, .2, 4, .2))
plot(pgm.individual.map1, ext = pgm.ext,
     col = c("#FBFCF8", "#FCF596", "#FBD28B", "#FF9C73", "#FF4545"), 
     breaks= seq(0, 1, by = .2))

par(mar = c(4, 2.5, 2, .2))
plot(stm.individual.map1, 
     col = c("#FBFCF8", "#FCF596", "#FBD28B", "#FF9C73", "#FF4545"), 
     breaks= seq(0, 1, by = .2), legend =F)

plot(stm.ext, add=T)

par(mar = c(4, .2, 2, .2))
plot(stm.individual.map1, ext = stm.ext,
     col = c("#FBFCF8", "#FCF596", "#FBD28B", "#FF9C73", "#FF4545"), 
     breaks= seq(0, 1, by = .2))


par(mar = c(2, 2.5, 4, .2))
plot(pgm.individual.map2, 
     col = c("#FBFCF8", "#FCF596", "#FBD28B", "#FF9C73", "#FF4545"), 
     breaks= seq(0, 1, by = .2), legend =F)

plot(pgm.ext, add=T)

par(mar = c(2, .2, 4, .2))
plot(pgm.individual.map2, ext = pgm.ext,
     col = c("#FBFCF8", "#FCF596", "#FBD28B", "#FF9C73", "#FF4545"), 
     breaks= seq(0, 1, by = .2))

par(mar = c(4, 2.5, 2, .2))
plot(stm.individual.map2, 
     col = c("#FBFCF8", "#FCF596", "#FBD28B", "#FF9C73", "#FF4545"), 
     breaks= seq(0, 1, by = .2), legend =F)

plot(stm.ext, add=T)

par(mar = c(4, .2, 2, .2))
plot(stm.individual.map2, ext = stm.ext,
     col = c("#FBFCF8", "#FCF596", "#FBD28B", "#FF9C73", "#FF4545"), 
     breaks= seq(0, 1, by = .2))

mtext("Aniba megaphylla", side = 3, line = -3, outer = T, font = 4)
mtext("Phlegopsis nigromaculata", side = 3, line = -38, outer = T, font = 4)


# Fig S6
layer.names <- c("Real 2010", "Real 2020", "Avoid deforestation", "Avoid degradation", 
                 "Restoration without avoid", "Avoid both", "Restoration and avoid deforestation",
                 "Restoration and avoid both", "Avoid deforestation PF only", "Avoid degradation PF only", 
                 "Avoid both PF only", "Restoration and avoid deforestation PF only",
                 "Restoration and avoid both PF only")

biodiversity.benefit.list <- list.files("models.output/biodiversity.benefits2/", pattern = ".tif", full.names = T, recursive = T)

pgm.biodiversity.benefit.list <- grep("PGM", biodiversity.benefit.list, value = T)

pgm.biodiversity.benefit.total <- stack(pgm.biodiversity.benefit.list)
pgm.biodiversity.benefit.total <- pgm.biodiversity.benefit.total[[c(1,8,4,6,13,2,11,9,5,7,3,12,10)]]
pgm.biodiversity.benefit.total <- mask(pgm.biodiversity.benefit.total, pgm.shp)

rm(pgm.biodiversity.benefit.list)

for (i in c(1:13)) {
  
  pgm.biodiversity.benefit.total[[i]][pgm.lulc[[i]][] == 0] <- 0
  names(pgm.biodiversity.benefit.total[[i]]) <- layer.names[i]
  
  cat("\n>working on layer", i, "now<\n")
}

gc()



stm.biodiversity.benefit.list <- grep("STM", biodiversity.benefit.list, value = T)

stm.biodiversity.benefit.total <- stack(stm.biodiversity.benefit.list)
stm.biodiversity.benefit.total <- stm.biodiversity.benefit.total[[c(1,8,4,6,13,2,11,9,5,7,3,12,10)]]
stm.biodiversity.benefit.total <- mask(stm.biodiversity.benefit.total, stm.shp)

rm(stm.biodiversity.benefit.list)

for (i in c(1:13)) {
  
  stm.biodiversity.benefit.total[[i]][stm.lulc[[i]][] == 0] <- 0
  names(stm.biodiversity.benefit.total[[i]]) <- layer.names[i]
  
  cat("\n>working on layer", i, "now<\n")
}

gc()



plot(stm.biodiversity.benefit.total[[2:13]], nr=3, 
     col = colorRampPalette(c("#FBFCF8", "#0AD1C8", "#14919B", "#213A57"))(length(seq(0, 410, by = 10))), 
     breaks= seq(0, 420, by = 10)) ## res = 1673 x 881



# Fig S8
carbon.benefit.list <- list.files("models.output/carbon.benefits/", pattern = ".tif", full.names = T, recursive = T)

pgm.carbon.benefit.list <- grep("PGM", carbon.benefit.list, value = T)

pgm.carbon.benefit.total <- stack(pgm.carbon.benefit.list)
pgm.carbon.benefit.total <- pgm.carbon.benefit.total[[c(1,8,4,6,13,2,11,9,5,7,3,12,10)]]
pgm.carbon.benefit.total <- mask(pgm.carbon.benefit.total, pgm.shp)

rm(pgm.carbon.benefit.list)

for (i in c(1:13)) {
  
  pgm.carbon.benefit.total[[i]][pgm.lulc[[i]][] == 0] <- 0
  names(pgm.carbon.benefit.total[[i]]) <- layer.names[i]
  
  cat("\n>working on layer", i, "now<\n")
}

gc()



stm.carbon.benefit.list <- grep("STM", carbon.benefit.list, value = T)

stm.carbon.benefit.total <- stack(stm.carbon.benefit.list)
stm.carbon.benefit.total <- stm.carbon.benefit.total[[c(1,8,4,6,13,2,11,9,5,7,3,12,10)]]
stm.carbon.benefit.total <- mask(stm.carbon.benefit.total, stm.shp)

rm(stm.carbon.benefit.list)

for (i in c(1:13)) {
  
  stm.carbon.benefit.total[[i]][stm.lulc[[i]][] == 0] <- 0
  names(stm.carbon.benefit.total[[i]]) <- layer.names[i]
  
  cat("\n>working on layer", i, "now<\n")
}

gc()



plot(stm.carbon.benefit.total[[2:13]], nr=3, 
     col = colorRampPalette(c("#FBFCF8", "#8fce00", "#374f00"))(length(seq(0, 200, by = 25))), 
     breaks= seq(0, 225, by = 25)) ## res = 1673 x 881



# Fig S9
# where is Paragominas and Santarem compared to the Brazilian Amazon? ==========
blm.states <- c("Rondônia", "Acre", "Amazonas", "Roraima", "Pará", "Amapá", "Tocantins", "Maranhão", "Mato Grosso")

#import data on deforestation from mapbiomas
#all municipalities in BLA
total.blm.deforestation.2010 <- readxl::read_xlsx("data/raw/mapbiomas_brasil_col9_state_municipality.xlsx", sheet = 2) %>% 
  filter(state %in% blm.states) %>%
  mutate(class_level_1 = ifelse(class_level_1 %in% c("3. Farming", "4. Non vegetated area"), "3. Deforest", class_level_1)) %>% 
  group_by(state, municipality, class_level_1) %>% 
  summarise(across('1985':'2020', sum)) %>% 
  dplyr::select(state, municipality, class_level_1, `2010`, `2020`) %>% 
  mutate(freq_2010 = `2010` / sum(`2010`),
         freq_2020 = `2020` / sum(`2020`))


#priority municipalities in deforestation arch
priority.mun.list <- c("Feijó","Manoel Urbano","Rio Branco","Sena Madureira","Tarauacá","Apuí","Boca do Acre",
                       "Canutama","Humaitá","Itapiranga","Lábrea","Manicoré","Maués","Novo Apurinã","Apiacás",
                       "Aripuanã","Bom Jesus do Araguaia","Cláudia","Colniza","Comodoro","Cotriguaçu","Feliz Natal",
                       "Gaúcha do Norte","Juara","Juína","Marcelândia","Nova Bandeirantes","Nova Maringá",
                       "Nova Ubiratã","Paranaíta","Paranatinga","Peixoto de Azevedo","Querência","Rondonlandia",
                       "São José do Xingú","União do Sul","Altamira","Anapu","Cumaru do Norte","Dom Eliseu","Itaituba",
                       "Itupiranga","Jacareaganga","Marabá","Medicilândia","Moju","Mojuí dos Campos","Novo Progresso",
                       "Novo Repartimento","Pacajá","Paragominas","Placas","Portel","Prainha","Rondon do Pará",
                       "Rurópolis","Santana do Araguaia","São Félix do Xingu","Senador José Porfírio","Trairão",
                       "Ulianópolis","Uruará","Buritis","Candeias Do Jamari","Cujubim","Machadinho D'Oeste",
                       "Nova Mamoré","Porto Velho","Mucajaí","Rorainópolis")

priority.mun.deforestation.2010 <- total.blm.deforestation.2010 %>% 
  filter(municipality %in% priority.mun.list)


#PGM and STM
pgm.deforestation.2010 <- total.blm.deforestation.2010 %>% 
  filter(municipality == "Paragominas" & class_level_1 == "3. Deforest") %>% 
  ungroup() %>% dplyr::select(freq_2010) %>% pull()

pgm.deforestation.2020 <- total.blm.deforestation.2010 %>% 
  filter(municipality == "Paragominas" & class_level_1 == "3. Deforest") %>% 
  ungroup() %>% dplyr::select(freq_2020) %>% pull()

stm.deforestation.2010 <- total.blm.deforestation.2010 %>% 
  filter(municipality == "Santarém" & class_level_1 == "3. Deforest") %>% 
  ungroup() %>% dplyr::select(freq_2010) %>% pull()

stm.deforestation.2020 <- total.blm.deforestation.2010 %>% 
  filter(municipality == "Santarém" & class_level_1 == "3. Deforest") %>% 
  ungroup() %>% dplyr::select(freq_2020) %>% pull()



# why/where benefit is higher? ====================================
## comparing benefit between avoid degradation and avoid deforestation

env <- c("TSDls", "edgedist", "UPFls")

pgm.env.2010 <- list.files("rasters/PGM/2010_real", pattern = ".tif", full.names = T, recursive = T)
pgm.env.2010 <- grep(paste(env, collapse = "|"), pgm.env.2010, value = T)
pgm.env.2010 <- stack(pgm.env.2010)
pgm.env.2010 <- mask(pgm.env.2010, pgm.shp)
pgm.env.2010.df <- as.data.frame(pgm.env.2010, xy = TRUE) %>%  
  mutate(Region = "PGM", Cell = row_number()) %>% drop_na()


stm.env.2010 <- list.files("rasters/STM/2010_real", pattern = ".tif", full.names = T, recursive = T)
stm.env.2010 <- grep(paste(env, collapse = "|"), stm.env.2010, value = T)
stm.env.2010 <- stack(stm.env.2010)
stm.env.2010 <- mask(stm.env.2010, stm.shp)
stm.env.2010.df <- as.data.frame(stm.env.2010, xy = TRUE) %>%  
  mutate(Region = "STM", Cell = row_number()) %>% drop_na()


env.2010 <- rbind(pgm.env.2010.df, stm.env.2010.df)


cell.deforest <- costs.principals %>% 
  filter(Scenario == "Avoid deforestation", area_change == "Direct") %>% 
  dplyr::select(Region, Cell, Cat, BBenefit, CBenefit) %>% 
  left_join(env.2010[,3:7])

cell.degrad <- costs.principals %>% 
  filter(Scenario == "Avoid degradation", area_change == "Direct") %>% 
  dplyr::select(Region, Cell, Cat, BBenefit, CBenefit) %>% 
  left_join(env.2010[,3:7])

cell.restor <- costs.principals %>% 
  filter(Scenario == "Restoration without avoid", area_change == "Direct") %>% 
  dplyr::select(Region, Cell, Cat, BBenefit, CBenefit) %>% 
  left_join(env.2010[,3:7])

cell.full <- costs.principals %>% 
  filter(Scenario == "Restoration and avoid both", area_change == "Direct") %>% 
  dplyr::select(Region, Cell, Cat, BBenefit, CBenefit) %>% 
  left_join(env.2010[,3:7])


ggarrange(
  ggplot() +
    geom_histogram(data=total.blm.deforestation.2010 %>% filter(class_level_1 == "3. Deforest"), 
                   aes(x=freq_2010, y=after_stat(count/nrow(total.blm.deforestation.2010)), fill = "All municipalities"), bins = 70) +
    geom_histogram(data=priority.mun.deforestation.2010 %>% filter(class_level_1 == "3. Deforest"), 
                   aes(x=freq_2010, y=after_stat(count/nrow(priority.mun.deforestation.2010)), fill = "Priority municipalities"), bins = 70, alpha = 0.55) +
    scale_fill_manual(values = c("All municipalities"="#e1d3cc", "Priority municipalities"="#e99561")) +
    geom_vline(xintercept = pgm.deforestation.2010, linetype = "dashed", color = "gray20", linewidth = 1) +
    geom_vline(xintercept = stm.deforestation.2010, linetype = "dashed", color = "gray20", linewidth = 1) +
    #scale_y_continuous(limits = c(0, 0.0000006)) +
    #scale_x_continuous(limits = c(0, 2000000)) +
    labs(x="Deforestation in the Amazon \nmunicipalities by 2010", y="") +
    theme_minimal() +
    theme(text = element_text(size = 16, family = "sans"),
          plot.title = element_text(hjust = 0.5),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 14),
          legend.title = element_blank(),
          legend.position = "top"),
  
  ggarrange(
    
    ggplot() +
      geom_histogram(data=env.2010, 
                     aes(x=edgedist, y=after_stat(count / nrow(env.2010))), fill = "#e1d3cc") +
      geom_freqpoly(data=cell.deforest, 
                    aes(x=edgedist, y=after_stat(count / nrow(cell.deforest)), color = "Deforestation"), linewidth = 1) +
      geom_freqpoly(data=cell.degrad, 
                    aes(x=edgedist, y=after_stat(count / nrow(cell.degrad)), color = "Degradation"), linewidth = 1, linetype = "dashed") +
      scale_x_continuous(limits = c(0, 7000)) +
      scale_y_continuous(limits = c(0, 0.4)) +
      scale_color_manual(values = c("Deforestation"="#294B29", "Degradation"="#50623A")) +
      labs(title = "", x = "Distance to the edge", y = "Proportion") +
      facet_wrap(~Region) +
      theme_minimal() +
      theme(text = element_text(size = 16, family = "sans"),
            plot.title = element_text(hjust = 0.5),
            axis.title = element_text(face="bold"),
            axis.text.x=element_text(size = 14),
            legend.title = element_blank()),
    
    ggplot() +
      geom_histogram(data=env.2010, 
                     aes(x=TSDls, y=after_stat(count / nrow(env.2010))), fill = "#e1d3cc") +
      geom_freqpoly(data=cell.deforest, 
                    aes(x=TSDls, y=after_stat(count / nrow(cell.deforest)), color = "Deforestation"), linewidth = 1) +
      geom_freqpoly(data=cell.degrad, 
                    aes(x=TSDls, y=after_stat(count / nrow(cell.degrad)), color = "Degradation"), linewidth = 1, linetype = "dashed") +
      scale_y_continuous(limits = c(0, 0.4)) +
      scale_color_manual(values = c("Deforestation"="#294B29", "Degradation"="#50623A")) +
      labs(title = "", x = "Time since degradation", y = "") +
      facet_wrap(~Region) +
      theme_minimal() +
      theme(text = element_text(size = 16, family = "sans"),
            plot.title = element_text(hjust = 0.5),
            axis.title = element_text(face="bold"),
            axis.text.x=element_text(size = 14),
            legend.title = element_blank()),
    
    nrow = 2, common.legend = T, legend = "bottom"),
  
  nrow = 2, heights = c(0.35, 0.65), legend = "top"
)  ## res = 881 x 1673












par(mfrow=c(1,2))
#plot(pgm.costs.total[[1]], main="Fire control", col = terrain.colors(9, rev = T), breaks= c(0,10,20,30,40,50,60,70,210)) ## res = 1673 x 881
#plot(pgm.costs.total[[4]], main="Restoration", col = terrain.colors(5, rev = T), breaks= c(0,10,100,500,800))
plot(pgm.costs.total[[3]], main="Opportunity farming", col = terrain.colors(9, rev = T), breaks= c(0,100,200,300,400,500,600,1000,3500))
#plot(pgm.costs.total[[2]], main="Opportunity logging", col = terrain.colors(9, rev = T), breaks= seq(0, 200, by = 20))
#plot(stm.costs.total[[1]], main="Fire control", col = terrain.colors(10, rev = T), breaks= c(0,10,20,30,40,50,60,70,100,300))
#plot(stm.costs.total[[4]], main="Restoration", col = terrain.colors(5, rev = T), breaks= c(0,10,100,500,800))
plot(stm.costs.total[[3]], main="Opportunity farming", col = terrain.colors(9, rev = T), breaks= c(0,100,200,300,400,500,600,1000,3000))
#plot(stm.costs.total[[2]], main="Opportunity logging", col = terrain.colors(9, rev = T), breaks= seq(0, 200, by = 20))




#
#




#=================================|end
























































#==============================| previous approach





# checking for impossible transition scenarios
#e.g., secondary forests becoming undegraded/degraded primary forest
# or non-forest degraded primary forest
#par(mfrow = c(2, 3))
#PGM 2010
#undegradded primary forest == 1
UPF2010.sk <- UPF2010
#UPF2010.sk <- raster("rasters/PGM/input/UPF2010_real.tif")
#UPF2010.sk[UPF2010.sk==1]<-1
#plot(UPF2010.sk, main="undegradded primary forest", legend=F)

#degradded primary forest == 10
DPF2010.sk <- DPF2010
#DPF2010.sk <- raster("rasters/PGM/input/DPF2010_real.tif")
DPF2010.sk[DPF2010.sk==1]<-10
#plot(DPF2010.sk, main="degradded primary forest", legend=F)

#secondary forest == 100
SF2010.sk <- SF2010
#SF2010.sk <- raster("rasters/PGM/input/SF2010_real.tif")
SF2010.sk[SF2010.sk==1]<-100
#plot(SF2010.sk, main="secondary forest", legend=F)

LULC2010 <- sum(UPF2010.sk, SF2010.sk, na.rm = T)
LULC2010 <- sum(LULC2010, DPF2010.sk, na.rm = T)
#plot(LULC2010)
#sort(unique(LULC2010[]))


#PGM 2020
#undegradded primary forest == 1
UPF2020.sk <- UPF2020
#UPF2020.sk <- raster("rasters/PGM/input/UPF2020_real.tif")
#UPF2020.sk[UPF2020.sk==1]<-1
#plot(UPF2020.sk, legend=F)

#degradded primary forest == 10
DPF2020.sk <- DPF2020
#DPF2020.sk <- raster("rasters/PGM/input/DPF2020_real.tif")
DPF2020.sk[DPF2020.sk==1]<-10
#plot(DPF2020.sk, legend=F)

#secondary forest == 100
SF2020.sk <- SF2020
#SF2020.sk <- raster("rasters/PGM/input/SF2020_real.tif")
SF2020.sk[SF2020.sk==1]<-100
#plot(SF2020.sk, legend=F)

LULC2020 <- sum(UPF2020.sk, SF2020.sk, na.rm = T)
LULC2020 <- sum(LULC2020, DPF2020.sk, na.rm = T)
#plot(LULC2020)
#sort(unique(LULC2020[]))



#set.seed(123)

# Create a data frame with the transition data
data_df <- data.frame(
  Period1 = factor(LULC2010[], levels = c(1, 10, 100, 0), labels = c("UPF", "DPF", "SF", "DF")),
  Period2 = factor(LULC2020[], levels = c(1, 10, 100, 0), labels = c("UPF", "DPF", "SF", "DF"))
)

# Create the alluvial plot
library(ggalluvial)
#> Loading required package: ggplot2
library(ggfittext)
library(scales)


#par(mfrow = c(1, 1))

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

data_df %>% drop_na() %>% sample_n(size = 50000, replace = T) %>%
  ggplot(aes(axis1 = Period1, axis2 = Period2)) +
  geom_flow(aes(fill = Period1), width = .15, curve_type = "quintic") +
  geom_stratum(width = .15) +
  scale_x_discrete(limits = c("Period1", "Period2"), 
                   breaks=c("Period1", "Period2"), 
                   labels=addline_format(c("2010 Real", "2020 Real")),
                   expand = c(.05, .05)) +
  scale_fill_manual(values = c("#263A29", "#65451F", "#83764F", "#F97B22")) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal()+
  theme(axis.text.y= element_blank(), legend.position = "none")






#data(vaccinations)
#vaccinations <- transform(vaccinations, freq = freq * 1000)
#
#
#ggplot(vaccinations,
#       aes(x = survey, stratum = response, alluvium = subject,
#           y = freq,
#           fill = response, label = freq)) +
#  scale_x_discrete(expand = c(.1, .1)) +
#  geom_flow() +
#  geom_stratum(alpha = .5) +
#  geom_fit_text(stat = "stratum", size = 10, min.size = 6, formatter = comma) +
#  theme(legend.position = "bottom") +
#  ggtitle("vaccination survey responses at three points in time")




# Load your four rasters
carbon.degradation <- pgm.conservact.avoiddegrad.carbbenefitmask
cost.degradation <- pgm.conservact.avoiddegrad.costmask

carbon.deforestation <- pgm.conservact.avoiddeforest.carbbenefitmask
cost.deforestation <- pgm.conservact.avoiddeforest.costmask

carbon.restoration <- pgm.conservact.restor_wo_avoid.carbbenefitmask
cost.restoration <- pgm.conservact.restor_wo_avoid.costmask

# Create an empty dataframe to store results
result_df <- data.frame(Budget = numeric(0),
                        Total_area = numeric(0),
                        Degradation_area = numeric(0), 
                        Deforestation_area = numeric(0), 
                        Restoration_area = numeric(0))


# Initialize areas reached
reached.areas <- raster("imagens/mapbiomas-brazil-collection-70-pgm-2010-100mpx.tif")
reached.areas[!is.na(reached.areas)] <- 0

ncell.degradation <- c()
ncell.deforestation <- c()
ncell.restoration <- c()

# Loop through each budget constraint
for (r in seq(0,100,2)) {
  # Initialize total carbon stocks and costs
  total.carbon.degradation <- 0
  total.cost.degradation <- 0
  
  total.carbon.deforestation <- 0
  total.cost.deforestation <- 0
  
  total.carbon.restoration <- 0
  total.cost.restoration <- 0

    
  # Iterate while budget is not exhausted
  while (round(sum(total.cost.deforestation, total.cost.restoration, total.cost.degradation), 0) <= r*1000000) {
    #excluding areas after conservation action
    carbon.degradation <- mask(carbon.degradation, reached.areas)
    cost.degradation <- mask(cost.degradation, reached.areas)
    
    carbon.deforestation <- mask(carbon.deforestation, reached.areas)
    cost.deforestation <- mask(cost.deforestation, reached.areas)
    
    carbon.restoration <- mask(carbon.restoration, reached.areas)
    cost.restoration <- mask(cost.restoration, reached.areas)
    
    # Calculate the remaining budget
    remaining.budget <- round(r*1000000 - (sum(total.cost.deforestation, total.cost.restoration, total.cost.degradation)),0)
    
    cat("\n\t> The remaining budget is", remaining.budget, "and this is the round", r)
    
    # Calculate the maximum carbon stock for each scenario
    max.carbon.degradation <- max(carbon.degradation[], na.rm = T)
    max.carbon.deforestation <- max(carbon.deforestation[], na.rm = T)
    max.carbon.restoration <- max(carbon.restoration[], na.rm = T)
        
    # Determine which scenario has the highest maximum carbon stock
    max.carbon <- max(max.carbon.deforestation, max.carbon.restoration, max.carbon.degradation)
    
    # Select the cell with the highest carbon stock based on the scenario
    if (max.carbon == max.carbon.degradation) {
      cell.to.select <- which(values(carbon.degradation)==max.carbon)
      ncell.degradation <- c(ncell.degradation, cell.to.select)
      total.cost.degradation <- sum(total.cost.degradation, cost.degradation[cell.to.select], na.rm = T)
    } else if (max.carbon == max.carbon.deforestation) {
      cell.to.select <- which(values(carbon.deforestation)==max.carbon)
      ncell.deforestation <- c(ncell.deforestation, cell.to.select)
      total.cost.deforestation <- sum(total.cost.deforestation, cost.deforestation[cell.to.select], na.rm = T)
    } else {
      cell.to.select <- which(values(carbon.restoration)==max.carbon)
      ncell.restoration <- c(ncell.restoration, cell.to.select)
      total.cost.restoration <- sum(total.cost.restoration, cost.restoration[cell.to.select], na.rm = T)
    }
    
    # Mark the selected cell as reached
    reached.areas[cell.to.select] <- NA
  }
  
  
  # Append results to the dataframe
  result_df <- rbind(result_df, data.frame(Budget = r*1000000,
                                           Total_area = length(which(is.na(values(reached.areas)))),
                                           Degradation_area = length(ncell.degradation), 
                                           Deforestation_area = length(ncell.deforestation), 
                                           Restoration_area = length(ncell.restoration)))
}







# Create a plot
result_df %>% mutate(Total = Total_area - lag(Total_area, default = first(Total_area)),
                     Degradation = Degradation_area - lag(Degradation_area, default = first(Degradation_area)),
                     Deforestation = Deforestation_area - lag(Deforestation_area, default = first(Deforestation_area)),
                     Restoration = Restoration_area - lag(Restoration_area, default = first(Restoration_area)),
                     Proportion_Degradation = Degradation/Total,
                     Proportion_Deforestation = Deforestation/Total,
                     Proportion_Restoration = Restoration/Total) %>% 
  ggplot(aes(x = Budget)) +
  geom_line(aes(y = Proportion_Degradation, color = "Avoid Degradation"), size = 1) +
  geom_line(aes(y = Proportion_Deforestation, color = "Avoid Deforestation"), size = 1) +
  geom_line(aes(y = Proportion_Restoration, color = "Passive Restoration"), size = 1) +
  labs(x = "Budget (Brazilian Reais)", y = "Proportion of Area") +
  scale_color_manual(values = c("Avoid Deforestation" = "blue", "Passive Restoration" = "green", "Avoid Degradation" = "red")) +
  theme_minimal() +
  theme(legend.title = element_blank())






par(mfrow = c(1, 3))







































