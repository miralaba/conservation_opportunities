
#' @title Cost-benefit of conservation actions in Amazon
#' @description script to build exploratory variables from 
#' land use - land cover, secondary forest, edge, 
#' degradation (logging and fire), temperature, precipitation, 
#' elevation and distances to road and water body
#' in Santarem, Belterra and Mojui dos Campos municipalities - PA;
#' this set of exploratory variables is used in fit 
#' regional species distribution models 

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


# creating directories =========================================================
dir.create("rasters/STM/input/LULC", recursive = T)



# importing raw rasters ========================================================
## standard projection
std.proj <- "+proj=longlat +datum=WGS84 +units=m +no_defs"

## shapefile paragominas
stm.shp <- readOGR(dsn = "shapes", layer = "Santarem")
#proj4string(stm.shp) <- CRS(std.proj)
stm.shp <- spTransform(stm.shp, crs(std.proj))
stm.shp2 <- gBuffer(stm.shp, width = 0.02700003)

# 100m resolution raster
#stm.lulc.100 <- raster("rasters/STM/raw/stm-lulc-mapbiomas-brazil-collection-80-2022-100res.tif")
#
#



## land use land cover from mapbiomas collection 7, 30m res [2010 and 2020]

stm.lulc <- stack(c("rasters/STM/raw/stm-lulc-mapbiomas-brazil-collection-80-2010.tif",
                    "rasters/STM/raw/stm-lulc-mapbiomas-brazil-collection-80-2020.tif"))
names(stm.lulc) <- c("stm.lulc.2010real", "stm.lulc.2020real")

# resample to 1ha resolution
#stm.lulc <- resample(stm.lulc, stm.lulc.100, method='ngb')
#checking
#st_crs(stm.lulc)==st_crs(stm.shp)
#sort(unique(values(stm.lulc[["stm.lulc.2020real"]])))

# land ues land cover pixel values and codes
# 0  == NA
# 3  == Forest Formation      == Forest
# 4  == Savanna Formation     == Forest
# 6  == Floodable Forest      == Forest
# 11 == Wetland               == Non Forest Natural Formation
# 12 == Grassland             == Non Forest Natural Formation
# 15 == Pasture               == Farming
# 24 == Urban Area            == Non vegetated area
# 33 == River, Lake and Ocean == Water
# 39 == Soybean               == Farming
# 41 == Other Temporary Crops == Farming

#stm.lulc.2020.df <- as.data.frame(stm.lulc[["stm.lulc.2020real"]], xy = TRUE)
#breakpoints <- sort(unique(stm.lulc.2020.df$stm.lulc.2020real))
#labels.legend <- c("Non Observed", "Forest Formation", "Savanna Formation", 
#                   "Floodable Forest", "Wetland", "Grassland",
#                   "Pasture", "Urban Area", "Water", 
#                   "Soybean", "Other Temporary Crops")
#mapbiomas.legend <- c("#ffffff", "#1f8d49", "#7dc975", 
#                      "#026975", "#519799", "#d6bc74", 
#                      "#edde8e", "#d4271e", "#2532e4", 
#                      "#f5b3c8", "#f54ca9")
#
#ggplot() +
#  geom_raster(data = stm.lulc.2020.df , 
#               aes(x = x, y = y, fill = factor(stm.lulc.2020real))) + 
#  scale_fill_manual(breaks = breakpoints, 
#                    values = mapbiomas.legend, 
#                    labels = labels.legend, 
#                    name = "LULC Classes") +
#  theme_void()


### isolating forest class pixels
stm.lulc.2010.forest.class <- stm.lulc[["stm.lulc.2010real"]]
stm.lulc.2010.forest.class[stm.lulc.2010.forest.class==3] <- 1
stm.lulc.2010.forest.class[stm.lulc.2010.forest.class>1] <- 0


stm.lulc.2010.forest.mask <- stm.lulc.2010.forest.class
stm.lulc.2010.forest.mask[stm.lulc.2010.forest.mask==0] <- NA


stm.lulc.2020.forest.class <- stm.lulc[["stm.lulc.2020real"]]
stm.lulc.2020.forest.class[stm.lulc.2020.forest.class==3] <- 1
stm.lulc.2020.forest.class[stm.lulc.2020.forest.class>1] <- 0


stm.lulc.2020.forest.mask <- stm.lulc.2020.forest.class
stm.lulc.2020.forest.mask[stm.lulc.2020.forest.mask==0] <- NA
#
#



### isolating deforestation class pixels (crops, pasture)
deforestation.class.list <- c(15,39,41,48)

stm.lulc.2010.deforestation.class <- stm.lulc[["stm.lulc.2010real"]]
stm.lulc.2010.deforestation.class[stm.lulc.2010.deforestation.class %in% deforestation.class.list] <- 1
stm.lulc.2010.deforestation.class[stm.lulc.2010.deforestation.class>1] <- 0


stm.lulc.2010.deforestation.mask <- stm.lulc.2010.deforestation.class
stm.lulc.2010.deforestation.mask[stm.lulc.2010.deforestation.mask==0] <- NA


stm.lulc.2020.deforestation.class <- stm.lulc[["stm.lulc.2020real"]]
stm.lulc.2020.deforestation.class[stm.lulc.2020.deforestation.class %in% deforestation.class.list] <- 1
stm.lulc.2020.deforestation.class[stm.lulc.2020.deforestation.class>1] <- 0


stm.lulc.2020.deforestation.mask <- stm.lulc.2020.deforestation.class
stm.lulc.2020.deforestation.mask[stm.lulc.2020.deforestation.mask==0] <- NA



### select pixels based on proximity to forest
stm.deforest <- stm.lulc[["stm.lulc.2010real"]]
stm.deforest[deforest %in% deforestation.class.list] <- NA
stm.deforest[deforest>1] <- 1

stm.deforest.dist <- distance(stm.deforest, doEdge=T)
#writeRaster(stm.deforest.dist, "rasters/STM/raw/stm-distance-to-forest.tif", format = "GTiff", overwrite = T)
#stm.deforest.dist <- raster("rasters/STM/raw/stm-distance-to-forest.tif")

stm.deforest.dist.copy <- stm.deforest.dist
#values(stm.deforest.dist)[values(stm.deforest.dist) == 0] <- NA
#values(stm.deforest.dist)[values(stm.deforest.dist) > 1000] <- NA
#values(stm.deforest.dist)[values(stm.deforest.dist) <= 1000] <- 1

rm(stm.deforest)
gc()
#
#





## time since degradation
#' @description mapbiomas fire is a time-series data from 1985 to present
#' degrad is the first monitoring system to detect degradation in brazil
#' it was functional from 2007 to 2016 but does not differentiate between classes of degradation.
#' deter is the current monitoring system, with data from 2016
#' see auxiliar.R script for details about combining these data sources

stm.degrad.pf <- stack(c("rasters/STM/raw/stm-2010-dpf-tsince0.tif",
                         "rasters/STM/raw/stm-2020-dpf-tsince0.tif"))
names(stm.degrad.pf) <- c("stm.degrad.2010real", "stm.degrad.2020real")


#checking
#st_crs(stm.degrad.pf)==st_crs(stm.shp)
#plot(stm.degrad.pf)
#range(values(stm.degrad.pf[["stm.degrad.2010real"]]), na.rm=T)

### non-degraded sites will be considered with 300 years following (BIB)
stm.degrad.pf[["stm.degrad.2010real"]][stm.degrad.pf[["stm.degrad.2010real"]]>25] <- 300
stm.degrad.pf[["stm.degrad.2020real"]][stm.degrad.pf[["stm.degrad.2020real"]]>35] <- 300

### isolating degraded primary forest class pixels
stm.degrad.2010.forest.class <- stm.degrad.pf[["stm.degrad.2010real"]]
#stm.degrad.2010.forest.class <- mask(stm.degrad.2010.forest.class, stm.lulc.2010.forest.mask)
stm.degrad.2010.forest.class[stm.degrad.2010.forest.class>25]<-NA
stm.degrad.2010.forest.class[stm.degrad.2010.forest.class<=25]<-1

stm.degrad.2010.mask <- stm.degrad.2010.forest.class

stm.degrad.2010.forest.class[is.na(stm.degrad.2010.forest.class)]<-0



stm.degrad.2020.forest.class <- stm.degrad.pf[["stm.degrad.2020real"]]
#stm.degrad.2020.forest.class <- mask(stm.degrad.2020.forest.class, stm.lulc.2020.forest.mask)
stm.degrad.2020.forest.class[stm.degrad.2020.forest.class>35]<-NA
stm.degrad.2020.forest.class[stm.degrad.2020.forest.class<=35]<-1

stm.degrad.2020.mask <- stm.degrad.2020.forest.class

stm.degrad.2020.forest.class[is.na(stm.degrad.2020.forest.class)]<-0


### accounting for repeated degradation
stm.repeateddegrad.2010.mask <- raster("rasters/STM/raw/stm-degfreq-1985_2010.tif")
stm.repeateddegrad.2010.mask[stm.repeateddegrad.2010.mask<2]<-NA
stm.repeateddegrad.2010.mask[stm.repeateddegrad.2010.mask>=2]<-1


stm.repeateddegrad.2020.mask <- raster("rasters/STM/raw/stm-degfreq-1985_2020.tif")
stm.repeateddegrad.2020.mask[stm.repeateddegrad.2020.mask<2]<-NA
stm.repeateddegrad.2020.mask[stm.repeateddegrad.2020.mask>=2]<-1

#
#





## secondary forest age  [2010 and 2020]
#' source: Silva Jr. et al 2020 [DOI: 10.1038/s41597-020-00600-4]

stm.sfage <- stack(c("rasters/STM/raw/stm-2010-sfage-mapbiomas-brazil-collection-60.tif",
                     "rasters/STM/raw/stm-2020-sfage-mapbiomas-brazil-collection-60.tif"))
names(stm.sfage) <- c("stm.sfage.2010real", "stm.sfage.2020real")

#checking
#st_crs(stm.sfage)==st_crs(stm.shp)
#plot(stm.sfage)
#range(values(stm.sfage[["stm.sfage.2010real"]]), na.rm = T)

# Conversion of rasters into same extent
#stm.sfage <- resample(stm.sfage, stm.lulc.100, method='ngb')

### excluding non-forest areas
stm.sfage[["stm.sfage.2010real"]] <- mask(stm.sfage[["stm.sfage.2010real"]], stm.lulc.2010.forest.mask)
stm.sfage[["stm.sfage.2010real"]][is.na(stm.sfage[["stm.sfage.2010real"]])] <- 0


stm.sfage[["stm.sfage.2020real"]] <- mask(stm.sfage[["stm.sfage.2020real"]], stm.lulc.2020.forest.mask)
stm.sfage[["stm.sfage.2020real"]][is.na(stm.sfage[["stm.sfage.2020real"]])] <- 0



### isolating secondary forest class pixels
stm.sfage.2010.all.class <- stm.sfage[["stm.sfage.2010real"]]
stm.sfage.2010.all.class[stm.sfage.2010.all.class>0] <- 1
stm.sfage.2010.all.class[stm.sfage.2010.all.class<1] <- 0

stm.sfage.2010.mask <- stm.sfage.2010.all.class
stm.sfage.2010.mask[stm.sfage.2010.mask==0] <- NA

stm.sfage.2020.all.class <- stm.sfage[["stm.sfage.2020real"]]
stm.sfage.2020.all.class[stm.sfage.2020.all.class>0] <- 1
stm.sfage.2020.all.class[stm.sfage.2020.all.class<1] <- 0

stm.sfage.2020.mask <- stm.sfage.2020.all.class
stm.sfage.2020.mask[stm.sfage.2020.mask==0] <- NA


## degraded secondary forest
stm.degrad.sf <- stack(c("rasters/STM/raw/stm-2010-dsf-tsince0.tif",
                         "rasters/STM/raw/stm-2020-dsf-tsince0.tif"))
names(stm.degrad.sf) <- c("stm.degradsf.2010real", "stm.degradsf.2020real")


### isolating degraded secondary forest class pixels
stm.dsf.2010.mask <- stm.degrad.sf[["stm.degradsf.2010real"]]
stm.dsf.2010.mask[stm.dsf.2010.mask>25]<-NA
stm.dsf.2010.mask[stm.dsf.2010.mask<=25]<-1


stm.dsf.2020.mask <- stm.degrad.sf[["stm.degradsf.2020real"]]
stm.dsf.2020.mask[stm.dsf.2020.mask>35]<-NA
stm.dsf.2020.mask[stm.dsf.2020.mask<=35]<-1



#
#



# candidate areas for restoration scenarios ==============|
#' @description [...] (provide a brief description)
#' 
#' see auxiliar.R script for details about the steps to select candidate areas for restoration
candidate.areas.final <- raster("rasters/STM/raw/stm_forest_cover_after_restoration_2010.tif")
candidate.areas.final <- sum(stm.lulc.2010.forest.mask, candidate.areas.final, na.rm = T)
candidate.areas.final[candidate.areas.final!=1] <- 0



#
#



# building real scenarios ======================================================
##Undegraded primary forest
###2010
UPF2010<-stm.lulc.2010.forest.class
UPF2010<-mask(UPF2010, stm.sfage.2010.mask, inverse=T)
UPF2010[is.na(UPF2010[])]<-0
UPF2010<-mask(UPF2010, stm.degrad.2010.mask, inverse=T)
UPF2010[is.na(UPF2010[])]<-0
#plot(UPF2010)
writeRaster(UPF2010, "rasters/STM/input/LULC/UPF2010_real.tif", format="GTiff", overwrite=T)
#UPF2010 <- raster("rasters/STM/input/LULC/UPF2010_real.tif")

###2020
UPF2020<-stm.lulc.2020.forest.class
UPF2020<-mask(UPF2020, stm.sfage.2020.mask, inverse=T)
UPF2020[is.na(UPF2020[])]<-0
UPF2020<-mask(UPF2020, stm.degrad.2020.mask, inverse=T)
UPF2020[is.na(UPF2020[])]<-0
#plot(UPF2020)
writeRaster(UPF2020, "rasters/STM/input/LULC/UPF2020_real.tif", format="GTiff", overwrite=T)
#UPF2020 <- raster("rasters/STM/input/LULC/UPF2020_real.tif")


##Degraded primary forest
###2010
DPF2010 <- stm.degrad.2010.forest.class
uDPF2010<-mask(DPF2010, stm.repeateddegrad.2010.mask, inverse=T)
uDPF2010[is.na(uDPF2010[])]<-0
#plot(uDPF2010)
writeRaster(uDPF2010, "rasters/STM/input/LULC/uDPF2010_real.tif", format="GTiff", overwrite=T)
#uDPF2010 <- raster("rasters/STM/input/LULC/uDPF2010_real.tif")

RDPF2010 <- mask(DPF2010, stm.repeateddegrad.2010.mask)
RDPF2010[is.na(RDPF2010[])]<-0
#plot(RDPF2010)
writeRaster(RDPF2010, "rasters/STM/input/LULC/RDPF2010_real.tif", format="GTiff", overwrite=T)
#RDPF2010 <- raster("rasters/STM/input/LULC/RDPF2010_real.tif")

###2020
DPF2020 <- stm.degrad.2020.forest.class
uDPF2020<-mask(DPF2020, stm.repeateddegrad.2020.mask, inverse=T)
uDPF2020[is.na(uDPF2020[])]<-0
#plot(uDPF2020)
writeRaster(uDPF2020, "rasters/STM/input/LULC/uDPF2020_real.tif", format="GTiff", overwrite=T)
#uDPF2020 <- raster("rasters/STM/input/LULC/uDPF2020_real.tif")

RDPF2020 <- mask(DPF2020, stm.repeateddegrad.2020.mask)
RDPF2020[is.na(RDPF2020[])]<-0
#plot(RDPF2020)
writeRaster(RDPF2020, "rasters/STM/input/LULC/RDPF2020_real.tif", format="GTiff", overwrite=T)
#RDPF2020 <- raster("rasters/STM/input/LULC/RDPF2020_real.tif")



##Secondary forest
### 2010
SF2010 <- stm.sfage.2010.all.class
uSF2010 <- mask(SF2010, stm.dsf.2010.mask, inverse=T)
uSF2010[is.na(uSF2010[])]<-0
#plot(uSF2010)
writeRaster(uSF2010, "rasters/STM/input/LULC/uSF2010_real.tif", format="GTiff", overwrite=T)
#uSF2010 <- raster("rasters/STM/input/LULC/uSF2010_real.tif")

DSF2010 <- mask(SF2010, stm.dsf.2010.mask)
DSF2010[is.na(DSF2010[])]<-0
#plot(DSF2010)
writeRaster(DSF2010, "rasters/STM/input/LULC/DSF2010_real.tif", format="GTiff", overwrite=T)
#DSF2010 <- raster("rasters/STM/input/LULC/DSF2010_real.tif")

### 2020
SF2020 <- stm.sfage.2020.all.class
uSF2020 <- mask(SF2020, stm.dsf.2020.mask, inverse=T)
uSF2020[is.na(uSF2020[])]<-0
#plot(uSF2020)
writeRaster(uSF2020, "rasters/STM/input/LULC/uSF2020_real.tif", format="GTiff", overwrite=T)
#uSF2020 <- raster("rasters/STM/input/LULC/uSF2020_real.tif")

DSF2020 <- mask(SF2020, stm.dsf.2020.mask)
DSF2020[is.na(DSF2020[])]<-0
#plot(DSF2020)
writeRaster(DSF2020, "rasters/STM/input/LULC/DSF2020_real.tif", format="GTiff", overwrite=T)
#DSF2020 <- raster("rasters/STM/input/LULC/DSF2020_real.tif")



#checking
#par(mfrow = c(2, 3))
###STM 2010
#undegraded primary forest == 1
UPF2010.sk <- UPF2010
#UPF2010.sk <- raster("rasters/STM/input/UPF2010_real.tif")
#UPF2010.sk[UPF2010.sk==1]<-1
UPF2010.sk <- mask(UPF2010.sk, stm.shp)
UPF2010.sk[UPF2010.sk[]==0] <- 333
#rasterVis::levelplot(ratify(UPF2010.sk), main="undegraded primary forest", col.regions=c("white","darkgreen"), att='ID', colorkey=F)

#degraded primary forest == 10
uDPF2010.sk <- uDPF2010
#uDPF2010.sk <- raster("rasters/STM/input/uDPF2010_real.tif")
uDPF2010.sk[uDPF2010.sk==1]<-10
uDPF2010.sk <- mask(uDPF2010.sk, stm.shp)
uDPF2010.sk[uDPF2010.sk[]==0] <- 333
#plot(uDPF2010.sk, main="degraded primary forest", legend=F)

#repeated degraded primary forest == 25
RDPF2010.sk <- RDPF2010
#RDPF2010.sk <- raster("rasters/STM/input/RDPF2010_real.tif")
RDPF2010.sk[RDPF2010.sk==1]<-25
RDPF2010.sk <- mask(RDPF2010.sk, stm.shp)
RDPF2010.sk[RDPF2010.sk[]==0] <- 333
#plot(RDPF2010.sk, main="degraded primary forest", legend=F)

#secondary forest == 100
uSF2010.sk <- uSF2010
#uSF2010.sk <- raster("rasters/STM/input/uSF2010_real.tif")
uSF2010.sk[uSF2010.sk==1]<-100
uSF2010.sk <- mask(uSF2010.sk, stm.shp)
uSF2010.sk[uSF2010.sk[]==0] <- 333
#plot(uSF2010.sk, main="secondary forest", legend=F)

#degraded secondary forest == 125
DSF2010.sk <- DSF2010
#DSF2010.sk <- raster("rasters/STM/input/DSF2010_real.tif")
DSF2010.sk[DSF2010.sk==1]<-125
DSF2010.sk <- mask(DSF2010.sk, stm.shp)
DSF2010.sk[DSF2010.sk[]==0] <- 333
#plot(DSF2010.sk, main="secondary forest", legend=F)

LULC2010 <- sum(UPF2010.sk, uSF2010.sk, na.rm = T)
LULC2010 <- sum(LULC2010, DSF2010.sk, na.rm = T)
LULC2010 <- sum(LULC2010, uDPF2010.sk, na.rm = T)
LULC2010 <- sum(LULC2010, RDPF2010.sk, na.rm = T)
LULC2010 <- mask(LULC2010, stm.shp)
#sort(unique(LULC2010[]))
LULC2010[LULC2010==1333]<-1
LULC2010[LULC2010==1342]<-10
LULC2010[LULC2010==1357]<-25
LULC2010[LULC2010==1432]<-100
LULC2010[LULC2010==1457]<-125
LULC2010[LULC2010==1665 ]<-0
#plot(LULC2010, main="forest cover 2010", legend=F)

###STM 2020
#undegraded primary forest == 1
UPF2020.sk <- UPF2020
#UPF2020.sk <- raster("rasters/STM/input/UPF2020_real.tif")
#UPF2020.sk[UPF2020.sk==1]<-1
UPF2020.sk <- mask(UPF2020.sk, stm.shp)
UPF2020.sk[UPF2020.sk[]==0] <- 333
#rasterVis::levelplot(ratify(UPF2020.sk), main="undegraded primary forest", col.regions=c("white","darkgreen"), att='ID', colorkey=F)

#degraded primary forest == 10
uDPF2020.sk <- uDPF2020
#uDPF2020.sk <- raster("rasters/STM/input/uDPF2020_real.tif")
uDPF2020.sk[uDPF2020.sk==1]<-10
uDPF2020.sk <- mask(uDPF2020.sk, stm.shp)
uDPF2020.sk[uDPF2020.sk[]==0] <- 333
#plot(uDPF2020.sk, main="degradded primary forest", legend=F)

#repeated degraded primary forest == 25
RDPF2020.sk <- RDPF2020
#RDPF2020.sk <- raster("rasters/STM/input/RDPF2020_real.tif")
RDPF2020.sk[RDPF2020.sk==1]<-25
RDPF2020.sk <- mask(RDPF2020.sk, stm.shp)
RDPF2020.sk[RDPF2020.sk[]==0] <- 333
#plot(RDPF2020.sk, main="degradded primary forest", legend=F)

#secondary forest == 100
uSF2020.sk <- uSF2020
#uSF2020.sk <- raster("rasters/STM/input/uSF2020_real.tif")
uSF2020.sk[uSF2020.sk==1]<-100
uSF2020.sk <- mask(uSF2020.sk, stm.shp)
uSF2020.sk[uSF2020.sk[]==0] <- 333
#plot(uSF2020.sk, main="secondary forest", legend=F)

#degraded secondary forest == 125
DSF2020.sk <- DSF2020
#DSF2020.sk <- raster("rasters/STM/input/DSF2020_real.tif")
DSF2020.sk[DSF2020.sk==1]<-125
DSF2020.sk <- mask(DSF2020.sk, stm.shp)
DSF2020.sk[DSF2020.sk[]==0] <- 333
#plot(DSF2020.sk, main="secondary forest", legend=F)

LULC2020 <- sum(UPF2020.sk, uSF2020.sk, na.rm = T)
LULC2020 <- sum(LULC2020, DSF2020.sk, na.rm = T)
LULC2020 <- sum(LULC2020, uDPF2020.sk, na.rm = T)
LULC2020 <- sum(LULC2020, RDPF2020.sk, na.rm = T)
LULC2020 <- mask(LULC2020, stm.shp)
#sort(unique(LULC2020[]))
LULC2020[LULC2020==1333]<-1
LULC2020[LULC2020==1342]<-10
LULC2020[LULC2020==1357]<-25
LULC2020[LULC2020==1432]<-100
LULC2020[LULC2020==1457]<-125
LULC2020[LULC2020==1665 ]<-0
#plot(LULC2020, main="forest cover 2020", legend=F)




# create a data frame with the transition data
data_df <- data.frame(
  Period1 = factor(LULC2010[], levels = c(1, 10, 25, 125, 100, 0), 
                   labels = c("UPF", "DPF", "RDPF", "DSF", "SF", "D")),
  Period2 = factor(LULC2020[], levels = c(1, 10, 25, 125, 100, 0), 
                   labels = c("UPF", "DPF", "RDPF", "DSF", "SF", "D"))
)

# Create the alluvial plot
#par(mfrow = c(1, 1))
# Loading required package:
library(ggalluvial)
library(ggfittext)
library(scales)

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

#change the axis2 to check lulc transition in each scenario
set.seed(1237)
data_df %>% drop_na() %>% sample_n(size = 100000, replace = T) %>%
  ggplot(aes(axis1 = Period1, axis2 = Period2)) +
  geom_flow(aes(fill = Period2), width = .15, curve_type = "quintic") +
  geom_stratum(width = .15) +
  scale_x_discrete(limits = c("Period1", "Period2"), 
                   breaks=c("Period1", "Period2"), 
                   labels=addline_format(c("2010 Real", "2020 real")),
                   expand = c(.05, .05)) +
  scale_fill_manual(values = c("#294B29", "#50623A", "#76453B", 
                               "#B19470", "#789461", "#F97B22")) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal()+
  theme(axis.text.y= element_blank(), legend.position = "none")


rm(list=ls()[ls() %in% grep("\\b.class\\b", ls(), value = T)])
rm(list=ls()[ls() %in% grep("\\b.sk\\b", ls(), value = T)])
rm(list=ls()[ls() %in% grep("LULC", ls(), value = T)])
gc()


#making adjustments against impossible situations
LULC2020.mask <- LULC2020
((length(LULC2020.mask[LULC2010[]==0 & LULC2020[]==1])+         #DF becoming UPF
    length(LULC2020.mask[LULC2010[]==0 & LULC2020[]==10])+      #DF becoming DPF
    length(LULC2020.mask[LULC2010[]==0 & LULC2020[]==25])+      #DF becoming RDPF
    length(LULC2020.mask[LULC2010[]==100 & LULC2020[]==1])+     #SF becoming UPF
    length(LULC2020.mask[LULC2010[]==100 & LULC2020[]==10])+    #SF becoming DPF
    length(LULC2020.mask[LULC2010[]==100 & LULC2020[]==25])+    #SF becoming RDPF
    length(LULC2020.mask[LULC2010[]==125 & LULC2020[]==1])+     #DSF becoming UPF
    length(LULC2020.mask[LULC2010[]==125 & LULC2020[]==10])+    #DSF becoming DPF
    length(LULC2020.mask[LULC2010[]==125 & LULC2020[]==25])+    #DSF becoming RDPF
    length(LULC2020.mask[LULC2010[]==125 & LULC2020[]==100])+   #DSF becoming SF
    length(LULC2020.mask[LULC2010[]==25 & LULC2020[]==1])+      #RDPF becoming UPF
    length(LULC2020.mask[LULC2010[]==25 & LULC2020[]==10])+     #RDPF becoming DPF
    length(LULC2020.mask[LULC2010[]==10 & LULC2020[]==1]))/     #DPF becoming UPF
    length(LULC2020.mask[!is.na(LULC2020.mask[])]))*100
#~0.4% of pixels

LULC2020.mask[] <- ifelse(LULC2010[]==0 & LULC2020[]==1, 1, NA)
LULC2020.mask[] <- ifelse(LULC2010[]==0 & LULC2020[]==10, 1, LULC2020.mask[])
LULC2020.mask[] <- ifelse(LULC2010[]==0 & LULC2020[]==25, 1, LULC2020.mask[])
LULC2020.mask[] <- ifelse(LULC2010[]==100 & LULC2020[]==1, 1, LULC2020.mask[])
LULC2020.mask[] <- ifelse(LULC2010[]==100 & LULC2020[]==10, 1, LULC2020.mask[])
LULC2020.mask[] <- ifelse(LULC2010[]==100 & LULC2020[]==25, 1, LULC2020.mask[])
LULC2020.mask[] <- ifelse(LULC2010[]==125 & LULC2020[]==1, 1, LULC2020.mask[])
LULC2020.mask[] <- ifelse(LULC2010[]==125 & LULC2020[]==10, 1, LULC2020.mask[])
LULC2020.mask[] <- ifelse(LULC2010[]==125 & LULC2020[]==25, 1, LULC2020.mask[])
LULC2020.mask[] <- ifelse(LULC2010[]==125 & LULC2020[]==100, 1, LULC2020.mask[])
LULC2020.mask[] <- ifelse(LULC2010[]==25 & LULC2020[]==1, 1, LULC2020.mask[])
LULC2020.mask[] <- ifelse(LULC2010[]==25 & LULC2020[]==10, 1, LULC2020.mask[])
LULC2020.mask[] <- ifelse(LULC2010[]==10 & LULC2020[]==1, 1, LULC2020.mask[])
plot(LULC2020.mask, col="black", legend=F)
plot(stm.shp, add=T)

#excluding impossible situations from scenarios
#SF2010.copy <- SF2010
SF2010 <- mask(SF2010, LULC2020.mask, inverse=T)
SF2010[is.na(SF2010[])]<-0

#uSF2010.copy <- uSF2010
uSF2010 <- mask(uSF2010, LULC2020.mask, inverse=T)
uSF2010[is.na(uSF2010[])]<-0

#DSF2010.copy <- DSF2010
DSF2010 <- mask(DSF2010, LULC2020.mask, inverse=T)
DSF2010[is.na(DSF2010[])]<-0

#SF2020.copy <- SF2020
SF2020 <- mask(SF2020, LULC2020.mask, inverse=T)
SF2020[is.na(SF2020[])]<-0

#uSF2020.copy <- uSF2020
uSF2020 <- mask(uSF2020, LULC2020.mask, inverse=T)
uSF2020[is.na(uSF2020[])]<-0

#DSF2020.copy <- DSF2020
DSF2020 <- mask(DSF2020, LULC2020.mask, inverse=T)
DSF2020[is.na(DSF2020[])]<-0

#DPF2010.copy <- DPF2010
DPF2010 <- mask(DPF2010, LULC2020.mask, inverse=T)
DPF2010[is.na(DPF2010[])]<-0

#uDPF2010.copy <- uDPF2010
uDPF2010 <- mask(uDPF2010, LULC2020.mask, inverse=T)
uDPF2010[is.na(uDPF2010[])]<-0

#RDPF2010.copy <- RDPF2010
RDPF2010 <- mask(RDPF2010, LULC2020.mask, inverse=T)
RDPF2010[is.na(RDPF2010[])]<-0

#DPF2020.copy <- DPF2020
DPF2020 <- mask(DPF2020, LULC2020.mask, inverse=T)
DPF2020[is.na(DPF2020[])]<-0

#uDPF2020.copy <- uDPF2020
uDPF2020 <- mask(uDPF2020, LULC2020.mask, inverse=T)
uDPF2020[is.na(uDPF2020[])]<-0

#RDPF2020.copy <- RDPF2020
RDPF2020 <- mask(RDPF2020, LULC2020.mask, inverse=T)
RDPF2020[is.na(RDPF2020[])]<-0

#UPF2010.copy <- UPF2010
UPF2010 <- mask(UPF2010, LULC2020.mask, inverse=T)
UPF2010[is.na(UPF2010[])]<-0

#UPF2020.copy <- UPF2020
UPF2020 <- mask(UPF2020, LULC2020.mask, inverse=T)
UPF2020[is.na(UPF2020[])]<-0

#candidate.areas.final.copy <- candidate.areas.final
candidate.areas.final <- mask(candidate.areas.final, LULC2020.mask, inverse=T)
candidate.areas.final[is.na(candidate.areas.final[])]<-0


#after the adjustments check again the lulc transiotion plot





# building contractual  scenarios ==============================================
##2020 avoid degradation (all)
###Undegraded primary forest
UPF2020_avoiddegrad <- UPF2020
UPF2020_avoiddegrad[] <- ifelse(UPF2010[]==1 & uDPF2020[]==1, 1, UPF2020_avoiddegrad[])
UPF2020_avoiddegrad[] <- ifelse(UPF2010[]==1 & RDPF2020[]==1, 1, UPF2020_avoiddegrad[])
#plot(UPF2020_avoiddegrad)
writeRaster(UPF2020_avoiddegrad, "rasters/STM/input/LULC/UPF2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#UPF2020_avoiddegrad <- raster("rasters/STM/input/LULC/UPF2020_avoiddegrad.tif")


###Degraded primary forest
uDPF2020_avoiddegrad <- uDPF2020
uDPF2020_avoiddegrad[] <- ifelse(uDPF2020[]==1 & UPF2010[]==1, 0, uDPF2020_avoiddegrad[])
uDPF2020_avoiddegrad[] <- ifelse(uDPF2010[]==1 & RDPF2020[]==1, 1, uDPF2020_avoiddegrad[])
#plot(uDPF2020_avoiddegrad)
writeRaster(uDPF2020_avoiddegrad, "rasters/STM/input/LULC/uDPF2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#uDPF2020_avoiddegrad <- raster("rasters/STM/input/LULC/uDPF2020_avoiddegrad.tif")

RDPF2020_avoiddegrad <- RDPF2020
RDPF2020_avoiddegrad[] <- ifelse(RDPF2020[]==1 & UPF2010[]==1, 0, RDPF2020_avoiddegrad[])
RDPF2020_avoiddegrad[] <- ifelse(RDPF2020[]==1 & uDPF2010[]==1, 0, RDPF2020_avoiddegrad[])
#plot(RDPF2020_avoiddegrad)
writeRaster(RDPF2020_avoiddegrad, "rasters/STM/input/LULC/RDPF2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#RDPF2020_avoiddegrad <- raster("rasters/STM/input/LULC/RDPF2020_avoiddegrad.tif")


###Secondary forest
uSF2020_avoiddegrad <- uSF2020
uSF2020_avoiddegrad[] <- ifelse(uSF2010[]==1 & DSF2020[]==1, 1, uSF2020_avoiddegrad[])
#plot(uSF2020_avoiddegrad)
writeRaster(uSF2020_avoiddegrad, "rasters/STM/input/LULC/uSF2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#uSF2020_avoiddegrad <- raster("rasters/STM/input/LULC/uSF2020_avoiddegrad.tif")

DSF2020_avoiddegrad <- DSF2020
DSF2020_avoiddegrad[] <- ifelse(DSF2020[]==1 & uSF2010[]==1, 0, DSF2020_avoiddegrad[])
#plot(DSF2020_avoiddegrad)
writeRaster(DSF2020_avoiddegrad, "rasters/STM/input/LULC/DSF2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#DSF2020_avoiddegrad <- raster("rasters/STM/input/LULC/DSF2020_avoiddegrad.tif")




##2020 avoid degradation (Primary forest only)
###Undegraded primary forest
UPF2020_avoiddegrad2 <- UPF2020
UPF2020_avoiddegrad2[] <- ifelse(UPF2010[]==1 & uDPF2020[]==1, 1, UPF2020_avoiddegrad2[])
UPF2020_avoiddegrad2[] <- ifelse(UPF2010[]==1 & RDPF2020[]==1, 1, UPF2020_avoiddegrad2[])
#plot(UPF2020_avoiddegrad2)
writeRaster(UPF2020_avoiddegrad2, "rasters/STM/input/LULC/UPF2020_avoiddegrad2.tif", format="GTiff", overwrite=T)
#UPF2020_avoiddegrad2 <- raster("rasters/STM/input/LULC/UPF2020_avoiddegrad2.tif")


###Degraded primary forest
uDPF2020_avoiddegrad2 <- uDPF2020
uDPF2020_avoiddegrad2[] <- ifelse(uDPF2020[]==1 & UPF2010[]==1, 0, uDPF2020_avoiddegrad2[])
uDPF2020_avoiddegrad2[] <- ifelse(uDPF2010[]==1 & RDPF2020[]==1, 1, uDPF2020_avoiddegrad2[])
#plot(uDPF2020_avoiddegrad2)
writeRaster(uDPF2020_avoiddegrad2, "rasters/STM/input/LULC/uDPF2020_avoiddegrad2.tif", format="GTiff", overwrite=T)
#uDPF2020_avoiddegrad2 <- raster("rasters/STM/input/LULC/uDPF2020_avoiddegrad2.tif")

RDPF2020_avoiddegrad2 <- RDPF2020
RDPF2020_avoiddegrad2[] <- ifelse(RDPF2020[]==1 & UPF2010[]==1, 0, RDPF2020_avoiddegrad2[])
RDPF2020_avoiddegrad2[] <- ifelse(RDPF2020[]==1 & uDPF2010[]==1, 0, RDPF2020_avoiddegrad2[])
#plot(RDPF2020_avoiddegrad2)
writeRaster(RDPF2020_avoiddegrad2, "rasters/STM/input/LULC/RDPF2020_avoiddegrad2.tif", format="GTiff", overwrite=T)
#RDPF2020_avoiddegrad2 <- raster("rasters/STM/input/LULC/RDPF2020_avoiddegrad2.tif")


###Secondary forest
uSF2020_avoiddegrad2 <- uSF2020
#plot(uSF2020_avoiddegrad2)
writeRaster(uSF2020_avoiddegrad2, "rasters/STM/input/LULC/uSF2020_avoiddegrad2.tif", format="GTiff", overwrite=T)
#uSF2020_avoiddegrad2 <- raster("rasters/STM/input/LULC/uSF2020_avoiddegrad2.tif")

DSF2020_avoiddegrad2 <- DSF2020
#plot(DSF2020_avoiddegrad2)
writeRaster(DSF2020_avoiddegrad2, "rasters/STM/input/LULC/DSF2020_avoiddegrad2.tif", format="GTiff", overwrite=T)
#DSF2020_avoiddegrad2 <- raster("rasters/STM/input/LULC/DSF2020_avoiddegrad2.tif")




##2020 avoid deforestation (all)
###Undegraded primary forest
UPF2020_avoiddeforest <- UPF2020
UPF2020_avoiddeforest[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, UPF2020_avoiddeforest[])
UPF2020_avoiddeforest[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, UPF2020_avoiddeforest[])
UPF2020_avoiddeforest[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, UPF2020_avoiddeforest[])
#plot(UPF2020_avoiddeforest)
writeRaster(UPF2020_avoiddeforest, "rasters/STM/input/LULC/UPF2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#UPF2020_avoiddeforest <- raster("rasters/STM/input/LULC/UPF2020_avoiddeforest.tif")


###Degraded primary forest
uDPF2020_avoiddeforest <- uDPF2020
uDPF2020_avoiddeforest[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, uDPF2020_avoiddeforest[])
uDPF2020_avoiddeforest[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, uDPF2020_avoiddeforest[])
uDPF2020_avoiddeforest[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, uDPF2020_avoiddeforest[])
#plot(uDPF2020_avoiddeforest)
writeRaster(uDPF2020_avoiddeforest, "rasters/STM/input/LULC/uDPF2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#uDPF2020_avoiddeforest <- raster("rasters/STM/input/LULC/uDPF2020_avoiddeforest.tif")

RDPF2020_avoiddeforest <- RDPF2020
RDPF2020_avoiddeforest[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, RDPF2020_avoiddeforest[])
RDPF2020_avoiddeforest[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, RDPF2020_avoiddeforest[])
RDPF2020_avoiddeforest[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, RDPF2020_avoiddeforest[])
#plot(RDPF2020_avoiddeforest)
writeRaster(RDPF2020_avoiddeforest, "rasters/STM/input/LULC/RDPF2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#RDPF2020_avoiddeforest <- raster("rasters/STM/input/LULC/RDPF2020_avoiddeforest.tif")


###Secondary forest
uSF2020_avoiddeforest <- uSF2020
uSF2020_avoiddeforest[] <- ifelse(uSF2020[]==1 & UPF2010[]==1, 0, uSF2020_avoiddeforest[])
uSF2020_avoiddeforest[] <- ifelse(uSF2020[]==1 & uDPF2010[]==1, 0, uSF2020_avoiddeforest[])
uSF2020_avoiddeforest[] <- ifelse(uSF2020[]==1 & RDPF2010[]==1, 0, uSF2020_avoiddeforest[])
uSF2020_avoiddeforest[] <- ifelse(DSF2020[]==0 & uSF2020[]==0 & uSF2010[]==1, 1, uSF2020_avoiddeforest[])
#plot(uSF2020_avoiddeforest)
writeRaster(uSF2020_avoiddeforest, "rasters/STM/input/LULC/uSF2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#uSF2020_avoiddeforest <- raster("rasters/STM/input/LULC/uSF2020_avoiddeforest.tif")

DSF2020_avoiddeforest <- DSF2020
DSF2020_avoiddeforest[] <- ifelse(DSF2020[]==1 & UPF2010[]==1, 0, DSF2020_avoiddeforest[])
DSF2020_avoiddeforest[] <- ifelse(DSF2020[]==1 & uDPF2010[]==1, 0, DSF2020_avoiddeforest[])
DSF2020_avoiddeforest[] <- ifelse(DSF2020[]==1 & RDPF2010[]==1, 0, DSF2020_avoiddeforest[])
DSF2020_avoiddeforest[] <- ifelse(uSF2020[]==0 & DSF2020[]==0 & DSF2010[]==1, 1, DSF2020_avoiddeforest[])
#plot(DSF2020_avoiddeforest)
writeRaster(DSF2020_avoiddeforest, "rasters/STM/input/LULC/DSF2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#DSF2020_avoiddeforest <- raster("rasters/STM/input/LULC/DSF2020_avoiddeforest.tif")




##2020 avoid deforestation (Primary forest only)
###Undegraded primary forest
UPF2020_avoiddeforest2 <- UPF2020
UPF2020_avoiddeforest2[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, UPF2020_avoiddeforest2[])
UPF2020_avoiddeforest2[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, UPF2020_avoiddeforest2[])
UPF2020_avoiddeforest2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, UPF2020_avoiddeforest2[])
#plot(UPF2020_avoiddeforest2)
writeRaster(UPF2020_avoiddeforest2, "rasters/STM/input/LULC/UPF2020_avoiddeforest2.tif", format="GTiff", overwrite=T)
#UPF2020_avoiddeforest2 <- raster("rasters/STM/input/LULC/UPF2020_avoiddeforest2.tif")


###Degraded primary forest
uDPF2020_avoiddeforest2 <- uDPF2020
uDPF2020_avoiddeforest2[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, uDPF2020_avoiddeforest2[])
uDPF2020_avoiddeforest2[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, uDPF2020_avoiddeforest2[])
uDPF2020_avoiddeforest2[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, uDPF2020_avoiddeforest2[])
#plot(uDPF2020_avoiddeforest2)
writeRaster(uDPF2020_avoiddeforest2, "rasters/STM/input/LULC/uDPF2020_avoiddeforest2.tif", format="GTiff", overwrite=T)
#uDPF2020_avoiddeforest2 <- raster("rasters/STM/input/LULC/uDPF2020_avoiddeforest2.tif")

RDPF2020_avoiddeforest2 <- RDPF2020
RDPF2020_avoiddeforest2[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, RDPF2020_avoiddeforest2[])
RDPF2020_avoiddeforest2[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, RDPF2020_avoiddeforest2[])
RDPF2020_avoiddeforest2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, RDPF2020_avoiddeforest2[])
#plot(RDPF2020_avoiddeforest2)
writeRaster(RDPF2020_avoiddeforest2, "rasters/STM/input/LULC/RDPF2020_avoiddeforest2.tif", format="GTiff", overwrite=T)
#RDPF2020_avoiddeforest2 <- raster("rasters/STM/input/LULC/RDPF2020_avoiddeforest2.tif")


###Secondary forest
uSF2020_avoiddeforest2 <- uSF2020
uSF2020_avoiddeforest2[] <- ifelse(uSF2020[]==1 & UPF2010[]==1, 0, uSF2020_avoiddeforest2[])
uSF2020_avoiddeforest2[] <- ifelse(uSF2020[]==1 & uDPF2010[]==1, 0, uSF2020_avoiddeforest2[])
uSF2020_avoiddeforest2[] <- ifelse(uSF2020[]==1 & RDPF2010[]==1, 0, uSF2020_avoiddeforest2[])
#plot(uSF2020_avoiddeforest2)
writeRaster(uSF2020_avoiddeforest2, "rasters/STM/input/LULC/uSF2020_avoiddeforest2.tif", format="GTiff", overwrite=T)
#uSF2020_avoiddeforest2 <- raster("rasters/STM/input/LULC/uSF2020_avoiddeforest2.tif")

DSF2020_avoiddeforest2 <- DSF2020
DSF2020_avoiddeforest2[] <- ifelse(DSF2020[]==1 & UPF2010[]==1, 0, DSF2020_avoiddeforest2[])
DSF2020_avoiddeforest2[] <- ifelse(DSF2020[]==1 & uDPF2010[]==1, 0, DSF2020_avoiddeforest2[])
DSF2020_avoiddeforest2[] <- ifelse(DSF2020[]==1 & RDPF2010[]==1, 0, DSF2020_avoiddeforest2[])
#plot(DSF2020_avoiddeforest2)
writeRaster(DSF2020_avoiddeforest2, "rasters/STM/input/LULC/DSF2020_avoiddeforest2.tif", format="GTiff", overwrite=T)
#DSF2020_avoiddeforest2 <- raster("rasters/STM/input/LULC/DSF2020_avoiddeforest2.tif")




##2020 restoration without avoid
###Secondary forest
uSF2020_restor_wo_avoid <- sum(uSF2020, candidate.areas.final)
uSF2020_restor_wo_avoid[] <-  ifelse(uSF2010[]==0 & uSF2020[]==1, 1, uSF2020_restor_wo_avoid[])
uSF2020_restor_wo_avoid[] <-  ifelse(uSF2020_restor_wo_avoid[]==1 & DSF2020[]==1, 0, uSF2020_restor_wo_avoid[])
uSF2020_restor_wo_avoid[] <-  ifelse(uSF2020_restor_wo_avoid[]==1 & UPF2020[]==1, 0, uSF2020_restor_wo_avoid[])
uSF2020_restor_wo_avoid[] <-  ifelse(uSF2020_restor_wo_avoid[]==1 & uDPF2020[]==1, 0, uSF2020_restor_wo_avoid[])
uSF2020_restor_wo_avoid[] <-  ifelse(uSF2020_restor_wo_avoid[]==1 & RDPF2020[]==1, 0, uSF2020_restor_wo_avoid[])
#plot(uSF2020_restor_wo_avoid)
writeRaster(uSF2020_restor_wo_avoid, "rasters/STM/input/LULC/uSF2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)
#uSF2020_restor_wo_avoid <- raster("rasters/STM/input/LULC/uSF2020_restor_wo_avoid.tif")

DSF2020_restor_wo_avoid <- DSF2020
#plot(DSF2020_restor_wo_avoid)
writeRaster(DSF2020_restor_wo_avoid, "rasters/STM/input/LULC/DSF2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)
#DSF2020_restor_wo_avoid <- raster("rasters/STM/input/LULC/DSF2020_restor_wo_avoid.tif")


###Undegraded primary forest
UPF2020_restor_wo_avoid <- UPF2020
#plot(UPF2020_restor_wo_avoid)
writeRaster(UPF2020_restor_wo_avoid, "rasters/STM/input/LULC/UPF2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)
#UPF2020_restor_wo_avoid <- raster("rasters/STM/input/LULC/UPF2020_restor_wo_avoid.tif")


###Degraded primary forest
uDPF2020_restor_wo_avoid <- uDPF2020
#plot(uDPF2020_restor_wo_avoid)
writeRaster(uDPF2020_restor_wo_avoid, "rasters/STM/input/LULC/uDPF2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)
#uDPF2020_restor_wo_avoid <- raster("rasters/STM/input/LULC/uDPF2020_restor_wo_avoid.tif")

RDPF2020_restor_wo_avoid <- RDPF2020
#plot(RDPF2020_restor_wo_avoid)
writeRaster(RDPF2020_restor_wo_avoid, "rasters/STM/input/LULC/RDPF2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)
#RDPF2020_restor_wo_avoid <- raster("rasters/STM/input/LULC/RDPF2020_restor_wo_avoid.tif")




##2020 avoid both (all)
###Undegraded primary forest
UPF2020_avoidboth <- UPF2020
UPF2020_avoidboth[] <- ifelse(UPF2010[]==1 & uDPF2020[]==1, 1, UPF2020_avoidboth[])
UPF2020_avoidboth[] <- ifelse(UPF2010[]==1 & RDPF2020[]==1, 1, UPF2020_avoidboth[])
UPF2020_avoidboth[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, UPF2020_avoidboth[])
UPF2020_avoidboth[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, UPF2020_avoidboth[])
UPF2020_avoidboth[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, UPF2020_avoidboth[])
#plot(UPF2020_avoidboth)
writeRaster(UPF2020_avoidboth, "rasters/STM/input/LULC/UPF2020_avoidboth.tif", format="GTiff", overwrite=T)
#UPF2020_avoidboth <- raster("rasters/STM/input/LULC/UPF2020_avoidboth.tif")


###Degraded primary forest
uDPF2020_avoidboth <- uDPF2020
uDPF2020_avoidboth[] <- ifelse(uDPF2020[]==1 & UPF2010[]==1, 0, uDPF2020_avoidboth[])
uDPF2020_avoidboth[] <- ifelse(uDPF2010[]==1 & RDPF2020[]==1, 1, uDPF2020_avoidboth[])
uDPF2020_avoidboth[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, uDPF2020_avoidboth[])
uDPF2020_avoidboth[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, uDPF2020_avoidboth[])
uDPF2020_avoidboth[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, uDPF2020_avoidboth[])
#plot(uDPF2020_avoidboth)
writeRaster(uDPF2020_avoidboth, "rasters/STM/input/LULC/uDPF2020_avoidboth.tif", format="GTiff", overwrite=T)
#uDPF2020_avoidboth <- raster("rasters/STM/input/LULC/uDPF2020_avoidboth.tif")

RDPF2020_avoidboth <- RDPF2020
RDPF2020_avoidboth[] <- ifelse(RDPF2020[]==1 & UPF2010[]==1, 0, RDPF2020_avoidboth[])
RDPF2020_avoidboth[] <- ifelse(RDPF2020[]==1 & uDPF2010[]==1, 0, RDPF2020_avoidboth[])
RDPF2020_avoidboth[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, RDPF2020_avoidboth[])
RDPF2020_avoidboth[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, RDPF2020_avoidboth[])
RDPF2020_avoidboth[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, RDPF2020_avoidboth[])
#plot(RDPF2020_avoidboth)
writeRaster(RDPF2020_avoidboth, "rasters/STM/input/LULC/RDPF2020_avoidboth.tif", format="GTiff", overwrite=T)
#RDPF2020_avoidboth <- raster("rasters/STM/input/LULC/RDPF2020_avoidboth.tif")


###Secondary forest
uSF2020_avoidboth <- uSF2020
uSF2020_avoidboth[] <- ifelse(uSF2010[]==1 & DSF2020[]==1, 1, uSF2020_avoidboth[])
uSF2020_avoidboth[] <- ifelse(uSF2020[]==1 & UPF2010[]==1, 0, uSF2020_avoidboth[])
uSF2020_avoidboth[] <- ifelse(uSF2020[]==1 & uDPF2010[]==1, 0, uSF2020_avoidboth[])
uSF2020_avoidboth[] <- ifelse(uSF2020[]==1 & RDPF2010[]==1, 0, uSF2020_avoidboth[])
uSF2020_avoidboth[] <- ifelse(DSF2020[]==0 & uSF2020[]==0 & uSF2010[]==1, 1, uSF2020_avoidboth[])
#plot(uSF2020_avoidboth)
writeRaster(uSF2020_avoidboth, "rasters/STM/input/LULC/uSF2020_avoidboth.tif", format="GTiff", overwrite=T)
#uSF2020_avoidboth <- raster("rasters/STM/input/LULC/uSF2020_avoidboth.tif")

DSF2020_avoidboth <- DSF2020
DSF2020_avoidboth[] <- ifelse(DSF2020[]==1 & uSF2010[]==1, 0, DSF2020_avoidboth[])
DSF2020_avoidboth[] <- ifelse(DSF2020[]==1 & UPF2010[]==1, 0, DSF2020_avoidboth[])
DSF2020_avoidboth[] <- ifelse(DSF2020[]==1 & uDPF2010[]==1, 0, DSF2020_avoidboth[])
DSF2020_avoidboth[] <- ifelse(DSF2020[]==1 & RDPF2010[]==1, 0, DSF2020_avoidboth[])
DSF2020_avoidboth[] <- ifelse(uSF2020[]==0 & DSF2020[]==0 & DSF2010[]==1, 1, DSF2020_avoidboth[])
#plot(DSF2020_avoidboth)
writeRaster(DSF2020_avoidboth, "rasters/STM/input/LULC/DSF2020_avoidboth.tif", format="GTiff", overwrite=T)
#DSF2020_avoidboth <- raster("rasters/STM/input/LULC/DSF2020_avoidboth.tif")




##2020 avoid both (Primary forest only)
###Undegraded primary forest
UPF2020_avoidboth2 <- UPF2020
UPF2020_avoidboth2[] <- ifelse(UPF2010[]==1 & uDPF2020[]==1, 1, UPF2020_avoidboth2[])
UPF2020_avoidboth2[] <- ifelse(UPF2010[]==1 & RDPF2020[]==1, 1, UPF2020_avoidboth2[])
UPF2020_avoidboth2[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, UPF2020_avoidboth2[])
UPF2020_avoidboth2[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, UPF2020_avoidboth2[])
UPF2020_avoidboth2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, UPF2020_avoidboth2[])
#plot(UPF2020_avoidboth2)
writeRaster(UPF2020_avoidboth2, "rasters/STM/input/LULC/UPF2020_avoidboth2.tif", format="GTiff", overwrite=T)
#UPF2020_avoidboth2 <- raster("rasters/STM/input/LULC/UPF2020_avoidboth2.tif")


###Degraded primary forest
uDPF2020_avoidboth2 <- uDPF2020
uDPF2020_avoidboth2[] <- ifelse(uDPF2020[]==1 & UPF2010[]==1, 0, uDPF2020_avoidboth2[])
uDPF2020_avoidboth2[] <- ifelse(uDPF2010[]==1 & RDPF2020[]==1, 1, uDPF2020_avoidboth2[])
uDPF2020_avoidboth2[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, uDPF2020_avoidboth2[])
uDPF2020_avoidboth2[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, uDPF2020_avoidboth2[])
uDPF2020_avoidboth2[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, uDPF2020_avoidboth2[])
#plot(uDPF2020_avoidboth2)
writeRaster(uDPF2020_avoidboth2, "rasters/STM/input/LULC/uDPF2020_avoidboth2.tif", format="GTiff", overwrite=T)
#uDPF2020_avoidboth2 <- raster("rasters/STM/input/LULC/uDPF2020_avoidboth2.tif")

RDPF2020_avoidboth2 <- RDPF2020
RDPF2020_avoidboth2[] <- ifelse(RDPF2020[]==1 & UPF2010[]==1, 0, RDPF2020_avoidboth2[])
RDPF2020_avoidboth2[] <- ifelse(RDPF2020[]==1 & uDPF2010[]==1, 0, RDPF2020_avoidboth2[])
RDPF2020_avoidboth2[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, RDPF2020_avoidboth2[])
RDPF2020_avoidboth2[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, RDPF2020_avoidboth2[])
RDPF2020_avoidboth2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, RDPF2020_avoidboth2[])
#plot(RDPF2020_avoidboth2)
writeRaster(RDPF2020_avoidboth2, "rasters/STM/input/LULC/RDPF2020_avoidboth2.tif", format="GTiff", overwrite=T)
#RDPF2020_avoidboth2 <- raster("rasters/STM/input/LULC/RDPF2020_avoidboth2.tif")


###Secondary forest
uSF2020_avoidboth2 <- uSF2020
uSF2020_avoidboth2[] <- ifelse(uSF2020[]==1 & UPF2010[]==1, 0, uSF2020_avoidboth2[])
uSF2020_avoidboth2[] <- ifelse(uSF2020[]==1 & uDPF2010[]==1, 0, uSF2020_avoidboth2[])
uSF2020_avoidboth2[] <- ifelse(uSF2020[]==1 & RDPF2010[]==1, 0, uSF2020_avoidboth2[])
#plot(uSF2020_avoidboth2)
writeRaster(uSF2020_avoidboth2, "rasters/STM/input/LULC/uSF2020_avoidboth2.tif", format="GTiff", overwrite=T)
#uSF2020_avoidboth2 <- raster("rasters/STM/input/LULC/uSF2020_avoidboth2.tif")

DSF2020_avoidboth2 <- DSF2020
DSF2020_avoidboth2[] <- ifelse(DSF2020[]==1 & UPF2010[]==1, 0, DSF2020_avoidboth2[])
DSF2020_avoidboth2[] <- ifelse(DSF2020[]==1 & uDPF2010[]==1, 0, DSF2020_avoidboth2[])
DSF2020_avoidboth2[] <- ifelse(DSF2020[]==1 & RDPF2010[]==1, 0, DSF2020_avoidboth2[])
#plot(DSF2020_avoidboth2)
writeRaster(DSF2020_avoidboth2, "rasters/STM/input/LULC/DSF2020_avoidboth2.tif", format="GTiff", overwrite=T)
#DSF2020_avoidboth2 <- raster("rasters/STM/input/LULC/DSF2020_avoidboth2.tif")




##2020 restoration and avoid deforestation (all)
###Secondary forest
uSF2020_restor_n_avoiddeforest <- sum(uSF2020, candidate.areas.final)
uSF2020_restor_n_avoiddeforest[] <-  ifelse(uSF2010[]==0 & uSF2020[]==1, 1, uSF2020_restor_n_avoiddeforest[])
uSF2020_restor_n_avoiddeforest[] <-  ifelse(uSF2020_restor_n_avoiddeforest[]==1 & DSF2020[]==1, 0, uSF2020_restor_n_avoiddeforest[])
uSF2020_restor_n_avoiddeforest[] <-  ifelse(uSF2020_restor_n_avoiddeforest[]==1 & UPF2020[]==1, 0, uSF2020_restor_n_avoiddeforest[])
uSF2020_restor_n_avoiddeforest[] <-  ifelse(uSF2020_restor_n_avoiddeforest[]==1 & uDPF2020[]==1, 0, uSF2020_restor_n_avoiddeforest[])
uSF2020_restor_n_avoiddeforest[] <-  ifelse(uSF2020_restor_n_avoiddeforest[]==1 & RDPF2020[]==1, 0, uSF2020_restor_n_avoiddeforest[])

uSF2020_restor_n_avoiddeforest[] <- ifelse(uSF2020[]==1 & UPF2010[]==1, 0, uSF2020_restor_n_avoiddeforest[])
uSF2020_restor_n_avoiddeforest[] <- ifelse(uSF2020[]==1 & uDPF2010[]==1, 0, uSF2020_restor_n_avoiddeforest[])
uSF2020_restor_n_avoiddeforest[] <- ifelse(uSF2020[]==1 & RDPF2010[]==1, 0, uSF2020_restor_n_avoiddeforest[])
uSF2020_restor_n_avoiddeforest[] <- ifelse(DSF2020[]==0 & uSF2020[]==0 & uSF2010[]==1, 1, uSF2020_restor_n_avoiddeforest[])
#plot(uSF2020_restor_n_avoiddeforest)
writeRaster(uSF2020_restor_n_avoiddeforest, "rasters/STM/input/LULC/uSF2020_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)
#uSF2020_restor_n_avoiddeforest <- raster("rasters/STM/input/LULC/uSF2020_restor_n_avoiddeforest.tif")

DSF2020_restor_n_avoiddeforest <- DSF2020
DSF2020_restor_n_avoiddeforest[] <- ifelse(DSF2020[]==1 & UPF2010[]==1, 0, DSF2020_restor_n_avoiddeforest[])
DSF2020_restor_n_avoiddeforest[] <- ifelse(DSF2020[]==1 & uDPF2010[]==1, 0, DSF2020_restor_n_avoiddeforest[])
DSF2020_restor_n_avoiddeforest[] <- ifelse(DSF2020[]==1 & RDPF2010[]==1, 0, DSF2020_restor_n_avoiddeforest[])
DSF2020_restor_n_avoiddeforest[] <- ifelse(uSF2020[]==0 & DSF2020[]==0 & DSF2010[]==1, 1, DSF2020_restor_n_avoiddeforest[])
#plot(DSF2020_restor_n_avoiddeforest)
writeRaster(DSF2020_restor_n_avoiddeforest, "rasters/STM/input/LULC/DSF2020_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)
#DSF2020_restor_n_avoiddeforest <- raster("rasters/STM/input/LULC/DSF2020_restor_n_avoiddeforest.tif")


###Undegraded primary forest
UPF2020_restor_n_avoiddeforest <- UPF2020
UPF2020_restor_n_avoiddeforest[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, UPF2020_restor_n_avoiddeforest[])
UPF2020_restor_n_avoiddeforest[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, UPF2020_restor_n_avoiddeforest[])
UPF2020_restor_n_avoiddeforest[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, UPF2020_restor_n_avoiddeforest[])
#plot(UPF2020_restor_n_avoiddeforest)
writeRaster(UPF2020_restor_n_avoiddeforest, "rasters/STM/input/LULC/UPF2020_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)
#UPF2020_restor_n_avoiddeforest <- raster("rasters/STM/input/LULC/UPF2020_restor_n_avoiddeforest.tif")


###Degraded primary forest
uDPF2020_restor_n_avoiddeforest <- uDPF2020
uDPF2020_restor_n_avoiddeforest[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, uDPF2020_restor_n_avoiddeforest[])
uDPF2020_restor_n_avoiddeforest[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, uDPF2020_restor_n_avoiddeforest[])
uDPF2020_restor_n_avoiddeforest[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, uDPF2020_restor_n_avoiddeforest[])
#plot(uDPF2020_restor_n_avoiddeforest)
writeRaster(uDPF2020_restor_n_avoiddeforest, "rasters/STM/input/LULC/uDPF2020_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)
#uDPF2020_restor_n_avoiddeforest <- raster("rasters/STM/input/LULC/uDPF2020_restor_n_avoiddeforest.tif")

RDPF2020_restor_n_avoiddeforest <- RDPF2020
RDPF2020_restor_n_avoiddeforest[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, RDPF2020_restor_n_avoiddeforest[])
RDPF2020_restor_n_avoiddeforest[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, RDPF2020_restor_n_avoiddeforest[])
RDPF2020_restor_n_avoiddeforest[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, RDPF2020_restor_n_avoiddeforest[])
#plot(RDPF2020_restor_n_avoiddeforest)
writeRaster(RDPF2020_restor_n_avoiddeforest, "rasters/STM/input/LULC/RDPF2020_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)
#RDPF2020_restor_n_avoiddeforest <- raster("rasters/STM/input/LULC/RDPF2020_restor_n_avoiddeforest.tif")




##2020 restoration and avoid deforestation (Primary forest only)
###Secondary forest
uSF2020_restor_n_avoiddeforest2 <- sum(uSF2020, candidate.areas.final)
uSF2020_restor_n_avoiddeforest2[] <-  ifelse(uSF2010[]==0 & uSF2020[]==1, 1, uSF2020_restor_n_avoiddeforest2[])
uSF2020_restor_n_avoiddeforest2[] <-  ifelse(uSF2020_restor_n_avoiddeforest2[]==1 & DSF2020[]==1, 0, uSF2020_restor_n_avoiddeforest2[])
uSF2020_restor_n_avoiddeforest2[] <-  ifelse(uSF2020_restor_n_avoiddeforest2[]==1 & UPF2020[]==1, 0, uSF2020_restor_n_avoiddeforest2[])
uSF2020_restor_n_avoiddeforest2[] <-  ifelse(uSF2020_restor_n_avoiddeforest2[]==1 & uDPF2020[]==1, 0, uSF2020_restor_n_avoiddeforest2[])
uSF2020_restor_n_avoiddeforest2[] <-  ifelse(uSF2020_restor_n_avoiddeforest2[]==1 & RDPF2020[]==1, 0, uSF2020_restor_n_avoiddeforest2[])

uSF2020_restor_n_avoiddeforest2[] <- ifelse(uSF2020[]==1 & UPF2010[]==1, 0, uSF2020_restor_n_avoiddeforest2[])
uSF2020_restor_n_avoiddeforest2[] <- ifelse(uSF2020[]==1 & uDPF2010[]==1, 0, uSF2020_restor_n_avoiddeforest2[])
uSF2020_restor_n_avoiddeforest2[] <- ifelse(uSF2020[]==1 & RDPF2010[]==1, 0, uSF2020_restor_n_avoiddeforest2[])
#plot(uSF2020_restor_n_avoiddeforest2)
writeRaster(uSF2020_restor_n_avoiddeforest2, "rasters/STM/input/LULC/uSF2020_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)
#uSF2020_restor_n_avoiddeforest2 <- raster("rasters/STM/input/LULC/uSF2020_restor_n_avoiddeforest2.tif")

DSF2020_restor_n_avoiddeforest2 <- DSF2020
DSF2020_restor_n_avoiddeforest2[] <- ifelse(DSF2020[]==1 & UPF2010[]==1, 0, DSF2020_restor_n_avoiddeforest2[])
DSF2020_restor_n_avoiddeforest2[] <- ifelse(DSF2020[]==1 & uDPF2010[]==1, 0, DSF2020_restor_n_avoiddeforest2[])
DSF2020_restor_n_avoiddeforest2[] <- ifelse(DSF2020[]==1 & RDPF2010[]==1, 0, DSF2020_restor_n_avoiddeforest2[])
#plot(DSF2020_restor_n_avoiddeforest2)
writeRaster(DSF2020_restor_n_avoiddeforest2, "rasters/STM/input/LULC/DSF2020_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)
#DSF2020_restor_n_avoiddeforest2 <- raster("rasters/STM/input/LULC/DSF2020_restor_n_avoiddeforest2.tif")


###Undegraded primary forest
UPF2020_restor_n_avoiddeforest2 <- UPF2020
UPF2020_restor_n_avoiddeforest2[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, UPF2020_restor_n_avoiddeforest2[])
UPF2020_restor_n_avoiddeforest2[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, UPF2020_restor_n_avoiddeforest2[])
UPF2020_restor_n_avoiddeforest2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, UPF2020_restor_n_avoiddeforest2[])
#plot(UPF2020_restor_n_avoiddeforest2)
writeRaster(UPF2020_restor_n_avoiddeforest2, "rasters/STM/input/LULC/UPF2020_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)
#UPF2020_restor_n_avoiddeforest2 <- raster("rasters/STM/input/LULC/UPF2020_restor_n_avoiddeforest2.tif")


###Degraded primary forest
uDPF2020_restor_n_avoiddeforest2 <- uDPF2020
uDPF2020_restor_n_avoiddeforest2[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, uDPF2020_restor_n_avoiddeforest2[])
uDPF2020_restor_n_avoiddeforest2[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, uDPF2020_restor_n_avoiddeforest2[])
uDPF2020_restor_n_avoiddeforest2[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, uDPF2020_restor_n_avoiddeforest2[])
#plot(uDPF2020_restor_n_avoiddeforest2)
writeRaster(uDPF2020_restor_n_avoiddeforest2, "rasters/STM/input/LULC/uDPF2020_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)
#uDPF2020_restor_n_avoiddeforest2 <- raster("rasters/STM/input/LULC/uDPF2020_restor_n_avoiddeforest2.tif")

RDPF2020_restor_n_avoiddeforest2 <- RDPF2020
RDPF2020_restor_n_avoiddeforest2[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, RDPF2020_restor_n_avoiddeforest2[])
RDPF2020_restor_n_avoiddeforest2[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, RDPF2020_restor_n_avoiddeforest2[])
RDPF2020_restor_n_avoiddeforest2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, RDPF2020_restor_n_avoiddeforest2[])
#plot(RDPF2020_restor_n_avoiddeforest2)
writeRaster(RDPF2020_restor_n_avoiddeforest2, "rasters/STM/input/LULC/RDPF2020_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)
#RDPF2020_restor_n_avoiddeforest2 <- raster("rasters/STM/input/LULC/RDPF2020_restor_n_avoiddeforest2.tif")




##2020 restoration and avoid both (all)
###Secondary forest
uSF2020_restor_n_avoidboth <- sum(uSF2020, candidate.areas.final)
uSF2020_restor_n_avoidboth[] <-  ifelse(uSF2010[]==0 & uSF2020[]==1, 1, uSF2020_restor_n_avoidboth[])
uSF2020_restor_n_avoidboth[] <-  ifelse(uSF2020_restor_n_avoidboth[]==1 & DSF2020[]==1, 0, uSF2020_restor_n_avoidboth[])
uSF2020_restor_n_avoidboth[] <-  ifelse(uSF2020_restor_n_avoidboth[]==1 & UPF2020[]==1, 0, uSF2020_restor_n_avoidboth[])
uSF2020_restor_n_avoidboth[] <-  ifelse(uSF2020_restor_n_avoidboth[]==1 & uDPF2020[]==1, 0, uSF2020_restor_n_avoidboth[])
uSF2020_restor_n_avoidboth[] <-  ifelse(uSF2020_restor_n_avoidboth[]==1 & RDPF2020[]==1, 0, uSF2020_restor_n_avoidboth[])

uSF2020_restor_n_avoidboth[] <- ifelse(uSF2010[]==1 & DSF2020[]==1, 1, uSF2020_restor_n_avoidboth[])
uSF2020_restor_n_avoidboth[] <- ifelse(uSF2020[]==1 & UPF2010[]==1, 0, uSF2020_restor_n_avoidboth[])
uSF2020_restor_n_avoidboth[] <- ifelse(uSF2020[]==1 & uDPF2010[]==1, 0, uSF2020_restor_n_avoidboth[])
uSF2020_restor_n_avoidboth[] <- ifelse(uSF2020[]==1 & RDPF2010[]==1, 0, uSF2020_restor_n_avoidboth[])
uSF2020_restor_n_avoidboth[] <- ifelse(DSF2020[]==0 & uSF2020[]==0 & uSF2010[]==1, 1, uSF2020_restor_n_avoidboth[])
#plot(uSF2020_restor_n_avoidboth)
writeRaster(uSF2020_restor_n_avoidboth, "rasters/STM/input/LULC/uSF2020_restor_n_avoidboth.tif", format="GTiff", overwrite=T)
#uSF2020_restor_n_avoidboth <- raster("rasters/STM/input/LULC/uSF2020_restor_n_avoidboth.tif")

DSF2020_restor_n_avoidboth <- DSF2020
DSF2020_restor_n_avoidboth[] <- ifelse(DSF2020[]==1 & uSF2010[]==1, 0, DSF2020_restor_n_avoidboth[])
DSF2020_restor_n_avoidboth[] <- ifelse(DSF2020[]==1 & UPF2010[]==1, 0, DSF2020_restor_n_avoidboth[])
DSF2020_restor_n_avoidboth[] <- ifelse(DSF2020[]==1 & uDPF2010[]==1, 0, DSF2020_restor_n_avoidboth[])
DSF2020_restor_n_avoidboth[] <- ifelse(DSF2020[]==1 & RDPF2010[]==1, 0, DSF2020_restor_n_avoidboth[])
DSF2020_restor_n_avoidboth[] <- ifelse(uSF2020[]==0 & DSF2020[]==0 & DSF2010[]==1, 1, DSF2020_restor_n_avoidboth[])
#plot(DSF2020_restor_n_avoidboth)
writeRaster(DSF2020_restor_n_avoidboth, "rasters/STM/input/LULC/DSF2020_restor_n_avoidboth.tif", format="GTiff", overwrite=T)
#DSF2020_restor_n_avoidboth <- raster("rasters/STM/input/LULC/DSF2020_restor_n_avoidboth.tif")


###Undegraded primary forest
UPF2020_restor_n_avoidboth <- UPF2020
UPF2020_restor_n_avoidboth[] <- ifelse(UPF2010[]==1 & uDPF2020[]==1, 1, UPF2020_restor_n_avoidboth[])
UPF2020_restor_n_avoidboth[] <- ifelse(UPF2010[]==1 & RDPF2020[]==1, 1, UPF2020_restor_n_avoidboth[])
UPF2020_restor_n_avoidboth[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, UPF2020_restor_n_avoidboth[])
UPF2020_restor_n_avoidboth[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, UPF2020_restor_n_avoidboth[])
UPF2020_restor_n_avoidboth[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, UPF2020_restor_n_avoidboth[])
#plot(UPF2020_restor_n_avoidboth)
writeRaster(UPF2020_restor_n_avoidboth, "rasters/STM/input/LULC/UPF2020_restor_n_avoidboth.tif", format="GTiff", overwrite=T)
#UPF2020_restor_n_avoidboth <- raster("rasters/STM/input/LULC/UPF2020_restor_n_avoidboth.tif")


###Degraded primary forest
uDPF2020_restor_n_avoidboth <- uDPF2020
uDPF2020_restor_n_avoidboth[] <- ifelse(uDPF2020[]==1 & UPF2010[]==1, 0, uDPF2020_restor_n_avoidboth[])
uDPF2020_restor_n_avoidboth[] <- ifelse(uDPF2010[]==1 & RDPF2020[]==1, 1, uDPF2020_restor_n_avoidboth[])
uDPF2020_restor_n_avoidboth[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, uDPF2020_restor_n_avoidboth[])
uDPF2020_restor_n_avoidboth[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, uDPF2020_restor_n_avoidboth[])
uDPF2020_restor_n_avoidboth[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, uDPF2020_restor_n_avoidboth[])
#plot(uDPF2020_restor_n_avoidboth)
writeRaster(uDPF2020_restor_n_avoidboth, "rasters/STM/input/LULC/uDPF2020_restor_n_avoidboth.tif", format="GTiff", overwrite=T)
#uDPF2020_restor_n_avoidboth <- raster("rasters/STM/input/LULC/uDPF2020_restor_n_avoidboth.tif")

RDPF2020_restor_n_avoidboth <- RDPF2020
RDPF2020_restor_n_avoidboth[] <- ifelse(RDPF2020[]==1 & UPF2010[]==1, 0, RDPF2020_restor_n_avoidboth[])
RDPF2020_restor_n_avoidboth[] <- ifelse(RDPF2020[]==1 & uDPF2010[]==1, 0, RDPF2020_restor_n_avoidboth[])
RDPF2020_restor_n_avoidboth[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, RDPF2020_restor_n_avoidboth[])
RDPF2020_restor_n_avoidboth[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, RDPF2020_restor_n_avoidboth[])
RDPF2020_restor_n_avoidboth[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, RDPF2020_restor_n_avoidboth[])
#plot(RDPF2020_restor_n_avoidboth)
writeRaster(RDPF2020_restor_n_avoidboth, "rasters/STM/input/LULC/RDPF2020_restor_n_avoidboth.tif", format="GTiff", overwrite=T)
#RDPF2020_restor_n_avoidboth <- raster("rasters/STM/input/LULC/RDPF2020_restor_n_avoidboth.tif")




##2020 restoration and avoid both (Primary forest only)
###Secondary forest
uSF2020_restor_n_avoidboth2 <- sum(uSF2020, candidate.areas.final)
uSF2020_restor_n_avoidboth2[] <-  ifelse(uSF2010[]==0 & uSF2020[]==1, 1, uSF2020_restor_n_avoidboth2[])
uSF2020_restor_n_avoidboth2[] <-  ifelse(uSF2020_restor_n_avoidboth2[]==1 & DSF2020[]==1, 0, uSF2020_restor_n_avoidboth2[])
uSF2020_restor_n_avoidboth2[] <-  ifelse(uSF2020_restor_n_avoidboth2[]==1 & UPF2020[]==1, 0, uSF2020_restor_n_avoidboth2[])
uSF2020_restor_n_avoidboth2[] <-  ifelse(uSF2020_restor_n_avoidboth2[]==1 & uDPF2020[]==1, 0, uSF2020_restor_n_avoidboth2[])
uSF2020_restor_n_avoidboth2[] <-  ifelse(uSF2020_restor_n_avoidboth2[]==1 & RDPF2020[]==1, 0, uSF2020_restor_n_avoidboth2[])

uSF2020_restor_n_avoidboth2[] <- ifelse(uSF2020[]==1 & UPF2010[]==1, 0, uSF2020_restor_n_avoidboth2[])
uSF2020_restor_n_avoidboth2[] <- ifelse(uSF2020[]==1 & uDPF2010[]==1, 0, uSF2020_restor_n_avoidboth2[])
uSF2020_restor_n_avoidboth2[] <- ifelse(uSF2020[]==1 & RDPF2010[]==1, 0, uSF2020_restor_n_avoidboth2[])
#plot(uSF2020_restor_n_avoidboth2)
writeRaster(uSF2020_restor_n_avoidboth2, "rasters/STM/input/LULC/uSF2020_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)
#uSF2020_restor_n_avoidboth2 <- raster("rasters/STM/input/LULC/uSF2020_restor_n_avoidboth2.tif")

DSF2020_restor_n_avoidboth2 <- DSF2020
DSF2020_restor_n_avoidboth2[] <- ifelse(DSF2020[]==1 & UPF2010[]==1, 0, DSF2020_restor_n_avoidboth2[])
DSF2020_restor_n_avoidboth2[] <- ifelse(DSF2020[]==1 & uDPF2010[]==1, 0, DSF2020_restor_n_avoidboth2[])
DSF2020_restor_n_avoidboth2[] <- ifelse(DSF2020[]==1 & RDPF2010[]==1, 0, DSF2020_restor_n_avoidboth2[])
#plot(DSF2020_restor_n_avoidboth2)
writeRaster(DSF2020_restor_n_avoidboth2, "rasters/STM/input/LULC/DSF2020_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)
#DSF2020_restor_n_avoidboth2 <- raster("rasters/STM/input/LULC/DSF2020_restor_n_avoidboth2.tif")


###Undegraded primary forest
UPF2020_restor_n_avoidboth2 <- UPF2020
UPF2020_restor_n_avoidboth2[] <- ifelse(UPF2010[]==1 & uDPF2020[]==1, 1, UPF2020_restor_n_avoidboth2[])
UPF2020_restor_n_avoidboth2[] <- ifelse(UPF2010[]==1 & RDPF2020[]==1, 1, UPF2020_restor_n_avoidboth2[])
UPF2020_restor_n_avoidboth2[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, UPF2020_restor_n_avoidboth2[])
UPF2020_restor_n_avoidboth2[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, UPF2020_restor_n_avoidboth2[])
UPF2020_restor_n_avoidboth2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, UPF2020_restor_n_avoidboth2[])
#plot(UPF2020_restor_n_avoidboth2)
writeRaster(UPF2020_restor_n_avoidboth2, "rasters/STM/input/LULC/UPF2020_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)
#UPF2020_restor_n_avoidboth2 <- raster("rasters/STM/input/LULC/UPF2020_restor_n_avoidboth2.tif")


###Degraded primary forest
uDPF2020_restor_n_avoidboth2 <- uDPF2020
uDPF2020_restor_n_avoidboth2[] <- ifelse(uDPF2020[]==1 & UPF2010[]==1, 0, uDPF2020_restor_n_avoidboth2[])
uDPF2020_restor_n_avoidboth2[] <- ifelse(uDPF2010[]==1 & RDPF2020[]==1, 1, uDPF2020_restor_n_avoidboth2[])
uDPF2020_restor_n_avoidboth2[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, uDPF2020_restor_n_avoidboth2[])
uDPF2020_restor_n_avoidboth2[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, uDPF2020_restor_n_avoidboth2[])
uDPF2020_restor_n_avoidboth2[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, uDPF2020_restor_n_avoidboth2[])
#plot(uDPF2020_restor_n_avoidboth2)
writeRaster(uDPF2020_restor_n_avoidboth2, "rasters/STM/input/LULC/uDPF2020_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)
#uDPF2020_restor_n_avoidboth2 <- raster("rasters/STM/input/LULC/uDPF2020_restor_n_avoidboth2.tif")

RDPF2020_restor_n_avoidboth2 <- RDPF2020
RDPF2020_restor_n_avoidboth2[] <- ifelse(RDPF2020[]==1 & UPF2010[]==1, 0, RDPF2020_restor_n_avoidboth2[])
RDPF2020_restor_n_avoidboth2[] <- ifelse(RDPF2020[]==1 & uDPF2010[]==1, 0, RDPF2020_restor_n_avoidboth2[])
RDPF2020_restor_n_avoidboth2[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, RDPF2020_restor_n_avoidboth2[])
RDPF2020_restor_n_avoidboth2[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, RDPF2020_restor_n_avoidboth2[])
RDPF2020_restor_n_avoidboth2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, RDPF2020_restor_n_avoidboth2[])
#plot(RDPF2020_restor_n_avoidboth2)
writeRaster(RDPF2020_restor_n_avoidboth2, "rasters/STM/input/LULC/RDPF2020_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)
#RDPF2020_restor_n_avoidboth2 <- raster("rasters/STM/input/LULC/RDPF2020_restor_n_avoidboth2.tif")
















###checking
#undegraded primary forest == 1
UPF2020_restor_n_avoidboth2.sk <- UPF2020_restor_n_avoidboth2
#UPF2020_restor_n_avoidboth2.sk <- raster("rasters/STM/input/UPF2020_restor_n_avoidboth2.tif")
#UPF2020_restor_n_avoidboth2.sk[UPF2020_restor_n_avoidboth2.sk==1]<-1
UPF2020_restor_n_avoidboth2.sk <- mask(UPF2020_restor_n_avoidboth2.sk, stm.shp)
UPF2020_restor_n_avoidboth2.sk[UPF2020_restor_n_avoidboth2.sk[]==0] <- 333
#plot(UPF2020_restor_n_avoidboth2.sk, main="undegraded primary forest", legend=F)

#degraded primary forest == 10
uDPF2020_restor_n_avoidboth2.sk <- uDPF2020_restor_n_avoidboth2
#DPF2020_restor_n_avoidboth2.sk <- raster("rasters/STM/input/uDPF2020_restor_n_avoidboth2.tif")
uDPF2020_restor_n_avoidboth2.sk[uDPF2020_restor_n_avoidboth2.sk==1]<-10
uDPF2020_restor_n_avoidboth2.sk <- mask(uDPF2020_restor_n_avoidboth2.sk, stm.shp)
uDPF2020_restor_n_avoidboth2.sk[uDPF2020_restor_n_avoidboth2.sk[]==0] <- 333
#plot(uDPF2020_restor_n_avoidboth2.sk, main="degraded primary forest", legend=F)

#repeated degraded primary forest == 25
RDPF2020_restor_n_avoidboth2.sk <- RDPF2020_restor_n_avoidboth2
#RDPF2020_restor_n_avoidboth2.sk <- raster("rasters/STM/input/RDPF2020_restor_n_avoidboth2.tif")
RDPF2020_restor_n_avoidboth2.sk[RDPF2020_restor_n_avoidboth2.sk==1]<-25
RDPF2020_restor_n_avoidboth2.sk <- mask(RDPF2020_restor_n_avoidboth2.sk, stm.shp)
RDPF2020_restor_n_avoidboth2.sk[RDPF2020_restor_n_avoidboth2.sk[]==0] <- 333
#plot(RDPF2020_restor_n_avoidboth2.sk, main="degraded primary forest", legend=F)

#secondary forest == 100
uSF2020_restor_n_avoidboth2.sk <- uSF2020_restor_n_avoidboth2
#uSF2020_restor_n_avoidboth2.sk <- raster("rasters/STM/input/uSF2020_restor_n_avoidboth2.tif")
uSF2020_restor_n_avoidboth2.sk[uSF2020_restor_n_avoidboth2.sk==1]<-100
uSF2020_restor_n_avoidboth2.sk <- mask(uSF2020_restor_n_avoidboth2.sk, stm.shp)
uSF2020_restor_n_avoidboth2.sk[uSF2020_restor_n_avoidboth2.sk[]==0] <- 333
#plot(uSF2020_restor_n_avoidboth2.sk, main="secondary forest", legend=F)

#degraded secondary forest == 125
DSF2020_restor_n_avoidboth2.sk <- DSF2020_restor_n_avoidboth2
#DSF2020_restor_n_avoidboth2.sk <- raster("rasters/STM/input/DSF2020_restor_n_avoidboth2.tif")
DSF2020_restor_n_avoidboth2.sk[DSF2020_restor_n_avoidboth2.sk==1]<-125
DSF2020_restor_n_avoidboth2.sk <- mask(DSF2020_restor_n_avoidboth2.sk, stm.shp)
DSF2020_restor_n_avoidboth2.sk[DSF2020_restor_n_avoidboth2.sk[]==0] <- 333
#plot(DSF2020_restor_n_avoidboth2.sk, main="secondary forest", legend=F)

LULC2020_restor_n_avoidboth2 <- sum(UPF2020_restor_n_avoidboth2.sk, uSF2020_restor_n_avoidboth2.sk, na.rm = T)
LULC2020_restor_n_avoidboth2 <- sum(LULC2020_restor_n_avoidboth2, DSF2020_restor_n_avoidboth2.sk, na.rm = T)
LULC2020_restor_n_avoidboth2 <- sum(LULC2020_restor_n_avoidboth2, uDPF2020_restor_n_avoidboth2.sk, na.rm = T)
LULC2020_restor_n_avoidboth2 <- sum(LULC2020_restor_n_avoidboth2, RDPF2020_restor_n_avoidboth2.sk, na.rm = T)
LULC2020_restor_n_avoidboth2 <- mask(LULC2020_restor_n_avoidboth2, stm.shp)
#sort(unique(LULC2020_restor_n_avoidboth2[]))
LULC2020_restor_n_avoidboth2[LULC2020_restor_n_avoidboth2==1333]<-1
LULC2020_restor_n_avoidboth2[LULC2020_restor_n_avoidboth2==1342]<-10
LULC2020_restor_n_avoidboth2[LULC2020_restor_n_avoidboth2==1357]<-25
LULC2020_restor_n_avoidboth2[LULC2020_restor_n_avoidboth2==1432]<-100
LULC2020_restor_n_avoidboth2[LULC2020_restor_n_avoidboth2==1457]<-125
LULC2020_restor_n_avoidboth2[LULC2020_restor_n_avoidboth2==1665 ]<-0
#plot(LULC2020_restor_n_avoidboth2, main="forest cover 2020", legend=F)





# Create a data frame with the transition data
data_df2 <- data.frame(
  Period1 = data_df$Period1,
  #Period2 = factor(LULC2020[], levels = c(1, 10, 25, 125, 100, 0), labels = c("UPF", "DPF", "RDPF", "DSF", "SF", "D"))
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
  Period2 = factor(LULC2020_restor_n_avoidboth2[], levels = c(1, 10, 25, 125, 100, 0), 
                   labels = c("UPF", "DPF", "RDPF", "DSF", "SF", "D"))
)


#change the axis2 to check lulc transition in each scenario
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


rm(list=ls()[ls() %in% grep("\\b.sk\\b", ls(), value = T)])
rm(list=ls()[ls() %in% grep("LULC", ls(), value = T)])
gc()


#
#




# setting variables ============================================================
## creating directories
dir.create("rasters/STM/2010_real", recursive = T)
dir.create("rasters/STM/2020_real", recursive = T)
dir.create("rasters/STM/2020_avoiddeforest", recursive = T)
dir.create("rasters/STM/2020_avoiddeforest2", recursive = T) #PF only
dir.create("rasters/STM/2020_avoiddegrad", recursive = T)
dir.create("rasters/STM/2020_avoiddegrad2", recursive = T) #PF only
dir.create("rasters/STM/2020_avoidboth", recursive = T)
dir.create("rasters/STM/2020_avoidboth2", recursive = T) #PF only
dir.create("rasters/STM/2020_restor_wo_avoid", recursive = T)
dir.create("rasters/STM/2020_restor_n_avoiddeforest", recursive = T)
dir.create("rasters/STM/2020_restor_n_avoiddeforest2", recursive = T) #PF only
dir.create("rasters/STM/2020_restor_n_avoidboth", recursive = T)
dir.create("rasters/STM/2020_restor_n_avoidboth2", recursive = T) #PF only
#dir.create("models.output/costs", recursive = T)



## [meantemp] annual average temperature from nasa earth observation
#' download and save global data
# urls <- c("https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755469&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-01"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755470&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-02"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755471&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-03"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755472&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-04"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755473&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-05"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755474&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-06"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755475&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-07"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755476&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-08"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755477&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-09"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755478&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-10"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755479&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-11"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755480&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-12"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1784090&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-01"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1785058&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-02"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1785890&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-03"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1786979&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-04"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1794500&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-05"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1795357&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-06"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1796358&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-07"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1797155&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-08"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1799174&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-09"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1799941&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-10"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1800680&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-11"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1801650&cs=rgb&format=TIFF&width=3600&height=1800") #"2020-12"
# 
# dir.create("rasters/STM/input/climate")


### scenario 2010
temp.list <- list.files("rasters/STM/raw/climate", "LSTD", full.names = T, recursive = T)

temp2010.list <- grep("2010", temp.list, value = T)
temp2010 <- stack(temp2010.list)
#plot(temp2010)

meantemp2010 <- mean(temp2010, na.rm=T)
stm.meantemp2010 <- crop(meantemp2010, extent(stm.lulc[[1]]))
#plot(stm.meantemp2010)
#plot(stm.shp, add=T)

stm.meantemp2010 <- resample(stm.meantemp2010, stm.lulc[[1]], method='bilinear')
stm.meantemp2010[] <- ifelse(stm.lulc[[1]][]==0, NA, stm.meantemp2010[])
#plot(stm.meantemp2010)

#saving
writeRaster(stm.meantemp2010, "rasters/STM/2010_real/meantemps.tif", format="GTiff", overwrite=T)
#

### scenario 2020
temp2020.list <- grep("2020", temp.list, value = T)
temp2020 <- stack(temp2020.list)
#plot(temp2020)

meantemp2020 <- mean(temp2020, na.rm=T)
stm.meantemp2020 <- crop(meantemp2020, extent(stm.lulc[[2]]))
#plot(stm.meantemp2020)
#plot(stm.shp, add=T)

stm.meantemp2020 <- resample(stm.meantemp2020, stm.lulc[[2]], method='bilinear')
stm.meantemp2020[] <- ifelse(stm.lulc[[2]][]==0, NA, stm.meantemp2020[])
#plot(stm.meantemp2020)

#saving
writeRaster(stm.meantemp2020, "rasters/STM/2020_real/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(stm.meantemp2020, "rasters/STM/2020_avoiddeforest/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(stm.meantemp2020, "rasters/STM/2020_avoiddeforest2/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(stm.meantemp2020, "rasters/STM/2020_avoiddegrad/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(stm.meantemp2020, "rasters/STM/2020_avoiddegrad2/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(stm.meantemp2020, "rasters/STM/2020_avoidboth/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(stm.meantemp2020, "rasters/STM/2020_avoidboth2/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(stm.meantemp2020, "rasters/STM/2020_restor_wo_avoid/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(stm.meantemp2020, "rasters/STM/2020_restor_n_avoiddeforest/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(stm.meantemp2020, "rasters/STM/2020_restor_n_avoiddeforest2/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(stm.meantemp2020, "rasters/STM/2020_restor_n_avoidboth/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(stm.meantemp2020, "rasters/STM/2020_restor_n_avoidboth2/meantemps.tif", format="GTiff", overwrite=T)
#

rm(list=ls()[ls() %in% c("temp.list", "temp2010.list", "meantemp2010", "temp2020.list", "meantemp2020")]) #keeping only raster stack
gc()



#
#




## [meanprecip] annual average precipitation from nasa earth observation
#' download and save global data
# urls <- c("https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843747&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-01"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843749&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-02"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843751&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-03"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843755&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-04"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843735&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-05"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843745&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-06"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843753&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-07"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843759&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-08"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843761&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-09"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843763&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-10"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843765&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-11"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843769&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-12"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843983&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-01"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843985&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-02"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843987&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-03"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843989&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-04"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843991&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-05"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843993&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-06"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843995&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-07"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843997&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-08"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843999&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-09"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1844001&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-10"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1844003&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-11"
#           "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1844005&cs=rgb&format=TIFF&width=3600&height=1800") #"2020-12"
# 


### scenario 2010
precip.list <- list.files("rasters/STM/raw/climate", "GPM", full.names = T, recursive = T)

precip2010.list <- grep("2010", precip.list, value = T)
precip2010 <- stack(precip2010.list)
#plot(precip2010)

meanprecip2010 <- mean(precip2010, na.rm=T)
stm.meanprecip2010 <- crop(meanprecip2010, extent(stm.lulc[[1]]))
#plot(stm.meanprecip2010)
#plot(stm.shp, add=T)

stm.meanprecip2010 <- resample(stm.meanprecip2010, stm.lulc[[1]], method='bilinear')
stm.meanprecip2010[] <- ifelse(stm.lulc[[1]][]==0, NA, stm.meanprecip2010[])
#plot(stm.meanprecip2010)

#saving
writeRaster(stm.meanprecip2010, "rasters/STM/2010_real/meanprecips.tif", format="GTiff", overwrite=T)
#

### scenario 2020
precip2020.list <- grep("2020", precip.list, value = T)
precip2020 <- stack(precip2020.list)
#plot(precip2020)


meanprecip2020 <- mean(precip2020, na.rm=T)
stm.meanprecip2020 <- crop(meanprecip2020, extent(stm.lulc[[2]]))
#plot(stm.meanprecip2020)
#plot(stm.shp, add=T)

stm.meanprecip2020 <- resample(stm.meanprecip2020, stm.lulc[[2]], method='bilinear')
stm.meanprecip2020[] <- ifelse(stm.lulc[[2]][]==0, NA, stm.meanprecip2020[])
#plot(stm.meanprecip2020)

#saving
writeRaster(stm.meanprecip2020, "rasters/STM/2020_real/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(stm.meanprecip2020, "rasters/STM/2020_avoiddeforest/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(stm.meanprecip2020, "rasters/STM/2020_avoiddeforest2/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(stm.meanprecip2020, "rasters/STM/2020_avoiddegrad/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(stm.meanprecip2020, "rasters/STM/2020_avoiddegrad2/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(stm.meanprecip2020, "rasters/STM/2020_avoidboth/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(stm.meanprecip2020, "rasters/STM/2020_avoidboth2/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(stm.meanprecip2020, "rasters/STM/2020_restor_wo_avoid/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(stm.meanprecip2020, "rasters/STM/2020_restor_n_avoiddeforest/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(stm.meanprecip2020, "rasters/STM/2020_restor_n_avoiddeforest2/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(stm.meanprecip2020, "rasters/STM/2020_restor_n_avoidboth/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(stm.meanprecip2020, "rasters/STM/2020_restor_n_avoidboth2/meanprecips.tif", format="GTiff", overwrite=T)
#

rm(list=ls()[ls() %in% c("precip.list", "precip2010.list", "meanprecip2010", "precip2020.list", "meanprecip2020")])
gc()



#
#




## [elevation]
elevation <- raster("rasters/STM/raw/stm-elevation.tif")
elevation <- projectRaster(elevation, crs = std.proj)
elevation <- resample(elevation, stm.lulc[[1]], method='bilinear')
elevation[] <- ifelse(stm.lulc[[1]][]==0, NA, elevation[])

#saving
writeRaster(elevation, "rasters/STM/2010_real/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/STM/2020_real/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/STM/2020_avoiddeforest/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/STM/2020_avoiddeforest2/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/STM/2020_avoiddegrad/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/STM/2020_avoiddegrad2/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/STM/2020_avoidboth/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/STM/2020_avoidboth2/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/STM/2020_restor_wo_avoid/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/STM/2020_restor_n_avoiddeforest/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/STM/2020_restor_n_avoiddeforest2/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/STM/2020_restor_n_avoidboth/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/STM/2020_restor_n_avoidboth2/elevation.tif", format="GTiff", overwrite=T)



#
#




## [distroad] distance to roads
dist.road <- raster("rasters/STM/raw/stm-distance-to-road.tif")
dist.road <- projectRaster(dist.road, crs = std.proj)
dist.road <- resample(dist.road, stm.lulc[[1]], method='bilinear')
dist.road[] <- ifelse(stm.lulc[[1]][]==0, NA, dist.road[])

#saving
writeRaster(dist.road, "rasters/STM/2010_real/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/STM/2020_real/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/STM/2020_avoiddeforest/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/STM/2020_avoiddeforest2/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/STM/2020_avoiddegrad/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/STM/2020_avoiddegrad2/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/STM/2020_avoidboth/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/STM/2020_avoidboth2/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/STM/2020_restor_wo_avoid/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/STM/2020_restor_n_avoiddeforest/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/STM/2020_restor_n_avoiddeforest2/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/STM/2020_restor_n_avoidboth/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/STM/2020_restor_n_avoidboth2/distroad.tif", format="GTiff", overwrite=T)



#
#




## Rivers
#' source: HydroSHEDS - https://www.hydrosheds.org/hydroatlas
stm.river <- readOGR(dsn = "rasters/STM/raw", layer = "stm_RiverATLAS_v10")
#head(stm.river@data)

### selecting variables and creating river width attribute
stm.river@data <- stm.river@data %>% dplyr::select(HYRIV_ID:LENGTH_KM, ria_ha_csu, ria_ha_usu) %>% 
  mutate(
    ril_m = LENGTH_KM * 1000, #converting to meters
    ria_m2 = ria_ha_csu * 10000, #converting to square meters
    riw_m = (ria_m2 / ril_m) #estimating the mean river segment width
  )

#checking
#st_crs(stm.river)==st_crs(stm.shp)
#plot(stm.lulc[["stm.lulc.2010real"]])
#plot(stm.river, add=T)


### [distriver] distance to rivers
stm.river.raster <- rasterize(stm.river, stm.lulc[[1]], field = 1, background = 0)
stm.river.raster.mask <- stm.river.raster
stm.river.raster.mask[stm.river.raster.mask==0] <- NA
#checking
#inv.stm.river
#plot(stm.river.raster.mask, col="black")

stm.dist.river <- distance(stm.river.raster.mask, doEdge=T)
#writeRaster(stm.dist.river, "rasters/STM/raw/stm-distance-to-river.tif", format = "GTiff", overwrite = T)
#stm.dist.river <- raster("rasters/STM/raw/stm-distance-to-river.tif")

stm.dist.river[] <- ifelse(stm.lulc[[1]][]==0, NA, stm.dist.river[])
##cheking
#anyNA(stm.dist.river[])
#plot(stm.dist.river)

#saving
writeRaster(stm.dist.river, "rasters/STM/2010_real/distriver.tif", format="GTiff", overwrite=T)
writeRaster(stm.dist.river, "rasters/STM/2020_real/distriver.tif", format="GTiff", overwrite=T)
writeRaster(stm.dist.river, "rasters/STM/2020_avoiddeforest/distriver.tif", format="GTiff", overwrite=T)
writeRaster(stm.dist.river, "rasters/STM/2020_avoiddeforest2/distriver.tif", format="GTiff", overwrite=T)
writeRaster(stm.dist.river, "rasters/STM/2020_avoiddegrad/distriver.tif", format="GTiff", overwrite=T)
writeRaster(stm.dist.river, "rasters/STM/2020_avoiddegrad2/distriver.tif", format="GTiff", overwrite=T)
writeRaster(stm.dist.river, "rasters/STM/2020_avoidboth/distriver.tif", format="GTiff", overwrite=T)
writeRaster(stm.dist.river, "rasters/STM/2020_avoidboth2/distriver.tif", format="GTiff", overwrite=T)
writeRaster(stm.dist.river, "rasters/STM/2020_restor_wo_avoid/distriver.tif", format="GTiff", overwrite=T)
writeRaster(stm.dist.river, "rasters/STM/2020_restor_n_avoiddeforest/distriver.tif", format="GTiff", overwrite=T)
writeRaster(stm.dist.river, "rasters/STM/2020_restor_n_avoiddeforest2/distriver.tif", format="GTiff", overwrite=T)
writeRaster(stm.dist.river, "rasters/STM/2020_restor_n_avoidboth/distriver.tif", format="GTiff", overwrite=T)
writeRaster(stm.dist.river, "rasters/STM/2020_restor_n_avoidboth2/distriver.tif", format="GTiff", overwrite=T)
#

rm(list=ls()[ls() %in% c("stm.river.raster", "stm.river.raster.mask")])
gc()



#
#




## [distmarket] distance to municipality nucleus
stm.munic.nucleus <- data.frame(ID = "stm", long = -54.7009, lat = -2.45063)

stm.munic.nucleus.coord <- SpatialPointsDataFrame(coords = stm.munic.nucleus[,c("long","lat")], 
                                                  data = stm.munic.nucleus, 
                                                  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
stm.munic.nucleus.coord <- spTransform(stm.munic.nucleus.coord, crs(std.proj))


distmarket <- distanceFromPoints(object = stm.lulc[[1]], xy = stm.munic.nucleus.coord)
distmarket[] <- ifelse(stm.lulc[[1]][]==0, NA, distmarket[])

#saving
writeRaster(distmarket, "rasters/STM/2010_real/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/STM/2020_real/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/STM/2020_avoiddeforest/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/STM/2020_avoiddeforest2/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/STM/2020_avoiddegrad/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/STM/2020_avoiddegrad2/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/STM/2020_avoidboth/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/STM/2020_avoidboth2/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/STM/2020_restor_wo_avoid/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/STM/2020_restor_n_avoiddeforest/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/STM/2020_restor_n_avoiddeforest2/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/STM/2020_restor_n_avoidboth/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/STM/2020_restor_n_avoidboth2/distmarket.tif", format="GTiff", overwrite=T)
#

rm(list=ls()[ls() %in% c("stm.munic.nucleus", "stm.munic.nucleus.coord")])
gc()



#
#




## import rural properties shapefiles and data from SISCAR
#' https://www.car.gov.br/publico/municipios/downloads
#' see the "restoration_candidate_stm.R" on the script folder
#' for details on reducing overlap between properties and 
stm.car <- readOGR(dsn = "rasters/STM/raw", layer = "stm_car_after_restoration_2010")
names(stm.car@data) <- c("COD_IMOVEL","NUM_AREA","COD_ESTADO","NOM_MUNICI","NUM_MODULO","TIPO_IMOVE","SITUACAO","CONDICAO_I",
                         "num_area_ha","num_modulo_new","num_area_flag","FOREST_COVER_2007","FOREST_COVER_2007_PP","FOREST_COVER_2010","FOREST_COVER_2010_PP",
                         "NEED_INCREMENT","APP_FOREST_COVER_INCREMENT","APP_FOREST_COVER_INCREMENT_PP","NEED_INCREMENT_AFTER_APP",
                         "ARL_FOREST_COVER_INCREMENT","ARL_FOREST_COVER_INCREMENT_PP","NEED_INCREMENT_AFTER_ARL")
stm.car <- spTransform(stm.car, crs(std.proj))

#checking
#nrow(stm.car)
#nrow(stm.car@data[stm.car@data$NEED_INCREMENT==1,])
#head(stm.car@data[,"num_area_ha"])

#note:
#according to Brazilian Forest Code
#small properties have less than or equal to 4 fiscal modules
#medium/big properties heave more than 4 fiscal modules
#in STM, 1 fiscal module = 55ha [or 550000m2]


## [propertysize]
# transforming the properties polygons into rasters
# cell value is the property size
stm.car.raster <- rasterize(stm.car, stm.lulc[[1]], field = "num_area_ha", fun = "mean", background = 0)
stm.car.raster[] <- ifelse(stm.lulc[[1]][]==0, NA, stm.car.raster[])
#checking
#st_crs(stm.car.raster)==st_crs(stm.shp)
#plot(stm.car.raster)
#plot(stm.car, add=T)

#saving
writeRaster(stm.car.raster, "rasters/STM/2010_real/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(stm.car.raster, "rasters/STM/2020_real/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(stm.car.raster, "rasters/STM/2020_avoiddeforest/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(stm.car.raster, "rasters/STM/2020_avoiddeforest2/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(stm.car.raster, "rasters/STM/2020_avoiddegrad/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(stm.car.raster, "rasters/STM/2020_avoiddegrad2/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(stm.car.raster, "rasters/STM/2020_avoidboth/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(stm.car.raster, "rasters/STM/2020_avoidboth2/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(stm.car.raster, "rasters/STM/2020_restor_wo_avoid/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(stm.car.raster, "rasters/STM/2020_restor_n_avoiddeforest/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(stm.car.raster, "rasters/STM/2020_restor_n_avoiddeforest2/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(stm.car.raster, "rasters/STM/2020_restor_n_avoidboth/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(stm.car.raster, "rasters/STM/2020_restor_n_avoidboth2/propertysize.tif", format="GTiff", overwrite=T)



#
#



## Scenario: 2010 Real =========================================================
### Undegraded primary forest
#UPF2010 <- raster("rasters/STM/input/LULC/UPF2010_real.tif")

### mean upf cover in local scale (90m)
UPF2010.px <- focal(UPF2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2010.px)<-"UPFpx"
UPF2010.px[is.nan(UPF2010.px)] <- 0
UPF2010.px[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2010.px[])
#saving
writeRaster(UPF2010.px, "rasters/STM/2010_real/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2010.ls <- focal(UPF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2010.ls)<-"UPFls"
UPF2010.ls[is.nan(UPF2010.ls)] <- 0
UPF2010.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2010.ls[])
#saving
writeRaster(UPF2010.ls, "rasters/STM/2010_real/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2010.px", "UPF2010.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2010 <- raster("rasters/STM/input/LULC/uDPF2010_real.tif")
#RDPF2010 <- raster("rasters/STM/input/LULC/RDPF2010_real.tif")
DPF2010 <- sum(uDPF2010, RDPF2010)

### mean dpf cover in local scale (90m)
DPF2010.px <- focal(DPF2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2010.px)<-"DPFpx"
DPF2010.px[is.nan(DPF2010.px)] <- 0
DPF2010.px[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2010.px[])
#saving
writeRaster(DPF2010.px, "rasters/STM/2010_real/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2010.ls <- focal(DPF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2010.ls)<-"DPFls"
DPF2010.ls[is.nan(DPF2010.ls)] <- 0
DPF2010.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2010.ls[])
#saving
writeRaster(DPF2010.ls, "rasters/STM/2010_real/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2010.px", "DPF2010.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2010.pf <- stm.degrad.pf[["stm.degrad.2010real"]]
TSD2010.pf[is.na(TSD2010.pf)] <- 0

TSD2010.sf <- stm.degrad.sf[["stm.degradsf.2010real"]]
TSD2010.sf[is.na(TSD2010.sf)] <- 0

TSD2010 <- sum(TSD2010.pf, TSD2010.sf)
writeRaster(TSD2010, "rasters/STM/input/TSD2010.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2010.px <- focal(TSD2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2010.px)<-"TSDpx"
TSD2010.px[is.nan(TSD2010.px)] <- 0
TSD2010.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2010.px[])
#saving
writeRaster(TSD2010.px, "rasters/STM/2010_real/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2010.ls <- focal(TSD2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2010.ls)<-"TSDls"
TSD2010.ls[is.nan(TSD2010.ls)] <- 0
TSD2010.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2010.ls[])
#saving
writeRaster(TSD2010.ls, "rasters/STM/2010_real/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2010.px", "TSD2010.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2010 <- raster("rasters/STM/input/LULC/uSF2010_real.tif")
#DSF2010 <- raster("rasters/STM/input/LULC/DSF2010_real.tif")
SF2010 <- sum(uSF2010, DSF2010)

### mean sf cover in local scale (90m)
SF2010.px <- focal(SF2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2010.px)<-"SFpx"
SF2010.px[is.nan(SF2010.px)] <- 0
SF2010.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2010.px[])
#saving
writeRaster(SF2010.px, "rasters/STM/2010_real/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2010.ls <- focal(SF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2010.ls)<-"SFls"
SF2010.ls[is.nan(SF2010.ls)] <- 0
SF2010.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2010.ls[])
#saving
writeRaster(SF2010.ls, "rasters/STM/2010_real/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2010.px", "SF2010.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2010 <- stm.sfage[["stm.sfage.2010real"]]
writeRaster(SFAge2010, "rasters/STM/input/SFAge2010.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2010.px <- focal(SFAge2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2010.px)<-"SFAgepx"
SFAge2010.px[is.nan(SFAge2010.px)] <- 0
SFAge2010.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2010.px[])
#saving
writeRaster(SFAge2010.px, "rasters/STM/2010_real/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2010.ls <- focal(SFAge2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2010.ls)<-"SFAgels"
SFAge2010.ls[is.nan(SFAge2010.ls)] <- 0
SFAge2010.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2010.ls[])
#saving
writeRaster(SFAge2010.ls, "rasters/STM/2010_real/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2010.px", "SFAge2010.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2010young <- SFAge2010
SF2010young[] <- ifelse(SF2010young[]>2, 1, 0)
writeRaster(SF2010young, "rasters/STM/input/SF2010young.tif", format="GTiff", overwrite=T)

TF2010 <- sum(UPF2010, DPF2010)
TF2010 <- sum(TF2010, SF2010young)
writeRaster(TF2010, "rasters/STM/input/TF2010.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2010.px <- focal(TF2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2010.px)<-"TFpx"
TF2010.px[is.nan(TF2010.px)] <- 0
TF2010.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2010.px[])
#saving
writeRaster(TF2010.px, "rasters/STM/2010_real/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2010.ls <- focal(TF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2010.ls)<-"TFls"
TF2010.ls[is.nan(TF2010.ls)] <- 0
TF2010.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2010.ls[])
#saving
writeRaster(TF2010.ls, "rasters/STM/2010_real/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2010young", "SFAge2010.px", "SFAge2010.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2010mature <- SFAge2010
SF2010mature[] <- ifelse(SF2010mature[]>5, 1, 0)
writeRaster(SF2010mature, "rasters/STM/input/SF2010mature.tif", format="GTiff", overwrite=T)

MF2010 <- sum(UPF2010, DPF2010)
MF2010 <- sum(MF2010, SF2010mature)
writeRaster(MF2010, "rasters/STM/input/MF2010.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2010.px <- focal(MF2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2010.px)<-"MFpx"
MF2010.px[is.nan(MF2010.px)] <- 0
MF2010.px[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2010.px[])
#saving
writeRaster(MF2010.px, "rasters/STM/2010_real/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2010.ls <- focal(MF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2010.ls)<-"MFls"
MF2010.ls[is.nan(MF2010.ls)] <- 0
MF2010.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2010.ls[])
#saving
writeRaster(MF2010.ls, "rasters/STM/2010_real/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2010mature", "SFAge2010.px", "SFAge2010.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2010 <- raster("rasters/STM/input/MF2010.tif")
inv.MF2010 <- MF2010
inv.MF2010[inv.MF2010==1]<-NA
#cheking
#inv.MF2010
#plot(inv.MF2010)

edge.dist.2010 <- distance(inv.MF2010, doEdge=T)
names(edge.dist.2010)<-"edgedist"
edge.dist.2010[is.nan(edge.dist.2010)] <- 0
edge.dist.2010[] <- ifelse(stm.lulc[[1]][]==0, NA, edge.dist.2010[])
#saving
writeRaster(edge.dist.2010, "rasters/STM/2010_real/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2010")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2010 <- edge.dist.2010
edge2010[] <- ifelse(edge2010[] < 200, 0, ifelse(edge2010[]>300, 0, 1))
writeRaster(edge2010, "rasters/STM/input/edge2010_real.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2010.px <- focal(edge2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2010.px)<-"edgepx"
edge2010.px[is.nan(edge2010.px)] <- 0
edge2010.px[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2010.px[])
#saving
writeRaster(edge2010.px, "rasters/STM/2010_real/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2010.ls <- focal(edge2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2010.ls)<-"edgels"
edge2010.ls[is.nan(edge2010.ls)] <- 0
edge2010.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2010.ls[])
#saving
writeRaster(edge2010.ls, "rasters/STM/2010_real/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2010.px", "edge2010.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: 2020 Real =========================================================
### Undegraded primary forest
#UPF2020 <- raster("rasters/STM/input/LULC/UPF2020_real.tif")

### mean upf cover in local scale (90m)
UPF2020.px <- focal(UPF2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020.px)<-"UPFpx"
UPF2020.px[is.nan(UPF2020.px)] <- 0
UPF2020.px[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020.px[])
#saving
writeRaster(UPF2020.px, "rasters/STM/2020_real/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020.ls <- focal(UPF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020.ls)<-"UPFls"
UPF2020.ls[is.nan(UPF2020.ls)] <- 0
UPF2020.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020.ls[])
#saving
writeRaster(UPF2020.ls, "rasters/STM/2020_real/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020.px", "UPF2020.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020 <- raster("rasters/STM/input/LULC/uDPF2020_real.tif")
#RDPF2020 <- raster("rasters/STM/input/LULC/RDPF2020_real.tif")
DPF2020 <- sum(uDPF2020, RDPF2020)

### mean dpf cover in local scale (90m)
DPF2020.px <- focal(DPF2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020.px)<-"DPFpx"
DPF2020.px[is.nan(DPF2020.px)] <- 0
DPF2020.px[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020.px[])
#saving
writeRaster(DPF2020.px, "rasters/STM/2020_real/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020.ls <- focal(DPF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020.ls)<-"DPFls"
DPF2020.ls[is.nan(DPF2020.ls)] <- 0
DPF2020.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020.ls[])
#saving
writeRaster(DPF2020.ls, "rasters/STM/2020_real/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020.px", "DPF2020.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020.pf <- stm.degrad.pf[["stm.degrad.2020real"]]
TSD2020.pf[is.na(TSD2020.pf)] <- 0

TSD2020.sf <- stm.degrad.sf[["stm.degradsf.2020real"]]
TSD2020.sf[is.na(TSD2020.sf)] <- 0

TSD2020 <- sum(TSD2020.pf, TSD2020.sf)
writeRaster(TSD2020, "rasters/STM/input/TSD2020.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020.px <- focal(TSD2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020.px)<-"TSDpx"
TSD2020.px[is.nan(TSD2020.px)] <- 0
TSD2020.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020.px[])
#saving
writeRaster(TSD2020.px, "rasters/STM/2020_real/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020.ls <- focal(TSD2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020.ls)<-"TSDls"
TSD2020.ls[is.nan(TSD2020.ls)] <- 0
TSD2020.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020.ls[])
#saving
writeRaster(TSD2020.ls, "rasters/STM/2020_real/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020.px", "TSD2020.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020 <- raster("rasters/STM/input/LULC/uSF2020_real.tif")
#DSF2020 <- raster("rasters/STM/input/LULC/DSF2020_real.tif")
SF2020 <- sum(uSF2020, DSF2020)

### mean sf cover in local scale (90m)
SF2020.px <- focal(SF2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020.px)<-"SFpx"
SF2020.px[is.nan(SF2020.px)] <- 0
SF2020.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020.px[])
#saving
writeRaster(SF2020.px, "rasters/STM/2020_real/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020.ls <- focal(SF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020.ls)<-"SFls"
SF2020.ls[is.nan(SF2020.ls)] <- 0
SF2020.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020.ls[])
#saving
writeRaster(SF2020.ls, "rasters/STM/2020_real/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020.px", "SF2020.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020 <- stm.sfage[["stm.sfage.2020real"]]
writeRaster(SFAge2020, "rasters/STM/input/SFAge2020.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020.px <- focal(SFAge2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020.px)<-"SFAgepx"
SFAge2020.px[is.nan(SFAge2020.px)] <- 0
SFAge2020.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020.px[])
#saving
writeRaster(SFAge2020.px, "rasters/STM/2020_real/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020.ls <- focal(SFAge2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020.ls)<-"SFAgels"
SFAge2020.ls[is.nan(SFAge2020.ls)] <- 0
SFAge2020.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020.ls[])
#saving
writeRaster(SFAge2020.ls, "rasters/STM/2020_real/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020.px", "SFAge2020.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young <- SFAge2020
SF2020young[] <- ifelse(SF2020young[]>2, 1, 0)
writeRaster(SF2020young, "rasters/STM/input/SF2020young.tif", format="GTiff", overwrite=T)

TF2020 <- sum(UPF2020, DPF2020)
TF2020 <- sum(TF2020, SF2020young)
writeRaster(TF2020, "rasters/STM/input/TF2020.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020.px <- focal(TF2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020.px)<-"TFpx"
TF2020.px[is.nan(TF2020.px)] <- 0
TF2020.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020.px[])
#saving
writeRaster(TF2020.px, "rasters/STM/2020_real/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020.ls <- focal(TF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020.ls)<-"TFls"
TF2020.ls[is.nan(TF2020.ls)] <- 0
TF2020.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020.ls[])
#saving
writeRaster(TF2020.ls, "rasters/STM/2020_real/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020young", "SFAge2020.px", "SFAge2020.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature <- SFAge2020
SF2020mature[] <- ifelse(SF2020mature[]>5, 1, 0)
writeRaster(SF2020mature, "rasters/STM/input/SF2020mature.tif", format="GTiff", overwrite=T)

MF2020 <- sum(UPF2020, DPF2020)
MF2020 <- sum(MF2020, SF2020mature)
writeRaster(MF2020, "rasters/STM/input/MF2020.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020.px <- focal(MF2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020.px)<-"MFpx"
MF2020.px[is.nan(MF2020.px)] <- 0
MF2020.px[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020.px[])
#saving
writeRaster(MF2020.px, "rasters/STM/2020_real/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020.ls <- focal(MF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020.ls)<-"MFls"
MF2020.ls[is.nan(MF2020.ls)] <- 0
MF2020.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020.ls[])
#saving
writeRaster(MF2020.ls, "rasters/STM/2020_real/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020mature", "SFAge2020.px", "SFAge2020.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020 <- raster("rasters/STM/input/MF2020.tif")
inv.MF2020 <- MF2020
inv.MF2020[inv.MF2020==1]<-NA
#cheking
#inv.MF2020
#plot(inv.MF2020)

edge.dist.2020 <- distance(inv.MF2020, doEdge=T)
names(edge.dist.2020)<-"edgedist"
edge.dist.2020[is.nan(edge.dist.2020)] <- 0
edge.dist.2020[] <- ifelse(stm.lulc[[1]][]==0, NA, edge.dist.2020[])
#saving
writeRaster(edge.dist.2020, "rasters/STM/2020_real/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020 <- edge.dist.2020
edge2020[] <- ifelse(edge2020[] < 200, 0, ifelse(edge2020[]>300, 0, 1))
writeRaster(edge2020, "rasters/STM/input/edge2020_real.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020.px <- focal(edge2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020.px)<-"edgepx"
edge2020.px[is.nan(edge2020.px)] <- 0
edge2020.px[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020.px[])
#saving
writeRaster(edge2020.px, "rasters/STM/2020_real/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020.ls <- focal(edge2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020.ls)<-"edgels"
edge2020.ls[is.nan(edge2020.ls)] <- 0
edge2020.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020.ls[])
#saving
writeRaster(edge2020.ls, "rasters/STM/2020_real/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020.px", "edge2020.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Avoid degradation (all) ===========================================
### Undegraded primary forest
#UPF2020_avoiddegrad <- raster("rasters/STM/input/LULC/UPF2020_avoiddegrad.tif")

### mean upf cover in local scale (90m)
UPF2020_avoiddegrad.px <- focal(UPF2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_avoiddegrad.px)<-"UPFpx"
UPF2020_avoiddegrad.px[is.nan(UPF2020_avoiddegrad.px)] <- 0
UPF2020_avoiddegrad.px[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_avoiddegrad.px[])
#saving
writeRaster(UPF2020_avoiddegrad.px, "rasters/STM/2020_avoiddegrad/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_avoiddegrad.ls <- focal(UPF2020_avoiddegrad, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_avoiddegrad.ls)<-"UPFls"
UPF2020_avoiddegrad.ls[is.nan(UPF2020_avoiddegrad.ls)] <- 0
UPF2020_avoiddegrad.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_avoiddegrad.ls[])
#saving
writeRaster(UPF2020_avoiddegrad.ls, "rasters/STM/2020_avoiddegrad/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_avoiddegrad.px", "UPF2020_avoiddegrad.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_avoiddegrad <- raster("rasters/STM/input/LULC/uDPF2020_avoiddegrad.tif")
#RDPF2020_avoiddegrad <- raster("rasters/STM/input/LULC/RDPF2020_avoiddegrad.tif")
DPF2020_avoiddegrad <- sum(uDPF2020_avoiddegrad, RDPF2020_avoiddegrad)

### mean dpf cover in local scale (90m)
DPF2020_avoiddegrad.px <- focal(DPF2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_avoiddegrad.px)<-"DPFpx"
DPF2020_avoiddegrad.px[is.nan(DPF2020_avoiddegrad.px)] <- 0
DPF2020_avoiddegrad.px[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_avoiddegrad.px[])
#saving
writeRaster(DPF2020_avoiddegrad.px, "rasters/STM/2020_avoiddegrad/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_avoiddegrad.ls <- focal(DPF2020_avoiddegrad, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_avoiddegrad.ls)<-"DPFls"
DPF2020_avoiddegrad.ls[is.nan(DPF2020_avoiddegrad.ls)] <- 0
DPF2020_avoiddegrad.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_avoiddegrad.ls[])
#saving
writeRaster(DPF2020_avoiddegrad.ls, "rasters/STM/2020_avoiddegrad/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_avoiddegrad.px", "DPF2020_avoiddegrad.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_avoiddegrad.pf <- calc(stm.degrad.pf[["stm.degrad.2010real"]], fun=function(x){ifelse(is.na(x), x, ifelse(x==0, x, ifelse(x==300, x, x+10)))})
TSD2020_avoiddegrad.pf[is.na(TSD2020_avoiddegrad.pf)] <- 0

TSD2020_avoiddegrad.sf <- calc(stm.degrad.sf[["stm.degradsf.2010real"]], fun=function(x){ifelse(is.na(x), x, ifelse(x==0, x, ifelse(x==300, x, x+10)))})
TSD2020_avoiddegrad.sf[is.na(TSD2020_avoiddegrad.sf)] <- 0

TSD2020_avoiddegrad <- sum(TSD2020_avoiddegrad.pf, TSD2020_avoiddegrad.sf)
writeRaster(TSD2020_avoiddegrad, "rasters/STM/input/TSD2020_avoiddegrad.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_avoiddegrad.px <- focal(TSD2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_avoiddegrad.px)<-"TSDpx"
TSD2020_avoiddegrad.px[is.nan(TSD2020_avoiddegrad.px)] <- 0
TSD2020_avoiddegrad.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_avoiddegrad.px[])
#saving
writeRaster(TSD2020_avoiddegrad.px, "rasters/STM/2020_avoiddegrad/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_avoiddegrad.ls <- focal(TSD2020_avoiddegrad, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_avoiddegrad.ls)<-"TSDls"
TSD2020_avoiddegrad.ls[is.nan(TSD2020_avoiddegrad.ls)] <- 0
TSD2020_avoiddegrad.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_avoiddegrad.ls[])
#saving
writeRaster(TSD2020_avoiddegrad.ls, "rasters/STM/2020_avoiddegrad/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_avoiddegrad.px", "TSD2020_avoiddegrad.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_avoiddegrad <- raster("rasters/STM/input/LULC/uSF2020_avoiddegrad.tif")
#DSF2020_avoiddegrad <- raster("rasters/STM/input/LULC/DSF2020_avoiddegrad.tif")
SF2020_avoiddegrad <- sum(uSF2020_avoiddegrad, DSF2020_avoiddegrad)

### mean sf cover in local scale (90m)
SF2020_avoiddegrad.px <- focal(SF2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_avoiddegrad.px)<-"SFpx"
SF2020_avoiddegrad.px[is.nan(SF2020_avoiddegrad.px)] <- 0
SF2020_avoiddegrad.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_avoiddegrad.px[])
#saving
writeRaster(SF2020_avoiddegrad.px, "rasters/STM/2020_avoiddegrad/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_avoiddegrad.ls <- focal(SF2020_avoiddegrad, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_avoiddegrad.ls)<-"SFls"
SF2020_avoiddegrad.ls[is.nan(SF2020_avoiddegrad.ls)] <- 0
SF2020_avoiddegrad.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_avoiddegrad.ls[])
#saving
writeRaster(SF2020_avoiddegrad.ls, "rasters/STM/2020_avoiddegrad/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddegrad.px", "SF2020_avoiddegrad.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_avoiddegrad <- stm.sfage[["stm.sfage.2020real"]]
writeRaster(SFAge2020_avoiddegrad, "rasters/STM/input/SFAge2020_avoiddegrad.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_avoiddegrad.px <- focal(SFAge2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_avoiddegrad.px)<-"SFAgepx"
SFAge2020_avoiddegrad.px[is.nan(SFAge2020_avoiddegrad.px)] <- 0
SFAge2020_avoiddegrad.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_avoiddegrad.px[])
#saving
writeRaster(SFAge2020_avoiddegrad.px, "rasters/STM/2020_avoiddegrad/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_avoiddegrad.ls <- focal(SFAge2020_avoiddegrad, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_avoiddegrad.ls)<-"SFAgels"
SFAge2020_avoiddegrad.ls[is.nan(SFAge2020_avoiddegrad.ls)] <- 0
SFAge2020_avoiddegrad.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_avoiddegrad.ls[])
#saving
writeRaster(SFAge2020_avoiddegrad.ls, "rasters/STM/2020_avoiddegrad/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_avoiddegrad.px", "SFAge2020_avoiddegrad.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young_avoiddegrad <- SFAge2020_avoiddegrad
SF2020young_avoiddegrad[] <- ifelse(SF2020young_avoiddegrad[]>2, 1, 0)
writeRaster(SF2020young_avoiddegrad, "rasters/STM/input/SF2020young_avoiddegrad.tif", format="GTiff", overwrite=T)

TF2020_avoiddegrad <- sum(UPF2020_avoiddegrad, DPF2020_avoiddegrad)
TF2020_avoiddegrad <- sum(TF2020_avoiddegrad, SF2020young_avoiddegrad)
writeRaster(TF2020_avoiddegrad, "rasters/STM/input/TF2020_avoiddegrad.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_avoiddegrad.px <- focal(TF2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_avoiddegrad.px)<-"TFpx"
TF2020_avoiddegrad.px[is.nan(TF2020_avoiddegrad.px)] <- 0
TF2020_avoiddegrad.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_avoiddegrad.px[])
#saving
writeRaster(TF2020_avoiddegrad.px, "rasters/STM/2020_avoiddegrad/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_avoiddegrad.ls <- focal(TF2020_avoiddegrad, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_avoiddegrad.ls)<-"TFls"
TF2020_avoiddegrad.ls[is.nan(TF2020_avoiddegrad.ls)] <- 0
TF2020_avoiddegrad.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_avoiddegrad.ls[])
#saving
writeRaster(TF2020_avoiddegrad.ls, "rasters/STM/2020_avoiddegrad/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddegradyoung", "SFAge2020_avoiddegrad.px", "SFAge2020_avoiddegrad.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature_avoiddegrad <- SFAge2020_avoiddegrad
SF2020mature_avoiddegrad[] <- ifelse(SF2020mature_avoiddegrad[]>5, 1, 0)
writeRaster(SF2020mature_avoiddegrad, "rasters/STM/input/SF2020mature_avoiddegrad.tif", format="GTiff", overwrite=T)

MF2020_avoiddegrad <- sum(UPF2020_avoiddegrad, DPF2020_avoiddegrad)
MF2020_avoiddegrad <- sum(MF2020_avoiddegrad, SF2020mature_avoiddegrad)
writeRaster(MF2020_avoiddegrad, "rasters/STM/input/MF2020_avoiddegrad.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_avoiddegrad.px <- focal(MF2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_avoiddegrad.px)<-"MFpx"
MF2020_avoiddegrad.px[is.nan(MF2020_avoiddegrad.px)] <- 0
MF2020_avoiddegrad.px[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_avoiddegrad.px[])
#saving
writeRaster(MF2020_avoiddegrad.px, "rasters/STM/2020_avoiddegrad/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_avoiddegrad.ls <- focal(MF2020_avoiddegrad, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_avoiddegrad.ls)<-"MFls"
MF2020_avoiddegrad.ls[is.nan(MF2020_avoiddegrad.ls)] <- 0
MF2020_avoiddegrad.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_avoiddegrad.ls[])
#saving
writeRaster(MF2020_avoiddegrad.ls, "rasters/STM/2020_avoiddegrad/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddegradmature", "SFAge2020_avoiddegrad.px", "SFAge2020_avoiddegrad.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_avoiddegrad <- raster("rasters/STM/input/MF2020_avoiddegrad.tif")
inv.MF2020_avoiddegrad <- MF2020_avoiddegrad
inv.MF2020_avoiddegrad[inv.MF2020_avoiddegrad==1]<-NA
#cheking
#inv.MF2020_avoiddegrad
#plot(inv.MF2020_avoiddegrad)

edge.dist.2020_avoiddegrad <- distance(inv.MF2020_avoiddegrad, doEdge=T)
names(edge.dist.2020_avoiddegrad)<-"edgedist"
edge.dist.2020_avoiddegrad[is.nan(edge.dist.2020_avoiddegrad)] <- 0
edge.dist.2020_avoiddegrad[] <- ifelse(stm.lulc[[1]][]==0, NA, edge.dist.2020_avoiddegrad[])
#saving
writeRaster(edge.dist.2020_avoiddegrad, "rasters/STM/2020_avoiddegrad/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_avoiddegrad")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_avoiddegrad <- edge.dist.2020_avoiddegrad
edge2020_avoiddegrad[] <- ifelse(edge2020_avoiddegrad[] < 200, 0, ifelse(edge2020_avoiddegrad[]>300, 0, 1))
writeRaster(edge2020_avoiddegrad, "rasters/STM/input/edge2020_avoiddegrad.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_avoiddegrad.px <- focal(edge2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_avoiddegrad.px)<-"edgepx"
edge2020_avoiddegrad.px[is.nan(edge2020_avoiddegrad.px)] <- 0
edge2020_avoiddegrad.px[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_avoiddegrad.px[])
#saving
writeRaster(edge2020_avoiddegrad.px, "rasters/STM/2020_avoiddegrad/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_avoiddegrad.ls <- focal(edge2020_avoiddegrad, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_avoiddegrad.ls)<-"edgels"
edge2020_avoiddegrad.ls[is.nan(edge2020_avoiddegrad.ls)] <- 0
edge2020_avoiddegrad.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_avoiddegrad.ls[])
#saving
writeRaster(edge2020_avoiddegrad.ls, "rasters/STM/2020_avoiddegrad/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_avoiddegrad.px", "edge2020_avoiddegrad.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Avoid degradation (Primary forests only) ==========================
### Undegraded primary forest
#UPF2020_avoiddegrad2 <- raster("rasters/STM/input/LULC/UPF2020_avoiddegrad2.tif")

### mean upf cover in local scale (90m)
UPF2020_avoiddegrad2.px <- focal(UPF2020_avoiddegrad2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_avoiddegrad2.px)<-"UPFpx"
UPF2020_avoiddegrad2.px[is.nan(UPF2020_avoiddegrad2.px)] <- 0
UPF2020_avoiddegrad2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_avoiddegrad2.px[])
#saving
writeRaster(UPF2020_avoiddegrad2.px, "rasters/STM/2020_avoiddegrad2/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_avoiddegrad2.ls <- focal(UPF2020_avoiddegrad2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_avoiddegrad2.ls)<-"UPFls"
UPF2020_avoiddegrad2.ls[is.nan(UPF2020_avoiddegrad2.ls)] <- 0
UPF2020_avoiddegrad2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_avoiddegrad2.ls[])
#saving
writeRaster(UPF2020_avoiddegrad2.ls, "rasters/STM/2020_avoiddegrad2/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_avoiddegrad2.px", "UPF2020_avoiddegrad2.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_avoiddegrad2 <- raster("rasters/STM/input/LULC/uDPF2020_avoiddegrad2.tif")
#RDPF2020_avoiddegrad2 <- raster("rasters/STM/input/LULC/RDPF2020_avoiddegrad2.tif")
DPF2020_avoiddegrad2 <- sum(uDPF2020_avoiddegrad2, RDPF2020_avoiddegrad2)

### mean dpf cover in local scale (90m)
DPF2020_avoiddegrad2.px <- focal(DPF2020_avoiddegrad2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_avoiddegrad2.px)<-"DPFpx"
DPF2020_avoiddegrad2.px[is.nan(DPF2020_avoiddegrad2.px)] <- 0
DPF2020_avoiddegrad2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_avoiddegrad2.px[])
#saving
writeRaster(DPF2020_avoiddegrad2.px, "rasters/STM/2020_avoiddegrad2/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_avoiddegrad2.ls <- focal(DPF2020_avoiddegrad2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_avoiddegrad2.ls)<-"DPFls"
DPF2020_avoiddegrad2.ls[is.nan(DPF2020_avoiddegrad2.ls)] <- 0
DPF2020_avoiddegrad2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_avoiddegrad2.ls[])
#saving
writeRaster(DPF2020_avoiddegrad2.ls, "rasters/STM/2020_avoiddegrad2/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_avoiddegrad2.px", "DPF2020_avoiddegrad2.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_avoiddegrad2.pf <- calc(stm.degrad.pf[["stm.degrad.2010real"]], fun=function(x){ifelse(is.na(x), x, ifelse(x==0, x, ifelse(x==300, x, x+10)))})
TSD2020_avoiddegrad2.pf[is.na(TSD2020_avoiddegrad2.pf)] <- 0

TSD2020_avoiddegrad2.sf <- stm.degrad.sf[["stm.degradsf.2020real"]]
TSD2020_avoiddegrad2.sf[is.na(TSD2020_avoiddegrad2.sf)] <- 0

TSD2020_avoiddegrad2 <- sum(TSD2020_avoiddegrad2.pf, TSD2020_avoiddegrad2.sf)
TSD2020_avoiddegrad2[TSD2020_avoiddegrad2>300] <- 0
writeRaster(TSD2020_avoiddegrad2, "rasters/STM/input/TSD2020_avoiddegrad2.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_avoiddegrad2.px <- focal(TSD2020_avoiddegrad2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_avoiddegrad2.px)<-"TSDpx"
TSD2020_avoiddegrad2.px[is.nan(TSD2020_avoiddegrad2.px)] <- 0
TSD2020_avoiddegrad2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_avoiddegrad2.px[])
#saving
writeRaster(TSD2020_avoiddegrad2.px, "rasters/STM/2020_avoiddegrad2/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_avoiddegrad2.ls <- focal(TSD2020_avoiddegrad2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_avoiddegrad2.ls)<-"TSDls"
TSD2020_avoiddegrad2.ls[is.nan(TSD2020_avoiddegrad2.ls)] <- 0
TSD2020_avoiddegrad2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_avoiddegrad2.ls[])
#saving
writeRaster(TSD2020_avoiddegrad2.ls, "rasters/STM/2020_avoiddegrad2/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_avoiddegrad2.px", "TSD2020_avoiddegrad2.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_avoiddegrad2 <- raster("rasters/STM/input/LULC/uSF2020_avoiddegrad2.tif")
#DSF2020_avoiddegrad2 <- raster("rasters/STM/input/LULC/DSF2020_avoiddegrad2.tif")
SF2020_avoiddegrad2 <- sum(uSF2020_avoiddegrad2, DSF2020_avoiddegrad2)

### mean sf cover in local scale (90m)
SF2020_avoiddegrad2.px <- focal(SF2020_avoiddegrad2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_avoiddegrad2.px)<-"SFpx"
SF2020_avoiddegrad2.px[is.nan(SF2020_avoiddegrad2.px)] <- 0
SF2020_avoiddegrad2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_avoiddegrad2.px[])
#saving
writeRaster(SF2020_avoiddegrad2.px, "rasters/STM/2020_avoiddegrad2/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_avoiddegrad2.ls <- focal(SF2020_avoiddegrad2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_avoiddegrad2.ls)<-"SFls"
SF2020_avoiddegrad2.ls[is.nan(SF2020_avoiddegrad2.ls)] <- 0
SF2020_avoiddegrad2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_avoiddegrad2.ls[])
#saving
writeRaster(SF2020_avoiddegrad2.ls, "rasters/STM/2020_avoiddegrad2/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddegrad2.px", "SF2020_avoiddegrad2.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_avoiddegrad2 <- stm.sfage[["stm.sfage.2020real"]]
writeRaster(SFAge2020_avoiddegrad2, "rasters/STM/input/SFAge2020_avoiddegrad2.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_avoiddegrad2.px <- focal(SFAge2020_avoiddegrad2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_avoiddegrad2.px)<-"SFAgepx"
SFAge2020_avoiddegrad2.px[is.nan(SFAge2020_avoiddegrad2.px)] <- 0
SFAge2020_avoiddegrad2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_avoiddegrad2.px[])
#saving
writeRaster(SFAge2020_avoiddegrad2.px, "rasters/STM/2020_avoiddegrad2/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_avoiddegrad2.ls <- focal(SFAge2020_avoiddegrad2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_avoiddegrad2.ls)<-"SFAgels"
SFAge2020_avoiddegrad2.ls[is.nan(SFAge2020_avoiddegrad2.ls)] <- 0
SFAge2020_avoiddegrad2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_avoiddegrad2.ls[])
#saving
writeRaster(SFAge2020_avoiddegrad2.ls, "rasters/STM/2020_avoiddegrad2/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_avoiddegrad2.px", "SFAge2020_avoiddegrad2.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young_avoiddegrad2 <- SFAge2020_avoiddegrad2
SF2020young_avoiddegrad2[] <- ifelse(SF2020young_avoiddegrad2[]>2, 1, 0)
writeRaster(SF2020young_avoiddegrad2, "rasters/STM/input/SF2020young_avoiddegrad2.tif", format="GTiff", overwrite=T)

TF2020_avoiddegrad2 <- sum(UPF2020_avoiddegrad2, DPF2020_avoiddegrad2)
TF2020_avoiddegrad2 <- sum(TF2020_avoiddegrad2, SF2020young_avoiddegrad2)
writeRaster(TF2020_avoiddegrad2, "rasters/STM/input/TF2020_avoiddegrad2.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_avoiddegrad2.px <- focal(TF2020_avoiddegrad2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_avoiddegrad2.px)<-"TFpx"
TF2020_avoiddegrad2.px[is.nan(TF2020_avoiddegrad2.px)] <- 0
TF2020_avoiddegrad2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_avoiddegrad2.px[])
#saving
writeRaster(TF2020_avoiddegrad2.px, "rasters/STM/2020_avoiddegrad2/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_avoiddegrad2.ls <- focal(TF2020_avoiddegrad2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_avoiddegrad2.ls)<-"TFls"
TF2020_avoiddegrad2.ls[is.nan(TF2020_avoiddegrad2.ls)] <- 0
TF2020_avoiddegrad2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_avoiddegrad2.ls[])
#saving
writeRaster(TF2020_avoiddegrad2.ls, "rasters/STM/2020_avoiddegrad2/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddegrad2young", "SFAge2020_avoiddegrad2.px", "SFAge2020_avoiddegrad2.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature_avoiddegrad2 <- SFAge2020_avoiddegrad2
SF2020mature_avoiddegrad2[] <- ifelse(SF2020mature_avoiddegrad2[]>5, 1, 0)
writeRaster(SF2020mature_avoiddegrad2, "rasters/STM/input/SF2020mature_avoiddegrad2.tif", format="GTiff", overwrite=T)

MF2020_avoiddegrad2 <- sum(UPF2020_avoiddegrad2, DPF2020_avoiddegrad2)
MF2020_avoiddegrad2 <- sum(MF2020_avoiddegrad2, SF2020mature_avoiddegrad2)
writeRaster(MF2020_avoiddegrad2, "rasters/STM/input/MF2020_avoiddegrad2.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_avoiddegrad2.px <- focal(MF2020_avoiddegrad2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_avoiddegrad2.px)<-"MFpx"
MF2020_avoiddegrad2.px[is.nan(MF2020_avoiddegrad2.px)] <- 0
MF2020_avoiddegrad2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_avoiddegrad2.px[])
#saving
writeRaster(MF2020_avoiddegrad2.px, "rasters/STM/2020_avoiddegrad2/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_avoiddegrad2.ls <- focal(MF2020_avoiddegrad2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_avoiddegrad2.ls)<-"MFls"
MF2020_avoiddegrad2.ls[is.nan(MF2020_avoiddegrad2.ls)] <- 0
MF2020_avoiddegrad2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_avoiddegrad2.ls[])
#saving
writeRaster(MF2020_avoiddegrad2.ls, "rasters/STM/2020_avoiddegrad2/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddegrad2mature", "SFAge2020_avoiddegrad2.px", "SFAge2020_avoiddegrad2.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_avoiddegrad2 <- raster("rasters/STM/input/MF2020_avoiddegrad2.tif")
inv.MF2020_avoiddegrad2 <- MF2020_avoiddegrad2
inv.MF2020_avoiddegrad2[inv.MF2020_avoiddegrad2==1]<-NA
#cheking
#inv.MF2020_avoiddegrad2
#plot(inv.MF2020_avoiddegrad2)

edge.dist.2020_avoiddegrad2 <- distance(inv.MF2020_avoiddegrad2, doEdge=T)
names(edge.dist.2020_avoiddegrad2)<-"edgedist"
edge.dist.2020_avoiddegrad2[is.nan(edge.dist.2020_avoiddegrad2)] <- 0
edge.dist.2020_avoiddegrad2[] <- ifelse(stm.lulc[[1]][]==0, NA, edge.dist.2020_avoiddegrad2[])
#saving
writeRaster(edge.dist.2020_avoiddegrad2, "rasters/STM/2020_avoiddegrad2/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_avoiddegrad2")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_avoiddegrad2 <- edge.dist.2020_avoiddegrad2
edge2020_avoiddegrad2[] <- ifelse(edge2020_avoiddegrad2[] < 200, 0, ifelse(edge2020_avoiddegrad2[]>300, 0, 1))
writeRaster(edge2020_avoiddegrad2, "rasters/STM/input/edge2020_avoiddegrad2.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_avoiddegrad2.px <- focal(edge2020_avoiddegrad2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_avoiddegrad2.px)<-"edgepx"
edge2020_avoiddegrad2.px[is.nan(edge2020_avoiddegrad2.px)] <- 0
edge2020_avoiddegrad2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_avoiddegrad2.px[])
#saving
writeRaster(edge2020_avoiddegrad2.px, "rasters/STM/2020_avoiddegrad2/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_avoiddegrad2.ls <- focal(edge2020_avoiddegrad2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_avoiddegrad2.ls)<-"edgels"
edge2020_avoiddegrad2.ls[is.nan(edge2020_avoiddegrad2.ls)] <- 0
edge2020_avoiddegrad2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_avoiddegrad2.ls[])
#saving
writeRaster(edge2020_avoiddegrad2.ls, "rasters/STM/2020_avoiddegrad2/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_avoiddegrad2.px", "edge2020_avoiddegrad2.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Avoid deforestation (all) =========================================
### Undegraded primary forest
#UPF2020_avoiddeforest <- raster("rasters/STM/input/LULC/UPF2020_avoiddeforest.tif")

### mean upf cover in local scale (90m)
UPF2020_avoiddeforest.px <- focal(UPF2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_avoiddeforest.px)<-"UPFpx"
UPF2020_avoiddeforest.px[is.nan(UPF2020_avoiddeforest.px)] <- 0
UPF2020_avoiddeforest.px[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_avoiddeforest.px[])
#saving
writeRaster(UPF2020_avoiddeforest.px, "rasters/STM/2020_avoiddeforest/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_avoiddeforest.ls <- focal(UPF2020_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_avoiddeforest.ls)<-"UPFls"
UPF2020_avoiddeforest.ls[is.nan(UPF2020_avoiddeforest.ls)] <- 0
UPF2020_avoiddeforest.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_avoiddeforest.ls[])
#saving
writeRaster(UPF2020_avoiddeforest.ls, "rasters/STM/2020_avoiddeforest/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_avoiddeforest.px", "UPF2020_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_avoiddeforest <- raster("rasters/STM/input/LULC/uDPF2020_avoiddeforest.tif")
#RDPF2020_avoiddeforest <- raster("rasters/STM/input/LULC/RDPF2020_avoiddeforest.tif")
DPF2020_avoiddeforest <- sum(uDPF2020_avoiddeforest, RDPF2020_avoiddeforest)

### mean dpf cover in local scale (90m)
DPF2020_avoiddeforest.px <- focal(DPF2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_avoiddeforest.px)<-"DPFpx"
DPF2020_avoiddeforest.px[is.nan(DPF2020_avoiddeforest.px)] <- 0
DPF2020_avoiddeforest.px[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_avoiddeforest.px[])
#saving
writeRaster(DPF2020_avoiddeforest.px, "rasters/STM/2020_avoiddeforest/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_avoiddeforest.ls <- focal(DPF2020_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_avoiddeforest.ls)<-"DPFls"
DPF2020_avoiddeforest.ls[is.nan(DPF2020_avoiddeforest.ls)] <- 0
DPF2020_avoiddeforest.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_avoiddeforest.ls[])
#saving
writeRaster(DPF2020_avoiddeforest.ls, "rasters/STM/2020_avoiddeforest/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_avoiddeforest.px", "DPF2020_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_avoiddeforest <- TSD2020
writeRaster(TSD2020_avoiddeforest, "rasters/STM/input/TSD2020_avoiddeforest.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_avoiddeforest.px <- focal(TSD2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_avoiddeforest.px)<-"TSDpx"
TSD2020_avoiddeforest.px[is.nan(TSD2020_avoiddeforest.px)] <- 0
TSD2020_avoiddeforest.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_avoiddeforest.px[])
#saving
writeRaster(TSD2020_avoiddeforest.px, "rasters/STM/2020_avoiddeforest/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_avoiddeforest.ls <- focal(TSD2020_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_avoiddeforest.ls)<-"TSDls"
TSD2020_avoiddeforest.ls[is.nan(TSD2020_avoiddeforest.ls)] <- 0
TSD2020_avoiddeforest.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_avoiddeforest.ls[])
#saving
writeRaster(TSD2020_avoiddeforest.ls, "rasters/STM/2020_avoiddeforest/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_avoiddeforest.px", "TSD2020_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_avoiddeforest <- raster("rasters/STM/input/LULC/uSF2020_avoiddeforest.tif")
#DSF2020_avoiddeforest <- raster("rasters/STM/input/LULC/DSF2020_avoiddeforest.tif")
SF2020_avoiddeforest <- sum(uSF2020_avoiddeforest, DSF2020_avoiddeforest)

### mean sf cover in local scale (90m)
SF2020_avoiddeforest.px <- focal(SF2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_avoiddeforest.px)<-"SFpx"
SF2020_avoiddeforest.px[is.nan(SF2020_avoiddeforest.px)] <- 0
SF2020_avoiddeforest.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_avoiddeforest.px[])
#saving
writeRaster(SF2020_avoiddeforest.px, "rasters/STM/2020_avoiddeforest/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_avoiddeforest.ls <- focal(SF2020_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_avoiddeforest.ls)<-"SFls"
SF2020_avoiddeforest.ls[is.nan(SF2020_avoiddeforest.ls)] <- 0
SF2020_avoiddeforest.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_avoiddeforest.ls[])
#saving
writeRaster(SF2020_avoiddeforest.ls, "rasters/STM/2020_avoiddeforest/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddeforest.px", "SF2020_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_avoiddeforest <- calc(stm.sfage[["stm.sfage.2010real"]], fun=function(x){ifelse(is.na(x), x, ifelse(x==0, x, x+10))})
writeRaster(SFAge2020_avoiddeforest, "rasters/STM/input/SFAge2020_avoiddeforest.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_avoiddeforest.px <- focal(SFAge2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_avoiddeforest.px)<-"SFAgepx"
SFAge2020_avoiddeforest.px[is.nan(SFAge2020_avoiddeforest.px)] <- 0
SFAge2020_avoiddeforest.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_avoiddeforest.px[])
#saving
writeRaster(SFAge2020_avoiddeforest.px, "rasters/STM/2020_avoiddeforest/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_avoiddeforest.ls <- focal(SFAge2020_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_avoiddeforest.ls)<-"SFAgels"
SFAge2020_avoiddeforest.ls[is.nan(SFAge2020_avoiddeforest.ls)] <- 0
SFAge2020_avoiddeforest.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_avoiddeforest.ls[])
#saving
writeRaster(SFAge2020_avoiddeforest.ls, "rasters/STM/2020_avoiddeforest/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_avoiddeforest.px", "SFAge2020_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young_avoiddeforest <- SFAge2020_avoiddeforest
SF2020young_avoiddeforest[] <- ifelse(SF2020young_avoiddeforest[]>2, 1, 0)
writeRaster(SF2020young_avoiddeforest, "rasters/STM/input/SF2020young_avoiddeforest.tif", format="GTiff", overwrite=T)

TF2020_avoiddeforest <- sum(UPF2020_avoiddeforest, DPF2020_avoiddeforest)
TF2020_avoiddeforest <- sum(TF2020_avoiddeforest, SF2020young_avoiddeforest)
writeRaster(TF2020_avoiddeforest, "rasters/STM/input/TF2020_avoiddeforest.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_avoiddeforest.px <- focal(TF2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_avoiddeforest.px)<-"TFpx"
TF2020_avoiddeforest.px[is.nan(TF2020_avoiddeforest.px)] <- 0
TF2020_avoiddeforest.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_avoiddeforest.px[])
#saving
writeRaster(TF2020_avoiddeforest.px, "rasters/STM/2020_avoiddeforest/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_avoiddeforest.ls <- focal(TF2020_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_avoiddeforest.ls)<-"TFls"
TF2020_avoiddeforest.ls[is.nan(TF2020_avoiddeforest.ls)] <- 0
TF2020_avoiddeforest.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_avoiddeforest.ls[])
#saving
writeRaster(TF2020_avoiddeforest.ls, "rasters/STM/2020_avoiddeforest/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddeforestyoung", "SFAge2020_avoiddeforest.px", "SFAge2020_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature_avoiddeforest <- SFAge2020_avoiddeforest
SF2020mature_avoiddeforest[] <- ifelse(SF2020mature_avoiddeforest[]>5, 1, 0)
writeRaster(SF2020mature_avoiddeforest, "rasters/STM/input/SF2020mature_avoiddeforest.tif", format="GTiff", overwrite=T)

MF2020_avoiddeforest <- sum(UPF2020_avoiddeforest, DPF2020_avoiddeforest)
MF2020_avoiddeforest <- sum(MF2020_avoiddeforest, SF2020mature_avoiddeforest)
writeRaster(MF2020_avoiddeforest, "rasters/STM/input/MF2020_avoiddeforest.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_avoiddeforest.px <- focal(MF2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_avoiddeforest.px)<-"MFpx"
MF2020_avoiddeforest.px[is.nan(MF2020_avoiddeforest.px)] <- 0
MF2020_avoiddeforest.px[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_avoiddeforest.px[])
#saving
writeRaster(MF2020_avoiddeforest.px, "rasters/STM/2020_avoiddeforest/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_avoiddeforest.ls <- focal(MF2020_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_avoiddeforest.ls)<-"MFls"
MF2020_avoiddeforest.ls[is.nan(MF2020_avoiddeforest.ls)] <- 0
MF2020_avoiddeforest.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_avoiddeforest.ls[])
#saving
writeRaster(MF2020_avoiddeforest.ls, "rasters/STM/2020_avoiddeforest/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddeforestmature", "SFAge2020_avoiddeforest.px", "SFAge2020_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_avoiddeforest <- raster("rasters/STM/input/MF2020_avoiddeforest.tif")
inv.MF2020_avoiddeforest <- MF2020_avoiddeforest
inv.MF2020_avoiddeforest[inv.MF2020_avoiddeforest==1]<-NA
#cheking
#inv.MF2020_avoiddeforest
#plot(inv.MF2020_avoiddeforest)

edge.dist.2020_avoiddeforest <- distance(inv.MF2020_avoiddeforest, doEdge=T)
names(edge.dist.2020_avoiddeforest)<-"edgedist"
edge.dist.2020_avoiddeforest[is.nan(edge.dist.2020_avoiddeforest)] <- 0
edge.dist.2020_avoiddeforest[] <- ifelse(stm.lulc[[1]][]==0, NA, edge.dist.2020_avoiddeforest[])
#saving
writeRaster(edge.dist.2020_avoiddeforest, "rasters/STM/2020_avoiddeforest/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_avoiddeforest")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_avoiddeforest <- edge.dist.2020_avoiddeforest
edge2020_avoiddeforest[] <- ifelse(edge2020_avoiddeforest[] < 200, 0, ifelse(edge2020_avoiddeforest[]>300, 0, 1))
writeRaster(edge2020_avoiddeforest, "rasters/STM/input/edge2020_avoiddeforest.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_avoiddeforest.px <- focal(edge2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_avoiddeforest.px)<-"edgepx"
edge2020_avoiddeforest.px[is.nan(edge2020_avoiddeforest.px)] <- 0
edge2020_avoiddeforest.px[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_avoiddeforest.px[])
#saving
writeRaster(edge2020_avoiddeforest.px, "rasters/STM/2020_avoiddeforest/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_avoiddeforest.ls <- focal(edge2020_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_avoiddeforest.ls)<-"edgels"
edge2020_avoiddeforest.ls[is.nan(edge2020_avoiddeforest.ls)] <- 0
edge2020_avoiddeforest.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_avoiddeforest.ls[])
#saving
writeRaster(edge2020_avoiddeforest.ls, "rasters/STM/2020_avoiddeforest/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_avoiddeforest.px", "edge2020_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Avoid deforestation (Primary forests only) ========================
### Undegraded primary forest
#UPF2020_avoiddeforest2 <- raster("rasters/STM/input/LULC/UPF2020_avoiddeforest2.tif")

### mean upf cover in local scale (90m)
UPF2020_avoiddeforest2.px <- focal(UPF2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_avoiddeforest2.px)<-"UPFpx"
UPF2020_avoiddeforest2.px[is.nan(UPF2020_avoiddeforest2.px)] <- 0
UPF2020_avoiddeforest2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_avoiddeforest2.px[])
#saving
writeRaster(UPF2020_avoiddeforest2.px, "rasters/STM/2020_avoiddeforest2/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_avoiddeforest2.ls <- focal(UPF2020_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_avoiddeforest2.ls)<-"UPFls"
UPF2020_avoiddeforest2.ls[is.nan(UPF2020_avoiddeforest2.ls)] <- 0
UPF2020_avoiddeforest2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_avoiddeforest2.ls[])
#saving
writeRaster(UPF2020_avoiddeforest2.ls, "rasters/STM/2020_avoiddeforest2/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_avoiddeforest2.px", "UPF2020_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_avoiddeforest2 <- raster("rasters/STM/input/LULC/uDPF2020_avoiddeforest2.tif")
#RDPF2020_avoiddeforest2 <- raster("rasters/STM/input/LULC/RDPF2020_avoiddeforest2.tif")
DPF2020_avoiddeforest2 <- sum(uDPF2020_avoiddeforest2, RDPF2020_avoiddeforest2)

### mean dpf cover in local scale (90m)
DPF2020_avoiddeforest2.px <- focal(DPF2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_avoiddeforest2.px)<-"DPFpx"
DPF2020_avoiddeforest2.px[is.nan(DPF2020_avoiddeforest2.px)] <- 0
DPF2020_avoiddeforest2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_avoiddeforest2.px[])
#saving
writeRaster(DPF2020_avoiddeforest2.px, "rasters/STM/2020_avoiddeforest2/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_avoiddeforest2.ls <- focal(DPF2020_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_avoiddeforest2.ls)<-"DPFls"
DPF2020_avoiddeforest2.ls[is.nan(DPF2020_avoiddeforest2.ls)] <- 0
DPF2020_avoiddeforest2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_avoiddeforest2.ls[])
#saving
writeRaster(DPF2020_avoiddeforest2.ls, "rasters/STM/2020_avoiddeforest2/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_avoiddeforest2.px", "DPF2020_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_avoiddeforest2 <- TSD2020
writeRaster(TSD2020_avoiddeforest2, "rasters/STM/input/TSD2020_avoiddeforest2.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_avoiddeforest2.px <- focal(TSD2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_avoiddeforest2.px)<-"TSDpx"
TSD2020_avoiddeforest2.px[is.nan(TSD2020_avoiddeforest2.px)] <- 0
TSD2020_avoiddeforest2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_avoiddeforest2.px[])
#saving
writeRaster(TSD2020_avoiddeforest2.px, "rasters/STM/2020_avoiddeforest2/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_avoiddeforest2.ls <- focal(TSD2020_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_avoiddeforest2.ls)<-"TSDls"
TSD2020_avoiddeforest2.ls[is.nan(TSD2020_avoiddeforest2.ls)] <- 0
TSD2020_avoiddeforest2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_avoiddeforest2.ls[])
#saving
writeRaster(TSD2020_avoiddeforest2.ls, "rasters/STM/2020_avoiddeforest2/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_avoiddeforest2.px", "TSD2020_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_avoiddeforest2 <- raster("rasters/STM/input/LULC/uSF2020_avoiddeforest2.tif")
#DSF2020_avoiddeforest2 <- raster("rasters/STM/input/LULC/DSF2020_avoiddeforest2.tif")
SF2020_avoiddeforest2 <- sum(uSF2020_avoiddeforest2, DSF2020_avoiddeforest2)

### mean sf cover in local scale (90m)
SF2020_avoiddeforest2.px <- focal(SF2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_avoiddeforest2.px)<-"SFpx"
SF2020_avoiddeforest2.px[is.nan(SF2020_avoiddeforest2.px)] <- 0
SF2020_avoiddeforest2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_avoiddeforest2.px[])
#saving
writeRaster(SF2020_avoiddeforest2.px, "rasters/STM/2020_avoiddeforest2/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_avoiddeforest2.ls <- focal(SF2020_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_avoiddeforest2.ls)<-"SFls"
SF2020_avoiddeforest2.ls[is.nan(SF2020_avoiddeforest2.ls)] <- 0
SF2020_avoiddeforest2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_avoiddeforest2.ls[])
#saving
writeRaster(SF2020_avoiddeforest2.ls, "rasters/STM/2020_avoiddeforest2/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddeforest2.px", "SF2020_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_avoiddeforest2 <- calc(stm.sfage[["stm.sfage.2010real"]], fun=function(x){ifelse(is.na(x), x, ifelse(x==0, x, x+10))})
SFAge2020_avoiddeforest2[] <- ifelse(SF2020_avoiddeforest2[]==0, 0, SFAge2020_avoiddeforest2[])
writeRaster(SFAge2020_avoiddeforest2, "rasters/STM/input/SFAge2020_avoiddeforest2.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_avoiddeforest2.px <- focal(SFAge2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_avoiddeforest2.px)<-"SFAgepx"
SFAge2020_avoiddeforest2.px[is.nan(SFAge2020_avoiddeforest2.px)] <- 0
SFAge2020_avoiddeforest2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_avoiddeforest2.px[])
#saving
writeRaster(SFAge2020_avoiddeforest2.px, "rasters/STM/2020_avoiddeforest2/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_avoiddeforest2.ls <- focal(SFAge2020_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_avoiddeforest2.ls)<-"SFAgels"
SFAge2020_avoiddeforest2.ls[is.nan(SFAge2020_avoiddeforest2.ls)] <- 0
SFAge2020_avoiddeforest2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_avoiddeforest2.ls[])
#saving
writeRaster(SFAge2020_avoiddeforest2.ls, "rasters/STM/2020_avoiddeforest2/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_avoiddeforest2.px", "SFAge2020_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young_avoiddeforest2 <- SFAge2020_avoiddeforest2
SF2020young_avoiddeforest2[] <- ifelse(SF2020young_avoiddeforest2[]>2, 1, 0)
writeRaster(SF2020young_avoiddeforest2, "rasters/STM/input/SF2020young_avoiddeforest2.tif", format="GTiff", overwrite=T)

TF2020_avoiddeforest2 <- sum(UPF2020_avoiddeforest2, DPF2020_avoiddeforest2)
TF2020_avoiddeforest2 <- sum(TF2020_avoiddeforest2, SF2020young_avoiddeforest2)
writeRaster(TF2020_avoiddeforest2, "rasters/STM/input/TF2020_avoiddeforest2.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_avoiddeforest2.px <- focal(TF2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_avoiddeforest2.px)<-"TFpx"
TF2020_avoiddeforest2.px[is.nan(TF2020_avoiddeforest2.px)] <- 0
TF2020_avoiddeforest2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_avoiddeforest2.px[])
#saving
writeRaster(TF2020_avoiddeforest2.px, "rasters/STM/2020_avoiddeforest2/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_avoiddeforest2.ls <- focal(TF2020_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_avoiddeforest2.ls)<-"TFls"
TF2020_avoiddeforest2.ls[is.nan(TF2020_avoiddeforest2.ls)] <- 0
TF2020_avoiddeforest2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_avoiddeforest2.ls[])
#saving
writeRaster(TF2020_avoiddeforest2.ls, "rasters/STM/2020_avoiddeforest2/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddeforest2young", "SFAge2020_avoiddeforest2.px", "SFAge2020_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature_avoiddeforest2 <- SFAge2020_avoiddeforest2
SF2020mature_avoiddeforest2[] <- ifelse(SF2020mature_avoiddeforest2[]>5, 1, 0)
writeRaster(SF2020mature_avoiddeforest2, "rasters/STM/input/SF2020mature_avoiddeforest2.tif", format="GTiff", overwrite=T)

MF2020_avoiddeforest2 <- sum(UPF2020_avoiddeforest2, DPF2020_avoiddeforest2)
MF2020_avoiddeforest2 <- sum(MF2020_avoiddeforest2, SF2020mature_avoiddeforest2)
writeRaster(MF2020_avoiddeforest2, "rasters/STM/input/MF2020_avoiddeforest2.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_avoiddeforest2.px <- focal(MF2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_avoiddeforest2.px)<-"MFpx"
MF2020_avoiddeforest2.px[is.nan(MF2020_avoiddeforest2.px)] <- 0
MF2020_avoiddeforest2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_avoiddeforest2.px[])
#saving
writeRaster(MF2020_avoiddeforest2.px, "rasters/STM/2020_avoiddeforest2/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_avoiddeforest2.ls <- focal(MF2020_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_avoiddeforest2.ls)<-"MFls"
MF2020_avoiddeforest2.ls[is.nan(MF2020_avoiddeforest2.ls)] <- 0
MF2020_avoiddeforest2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_avoiddeforest2.ls[])
#saving
writeRaster(MF2020_avoiddeforest2.ls, "rasters/STM/2020_avoiddeforest2/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddeforest2mature", "SFAge2020_avoiddeforest2.px", "SFAge2020_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_avoiddeforest2 <- raster("rasters/STM/input/MF2020_avoiddeforest2.tif")
inv.MF2020_avoiddeforest2 <- MF2020_avoiddeforest2
inv.MF2020_avoiddeforest2[inv.MF2020_avoiddeforest2==1]<-NA
#cheking
#inv.MF2020_avoiddeforest2
#plot(inv.MF2020_avoiddeforest2)

edge.dist.2020_avoiddeforest2 <- distance(inv.MF2020_avoiddeforest2, doEdge=T)
names(edge.dist.2020_avoiddeforest2)<-"edgedist"
edge.dist.2020_avoiddeforest2[is.nan(edge.dist.2020_avoiddeforest2)] <- 0
edge.dist.2020_avoiddeforest2[] <- ifelse(stm.lulc[[1]][]==0, NA, edge.dist.2020_avoiddeforest2[])
#saving
writeRaster(edge.dist.2020_avoiddeforest2, "rasters/STM/2020_avoiddeforest2/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_avoiddeforest2")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_avoiddeforest2 <- edge.dist.2020_avoiddeforest2
edge2020_avoiddeforest2[] <- ifelse(edge2020_avoiddeforest2[] < 200, 0, ifelse(edge2020_avoiddeforest2[]>300, 0, 1))
writeRaster(edge2020_avoiddeforest2, "rasters/STM/input/edge2020_avoiddeforest2.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_avoiddeforest2.px <- focal(edge2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_avoiddeforest2.px)<-"edgepx"
edge2020_avoiddeforest2.px[is.nan(edge2020_avoiddeforest2.px)] <- 0
edge2020_avoiddeforest2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_avoiddeforest2.px[])
#saving
writeRaster(edge2020_avoiddeforest2.px, "rasters/STM/2020_avoiddeforest2/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_avoiddeforest2.ls <- focal(edge2020_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_avoiddeforest2.ls)<-"edgels"
edge2020_avoiddeforest2.ls[is.nan(edge2020_avoiddeforest2.ls)] <- 0
edge2020_avoiddeforest2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_avoiddeforest2.ls[])
#saving
writeRaster(edge2020_avoiddeforest2.ls, "rasters/STM/2020_avoiddeforest2/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_avoiddeforest2.px", "edge2020_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Restoration without avoiding ======================================
### Undegraded primary forest
#UPF2020_restor_wo_avoid <- raster("rasters/STM/input/LULC/UPF2020_restor_wo_avoid.tif")

### mean upf cover in local scale (90m)
UPF2020_restor_wo_avoid.px <- focal(UPF2020_restor_wo_avoid, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_restor_wo_avoid.px)<-"UPFpx"
UPF2020_restor_wo_avoid.px[is.nan(UPF2020_restor_wo_avoid.px)] <- 0
UPF2020_restor_wo_avoid.px[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_restor_wo_avoid.px[])
#saving
writeRaster(UPF2020_restor_wo_avoid.px, "rasters/STM/2020_restor_wo_avoid/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_restor_wo_avoid.ls <- focal(UPF2020_restor_wo_avoid, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_restor_wo_avoid.ls)<-"UPFls"
UPF2020_restor_wo_avoid.ls[is.nan(UPF2020_restor_wo_avoid.ls)] <- 0
UPF2020_restor_wo_avoid.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_restor_wo_avoid.ls[])
#saving
writeRaster(UPF2020_restor_wo_avoid.ls, "rasters/STM/2020_restor_wo_avoid/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_restor_wo_avoid.px", "UPF2020_restor_wo_avoid.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_restor_wo_avoid <- raster("rasters/STM/input/LULC/uDPF2020_restor_wo_avoid.tif")
#RDPF2020_restor_wo_avoid <- raster("rasters/STM/input/LULC/RDPF2020_restor_wo_avoid.tif")
DPF2020_restor_wo_avoid <- sum(uDPF2020_restor_wo_avoid, RDPF2020_restor_wo_avoid)

### mean dpf cover in local scale (90m)
DPF2020_restor_wo_avoid.px <- focal(DPF2020_restor_wo_avoid, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_restor_wo_avoid.px)<-"DPFpx"
DPF2020_restor_wo_avoid.px[is.nan(DPF2020_restor_wo_avoid.px)] <- 0
DPF2020_restor_wo_avoid.px[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_restor_wo_avoid.px[])
#saving
writeRaster(DPF2020_restor_wo_avoid.px, "rasters/STM/2020_restor_wo_avoid/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_restor_wo_avoid.ls <- focal(DPF2020_restor_wo_avoid, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_restor_wo_avoid.ls)<-"DPFls"
DPF2020_restor_wo_avoid.ls[is.nan(DPF2020_restor_wo_avoid.ls)] <- 0
DPF2020_restor_wo_avoid.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_restor_wo_avoid.ls[])
#saving
writeRaster(DPF2020_restor_wo_avoid.ls, "rasters/STM/2020_restor_wo_avoid/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_restor_wo_avoid.px", "DPF2020_restor_wo_avoid.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_restor_wo_avoid <- TSD2020
writeRaster(TSD2020_restor_wo_avoid, "rasters/STM/input/TSD2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_restor_wo_avoid.px <- focal(TSD2020_restor_wo_avoid, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_restor_wo_avoid.px)<-"TSDpx"
TSD2020_restor_wo_avoid.px[is.nan(TSD2020_restor_wo_avoid.px)] <- 0
TSD2020_restor_wo_avoid.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_restor_wo_avoid.px[])
#saving
writeRaster(TSD2020_restor_wo_avoid.px, "rasters/STM/2020_restor_wo_avoid/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_restor_wo_avoid.ls <- focal(TSD2020_restor_wo_avoid, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_restor_wo_avoid.ls)<-"TSDls"
TSD2020_restor_wo_avoid.ls[is.nan(TSD2020_restor_wo_avoid.ls)] <- 0
TSD2020_restor_wo_avoid.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_restor_wo_avoid.ls[])
#saving
writeRaster(TSD2020_restor_wo_avoid.ls, "rasters/STM/2020_restor_wo_avoid/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_restor_wo_avoid.px", "TSD2020_restor_wo_avoid.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_restor_wo_avoid <- raster("rasters/STM/input/LULC/uSF2020_restor_wo_avoid.tif")
#DSF2020_restor_wo_avoid <- raster("rasters/STM/input/LULC/DSF2020_restor_wo_avoid.tif")
SF2020_restor_wo_avoid <- sum(uSF2020_restor_wo_avoid, DSF2020_restor_wo_avoid)

### mean sf cover in local scale (90m)
SF2020_restor_wo_avoid.px <- focal(SF2020_restor_wo_avoid, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_restor_wo_avoid.px)<-"SFpx"
SF2020_restor_wo_avoid.px[is.nan(SF2020_restor_wo_avoid.px)] <- 0
SF2020_restor_wo_avoid.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_restor_wo_avoid.px[])
#saving
writeRaster(SF2020_restor_wo_avoid.px, "rasters/STM/2020_restor_wo_avoid/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_restor_wo_avoid.ls <- focal(SF2020_restor_wo_avoid, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_restor_wo_avoid.ls)<-"SFls"
SF2020_restor_wo_avoid.ls[is.nan(SF2020_restor_wo_avoid.ls)] <- 0
SF2020_restor_wo_avoid.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_restor_wo_avoid.ls[])
#saving
writeRaster(SF2020_restor_wo_avoid.ls, "rasters/STM/2020_restor_wo_avoid/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_wo_avoid.px", "SF2020_restor_wo_avoid.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_restor_wo_avoid <- stm.sfage[["stm.sfage.2020real"]]
SFAge2020_restor_wo_avoid[] <- ifelse(SFAge2020_restor_wo_avoid[]==0 & SF2020_restor_wo_avoid[]==1, 10, SFAge2020_restor_wo_avoid[])
writeRaster(SFAge2020_restor_wo_avoid, "rasters/STM/input/SFAge2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_restor_wo_avoid.px <- focal(SFAge2020_restor_wo_avoid, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_restor_wo_avoid.px)<-"SFAgepx"
SFAge2020_restor_wo_avoid.px[is.nan(SFAge2020_restor_wo_avoid.px)] <- 0
SFAge2020_restor_wo_avoid.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_restor_wo_avoid.px[])
#saving
writeRaster(SFAge2020_restor_wo_avoid.px, "rasters/STM/2020_restor_wo_avoid/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_restor_wo_avoid.ls <- focal(SFAge2020_restor_wo_avoid, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_restor_wo_avoid.ls)<-"SFAgels"
SFAge2020_restor_wo_avoid.ls[is.nan(SFAge2020_restor_wo_avoid.ls)] <- 0
SFAge2020_restor_wo_avoid.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_restor_wo_avoid.ls[])
#saving
writeRaster(SFAge2020_restor_wo_avoid.ls, "rasters/STM/2020_restor_wo_avoid/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_restor_wo_avoid.px", "SFAge2020_restor_wo_avoid.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020_restor_wo_avoidyoung <- SFAge2020_restor_wo_avoid
SF2020_restor_wo_avoidyoung[] <- ifelse(SF2020_restor_wo_avoidyoung[]>2, 1, 0)
writeRaster(SF2020_restor_wo_avoidyoung, "rasters/STM/input/SF2020_restor_wo_avoidyoung.tif", format="GTiff", overwrite=T)

TF2020_restor_wo_avoid <- sum(UPF2020_restor_wo_avoid, DPF2020_restor_wo_avoid)
TF2020_restor_wo_avoid <- sum(TF2020_restor_wo_avoid, SF2020_restor_wo_avoidyoung)
writeRaster(TF2020_restor_wo_avoid, "rasters/STM/input/TF2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_restor_wo_avoid.px <- focal(TF2020_restor_wo_avoid, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_restor_wo_avoid.px)<-"TFpx"
TF2020_restor_wo_avoid.px[is.nan(TF2020_restor_wo_avoid.px)] <- 0
TF2020_restor_wo_avoid.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_restor_wo_avoid.px[])
#saving
writeRaster(TF2020_restor_wo_avoid.px, "rasters/STM/2020_restor_wo_avoid/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_restor_wo_avoid.ls <- focal(TF2020_restor_wo_avoid, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_restor_wo_avoid.ls)<-"TFls"
TF2020_restor_wo_avoid.ls[is.nan(TF2020_restor_wo_avoid.ls)] <- 0
TF2020_restor_wo_avoid.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_restor_wo_avoid.ls[])
#saving
writeRaster(TF2020_restor_wo_avoid.ls, "rasters/STM/2020_restor_wo_avoid/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_wo_avoidyoung", "SFAge2020_restor_wo_avoid.px", "SFAge2020_restor_wo_avoid.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020_restor_wo_avoidmature <- SFAge2020_restor_wo_avoid
SF2020_restor_wo_avoidmature[] <- ifelse(SF2020_restor_wo_avoidmature[]>5, 1, 0)
writeRaster(SF2020_restor_wo_avoidmature, "rasters/STM/input/SF2020_restor_wo_avoidmature.tif", format="GTiff", overwrite=T)

MF2020_restor_wo_avoid <- sum(UPF2020_restor_wo_avoid, DPF2020_restor_wo_avoid)
MF2020_restor_wo_avoid <- sum(MF2020_restor_wo_avoid, SF2020_restor_wo_avoidmature)
writeRaster(MF2020_restor_wo_avoid, "rasters/STM/input/MF2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_restor_wo_avoid.px <- focal(MF2020_restor_wo_avoid, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_restor_wo_avoid.px)<-"MFpx"
MF2020_restor_wo_avoid.px[is.nan(MF2020_restor_wo_avoid.px)] <- 0
MF2020_restor_wo_avoid.px[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_restor_wo_avoid.px[])
#saving
writeRaster(MF2020_restor_wo_avoid.px, "rasters/STM/2020_restor_wo_avoid/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_restor_wo_avoid.ls <- focal(MF2020_restor_wo_avoid, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_restor_wo_avoid.ls)<-"MFls"
MF2020_restor_wo_avoid.ls[is.nan(MF2020_restor_wo_avoid.ls)] <- 0
MF2020_restor_wo_avoid.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_restor_wo_avoid.ls[])
#saving
writeRaster(MF2020_restor_wo_avoid.ls, "rasters/STM/2020_restor_wo_avoid/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_wo_avoidmature", "SFAge2020_restor_wo_avoid.px", "SFAge2020_restor_wo_avoid.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_restor_wo_avoid <- raster("rasters/STM/input/MF2020_restor_wo_avoid.tif")
inv.MF2020_restor_wo_avoid <- MF2020_restor_wo_avoid
inv.MF2020_restor_wo_avoid[inv.MF2020_restor_wo_avoid==1]<-NA
#cheking
#inv.MF2020_restor_wo_avoid
#plot(inv.MF2020_restor_wo_avoid)

edge.dist.2020_restor_wo_avoid <- distance(inv.MF2020_restor_wo_avoid, doEdge=T)
names(edge.dist.2020_restor_wo_avoid)<-"edgedist"
edge.dist.2020_restor_wo_avoid[is.nan(edge.dist.2020_restor_wo_avoid)] <- 0
edge.dist.2020_restor_wo_avoid[] <- ifelse(stm.lulc[[1]][]==0, NA, edge.dist.2020_restor_wo_avoid[])
#saving
writeRaster(edge.dist.2020_restor_wo_avoid, "rasters/STM/2020_restor_wo_avoid/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_restor_wo_avoid")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_restor_wo_avoid <- edge.dist.2020_restor_wo_avoid
edge2020_restor_wo_avoid[] <- ifelse(edge2020_restor_wo_avoid[] < 200, 0, ifelse(edge2020_restor_wo_avoid[]>300, 0, 1))
writeRaster(edge2020_restor_wo_avoid, "rasters/STM/input/edge2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_restor_wo_avoid.px <- focal(edge2020_restor_wo_avoid, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_restor_wo_avoid.px)<-"edgepx"
edge2020_restor_wo_avoid.px[is.nan(edge2020_restor_wo_avoid.px)] <- 0
edge2020_restor_wo_avoid.px[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_restor_wo_avoid.px[])
#saving
writeRaster(edge2020_restor_wo_avoid.px, "rasters/STM/2020_restor_wo_avoid/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_restor_wo_avoid.ls <- focal(edge2020_restor_wo_avoid, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_restor_wo_avoid.ls)<-"edgels"
edge2020_restor_wo_avoid.ls[is.nan(edge2020_restor_wo_avoid.ls)] <- 0
edge2020_restor_wo_avoid.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_restor_wo_avoid.ls[])
#saving
writeRaster(edge2020_restor_wo_avoid.ls, "rasters/STM/2020_restor_wo_avoid/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_restor_wo_avoid.px", "edge2020_restor_wo_avoid.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Avoid both (all) ==================================================
### Undegraded primary forest
#UPF2020_avoidboth <- raster("rasters/STM/input/LULC/UPF2020_avoidboth.tif")

### mean upf cover in local scale (90m)
UPF2020_avoidboth.px <- focal(UPF2020_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_avoidboth.px)<-"UPFpx"
UPF2020_avoidboth.px[is.nan(UPF2020_avoidboth.px)] <- 0
UPF2020_avoidboth.px[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_avoidboth.px[])
#saving
writeRaster(UPF2020_avoidboth.px, "rasters/STM/2020_avoidboth/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_avoidboth.ls <- focal(UPF2020_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_avoidboth.ls)<-"UPFls"
UPF2020_avoidboth.ls[is.nan(UPF2020_avoidboth.ls)] <- 0
UPF2020_avoidboth.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_avoidboth.ls[])
#saving
writeRaster(UPF2020_avoidboth.ls, "rasters/STM/2020_avoidboth/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_avoidboth.px", "UPF2020_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_avoidboth <- raster("rasters/STM/input/LULC/uDPF2020_avoidboth.tif")
#RDPF2020_avoidboth <- raster("rasters/STM/input/LULC/RDPF2020_avoidboth.tif")
DPF2020_avoidboth <- sum(uDPF2020_avoidboth, RDPF2020_avoidboth)

### mean dpf cover in local scale (90m)
DPF2020_avoidboth.px <- focal(DPF2020_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_avoidboth.px)<-"DPFpx"
DPF2020_avoidboth.px[is.nan(DPF2020_avoidboth.px)] <- 0
DPF2020_avoidboth.px[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_avoidboth.px[])
#saving
writeRaster(DPF2020_avoidboth.px, "rasters/STM/2020_avoidboth/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_avoidboth.ls <- focal(DPF2020_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_avoidboth.ls)<-"DPFls"
DPF2020_avoidboth.ls[is.nan(DPF2020_avoidboth.ls)] <- 0
DPF2020_avoidboth.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_avoidboth.ls[])
#saving
writeRaster(DPF2020_avoidboth.ls, "rasters/STM/2020_avoidboth/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_avoidboth.px", "DPF2020_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_avoidboth <- TSD2020_avoiddegrad
writeRaster(TSD2020_avoidboth, "rasters/STM/input/TSD2020_avoidboth.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_avoidboth.px <- focal(TSD2020_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_avoidboth.px)<-"TSDpx"
TSD2020_avoidboth.px[is.nan(TSD2020_avoidboth.px)] <- 0
TSD2020_avoidboth.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_avoidboth.px[])
#saving
writeRaster(TSD2020_avoidboth.px, "rasters/STM/2020_avoidboth/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_avoidboth.ls <- focal(TSD2020_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_avoidboth.ls)<-"TSDls"
TSD2020_avoidboth.ls[is.nan(TSD2020_avoidboth.ls)] <- 0
TSD2020_avoidboth.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_avoidboth.ls[])
#saving
writeRaster(TSD2020_avoidboth.ls, "rasters/STM/2020_avoidboth/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_avoidboth.px", "TSD2020_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_avoidboth <- raster("rasters/STM/input/LULC/uSF2020_avoidboth.tif")
#DSF2020_avoidboth <- raster("rasters/STM/input/LULC/DSF2020_avoidboth.tif")
SF2020_avoidboth <- sum(uSF2020_avoidboth, DSF2020_avoidboth)

### mean sf cover in local scale (90m)
SF2020_avoidboth.px <- focal(SF2020_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_avoidboth.px)<-"SFpx"
SF2020_avoidboth.px[is.nan(SF2020_avoidboth.px)] <- 0
SF2020_avoidboth.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_avoidboth.px[])
#saving
writeRaster(SF2020_avoidboth.px, "rasters/STM/2020_avoidboth/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_avoidboth.ls <- focal(SF2020_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_avoidboth.ls)<-"SFls"
SF2020_avoidboth.ls[is.nan(SF2020_avoidboth.ls)] <- 0
SF2020_avoidboth.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_avoidboth.ls[])
#saving
writeRaster(SF2020_avoidboth.ls, "rasters/STM/2020_avoidboth/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoidboth.px", "SF2020_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_avoidboth <- SFAge2020_avoiddeforest
writeRaster(SFAge2020_avoidboth, "rasters/STM/input/SFAge2020_avoidboth.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_avoidboth.px <- focal(SFAge2020_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_avoidboth.px)<-"SFAgepx"
SFAge2020_avoidboth.px[is.nan(SFAge2020_avoidboth.px)] <- 0
SFAge2020_avoidboth.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_avoidboth.px[])
#saving
writeRaster(SFAge2020_avoidboth.px, "rasters/STM/2020_avoidboth/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_avoidboth.ls <- focal(SFAge2020_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_avoidboth.ls)<-"SFAgels"
SFAge2020_avoidboth.ls[is.nan(SFAge2020_avoidboth.ls)] <- 0
SFAge2020_avoidboth.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_avoidboth.ls[])
#saving
writeRaster(SFAge2020_avoidboth.ls, "rasters/STM/2020_avoidboth/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_avoidboth.px", "SFAge2020_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young_avoidboth <- SFAge2020_avoidboth
SF2020young_avoidboth[] <- ifelse(SF2020young_avoidboth[]>2, 1, 0)
writeRaster(SF2020young_avoidboth, "rasters/STM/input/SF2020young_avoidboth.tif", format="GTiff", overwrite=T)

TF2020_avoidboth <- sum(UPF2020_avoidboth, DPF2020_avoidboth)
TF2020_avoidboth <- sum(TF2020_avoidboth, SF2020young_avoidboth)
writeRaster(TF2020_avoidboth, "rasters/STM/input/TF2020_avoidboth.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_avoidboth.px <- focal(TF2020_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_avoidboth.px)<-"TFpx"
TF2020_avoidboth.px[is.nan(TF2020_avoidboth.px)] <- 0
TF2020_avoidboth.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_avoidboth.px[])
#saving
writeRaster(TF2020_avoidboth.px, "rasters/STM/2020_avoidboth/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_avoidboth.ls <- focal(TF2020_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_avoidboth.ls)<-"TFls"
TF2020_avoidboth.ls[is.nan(TF2020_avoidboth.ls)] <- 0
TF2020_avoidboth.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_avoidboth.ls[])
#saving
writeRaster(TF2020_avoidboth.ls, "rasters/STM/2020_avoidboth/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoidbothyoung", "SFAge2020_avoidboth.px", "SFAge2020_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature_avoidboth <- SFAge2020_avoidboth
SF2020mature_avoidboth[] <- ifelse(SF2020mature_avoidboth[]>5, 1, 0)
writeRaster(SF2020mature_avoidboth, "rasters/STM/input/SF2020mature_avoidboth.tif", format="GTiff", overwrite=T)

MF2020_avoidboth <- sum(UPF2020_avoidboth, DPF2020_avoidboth)
MF2020_avoidboth <- sum(MF2020_avoidboth, SF2020mature_avoidboth)
writeRaster(MF2020_avoidboth, "rasters/STM/input/MF2020_avoidboth.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_avoidboth.px <- focal(MF2020_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_avoidboth.px)<-"MFpx"
MF2020_avoidboth.px[is.nan(MF2020_avoidboth.px)] <- 0
MF2020_avoidboth.px[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_avoidboth.px[])
#saving
writeRaster(MF2020_avoidboth.px, "rasters/STM/2020_avoidboth/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_avoidboth.ls <- focal(MF2020_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_avoidboth.ls)<-"MFls"
MF2020_avoidboth.ls[is.nan(MF2020_avoidboth.ls)] <- 0
MF2020_avoidboth.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_avoidboth.ls[])
#saving
writeRaster(MF2020_avoidboth.ls, "rasters/STM/2020_avoidboth/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoidbothmature", "SFAge2020_avoidboth.px", "SFAge2020_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_avoidboth <- raster("rasters/STM/input/MF2020_avoidboth.tif")
inv.MF2020_avoidboth <- MF2020_avoidboth
inv.MF2020_avoidboth[inv.MF2020_avoidboth==1]<-NA
#cheking
#inv.MF2020_avoidboth
#plot(inv.MF2020_avoidboth)

edge.dist.2020_avoidboth <- distance(inv.MF2020_avoidboth, doEdge=T)
names(edge.dist.2020_avoidboth)<-"edgedist"
edge.dist.2020_avoidboth[is.nan(edge.dist.2020_avoidboth)] <- 0
edge.dist.2020_avoidboth[] <- ifelse(stm.lulc[[1]][]==0, NA, edge.dist.2020_avoidboth[])
#saving
writeRaster(edge.dist.2020_avoidboth, "rasters/STM/2020_avoidboth/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_avoidboth")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_avoidboth <- edge.dist.2020_avoidboth
edge2020_avoidboth[] <- ifelse(edge2020_avoidboth[] < 200, 0, ifelse(edge2020_avoidboth[]>300, 0, 1))
writeRaster(edge2020_avoidboth, "rasters/STM/input/edge2020_avoidboth.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_avoidboth.px <- focal(edge2020_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_avoidboth.px)<-"edgepx"
edge2020_avoidboth.px[is.nan(edge2020_avoidboth.px)] <- 0
edge2020_avoidboth.px[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_avoidboth.px[])
#saving
writeRaster(edge2020_avoidboth.px, "rasters/STM/2020_avoidboth/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_avoidboth.ls <- focal(edge2020_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_avoidboth.ls)<-"edgels"
edge2020_avoidboth.ls[is.nan(edge2020_avoidboth.ls)] <- 0
edge2020_avoidboth.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_avoidboth.ls[])
#saving
writeRaster(edge2020_avoidboth.ls, "rasters/STM/2020_avoidboth/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_avoidboth.px", "edge2020_avoidboth.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Avoid both (Primary forests only) =================================
### Undegraded primary forest
#UPF2020_avoidboth2 <- raster("rasters/STM/input/LULC/UPF2020_avoidboth2.tif")

### mean upf cover in local scale (90m)
UPF2020_avoidboth2.px <- focal(UPF2020_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_avoidboth2.px)<-"UPFpx"
UPF2020_avoidboth2.px[is.nan(UPF2020_avoidboth2.px)] <- 0
UPF2020_avoidboth2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_avoidboth2.px[])
#saving
writeRaster(UPF2020_avoidboth2.px, "rasters/STM/2020_avoidboth2/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_avoidboth2.ls <- focal(UPF2020_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_avoidboth2.ls)<-"UPFls"
UPF2020_avoidboth2.ls[is.nan(UPF2020_avoidboth2.ls)] <- 0
UPF2020_avoidboth2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_avoidboth2.ls[])
#saving
writeRaster(UPF2020_avoidboth2.ls, "rasters/STM/2020_avoidboth2/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_avoidboth2.px", "UPF2020_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_avoidboth2 <- raster("rasters/STM/input/LULC/uDPF2020_avoidboth2.tif")
#RDPF2020_avoidboth2 <- raster("rasters/STM/input/LULC/RDPF2020_avoidboth2.tif")
DPF2020_avoidboth2 <- sum(uDPF2020_avoidboth2, RDPF2020_avoidboth2)

### mean dpf cover in local scale (90m)
DPF2020_avoidboth2.px <- focal(DPF2020_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_avoidboth2.px)<-"DPFpx"
DPF2020_avoidboth2.px[is.nan(DPF2020_avoidboth2.px)] <- 0
DPF2020_avoidboth2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_avoidboth2.px[])
#saving
writeRaster(DPF2020_avoidboth2.px, "rasters/STM/2020_avoidboth2/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_avoidboth2.ls <- focal(DPF2020_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_avoidboth2.ls)<-"DPFls"
DPF2020_avoidboth2.ls[is.nan(DPF2020_avoidboth2.ls)] <- 0
DPF2020_avoidboth2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_avoidboth2.ls[])
#saving
writeRaster(DPF2020_avoidboth2.ls, "rasters/STM/2020_avoidboth2/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_avoidboth2.px", "DPF2020_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_avoidboth2 <- TSD2020_avoiddegrad2
writeRaster(TSD2020_avoidboth2, "rasters/STM/input/TSD2020_avoidboth2.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_avoidboth2.px <- focal(TSD2020_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_avoidboth2.px)<-"TSDpx"
TSD2020_avoidboth2.px[is.nan(TSD2020_avoidboth2.px)] <- 0
TSD2020_avoidboth2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_avoidboth2.px[])
#saving
writeRaster(TSD2020_avoidboth2.px, "rasters/STM/2020_avoidboth2/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_avoidboth2.ls <- focal(TSD2020_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_avoidboth2.ls)<-"TSDls"
TSD2020_avoidboth2.ls[is.nan(TSD2020_avoidboth2.ls)] <- 0
TSD2020_avoidboth2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_avoidboth2.ls[])
#saving
writeRaster(TSD2020_avoidboth2.ls, "rasters/STM/2020_avoidboth2/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_avoidboth2.px", "TSD2020_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_avoidboth2 <- raster("rasters/STM/input/LULC/uSF2020_avoidboth2.tif")
#DSF2020_avoidboth2 <- raster("rasters/STM/input/LULC/DSF2020_avoidboth2.tif")
SF2020_avoidboth2 <- sum(uSF2020_avoidboth2, DSF2020_avoidboth2)

### mean sf cover in local scale (90m)
SF2020_avoidboth2.px <- focal(SF2020_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_avoidboth2.px)<-"SFpx"
SF2020_avoidboth2.px[is.nan(SF2020_avoidboth2.px)] <- 0
SF2020_avoidboth2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_avoidboth2.px[])
#saving
writeRaster(SF2020_avoidboth2.px, "rasters/STM/2020_avoidboth2/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_avoidboth2.ls <- focal(SF2020_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_avoidboth2.ls)<-"SFls"
SF2020_avoidboth2.ls[is.nan(SF2020_avoidboth2.ls)] <- 0
SF2020_avoidboth2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_avoidboth2.ls[])
#saving
writeRaster(SF2020_avoidboth2.ls, "rasters/STM/2020_avoidboth2/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoidboth2.px", "SF2020_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_avoidboth2 <- SFAge2020_avoiddeforest2
writeRaster(SFAge2020_avoidboth2, "rasters/STM/input/SFAge2020_avoidboth2.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_avoidboth2.px <- focal(SFAge2020_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_avoidboth2.px)<-"SFAgepx"
SFAge2020_avoidboth2.px[is.nan(SFAge2020_avoidboth2.px)] <- 0
SFAge2020_avoidboth2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_avoidboth2.px[])
#saving
writeRaster(SFAge2020_avoidboth2.px, "rasters/STM/2020_avoidboth2/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_avoidboth2.ls <- focal(SFAge2020_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_avoidboth2.ls)<-"SFAgels"
SFAge2020_avoidboth2.ls[is.nan(SFAge2020_avoidboth2.ls)] <- 0
SFAge2020_avoidboth2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_avoidboth2.ls[])
#saving
writeRaster(SFAge2020_avoidboth2.ls, "rasters/STM/2020_avoidboth2/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_avoidboth2.px", "SFAge2020_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young_avoidboth2 <- SFAge2020_avoidboth2
SF2020young_avoidboth2[] <- ifelse(SF2020young_avoidboth2[]>2, 1, 0)
writeRaster(SF2020young_avoidboth2, "rasters/STM/input/SF2020young_avoidboth2.tif", format="GTiff", overwrite=T)

TF2020_avoidboth2 <- sum(UPF2020_avoidboth2, DPF2020_avoidboth2)
TF2020_avoidboth2 <- sum(TF2020_avoidboth2, SF2020young_avoidboth2)
writeRaster(TF2020_avoidboth2, "rasters/STM/input/TF2020_avoidboth2.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_avoidboth2.px <- focal(TF2020_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_avoidboth2.px)<-"TFpx"
TF2020_avoidboth2.px[is.nan(TF2020_avoidboth2.px)] <- 0
TF2020_avoidboth2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_avoidboth2.px[])
#saving
writeRaster(TF2020_avoidboth2.px, "rasters/STM/2020_avoidboth2/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_avoidboth2.ls <- focal(TF2020_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_avoidboth2.ls)<-"TFls"
TF2020_avoidboth2.ls[is.nan(TF2020_avoidboth2.ls)] <- 0
TF2020_avoidboth2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_avoidboth2.ls[])
#saving
writeRaster(TF2020_avoidboth2.ls, "rasters/STM/2020_avoidboth2/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoidboth2young", "SFAge2020_avoidboth2.px", "SFAge2020_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature_avoidboth2 <- SFAge2020_avoidboth2
SF2020mature_avoidboth2[] <- ifelse(SF2020mature_avoidboth2[]>5, 1, 0)
writeRaster(SF2020mature_avoidboth2, "rasters/STM/input/SF2020mature_avoidboth2.tif", format="GTiff", overwrite=T)

MF2020_avoidboth2 <- sum(UPF2020_avoidboth2, DPF2020_avoidboth2)
MF2020_avoidboth2 <- sum(MF2020_avoidboth2, SF2020mature_avoidboth2)
writeRaster(MF2020_avoidboth2, "rasters/STM/input/MF2020_avoidboth2.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_avoidboth2.px <- focal(MF2020_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_avoidboth2.px)<-"MFpx"
MF2020_avoidboth2.px[is.nan(MF2020_avoidboth2.px)] <- 0
MF2020_avoidboth2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_avoidboth2.px[])
#saving
writeRaster(MF2020_avoidboth2.px, "rasters/STM/2020_avoidboth2/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_avoidboth2.ls <- focal(MF2020_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_avoidboth2.ls)<-"MFls"
MF2020_avoidboth2.ls[is.nan(MF2020_avoidboth2.ls)] <- 0
MF2020_avoidboth2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_avoidboth2.ls[])
#saving
writeRaster(MF2020_avoidboth2.ls, "rasters/STM/2020_avoidboth2/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoidboth2mature", "SFAge2020_avoidboth2.px", "SFAge2020_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_avoidboth2 <- raster("rasters/STM/input/MF2020_avoidboth2.tif")
inv.MF2020_avoidboth2 <- MF2020_avoidboth2
inv.MF2020_avoidboth2[inv.MF2020_avoidboth2==1]<-NA
#cheking
#inv.MF2020_avoidboth2
#plot(inv.MF2020_avoidboth2)

edge.dist.2020_avoidboth2 <- distance(inv.MF2020_avoidboth2, doEdge=T)
names(edge.dist.2020_avoidboth2)<-"edgedist"
edge.dist.2020_avoidboth2[is.nan(edge.dist.2020_avoidboth2)] <- 0
edge.dist.2020_avoidboth2[] <- ifelse(stm.lulc[[1]][]==0, NA, edge.dist.2020_avoidboth2[])
#saving
writeRaster(edge.dist.2020_avoidboth2, "rasters/STM/2020_avoidboth2/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_avoidboth2")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_avoidboth2 <- edge.dist.2020_avoidboth2
edge2020_avoidboth2[] <- ifelse(edge2020_avoidboth2[] < 200, 0, ifelse(edge2020_avoidboth2[]>300, 0, 1))
writeRaster(edge2020_avoidboth2, "rasters/STM/input/edge2020_avoidboth2.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_avoidboth2.px <- focal(edge2020_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_avoidboth2.px)<-"edgepx"
edge2020_avoidboth2.px[is.nan(edge2020_avoidboth2.px)] <- 0
edge2020_avoidboth2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_avoidboth2.px[])
#saving
writeRaster(edge2020_avoidboth2.px, "rasters/STM/2020_avoidboth2/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_avoidboth2.ls <- focal(edge2020_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_avoidboth2.ls)<-"edgels"
edge2020_avoidboth2.ls[is.nan(edge2020_avoidboth2.ls)] <- 0
edge2020_avoidboth2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_avoidboth2.ls[])
#saving
writeRaster(edge2020_avoidboth2.ls, "rasters/STM/2020_avoidboth2/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_avoidboth2.px", "edge2020_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Restoration and avoid deforestation (all) =========================
### Undegraded primary forest
#UPF2020_restor_n_avoiddeforest <- raster("rasters/STM/input/LULC/UPF2020_restor_n_avoiddeforest.tif")

### mean upf cover in local scale (90m)
UPF2020_restor_n_avoiddeforest.px <- focal(UPF2020_restor_n_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_restor_n_avoiddeforest.px)<-"UPFpx"
UPF2020_restor_n_avoiddeforest.px[is.nan(UPF2020_restor_n_avoiddeforest.px)] <- 0
UPF2020_restor_n_avoiddeforest.px[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_restor_n_avoiddeforest.px[])
#saving
writeRaster(UPF2020_restor_n_avoiddeforest.px, "rasters/STM/2020_restor_n_avoiddeforest/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_restor_n_avoiddeforest.ls <- focal(UPF2020_restor_n_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_restor_n_avoiddeforest.ls)<-"UPFls"
UPF2020_restor_n_avoiddeforest.ls[is.nan(UPF2020_restor_n_avoiddeforest.ls)] <- 0
UPF2020_restor_n_avoiddeforest.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_restor_n_avoiddeforest.ls[])
#saving
writeRaster(UPF2020_restor_n_avoiddeforest.ls, "rasters/STM/2020_restor_n_avoiddeforest/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_restor_n_avoiddeforest.px", "UPF2020_restor_n_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_restor_n_avoiddeforest <- raster("rasters/STM/input/LULC/uDPF2020_restor_n_avoiddeforest.tif")
#RDPF2020_restor_n_avoiddeforest <- raster("rasters/STM/input/LULC/RDPF2020_restor_n_avoiddeforest.tif")
DPF2020_restor_n_avoiddeforest <- sum(uDPF2020_restor_n_avoiddeforest, RDPF2020_restor_n_avoiddeforest)

### mean dpf cover in local scale (90m)
DPF2020_restor_n_avoiddeforest.px <- focal(DPF2020_restor_n_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_restor_n_avoiddeforest.px)<-"DPFpx"
DPF2020_restor_n_avoiddeforest.px[is.nan(DPF2020_restor_n_avoiddeforest.px)] <- 0
DPF2020_restor_n_avoiddeforest.px[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_restor_n_avoiddeforest.px[])
#saving
writeRaster(DPF2020_restor_n_avoiddeforest.px, "rasters/STM/2020_restor_n_avoiddeforest/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_restor_n_avoiddeforest.ls <- focal(DPF2020_restor_n_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_restor_n_avoiddeforest.ls)<-"DPFls"
DPF2020_restor_n_avoiddeforest.ls[is.nan(DPF2020_restor_n_avoiddeforest.ls)] <- 0
DPF2020_restor_n_avoiddeforest.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_restor_n_avoiddeforest.ls[])
#saving
writeRaster(DPF2020_restor_n_avoiddeforest.ls, "rasters/STM/2020_restor_n_avoiddeforest/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_restor_n_avoiddeforest.px", "DPF2020_restor_n_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_restor_n_avoiddeforest <- TSD2020_avoiddeforest
writeRaster(TSD2020_restor_n_avoiddeforest, "rasters/STM/input/TSD2020_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_restor_n_avoiddeforest.px <- focal(TSD2020_restor_n_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_restor_n_avoiddeforest.px)<-"TSDpx"
TSD2020_restor_n_avoiddeforest.px[is.nan(TSD2020_restor_n_avoiddeforest.px)] <- 0
TSD2020_restor_n_avoiddeforest.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_restor_n_avoiddeforest.px[])
#saving
writeRaster(TSD2020_restor_n_avoiddeforest.px, "rasters/STM/2020_restor_n_avoiddeforest/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_restor_n_avoiddeforest.ls <- focal(TSD2020_restor_n_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_restor_n_avoiddeforest.ls)<-"TSDls"
TSD2020_restor_n_avoiddeforest.ls[is.nan(TSD2020_restor_n_avoiddeforest.ls)] <- 0
TSD2020_restor_n_avoiddeforest.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_restor_n_avoiddeforest.ls[])
#saving
writeRaster(TSD2020_restor_n_avoiddeforest.ls, "rasters/STM/2020_restor_n_avoiddeforest/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_restor_n_avoiddeforest.px", "TSD2020_restor_n_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_restor_n_avoiddeforest <- raster("rasters/STM/input/LULC/uSF2020_restor_n_avoiddeforest.tif")
#DSF2020_restor_n_avoiddeforest <- raster("rasters/STM/input/LULC/DSF2020_restor_n_avoiddeforest.tif")
SF2020_restor_n_avoiddeforest <- sum(uSF2020_restor_n_avoiddeforest, DSF2020_restor_n_avoiddeforest)

### mean sf cover in local scale (90m)
SF2020_restor_n_avoiddeforest.px <- focal(SF2020_restor_n_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_restor_n_avoiddeforest.px)<-"SFpx"
SF2020_restor_n_avoiddeforest.px[is.nan(SF2020_restor_n_avoiddeforest.px)] <- 0
SF2020_restor_n_avoiddeforest.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_restor_n_avoiddeforest.px[])
#saving
writeRaster(SF2020_restor_n_avoiddeforest.px, "rasters/STM/2020_restor_n_avoiddeforest/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_restor_n_avoiddeforest.ls <- focal(SF2020_restor_n_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_restor_n_avoiddeforest.ls)<-"SFls"
SF2020_restor_n_avoiddeforest.ls[is.nan(SF2020_restor_n_avoiddeforest.ls)] <- 0
SF2020_restor_n_avoiddeforest.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_restor_n_avoiddeforest.ls[])
#saving
writeRaster(SF2020_restor_n_avoiddeforest.ls, "rasters/STM/2020_restor_n_avoiddeforest/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoiddeforest.px", "SF2020_restor_n_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_restor_n_avoiddeforest <- SFAge2020_avoiddeforest
SFAge2020_restor_n_avoiddeforest[] <- ifelse(SFAge2020_restor_n_avoiddeforest[]==0 & SF2020_restor_n_avoiddeforest[]==1, 10, SFAge2020_restor_n_avoiddeforest[])
writeRaster(SFAge2020_restor_n_avoiddeforest, "rasters/STM/input/SFAge2020_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_restor_n_avoiddeforest.px <- focal(SFAge2020_restor_n_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_restor_n_avoiddeforest.px)<-"SFAgepx"
SFAge2020_restor_n_avoiddeforest.px[is.nan(SFAge2020_restor_n_avoiddeforest.px)] <- 0
SFAge2020_restor_n_avoiddeforest.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_restor_n_avoiddeforest.px[])
#saving
writeRaster(SFAge2020_restor_n_avoiddeforest.px, "rasters/STM/2020_restor_n_avoiddeforest/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_restor_n_avoiddeforest.ls <- focal(SFAge2020_restor_n_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_restor_n_avoiddeforest.ls)<-"SFAgels"
SFAge2020_restor_n_avoiddeforest.ls[is.nan(SFAge2020_restor_n_avoiddeforest.ls)] <- 0
SFAge2020_restor_n_avoiddeforest.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_restor_n_avoiddeforest.ls[])
#saving
writeRaster(SFAge2020_restor_n_avoiddeforest.ls, "rasters/STM/2020_restor_n_avoiddeforest/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_restor_n_avoiddeforest.px", "SFAge2020_restor_n_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young_restor_n_avoiddeforest <- SFAge2020_restor_n_avoiddeforest
SF2020young_restor_n_avoiddeforest[] <- ifelse(SF2020young_restor_n_avoiddeforest[]>2, 1, 0)
writeRaster(SF2020young_restor_n_avoiddeforest, "rasters/STM/input/SF2020young_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)

TF2020_restor_n_avoiddeforest <- sum(UPF2020_restor_n_avoiddeforest, DPF2020_restor_n_avoiddeforest)
TF2020_restor_n_avoiddeforest <- sum(TF2020_restor_n_avoiddeforest, SF2020young_restor_n_avoiddeforest)
writeRaster(TF2020_restor_n_avoiddeforest, "rasters/STM/input/TF2020_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_restor_n_avoiddeforest.px <- focal(TF2020_restor_n_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_restor_n_avoiddeforest.px)<-"TFpx"
TF2020_restor_n_avoiddeforest.px[is.nan(TF2020_restor_n_avoiddeforest.px)] <- 0
TF2020_restor_n_avoiddeforest.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_restor_n_avoiddeforest.px[])
#saving
writeRaster(TF2020_restor_n_avoiddeforest.px, "rasters/STM/2020_restor_n_avoiddeforest/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_restor_n_avoiddeforest.ls <- focal(TF2020_restor_n_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_restor_n_avoiddeforest.ls)<-"TFls"
TF2020_restor_n_avoiddeforest.ls[is.nan(TF2020_restor_n_avoiddeforest.ls)] <- 0
TF2020_restor_n_avoiddeforest.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_restor_n_avoiddeforest.ls[])
#saving
writeRaster(TF2020_restor_n_avoiddeforest.ls, "rasters/STM/2020_restor_n_avoiddeforest/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoiddeforestyoung", "SFAge2020_restor_n_avoiddeforest.px", "SFAge2020_restor_n_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature_restor_n_avoiddeforest <- SFAge2020_restor_n_avoiddeforest
SF2020mature_restor_n_avoiddeforest[] <- ifelse(SF2020mature_restor_n_avoiddeforest[]>5, 1, 0)
writeRaster(SF2020mature_restor_n_avoiddeforest, "rasters/STM/input/SF2020mature_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)

MF2020_restor_n_avoiddeforest <- sum(UPF2020_restor_n_avoiddeforest, DPF2020_restor_n_avoiddeforest)
MF2020_restor_n_avoiddeforest <- sum(MF2020_restor_n_avoiddeforest, SF2020mature_restor_n_avoiddeforest)
writeRaster(MF2020_restor_n_avoiddeforest, "rasters/STM/input/MF2020_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_restor_n_avoiddeforest.px <- focal(MF2020_restor_n_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_restor_n_avoiddeforest.px)<-"MFpx"
MF2020_restor_n_avoiddeforest.px[is.nan(MF2020_restor_n_avoiddeforest.px)] <- 0
MF2020_restor_n_avoiddeforest.px[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_restor_n_avoiddeforest.px[])
#saving
writeRaster(MF2020_restor_n_avoiddeforest.px, "rasters/STM/2020_restor_n_avoiddeforest/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_restor_n_avoiddeforest.ls <- focal(MF2020_restor_n_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_restor_n_avoiddeforest.ls)<-"MFls"
MF2020_restor_n_avoiddeforest.ls[is.nan(MF2020_restor_n_avoiddeforest.ls)] <- 0
MF2020_restor_n_avoiddeforest.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_restor_n_avoiddeforest.ls[])
#saving
writeRaster(MF2020_restor_n_avoiddeforest.ls, "rasters/STM/2020_restor_n_avoiddeforest/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoiddeforestmature", "SFAge2020_restor_n_avoiddeforest.px", "SFAge2020_restor_n_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_restor_n_avoiddeforest <- raster("rasters/STM/input/MF2020_restor_n_avoiddeforest.tif")
inv.MF2020_restor_n_avoiddeforest <- MF2020_restor_n_avoiddeforest
inv.MF2020_restor_n_avoiddeforest[inv.MF2020_restor_n_avoiddeforest==1]<-NA
#cheking
#inv.MF2020_restor_n_avoiddeforest
#plot(inv.MF2020_restor_n_avoiddeforest)

edge.dist.2020_restor_n_avoiddeforest <- distance(inv.MF2020_restor_n_avoiddeforest, doEdge=T)
names(edge.dist.2020_restor_n_avoiddeforest)<-"edgedist"
edge.dist.2020_restor_n_avoiddeforest[is.nan(edge.dist.2020_restor_n_avoiddeforest)] <- 0
edge.dist.2020_restor_n_avoiddeforest[] <- ifelse(stm.lulc[[1]][]==0, NA, edge.dist.2020_restor_n_avoiddeforest[])
#saving
writeRaster(edge.dist.2020_restor_n_avoiddeforest, "rasters/STM/2020_restor_n_avoiddeforest/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_restor_n_avoiddeforest")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_restor_n_avoiddeforest <- edge.dist.2020_restor_n_avoiddeforest
edge2020_restor_n_avoiddeforest[] <- ifelse(edge2020_restor_n_avoiddeforest[] < 200, 0, ifelse(edge2020_restor_n_avoiddeforest[]>300, 0, 1))
writeRaster(edge2020_restor_n_avoiddeforest, "rasters/STM/input/edge2020_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_restor_n_avoiddeforest.px <- focal(edge2020_restor_n_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_restor_n_avoiddeforest.px)<-"edgepx"
edge2020_restor_n_avoiddeforest.px[is.nan(edge2020_restor_n_avoiddeforest.px)] <- 0
edge2020_restor_n_avoiddeforest.px[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_restor_n_avoiddeforest.px[])
#saving
writeRaster(edge2020_restor_n_avoiddeforest.px, "rasters/STM/2020_restor_n_avoiddeforest/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_restor_n_avoiddeforest.ls <- focal(edge2020_restor_n_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_restor_n_avoiddeforest.ls)<-"edgels"
edge2020_restor_n_avoiddeforest.ls[is.nan(edge2020_restor_n_avoiddeforest.ls)] <- 0
edge2020_restor_n_avoiddeforest.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_restor_n_avoiddeforest.ls[])
#saving
writeRaster(edge2020_restor_n_avoiddeforest.ls, "rasters/STM/2020_restor_n_avoiddeforest/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_restor_n_avoiddeforest.px", "edge2020_restor_n_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Restoration and avoid deforestation (Primary forests only) ========
### Undegraded primary forest
#UPF2020_restor_n_avoiddeforest2 <- raster("rasters/STM/input/LULC/UPF2020_restor_n_avoiddeforest2.tif")

### mean upf cover in local scale (90m)
UPF2020_restor_n_avoiddeforest2.px <- focal(UPF2020_restor_n_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_restor_n_avoiddeforest2.px)<-"UPFpx"
UPF2020_restor_n_avoiddeforest2.px[is.nan(UPF2020_restor_n_avoiddeforest2.px)] <- 0
UPF2020_restor_n_avoiddeforest2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_restor_n_avoiddeforest2.px[])
#saving
writeRaster(UPF2020_restor_n_avoiddeforest2.px, "rasters/STM/2020_restor_n_avoiddeforest2/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_restor_n_avoiddeforest2.ls <- focal(UPF2020_restor_n_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_restor_n_avoiddeforest2.ls)<-"UPFls"
UPF2020_restor_n_avoiddeforest2.ls[is.nan(UPF2020_restor_n_avoiddeforest2.ls)] <- 0
UPF2020_restor_n_avoiddeforest2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_restor_n_avoiddeforest2.ls[])
#saving
writeRaster(UPF2020_restor_n_avoiddeforest2.ls, "rasters/STM/2020_restor_n_avoiddeforest2/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_restor_n_avoiddeforest2.px", "UPF2020_restor_n_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_restor_n_avoiddeforest2 <- raster("rasters/STM/input/LULC/uDPF2020_restor_n_avoiddeforest2.tif")
#RDPF2020_restor_n_avoiddeforest2 <- raster("rasters/STM/input/LULC/RDPF2020_restor_n_avoiddeforest2.tif")
DPF2020_restor_n_avoiddeforest2 <- sum(uDPF2020_restor_n_avoiddeforest2, RDPF2020_restor_n_avoiddeforest2)

### mean dpf cover in local scale (90m)
DPF2020_restor_n_avoiddeforest2.px <- focal(DPF2020_restor_n_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_restor_n_avoiddeforest2.px)<-"DPFpx"
DPF2020_restor_n_avoiddeforest2.px[is.nan(DPF2020_restor_n_avoiddeforest2.px)] <- 0
DPF2020_restor_n_avoiddeforest2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_restor_n_avoiddeforest2.px[])
#saving
writeRaster(DPF2020_restor_n_avoiddeforest2.px, "rasters/STM/2020_restor_n_avoiddeforest2/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_restor_n_avoiddeforest2.ls <- focal(DPF2020_restor_n_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_restor_n_avoiddeforest2.ls)<-"DPFls"
DPF2020_restor_n_avoiddeforest2.ls[is.nan(DPF2020_restor_n_avoiddeforest2.ls)] <- 0
DPF2020_restor_n_avoiddeforest2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_restor_n_avoiddeforest2.ls[])
#saving
writeRaster(DPF2020_restor_n_avoiddeforest2.ls, "rasters/STM/2020_restor_n_avoiddeforest2/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_restor_n_avoiddeforest2.px", "DPF2020_restor_n_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_restor_n_avoiddeforest2 <- TSD2020_avoiddeforest2
writeRaster(TSD2020_restor_n_avoiddeforest2, "rasters/STM/input/TSD2020_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_restor_n_avoiddeforest2.px <- focal(TSD2020_restor_n_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_restor_n_avoiddeforest2.px)<-"TSDpx"
TSD2020_restor_n_avoiddeforest2.px[is.nan(TSD2020_restor_n_avoiddeforest2.px)] <- 0
TSD2020_restor_n_avoiddeforest2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_restor_n_avoiddeforest2.px[])
#saving
writeRaster(TSD2020_restor_n_avoiddeforest2.px, "rasters/STM/2020_restor_n_avoiddeforest2/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_restor_n_avoiddeforest2.ls <- focal(TSD2020_restor_n_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_restor_n_avoiddeforest2.ls)<-"TSDls"
TSD2020_restor_n_avoiddeforest2.ls[is.nan(TSD2020_restor_n_avoiddeforest2.ls)] <- 0
TSD2020_restor_n_avoiddeforest2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_restor_n_avoiddeforest2.ls[])
#saving
writeRaster(TSD2020_restor_n_avoiddeforest2.ls, "rasters/STM/2020_restor_n_avoiddeforest2/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_restor_n_avoiddeforest2.px", "TSD2020_restor_n_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_restor_n_avoiddeforest2 <- raster("rasters/STM/input/LULC/uSF2020_restor_n_avoiddeforest2.tif")
#DSF2020_restor_n_avoiddeforest2 <- raster("rasters/STM/input/LULC/DSF2020_restor_n_avoiddeforest2.tif")
SF2020_restor_n_avoiddeforest2 <- sum(uSF2020_restor_n_avoiddeforest2, DSF2020_restor_n_avoiddeforest2)

### mean sf cover in local scale (90m)
SF2020_restor_n_avoiddeforest2.px <- focal(SF2020_restor_n_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_restor_n_avoiddeforest2.px)<-"SFpx"
SF2020_restor_n_avoiddeforest2.px[is.nan(SF2020_restor_n_avoiddeforest2.px)] <- 0
SF2020_restor_n_avoiddeforest2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_restor_n_avoiddeforest2.px[])
#saving
writeRaster(SF2020_restor_n_avoiddeforest2.px, "rasters/STM/2020_restor_n_avoiddeforest2/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_restor_n_avoiddeforest2.ls <- focal(SF2020_restor_n_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_restor_n_avoiddeforest2.ls)<-"SFls"
SF2020_restor_n_avoiddeforest2.ls[is.nan(SF2020_restor_n_avoiddeforest2.ls)] <- 0
SF2020_restor_n_avoiddeforest2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_restor_n_avoiddeforest2.ls[])
#saving
writeRaster(SF2020_restor_n_avoiddeforest2.ls, "rasters/STM/2020_restor_n_avoiddeforest2/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoiddeforest2.px", "SF2020_restor_n_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_restor_n_avoiddeforest2 <- SFAge2020_avoiddeforest2
SFAge2020_restor_n_avoiddeforest2[] <- ifelse(SFAge2020_restor_n_avoiddeforest2[]==0 & SF2020_restor_n_avoiddeforest2[]==1, 10, SFAge2020_restor_n_avoiddeforest2[])
writeRaster(SFAge2020_restor_n_avoiddeforest2, "rasters/STM/input/SFAge2020_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_restor_n_avoiddeforest2.px <- focal(SFAge2020_restor_n_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_restor_n_avoiddeforest2.px)<-"SFAgepx"
SFAge2020_restor_n_avoiddeforest2.px[is.nan(SFAge2020_restor_n_avoiddeforest2.px)] <- 0
SFAge2020_restor_n_avoiddeforest2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_restor_n_avoiddeforest2.px[])
#saving
writeRaster(SFAge2020_restor_n_avoiddeforest2.px, "rasters/STM/2020_restor_n_avoiddeforest2/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_restor_n_avoiddeforest2.ls <- focal(SFAge2020_restor_n_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_restor_n_avoiddeforest2.ls)<-"SFAgels"
SFAge2020_restor_n_avoiddeforest2.ls[is.nan(SFAge2020_restor_n_avoiddeforest2.ls)] <- 0
SFAge2020_restor_n_avoiddeforest2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_restor_n_avoiddeforest2.ls[])
#saving
writeRaster(SFAge2020_restor_n_avoiddeforest2.ls, "rasters/STM/2020_restor_n_avoiddeforest2/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_restor_n_avoiddeforest2.px", "SFAge2020_restor_n_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young_restor_n_avoiddeforest2 <- SFAge2020_restor_n_avoiddeforest2
SF2020young_restor_n_avoiddeforest2[] <- ifelse(SF2020young_restor_n_avoiddeforest2[]>2, 1, 0)
writeRaster(SF2020young_restor_n_avoiddeforest2, "rasters/STM/input/SF2020young_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)

TF2020_restor_n_avoiddeforest2 <- sum(UPF2020_restor_n_avoiddeforest2, DPF2020_restor_n_avoiddeforest2)
TF2020_restor_n_avoiddeforest2 <- sum(TF2020_restor_n_avoiddeforest2, SF2020young_restor_n_avoiddeforest2)
writeRaster(TF2020_restor_n_avoiddeforest2, "rasters/STM/input/TF2020_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_restor_n_avoiddeforest2.px <- focal(TF2020_restor_n_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_restor_n_avoiddeforest2.px)<-"TFpx"
TF2020_restor_n_avoiddeforest2.px[is.nan(TF2020_restor_n_avoiddeforest2.px)] <- 0
TF2020_restor_n_avoiddeforest2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_restor_n_avoiddeforest2.px[])
#saving
writeRaster(TF2020_restor_n_avoiddeforest2.px, "rasters/STM/2020_restor_n_avoiddeforest2/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_restor_n_avoiddeforest2.ls <- focal(TF2020_restor_n_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_restor_n_avoiddeforest2.ls)<-"TFls"
TF2020_restor_n_avoiddeforest2.ls[is.nan(TF2020_restor_n_avoiddeforest2.ls)] <- 0
TF2020_restor_n_avoiddeforest2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_restor_n_avoiddeforest2.ls[])
#saving
writeRaster(TF2020_restor_n_avoiddeforest2.ls, "rasters/STM/2020_restor_n_avoiddeforest2/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoiddeforest2young", "SFAge2020_restor_n_avoiddeforest2.px", "SFAge2020_restor_n_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature_restor_n_avoiddeforest2 <- SFAge2020_restor_n_avoiddeforest2
SF2020mature_restor_n_avoiddeforest2[] <- ifelse(SF2020mature_restor_n_avoiddeforest2[]>5, 1, 0)
writeRaster(SF2020mature_restor_n_avoiddeforest2, "rasters/STM/input/SF2020mature_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)

MF2020_restor_n_avoiddeforest2 <- sum(UPF2020_restor_n_avoiddeforest2, DPF2020_restor_n_avoiddeforest2)
MF2020_restor_n_avoiddeforest2 <- sum(MF2020_restor_n_avoiddeforest2, SF2020mature_restor_n_avoiddeforest2)
writeRaster(MF2020_restor_n_avoiddeforest2, "rasters/STM/input/MF2020_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_restor_n_avoiddeforest2.px <- focal(MF2020_restor_n_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_restor_n_avoiddeforest2.px)<-"MFpx"
MF2020_restor_n_avoiddeforest2.px[is.nan(MF2020_restor_n_avoiddeforest2.px)] <- 0
MF2020_restor_n_avoiddeforest2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_restor_n_avoiddeforest2.px[])
#saving
writeRaster(MF2020_restor_n_avoiddeforest2.px, "rasters/STM/2020_restor_n_avoiddeforest2/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_restor_n_avoiddeforest2.ls <- focal(MF2020_restor_n_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_restor_n_avoiddeforest2.ls)<-"MFls"
MF2020_restor_n_avoiddeforest2.ls[is.nan(MF2020_restor_n_avoiddeforest2.ls)] <- 0
MF2020_restor_n_avoiddeforest2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_restor_n_avoiddeforest2.ls[])
#saving
writeRaster(MF2020_restor_n_avoiddeforest2.ls, "rasters/STM/2020_restor_n_avoiddeforest2/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoiddeforest2mature", "SFAge2020_restor_n_avoiddeforest2.px", "SFAge2020_restor_n_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_restor_n_avoiddeforest2 <- raster("rasters/STM/input/MF2020_restor_n_avoiddeforest2.tif")
inv.MF2020_restor_n_avoiddeforest2 <- MF2020_restor_n_avoiddeforest2
inv.MF2020_restor_n_avoiddeforest2[inv.MF2020_restor_n_avoiddeforest2==1]<-NA
#cheking
#inv.MF2020_restor_n_avoiddeforest2
#plot(inv.MF2020_restor_n_avoiddeforest2)

edge.dist.2020_restor_n_avoiddeforest2 <- distance(inv.MF2020_restor_n_avoiddeforest2, doEdge=T)
names(edge.dist.2020_restor_n_avoiddeforest2)<-"edgedist"
edge.dist.2020_restor_n_avoiddeforest2[is.nan(edge.dist.2020_restor_n_avoiddeforest2)] <- 0
edge.dist.2020_restor_n_avoiddeforest2[] <- ifelse(stm.lulc[[1]][]==0, NA, edge.dist.2020_restor_n_avoiddeforest2[])
#saving
writeRaster(edge.dist.2020_restor_n_avoiddeforest2, "rasters/STM/2020_restor_n_avoiddeforest2/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_restor_n_avoiddeforest2")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_restor_n_avoiddeforest2 <- edge.dist.2020_restor_n_avoiddeforest2
edge2020_restor_n_avoiddeforest2[] <- ifelse(edge2020_restor_n_avoiddeforest2[] < 200, 0, ifelse(edge2020_restor_n_avoiddeforest2[]>300, 0, 1))
writeRaster(edge2020_restor_n_avoiddeforest2, "rasters/STM/input/edge2020_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_restor_n_avoiddeforest2.px <- focal(edge2020_restor_n_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_restor_n_avoiddeforest2.px)<-"edgepx"
edge2020_restor_n_avoiddeforest2.px[is.nan(edge2020_restor_n_avoiddeforest2.px)] <- 0
edge2020_restor_n_avoiddeforest2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_restor_n_avoiddeforest2.px[])
#saving
writeRaster(edge2020_restor_n_avoiddeforest2.px, "rasters/STM/2020_restor_n_avoiddeforest2/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_restor_n_avoiddeforest2.ls <- focal(edge2020_restor_n_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_restor_n_avoiddeforest2.ls)<-"edgels"
edge2020_restor_n_avoiddeforest2.ls[is.nan(edge2020_restor_n_avoiddeforest2.ls)] <- 0
edge2020_restor_n_avoiddeforest2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_restor_n_avoiddeforest2.ls[])
#saving
writeRaster(edge2020_restor_n_avoiddeforest2.ls, "rasters/STM/2020_restor_n_avoiddeforest2/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_restor_n_avoiddeforest2.px", "edge2020_restor_n_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Restoration and avoid both (all) ==================================
### Undegraded primary forest
#UPF2020_restor_n_avoidboth <- raster("rasters/STM/input/LULC/UPF2020_restor_n_avoidboth.tif")

### mean upf cover in local scale (90m)
UPF2020_restor_n_avoidboth.px <- focal(UPF2020_restor_n_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_restor_n_avoidboth.px)<-"UPFpx"
UPF2020_restor_n_avoidboth.px[is.nan(UPF2020_restor_n_avoidboth.px)] <- 0
UPF2020_restor_n_avoidboth.px[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_restor_n_avoidboth.px[])
#saving
writeRaster(UPF2020_restor_n_avoidboth.px, "rasters/STM/2020_restor_n_avoidboth/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_restor_n_avoidboth.ls <- focal(UPF2020_restor_n_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_restor_n_avoidboth.ls)<-"UPFls"
UPF2020_restor_n_avoidboth.ls[is.nan(UPF2020_restor_n_avoidboth.ls)] <- 0
UPF2020_restor_n_avoidboth.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_restor_n_avoidboth.ls[])
#saving
writeRaster(UPF2020_restor_n_avoidboth.ls, "rasters/STM/2020_restor_n_avoidboth/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_restor_n_avoidboth.px", "UPF2020_restor_n_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_restor_n_avoidboth <- raster("rasters/STM/input/LULC/uDPF2020_restor_n_avoidboth.tif")
#RDPF2020_restor_n_avoidboth <- raster("rasters/STM/input/LULC/RDPF2020_restor_n_avoidboth.tif")
DPF2020_restor_n_avoidboth <- sum(uDPF2020_restor_n_avoidboth, RDPF2020_restor_n_avoidboth)

### mean dpf cover in local scale (90m)
DPF2020_restor_n_avoidboth.px <- focal(DPF2020_restor_n_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_restor_n_avoidboth.px)<-"DPFpx"
DPF2020_restor_n_avoidboth.px[is.nan(DPF2020_restor_n_avoidboth.px)] <- 0
DPF2020_restor_n_avoidboth.px[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_restor_n_avoidboth.px[])
#saving
writeRaster(DPF2020_restor_n_avoidboth.px, "rasters/STM/2020_restor_n_avoidboth/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_restor_n_avoidboth.ls <- focal(DPF2020_restor_n_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_restor_n_avoidboth.ls)<-"DPFls"
DPF2020_restor_n_avoidboth.ls[is.nan(DPF2020_restor_n_avoidboth.ls)] <- 0
DPF2020_restor_n_avoidboth.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_restor_n_avoidboth.ls[])
#saving
writeRaster(DPF2020_restor_n_avoidboth.ls, "rasters/STM/2020_restor_n_avoidboth/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_restor_n_avoidboth.px", "DPF2020_restor_n_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_restor_n_avoidboth <- TSD2020_avoidboth
writeRaster(TSD2020_restor_n_avoidboth, "rasters/STM/input/TSD2020_restor_n_avoidboth.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_restor_n_avoidboth.px <- focal(TSD2020_restor_n_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_restor_n_avoidboth.px)<-"TSDpx"
TSD2020_restor_n_avoidboth.px[is.nan(TSD2020_restor_n_avoidboth.px)] <- 0
TSD2020_restor_n_avoidboth.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_restor_n_avoidboth.px[])
#saving
writeRaster(TSD2020_restor_n_avoidboth.px, "rasters/STM/2020_restor_n_avoidboth/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_restor_n_avoidboth.ls <- focal(TSD2020_restor_n_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_restor_n_avoidboth.ls)<-"TSDls"
TSD2020_restor_n_avoidboth.ls[is.nan(TSD2020_restor_n_avoidboth.ls)] <- 0
TSD2020_restor_n_avoidboth.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_restor_n_avoidboth.ls[])
#saving
writeRaster(TSD2020_restor_n_avoidboth.ls, "rasters/STM/2020_restor_n_avoidboth/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_restor_n_avoidboth.px", "TSD2020_restor_n_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_restor_n_avoidboth <- raster("rasters/STM/input/LULC/uSF2020_restor_n_avoidboth.tif")
#DSF2020_restor_n_avoidboth <- raster("rasters/STM/input/LULC/DSF2020_restor_n_avoidboth.tif")
SF2020_restor_n_avoidboth <- sum(uSF2020_restor_n_avoidboth, DSF2020_restor_n_avoidboth)

### mean sf cover in local scale (90m)
SF2020_restor_n_avoidboth.px <- focal(SF2020_restor_n_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_restor_n_avoidboth.px)<-"SFpx"
SF2020_restor_n_avoidboth.px[is.nan(SF2020_restor_n_avoidboth.px)] <- 0
SF2020_restor_n_avoidboth.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_restor_n_avoidboth.px[])
#saving
writeRaster(SF2020_restor_n_avoidboth.px, "rasters/STM/2020_restor_n_avoidboth/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_restor_n_avoidboth.ls <- focal(SF2020_restor_n_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_restor_n_avoidboth.ls)<-"SFls"
SF2020_restor_n_avoidboth.ls[is.nan(SF2020_restor_n_avoidboth.ls)] <- 0
SF2020_restor_n_avoidboth.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_restor_n_avoidboth.ls[])
#saving
writeRaster(SF2020_restor_n_avoidboth.ls, "rasters/STM/2020_restor_n_avoidboth/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoidboth.px", "SF2020_restor_n_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_restor_n_avoidboth <- SFAge2020_avoidboth
SFAge2020_restor_n_avoidboth[] <- ifelse(SFAge2020_restor_n_avoidboth[]==0 & SF2020_restor_n_avoidboth[]==1, 10, SFAge2020_restor_n_avoidboth[])
writeRaster(SFAge2020_restor_n_avoidboth, "rasters/STM/input/SFAge2020_restor_n_avoidboth.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_restor_n_avoidboth.px <- focal(SFAge2020_restor_n_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_restor_n_avoidboth.px)<-"SFAgepx"
SFAge2020_restor_n_avoidboth.px[is.nan(SFAge2020_restor_n_avoidboth.px)] <- 0
SFAge2020_restor_n_avoidboth.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_restor_n_avoidboth.px[])
#saving
writeRaster(SFAge2020_restor_n_avoidboth.px, "rasters/STM/2020_restor_n_avoidboth/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_restor_n_avoidboth.ls <- focal(SFAge2020_restor_n_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_restor_n_avoidboth.ls)<-"SFAgels"
SFAge2020_restor_n_avoidboth.ls[is.nan(SFAge2020_restor_n_avoidboth.ls)] <- 0
SFAge2020_restor_n_avoidboth.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_restor_n_avoidboth.ls[])
#saving
writeRaster(SFAge2020_restor_n_avoidboth.ls, "rasters/STM/2020_restor_n_avoidboth/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_restor_n_avoidboth.px", "SFAge2020_restor_n_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young_restor_n_avoidboth <- SFAge2020_restor_n_avoidboth
SF2020young_restor_n_avoidboth[] <- ifelse(SF2020young_restor_n_avoidboth[]>2, 1, 0)
writeRaster(SF2020young_restor_n_avoidboth, "rasters/STM/input/SF2020young_restor_n_avoidboth.tif", format="GTiff", overwrite=T)

TF2020_restor_n_avoidboth <- sum(UPF2020_restor_n_avoidboth, DPF2020_restor_n_avoidboth)
TF2020_restor_n_avoidboth <- sum(TF2020_restor_n_avoidboth, SF2020young_restor_n_avoidboth)
writeRaster(TF2020_restor_n_avoidboth, "rasters/STM/input/TF2020_restor_n_avoidboth.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_restor_n_avoidboth.px <- focal(TF2020_restor_n_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_restor_n_avoidboth.px)<-"TFpx"
TF2020_restor_n_avoidboth.px[is.nan(TF2020_restor_n_avoidboth.px)] <- 0
TF2020_restor_n_avoidboth.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_restor_n_avoidboth.px[])
#saving
writeRaster(TF2020_restor_n_avoidboth.px, "rasters/STM/2020_restor_n_avoidboth/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_restor_n_avoidboth.ls <- focal(TF2020_restor_n_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_restor_n_avoidboth.ls)<-"TFls"
TF2020_restor_n_avoidboth.ls[is.nan(TF2020_restor_n_avoidboth.ls)] <- 0
TF2020_restor_n_avoidboth.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_restor_n_avoidboth.ls[])
#saving
writeRaster(TF2020_restor_n_avoidboth.ls, "rasters/STM/2020_restor_n_avoidboth/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoidbothyoung", "SFAge2020_restor_n_avoidboth.px", "SFAge2020_restor_n_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature_restor_n_avoidboth <- SFAge2020_restor_n_avoidboth
SF2020mature_restor_n_avoidboth[] <- ifelse(SF2020mature_restor_n_avoidboth[]>5, 1, 0)
writeRaster(SF2020mature_restor_n_avoidboth, "rasters/STM/input/SF2020mature_restor_n_avoidboth.tif", format="GTiff", overwrite=T)

MF2020_restor_n_avoidboth <- sum(UPF2020_restor_n_avoidboth, DPF2020_restor_n_avoidboth)
MF2020_restor_n_avoidboth <- sum(MF2020_restor_n_avoidboth, SF2020mature_restor_n_avoidboth)
writeRaster(MF2020_restor_n_avoidboth, "rasters/STM/input/MF2020_restor_n_avoidboth.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_restor_n_avoidboth.px <- focal(MF2020_restor_n_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_restor_n_avoidboth.px)<-"MFpx"
MF2020_restor_n_avoidboth.px[is.nan(MF2020_restor_n_avoidboth.px)] <- 0
MF2020_restor_n_avoidboth.px[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_restor_n_avoidboth.px[])
#saving
writeRaster(MF2020_restor_n_avoidboth.px, "rasters/STM/2020_restor_n_avoidboth/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_restor_n_avoidboth.ls <- focal(MF2020_restor_n_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_restor_n_avoidboth.ls)<-"MFls"
MF2020_restor_n_avoidboth.ls[is.nan(MF2020_restor_n_avoidboth.ls)] <- 0
MF2020_restor_n_avoidboth.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_restor_n_avoidboth.ls[])
#saving
writeRaster(MF2020_restor_n_avoidboth.ls, "rasters/STM/2020_restor_n_avoidboth/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoidbothmature", "SFAge2020_restor_n_avoidboth.px", "SFAge2020_restor_n_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_restor_n_avoidboth <- raster("rasters/STM/input/MF2020_restor_n_avoidboth.tif")
inv.MF2020_restor_n_avoidboth <- MF2020_restor_n_avoidboth
inv.MF2020_restor_n_avoidboth[inv.MF2020_restor_n_avoidboth==1]<-NA
#cheking
#inv.MF2020_restor_n_avoidboth
#plot(inv.MF2020_restor_n_avoidboth)

edge.dist.2020_restor_n_avoidboth <- distance(inv.MF2020_restor_n_avoidboth, doEdge=T)
names(edge.dist.2020_restor_n_avoidboth)<-"edgedist"
edge.dist.2020_restor_n_avoidboth[is.nan(edge.dist.2020_restor_n_avoidboth)] <- 0
edge.dist.2020_restor_n_avoidboth[] <- ifelse(stm.lulc[[1]][]==0, NA, edge.dist.2020_restor_n_avoidboth[])
#saving
writeRaster(edge.dist.2020_restor_n_avoidboth, "rasters/STM/2020_restor_n_avoidboth/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_restor_n_avoidboth")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_restor_n_avoidboth <- edge.dist.2020_restor_n_avoidboth
edge2020_restor_n_avoidboth[] <- ifelse(edge2020_restor_n_avoidboth[] < 200, 0, ifelse(edge2020_restor_n_avoidboth[]>300, 0, 1))
writeRaster(edge2020_restor_n_avoidboth, "rasters/STM/input/edge2020_restor_n_avoidboth.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_restor_n_avoidboth.px <- focal(edge2020_restor_n_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_restor_n_avoidboth.px)<-"edgepx"
edge2020_restor_n_avoidboth.px[is.nan(edge2020_restor_n_avoidboth.px)] <- 0
edge2020_restor_n_avoidboth.px[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_restor_n_avoidboth.px[])
#saving
writeRaster(edge2020_restor_n_avoidboth.px, "rasters/STM/2020_restor_n_avoidboth/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_restor_n_avoidboth.ls <- focal(edge2020_restor_n_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_restor_n_avoidboth.ls)<-"edgels"
edge2020_restor_n_avoidboth.ls[is.nan(edge2020_restor_n_avoidboth.ls)] <- 0
edge2020_restor_n_avoidboth.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_restor_n_avoidboth.ls[])
#saving
writeRaster(edge2020_restor_n_avoidboth.ls, "rasters/STM/2020_restor_n_avoidboth/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_restor_n_avoidboth.px", "edge2020_restor_n_avoidboth.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Restoration and avoid both (Primary forests only) =================
### Undegraded primary forest
#UPF2020_restor_n_avoidboth2 <- raster("rasters/STM/input/LULC/UPF2020_restor_n_avoidboth2.tif")

### mean upf cover in local scale (90m)
UPF2020_restor_n_avoidboth2.px <- focal(UPF2020_restor_n_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_restor_n_avoidboth2.px)<-"UPFpx"
UPF2020_restor_n_avoidboth2.px[is.nan(UPF2020_restor_n_avoidboth2.px)] <- 0
UPF2020_restor_n_avoidboth2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_restor_n_avoidboth2.px[])
#saving
writeRaster(UPF2020_restor_n_avoidboth2.px, "rasters/STM/2020_restor_n_avoidboth2/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_restor_n_avoidboth2.ls <- focal(UPF2020_restor_n_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_restor_n_avoidboth2.ls)<-"UPFls"
UPF2020_restor_n_avoidboth2.ls[is.nan(UPF2020_restor_n_avoidboth2.ls)] <- 0
UPF2020_restor_n_avoidboth2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, UPF2020_restor_n_avoidboth2.ls[])
#saving
writeRaster(UPF2020_restor_n_avoidboth2.ls, "rasters/STM/2020_restor_n_avoidboth2/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_restor_n_avoidboth2.px", "UPF2020_restor_n_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_restor_n_avoidboth2 <- raster("rasters/STM/input/LULC/uDPF2020_restor_n_avoidboth2.tif")
#RDPF2020_restor_n_avoidboth2 <- raster("rasters/STM/input/LULC/RDPF2020_restor_n_avoidboth2.tif")
DPF2020_restor_n_avoidboth2 <- sum(uDPF2020_restor_n_avoidboth2, RDPF2020_restor_n_avoidboth2)

### mean dpf cover in local scale (90m)
DPF2020_restor_n_avoidboth2.px <- focal(DPF2020_restor_n_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_restor_n_avoidboth2.px)<-"DPFpx"
DPF2020_restor_n_avoidboth2.px[is.nan(DPF2020_restor_n_avoidboth2.px)] <- 0
DPF2020_restor_n_avoidboth2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_restor_n_avoidboth2.px[])
#saving
writeRaster(DPF2020_restor_n_avoidboth2.px, "rasters/STM/2020_restor_n_avoidboth2/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_restor_n_avoidboth2.ls <- focal(DPF2020_restor_n_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_restor_n_avoidboth2.ls)<-"DPFls"
DPF2020_restor_n_avoidboth2.ls[is.nan(DPF2020_restor_n_avoidboth2.ls)] <- 0
DPF2020_restor_n_avoidboth2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, DPF2020_restor_n_avoidboth2.ls[])
#saving
writeRaster(DPF2020_restor_n_avoidboth2.ls, "rasters/STM/2020_restor_n_avoidboth2/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_restor_n_avoidboth2.px", "DPF2020_restor_n_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_restor_n_avoidboth2 <- TSD2020_avoidboth2
writeRaster(TSD2020_restor_n_avoidboth2, "rasters/STM/input/TSD2020_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_restor_n_avoidboth2.px <- focal(TSD2020_restor_n_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_restor_n_avoidboth2.px)<-"TSDpx"
TSD2020_restor_n_avoidboth2.px[is.nan(TSD2020_restor_n_avoidboth2.px)] <- 0
TSD2020_restor_n_avoidboth2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_restor_n_avoidboth2.px[])
#saving
writeRaster(TSD2020_restor_n_avoidboth2.px, "rasters/STM/2020_restor_n_avoidboth2/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_restor_n_avoidboth2.ls <- focal(TSD2020_restor_n_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_restor_n_avoidboth2.ls)<-"TSDls"
TSD2020_restor_n_avoidboth2.ls[is.nan(TSD2020_restor_n_avoidboth2.ls)] <- 0
TSD2020_restor_n_avoidboth2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TSD2020_restor_n_avoidboth2.ls[])
#saving
writeRaster(TSD2020_restor_n_avoidboth2.ls, "rasters/STM/2020_restor_n_avoidboth2/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_restor_n_avoidboth2.px", "TSD2020_restor_n_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_restor_n_avoidboth2 <- raster("rasters/STM/input/LULC/uSF2020_restor_n_avoidboth2.tif")
#DSF2020_restor_n_avoidboth2 <- raster("rasters/STM/input/LULC/DSF2020_restor_n_avoidboth2.tif")
SF2020_restor_n_avoidboth2 <- sum(uSF2020_restor_n_avoidboth2, DSF2020_restor_n_avoidboth2)

### mean sf cover in local scale (90m)
SF2020_restor_n_avoidboth2.px <- focal(SF2020_restor_n_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_restor_n_avoidboth2.px)<-"SFpx"
SF2020_restor_n_avoidboth2.px[is.nan(SF2020_restor_n_avoidboth2.px)] <- 0
SF2020_restor_n_avoidboth2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_restor_n_avoidboth2.px[])
#saving
writeRaster(SF2020_restor_n_avoidboth2.px, "rasters/STM/2020_restor_n_avoidboth2/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_restor_n_avoidboth2.ls <- focal(SF2020_restor_n_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_restor_n_avoidboth2.ls)<-"SFls"
SF2020_restor_n_avoidboth2.ls[is.nan(SF2020_restor_n_avoidboth2.ls)] <- 0
SF2020_restor_n_avoidboth2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SF2020_restor_n_avoidboth2.ls[])
#saving
writeRaster(SF2020_restor_n_avoidboth2.ls, "rasters/STM/2020_restor_n_avoidboth2/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoidboth2.px", "SF2020_restor_n_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_restor_n_avoidboth2 <- SFAge2020_avoidboth2
SFAge2020_restor_n_avoidboth2[] <- ifelse(SFAge2020_restor_n_avoidboth2[]==0 & SF2020_restor_n_avoidboth2[]==1, 10, SFAge2020_restor_n_avoidboth2[])
writeRaster(SFAge2020_restor_n_avoidboth2, "rasters/STM/input/SFAge2020_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_restor_n_avoidboth2.px <- focal(SFAge2020_restor_n_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_restor_n_avoidboth2.px)<-"SFAgepx"
SFAge2020_restor_n_avoidboth2.px[is.nan(SFAge2020_restor_n_avoidboth2.px)] <- 0
SFAge2020_restor_n_avoidboth2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_restor_n_avoidboth2.px[])
#saving
writeRaster(SFAge2020_restor_n_avoidboth2.px, "rasters/STM/2020_restor_n_avoidboth2/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_restor_n_avoidboth2.ls <- focal(SFAge2020_restor_n_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_restor_n_avoidboth2.ls)<-"SFAgels"
SFAge2020_restor_n_avoidboth2.ls[is.nan(SFAge2020_restor_n_avoidboth2.ls)] <- 0
SFAge2020_restor_n_avoidboth2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, SFAge2020_restor_n_avoidboth2.ls[])
#saving
writeRaster(SFAge2020_restor_n_avoidboth2.ls, "rasters/STM/2020_restor_n_avoidboth2/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_restor_n_avoidboth2.px", "SFAge2020_restor_n_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young_restor_n_avoidboth2 <- SFAge2020_restor_n_avoidboth2
SF2020young_restor_n_avoidboth2[] <- ifelse(SF2020young_restor_n_avoidboth2[]>2, 1, 0)
writeRaster(SF2020young_restor_n_avoidboth2, "rasters/STM/input/SF2020young_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)

TF2020_restor_n_avoidboth2 <- sum(UPF2020_restor_n_avoidboth2, DPF2020_restor_n_avoidboth2)
TF2020_restor_n_avoidboth2 <- sum(TF2020_restor_n_avoidboth2, SF2020young_restor_n_avoidboth2)
writeRaster(TF2020_restor_n_avoidboth2, "rasters/STM/input/TF2020_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_restor_n_avoidboth2.px <- focal(TF2020_restor_n_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_restor_n_avoidboth2.px)<-"TFpx"
TF2020_restor_n_avoidboth2.px[is.nan(TF2020_restor_n_avoidboth2.px)] <- 0
TF2020_restor_n_avoidboth2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_restor_n_avoidboth2.px[])
#saving
writeRaster(TF2020_restor_n_avoidboth2.px, "rasters/STM/2020_restor_n_avoidboth2/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_restor_n_avoidboth2.ls <- focal(TF2020_restor_n_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_restor_n_avoidboth2.ls)<-"TFls"
TF2020_restor_n_avoidboth2.ls[is.nan(TF2020_restor_n_avoidboth2.ls)] <- 0
TF2020_restor_n_avoidboth2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, TF2020_restor_n_avoidboth2.ls[])
#saving
writeRaster(TF2020_restor_n_avoidboth2.ls, "rasters/STM/2020_restor_n_avoidboth2/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoidboth2young", "SFAge2020_restor_n_avoidboth2.px", "SFAge2020_restor_n_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature_restor_n_avoidboth2 <- SFAge2020_restor_n_avoidboth2
SF2020mature_restor_n_avoidboth2[] <- ifelse(SF2020mature_restor_n_avoidboth2[]>5, 1, 0)
writeRaster(SF2020mature_restor_n_avoidboth2, "rasters/STM/input/SF2020mature_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)

MF2020_restor_n_avoidboth2 <- sum(UPF2020_restor_n_avoidboth2, DPF2020_restor_n_avoidboth2)
MF2020_restor_n_avoidboth2 <- sum(MF2020_restor_n_avoidboth2, SF2020mature_restor_n_avoidboth2)
writeRaster(MF2020_restor_n_avoidboth2, "rasters/STM/input/MF2020_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_restor_n_avoidboth2.px <- focal(MF2020_restor_n_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_restor_n_avoidboth2.px)<-"MFpx"
MF2020_restor_n_avoidboth2.px[is.nan(MF2020_restor_n_avoidboth2.px)] <- 0
MF2020_restor_n_avoidboth2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_restor_n_avoidboth2.px[])
#saving
writeRaster(MF2020_restor_n_avoidboth2.px, "rasters/STM/2020_restor_n_avoidboth2/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_restor_n_avoidboth2.ls <- focal(MF2020_restor_n_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_restor_n_avoidboth2.ls)<-"MFls"
MF2020_restor_n_avoidboth2.ls[is.nan(MF2020_restor_n_avoidboth2.ls)] <- 0
MF2020_restor_n_avoidboth2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, MF2020_restor_n_avoidboth2.ls[])
#saving
writeRaster(MF2020_restor_n_avoidboth2.ls, "rasters/STM/2020_restor_n_avoidboth2/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoidboth2mature", "SFAge2020_restor_n_avoidboth2.px", "SFAge2020_restor_n_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_restor_n_avoidboth2 <- raster("rasters/STM/input/MF2020_restor_n_avoidboth2.tif")
inv.MF2020_restor_n_avoidboth2 <- MF2020_restor_n_avoidboth2
inv.MF2020_restor_n_avoidboth2[inv.MF2020_restor_n_avoidboth2==1]<-NA
#cheking
#inv.MF2020_restor_n_avoidboth2
#plot(inv.MF2020_restor_n_avoidboth2)

edge.dist.2020_restor_n_avoidboth2 <- distance(inv.MF2020_restor_n_avoidboth2, doEdge=T)
names(edge.dist.2020_restor_n_avoidboth2)<-"edgedist"
edge.dist.2020_restor_n_avoidboth2[is.nan(edge.dist.2020_restor_n_avoidboth2)] <- 0
edge.dist.2020_restor_n_avoidboth2[] <- ifelse(stm.lulc[[1]][]==0, NA, edge.dist.2020_restor_n_avoidboth2[])
#saving
writeRaster(edge.dist.2020_restor_n_avoidboth2, "rasters/STM/2020_restor_n_avoidboth2/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_restor_n_avoidboth2")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_restor_n_avoidboth2 <- edge.dist.2020_restor_n_avoidboth2
edge2020_restor_n_avoidboth2[] <- ifelse(edge2020_restor_n_avoidboth2[] < 200, 0, ifelse(edge2020_restor_n_avoidboth2[]>300, 0, 1))
writeRaster(edge2020_restor_n_avoidboth2, "rasters/STM/input/edge2020_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_restor_n_avoidboth2.px <- focal(edge2020_restor_n_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_restor_n_avoidboth2.px)<-"edgepx"
edge2020_restor_n_avoidboth2.px[is.nan(edge2020_restor_n_avoidboth2.px)] <- 0
edge2020_restor_n_avoidboth2.px[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_restor_n_avoidboth2.px[])
#saving
writeRaster(edge2020_restor_n_avoidboth2.px, "rasters/STM/2020_restor_n_avoidboth2/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_restor_n_avoidboth2.ls <- focal(edge2020_restor_n_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_restor_n_avoidboth2.ls)<-"edgels"
edge2020_restor_n_avoidboth2.ls[is.nan(edge2020_restor_n_avoidboth2.ls)] <- 0
edge2020_restor_n_avoidboth2.ls[] <- ifelse(stm.lulc[[1]][]==0, NA, edge2020_restor_n_avoidboth2.ls[])
#saving
writeRaster(edge2020_restor_n_avoidboth2.ls, "rasters/STM/2020_restor_n_avoidboth2/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_restor_n_avoidboth2.px", "edge2020_restor_n_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#



## Fixed costs =================================================================
## avoid degradation costs -- adding fire control costs
#' creating and maintaining (i.e., every four years on average) 
#' 6m wide fire breaks at the edge of forested areas

#calculating forest perimeter in each property
#adding variable for forest cover
stm.car@data$FOREST_PERIMETER <- NA
stm.car@data$ncell <- NA

j=nrow(stm.car@data)
for (i in stm.car$COD_IMOVEL) {
  
  rural.property <- stm.car[stm.car$COD_IMOVEL==i,]
  rural.property.edge <- crop(stm.deforest.dist.copy, extent(rural.property))
  rural.property.edge <- mask(rural.property.edge, rural.property)
  rural.property.edge[rural.property.edge > 100]<-NA
  rural.property.edge[rural.property.edge < 1]<-NA
  rural.property.edge[rural.property.edge<=100]<-1
  stm.car@data[stm.car$COD_IMOVEL==i,"ncell"] <- ncell(rural.property.edge)
  if(all(is.na(values(rural.property.edge)))) next
  #convert raster to polygons
  rural.property.edge.shp <- as_Spatial(st_as_sf(st_as_stars(rural.property.edge),
                                                 as_points = FALSE, merge = TRUE))
  
  #cheking & adjustments
  #st_crs(forest.cover.shp)==st_crs(stm.shp)
  #gIsValid(forest.cover.shp)
  #FALSE here means that you'll need to run the buffer routine:
  #forest.cover.shp <- rgeos::gBuffer(forest.cover.shp, byid = TRUE, width = 0)
  
  #estimating the perimeter
  stm.car@data[stm.car$COD_IMOVEL==i,"FOREST_PERIMETER"] <- sum(st_length(st_cast(st_as_sf(rural.property.edge.shp),"MULTILINESTRING")), na.rm=T)/2
  
  
  j=j-1
  cat("\n>", j, "out of", nrow(stm.car@data), "properties left<\n")
  
}


#' @description fire breaks could be cleared at rate of 33.333 meters per day, costing R$100 per day according to IPAM
#' source: https://www.terrabrasilis.org.br/ecotecadigital/pdf/tecnicas-de-prevencao-de-fogo-acidental-metodo-bom-manejo-de-fogo-para-areas-de-agricultura-familiar.pdf
#' cost of fire control was (P/33.33333) x 100, where P is the perimeter in meters of forested area in the property
stm.car@data$cost <- (as.numeric((stm.car@data$FOREST_PERIMETER/33.33333) * 100))/stm.car@data$ncell

#convert to raster
avoid.degrad.cost <- rasterize(stm.car, stm.lulc.2010.forest.class, field = "cost", fun = mean)
avoid.degrad.cost[is.na(avoid.degrad.cost)] <- 0
avoid.degrad.cost <- mask(avoid.degrad.cost, stm.shp)
#plot(avoid.degrad.cost)

#saving
writeRaster(avoid.degrad.cost, "models.output/opportunity.costs/STM_2010_real_base_firecontrol.tif", format="GTiff", overwrite=T)



#
#




## passive restoration cost: average per hectare cost associated with:
#' natural regeneration without fence
#' 48.87 ± 0.7 US$/ha -- US$1.00 ~ R$3.87 -- R$189.13
#' natural regeneration with fence
#' 344.07 ± 156 US$/ha -- US$1.00 ~ R$3.87 -- R$1331.55
#' active restoration
#' 2041.27 ± 728 US$/ha -- US$1.00 ~ R$3.87 -- R$7899.71
#' source: Bracalion et al. 2019 https://doi.org/10.1016/j.biocon.2019.108274


#natural regeneration without fence
agro.class <- c(39,41,48)

restor.cost1 <- stm.lulc[["stm.lulc.2010real"]]

values(restor.cost1)[values(restor.cost1) %in% agro.class] <- 1
values(restor.cost1)[values(restor.cost1) > 1] <- 0
names(restor.cost1) <- "restoration.no.fences"


restor.cost1.deforest.dist <- deforest.dist.copy
values(restor.cost1.deforest.dist)[values(restor.cost1.deforest.dist) == 0] <- NA
values(restor.cost1.deforest.dist)[values(restor.cost1.deforest.dist) > 500] <- NA
values(restor.cost1.deforest.dist)[values(restor.cost1.deforest.dist) <= 500] <- 1

restor.cost1 <- mask(restor.cost1, restor.cost1.deforest.dist)
restor.cost1[is.na(restor.cost1)] <- 0
restor.cost1[restor.cost1==1] <- 189.13

#natural regeneration with fence
restor.cost2 <- stm.lulc[["stm.lulc.2010real"]]

restor.cost2[restor.cost2 == 15] <- 1
restor.cost2[restor.cost2 > 1] <- 0
names(restor.cost2) <- "restoration.fences"


restor.cost2 <- mask(restor.cost2, restor.cost1.deforest.dist)
restor.cost2[is.na(restor.cost2)] <- 0
restor.cost2[restor.cost2==1] <- 1331.55

#active restoration
restor.cost3 <- stm.lulc[["stm.lulc.2010real"]]

values(restor.cost3)[values(restor.cost3) %in% deforestation.class.list] <- 1
values(restor.cost3)[values(restor.cost3) > 1] <- 0
names(restor.cost3) <- "restoration.active"


restor.cost3.deforest.dist <- deforest.dist.copy
values(restor.cost3.deforest.dist)[values(restor.cost3.deforest.dist) <= 500] <- NA
values(restor.cost3.deforest.dist)[values(restor.cost3.deforest.dist) > 500] <- 1

restor.cost3 <- mask(restor.cost3, restor.cost3.deforest.dist)
restor.cost3[is.na(restor.cost3)] <- 0
restor.cost3[restor.cost3==1] <- 7899.71

#restoration cost layer
restor.cost.final <- sum(restor.cost1, restor.cost2, restor.cost3)
restor.cost.final <- mask(restor.cost.final, candidate.areas.final)

writeRaster(restor.cost.final, paste0("models.output/opportunity.costs/STM_2010_real_base_passiverestoration.tif"), format = "GTiff", overwrite = T)



rm(list=ls()[!ls() %in% c("...")]) #keeping only raster stack
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
#




#=================================|end