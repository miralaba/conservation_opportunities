
#' @title Cost-benefit of conservation actions in Amazon
#' @description script to build exploratory variables from 
#' land use - land cover, secondary forest, edge, 
#' degradation (logging and fire), temperature, precipitation, 
#' elevation and distances to road and water body
#' in Paragominas municipality - PA;
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
dir.create("rasters/PGM/input/LULC", recursive = T)



# importing raw rasters ========================================================
## standard projection
std.proj <- "+proj=longlat +datum=WGS84 +units=m +no_defs"

## shapefile paragominas
pgm.shp <- readOGR(dsn = "shapes", layer = "Paragominas_Mask_R3")
proj4string(pgm.shp) <- CRS(std.proj)
pgm.shp <- spTransform(pgm.shp, crs(std.proj))
pgm.shp2 <- gBuffer(pgm.shp, width = 0.02700003)

# 100m resolution raster
#pgm.lulc.100 <- raster("rasters/PGM/raw/pgm-lulc-mapbiomas-brazil-collection-80-2022-100res.tif")
#
#



## land use land cover from mapbiomas collection 7, 30m res [2010 and 2020]

pgm.lulc <- stack(c("rasters/PGM/raw/pgm-lulc-mapbiomas-brazil-collection-80-2010.tif",
                    "rasters/PGM/raw/pgm-lulc-mapbiomas-brazil-collection-80-2020.tif"))
names(pgm.lulc) <- c("pgm.lulc.2010real", "pgm.lulc.2020real")

# resample to 1ha resolution
#pgm.lulc <- resample(pgm.lulc, pgm.lulc.100, method='ngb')
#checking
#st_crs(pgm.lulc)==st_crs(pgm.shp)
#sort(unique(values(pgm.lulc[["pgm.lulc.2020real"]])))

# land ues land cover pixel values and codes
# 0  == NA
# 3  == Forest Formation      == Forest
# 4  == Savanna Formation     == Forest
# 6  == Floodable Forest      == Forest
# 9  == Forest Plantation     == Farming
# 11 == Wetland               == Non Forest Natural Formation
# 12 == Grassland             == Non Forest Natural Formation
# 15 == Pasture               == Farming
# 24 == Urban Area            == Non vegetated area
# 30 == Mining                == Non vegetated area
# 33 == River, Lake and Ocean == Water
# 35 == Palm Oil              == Farming
# 39 == Soybean               == Farming
# 41 == Other Temporary Crops == Farming

#pgm.lulc.2020.df <- as.data.frame(pgm.lulc[["pgm.lulc.2020real"]], xy = TRUE)
#breakpoints <- sort(unique(pgm.lulc.2020.df$pgm.lulc.2020real))
#labels.legend <- c("Non Observed", "Forest Formation", "Savanna Formation", 
#                   "Floodable Forest", "Forest Plantation", "Wetland", "Grassland",
#                   "Pasture", "Urban Area", "Mining", "Water", "Palm Oil", 
#                   "Soybean", "Other Temporary Crops")
#mapbiomas.legend <- c("#ffffff", "#1f8d49", "#7dc975", 
#                      "#026975", "#7a5900", "#519799", "#d6bc74", 
#                      "#edde8e", "#d4271e", "#9c0027", "#2532e4", "#9065d0", 
#                      "#f5b3c8", "#f54ca9")
#
#ggplot() +
#  geom_raster(data = pgm.lulc.2020.df , 
#               aes(x = x, y = y, fill = factor(pgm.lulc.2020real))) + 
#  scale_fill_manual(breaks = breakpoints, 
#                    values = mapbiomas.legend, 
#                    labels = labels.legend, 
#                    name = "LULC Classes") +
#  theme_void()


### isolating forest class pixels
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



### isolating deforestation class pixels (crops, pasture)
deforestation.class.list <- c(15,39,41,48)

pgm.lulc.2010.deforestation.class <- pgm.lulc[["pgm.lulc.2010real"]]
pgm.lulc.2010.deforestation.class[pgm.lulc.2010.deforestation.class %in% deforestation.class.list] <- 1
pgm.lulc.2010.deforestation.class[pgm.lulc.2010.deforestation.class>1] <- 0


pgm.lulc.2010.deforestation.mask <- pgm.lulc.2010.deforestation.class
pgm.lulc.2010.deforestation.mask[pgm.lulc.2010.deforestation.mask==0] <- NA


pgm.lulc.2020.deforestation.class <- pgm.lulc[["pgm.lulc.2020real"]]
pgm.lulc.2020.deforestation.class[pgm.lulc.2020.deforestation.class %in% deforestation.class.list] <- 1
pgm.lulc.2020.deforestation.class[pgm.lulc.2020.deforestation.class>1] <- 0


pgm.lulc.2020.deforestation.mask <- pgm.lulc.2020.deforestation.class
pgm.lulc.2020.deforestation.mask[pgm.lulc.2020.deforestation.mask==0] <- NA



### select pixels based on proximity to forest
pgm.deforest <- pgm.lulc[["pgm.lulc.2010real"]]
pgm.deforest[pgm.deforest %in% deforestation.class.list] <- NA
pgm.deforest[pgm.deforest>1] <- 1

pgm.deforest.dist <- distance(pgm.deforest, doEdge=T)
#writeRaster(pgm.deforest.dist, "rasters/PGM/raw/pgm-distance-to-forest.tif", format = "GTiff", overwrite = T)
#pgm.deforest.dist <- raster("rasters/PGM/raw/pgm-distance-to-forest.tif")

pgm.deforest.dist.copy <- pgm.deforest.dist
#values(pgm.deforest.dist)[values(pgm.deforest.dist) == 0] <- NA
#values(pgm.deforest.dist)[values(pgm.deforest.dist) > 1000] <- NA
#values(pgm.deforest.dist)[values(pgm.deforest.dist) <= 1000] <- 1

rm(pgm.deforest)
gc()
#
#





## time since degradation
#' @description mapbiomas fire is a time-series data from 1985 to present
#' degrad is the first monitoring system to detect degradation in brazil
#' it was functional from 2007 to 2016 but does not differentiate between classes of degradation.
#' deter is the current monitoring system, with data from 2016
#' see auxiliar.R script for details about combining these data sources

pgm.degrad.pf <- stack(c("rasters/PGM/raw/pgm-2010-dpf-tsince0.tif",
                         "rasters/PGM/raw/pgm-2020-dpf-tsince0.tif"))
names(pgm.degrad.pf) <- c("pgm.degrad.2010real", "pgm.degrad.2020real")


#checking
#st_crs(pgm.degrad.pf)==st_crs(pgm.shp)
#plot(pgm.degrad.pf)
#range(values(pgm.degrad.pf[["pgm.degrad.2010real"]]), na.rm=T)

### non-degraded sites will be considered with 300 years following (Vieira et al. 2005 - https://doi.org/10.1073/pnas.0505966102)
pgm.degrad.pf[["pgm.degrad.2010real"]][pgm.degrad.pf[["pgm.degrad.2010real"]]>25] <- 300
pgm.degrad.pf[["pgm.degrad.2020real"]][pgm.degrad.pf[["pgm.degrad.2020real"]]>35] <- 300

### isolating degraded primary forest class pixels
pgm.degrad.2010.forest.class <- pgm.degrad.pf[["pgm.degrad.2010real"]]
#pgm.degrad.2010.forest.class <- mask(pgm.degrad.2010.forest.class, pgm.lulc.2010.forest.mask)
pgm.degrad.2010.forest.class[pgm.degrad.2010.forest.class>25]<-NA
pgm.degrad.2010.forest.class[pgm.degrad.2010.forest.class<=25]<-1

pgm.degrad.2010.mask <- pgm.degrad.2010.forest.class

pgm.degrad.2010.forest.class[is.na(pgm.degrad.2010.forest.class)]<-0



pgm.degrad.2020.forest.class <- pgm.degrad.pf[["pgm.degrad.2020real"]]
#pgm.degrad.2020.forest.class <- mask(pgm.degrad.2020.forest.class, pgm.lulc.2020.forest.mask)
pgm.degrad.2020.forest.class[pgm.degrad.2020.forest.class>35]<-NA
pgm.degrad.2020.forest.class[pgm.degrad.2020.forest.class<=35]<-1

pgm.degrad.2020.mask <- pgm.degrad.2020.forest.class

pgm.degrad.2020.forest.class[is.na(pgm.degrad.2020.forest.class)]<-0


### accounting for repeated degradation
pgm.repeateddegrad.2010.mask <- raster("rasters/PGM/raw/pgm-degfreq-1985_2010.tif")
pgm.repeateddegrad.2010.mask[pgm.repeateddegrad.2010.mask<2]<-NA
pgm.repeateddegrad.2010.mask[pgm.repeateddegrad.2010.mask>=2]<-1


pgm.repeateddegrad.2020.mask <- raster("rasters/PGM/raw/pgm-degfreq-1985_2020.tif")
pgm.repeateddegrad.2020.mask[pgm.repeateddegrad.2020.mask<2]<-NA
pgm.repeateddegrad.2020.mask[pgm.repeateddegrad.2020.mask>=2]<-1

#
#





## secondary forest age  [2010 and 2020]
#' source: Silva Jr. et al 2020 [DOI: 10.1038/s41597-020-00600-4]

pgm.sfage <- stack(c("rasters/PGM/raw/pgm-2010-sfage-mapbiomas-brazil-collection-60.tif",
                     "rasters/PGM/raw/pgm-2020-sfage-mapbiomas-brazil-collection-60.tif"))
names(pgm.sfage) <- c("pgm.sfage.2010real", "pgm.sfage.2020real")

#checking
#st_crs(pgm.sfage)==st_crs(pgm.shp)
#plot(pgm.sfage)
#range(values(pgm.sfage[["pgm.sfage.2010real"]]), na.rm = T)

# Conversion of rasters into same extent
#pgm.sfage <- resample(pgm.sfage, pgm.lulc.100, method='ngb')

### excluding non-forest areas
pgm.sfage[["pgm.sfage.2010real"]] <- mask(pgm.sfage[["pgm.sfage.2010real"]], pgm.lulc.2010.forest.mask)
pgm.sfage[["pgm.sfage.2010real"]][is.na(pgm.sfage[["pgm.sfage.2010real"]])] <- 0


pgm.sfage[["pgm.sfage.2020real"]] <- mask(pgm.sfage[["pgm.sfage.2020real"]], pgm.lulc.2020.forest.mask)
pgm.sfage[["pgm.sfage.2020real"]][is.na(pgm.sfage[["pgm.sfage.2020real"]])] <- 0



### isolating secondary forest class pixels
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


## degraded secondary forest
pgm.degrad.sf <- stack(c("rasters/PGM/raw/pgm-2010-dsf-tsince0.tif",
                         "rasters/PGM/raw/pgm-2020-dsf-tsince0.tif"))
names(pgm.degrad.sf) <- c("pgm.degradsf.2010real", "pgm.degradsf.2020real")


### isolating degraded secondary forest class pixels
pgm.dsf.2010.mask <- pgm.degrad.sf[["pgm.degradsf.2010real"]]
pgm.dsf.2010.mask[pgm.dsf.2010.mask>25]<-NA
pgm.dsf.2010.mask[pgm.dsf.2010.mask<=25]<-1


pgm.dsf.2020.mask <- pgm.degrad.sf[["pgm.degradsf.2020real"]]
pgm.dsf.2020.mask[pgm.dsf.2020.mask>35]<-NA
pgm.dsf.2020.mask[pgm.dsf.2020.mask<=35]<-1



#
#



# candidate areas for restoration scenarios ==============|
#' @description [...] (provide a brief description)
#' 
#' see auxiliar.R script for details about the steps to select candidate areas for restoration
candidate.areas.final <- raster("rasters/PGM/raw/pgm_forest_cover_after_restoration_2010.tif")
candidate.areas.final <- sum(pgm.lulc.2010.forest.mask, candidate.areas.final, na.rm = T)
candidate.areas.final[candidate.areas.final!=1] <- 0



#
#



# building real scenarios ======================================================
##Undegraded primary forest
###2010
UPF2010<-pgm.lulc.2010.forest.class
UPF2010<-mask(UPF2010, pgm.sfage.2010.mask, inverse=T)
UPF2010[is.na(UPF2010[])]<-0
UPF2010<-mask(UPF2010, pgm.degrad.2010.mask, inverse=T)
UPF2010[is.na(UPF2010[])]<-0
#plot(UPF2010)
writeRaster(UPF2010, "rasters/PGM/input/LULC/UPF2010_real.tif", format="GTiff", overwrite=T)
#UPF2010 <- raster("rasters/PGM/input/LULC/UPF2010_real.tif")

###2020
UPF2020<-pgm.lulc.2020.forest.class
UPF2020<-mask(UPF2020, pgm.sfage.2020.mask, inverse=T)
UPF2020[is.na(UPF2020[])]<-0
UPF2020<-mask(UPF2020, pgm.degrad.2020.mask, inverse=T)
UPF2020[is.na(UPF2020[])]<-0
#plot(UPF2020)
writeRaster(UPF2020, "rasters/PGM/input/LULC/UPF2020_real.tif", format="GTiff", overwrite=T)
#UPF2020 <- raster("rasters/PGM/input/LULC/UPF2020_real.tif")


##Degraded primary forest
###2010
DPF2010 <- pgm.degrad.2010.forest.class
uDPF2010<-mask(DPF2010, pgm.repeateddegrad.2010.mask, inverse=T)
uDPF2010[is.na(uDPF2010[])]<-0
#plot(uDPF2010)
writeRaster(uDPF2010, "rasters/PGM/input/LULC/uDPF2010_real.tif", format="GTiff", overwrite=T)
#uDPF2010 <- raster("rasters/PGM/input/LULC/uDPF2010_real.tif")

RDPF2010 <- mask(DPF2010, pgm.repeateddegrad.2010.mask)
RDPF2010[is.na(RDPF2010[])]<-0
#plot(RDPF2010)
writeRaster(RDPF2010, "rasters/PGM/input/LULC/RDPF2010_real.tif", format="GTiff", overwrite=T)
#RDPF2010 <- raster("rasters/PGM/input/LULC/RDPF2010_real.tif")

###2020
DPF2020 <- pgm.degrad.2020.forest.class
uDPF2020<-mask(DPF2020, pgm.repeateddegrad.2020.mask, inverse=T)
uDPF2020[is.na(uDPF2020[])]<-0
#plot(uDPF2020)
writeRaster(uDPF2020, "rasters/PGM/input/LULC/uDPF2020_real.tif", format="GTiff", overwrite=T)
#uDPF2020 <- raster("rasters/PGM/input/LULC/uDPF2020_real.tif")

RDPF2020 <- mask(DPF2020, pgm.repeateddegrad.2020.mask)
RDPF2020[is.na(RDPF2020[])]<-0
#plot(RDPF2020)
writeRaster(RDPF2020, "rasters/PGM/input/LULC/RDPF2020_real.tif", format="GTiff", overwrite=T)
#RDPF2020 <- raster("rasters/PGM/input/LULC/RDPF2020_real.tif")



##Secondary forest
### 2010
SF2010 <- pgm.sfage.2010.all.class
uSF2010 <- mask(SF2010, pgm.dsf.2010.mask, inverse=T)
uSF2010[is.na(uSF2010[])]<-0
#plot(uSF2010)
writeRaster(uSF2010, "rasters/PGM/input/LULC/uSF2010_real.tif", format="GTiff", overwrite=T)
#uSF2010 <- raster("rasters/PGM/input/LULC/uSF2010_real.tif")

DSF2010 <- mask(SF2010, pgm.dsf.2010.mask)
DSF2010[is.na(DSF2010[])]<-0
#plot(DSF2010)
writeRaster(DSF2010, "rasters/PGM/input/LULC/DSF2010_real.tif", format="GTiff", overwrite=T)
#DSF2010 <- raster("rasters/PGM/input/LULC/DSF2010_real.tif")

### 2020
SF2020 <- pgm.sfage.2020.all.class
uSF2020 <- mask(SF2020, pgm.dsf.2020.mask, inverse=T)
uSF2020[is.na(uSF2020[])]<-0
#plot(uSF2020)
writeRaster(uSF2020, "rasters/PGM/input/LULC/uSF2020_real.tif", format="GTiff", overwrite=T)
#uSF2020 <- raster("rasters/PGM/input/LULC/uSF2020_real.tif")

DSF2020 <- mask(SF2020, pgm.dsf.2020.mask)
DSF2020[is.na(DSF2020[])]<-0
#plot(DSF2020)
writeRaster(DSF2020, "rasters/PGM/input/LULC/DSF2020_real.tif", format="GTiff", overwrite=T)
#DSF2020 <- raster("rasters/PGM/input/LULC/DSF2020_real.tif")



#checking
#par(mfrow = c(2, 3))
###PGM 2010
#undegraded primary forest == 1
UPF2010.sk <- UPF2010
#UPF2010.sk <- raster("rasters/PGM/input/UPF2010_real.tif")
#UPF2010.sk[UPF2010.sk==1]<-1
UPF2010.sk <- mask(UPF2010.sk, pgm.shp)
UPF2010.sk[UPF2010.sk[]==0] <- 333
#rasterVis::levelplot(ratify(UPF2010.sk), main="undegraded primary forest", col.regions=c("white","darkgreen"), att='ID', colorkey=F)

#degraded primary forest == 10
uDPF2010.sk <- uDPF2010
#uDPF2010.sk <- raster("rasters/PGM/input/uDPF2010_real.tif")
uDPF2010.sk[uDPF2010.sk==1]<-10
uDPF2010.sk <- mask(uDPF2010.sk, pgm.shp)
uDPF2010.sk[uDPF2010.sk[]==0] <- 333
#plot(uDPF2010.sk, main="degraded primary forest", legend=F)

#repeated degraded primary forest == 25
RDPF2010.sk <- RDPF2010
#RDPF2010.sk <- raster("rasters/PGM/input/RDPF2010_real.tif")
RDPF2010.sk[RDPF2010.sk==1]<-25
RDPF2010.sk <- mask(RDPF2010.sk, pgm.shp)
RDPF2010.sk[RDPF2010.sk[]==0] <- 333
#plot(RDPF2010.sk, main="degraded primary forest", legend=F)

#secondary forest == 100
uSF2010.sk <- uSF2010
#uSF2010.sk <- raster("rasters/PGM/input/uSF2010_real.tif")
uSF2010.sk[uSF2010.sk==1]<-100
uSF2010.sk <- mask(uSF2010.sk, pgm.shp)
uSF2010.sk[uSF2010.sk[]==0] <- 333
#plot(uSF2010.sk, main="secondary forest", legend=F)

#degraded secondary forest == 125
DSF2010.sk <- DSF2010
#DSF2010.sk <- raster("rasters/PGM/input/DSF2010_real.tif")
DSF2010.sk[DSF2010.sk==1]<-125
DSF2010.sk <- mask(DSF2010.sk, pgm.shp)
DSF2010.sk[DSF2010.sk[]==0] <- 333
#plot(DSF2010.sk, main="secondary forest", legend=F)

LULC2010 <- sum(UPF2010.sk, uSF2010.sk, na.rm = T)
LULC2010 <- sum(LULC2010, DSF2010.sk, na.rm = T)
LULC2010 <- sum(LULC2010, uDPF2010.sk, na.rm = T)
LULC2010 <- sum(LULC2010, RDPF2010.sk, na.rm = T)
LULC2010 <- mask(LULC2010, pgm.shp)
#sort(unique(LULC2010[]))
LULC2010[LULC2010==1333]<-1
LULC2010[LULC2010==1342]<-10
LULC2010[LULC2010==1357]<-25
LULC2010[LULC2010==1432]<-100
LULC2010[LULC2010==1457]<-125
LULC2010[LULC2010==1665 ]<-0
#plot(LULC2020, main="forest cover 2020", legend=F)

###PGM 2020
#undegraded primary forest == 1
UPF2020.sk <- UPF2020
#UPF2020.sk <- raster("rasters/PGM/input/UPF2020_real.tif")
#UPF2020.sk[UPF2020.sk==1]<-1
UPF2020.sk <- mask(UPF2020.sk, pgm.shp)
UPF2020.sk[UPF2020.sk[]==0] <- 333
#rasterVis::levelplot(ratify(UPF2020.sk), main="undegraded primary forest", col.regions=c("white","darkgreen"), att='ID', colorkey=F)

#degraded primary forest == 10
uDPF2020.sk <- uDPF2020
#uDPF2020.sk <- raster("rasters/PGM/input/uDPF2020_real.tif")
uDPF2020.sk[uDPF2020.sk==1]<-10
uDPF2020.sk <- mask(uDPF2020.sk, pgm.shp)
uDPF2020.sk[uDPF2020.sk[]==0] <- 333
#plot(uDPF2020.sk, main="degradded primary forest", legend=F)

#repeated degraded primary forest == 25
RDPF2020.sk <- RDPF2020
#RDPF2020.sk <- raster("rasters/PGM/input/RDPF2020_real.tif")
RDPF2020.sk[RDPF2020.sk==1]<-25
RDPF2020.sk <- mask(RDPF2020.sk, pgm.shp)
RDPF2020.sk[RDPF2020.sk[]==0] <- 333
#plot(RDPF2020.sk, main="degradded primary forest", legend=F)

#secondary forest == 100
uSF2020.sk <- uSF2020
#uSF2020.sk <- raster("rasters/PGM/input/uSF2020_real.tif")
uSF2020.sk[uSF2020.sk==1]<-100
uSF2020.sk <- mask(uSF2020.sk, pgm.shp)
uSF2020.sk[uSF2020.sk[]==0] <- 333
#plot(uSF2020.sk, main="secondary forest", legend=F)

#degraded secondary forest == 125
DSF2020.sk <- DSF2020
#DSF2020.sk <- raster("rasters/PGM/input/DSF2020_real.tif")
DSF2020.sk[DSF2020.sk==1]<-125
DSF2020.sk <- mask(DSF2020.sk, pgm.shp)
DSF2020.sk[DSF2020.sk[]==0] <- 333
#plot(DSF2020.sk, main="secondary forest", legend=F)

LULC2020 <- sum(UPF2020.sk, uSF2020.sk, na.rm = T)
LULC2020 <- sum(LULC2020, DSF2020.sk, na.rm = T)
LULC2020 <- sum(LULC2020, uDPF2020.sk, na.rm = T)
LULC2020 <- sum(LULC2020, RDPF2020.sk, na.rm = T)
LULC2020 <- mask(LULC2020, pgm.shp)
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
#~1.3% of pixels

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
plot(pgm.shp, add=T)

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

#isolating areas that change
#2010 Real
LULC2010 <- sum(UPF2010, uSF2010, na.rm = T)
LULC2010 <- sum(LULC2010, DSF2010, na.rm = T)
LULC2010 <- sum(LULC2010, uDPF2010, na.rm = T)
LULC2010 <- sum(LULC2010, RDPF2010, na.rm = T)


#2020 Real
LULC2020 <- sum(uSF2020, DSF2020, na.rm = T)
LULC2020 <- sum(LULC2020, uDPF2020, na.rm = T)
LULC2020 <- sum(LULC2020, RDPF2020, na.rm = T)
LULC2020[] <- ifelse(UPF2010[]==1 & UPF2020[]==0, 1, LULC2020[])
LULC2020[] <- ifelse(LULC2020[]==1, 1, NA)
#length(which(LULC2020[]==1))

LULC2010[] <- ifelse(LULC2020[]==1, 1, LULC2010[])
#length(which(LULC2010[]==1))

writeRaster(LULC2010, "rasters/PGM/input/LULC/area_change_2010_real_bin.tif", format="GTiff", overwrite=T)
#LULC2010.bin <- raster("rasters/PGM/input/LULC/area_change_2010_real_bin.tif")

writeRaster(LULC2020, "rasters/PGM/input/LULC/area_change_2020_real_bin.tif", format="GTiff", overwrite=T)
#LULC2020.bin <- raster("rasters/PGM/input/LULC/area_change_2020_real_bin.tif")


LULC2010[] <- ifelse(LULC2010[]==1 & uDPF2010[]==1, 10, LULC2010[])
LULC2010[] <- ifelse(LULC2010[]==1 & RDPF2010[]==1, 25, LULC2010[])
LULC2010[] <- ifelse(LULC2010[]==1 & uSF2010[]==1, 100, LULC2010[])
LULC2010[] <- ifelse(LULC2010[]==1 & DSF2010[]==1, 125, LULC2010[])

writeRaster(LULC2010, "rasters/PGM/input/LULC/area_change_2010_real.tif", format="GTiff", overwrite=T)
#LULC2010 <- raster("rasters/PGM/input/LULC/area_change_2010_real.tif")


LULC2020[] <- ifelse(LULC2020[]==1 & uDPF2020[]==1, 10, LULC2020[])
LULC2020[] <- ifelse(LULC2020[]==1 & RDPF2020[]==1, 25, LULC2020[])
LULC2020[] <- ifelse(LULC2020[]==1 & uSF2020[]==1, 100, LULC2020[])
LULC2020[] <- ifelse(LULC2020[]==1 & DSF2020[]==1, 125, LULC2020[])

writeRaster(LULC2020, "rasters/PGM/input/LULC/area_change_2020_real.tif", format="GTiff", overwrite=T)
#LULC2020 <- raster("rasters/PGM/input/LULC/area_change_2020_real.tif")



# building contractual  scenarios ==============================================
##2020 avoid degradation (all)
###Undegraded primary forest
UPF2020_avoiddegrad <- UPF2020
UPF2020_avoiddegrad[] <- ifelse(UPF2010[]==1 & uDPF2020[]==1, 1, UPF2020_avoiddegrad[])
UPF2020_avoiddegrad[] <- ifelse(UPF2010[]==1 & RDPF2020[]==1, 1, UPF2020_avoiddegrad[])
#plot(UPF2020_avoiddegrad)
writeRaster(UPF2020_avoiddegrad, "rasters/PGM/input/LULC/UPF2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#UPF2020_avoiddegrad <- raster("rasters/PGM/input/LULC/UPF2020_avoiddegrad.tif")


###Degraded primary forest
uDPF2020_avoiddegrad <- uDPF2020
uDPF2020_avoiddegrad[] <- ifelse(uDPF2020[]==1 & UPF2010[]==1, 0, uDPF2020_avoiddegrad[])
uDPF2020_avoiddegrad[] <- ifelse(uDPF2010[]==1 & RDPF2020[]==1, 1, uDPF2020_avoiddegrad[])
#plot(uDPF2020_avoiddegrad)
writeRaster(uDPF2020_avoiddegrad, "rasters/PGM/input/LULC/uDPF2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#uDPF2020_avoiddegrad <- raster("rasters/PGM/input/LULC/uDPF2020_avoiddegrad.tif")

RDPF2020_avoiddegrad <- RDPF2020
RDPF2020_avoiddegrad[] <- ifelse(RDPF2020[]==1 & UPF2010[]==1, 0, RDPF2020_avoiddegrad[])
RDPF2020_avoiddegrad[] <- ifelse(RDPF2020[]==1 & uDPF2010[]==1, 0, RDPF2020_avoiddegrad[])
#plot(RDPF2020_avoiddegrad)
writeRaster(RDPF2020_avoiddegrad, "rasters/PGM/input/LULC/RDPF2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#RDPF2020_avoiddegrad <- raster("rasters/PGM/input/LULC/RDPF2020_avoiddegrad.tif")


###Secondary forest
uSF2020_avoiddegrad <- uSF2020
uSF2020_avoiddegrad[] <- ifelse(uSF2010[]==1 & DSF2020[]==1, 1, uSF2020_avoiddegrad[])
#plot(uSF2020_avoiddegrad)
writeRaster(uSF2020_avoiddegrad, "rasters/PGM/input/LULC/uSF2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#uSF2020_avoiddegrad <- raster("rasters/PGM/input/LULC/uSF2020_avoiddegrad.tif")

DSF2020_avoiddegrad <- DSF2020
DSF2020_avoiddegrad[] <- ifelse(DSF2020[]==1 & uSF2010[]==1, 0, DSF2020_avoiddegrad[])
#plot(DSF2020_avoiddegrad)
writeRaster(DSF2020_avoiddegrad, "rasters/PGM/input/LULC/DSF2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#DSF2020_avoiddegrad <- raster("rasters/PGM/input/LULC/DSF2020_avoiddegrad.tif")

#isolating areas that change
LULC2020_avoiddegrad <- LULC2020.bin
LULC2020_avoiddegrad[] <- ifelse(UPF2010[]==1 & uDPF2020[]==1, 1, NA)
LULC2020_avoiddegrad[] <- ifelse(UPF2010[]==1 & RDPF2020[]==1, 1, LULC2020_avoiddegrad[])
LULC2020_avoiddegrad[] <- ifelse(uDPF2010[]==1 & RDPF2020[]==1, 1, LULC2020_avoiddegrad[])
LULC2020_avoiddegrad[] <- ifelse(uSF2010[]==1 & DSF2020[]==1, 1, LULC2020_avoiddegrad[])
#length(which(LULC2020_avoiddegrad[]==1))
writeRaster(LULC2020_avoiddegrad, "rasters/PGM/input/LULC/area_change_2020_avoiddegrad_bin.tif", format="GTiff", overwrite=T)
#LULC2020_avoiddegrad.bin <- raster("rasters/PGM/input/LULC/area_change_2020_avoiddegrad_bin.tif")


LULC2020_avoiddegrad[] <- ifelse(LULC2020_avoiddegrad[]==1 & uDPF2020_avoiddegrad[]==1, 10, LULC2020_avoiddegrad[])
LULC2020_avoiddegrad[] <- ifelse(LULC2020_avoiddegrad[]==1 & RDPF2020_avoiddegrad[]==1, 25, LULC2020_avoiddegrad[])
LULC2020_avoiddegrad[] <- ifelse(LULC2020_avoiddegrad[]==1 & uSF2020_avoiddegrad[]==1, 100, LULC2020_avoiddegrad[])
LULC2020_avoiddegrad[] <- ifelse(LULC2020_avoiddegrad[]==1 & DSF2020_avoiddegrad[]==1, 125, LULC2020_avoiddegrad[])

writeRaster(LULC2020_avoiddegrad, "rasters/PGM/input/LULC/area_change_2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#LULC2020_avoiddegrad <- raster("rasters/PGM/input/LULC/area_change_2020_avoiddegrad.tif")



##2020 avoid degradation (Primary forest only)
###Undegraded primary forest
UPF2020_avoiddegrad2 <- UPF2020
UPF2020_avoiddegrad2[] <- ifelse(UPF2010[]==1 & uDPF2020[]==1, 1, UPF2020_avoiddegrad2[])
UPF2020_avoiddegrad2[] <- ifelse(UPF2010[]==1 & RDPF2020[]==1, 1, UPF2020_avoiddegrad2[])
#plot(UPF2020_avoiddegrad2)
writeRaster(UPF2020_avoiddegrad2, "rasters/PGM/input/LULC/UPF2020_avoiddegrad2.tif", format="GTiff", overwrite=T)
#UPF2020_avoiddegrad2 <- raster("rasters/PGM/input/LULC/UPF2020_avoiddegrad2.tif")


###Degraded primary forest
uDPF2020_avoiddegrad2 <- uDPF2020
uDPF2020_avoiddegrad2[] <- ifelse(uDPF2020[]==1 & UPF2010[]==1, 0, uDPF2020_avoiddegrad2[])
uDPF2020_avoiddegrad2[] <- ifelse(uDPF2010[]==1 & RDPF2020[]==1, 1, uDPF2020_avoiddegrad2[])
#plot(uDPF2020_avoiddegrad2)
writeRaster(uDPF2020_avoiddegrad2, "rasters/PGM/input/LULC/uDPF2020_avoiddegrad2.tif", format="GTiff", overwrite=T)
#uDPF2020_avoiddegrad2 <- raster("rasters/PGM/input/LULC/uDPF2020_avoiddegrad2.tif")

RDPF2020_avoiddegrad2 <- RDPF2020
RDPF2020_avoiddegrad2[] <- ifelse(RDPF2020[]==1 & UPF2010[]==1, 0, RDPF2020_avoiddegrad2[])
RDPF2020_avoiddegrad2[] <- ifelse(RDPF2020[]==1 & uDPF2010[]==1, 0, RDPF2020_avoiddegrad2[])
#plot(RDPF2020_avoiddegrad2)
writeRaster(RDPF2020_avoiddegrad2, "rasters/PGM/input/LULC/RDPF2020_avoiddegrad2.tif", format="GTiff", overwrite=T)
#RDPF2020_avoiddegrad2 <- raster("rasters/PGM/input/LULC/RDPF2020_avoiddegrad2.tif")


###Secondary forest
uSF2020_avoiddegrad2 <- uSF2020
#plot(uSF2020_avoiddegrad2)
writeRaster(uSF2020_avoiddegrad2, "rasters/PGM/input/LULC/uSF2020_avoiddegrad2.tif", format="GTiff", overwrite=T)
#uSF2020_avoiddegrad2 <- raster("rasters/PGM/input/LULC/uSF2020_avoiddegrad2.tif")

DSF2020_avoiddegrad2 <- DSF2020
#plot(DSF2020_avoiddegrad2)
writeRaster(DSF2020_avoiddegrad2, "rasters/PGM/input/LULC/DSF2020_avoiddegrad2.tif", format="GTiff", overwrite=T)
#DSF2020_avoiddegrad2 <- raster("rasters/PGM/input/LULC/DSF2020_avoiddegrad2.tif")

#isolating areas that change
LULC2020_avoiddegrad2 <- LULC2020.bin
LULC2020_avoiddegrad2[] <- ifelse(UPF2010[]==1 & uDPF2020[]==1, 1, NA)
LULC2020_avoiddegrad2[] <- ifelse(UPF2010[]==1 & RDPF2020[]==1, 1, LULC2020_avoiddegrad2[])
LULC2020_avoiddegrad2[] <- ifelse(uDPF2010[]==1 & RDPF2020[]==1, 1, LULC2020_avoiddegrad2[])
#length(which(LULC2020_avoiddegrad2[]==1))
writeRaster(LULC2020_avoiddegrad2, "rasters/PGM/input/LULC/area_change_2020_avoiddegrad2_bin.tif", format="GTiff", overwrite=T)
#LULC2020_avoiddegrad2.bin <- raster("rasters/PGM/input/LULC/area_change_2020_avoiddegrad_bin.tif")


LULC2020_avoiddegrad2[] <- ifelse(LULC2020_avoiddegrad2[]==1 & uDPF2020_avoiddegrad2[]==1, 10, LULC2020_avoiddegrad2[])
LULC2020_avoiddegrad2[] <- ifelse(LULC2020_avoiddegrad2[]==1 & RDPF2020_avoiddegrad2[]==1, 25, LULC2020_avoiddegrad2[])
LULC2020_avoiddegrad2[] <- ifelse(LULC2020_avoiddegrad2[]==1 & uSF2020_avoiddegrad2[]==1, 100, LULC2020_avoiddegrad2[])
LULC2020_avoiddegrad2[] <- ifelse(LULC2020_avoiddegrad2[]==1 & DSF2020_avoiddegrad2[]==1, 125, LULC2020_avoiddegrad2[])

writeRaster(LULC2020_avoiddegrad2, "rasters/PGM/input/LULC/area_change_2020_avoiddegrad2.tif", format="GTiff", overwrite=T)
#LULC2020_avoiddegrad <- raster("rasters/PGM/input/LULC/area_change_2020_avoiddegrad2.tif")



##2020 avoid deforestation (all)
###Undegraded primary forest
UPF2020_avoiddeforest <- UPF2020
UPF2020_avoiddeforest[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, UPF2020_avoiddeforest[])
UPF2020_avoiddeforest[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, UPF2020_avoiddeforest[])
UPF2020_avoiddeforest[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, UPF2020_avoiddeforest[])
#plot(UPF2020_avoiddeforest)
writeRaster(UPF2020_avoiddeforest, "rasters/PGM/input/LULC/UPF2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#UPF2020_avoiddeforest <- raster("rasters/PGM/input/LULC/UPF2020_avoiddeforest.tif")


###Degraded primary forest
uDPF2020_avoiddeforest <- uDPF2020
uDPF2020_avoiddeforest[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, uDPF2020_avoiddeforest[])
uDPF2020_avoiddeforest[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, uDPF2020_avoiddeforest[])
uDPF2020_avoiddeforest[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, uDPF2020_avoiddeforest[])
#plot(uDPF2020_avoiddeforest)
writeRaster(uDPF2020_avoiddeforest, "rasters/PGM/input/LULC/uDPF2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#uDPF2020_avoiddeforest <- raster("rasters/PGM/input/LULC/uDPF2020_avoiddeforest.tif")

RDPF2020_avoiddeforest <- RDPF2020
RDPF2020_avoiddeforest[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, RDPF2020_avoiddeforest[])
RDPF2020_avoiddeforest[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, RDPF2020_avoiddeforest[])
RDPF2020_avoiddeforest[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, RDPF2020_avoiddeforest[])
#plot(RDPF2020_avoiddeforest)
writeRaster(RDPF2020_avoiddeforest, "rasters/PGM/input/LULC/RDPF2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#RDPF2020_avoiddeforest <- raster("rasters/PGM/input/LULC/RDPF2020_avoiddeforest.tif")


###Secondary forest
uSF2020_avoiddeforest <- uSF2020
uSF2020_avoiddeforest[] <- ifelse(uSF2020[]==1 & UPF2010[]==1, 0, uSF2020_avoiddeforest[])
uSF2020_avoiddeforest[] <- ifelse(uSF2020[]==1 & uDPF2010[]==1, 0, uSF2020_avoiddeforest[])
uSF2020_avoiddeforest[] <- ifelse(uSF2020[]==1 & RDPF2010[]==1, 0, uSF2020_avoiddeforest[])
uSF2020_avoiddeforest[] <- ifelse(DSF2020[]==0 & uSF2020[]==0 & uSF2010[]==1, 1, uSF2020_avoiddeforest[])
#plot(uSF2020_avoiddeforest)
writeRaster(uSF2020_avoiddeforest, "rasters/PGM/input/LULC/uSF2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#uSF2020_avoiddeforest <- raster("rasters/PGM/input/LULC/uSF2020_avoiddeforest.tif")

DSF2020_avoiddeforest <- DSF2020
DSF2020_avoiddeforest[] <- ifelse(DSF2020[]==1 & UPF2010[]==1, 0, DSF2020_avoiddeforest[])
DSF2020_avoiddeforest[] <- ifelse(DSF2020[]==1 & uDPF2010[]==1, 0, DSF2020_avoiddeforest[])
DSF2020_avoiddeforest[] <- ifelse(DSF2020[]==1 & RDPF2010[]==1, 0, DSF2020_avoiddeforest[])
DSF2020_avoiddeforest[] <- ifelse(uSF2020[]==0 & DSF2020[]==0 & DSF2010[]==1, 1, DSF2020_avoiddeforest[])
#plot(DSF2020_avoiddeforest)
writeRaster(DSF2020_avoiddeforest, "rasters/PGM/input/LULC/DSF2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#DSF2020_avoiddeforest <- raster("rasters/PGM/input/LULC/DSF2020_avoiddeforest.tif")

#isolating areas that change
LULC2020_avoiddeforest <- LULC2020.bin
LULC2020_avoiddeforest[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, NA)
LULC2020_avoiddeforest[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, LULC2020_avoiddeforest[])
LULC2020_avoiddeforest[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, LULC2020_avoiddeforest[])
LULC2020_avoiddeforest[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, LULC2020_avoiddeforest[])
LULC2020_avoiddeforest[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, LULC2020_avoiddeforest[])
LULC2020_avoiddeforest[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, LULC2020_avoiddeforest[])
LULC2020_avoiddeforest[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, LULC2020_avoiddeforest[])
LULC2020_avoiddeforest[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, LULC2020_avoiddeforest[])
LULC2020_avoiddeforest[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, LULC2020_avoiddeforest[])
LULC2020_avoiddeforest[] <- ifelse(DSF2020[]==0 & uSF2020[]==0 & uSF2010[]==1, 1, LULC2020_avoiddeforest[])
LULC2020_avoiddeforest[] <- ifelse(uSF2020[]==0 & DSF2020[]==0 & DSF2010[]==1, 1, LULC2020_avoiddeforest[])
#length(which(LULC2020_avoiddeforest[]==1))
writeRaster(LULC2020_avoiddeforest, "rasters/PGM/input/LULC/area_change_2020_avoiddeforest_bin.tif", format="GTiff", overwrite=T)
#LULC2020_avoiddeforest.bin <- raster("rasters/PGM/input/LULC/area_change_2020_avoiddeforest_bin.tif")


LULC2020_avoiddeforest[] <- ifelse(LULC2020_avoiddeforest[]==1 & uDPF2020_avoiddeforest[]==1, 10, LULC2020_avoiddeforest[])
LULC2020_avoiddeforest[] <- ifelse(LULC2020_avoiddeforest[]==1 & RDPF2020_avoiddeforest[]==1, 25, LULC2020_avoiddeforest[])
LULC2020_avoiddeforest[] <- ifelse(LULC2020_avoiddeforest[]==1 & uSF2020_avoiddeforest[]==1, 100, LULC2020_avoiddeforest[])
LULC2020_avoiddeforest[] <- ifelse(LULC2020_avoiddeforest[]==1 & DSF2020_avoiddeforest[]==1, 125, LULC2020_avoiddeforest[])

writeRaster(LULC2020_avoiddeforest, "rasters/PGM/input/LULC/area_change_2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#LULC2020_avoiddeforest <- raster("rasters/PGM/input/LULC/area_change_2020_avoiddeforest.tif")



##2020 avoid deforestation (Primary forest only)
###Undegraded primary forest
UPF2020_avoiddeforest2 <- UPF2020
UPF2020_avoiddeforest2[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, UPF2020_avoiddeforest2[])
UPF2020_avoiddeforest2[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, UPF2020_avoiddeforest2[])
UPF2020_avoiddeforest2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, UPF2020_avoiddeforest2[])
#plot(UPF2020_avoiddeforest2)
writeRaster(UPF2020_avoiddeforest2, "rasters/PGM/input/LULC/UPF2020_avoiddeforest2.tif", format="GTiff", overwrite=T)
#UPF2020_avoiddeforest2 <- raster("rasters/PGM/input/LULC/UPF2020_avoiddeforest2.tif")


###Degraded primary forest
uDPF2020_avoiddeforest2 <- uDPF2020
uDPF2020_avoiddeforest2[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, uDPF2020_avoiddeforest2[])
uDPF2020_avoiddeforest2[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, uDPF2020_avoiddeforest2[])
uDPF2020_avoiddeforest2[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, uDPF2020_avoiddeforest2[])
#plot(uDPF2020_avoiddeforest2)
writeRaster(uDPF2020_avoiddeforest2, "rasters/PGM/input/LULC/uDPF2020_avoiddeforest2.tif", format="GTiff", overwrite=T)
#uDPF2020_avoiddeforest2 <- raster("rasters/PGM/input/LULC/uDPF2020_avoiddeforest2.tif")

RDPF2020_avoiddeforest2 <- RDPF2020
RDPF2020_avoiddeforest2[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, RDPF2020_avoiddeforest2[])
RDPF2020_avoiddeforest2[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, RDPF2020_avoiddeforest2[])
RDPF2020_avoiddeforest2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, RDPF2020_avoiddeforest2[])
#plot(RDPF2020_avoiddeforest2)
writeRaster(RDPF2020_avoiddeforest2, "rasters/PGM/input/LULC/RDPF2020_avoiddeforest2.tif", format="GTiff", overwrite=T)
#RDPF2020_avoiddeforest2 <- raster("rasters/PGM/input/LULC/RDPF2020_avoiddeforest2.tif")


###Secondary forest
uSF2020_avoiddeforest2 <- uSF2020
uSF2020_avoiddeforest2[] <- ifelse(uSF2020[]==1 & UPF2010[]==1, 0, uSF2020_avoiddeforest2[])
uSF2020_avoiddeforest2[] <- ifelse(uSF2020[]==1 & uDPF2010[]==1, 0, uSF2020_avoiddeforest2[])
uSF2020_avoiddeforest2[] <- ifelse(uSF2020[]==1 & RDPF2010[]==1, 0, uSF2020_avoiddeforest2[])
#plot(uSF2020_avoiddeforest2)
writeRaster(uSF2020_avoiddeforest2, "rasters/PGM/input/LULC/uSF2020_avoiddeforest2.tif", format="GTiff", overwrite=T)
#uSF2020_avoiddeforest2 <- raster("rasters/PGM/input/LULC/uSF2020_avoiddeforest2.tif")

DSF2020_avoiddeforest2 <- DSF2020
DSF2020_avoiddeforest2[] <- ifelse(DSF2020[]==1 & UPF2010[]==1, 0, DSF2020_avoiddeforest2[])
DSF2020_avoiddeforest2[] <- ifelse(DSF2020[]==1 & uDPF2010[]==1, 0, DSF2020_avoiddeforest2[])
DSF2020_avoiddeforest2[] <- ifelse(DSF2020[]==1 & RDPF2010[]==1, 0, DSF2020_avoiddeforest2[])
#plot(DSF2020_avoiddeforest2)
writeRaster(DSF2020_avoiddeforest2, "rasters/PGM/input/LULC/DSF2020_avoiddeforest2.tif", format="GTiff", overwrite=T)
#DSF2020_avoiddeforest2 <- raster("rasters/PGM/input/LULC/DSF2020_avoiddeforest2.tif")

#isolating areas that change
LULC2020_avoiddeforest2 <- LULC2020.bin
LULC2020_avoiddeforest2[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, NA)
LULC2020_avoiddeforest2[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, LULC2020_avoiddeforest2[])
LULC2020_avoiddeforest2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, LULC2020_avoiddeforest2[])
LULC2020_avoiddeforest2[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, LULC2020_avoiddeforest2[])
LULC2020_avoiddeforest2[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, LULC2020_avoiddeforest2[])
LULC2020_avoiddeforest2[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, LULC2020_avoiddeforest2[])
LULC2020_avoiddeforest2[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, LULC2020_avoiddeforest2[])
LULC2020_avoiddeforest2[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, LULC2020_avoiddeforest2[])
LULC2020_avoiddeforest2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, LULC2020_avoiddeforest2[])
#length(which(LULC2020_avoiddeforest2[]==1))
writeRaster(LULC2020_avoiddeforest2, "rasters/PGM/input/LULC/area_change_2020_avoiddeforest2_bin.tif", format="GTiff", overwrite=T)
#LULC2020_avoiddeforest2.bin <- raster("rasters/PGM/input/LULC/area_change_2020_avoiddeforest2_bin.tif")


LULC2020_avoiddeforest2[] <- ifelse(LULC2020_avoiddeforest2[]==1 & uDPF2020_avoiddeforest2[]==1, 10, LULC2020_avoiddeforest2[])
LULC2020_avoiddeforest2[] <- ifelse(LULC2020_avoiddeforest2[]==1 & RDPF2020_avoiddeforest2[]==1, 25, LULC2020_avoiddeforest2[])
LULC2020_avoiddeforest2[] <- ifelse(LULC2020_avoiddeforest2[]==1 & uSF2020_avoiddeforest2[]==1, 100, LULC2020_avoiddeforest2[])
LULC2020_avoiddeforest2[] <- ifelse(LULC2020_avoiddeforest2[]==1 & DSF2020_avoiddeforest2[]==1, 125, LULC2020_avoiddeforest2[])

writeRaster(LULC2020_avoiddeforest2, "rasters/PGM/input/LULC/area_change_2020_avoiddeforest2.tif", format="GTiff", overwrite=T)
#LULC2020_avoiddeforest2 <- raster("rasters/PGM/input/LULC/area_change_2020_avoiddeforest2.tif")



##2020 restoration without avoid
###Secondary forest
uSF2020_restor_wo_avoid <- sum(uSF2020, candidate.areas.final)
uSF2020_restor_wo_avoid[] <-  ifelse(uSF2020_restor_wo_avoid[]!=0, 1, uSF2020_restor_wo_avoid[])
uSF2020_restor_wo_avoid[] <-  ifelse(uSF2020_restor_wo_avoid[]==1 & DSF2020[]==1, 0, uSF2020_restor_wo_avoid[])
uSF2020_restor_wo_avoid[] <-  ifelse(uSF2020_restor_wo_avoid[]==1 & UPF2020[]==1, 0, uSF2020_restor_wo_avoid[])
uSF2020_restor_wo_avoid[] <-  ifelse(uSF2020_restor_wo_avoid[]==1 & uDPF2020[]==1, 0, uSF2020_restor_wo_avoid[])
uSF2020_restor_wo_avoid[] <-  ifelse(uSF2020_restor_wo_avoid[]==1 & RDPF2020[]==1, 0, uSF2020_restor_wo_avoid[])
#plot(uSF2020_restor_wo_avoid)
writeRaster(uSF2020_restor_wo_avoid, "rasters/PGM/input/LULC/uSF2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)
#uSF2020_restor_wo_avoid <- raster("rasters/PGM/input/LULC/uSF2020_restor_wo_avoid.tif")

DSF2020_restor_wo_avoid <- DSF2020
#plot(DSF2020_restor_wo_avoid)
writeRaster(DSF2020_restor_wo_avoid, "rasters/PGM/input/LULC/DSF2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)
#DSF2020_restor_wo_avoid <- raster("rasters/PGM/input/LULC/DSF2020_restor_wo_avoid.tif")


###Undegraded primary forest
UPF2020_restor_wo_avoid <- UPF2020
#plot(UPF2020_restor_wo_avoid)
writeRaster(UPF2020_restor_wo_avoid, "rasters/PGM/input/LULC/UPF2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)
#UPF2020_restor_wo_avoid <- raster("rasters/PGM/input/LULC/UPF2020_restor_wo_avoid.tif")


###Degraded primary forest
uDPF2020_restor_wo_avoid <- uDPF2020
#plot(uDPF2020_restor_wo_avoid)
writeRaster(uDPF2020_restor_wo_avoid, "rasters/PGM/input/LULC/uDPF2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)
#uDPF2020_restor_wo_avoid <- raster("rasters/PGM/input/LULC/uDPF2020_restor_wo_avoid.tif")

RDPF2020_restor_wo_avoid <- RDPF2020
#plot(RDPF2020_restor_wo_avoid)
writeRaster(RDPF2020_restor_wo_avoid, "rasters/PGM/input/LULC/RDPF2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)
#RDPF2020_restor_wo_avoid <- raster("rasters/PGM/input/LULC/RDPF2020_restor_wo_avoid.tif")

#isolating areas that change
LULC2020_restor_wo_avoid <- candidate.areas.final
LULC2020_restor_wo_avoid[] <- ifelse(LULC2020_restor_wo_avoid[]==0, NA, LULC2020_restor_wo_avoid[])
#length(which(LULC2020_restor_wo_avoid[]==1))
writeRaster(LULC2020_restor_wo_avoid, "rasters/PGM/input/LULC/area_change_2020_restor_wo_avoid_bin.tif", format="GTiff", overwrite=T)
#LULC2020_restor_wo_avoid.bin <- raster("rasters/PGM/input/LULC/area_change_2020_restor_wo_avoid_bin.tif")


LULC2020_restor_wo_avoid[LULC2020_restor_wo_avoid[]==1] <- 100

writeRaster(LULC2020_restor_wo_avoid, "rasters/PGM/input/LULC/area_change_2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)
#LULC2020_restor_wo_avoid <- raster("rasters/PGM/input/LULC/area_change_2020_restor_wo_avoid.tif")



##2020 avoid both (all)
###Undegraded primary forest
UPF2020_avoidboth <- UPF2020
UPF2020_avoidboth[] <- ifelse(UPF2010[]==1 & uDPF2020[]==1, 1, UPF2020_avoidboth[])
UPF2020_avoidboth[] <- ifelse(UPF2010[]==1 & RDPF2020[]==1, 1, UPF2020_avoidboth[])
UPF2020_avoidboth[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, UPF2020_avoidboth[])
UPF2020_avoidboth[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, UPF2020_avoidboth[])
UPF2020_avoidboth[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, UPF2020_avoidboth[])
#plot(UPF2020_avoidboth)
writeRaster(UPF2020_avoidboth, "rasters/PGM/input/LULC/UPF2020_avoidboth.tif", format="GTiff", overwrite=T)
#UPF2020_avoidboth <- raster("rasters/PGM/input/LULC/UPF2020_avoidboth.tif")


###Degraded primary forest
uDPF2020_avoidboth <- uDPF2020
uDPF2020_avoidboth[] <- ifelse(uDPF2020[]==1 & UPF2010[]==1, 0, uDPF2020_avoidboth[])
uDPF2020_avoidboth[] <- ifelse(uDPF2010[]==1 & RDPF2020[]==1, 1, uDPF2020_avoidboth[])
uDPF2020_avoidboth[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, uDPF2020_avoidboth[])
uDPF2020_avoidboth[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, uDPF2020_avoidboth[])
uDPF2020_avoidboth[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, uDPF2020_avoidboth[])
#plot(uDPF2020_avoidboth)
writeRaster(uDPF2020_avoidboth, "rasters/PGM/input/LULC/uDPF2020_avoidboth.tif", format="GTiff", overwrite=T)
#uDPF2020_avoidboth <- raster("rasters/PGM/input/LULC/uDPF2020_avoidboth.tif")

RDPF2020_avoidboth <- RDPF2020
RDPF2020_avoidboth[] <- ifelse(RDPF2020[]==1 & UPF2010[]==1, 0, RDPF2020_avoidboth[])
RDPF2020_avoidboth[] <- ifelse(RDPF2020[]==1 & uDPF2010[]==1, 0, RDPF2020_avoidboth[])
RDPF2020_avoidboth[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, RDPF2020_avoidboth[])
RDPF2020_avoidboth[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, RDPF2020_avoidboth[])
RDPF2020_avoidboth[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, RDPF2020_avoidboth[])
#plot(RDPF2020_avoidboth)
writeRaster(RDPF2020_avoidboth, "rasters/PGM/input/LULC/RDPF2020_avoidboth.tif", format="GTiff", overwrite=T)
#RDPF2020_avoidboth <- raster("rasters/PGM/input/LULC/RDPF2020_avoidboth.tif")


###Secondary forest
uSF2020_avoidboth <- uSF2020
uSF2020_avoidboth[] <- ifelse(uSF2010[]==1 & DSF2020[]==1, 1, uSF2020_avoidboth[])
uSF2020_avoidboth[] <- ifelse(uSF2020[]==1 & UPF2010[]==1, 0, uSF2020_avoidboth[])
uSF2020_avoidboth[] <- ifelse(uSF2020[]==1 & uDPF2010[]==1, 0, uSF2020_avoidboth[])
uSF2020_avoidboth[] <- ifelse(uSF2020[]==1 & RDPF2010[]==1, 0, uSF2020_avoidboth[])
uSF2020_avoidboth[] <- ifelse(DSF2020[]==0 & uSF2020[]==0 & uSF2010[]==1, 1, uSF2020_avoidboth[])
#plot(uSF2020_avoidboth)
writeRaster(uSF2020_avoidboth, "rasters/PGM/input/LULC/uSF2020_avoidboth.tif", format="GTiff", overwrite=T)
#uSF2020_avoidboth <- raster("rasters/PGM/input/LULC/uSF2020_avoidboth.tif")

DSF2020_avoidboth <- DSF2020
DSF2020_avoidboth[] <- ifelse(DSF2020[]==1 & uSF2010[]==1, 0, DSF2020_avoidboth[])
DSF2020_avoidboth[] <- ifelse(DSF2020[]==1 & UPF2010[]==1, 0, DSF2020_avoidboth[])
DSF2020_avoidboth[] <- ifelse(DSF2020[]==1 & uDPF2010[]==1, 0, DSF2020_avoidboth[])
DSF2020_avoidboth[] <- ifelse(DSF2020[]==1 & RDPF2010[]==1, 0, DSF2020_avoidboth[])
DSF2020_avoidboth[] <- ifelse(uSF2020[]==0 & DSF2020[]==0 & DSF2010[]==1, 1, DSF2020_avoidboth[])
#plot(DSF2020_avoidboth)
writeRaster(DSF2020_avoidboth, "rasters/PGM/input/LULC/DSF2020_avoidboth.tif", format="GTiff", overwrite=T)
#DSF2020_avoidboth <- raster("rasters/PGM/input/LULC/DSF2020_avoidboth.tif")

#isolating areas that change
LULC2020_avoidboth <- LULC2020.bin
LULC2020_avoidboth[] <- ifelse(UPF2010[]==1 & uDPF2020[]==1, 1, NA)
LULC2020_avoidboth[] <- ifelse(UPF2010[]==1 & RDPF2020[]==1, 1, LULC2020_avoidboth[])
LULC2020_avoidboth[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, LULC2020_avoidboth[])
LULC2020_avoidboth[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, LULC2020_avoidboth[])
LULC2020_avoidboth[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, LULC2020_avoidboth[])
LULC2020_avoidboth[] <- ifelse(uDPF2010[]==1 & RDPF2020[]==1, 1, LULC2020_avoidboth[])
LULC2020_avoidboth[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, LULC2020_avoidboth[])
LULC2020_avoidboth[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, LULC2020_avoidboth[])
LULC2020_avoidboth[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, LULC2020_avoidboth[])
LULC2020_avoidboth[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, LULC2020_avoidboth[])
LULC2020_avoidboth[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, LULC2020_avoidboth[])
LULC2020_avoidboth[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, LULC2020_avoidboth[])
LULC2020_avoidboth[] <- ifelse(uSF2010[]==1 & DSF2020[]==1, 1, LULC2020_avoidboth[])
LULC2020_avoidboth[] <- ifelse(DSF2020[]==0 & uSF2020[]==0 & uSF2010[]==1, 1, LULC2020_avoidboth[])
LULC2020_avoidboth[] <- ifelse(uSF2020[]==0 & DSF2020[]==0 & DSF2010[]==1, 1, LULC2020_avoidboth[])
#length(which(LULC2020_avoidboth[]==1))
writeRaster(LULC2020_avoidboth, "rasters/PGM/input/LULC/area_change_2020_avoidboth_bin.tif", format="GTiff", overwrite=T)
#LULC2020_avoidboth.bin <- raster("rasters/PGM/input/LULC/area_change_2020_avoidboth_bin.tif")


LULC2020_avoidboth[] <- ifelse(LULC2020_avoidboth[]==1 & uDPF2020_avoidboth[]==1, 10, LULC2020_avoidboth[])
LULC2020_avoidboth[] <- ifelse(LULC2020_avoidboth[]==1 & RDPF2020_avoidboth[]==1, 25, LULC2020_avoidboth[])
LULC2020_avoidboth[] <- ifelse(LULC2020_avoidboth[]==1 & uSF2020_avoidboth[]==1, 100, LULC2020_avoidboth[])
LULC2020_avoidboth[] <- ifelse(LULC2020_avoidboth[]==1 & DSF2020_avoidboth[]==1, 125, LULC2020_avoidboth[])

writeRaster(LULC2020_avoidboth, "rasters/PGM/input/LULC/area_change_2020_avoidboth.tif", format="GTiff", overwrite=T)
#LULC2020_avoidboth <- raster("rasters/PGM/input/LULC/area_change_2020_avoidboth.tif")



##2020 avoid both (Primary forest only)
###Undegraded primary forest
UPF2020_avoidboth2 <- UPF2020
UPF2020_avoidboth2[] <- ifelse(UPF2010[]==1 & uDPF2020[]==1, 1, UPF2020_avoidboth2[])
UPF2020_avoidboth2[] <- ifelse(UPF2010[]==1 & RDPF2020[]==1, 1, UPF2020_avoidboth2[])
UPF2020_avoidboth2[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, UPF2020_avoidboth2[])
UPF2020_avoidboth2[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, UPF2020_avoidboth2[])
UPF2020_avoidboth2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, UPF2020_avoidboth2[])
#plot(UPF2020_avoidboth2)
writeRaster(UPF2020_avoidboth2, "rasters/PGM/input/LULC/UPF2020_avoidboth2.tif", format="GTiff", overwrite=T)
#UPF2020_avoidboth2 <- raster("rasters/PGM/input/LULC/UPF2020_avoidboth2.tif")


###Degraded primary forest
uDPF2020_avoidboth2 <- uDPF2020
uDPF2020_avoidboth2[] <- ifelse(uDPF2020[]==1 & UPF2010[]==1, 0, uDPF2020_avoidboth2[])
uDPF2020_avoidboth2[] <- ifelse(uDPF2010[]==1 & RDPF2020[]==1, 1, uDPF2020_avoidboth2[])
uDPF2020_avoidboth2[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, uDPF2020_avoidboth2[])
uDPF2020_avoidboth2[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, uDPF2020_avoidboth2[])
uDPF2020_avoidboth2[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, uDPF2020_avoidboth2[])
#plot(uDPF2020_avoidboth2)
writeRaster(uDPF2020_avoidboth2, "rasters/PGM/input/LULC/uDPF2020_avoidboth2.tif", format="GTiff", overwrite=T)
#uDPF2020_avoidboth2 <- raster("rasters/PGM/input/LULC/uDPF2020_avoidboth2.tif")

RDPF2020_avoidboth2 <- RDPF2020
RDPF2020_avoidboth2[] <- ifelse(RDPF2020[]==1 & UPF2010[]==1, 0, RDPF2020_avoidboth2[])
RDPF2020_avoidboth2[] <- ifelse(RDPF2020[]==1 & uDPF2010[]==1, 0, RDPF2020_avoidboth2[])
RDPF2020_avoidboth2[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, RDPF2020_avoidboth2[])
RDPF2020_avoidboth2[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, RDPF2020_avoidboth2[])
RDPF2020_avoidboth2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, RDPF2020_avoidboth2[])
#plot(RDPF2020_avoidboth2)
writeRaster(RDPF2020_avoidboth2, "rasters/PGM/input/LULC/RDPF2020_avoidboth2.tif", format="GTiff", overwrite=T)
#RDPF2020_avoidboth2 <- raster("rasters/PGM/input/LULC/RDPF2020_avoidboth2.tif")


###Secondary forest
uSF2020_avoidboth2 <- uSF2020
uSF2020_avoidboth2[] <- ifelse(uSF2020[]==1 & UPF2010[]==1, 0, uSF2020_avoidboth2[])
uSF2020_avoidboth2[] <- ifelse(uSF2020[]==1 & uDPF2010[]==1, 0, uSF2020_avoidboth2[])
uSF2020_avoidboth2[] <- ifelse(uSF2020[]==1 & RDPF2010[]==1, 0, uSF2020_avoidboth2[])
#plot(uSF2020_avoidboth2)
writeRaster(uSF2020_avoidboth2, "rasters/PGM/input/LULC/uSF2020_avoidboth2.tif", format="GTiff", overwrite=T)
#uSF2020_avoidboth2 <- raster("rasters/PGM/input/LULC/uSF2020_avoidboth2.tif")

DSF2020_avoidboth2 <- DSF2020
DSF2020_avoidboth2[] <- ifelse(DSF2020[]==1 & UPF2010[]==1, 0, DSF2020_avoidboth2[])
DSF2020_avoidboth2[] <- ifelse(DSF2020[]==1 & uDPF2010[]==1, 0, DSF2020_avoidboth2[])
DSF2020_avoidboth2[] <- ifelse(DSF2020[]==1 & RDPF2010[]==1, 0, DSF2020_avoidboth2[])
#plot(DSF2020_avoidboth2)
writeRaster(DSF2020_avoidboth2, "rasters/PGM/input/LULC/DSF2020_avoidboth2.tif", format="GTiff", overwrite=T)
#DSF2020_avoidboth2 <- raster("rasters/PGM/input/LULC/DSF2020_avoidboth2.tif")

#isolating areas that change
LULC2020_avoidboth2 <- LULC2020.bin
LULC2020_avoidboth2[] <- ifelse(UPF2010[]==1 & uDPF2020[]==1, 1, NA)
LULC2020_avoidboth2[] <- ifelse(UPF2010[]==1 & RDPF2020[]==1, 1, LULC2020_avoidboth2[])
LULC2020_avoidboth2[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, LULC2020_avoidboth2[])
LULC2020_avoidboth2[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, LULC2020_avoidboth2[])
LULC2020_avoidboth2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, LULC2020_avoidboth2[])
LULC2020_avoidboth2[] <- ifelse(uDPF2010[]==1 & RDPF2020[]==1, 1, LULC2020_avoidboth2[])
LULC2020_avoidboth2[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, LULC2020_avoidboth2[])
LULC2020_avoidboth2[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, LULC2020_avoidboth2[])
LULC2020_avoidboth2[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, LULC2020_avoidboth2[])
LULC2020_avoidboth2[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, LULC2020_avoidboth2[])
LULC2020_avoidboth2[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, LULC2020_avoidboth2[])
LULC2020_avoidboth2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, LULC2020_avoidboth2[])
#length(which(LULC2020_avoidboth2[]==1))
writeRaster(LULC2020_avoidboth2, "rasters/PGM/input/LULC/area_change_2020_avoidboth2_bin.tif", format="GTiff", overwrite=T)
#LULC2020_avoidboth2.bin <- raster("rasters/PGM/input/LULC/area_change_2020_avoidboth2_bin.tif")


LULC2020_avoidboth2[] <- ifelse(LULC2020_avoidboth2[]==1 & uDPF2020_avoidboth2[]==1, 10, LULC2020_avoidboth2[])
LULC2020_avoidboth2[] <- ifelse(LULC2020_avoidboth2[]==1 & RDPF2020_avoidboth2[]==1, 25, LULC2020_avoidboth2[])
LULC2020_avoidboth2[] <- ifelse(LULC2020_avoidboth2[]==1 & uSF2020_avoidboth2[]==1, 100, LULC2020_avoidboth2[])
LULC2020_avoidboth2[] <- ifelse(LULC2020_avoidboth2[]==1 & DSF2020_avoidboth2[]==1, 125, LULC2020_avoidboth2[])

writeRaster(LULC2020_avoidboth2, "rasters/PGM/input/LULC/area_change_2020_avoidboth2.tif", format="GTiff", overwrite=T)
#LULC2020_avoidboth2 <- raster("rasters/PGM/input/LULC/area_change_2020_avoidboth2.tif")



##2020 restoration and avoid deforestation (all)
###Secondary forest
uSF2020_restor_n_avoiddeforest <- sum(uSF2020, candidate.areas.final)
uSF2020_restor_n_avoiddeforest[] <-  ifelse(uSF2020_restor_n_avoiddeforest[]!=0, 1, uSF2020_restor_n_avoiddeforest[])

uSF2020_restor_n_avoiddeforest[] <- ifelse(uSF2020[]==1 & UPF2010[]==1, 0, uSF2020_restor_n_avoiddeforest[])
uSF2020_restor_n_avoiddeforest[] <- ifelse(uSF2020[]==1 & uDPF2010[]==1, 0, uSF2020_restor_n_avoiddeforest[])
uSF2020_restor_n_avoiddeforest[] <- ifelse(uSF2020[]==1 & RDPF2010[]==1, 0, uSF2020_restor_n_avoiddeforest[])
uSF2020_restor_n_avoiddeforest[] <- ifelse(DSF2020[]==0 & uSF2020[]==0 & uSF2010[]==1, 1, uSF2020_restor_n_avoiddeforest[])

uSF2020_restor_n_avoiddeforest[] <-  ifelse(uSF2020_restor_n_avoiddeforest[]==1 & DSF2020[]==1, 0, uSF2020_restor_n_avoiddeforest[])
uSF2020_restor_n_avoiddeforest[] <-  ifelse(uSF2020_restor_n_avoiddeforest[]==1 & UPF2020[]==1, 0, uSF2020_restor_n_avoiddeforest[])
uSF2020_restor_n_avoiddeforest[] <-  ifelse(uSF2020_restor_n_avoiddeforest[]==1 & uDPF2020[]==1, 0, uSF2020_restor_n_avoiddeforest[])
uSF2020_restor_n_avoiddeforest[] <-  ifelse(uSF2020_restor_n_avoiddeforest[]==1 & RDPF2020[]==1, 0, uSF2020_restor_n_avoiddeforest[])

#plot(uSF2020_restor_n_avoiddeforest)
writeRaster(uSF2020_restor_n_avoiddeforest, "rasters/PGM/input/LULC/uSF2020_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)
#uSF2020_restor_n_avoiddeforest <- raster("rasters/PGM/input/LULC/uSF2020_restor_n_avoiddeforest.tif")

DSF2020_restor_n_avoiddeforest <- DSF2020
DSF2020_restor_n_avoiddeforest[] <- ifelse(DSF2020[]==1 & UPF2010[]==1, 0, DSF2020_restor_n_avoiddeforest[])
DSF2020_restor_n_avoiddeforest[] <- ifelse(DSF2020[]==1 & uDPF2010[]==1, 0, DSF2020_restor_n_avoiddeforest[])
DSF2020_restor_n_avoiddeforest[] <- ifelse(DSF2020[]==1 & RDPF2010[]==1, 0, DSF2020_restor_n_avoiddeforest[])
DSF2020_restor_n_avoiddeforest[] <- ifelse(uSF2020[]==0 & DSF2020[]==0 & DSF2010[]==1, 1, DSF2020_restor_n_avoiddeforest[])

DSF2020_restor_n_avoiddeforest[] <-  ifelse(DSF2020_restor_n_avoiddeforest[]==1 & uDPF2020[]==1, 0, DSF2020_restor_n_avoiddeforest[])
DSF2020_restor_n_avoiddeforest[] <-  ifelse(DSF2020_restor_n_avoiddeforest[]==1 & RDPF2020[]==1, 0, DSF2020_restor_n_avoiddeforest[])

#plot(DSF2020_restor_n_avoiddeforest)
writeRaster(DSF2020_restor_n_avoiddeforest, "rasters/PGM/input/LULC/DSF2020_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)
#DSF2020_restor_n_avoiddeforest <- raster("rasters/PGM/input/LULC/DSF2020_restor_n_avoiddeforest.tif")


###Undegraded primary forest
UPF2020_restor_n_avoiddeforest <- UPF2020
UPF2020_restor_n_avoiddeforest[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, UPF2020_restor_n_avoiddeforest[])
UPF2020_restor_n_avoiddeforest[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, UPF2020_restor_n_avoiddeforest[])
UPF2020_restor_n_avoiddeforest[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, UPF2020_restor_n_avoiddeforest[])
#plot(UPF2020_restor_n_avoiddeforest)
writeRaster(UPF2020_restor_n_avoiddeforest, "rasters/PGM/input/LULC/UPF2020_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)
#UPF2020_restor_n_avoiddeforest <- raster("rasters/PGM/input/LULC/UPF2020_restor_n_avoiddeforest.tif")


###Degraded primary forest
uDPF2020_restor_n_avoiddeforest <- uDPF2020
uDPF2020_restor_n_avoiddeforest[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, uDPF2020_restor_n_avoiddeforest[])
uDPF2020_restor_n_avoiddeforest[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, uDPF2020_restor_n_avoiddeforest[])
uDPF2020_restor_n_avoiddeforest[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, uDPF2020_restor_n_avoiddeforest[])
#plot(uDPF2020_restor_n_avoiddeforest)
writeRaster(uDPF2020_restor_n_avoiddeforest, "rasters/PGM/input/LULC/uDPF2020_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)
#uDPF2020_restor_n_avoiddeforest <- raster("rasters/PGM/input/LULC/uDPF2020_restor_n_avoiddeforest.tif")

RDPF2020_restor_n_avoiddeforest <- RDPF2020
RDPF2020_restor_n_avoiddeforest[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, RDPF2020_restor_n_avoiddeforest[])
RDPF2020_restor_n_avoiddeforest[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, RDPF2020_restor_n_avoiddeforest[])
RDPF2020_restor_n_avoiddeforest[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, RDPF2020_restor_n_avoiddeforest[])
#plot(RDPF2020_restor_n_avoiddeforest)
writeRaster(RDPF2020_restor_n_avoiddeforest, "rasters/PGM/input/LULC/RDPF2020_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)
#RDPF2020_restor_n_avoiddeforest <- raster("rasters/PGM/input/LULC/RDPF2020_restor_n_avoiddeforest.tif")

#isolating areas that change
LULC2020_restor_n_avoiddeforest <- candidate.areas.final
LULC2020_restor_n_avoiddeforest[] <- ifelse(LULC2020_restor_n_avoiddeforest[]==0, NA, LULC2020_restor_n_avoiddeforest[])
LULC2020_restor_n_avoiddeforest[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, LULC2020_restor_n_avoiddeforest[])
LULC2020_restor_n_avoiddeforest[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, LULC2020_restor_n_avoiddeforest[])
LULC2020_restor_n_avoiddeforest[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, LULC2020_restor_n_avoiddeforest[])
LULC2020_restor_n_avoiddeforest[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, LULC2020_restor_n_avoiddeforest[])
LULC2020_restor_n_avoiddeforest[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, LULC2020_restor_n_avoiddeforest[])
LULC2020_restor_n_avoiddeforest[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, LULC2020_restor_n_avoiddeforest[])
LULC2020_restor_n_avoiddeforest[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, LULC2020_restor_n_avoiddeforest[])
LULC2020_restor_n_avoiddeforest[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, LULC2020_restor_n_avoiddeforest[])
LULC2020_restor_n_avoiddeforest[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, LULC2020_restor_n_avoiddeforest[])
LULC2020_restor_n_avoiddeforest[] <- ifelse(DSF2020[]==0 & uSF2020[]==0 & uSF2010[]==1, 1, LULC2020_restor_n_avoiddeforest[])
LULC2020_restor_n_avoiddeforest[] <- ifelse(uSF2020[]==0 & DSF2020[]==0 & DSF2010[]==1, 1, LULC2020_restor_n_avoiddeforest[])
#length(which(LULC2020_restor_n_avoiddeforest[]==1))
writeRaster(LULC2020_restor_n_avoiddeforest, "rasters/PGM/input/LULC/area_change_2020_restor_n_avoiddeforest_bin.tif", format="GTiff", overwrite=T)
#LULC2020_restor_n_avoiddeforest.bin <- raster("rasters/PGM/input/LULC/area_change_2020_restor_n_avoiddeforest_bin.tif")


LULC2020_restor_n_avoiddeforest[] <- ifelse(LULC2020_restor_n_avoiddeforest[]==1 & uDPF2020_restor_n_avoiddeforest[]==1, 10, LULC2020_restor_n_avoiddeforest[])
LULC2020_restor_n_avoiddeforest[] <- ifelse(LULC2020_restor_n_avoiddeforest[]==1 & RDPF2020_restor_n_avoiddeforest[]==1, 25, LULC2020_restor_n_avoiddeforest[])
LULC2020_restor_n_avoiddeforest[] <- ifelse(LULC2020_restor_n_avoiddeforest[]==1 & uSF2020_restor_n_avoiddeforest[]==1, 100, LULC2020_restor_n_avoiddeforest[])
LULC2020_restor_n_avoiddeforest[] <- ifelse(LULC2020_restor_n_avoiddeforest[]==1 & DSF2020_restor_n_avoiddeforest[]==1, 125, LULC2020_restor_n_avoiddeforest[])

writeRaster(LULC2020_restor_n_avoiddeforest, "rasters/PGM/input/LULC/area_change_2020_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)
#LULC2020_restor_n_avoiddeforest <- raster("rasters/PGM/input/LULC/area_change_2020_restor_n_avoiddeforest.tif")



##2020 restoration and avoid deforestation (Primary forest only)
###Secondary forest
uSF2020_restor_n_avoiddeforest2 <- sum(uSF2020, candidate.areas.final)
uSF2020_restor_n_avoiddeforest2[] <-  ifelse(uSF2020_restor_n_avoiddeforest2[]!=0, 1, uSF2020_restor_n_avoiddeforest2[])

uSF2020_restor_n_avoiddeforest2[] <- ifelse(uSF2020[]==1 & UPF2010[]==1, 0, uSF2020_restor_n_avoiddeforest2[])
uSF2020_restor_n_avoiddeforest2[] <- ifelse(uSF2020[]==1 & uDPF2010[]==1, 0, uSF2020_restor_n_avoiddeforest2[])
uSF2020_restor_n_avoiddeforest2[] <- ifelse(uSF2020[]==1 & RDPF2010[]==1, 0, uSF2020_restor_n_avoiddeforest2[])

uSF2020_restor_n_avoiddeforest2[] <-  ifelse(uSF2020_restor_n_avoiddeforest2[]==1 & DSF2020[]==1, 0, uSF2020_restor_n_avoiddeforest2[])
uSF2020_restor_n_avoiddeforest2[] <-  ifelse(uSF2020_restor_n_avoiddeforest2[]==1 & UPF2020[]==1, 0, uSF2020_restor_n_avoiddeforest2[])
uSF2020_restor_n_avoiddeforest2[] <-  ifelse(uSF2020_restor_n_avoiddeforest2[]==1 & uDPF2020[]==1, 0, uSF2020_restor_n_avoiddeforest2[])
uSF2020_restor_n_avoiddeforest2[] <-  ifelse(uSF2020_restor_n_avoiddeforest2[]==1 & RDPF2020[]==1, 0, uSF2020_restor_n_avoiddeforest2[])
#plot(uSF2020_restor_n_avoiddeforest2)
writeRaster(uSF2020_restor_n_avoiddeforest2, "rasters/PGM/input/LULC/uSF2020_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)
#uSF2020_restor_n_avoiddeforest2 <- raster("rasters/PGM/input/LULC/uSF2020_restor_n_avoiddeforest2.tif")

DSF2020_restor_n_avoiddeforest2 <- DSF2020
DSF2020_restor_n_avoiddeforest2[] <- ifelse(DSF2020[]==1 & UPF2010[]==1, 0, DSF2020_restor_n_avoiddeforest2[])
DSF2020_restor_n_avoiddeforest2[] <- ifelse(DSF2020[]==1 & uDPF2010[]==1, 0, DSF2020_restor_n_avoiddeforest2[])
DSF2020_restor_n_avoiddeforest2[] <- ifelse(DSF2020[]==1 & RDPF2010[]==1, 0, DSF2020_restor_n_avoiddeforest2[])

DSF2020_restor_n_avoiddeforest2[] <-  ifelse(DSF2020_restor_n_avoiddeforest2[]==1 & uDPF2020[]==1, 0, DSF2020_restor_n_avoiddeforest2[])
DSF2020_restor_n_avoiddeforest2[] <-  ifelse(DSF2020_restor_n_avoiddeforest2[]==1 & RDPF2020[]==1, 0, DSF2020_restor_n_avoiddeforest2[])

#plot(DSF2020_restor_n_avoiddeforest2)
writeRaster(DSF2020_restor_n_avoiddeforest2, "rasters/PGM/input/LULC/DSF2020_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)
#DSF2020_restor_n_avoiddeforest2 <- raster("rasters/PGM/input/LULC/DSF2020_restor_n_avoiddeforest2.tif")


###Undegraded primary forest
UPF2020_restor_n_avoiddeforest2 <- UPF2020
UPF2020_restor_n_avoiddeforest2[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, UPF2020_restor_n_avoiddeforest2[])
UPF2020_restor_n_avoiddeforest2[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, UPF2020_restor_n_avoiddeforest2[])
UPF2020_restor_n_avoiddeforest2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, UPF2020_restor_n_avoiddeforest2[])
#plot(UPF2020_restor_n_avoiddeforest2)
writeRaster(UPF2020_restor_n_avoiddeforest2, "rasters/PGM/input/LULC/UPF2020_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)
#UPF2020_restor_n_avoiddeforest2 <- raster("rasters/PGM/input/LULC/UPF2020_restor_n_avoiddeforest2.tif")


###Degraded primary forest
uDPF2020_restor_n_avoiddeforest2 <- uDPF2020
uDPF2020_restor_n_avoiddeforest2[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, uDPF2020_restor_n_avoiddeforest2[])
uDPF2020_restor_n_avoiddeforest2[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, uDPF2020_restor_n_avoiddeforest2[])
uDPF2020_restor_n_avoiddeforest2[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, uDPF2020_restor_n_avoiddeforest2[])
#plot(uDPF2020_restor_n_avoiddeforest2)
writeRaster(uDPF2020_restor_n_avoiddeforest2, "rasters/PGM/input/LULC/uDPF2020_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)
#uDPF2020_restor_n_avoiddeforest2 <- raster("rasters/PGM/input/LULC/uDPF2020_restor_n_avoiddeforest2.tif")

RDPF2020_restor_n_avoiddeforest2 <- RDPF2020
RDPF2020_restor_n_avoiddeforest2[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, RDPF2020_restor_n_avoiddeforest2[])
RDPF2020_restor_n_avoiddeforest2[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, RDPF2020_restor_n_avoiddeforest2[])
RDPF2020_restor_n_avoiddeforest2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, RDPF2020_restor_n_avoiddeforest2[])
#plot(RDPF2020_restor_n_avoiddeforest2)
writeRaster(RDPF2020_restor_n_avoiddeforest2, "rasters/PGM/input/LULC/RDPF2020_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)
#RDPF2020_restor_n_avoiddeforest2 <- raster("rasters/PGM/input/LULC/RDPF2020_restor_n_avoiddeforest2.tif")

#isolating areas that change
LULC2020_restor_n_avoiddeforest2 <- candidate.areas.final
LULC2020_restor_n_avoiddeforest2[] <- ifelse(LULC2020_restor_n_avoiddeforest2[]==0, NA, LULC2020_restor_n_avoiddeforest2[])
LULC2020_restor_n_avoiddeforest2[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, LULC2020_restor_n_avoiddeforest2[])
LULC2020_restor_n_avoiddeforest2[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, LULC2020_restor_n_avoiddeforest2[])
LULC2020_restor_n_avoiddeforest2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, LULC2020_restor_n_avoiddeforest2[])
LULC2020_restor_n_avoiddeforest2[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, LULC2020_restor_n_avoiddeforest2[])
LULC2020_restor_n_avoiddeforest2[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, LULC2020_restor_n_avoiddeforest2[])
LULC2020_restor_n_avoiddeforest2[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, LULC2020_restor_n_avoiddeforest2[])
LULC2020_restor_n_avoiddeforest2[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, LULC2020_restor_n_avoiddeforest2[])
LULC2020_restor_n_avoiddeforest2[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, LULC2020_restor_n_avoiddeforest2[])
LULC2020_restor_n_avoiddeforest2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, LULC2020_restor_n_avoiddeforest2[])
#length(which(LULC2020_restor_n_avoiddeforest2[]==1))
writeRaster(LULC2020_restor_n_avoiddeforest2, "rasters/PGM/input/LULC/area_change_2020_restor_n_avoiddeforest2_bin.tif", format="GTiff", overwrite=T)
#LULC2020_restor_n_avoiddeforest2.bin <- raster("rasters/PGM/input/LULC/area_change_2020_restor_n_avoiddeforest2_bin.tif")


LULC2020_restor_n_avoiddeforest2[] <- ifelse(LULC2020_restor_n_avoiddeforest2[]==1 & uDPF2020_restor_n_avoiddeforest2[]==1, 10, LULC2020_restor_n_avoiddeforest2[])
LULC2020_restor_n_avoiddeforest2[] <- ifelse(LULC2020_restor_n_avoiddeforest2[]==1 & RDPF2020_restor_n_avoiddeforest2[]==1, 25, LULC2020_restor_n_avoiddeforest2[])
LULC2020_restor_n_avoiddeforest2[] <- ifelse(LULC2020_restor_n_avoiddeforest2[]==1 & uSF2020_restor_n_avoiddeforest2[]==1, 100, LULC2020_restor_n_avoiddeforest2[])
LULC2020_restor_n_avoiddeforest2[] <- ifelse(LULC2020_restor_n_avoiddeforest2[]==1 & DSF2020_restor_n_avoiddeforest2[]==1, 125, LULC2020_restor_n_avoiddeforest2[])

writeRaster(LULC2020_restor_n_avoiddeforest2, "rasters/PGM/input/LULC/area_change_2020_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)
#LULC2020_restor_n_avoiddeforest2 <- raster("rasters/PGM/input/LULC/area_change_2020_restor_n_avoiddeforest2.tif")



##2020 restoration and avoid both (all)
###Secondary forest
uSF2020_restor_n_avoidboth <- sum(uSF2020, candidate.areas.final)
uSF2020_restor_n_avoidboth[] <-  ifelse(uSF2020_restor_n_avoidboth[]!=0, 1, uSF2020_restor_n_avoidboth[])

uSF2020_restor_n_avoidboth[] <- ifelse(uSF2010[]==1 & DSF2020[]==1, 1, uSF2020_restor_n_avoidboth[])
uSF2020_restor_n_avoidboth[] <- ifelse(uSF2020[]==1 & UPF2010[]==1, 0, uSF2020_restor_n_avoidboth[])
uSF2020_restor_n_avoidboth[] <- ifelse(uSF2020[]==1 & uDPF2010[]==1, 0, uSF2020_restor_n_avoidboth[])
uSF2020_restor_n_avoidboth[] <- ifelse(uSF2020[]==1 & RDPF2010[]==1, 0, uSF2020_restor_n_avoidboth[])
uSF2020_restor_n_avoidboth[] <- ifelse(DSF2020[]==0 & uSF2020[]==0 & uSF2010[]==1, 1, uSF2020_restor_n_avoidboth[])

uSF2020_restor_n_avoidboth[] <-  ifelse(uSF2020_restor_n_avoidboth[]==1 & DSF2020[]==1, 0, uSF2020_restor_n_avoidboth[])
uSF2020_restor_n_avoidboth[] <-  ifelse(uSF2020_restor_n_avoidboth[]==1 & UPF2020[]==1, 0, uSF2020_restor_n_avoidboth[])
uSF2020_restor_n_avoidboth[] <-  ifelse(uSF2020_restor_n_avoidboth[]==1 & uDPF2020[]==1, 0, uSF2020_restor_n_avoidboth[])
uSF2020_restor_n_avoidboth[] <-  ifelse(uSF2020_restor_n_avoidboth[]==1 & RDPF2020[]==1, 0, uSF2020_restor_n_avoidboth[])
#plot(uSF2020_restor_n_avoidboth)
writeRaster(uSF2020_restor_n_avoidboth, "rasters/PGM/input/LULC/uSF2020_restor_n_avoidboth.tif", format="GTiff", overwrite=T)
#uSF2020_restor_n_avoidboth <- raster("rasters/PGM/input/LULC/uSF2020_restor_n_avoidboth.tif")

DSF2020_restor_n_avoidboth <- DSF2020
DSF2020_restor_n_avoidboth[] <- ifelse(DSF2020[]==1 & uSF2010[]==1, 0, DSF2020_restor_n_avoidboth[])
DSF2020_restor_n_avoidboth[] <- ifelse(DSF2020[]==1 & UPF2010[]==1, 0, DSF2020_restor_n_avoidboth[])
DSF2020_restor_n_avoidboth[] <- ifelse(DSF2020[]==1 & uDPF2010[]==1, 0, DSF2020_restor_n_avoidboth[])
DSF2020_restor_n_avoidboth[] <- ifelse(DSF2020[]==1 & RDPF2010[]==1, 0, DSF2020_restor_n_avoidboth[])
DSF2020_restor_n_avoidboth[] <- ifelse(uSF2020[]==0 & DSF2020[]==0 & DSF2010[]==1, 1, DSF2020_restor_n_avoidboth[])

DSF2020_restor_n_avoidboth[] <-  ifelse(DSF2020_restor_n_avoidboth[]==1 & uDPF2020[]==1, 0, DSF2020_restor_n_avoidboth[])
DSF2020_restor_n_avoidboth[] <-  ifelse(DSF2020_restor_n_avoidboth[]==1 & RDPF2020[]==1, 0, DSF2020_restor_n_avoidboth[])

#plot(DSF2020_restor_n_avoidboth)
writeRaster(DSF2020_restor_n_avoidboth, "rasters/PGM/input/LULC/DSF2020_restor_n_avoidboth.tif", format="GTiff", overwrite=T)
#DSF2020_restor_n_avoidboth <- raster("rasters/PGM/input/LULC/DSF2020_restor_n_avoidboth.tif")


###Undegraded primary forest
UPF2020_restor_n_avoidboth <- UPF2020
UPF2020_restor_n_avoidboth[] <- ifelse(UPF2010[]==1 & uDPF2020[]==1, 1, UPF2020_restor_n_avoidboth[])
UPF2020_restor_n_avoidboth[] <- ifelse(UPF2010[]==1 & RDPF2020[]==1, 1, UPF2020_restor_n_avoidboth[])
UPF2020_restor_n_avoidboth[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, UPF2020_restor_n_avoidboth[])
UPF2020_restor_n_avoidboth[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, UPF2020_restor_n_avoidboth[])
UPF2020_restor_n_avoidboth[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, UPF2020_restor_n_avoidboth[])
#plot(UPF2020_restor_n_avoidboth)
writeRaster(UPF2020_restor_n_avoidboth, "rasters/PGM/input/LULC/UPF2020_restor_n_avoidboth.tif", format="GTiff", overwrite=T)
#UPF2020_restor_n_avoidboth <- raster("rasters/PGM/input/LULC/UPF2020_restor_n_avoidboth.tif")


###Degraded primary forest
uDPF2020_restor_n_avoidboth <- uDPF2020
uDPF2020_restor_n_avoidboth[] <- ifelse(uDPF2020[]==1 & UPF2010[]==1, 0, uDPF2020_restor_n_avoidboth[])
uDPF2020_restor_n_avoidboth[] <- ifelse(uDPF2010[]==1 & RDPF2020[]==1, 1, uDPF2020_restor_n_avoidboth[])
uDPF2020_restor_n_avoidboth[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, uDPF2020_restor_n_avoidboth[])
uDPF2020_restor_n_avoidboth[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, uDPF2020_restor_n_avoidboth[])
uDPF2020_restor_n_avoidboth[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, uDPF2020_restor_n_avoidboth[])
#plot(uDPF2020_restor_n_avoidboth)
writeRaster(uDPF2020_restor_n_avoidboth, "rasters/PGM/input/LULC/uDPF2020_restor_n_avoidboth.tif", format="GTiff", overwrite=T)
#uDPF2020_restor_n_avoidboth <- raster("rasters/PGM/input/LULC/uDPF2020_restor_n_avoidboth.tif")

RDPF2020_restor_n_avoidboth <- RDPF2020
RDPF2020_restor_n_avoidboth[] <- ifelse(RDPF2020[]==1 & UPF2010[]==1, 0, RDPF2020_restor_n_avoidboth[])
RDPF2020_restor_n_avoidboth[] <- ifelse(RDPF2020[]==1 & uDPF2010[]==1, 0, RDPF2020_restor_n_avoidboth[])
RDPF2020_restor_n_avoidboth[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, RDPF2020_restor_n_avoidboth[])
RDPF2020_restor_n_avoidboth[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, RDPF2020_restor_n_avoidboth[])
RDPF2020_restor_n_avoidboth[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, RDPF2020_restor_n_avoidboth[])
#plot(RDPF2020_restor_n_avoidboth)
writeRaster(RDPF2020_restor_n_avoidboth, "rasters/PGM/input/LULC/RDPF2020_restor_n_avoidboth.tif", format="GTiff", overwrite=T)
#RDPF2020_restor_n_avoidboth <- raster("rasters/PGM/input/LULC/RDPF2020_restor_n_avoidboth.tif")

#isolating areas that change
LULC2020_restor_n_avoidboth <- candidate.areas.final
LULC2020_restor_n_avoidboth[] <- ifelse(LULC2020_restor_n_avoidboth[]==0, NA, LULC2020_restor_n_avoidboth[])
LULC2020_restor_n_avoidboth[] <- ifelse(UPF2010[]==1 & uDPF2020[]==1, 1, LULC2020_restor_n_avoidboth[])
LULC2020_restor_n_avoidboth[] <- ifelse(UPF2010[]==1 & RDPF2020[]==1, 1, LULC2020_restor_n_avoidboth[])
LULC2020_restor_n_avoidboth[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, LULC2020_restor_n_avoidboth[])
LULC2020_restor_n_avoidboth[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, LULC2020_restor_n_avoidboth[])
LULC2020_restor_n_avoidboth[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, LULC2020_restor_n_avoidboth[])
LULC2020_restor_n_avoidboth[] <- ifelse(uDPF2010[]==1 & RDPF2020[]==1, 1, LULC2020_restor_n_avoidboth[])
LULC2020_restor_n_avoidboth[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, LULC2020_restor_n_avoidboth[])
LULC2020_restor_n_avoidboth[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, LULC2020_restor_n_avoidboth[])
LULC2020_restor_n_avoidboth[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, LULC2020_restor_n_avoidboth[])
LULC2020_restor_n_avoidboth[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, LULC2020_restor_n_avoidboth[])
LULC2020_restor_n_avoidboth[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, LULC2020_restor_n_avoidboth[])
LULC2020_restor_n_avoidboth[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, LULC2020_restor_n_avoidboth[])
LULC2020_restor_n_avoidboth[] <- ifelse(uSF2010[]==1 & DSF2020[]==1, 1, LULC2020_restor_n_avoidboth[])
LULC2020_restor_n_avoidboth[] <- ifelse(DSF2020[]==0 & uSF2020[]==0 & uSF2010[]==1, 1, LULC2020_restor_n_avoidboth[])
LULC2020_restor_n_avoidboth[] <- ifelse(uSF2020[]==0 & DSF2020[]==0 & DSF2010[]==1, 1, LULC2020_restor_n_avoidboth[])
#length(which(LULC2020_restor_n_avoidboth[]==1))
writeRaster(LULC2020_restor_n_avoidboth, "rasters/PGM/input/LULC/area_change_2020_restor_n_avoidboth_bin.tif", format="GTiff", overwrite=T)
#LULC2020_restor_n_avoidboth.bin <- raster("rasters/PGM/input/LULC/area_change_2020_restor_n_avoidboth_bin.tif")


LULC2020_restor_n_avoidboth[] <- ifelse(LULC2020_restor_n_avoidboth[]==1 & uDPF2020_restor_n_avoidboth[]==1, 10, LULC2020_restor_n_avoidboth[])
LULC2020_restor_n_avoidboth[] <- ifelse(LULC2020_restor_n_avoidboth[]==1 & RDPF2020_restor_n_avoidboth[]==1, 25, LULC2020_restor_n_avoidboth[])
LULC2020_restor_n_avoidboth[] <- ifelse(LULC2020_restor_n_avoidboth[]==1 & uSF2020_restor_n_avoidboth[]==1, 100, LULC2020_restor_n_avoidboth[])
LULC2020_restor_n_avoidboth[] <- ifelse(LULC2020_restor_n_avoidboth[]==1 & DSF2020_restor_n_avoidboth[]==1, 125, LULC2020_restor_n_avoidboth[])

writeRaster(LULC2020_restor_n_avoidboth, "rasters/PGM/input/LULC/area_change_2020_restor_n_avoidboth.tif", format="GTiff", overwrite=T)
#LULC2020_restor_n_avoidboth <- raster("rasters/PGM/input/LULC/area_change_2020_restor_n_avoidboth.tif")



##2020 restoration and avoid both (Primary forest only)
###Secondary forest
uSF2020_restor_n_avoidboth2 <- sum(uSF2020, candidate.areas.final)
uSF2020_restor_n_avoidboth2[] <-  ifelse(uSF2020_restor_n_avoidboth2[]!=0, 1, uSF2020_restor_n_avoidboth2[])

uSF2020_restor_n_avoidboth2[] <- ifelse(uSF2020[]==1 & UPF2010[]==1, 0, uSF2020_restor_n_avoidboth2[])
uSF2020_restor_n_avoidboth2[] <- ifelse(uSF2020[]==1 & uDPF2010[]==1, 0, uSF2020_restor_n_avoidboth2[])
uSF2020_restor_n_avoidboth2[] <- ifelse(uSF2020[]==1 & RDPF2010[]==1, 0, uSF2020_restor_n_avoidboth2[])

uSF2020_restor_n_avoidboth2[] <-  ifelse(uSF2020_restor_n_avoidboth2[]==1 & DSF2020[]==1, 0, uSF2020_restor_n_avoidboth2[])
uSF2020_restor_n_avoidboth2[] <-  ifelse(uSF2020_restor_n_avoidboth2[]==1 & UPF2020[]==1, 0, uSF2020_restor_n_avoidboth2[])
uSF2020_restor_n_avoidboth2[] <-  ifelse(uSF2020_restor_n_avoidboth2[]==1 & uDPF2020[]==1, 0, uSF2020_restor_n_avoidboth2[])
uSF2020_restor_n_avoidboth2[] <-  ifelse(uSF2020_restor_n_avoidboth2[]==1 & RDPF2020[]==1, 0, uSF2020_restor_n_avoidboth2[])
#plot(uSF2020_restor_n_avoidboth2)
writeRaster(uSF2020_restor_n_avoidboth2, "rasters/PGM/input/LULC/uSF2020_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)
#uSF2020_restor_n_avoidboth2 <- raster("rasters/PGM/input/LULC/uSF2020_restor_n_avoidboth2.tif")

DSF2020_restor_n_avoidboth2 <- DSF2020
DSF2020_restor_n_avoidboth2[] <- ifelse(DSF2020[]==1 & UPF2010[]==1, 0, DSF2020_restor_n_avoidboth2[])
DSF2020_restor_n_avoidboth2[] <- ifelse(DSF2020[]==1 & uDPF2010[]==1, 0, DSF2020_restor_n_avoidboth2[])
DSF2020_restor_n_avoidboth2[] <- ifelse(DSF2020[]==1 & RDPF2010[]==1, 0, DSF2020_restor_n_avoidboth2[])

DSF2020_restor_n_avoidboth2[] <-  ifelse(DSF2020_restor_n_avoidboth2[]==1 & uDPF2020[]==1, 0, DSF2020_restor_n_avoidboth2[])
DSF2020_restor_n_avoidboth2[] <-  ifelse(DSF2020_restor_n_avoidboth2[]==1 & RDPF2020[]==1, 0, DSF2020_restor_n_avoidboth2[])

#plot(DSF2020_restor_n_avoidboth2)
writeRaster(DSF2020_restor_n_avoidboth2, "rasters/PGM/input/LULC/DSF2020_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)
#DSF2020_restor_n_avoidboth2 <- raster("rasters/PGM/input/LULC/DSF2020_restor_n_avoidboth2.tif")


###Undegraded primary forest
UPF2020_restor_n_avoidboth2 <- UPF2020
UPF2020_restor_n_avoidboth2[] <- ifelse(UPF2010[]==1 & uDPF2020[]==1, 1, UPF2020_restor_n_avoidboth2[])
UPF2020_restor_n_avoidboth2[] <- ifelse(UPF2010[]==1 & RDPF2020[]==1, 1, UPF2020_restor_n_avoidboth2[])
UPF2020_restor_n_avoidboth2[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, UPF2020_restor_n_avoidboth2[])
UPF2020_restor_n_avoidboth2[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, UPF2020_restor_n_avoidboth2[])
UPF2020_restor_n_avoidboth2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, UPF2020_restor_n_avoidboth2[])
#plot(UPF2020_restor_n_avoidboth2)
writeRaster(UPF2020_restor_n_avoidboth2, "rasters/PGM/input/LULC/UPF2020_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)
#UPF2020_restor_n_avoidboth2 <- raster("rasters/PGM/input/LULC/UPF2020_restor_n_avoidboth2.tif")


###Degraded primary forest
uDPF2020_restor_n_avoidboth2 <- uDPF2020
uDPF2020_restor_n_avoidboth2[] <- ifelse(uDPF2020[]==1 & UPF2010[]==1, 0, uDPF2020_restor_n_avoidboth2[])
uDPF2020_restor_n_avoidboth2[] <- ifelse(uDPF2010[]==1 & RDPF2020[]==1, 1, uDPF2020_restor_n_avoidboth2[])
uDPF2020_restor_n_avoidboth2[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, uDPF2020_restor_n_avoidboth2[])
uDPF2020_restor_n_avoidboth2[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, uDPF2020_restor_n_avoidboth2[])
uDPF2020_restor_n_avoidboth2[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, uDPF2020_restor_n_avoidboth2[])
#plot(uDPF2020_restor_n_avoidboth2)
writeRaster(uDPF2020_restor_n_avoidboth2, "rasters/PGM/input/LULC/uDPF2020_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)
#uDPF2020_restor_n_avoidboth2 <- raster("rasters/PGM/input/LULC/uDPF2020_restor_n_avoidboth2.tif")

RDPF2020_restor_n_avoidboth2 <- RDPF2020
RDPF2020_restor_n_avoidboth2[] <- ifelse(RDPF2020[]==1 & UPF2010[]==1, 0, RDPF2020_restor_n_avoidboth2[])
RDPF2020_restor_n_avoidboth2[] <- ifelse(RDPF2020[]==1 & uDPF2010[]==1, 0, RDPF2020_restor_n_avoidboth2[])
RDPF2020_restor_n_avoidboth2[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, RDPF2020_restor_n_avoidboth2[])
RDPF2020_restor_n_avoidboth2[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, RDPF2020_restor_n_avoidboth2[])
RDPF2020_restor_n_avoidboth2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, RDPF2020_restor_n_avoidboth2[])
#plot(RDPF2020_restor_n_avoidboth2)
writeRaster(RDPF2020_restor_n_avoidboth2, "rasters/PGM/input/LULC/RDPF2020_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)
#RDPF2020_restor_n_avoidboth2 <- raster("rasters/PGM/input/LULC/RDPF2020_restor_n_avoidboth2.tif")

#isolating areas that change
LULC2020_restor_n_avoidboth2 <- candidate.areas.final
LULC2020_restor_n_avoidboth2[] <- ifelse(LULC2020_restor_n_avoidboth2[]==0, NA, LULC2020_restor_n_avoidboth2[])
LULC2020_restor_n_avoidboth2[] <- ifelse(UPF2010[]==1 & uDPF2020[]==1, 1, LULC2020_restor_n_avoidboth2[])
LULC2020_restor_n_avoidboth2[] <- ifelse(UPF2010[]==1 & RDPF2020[]==1, 1, LULC2020_restor_n_avoidboth2[])
LULC2020_restor_n_avoidboth2[] <- ifelse(UPF2010[]==1 & DSF2020[]==1, 1, LULC2020_restor_n_avoidboth2[])
LULC2020_restor_n_avoidboth2[] <- ifelse(UPF2010[]==1 & uSF2020[]==1, 1, LULC2020_restor_n_avoidboth2[])
LULC2020_restor_n_avoidboth2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & UPF2020[]==0 & UPF2010[]==1, 1, LULC2020_restor_n_avoidboth2[])
LULC2020_restor_n_avoidboth2[] <- ifelse(uDPF2010[]==1 & RDPF2020[]==1, 1, LULC2020_restor_n_avoidboth2[])
LULC2020_restor_n_avoidboth2[] <- ifelse(uDPF2010[]==1 & DSF2020[]==1, 1, LULC2020_restor_n_avoidboth2[])
LULC2020_restor_n_avoidboth2[] <- ifelse(uDPF2010[]==1 & uSF2020[]==1, 1, LULC2020_restor_n_avoidboth2[])
LULC2020_restor_n_avoidboth2[] <- ifelse(RDPF2020[]==0 & uDPF2020[]==0 & uDPF2010[]==1, 1, LULC2020_restor_n_avoidboth2[])
LULC2020_restor_n_avoidboth2[] <- ifelse(RDPF2010[]==1 & DSF2020[]==1, 1, LULC2020_restor_n_avoidboth2[])
LULC2020_restor_n_avoidboth2[] <- ifelse(RDPF2010[]==1 & uSF2020[]==1, 1, LULC2020_restor_n_avoidboth2[])
LULC2020_restor_n_avoidboth2[] <- ifelse(uDPF2020[]==0 & RDPF2020[]==0 & RDPF2010[]==1, 1, LULC2020_restor_n_avoidboth2[])
#length(which(LULC2020_restor_n_avoidboth2[]==1))
writeRaster(LULC2020_restor_n_avoidboth2, "rasters/PGM/input/LULC/area_change_2020_restor_n_avoidboth2_bin.tif", format="GTiff", overwrite=T)
#LULC2020_restor_n_avoidboth2.bin <- raster("rasters/PGM/input/LULC/area_change_2020_restor_n_avoidboth2_bin.tif")


LULC2020_restor_n_avoidboth2[] <- ifelse(LULC2020_restor_n_avoidboth2[]==1 & uDPF2020_restor_n_avoidboth2[]==1, 10, LULC2020_restor_n_avoidboth2[])
LULC2020_restor_n_avoidboth2[] <- ifelse(LULC2020_restor_n_avoidboth2[]==1 & RDPF2020_restor_n_avoidboth2[]==1, 25, LULC2020_restor_n_avoidboth2[])
LULC2020_restor_n_avoidboth2[] <- ifelse(LULC2020_restor_n_avoidboth2[]==1 & uSF2020_restor_n_avoidboth2[]==1, 100, LULC2020_restor_n_avoidboth2[])
LULC2020_restor_n_avoidboth2[] <- ifelse(LULC2020_restor_n_avoidboth2[]==1 & DSF2020_restor_n_avoidboth2[]==1, 125, LULC2020_restor_n_avoidboth2[])

writeRaster(LULC2020_restor_n_avoidboth2, "rasters/PGM/input/LULC/area_change_2020_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)
#LULC2020_restor_n_avoidboth2 <- raster("rasters/PGM/input/LULC/area_change_2020_restor_n_avoidboth2.tif")




















###checking
#undegraded primary forest == 1
UPF2020_restor_n_avoidboth2.sk <- UPF2020_restor_n_avoidboth2
#UPF2020_restor_n_avoidboth2.sk <- raster("rasters/PGM/input/UPF2020_restor_n_avoidboth2.tif")
#UPF2020_restor_n_avoidboth2.sk[UPF2020_restor_n_avoidboth2.sk==1]<-1
UPF2020_restor_n_avoidboth2.sk <- mask(UPF2020_restor_n_avoidboth2.sk, pgm.shp)
UPF2020_restor_n_avoidboth2.sk[UPF2020_restor_n_avoidboth2.sk[]==0] <- 333
#plot(UPF2020_restor_n_avoidboth2.sk, main="undegraded primary forest", legend=F)

#degraded primary forest == 10
uDPF2020_restor_n_avoidboth2.sk <- uDPF2020_restor_n_avoidboth2
#DPF2020_restor_n_avoidboth2.sk <- raster("rasters/PGM/input/uDPF2020_restor_n_avoidboth2.tif")
uDPF2020_restor_n_avoidboth2.sk[uDPF2020_restor_n_avoidboth2.sk==1]<-10
uDPF2020_restor_n_avoidboth2.sk <- mask(uDPF2020_restor_n_avoidboth2.sk, pgm.shp)
uDPF2020_restor_n_avoidboth2.sk[uDPF2020_restor_n_avoidboth2.sk[]==0] <- 333
#plot(uDPF2020_restor_n_avoidboth2.sk, main="degraded primary forest", legend=F)

#repeated degraded primary forest == 25
RDPF2020_restor_n_avoidboth2.sk <- RDPF2020_restor_n_avoidboth2
#RDPF2020_restor_n_avoidboth2.sk <- raster("rasters/PGM/input/RDPF2020_restor_n_avoidboth2.tif")
RDPF2020_restor_n_avoidboth2.sk[RDPF2020_restor_n_avoidboth2.sk==1]<-25
RDPF2020_restor_n_avoidboth2.sk <- mask(RDPF2020_restor_n_avoidboth2.sk, pgm.shp)
RDPF2020_restor_n_avoidboth2.sk[RDPF2020_restor_n_avoidboth2.sk[]==0] <- 333
#plot(RDPF2020_restor_n_avoidboth2.sk, main="degraded primary forest", legend=F)

#secondary forest == 100
uSF2020_restor_n_avoidboth2.sk <- uSF2020_restor_n_avoidboth2
#uSF2020_restor_n_avoidboth2.sk <- raster("rasters/PGM/input/uSF2020_restor_n_avoidboth2.tif")
uSF2020_restor_n_avoidboth2.sk[uSF2020_restor_n_avoidboth2.sk==1]<-100
uSF2020_restor_n_avoidboth2.sk <- mask(uSF2020_restor_n_avoidboth2.sk, pgm.shp)
uSF2020_restor_n_avoidboth2.sk[uSF2020_restor_n_avoidboth2.sk[]==0] <- 333
#plot(uSF2020_restor_n_avoidboth2.sk, main="secondary forest", legend=F)

#degraded secondary forest == 125
DSF2020_restor_n_avoidboth2.sk <- DSF2020_restor_n_avoidboth2
#DSF2020_restor_n_avoidboth2.sk <- raster("rasters/PGM/input/DSF2020_restor_n_avoidboth2.tif")
DSF2020_restor_n_avoidboth2.sk[DSF2020_restor_n_avoidboth2.sk==1]<-125
DSF2020_restor_n_avoidboth2.sk <- mask(DSF2020_restor_n_avoidboth2.sk, pgm.shp)
DSF2020_restor_n_avoidboth2.sk[DSF2020_restor_n_avoidboth2.sk[]==0] <- 333
#plot(DSF2020_restor_n_avoidboth2.sk, main="secondary forest", legend=F)

LULC2020_restor_n_avoidboth2 <- sum(UPF2020_restor_n_avoidboth2.sk, uSF2020_restor_n_avoidboth2.sk, na.rm = T)
LULC2020_restor_n_avoidboth2 <- sum(LULC2020_restor_n_avoidboth2, DSF2020_restor_n_avoidboth2.sk, na.rm = T)
LULC2020_restor_n_avoidboth2 <- sum(LULC2020_restor_n_avoidboth2, uDPF2020_restor_n_avoidboth2.sk, na.rm = T)
LULC2020_restor_n_avoidboth2 <- sum(LULC2020_restor_n_avoidboth2, RDPF2020_restor_n_avoidboth2.sk, na.rm = T)
LULC2020_restor_n_avoidboth2 <- mask(LULC2020_restor_n_avoidboth2, pgm.shp)
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
dir.create("rasters/PGM/2010_real", recursive = T)
dir.create("rasters/PGM/2020_real", recursive = T)
dir.create("rasters/PGM/2020_avoiddeforest", recursive = T)
dir.create("rasters/PGM/2020_avoiddeforest2", recursive = T) #PF only
dir.create("rasters/PGM/2020_avoiddegrad", recursive = T)
dir.create("rasters/PGM/2020_avoiddegrad2", recursive = T) #PF only
dir.create("rasters/PGM/2020_avoidboth", recursive = T)
dir.create("rasters/PGM/2020_avoidboth2", recursive = T) #PF only
dir.create("rasters/PGM/2020_restor_wo_avoid", recursive = T)
dir.create("rasters/PGM/2020_restor_n_avoiddeforest", recursive = T)
dir.create("rasters/PGM/2020_restor_n_avoiddeforest2", recursive = T) #PF only
dir.create("rasters/PGM/2020_restor_n_avoidboth", recursive = T)
dir.create("rasters/PGM/2020_restor_n_avoidboth2", recursive = T) #PF only
dir.create("models.output/costs", recursive = T)



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
# dir.create("rasters/PGM/input/climate")


### scenario 2010
temp.list <- list.files("rasters/PGM/raw/climate", "LSTD", full.names = T, recursive = T)

temp2010.list <- grep("2010", temp.list, value = T)
temp2010 <- stack(temp2010.list)
#plot(temp2010)

meantemp2010 <- mean(temp2010, na.rm=T)
pgm.meantemp2010 <- crop(meantemp2010, extent(pgm.lulc[[1]]))
#plot(pgm.meantemp2010)
#plot(pgm.shp, add=T)

pgm.meantemp2010 <- resample(pgm.meantemp2010, pgm.lulc[[1]], method='bilinear')
pgm.meantemp2010[] <- ifelse(pgm.lulc[[1]][]==0, NA, pgm.meantemp2010[])
#plot(pgm.meantemp2010)

#saving
writeRaster(pgm.meantemp2010, "rasters/PGM/2010_real/meantemps.tif", format="GTiff", overwrite=T)
#

### scenario 2020
temp2020.list <- grep("2020", temp.list, value = T)
temp2020 <- stack(temp2020.list)
#plot(temp2020)

meantemp2020 <- mean(temp2020, na.rm=T)
pgm.meantemp2020 <- crop(meantemp2020, extent(pgm.lulc[[2]]))
#plot(pgm.meantemp2020)
#plot(pgm.shp, add=T)

pgm.meantemp2020 <- resample(pgm.meantemp2020, pgm.lulc[[2]], method='bilinear')
pgm.meantemp2020[] <- ifelse(pgm.lulc[[2]][]==0, NA, pgm.meantemp2020[])
#plot(pgm.meantemp2020)

#saving
writeRaster(pgm.meantemp2020, "rasters/PGM/2020_real/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meantemp2020, "rasters/PGM/2020_avoiddeforest/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meantemp2020, "rasters/PGM/2020_avoiddeforest2/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meantemp2020, "rasters/PGM/2020_avoiddegrad/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meantemp2020, "rasters/PGM/2020_avoiddegrad2/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meantemp2020, "rasters/PGM/2020_avoidboth/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meantemp2020, "rasters/PGM/2020_avoidboth2/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meantemp2020, "rasters/PGM/2020_restor_wo_avoid/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meantemp2020, "rasters/PGM/2020_restor_n_avoiddeforest/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meantemp2020, "rasters/PGM/2020_restor_n_avoiddeforest2/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meantemp2020, "rasters/PGM/2020_restor_n_avoidboth/meantemps.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meantemp2020, "rasters/PGM/2020_restor_n_avoidboth2/meantemps.tif", format="GTiff", overwrite=T)
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
precip.list <- list.files("rasters/PGM/raw/climate", "GPM", full.names = T, recursive = T)

precip2010.list <- grep("2010", precip.list, value = T)
precip2010 <- stack(precip2010.list)
#plot(precip2010)

meanprecip2010 <- mean(precip2010, na.rm=T)
pgm.meanprecip2010 <- crop(meanprecip2010, extent(pgm.lulc[[1]]))
#plot(pgm.meanprecip2010)
#plot(pgm.shp, add=T)

pgm.meanprecip2010 <- resample(pgm.meanprecip2010, pgm.lulc[[1]], method='bilinear')
pgm.meanprecip2010[] <- ifelse(pgm.lulc[[1]][]==0, NA, pgm.meanprecip2010[])
#plot(pgm.meanprecip2010)

#saving
writeRaster(pgm.meanprecip2010, "rasters/PGM/2010_real/meanprecips.tif", format="GTiff", overwrite=T)
#

### scenario 2020
precip2020.list <- grep("2020", precip.list, value = T)
precip2020 <- stack(precip2020.list)
#plot(precip2020)


meanprecip2020 <- mean(precip2020, na.rm=T)
pgm.meanprecip2020 <- crop(meanprecip2020, extent(pgm.lulc[[2]]))
#plot(pgm.meanprecip2020)
#plot(pgm.shp, add=T)

pgm.meanprecip2020 <- resample(pgm.meanprecip2020, pgm.lulc[[2]], method='bilinear')
pgm.meanprecip2020[] <- ifelse(pgm.lulc[[2]][]==0, NA, pgm.meanprecip2020[])
#plot(pgm.meanprecip2020)

#saving
writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_real/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_avoiddeforest/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_avoiddeforest2/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_avoiddegrad/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_avoiddegrad2/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_avoidboth/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_avoidboth2/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_restor_wo_avoid/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_restor_n_avoiddeforest/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_restor_n_avoiddeforest2/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_restor_n_avoidboth/meanprecips.tif", format="GTiff", overwrite=T)
writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_restor_n_avoidboth2/meanprecips.tif", format="GTiff", overwrite=T)
#

rm(list=ls()[ls() %in% c("precip.list", "precip2010.list", "meanprecip2010", "precip2020.list", "meanprecip2020")])
gc()



#
#




## [elevation]
elevation <- raster("rasters/PGM/raw/pgm-elevation.tif")
elevation <- projectRaster(elevation, crs = std.proj)
elevation <- resample(elevation, pgm.lulc[[1]], method='bilinear')
elevation[] <- ifelse(pgm.lulc[[1]][]==0, NA, elevation[])

#saving
writeRaster(elevation, "rasters/PGM/2010_real/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/PGM/2020_real/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/PGM/2020_avoiddeforest/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/PGM/2020_avoiddeforest2/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/PGM/2020_avoiddegrad/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/PGM/2020_avoiddegrad2/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/PGM/2020_avoidboth/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/PGM/2020_avoidboth2/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/PGM/2020_restor_wo_avoid/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/PGM/2020_restor_n_avoiddeforest/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/PGM/2020_restor_n_avoiddeforest2/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/PGM/2020_restor_n_avoidboth/elevation.tif", format="GTiff", overwrite=T)
writeRaster(elevation, "rasters/PGM/2020_restor_n_avoidboth2/elevation.tif", format="GTiff", overwrite=T)



#
#




## [distroad] distance to roads
dist.road <- raster("rasters/PGM/raw/pgm-distance-to-road.tif")
dist.road <- projectRaster(dist.road, crs = std.proj)
dist.road <- resample(dist.road, pgm.lulc[[1]], method='bilinear')
dist.road[] <- ifelse(pgm.lulc[[1]][]==0, NA, dist.road[])

#saving
writeRaster(dist.road, "rasters/PGM/2010_real/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/PGM/2020_real/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/PGM/2020_avoiddeforest/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/PGM/2020_avoiddeforest2/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/PGM/2020_avoiddegrad/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/PGM/2020_avoiddegrad2/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/PGM/2020_avoidboth/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/PGM/2020_avoidboth2/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/PGM/2020_restor_wo_avoid/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/PGM/2020_restor_n_avoiddeforest/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/PGM/2020_restor_n_avoiddeforest2/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/PGM/2020_restor_n_avoidboth/distroad.tif", format="GTiff", overwrite=T)
writeRaster(dist.road, "rasters/PGM/2020_restor_n_avoidboth2/distroad.tif", format="GTiff", overwrite=T)



#
#




## Rivers
#' source: HydroSHEDS - https://www.hydrosheds.org/hydroatlas
pgm.river <- readOGR(dsn = "rasters/PGM/raw", layer = "pgm_RiverATLAS_v10")
#head(pgm.river@data)

### selecting variables and creating river width attribute
pgm.river@data <- pgm.river@data %>% dplyr::select(HYRIV_ID:LENGTH_KM, ria_ha_csu, ria_ha_usu) %>% 
  mutate(
    ril_m = LENGTH_KM * 1000, #converting to meters
    ria_m2 = ria_ha_csu * 10000, #converting to square meters
    riw_m = (ria_m2 / ril_m) #estimating the mean river segment width
    )

#checking
#st_crs(pgm.river)==st_crs(pgm.shp)
#plot(pgm.lulc[["pgm.lulc.2010real"]])
#plot(pgm.river, add=T)


### [distriver] distance to rivers
pgm.river.raster <- rasterize(pgm.river, pgm.lulc[[1]], field = 1, background = 0)
pgm.river.raster.mask <- pgm.river.raster
pgm.river.raster.mask[pgm.river.raster.mask==0] <- NA
#checking
#inv.pgm.river
#plot(pgm.river.raster.mask, col="black")

pgm.dist.river <- distance(pgm.river.raster.mask, doEdge=T)
#writeRaster(pgm.dist.river, "rasters/PGM/raw/pgm-distance-to-river.tif", format = "GTiff", overwrite = T)
#pgm.dist.river <- raster("rasters/PGM/raw/pgm-distance-to-river.tif")

pgm.dist.river[] <- ifelse(pgm.lulc[[1]][]==0, NA, pgm.dist.river[])
##cheking
#anyNA(pgm.dist.river[])
#plot(pgm.dist.river)

#saving
writeRaster(pgm.dist.river, "rasters/PGM/2010_real/distriver.tif", format="GTiff", overwrite=T)
writeRaster(pgm.dist.river, "rasters/PGM/2020_real/distriver.tif", format="GTiff", overwrite=T)
writeRaster(pgm.dist.river, "rasters/PGM/2020_avoiddeforest/distriver.tif", format="GTiff", overwrite=T)
writeRaster(pgm.dist.river, "rasters/PGM/2020_avoiddeforest2/distriver.tif", format="GTiff", overwrite=T)
writeRaster(pgm.dist.river, "rasters/PGM/2020_avoiddegrad/distriver.tif", format="GTiff", overwrite=T)
writeRaster(pgm.dist.river, "rasters/PGM/2020_avoiddegrad2/distriver.tif", format="GTiff", overwrite=T)
writeRaster(pgm.dist.river, "rasters/PGM/2020_avoidboth/distriver.tif", format="GTiff", overwrite=T)
writeRaster(pgm.dist.river, "rasters/PGM/2020_avoidboth2/distriver.tif", format="GTiff", overwrite=T)
writeRaster(pgm.dist.river, "rasters/PGM/2020_restor_wo_avoid/distriver.tif", format="GTiff", overwrite=T)
writeRaster(pgm.dist.river, "rasters/PGM/2020_restor_n_avoiddeforest/distriver.tif", format="GTiff", overwrite=T)
writeRaster(pgm.dist.river, "rasters/PGM/2020_restor_n_avoiddeforest2/distriver.tif", format="GTiff", overwrite=T)
writeRaster(pgm.dist.river, "rasters/PGM/2020_restor_n_avoidboth/distriver.tif", format="GTiff", overwrite=T)
writeRaster(pgm.dist.river, "rasters/PGM/2020_restor_n_avoidboth2/distriver.tif", format="GTiff", overwrite=T)
#

rm(list=ls()[ls() %in% c("pgm.river.raster", "pgm.river.raster.mask")])
gc()



#
#




## [distmarket] distance to municipality nucleus
pgm.munic.nucleus <- data.frame(ID = "pgm", long = -47.35311, lat = -3.00249)

pgm.munic.nucleus.coord <- SpatialPointsDataFrame(coords = pgm.munic.nucleus[,c("long","lat")], 
                                                  data = pgm.munic.nucleus, 
                                                  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
pgm.munic.nucleus.coord <- spTransform(pgm.munic.nucleus.coord, crs(std.proj))


distmarket <- distanceFromPoints(object = pgm.lulc[[1]], xy = pgm.munic.nucleus.coord)
distmarket[] <- ifelse(pgm.lulc[[1]][]==0, NA, distmarket[])

#saving
writeRaster(distmarket, "rasters/PGM/2010_real/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/PGM/2020_real/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/PGM/2020_avoiddeforest/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/PGM/2020_avoiddeforest2/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/PGM/2020_avoiddegrad/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/PGM/2020_avoiddegrad2/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/PGM/2020_avoidboth/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/PGM/2020_avoidboth2/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/PGM/2020_restor_wo_avoid/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/PGM/2020_restor_n_avoiddeforest/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/PGM/2020_restor_n_avoiddeforest2/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/PGM/2020_restor_n_avoidboth/distmarket.tif", format="GTiff", overwrite=T)
writeRaster(distmarket, "rasters/PGM/2020_restor_n_avoidboth2/distmarket.tif", format="GTiff", overwrite=T)
#

rm(list=ls()[ls() %in% c("pgm.munic.nucleus", "pgm.munic.nucleus.coord")])
gc()



#
#



## Scenario: 2010 Real =========================================================
### Undegraded primary forest
#UPF2010 <- raster("rasters/PGM/input/LULC/UPF2010_real.tif")

### mean upf cover in local scale (90m)
UPF2010.px <- focal(UPF2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2010.px)<-"UPFpx"
UPF2010.px[is.nan(UPF2010.px)] <- 0
UPF2010.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2010.px[])
#saving
writeRaster(UPF2010.px, "rasters/PGM/2010_real/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2010.ls <- focal(UPF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2010.ls)<-"UPFls"
UPF2010.ls[is.nan(UPF2010.ls)] <- 0
UPF2010.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2010.ls[])
#saving
writeRaster(UPF2010.ls, "rasters/PGM/2010_real/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2010.px", "UPF2010.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2010 <- raster("rasters/PGM/input/LULC/uDPF2010_real.tif")
#RDPF2010 <- raster("rasters/PGM/input/LULC/RDPF2010_real.tif")
DPF2010 <- sum(uDPF2010, RDPF2010)

### mean dpf cover in local scale (90m)
DPF2010.px <- focal(DPF2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2010.px)<-"DPFpx"
DPF2010.px[is.nan(DPF2010.px)] <- 0
DPF2010.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2010.px[])
#saving
writeRaster(DPF2010.px, "rasters/PGM/2010_real/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2010.ls <- focal(DPF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2010.ls)<-"DPFls"
DPF2010.ls[is.nan(DPF2010.ls)] <- 0
DPF2010.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2010.ls[])
#saving
writeRaster(DPF2010.ls, "rasters/PGM/2010_real/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2010.px", "DPF2010.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2010.pf <- pgm.degrad.pf[["pgm.degrad.2010real"]]
TSD2010.pf[is.na(TSD2010.pf)] <- 0

TSD2010.sf <- pgm.degrad.sf[["pgm.degradsf.2010real"]]
TSD2010.sf[is.na(TSD2010.sf)] <- 0

TSD2010 <- sum(TSD2010.pf, TSD2010.sf)
writeRaster(TSD2010, "rasters/PGM/input/TSD2010.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2010.px <- focal(TSD2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2010.px)<-"TSDpx"
TSD2010.px[is.nan(TSD2010.px)] <- 0
TSD2010.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2010.px[])
#saving
writeRaster(TSD2010.px, "rasters/PGM/2010_real/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2010.ls <- focal(TSD2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2010.ls)<-"TSDls"
TSD2010.ls[is.nan(TSD2010.ls)] <- 0
TSD2010.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2010.ls[])
#saving
writeRaster(TSD2010.ls, "rasters/PGM/2010_real/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2010.px", "TSD2010.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2010 <- raster("rasters/PGM/input/LULC/uSF2010_real.tif")
#DSF2010 <- raster("rasters/PGM/input/LULC/DSF2010_real.tif")
SF2010 <- sum(uSF2010, DSF2010)

### mean sf cover in local scale (90m)
SF2010.px <- focal(SF2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2010.px)<-"SFpx"
SF2010.px[is.nan(SF2010.px)] <- 0
SF2010.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2010.px[])
#saving
writeRaster(SF2010.px, "rasters/PGM/2010_real/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2010.ls <- focal(SF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2010.ls)<-"SFls"
SF2010.ls[is.nan(SF2010.ls)] <- 0
SF2010.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2010.ls[])
#saving
writeRaster(SF2010.ls, "rasters/PGM/2010_real/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2010.px", "SF2010.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2010 <- pgm.sfage[["pgm.sfage.2010real"]]
writeRaster(SFAge2010, "rasters/PGM/input/SFAge2010.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2010.px <- focal(SFAge2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2010.px)<-"SFAgepx"
SFAge2010.px[is.nan(SFAge2010.px)] <- 0
SFAge2010.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2010.px[])
#saving
writeRaster(SFAge2010.px, "rasters/PGM/2010_real/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2010.ls <- focal(SFAge2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2010.ls)<-"SFAgels"
SFAge2010.ls[is.nan(SFAge2010.ls)] <- 0
SFAge2010.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2010.ls[])
#saving
writeRaster(SFAge2010.ls, "rasters/PGM/2010_real/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2010.px", "SFAge2010.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2010young <- SFAge2010
SF2010young[] <- ifelse(SF2010young[]>2, 1, 0)
writeRaster(SF2010young, "rasters/PGM/input/SF2010young.tif", format="GTiff", overwrite=T)

TF2010 <- sum(UPF2010, DPF2010)
TF2010 <- sum(TF2010, SF2010young)
writeRaster(TF2010, "rasters/PGM/input/TF2010.tif", format="GTiff", overwrite=T)
  
### mean tf cover in local scale (90m)
TF2010.px <- focal(TF2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2010.px)<-"TFpx"
TF2010.px[is.nan(TF2010.px)] <- 0
TF2010.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2010.px[])
#saving
writeRaster(TF2010.px, "rasters/PGM/2010_real/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2010.ls <- focal(TF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2010.ls)<-"TFls"
TF2010.ls[is.nan(TF2010.ls)] <- 0
TF2010.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2010.ls[])
#saving
writeRaster(TF2010.ls, "rasters/PGM/2010_real/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2010young", "TF2010.px", "TF2010.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2010mature <- SFAge2010
SF2010mature[] <- ifelse(SF2010mature[]>5, 1, 0)
writeRaster(SF2010mature, "rasters/PGM/input/SF2010mature.tif", format="GTiff", overwrite=T)

MF2010 <- sum(UPF2010, DPF2010)
MF2010 <- sum(MF2010, SF2010mature)
writeRaster(MF2010, "rasters/PGM/input/MF2010.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2010.px <- focal(MF2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2010.px)<-"MFpx"
MF2010.px[is.nan(MF2010.px)] <- 0
MF2010.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2010.px[])
#saving
writeRaster(MF2010.px, "rasters/PGM/2010_real/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2010.ls <- focal(MF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2010.ls)<-"MFls"
MF2010.ls[is.nan(MF2010.ls)] <- 0
MF2010.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2010.ls[])
#saving
writeRaster(MF2010.ls, "rasters/PGM/2010_real/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2010mature", "MF2010.px", "MF2010.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2010 <- raster("rasters/PGM/input/MF2010.tif")
inv.MF2010 <- MF2010
inv.MF2010[inv.MF2010==1]<-NA
#cheking
#inv.MF2010
#plot(inv.MF2010)

edge.dist.2010 <- distance(inv.MF2010, doEdge=T)
names(edge.dist.2010)<-"edgedist"
edge.dist.2010[is.nan(edge.dist.2010)] <- 0
edge.dist.2010[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge.dist.2010[])
#saving
writeRaster(edge.dist.2010, "rasters/PGM/2010_real/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2010")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2010 <- edge.dist.2010
edge2010[] <- ifelse(edge2010[] < 200, 0, ifelse(edge2010[]>300, 0, 1))
writeRaster(edge2010, "rasters/PGM/input/edge2010_real.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2010.px <- focal(edge2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2010.px)<-"edgepx"
edge2010.px[is.nan(edge2010.px)] <- 0
edge2010.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2010.px[])
#saving
writeRaster(edge2010.px, "rasters/PGM/2010_real/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2010.ls <- focal(edge2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2010.ls)<-"edgels"
edge2010.ls[is.nan(edge2010.ls)] <- 0
edge2010.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2010.ls[])
#saving
writeRaster(edge2010.ls, "rasters/PGM/2010_real/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2010.px", "edge2010.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: 2020 Real =========================================================
### Undegraded primary forest
#UPF2020 <- raster("rasters/PGM/input/LULC/UPF2020_real.tif")

### mean upf cover in local scale (90m)
UPF2020.px <- focal(UPF2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020.px)<-"UPFpx"
UPF2020.px[is.nan(UPF2020.px)] <- 0
UPF2020.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020.px[])
#saving
writeRaster(UPF2020.px, "rasters/PGM/2020_real/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020.ls <- focal(UPF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020.ls)<-"UPFls"
UPF2020.ls[is.nan(UPF2020.ls)] <- 0
UPF2020.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020.ls[])
#saving
writeRaster(UPF2020.ls, "rasters/PGM/2020_real/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020.px", "UPF2020.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020 <- raster("rasters/PGM/input/LULC/uDPF2020_real.tif")
#RDPF2020 <- raster("rasters/PGM/input/LULC/RDPF2020_real.tif")
DPF2020 <- sum(uDPF2020, RDPF2020)

### mean dpf cover in local scale (90m)
DPF2020.px <- focal(DPF2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020.px)<-"DPFpx"
DPF2020.px[is.nan(DPF2020.px)] <- 0
DPF2020.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020.px[])
#saving
writeRaster(DPF2020.px, "rasters/PGM/2020_real/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020.ls <- focal(DPF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020.ls)<-"DPFls"
DPF2020.ls[is.nan(DPF2020.ls)] <- 0
DPF2020.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020.ls[])
#saving
writeRaster(DPF2020.ls, "rasters/PGM/2020_real/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020.px", "DPF2020.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020.pf <- pgm.degrad.pf[["pgm.degrad.2020real"]]
TSD2020.pf[is.na(TSD2020.pf)] <- 0

TSD2020.sf <- pgm.degrad.sf[["pgm.degradsf.2020real"]]
TSD2020.sf[is.na(TSD2020.sf)] <- 0

TSD2020 <- sum(TSD2020.pf, TSD2020.sf)
writeRaster(TSD2020, "rasters/PGM/input/TSD2020.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020.px <- focal(TSD2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020.px)<-"TSDpx"
TSD2020.px[is.nan(TSD2020.px)] <- 0
TSD2020.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020.px[])
#saving
writeRaster(TSD2020.px, "rasters/PGM/2020_real/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020.ls <- focal(TSD2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020.ls)<-"TSDls"
TSD2020.ls[is.nan(TSD2020.ls)] <- 0
TSD2020.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020.ls[])
#saving
writeRaster(TSD2020.ls, "rasters/PGM/2020_real/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020.px", "TSD2020.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020 <- raster("rasters/PGM/input/LULC/uSF2020_real.tif")
#DSF2020 <- raster("rasters/PGM/input/LULC/DSF2020_real.tif")
SF2020 <- sum(uSF2020, DSF2020)

### mean sf cover in local scale (90m)
SF2020.px <- focal(SF2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020.px)<-"SFpx"
SF2020.px[is.nan(SF2020.px)] <- 0
SF2020.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020.px[])
#saving
writeRaster(SF2020.px, "rasters/PGM/2020_real/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020.ls <- focal(SF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020.ls)<-"SFls"
SF2020.ls[is.nan(SF2020.ls)] <- 0
SF2020.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020.ls[])
#saving
writeRaster(SF2020.ls, "rasters/PGM/2020_real/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020.px", "SF2020.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020 <- pgm.sfage[["pgm.sfage.2020real"]]
writeRaster(SFAge2020, "rasters/PGM/input/SFAge2020.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020.px <- focal(SFAge2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020.px)<-"SFAgepx"
SFAge2020.px[is.nan(SFAge2020.px)] <- 0
SFAge2020.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020.px[])
#saving
writeRaster(SFAge2020.px, "rasters/PGM/2020_real/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020.ls <- focal(SFAge2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020.ls)<-"SFAgels"
SFAge2020.ls[is.nan(SFAge2020.ls)] <- 0
SFAge2020.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020.ls[])
#saving
writeRaster(SFAge2020.ls, "rasters/PGM/2020_real/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020.px", "SFAge2020.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young <- SFAge2020
SF2020young[] <- ifelse(SF2020young[]>2, 1, 0)
writeRaster(SF2020young, "rasters/PGM/input/SF2020young.tif", format="GTiff", overwrite=T)

TF2020 <- sum(UPF2020, DPF2020)
TF2020 <- sum(TF2020, SF2020young)
writeRaster(TF2020, "rasters/PGM/input/TF2020.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020.px <- focal(TF2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020.px)<-"TFpx"
TF2020.px[is.nan(TF2020.px)] <- 0
TF2020.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020.px[])
#saving
writeRaster(TF2020.px, "rasters/PGM/2020_real/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020.ls <- focal(TF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020.ls)<-"TFls"
TF2020.ls[is.nan(TF2020.ls)] <- 0
TF2020.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020.ls[])
#saving
writeRaster(TF2020.ls, "rasters/PGM/2020_real/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020young", "TF2020.px", "TF2020.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature <- SFAge2020
SF2020mature[] <- ifelse(SF2020mature[]>5, 1, 0)
writeRaster(SF2020mature, "rasters/PGM/input/SF2020mature.tif", format="GTiff", overwrite=T)

MF2020 <- sum(UPF2020, DPF2020)
MF2020 <- sum(MF2020, SF2020mature)
writeRaster(MF2020, "rasters/PGM/input/MF2020.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020.px <- focal(MF2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020.px)<-"MFpx"
MF2020.px[is.nan(MF2020.px)] <- 0
MF2020.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020.px[])
#saving
writeRaster(MF2020.px, "rasters/PGM/2020_real/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020.ls <- focal(MF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020.ls)<-"MFls"
MF2020.ls[is.nan(MF2020.ls)] <- 0
MF2020.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020.ls[])
#saving
writeRaster(MF2020.ls, "rasters/PGM/2020_real/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020mature", "MF2020.px", "MF2020.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020 <- raster("rasters/PGM/input/MF2020.tif")
inv.MF2020 <- MF2020
inv.MF2020[inv.MF2020==1]<-NA
#cheking
#inv.MF2020
#plot(inv.MF2020)

edge.dist.2020 <- distance(inv.MF2020, doEdge=T)
names(edge.dist.2020)<-"edgedist"
edge.dist.2020[is.nan(edge.dist.2020)] <- 0
edge.dist.2020[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge.dist.2020[])
#saving
writeRaster(edge.dist.2020, "rasters/PGM/2020_real/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020 <- edge.dist.2020
edge2020[] <- ifelse(edge2020[] < 200, 0, ifelse(edge2020[]>300, 0, 1))
writeRaster(edge2020, "rasters/PGM/input/edge2020_real.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020.px <- focal(edge2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020.px)<-"edgepx"
edge2020.px[is.nan(edge2020.px)] <- 0
edge2020.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020.px[])
#saving
writeRaster(edge2020.px, "rasters/PGM/2020_real/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020.ls <- focal(edge2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020.ls)<-"edgels"
edge2020.ls[is.nan(edge2020.ls)] <- 0
edge2020.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020.ls[])
#saving
writeRaster(edge2020.ls, "rasters/PGM/2020_real/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020.px", "edge2020.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Avoid degradation (all) ===========================================
### Undegraded primary forest
#UPF2020_avoiddegrad <- raster("rasters/PGM/input/LULC/UPF2020_avoiddegrad.tif")

### mean upf cover in local scale (90m)
UPF2020_avoiddegrad.px <- focal(UPF2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_avoiddegrad.px)<-"UPFpx"
UPF2020_avoiddegrad.px[is.nan(UPF2020_avoiddegrad.px)] <- 0
UPF2020_avoiddegrad.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_avoiddegrad.px[])
#saving
writeRaster(UPF2020_avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_avoiddegrad.ls <- focal(UPF2020_avoiddegrad, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_avoiddegrad.ls)<-"UPFls"
UPF2020_avoiddegrad.ls[is.nan(UPF2020_avoiddegrad.ls)] <- 0
UPF2020_avoiddegrad.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_avoiddegrad.ls[])
#saving
writeRaster(UPF2020_avoiddegrad.ls, "rasters/PGM/2020_avoiddegrad/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_avoiddegrad.px", "UPF2020_avoiddegrad.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_avoiddegrad <- raster("rasters/PGM/input/LULC/uDPF2020_avoiddegrad.tif")
#RDPF2020_avoiddegrad <- raster("rasters/PGM/input/LULC/RDPF2020_avoiddegrad.tif")
DPF2020_avoiddegrad <- sum(uDPF2020_avoiddegrad, RDPF2020_avoiddegrad)

### mean dpf cover in local scale (90m)
DPF2020_avoiddegrad.px <- focal(DPF2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_avoiddegrad.px)<-"DPFpx"
DPF2020_avoiddegrad.px[is.nan(DPF2020_avoiddegrad.px)] <- 0
DPF2020_avoiddegrad.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_avoiddegrad.px[])
#saving
writeRaster(DPF2020_avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_avoiddegrad.ls <- focal(DPF2020_avoiddegrad, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_avoiddegrad.ls)<-"DPFls"
DPF2020_avoiddegrad.ls[is.nan(DPF2020_avoiddegrad.ls)] <- 0
DPF2020_avoiddegrad.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_avoiddegrad.ls[])
#saving
writeRaster(DPF2020_avoiddegrad.ls, "rasters/PGM/2020_avoiddegrad/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_avoiddegrad.px", "DPF2020_avoiddegrad.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_avoiddegrad.pf <- calc(pgm.degrad.pf[["pgm.degrad.2010real"]], fun=function(x){ifelse(is.na(x), x, ifelse(x==0, x, ifelse(x==300, x, x+10)))})
TSD2020_avoiddegrad.pf[is.na(TSD2020_avoiddegrad.pf)] <- 0

TSD2020_avoiddegrad.sf <- calc(pgm.degrad.sf[["pgm.degradsf.2010real"]], fun=function(x){ifelse(is.na(x), x, ifelse(x==0, x, ifelse(x==300, x, x+10)))})
TSD2020_avoiddegrad.sf[is.na(TSD2020_avoiddegrad.sf)] <- 0

TSD2020_avoiddegrad <- sum(TSD2020_avoiddegrad.pf, TSD2020_avoiddegrad.sf)
writeRaster(TSD2020_avoiddegrad, "rasters/PGM/input/TSD2020_avoiddegrad.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_avoiddegrad.px <- focal(TSD2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_avoiddegrad.px)<-"TSDpx"
TSD2020_avoiddegrad.px[is.nan(TSD2020_avoiddegrad.px)] <- 0
TSD2020_avoiddegrad.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_avoiddegrad.px[])
#saving
writeRaster(TSD2020_avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_avoiddegrad.ls <- focal(TSD2020_avoiddegrad, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_avoiddegrad.ls)<-"TSDls"
TSD2020_avoiddegrad.ls[is.nan(TSD2020_avoiddegrad.ls)] <- 0
TSD2020_avoiddegrad.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_avoiddegrad.ls[])
#saving
writeRaster(TSD2020_avoiddegrad.ls, "rasters/PGM/2020_avoiddegrad/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_avoiddegrad.px", "TSD2020_avoiddegrad.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_avoiddegrad <- raster("rasters/PGM/input/LULC/uSF2020_avoiddegrad.tif")
#DSF2020_avoiddegrad <- raster("rasters/PGM/input/LULC/DSF2020_avoiddegrad.tif")
SF2020_avoiddegrad <- sum(uSF2020_avoiddegrad, DSF2020_avoiddegrad)

### mean sf cover in local scale (90m)
SF2020_avoiddegrad.px <- focal(SF2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_avoiddegrad.px)<-"SFpx"
SF2020_avoiddegrad.px[is.nan(SF2020_avoiddegrad.px)] <- 0
SF2020_avoiddegrad.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_avoiddegrad.px[])
#saving
writeRaster(SF2020_avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_avoiddegrad.ls <- focal(SF2020_avoiddegrad, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_avoiddegrad.ls)<-"SFls"
SF2020_avoiddegrad.ls[is.nan(SF2020_avoiddegrad.ls)] <- 0
SF2020_avoiddegrad.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_avoiddegrad.ls[])
#saving
writeRaster(SF2020_avoiddegrad.ls, "rasters/PGM/2020_avoiddegrad/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddegrad.px", "SF2020_avoiddegrad.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_avoiddegrad <- pgm.sfage[["pgm.sfage.2020real"]]
writeRaster(SFAge2020_avoiddegrad, "rasters/PGM/input/SFAge2020_avoiddegrad.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_avoiddegrad.px <- focal(SFAge2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_avoiddegrad.px)<-"SFAgepx"
SFAge2020_avoiddegrad.px[is.nan(SFAge2020_avoiddegrad.px)] <- 0
SFAge2020_avoiddegrad.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_avoiddegrad.px[])
#saving
writeRaster(SFAge2020_avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_avoiddegrad.ls <- focal(SFAge2020_avoiddegrad, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_avoiddegrad.ls)<-"SFAgels"
SFAge2020_avoiddegrad.ls[is.nan(SFAge2020_avoiddegrad.ls)] <- 0
SFAge2020_avoiddegrad.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_avoiddegrad.ls[])
#saving
writeRaster(SFAge2020_avoiddegrad.ls, "rasters/PGM/2020_avoiddegrad/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_avoiddegrad.px", "SFAge2020_avoiddegrad.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young_avoiddegrad <- SFAge2020_avoiddegrad
SF2020young_avoiddegrad[] <- ifelse(SF2020young_avoiddegrad[]>2, 1, 0)
writeRaster(SF2020young_avoiddegrad, "rasters/PGM/input/SF2020young_avoiddegrad.tif", format="GTiff", overwrite=T)

TF2020_avoiddegrad <- sum(UPF2020_avoiddegrad, DPF2020_avoiddegrad)
TF2020_avoiddegrad <- sum(TF2020_avoiddegrad, SF2020young_avoiddegrad)
writeRaster(TF2020_avoiddegrad, "rasters/PGM/input/TF2020_avoiddegrad.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_avoiddegrad.px <- focal(TF2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_avoiddegrad.px)<-"TFpx"
TF2020_avoiddegrad.px[is.nan(TF2020_avoiddegrad.px)] <- 0
TF2020_avoiddegrad.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_avoiddegrad.px[])
#saving
writeRaster(TF2020_avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_avoiddegrad.ls <- focal(TF2020_avoiddegrad, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_avoiddegrad.ls)<-"TFls"
TF2020_avoiddegrad.ls[is.nan(TF2020_avoiddegrad.ls)] <- 0
TF2020_avoiddegrad.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_avoiddegrad.ls[])
#saving
writeRaster(TF2020_avoiddegrad.ls, "rasters/PGM/2020_avoiddegrad/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddegradyoung", "TF2020_avoiddegrad.px", "TF2020_avoiddegrad.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature_avoiddegrad <- SFAge2020_avoiddegrad
SF2020mature_avoiddegrad[] <- ifelse(SF2020mature_avoiddegrad[]>5, 1, 0)
writeRaster(SF2020mature_avoiddegrad, "rasters/PGM/input/SF2020mature_avoiddegrad.tif", format="GTiff", overwrite=T)

MF2020_avoiddegrad <- sum(UPF2020_avoiddegrad, DPF2020_avoiddegrad)
MF2020_avoiddegrad <- sum(MF2020_avoiddegrad, SF2020mature_avoiddegrad)
writeRaster(MF2020_avoiddegrad, "rasters/PGM/input/MF2020_avoiddegrad.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_avoiddegrad.px <- focal(MF2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_avoiddegrad.px)<-"MFpx"
MF2020_avoiddegrad.px[is.nan(MF2020_avoiddegrad.px)] <- 0
MF2020_avoiddegrad.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_avoiddegrad.px[])
#saving
writeRaster(MF2020_avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_avoiddegrad.ls <- focal(MF2020_avoiddegrad, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_avoiddegrad.ls)<-"MFls"
MF2020_avoiddegrad.ls[is.nan(MF2020_avoiddegrad.ls)] <- 0
MF2020_avoiddegrad.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_avoiddegrad.ls[])
#saving
writeRaster(MF2020_avoiddegrad.ls, "rasters/PGM/2020_avoiddegrad/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddegradmature", "MF2020_avoiddegrad.px", "MF2020_avoiddegrad.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_avoiddegrad <- raster("rasters/PGM/input/MF2020_avoiddegrad.tif")
inv.MF2020_avoiddegrad <- MF2020_avoiddegrad
inv.MF2020_avoiddegrad[inv.MF2020_avoiddegrad==1]<-NA
#cheking
#inv.MF2020_avoiddegrad
#plot(inv.MF2020_avoiddegrad)

edge.dist.2020_avoiddegrad <- distance(inv.MF2020_avoiddegrad, doEdge=T)
names(edge.dist.2020_avoiddegrad)<-"edgedist"
edge.dist.2020_avoiddegrad[is.nan(edge.dist.2020_avoiddegrad)] <- 0
edge.dist.2020_avoiddegrad[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge.dist.2020_avoiddegrad[])
#saving
writeRaster(edge.dist.2020_avoiddegrad, "rasters/PGM/2020_avoiddegrad/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_avoiddegrad")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_avoiddegrad <- edge.dist.2020_avoiddegrad
edge2020_avoiddegrad[] <- ifelse(edge2020_avoiddegrad[] < 200, 0, ifelse(edge2020_avoiddegrad[]>300, 0, 1))
writeRaster(edge2020_avoiddegrad, "rasters/PGM/input/edge2020_avoiddegrad.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_avoiddegrad.px <- focal(edge2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_avoiddegrad.px)<-"edgepx"
edge2020_avoiddegrad.px[is.nan(edge2020_avoiddegrad.px)] <- 0
edge2020_avoiddegrad.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_avoiddegrad.px[])
#saving
writeRaster(edge2020_avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_avoiddegrad.ls <- focal(edge2020_avoiddegrad, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_avoiddegrad.ls)<-"edgels"
edge2020_avoiddegrad.ls[is.nan(edge2020_avoiddegrad.ls)] <- 0
edge2020_avoiddegrad.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_avoiddegrad.ls[])
#saving
writeRaster(edge2020_avoiddegrad.ls, "rasters/PGM/2020_avoiddegrad/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_avoiddegrad.px", "edge2020_avoiddegrad.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Avoid degradation (Primary forests only) ==========================
### Undegraded primary forest
#UPF2020_avoiddegrad2 <- raster("rasters/PGM/input/LULC/UPF2020_avoiddegrad2.tif")

### mean upf cover in local scale (90m)
UPF2020_avoiddegrad2.px <- focal(UPF2020_avoiddegrad2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_avoiddegrad2.px)<-"UPFpx"
UPF2020_avoiddegrad2.px[is.nan(UPF2020_avoiddegrad2.px)] <- 0
UPF2020_avoiddegrad2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_avoiddegrad2.px[])
#saving
writeRaster(UPF2020_avoiddegrad2.px, "rasters/PGM/2020_avoiddegrad2/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_avoiddegrad2.ls <- focal(UPF2020_avoiddegrad2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_avoiddegrad2.ls)<-"UPFls"
UPF2020_avoiddegrad2.ls[is.nan(UPF2020_avoiddegrad2.ls)] <- 0
UPF2020_avoiddegrad2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_avoiddegrad2.ls[])
#saving
writeRaster(UPF2020_avoiddegrad2.ls, "rasters/PGM/2020_avoiddegrad2/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_avoiddegrad2.px", "UPF2020_avoiddegrad2.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_avoiddegrad2 <- raster("rasters/PGM/input/LULC/uDPF2020_avoiddegrad2.tif")
#RDPF2020_avoiddegrad2 <- raster("rasters/PGM/input/LULC/RDPF2020_avoiddegrad2.tif")
DPF2020_avoiddegrad2 <- sum(uDPF2020_avoiddegrad2, RDPF2020_avoiddegrad2)

### mean dpf cover in local scale (90m)
DPF2020_avoiddegrad2.px <- focal(DPF2020_avoiddegrad2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_avoiddegrad2.px)<-"DPFpx"
DPF2020_avoiddegrad2.px[is.nan(DPF2020_avoiddegrad2.px)] <- 0
DPF2020_avoiddegrad2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_avoiddegrad2.px[])
#saving
writeRaster(DPF2020_avoiddegrad2.px, "rasters/PGM/2020_avoiddegrad2/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_avoiddegrad2.ls <- focal(DPF2020_avoiddegrad2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_avoiddegrad2.ls)<-"DPFls"
DPF2020_avoiddegrad2.ls[is.nan(DPF2020_avoiddegrad2.ls)] <- 0
DPF2020_avoiddegrad2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_avoiddegrad2.ls[])
#saving
writeRaster(DPF2020_avoiddegrad2.ls, "rasters/PGM/2020_avoiddegrad2/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_avoiddegrad2.px", "DPF2020_avoiddegrad2.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_avoiddegrad2.pf <- calc(pgm.degrad.pf[["pgm.degrad.2010real"]], fun=function(x){ifelse(is.na(x), x, ifelse(x==0, x, ifelse(x==300, x, x+10)))})
TSD2020_avoiddegrad2.pf[is.na(TSD2020_avoiddegrad2.pf)] <- 0

TSD2020_avoiddegrad2.sf <- pgm.degrad.sf[["pgm.degradsf.2020real"]]
TSD2020_avoiddegrad2.sf[is.na(TSD2020_avoiddegrad2.sf)] <- 0

TSD2020_avoiddegrad2 <- sum(TSD2020_avoiddegrad2.pf, TSD2020_avoiddegrad2.sf)
TSD2020_avoiddegrad2[TSD2020_avoiddegrad2>300] <- 0
writeRaster(TSD2020_avoiddegrad2, "rasters/PGM/input/TSD2020_avoiddegrad2.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_avoiddegrad2.px <- focal(TSD2020_avoiddegrad2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_avoiddegrad2.px)<-"TSDpx"
TSD2020_avoiddegrad2.px[is.nan(TSD2020_avoiddegrad2.px)] <- 0
TSD2020_avoiddegrad2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_avoiddegrad2.px[])
#saving
writeRaster(TSD2020_avoiddegrad2.px, "rasters/PGM/2020_avoiddegrad2/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_avoiddegrad2.ls <- focal(TSD2020_avoiddegrad2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_avoiddegrad2.ls)<-"TSDls"
TSD2020_avoiddegrad2.ls[is.nan(TSD2020_avoiddegrad2.ls)] <- 0
TSD2020_avoiddegrad2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_avoiddegrad2.ls[])
#saving
writeRaster(TSD2020_avoiddegrad2.ls, "rasters/PGM/2020_avoiddegrad2/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_avoiddegrad2.px", "TSD2020_avoiddegrad2.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_avoiddegrad2 <- raster("rasters/PGM/input/LULC/uSF2020_avoiddegrad2.tif")
#DSF2020_avoiddegrad2 <- raster("rasters/PGM/input/LULC/DSF2020_avoiddegrad2.tif")
SF2020_avoiddegrad2 <- sum(uSF2020_avoiddegrad2, DSF2020_avoiddegrad2)

### mean sf cover in local scale (90m)
SF2020_avoiddegrad2.px <- focal(SF2020_avoiddegrad2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_avoiddegrad2.px)<-"SFpx"
SF2020_avoiddegrad2.px[is.nan(SF2020_avoiddegrad2.px)] <- 0
SF2020_avoiddegrad2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_avoiddegrad2.px[])
#saving
writeRaster(SF2020_avoiddegrad2.px, "rasters/PGM/2020_avoiddegrad2/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_avoiddegrad2.ls <- focal(SF2020_avoiddegrad2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_avoiddegrad2.ls)<-"SFls"
SF2020_avoiddegrad2.ls[is.nan(SF2020_avoiddegrad2.ls)] <- 0
SF2020_avoiddegrad2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_avoiddegrad2.ls[])
#saving
writeRaster(SF2020_avoiddegrad2.ls, "rasters/PGM/2020_avoiddegrad2/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddegrad2.px", "SF2020_avoiddegrad2.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_avoiddegrad2 <- pgm.sfage[["pgm.sfage.2020real"]]
writeRaster(SFAge2020_avoiddegrad2, "rasters/PGM/input/SFAge2020_avoiddegrad2.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_avoiddegrad2.px <- focal(SFAge2020_avoiddegrad2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_avoiddegrad2.px)<-"SFAgepx"
SFAge2020_avoiddegrad2.px[is.nan(SFAge2020_avoiddegrad2.px)] <- 0
SFAge2020_avoiddegrad2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_avoiddegrad2.px[])
#saving
writeRaster(SFAge2020_avoiddegrad2.px, "rasters/PGM/2020_avoiddegrad2/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_avoiddegrad2.ls <- focal(SFAge2020_avoiddegrad2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_avoiddegrad2.ls)<-"SFAgels"
SFAge2020_avoiddegrad2.ls[is.nan(SFAge2020_avoiddegrad2.ls)] <- 0
SFAge2020_avoiddegrad2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_avoiddegrad2.ls[])
#saving
writeRaster(SFAge2020_avoiddegrad2.ls, "rasters/PGM/2020_avoiddegrad2/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_avoiddegrad2.px", "SFAge2020_avoiddegrad2.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young_avoiddegrad2 <- SFAge2020_avoiddegrad2
SF2020young_avoiddegrad2[] <- ifelse(SF2020young_avoiddegrad2[]>2, 1, 0)
writeRaster(SF2020young_avoiddegrad2, "rasters/PGM/input/SF2020young_avoiddegrad2.tif", format="GTiff", overwrite=T)

TF2020_avoiddegrad2 <- sum(UPF2020_avoiddegrad2, DPF2020_avoiddegrad2)
TF2020_avoiddegrad2 <- sum(TF2020_avoiddegrad2, SF2020young_avoiddegrad2)
writeRaster(TF2020_avoiddegrad2, "rasters/PGM/input/TF2020_avoiddegrad2.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_avoiddegrad2.px <- focal(TF2020_avoiddegrad2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_avoiddegrad2.px)<-"TFpx"
TF2020_avoiddegrad2.px[is.nan(TF2020_avoiddegrad2.px)] <- 0
TF2020_avoiddegrad2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_avoiddegrad2.px[])
#saving
writeRaster(TF2020_avoiddegrad2.px, "rasters/PGM/2020_avoiddegrad2/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_avoiddegrad2.ls <- focal(TF2020_avoiddegrad2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_avoiddegrad2.ls)<-"TFls"
TF2020_avoiddegrad2.ls[is.nan(TF2020_avoiddegrad2.ls)] <- 0
TF2020_avoiddegrad2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_avoiddegrad2.ls[])
#saving
writeRaster(TF2020_avoiddegrad2.ls, "rasters/PGM/2020_avoiddegrad2/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddegrad2young", "TF2020_avoiddegrad2.px", "TF2020_avoiddegrad2.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature_avoiddegrad2 <- SFAge2020_avoiddegrad2
SF2020mature_avoiddegrad2[] <- ifelse(SF2020mature_avoiddegrad2[]>5, 1, 0)
writeRaster(SF2020mature_avoiddegrad2, "rasters/PGM/input/SF2020mature_avoiddegrad2.tif", format="GTiff", overwrite=T)

MF2020_avoiddegrad2 <- sum(UPF2020_avoiddegrad2, DPF2020_avoiddegrad2)
MF2020_avoiddegrad2 <- sum(MF2020_avoiddegrad2, SF2020mature_avoiddegrad2)
writeRaster(MF2020_avoiddegrad2, "rasters/PGM/input/MF2020_avoiddegrad2.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_avoiddegrad2.px <- focal(MF2020_avoiddegrad2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_avoiddegrad2.px)<-"MFpx"
MF2020_avoiddegrad2.px[is.nan(MF2020_avoiddegrad2.px)] <- 0
MF2020_avoiddegrad2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_avoiddegrad2.px[])
#saving
writeRaster(MF2020_avoiddegrad2.px, "rasters/PGM/2020_avoiddegrad2/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_avoiddegrad2.ls <- focal(MF2020_avoiddegrad2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_avoiddegrad2.ls)<-"MFls"
MF2020_avoiddegrad2.ls[is.nan(MF2020_avoiddegrad2.ls)] <- 0
MF2020_avoiddegrad2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_avoiddegrad2.ls[])
#saving
writeRaster(MF2020_avoiddegrad2.ls, "rasters/PGM/2020_avoiddegrad2/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddegrad2mature", "MF2020_avoiddegrad2.px", "MF2020_avoiddegrad2.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_avoiddegrad2 <- raster("rasters/PGM/input/MF2020_avoiddegrad2.tif")
inv.MF2020_avoiddegrad2 <- MF2020_avoiddegrad2
inv.MF2020_avoiddegrad2[inv.MF2020_avoiddegrad2==1]<-NA
#cheking
#inv.MF2020_avoiddegrad2
#plot(inv.MF2020_avoiddegrad2)

edge.dist.2020_avoiddegrad2 <- distance(inv.MF2020_avoiddegrad2, doEdge=T)
names(edge.dist.2020_avoiddegrad2)<-"edgedist"
edge.dist.2020_avoiddegrad2[is.nan(edge.dist.2020_avoiddegrad2)] <- 0
edge.dist.2020_avoiddegrad2[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge.dist.2020_avoiddegrad2[])
#saving
writeRaster(edge.dist.2020_avoiddegrad2, "rasters/PGM/2020_avoiddegrad2/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_avoiddegrad2")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_avoiddegrad2 <- edge.dist.2020_avoiddegrad2
edge2020_avoiddegrad2[] <- ifelse(edge2020_avoiddegrad2[] < 200, 0, ifelse(edge2020_avoiddegrad2[]>300, 0, 1))
writeRaster(edge2020_avoiddegrad2, "rasters/PGM/input/edge2020_avoiddegrad2.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_avoiddegrad2.px <- focal(edge2020_avoiddegrad2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_avoiddegrad2.px)<-"edgepx"
edge2020_avoiddegrad2.px[is.nan(edge2020_avoiddegrad2.px)] <- 0
edge2020_avoiddegrad2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_avoiddegrad2.px[])
#saving
writeRaster(edge2020_avoiddegrad2.px, "rasters/PGM/2020_avoiddegrad2/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_avoiddegrad2.ls <- focal(edge2020_avoiddegrad2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_avoiddegrad2.ls)<-"edgels"
edge2020_avoiddegrad2.ls[is.nan(edge2020_avoiddegrad2.ls)] <- 0
edge2020_avoiddegrad2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_avoiddegrad2.ls[])
#saving
writeRaster(edge2020_avoiddegrad2.ls, "rasters/PGM/2020_avoiddegrad2/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_avoiddegrad2.px", "edge2020_avoiddegrad2.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Avoid deforestation (all) =========================================
### Undegraded primary forest
#UPF2020_avoiddeforest <- raster("rasters/PGM/input/LULC/UPF2020_avoiddeforest.tif")

### mean upf cover in local scale (90m)
UPF2020_avoiddeforest.px <- focal(UPF2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_avoiddeforest.px)<-"UPFpx"
UPF2020_avoiddeforest.px[is.nan(UPF2020_avoiddeforest.px)] <- 0
UPF2020_avoiddeforest.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_avoiddeforest.px[])
#saving
writeRaster(UPF2020_avoiddeforest.px, "rasters/PGM/2020_avoiddeforest/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_avoiddeforest.ls <- focal(UPF2020_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_avoiddeforest.ls)<-"UPFls"
UPF2020_avoiddeforest.ls[is.nan(UPF2020_avoiddeforest.ls)] <- 0
UPF2020_avoiddeforest.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_avoiddeforest.ls[])
#saving
writeRaster(UPF2020_avoiddeforest.ls, "rasters/PGM/2020_avoiddeforest/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_avoiddeforest.px", "UPF2020_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_avoiddeforest <- raster("rasters/PGM/input/LULC/uDPF2020_avoiddeforest.tif")
#RDPF2020_avoiddeforest <- raster("rasters/PGM/input/LULC/RDPF2020_avoiddeforest.tif")
DPF2020_avoiddeforest <- sum(uDPF2020_avoiddeforest, RDPF2020_avoiddeforest)

### mean dpf cover in local scale (90m)
DPF2020_avoiddeforest.px <- focal(DPF2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_avoiddeforest.px)<-"DPFpx"
DPF2020_avoiddeforest.px[is.nan(DPF2020_avoiddeforest.px)] <- 0
DPF2020_avoiddeforest.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_avoiddeforest.px[])
#saving
writeRaster(DPF2020_avoiddeforest.px, "rasters/PGM/2020_avoiddeforest/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_avoiddeforest.ls <- focal(DPF2020_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_avoiddeforest.ls)<-"DPFls"
DPF2020_avoiddeforest.ls[is.nan(DPF2020_avoiddeforest.ls)] <- 0
DPF2020_avoiddeforest.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_avoiddeforest.ls[])
#saving
writeRaster(DPF2020_avoiddeforest.ls, "rasters/PGM/2020_avoiddeforest/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_avoiddeforest.px", "DPF2020_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_avoiddeforest <- TSD2020
writeRaster(TSD2020_avoiddeforest, "rasters/PGM/input/TSD2020_avoiddeforest.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_avoiddeforest.px <- focal(TSD2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_avoiddeforest.px)<-"TSDpx"
TSD2020_avoiddeforest.px[is.nan(TSD2020_avoiddeforest.px)] <- 0
TSD2020_avoiddeforest.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_avoiddeforest.px[])
#saving
writeRaster(TSD2020_avoiddeforest.px, "rasters/PGM/2020_avoiddeforest/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_avoiddeforest.ls <- focal(TSD2020_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_avoiddeforest.ls)<-"TSDls"
TSD2020_avoiddeforest.ls[is.nan(TSD2020_avoiddeforest.ls)] <- 0
TSD2020_avoiddeforest.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_avoiddeforest.ls[])
#saving
writeRaster(TSD2020_avoiddeforest.ls, "rasters/PGM/2020_avoiddeforest/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_avoiddeforest.px", "TSD2020_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_avoiddeforest <- raster("rasters/PGM/input/LULC/uSF2020_avoiddeforest.tif")
#DSF2020_avoiddeforest <- raster("rasters/PGM/input/LULC/DSF2020_avoiddeforest.tif")
SF2020_avoiddeforest <- sum(uSF2020_avoiddeforest, DSF2020_avoiddeforest)

### mean sf cover in local scale (90m)
SF2020_avoiddeforest.px <- focal(SF2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_avoiddeforest.px)<-"SFpx"
SF2020_avoiddeforest.px[is.nan(SF2020_avoiddeforest.px)] <- 0
SF2020_avoiddeforest.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_avoiddeforest.px[])
#saving
writeRaster(SF2020_avoiddeforest.px, "rasters/PGM/2020_avoiddeforest/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_avoiddeforest.ls <- focal(SF2020_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_avoiddeforest.ls)<-"SFls"
SF2020_avoiddeforest.ls[is.nan(SF2020_avoiddeforest.ls)] <- 0
SF2020_avoiddeforest.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_avoiddeforest.ls[])
#saving
writeRaster(SF2020_avoiddeforest.ls, "rasters/PGM/2020_avoiddeforest/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddeforest.px", "SF2020_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_avoiddeforest <- calc(pgm.sfage[["pgm.sfage.2010real"]], fun=function(x){ifelse(is.na(x), x, ifelse(x==0, x, x+10))})
writeRaster(SFAge2020_avoiddeforest, "rasters/PGM/input/SFAge2020_avoiddeforest.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_avoiddeforest.px <- focal(SFAge2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_avoiddeforest.px)<-"SFAgepx"
SFAge2020_avoiddeforest.px[is.nan(SFAge2020_avoiddeforest.px)] <- 0
SFAge2020_avoiddeforest.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_avoiddeforest.px[])
#saving
writeRaster(SFAge2020_avoiddeforest.px, "rasters/PGM/2020_avoiddeforest/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_avoiddeforest.ls <- focal(SFAge2020_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_avoiddeforest.ls)<-"SFAgels"
SFAge2020_avoiddeforest.ls[is.nan(SFAge2020_avoiddeforest.ls)] <- 0
SFAge2020_avoiddeforest.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_avoiddeforest.ls[])
#saving
writeRaster(SFAge2020_avoiddeforest.ls, "rasters/PGM/2020_avoiddeforest/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_avoiddeforest.px", "SFAge2020_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young_avoiddeforest <- SFAge2020_avoiddeforest
SF2020young_avoiddeforest[] <- ifelse(SF2020young_avoiddeforest[]>2, 1, 0)
writeRaster(SF2020young_avoiddeforest, "rasters/PGM/input/SF2020young_avoiddeforest.tif", format="GTiff", overwrite=T)

TF2020_avoiddeforest <- sum(UPF2020_avoiddeforest, DPF2020_avoiddeforest)
TF2020_avoiddeforest <- sum(TF2020_avoiddeforest, SF2020young_avoiddeforest)
writeRaster(TF2020_avoiddeforest, "rasters/PGM/input/TF2020_avoiddeforest.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_avoiddeforest.px <- focal(TF2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_avoiddeforest.px)<-"TFpx"
TF2020_avoiddeforest.px[is.nan(TF2020_avoiddeforest.px)] <- 0
TF2020_avoiddeforest.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_avoiddeforest.px[])
#saving
writeRaster(TF2020_avoiddeforest.px, "rasters/PGM/2020_avoiddeforest/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_avoiddeforest.ls <- focal(TF2020_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_avoiddeforest.ls)<-"TFls"
TF2020_avoiddeforest.ls[is.nan(TF2020_avoiddeforest.ls)] <- 0
TF2020_avoiddeforest.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_avoiddeforest.ls[])
#saving
writeRaster(TF2020_avoiddeforest.ls, "rasters/PGM/2020_avoiddeforest/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddeforestyoung", "TF2020_avoiddeforest.px", "TF2020_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature_avoiddeforest <- SFAge2020_avoiddeforest
SF2020mature_avoiddeforest[] <- ifelse(SF2020mature_avoiddeforest[]>5, 1, 0)
writeRaster(SF2020mature_avoiddeforest, "rasters/PGM/input/SF2020mature_avoiddeforest.tif", format="GTiff", overwrite=T)

MF2020_avoiddeforest <- sum(UPF2020_avoiddeforest, DPF2020_avoiddeforest)
MF2020_avoiddeforest <- sum(MF2020_avoiddeforest, SF2020mature_avoiddeforest)
writeRaster(MF2020_avoiddeforest, "rasters/PGM/input/MF2020_avoiddeforest.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_avoiddeforest.px <- focal(MF2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_avoiddeforest.px)<-"MFpx"
MF2020_avoiddeforest.px[is.nan(MF2020_avoiddeforest.px)] <- 0
MF2020_avoiddeforest.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_avoiddeforest.px[])
#saving
writeRaster(MF2020_avoiddeforest.px, "rasters/PGM/2020_avoiddeforest/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_avoiddeforest.ls <- focal(MF2020_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_avoiddeforest.ls)<-"MFls"
MF2020_avoiddeforest.ls[is.nan(MF2020_avoiddeforest.ls)] <- 0
MF2020_avoiddeforest.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_avoiddeforest.ls[])
#saving
writeRaster(MF2020_avoiddeforest.ls, "rasters/PGM/2020_avoiddeforest/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddeforestmature", "MF2020_avoiddeforest.px", "MF2020_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_avoiddeforest <- raster("rasters/PGM/input/MF2020_avoiddeforest.tif")
inv.MF2020_avoiddeforest <- MF2020_avoiddeforest
inv.MF2020_avoiddeforest[inv.MF2020_avoiddeforest==1]<-NA
#cheking
#inv.MF2020_avoiddeforest
#plot(inv.MF2020_avoiddeforest)

edge.dist.2020_avoiddeforest <- distance(inv.MF2020_avoiddeforest, doEdge=T)
names(edge.dist.2020_avoiddeforest)<-"edgedist"
edge.dist.2020_avoiddeforest[is.nan(edge.dist.2020_avoiddeforest)] <- 0
edge.dist.2020_avoiddeforest[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge.dist.2020_avoiddeforest[])
#saving
writeRaster(edge.dist.2020_avoiddeforest, "rasters/PGM/2020_avoiddeforest/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_avoiddeforest")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_avoiddeforest <- edge.dist.2020_avoiddeforest
edge2020_avoiddeforest[] <- ifelse(edge2020_avoiddeforest[] < 200, 0, ifelse(edge2020_avoiddeforest[]>300, 0, 1))
writeRaster(edge2020_avoiddeforest, "rasters/PGM/input/edge2020_avoiddeforest.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_avoiddeforest.px <- focal(edge2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_avoiddeforest.px)<-"edgepx"
edge2020_avoiddeforest.px[is.nan(edge2020_avoiddeforest.px)] <- 0
edge2020_avoiddeforest.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_avoiddeforest.px[])
#saving
writeRaster(edge2020_avoiddeforest.px, "rasters/PGM/2020_avoiddeforest/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_avoiddeforest.ls <- focal(edge2020_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_avoiddeforest.ls)<-"edgels"
edge2020_avoiddeforest.ls[is.nan(edge2020_avoiddeforest.ls)] <- 0
edge2020_avoiddeforest.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_avoiddeforest.ls[])
#saving
writeRaster(edge2020_avoiddeforest.ls, "rasters/PGM/2020_avoiddeforest/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_avoiddeforest.px", "edge2020_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Avoid deforestation (Primary forests only) ========================
### Undegraded primary forest
#UPF2020_avoiddeforest2 <- raster("rasters/PGM/input/LULC/UPF2020_avoiddeforest2.tif")

### mean upf cover in local scale (90m)
UPF2020_avoiddeforest2.px <- focal(UPF2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_avoiddeforest2.px)<-"UPFpx"
UPF2020_avoiddeforest2.px[is.nan(UPF2020_avoiddeforest2.px)] <- 0
UPF2020_avoiddeforest2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_avoiddeforest2.px[])
#saving
writeRaster(UPF2020_avoiddeforest2.px, "rasters/PGM/2020_avoiddeforest2/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_avoiddeforest2.ls <- focal(UPF2020_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_avoiddeforest2.ls)<-"UPFls"
UPF2020_avoiddeforest2.ls[is.nan(UPF2020_avoiddeforest2.ls)] <- 0
UPF2020_avoiddeforest2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_avoiddeforest2.ls[])
#saving
writeRaster(UPF2020_avoiddeforest2.ls, "rasters/PGM/2020_avoiddeforest2/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_avoiddeforest2.px", "UPF2020_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_avoiddeforest2 <- raster("rasters/PGM/input/LULC/uDPF2020_avoiddeforest2.tif")
#RDPF2020_avoiddeforest2 <- raster("rasters/PGM/input/LULC/RDPF2020_avoiddeforest2.tif")
DPF2020_avoiddeforest2 <- sum(uDPF2020_avoiddeforest2, RDPF2020_avoiddeforest2)

### mean dpf cover in local scale (90m)
DPF2020_avoiddeforest2.px <- focal(DPF2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_avoiddeforest2.px)<-"DPFpx"
DPF2020_avoiddeforest2.px[is.nan(DPF2020_avoiddeforest2.px)] <- 0
DPF2020_avoiddeforest2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_avoiddeforest2.px[])
#saving
writeRaster(DPF2020_avoiddeforest2.px, "rasters/PGM/2020_avoiddeforest2/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_avoiddeforest2.ls <- focal(DPF2020_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_avoiddeforest2.ls)<-"DPFls"
DPF2020_avoiddeforest2.ls[is.nan(DPF2020_avoiddeforest2.ls)] <- 0
DPF2020_avoiddeforest2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_avoiddeforest2.ls[])
#saving
writeRaster(DPF2020_avoiddeforest2.ls, "rasters/PGM/2020_avoiddeforest2/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_avoiddeforest2.px", "DPF2020_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_avoiddeforest2 <- TSD2020
writeRaster(TSD2020_avoiddeforest2, "rasters/PGM/input/TSD2020_avoiddeforest2.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_avoiddeforest2.px <- focal(TSD2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_avoiddeforest2.px)<-"TSDpx"
TSD2020_avoiddeforest2.px[is.nan(TSD2020_avoiddeforest2.px)] <- 0
TSD2020_avoiddeforest2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_avoiddeforest2.px[])
#saving
writeRaster(TSD2020_avoiddeforest2.px, "rasters/PGM/2020_avoiddeforest2/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_avoiddeforest2.ls <- focal(TSD2020_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_avoiddeforest2.ls)<-"TSDls"
TSD2020_avoiddeforest2.ls[is.nan(TSD2020_avoiddeforest2.ls)] <- 0
TSD2020_avoiddeforest2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_avoiddeforest2.ls[])
#saving
writeRaster(TSD2020_avoiddeforest2.ls, "rasters/PGM/2020_avoiddeforest2/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_avoiddeforest2.px", "TSD2020_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_avoiddeforest2 <- raster("rasters/PGM/input/LULC/uSF2020_avoiddeforest2.tif")
#DSF2020_avoiddeforest2 <- raster("rasters/PGM/input/LULC/DSF2020_avoiddeforest2.tif")
SF2020_avoiddeforest2 <- sum(uSF2020_avoiddeforest2, DSF2020_avoiddeforest2)

### mean sf cover in local scale (90m)
SF2020_avoiddeforest2.px <- focal(SF2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_avoiddeforest2.px)<-"SFpx"
SF2020_avoiddeforest2.px[is.nan(SF2020_avoiddeforest2.px)] <- 0
SF2020_avoiddeforest2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_avoiddeforest2.px[])
#saving
writeRaster(SF2020_avoiddeforest2.px, "rasters/PGM/2020_avoiddeforest2/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_avoiddeforest2.ls <- focal(SF2020_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_avoiddeforest2.ls)<-"SFls"
SF2020_avoiddeforest2.ls[is.nan(SF2020_avoiddeforest2.ls)] <- 0
SF2020_avoiddeforest2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_avoiddeforest2.ls[])
#saving
writeRaster(SF2020_avoiddeforest2.ls, "rasters/PGM/2020_avoiddeforest2/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddeforest2.px", "SF2020_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_avoiddeforest2 <- calc(pgm.sfage[["pgm.sfage.2010real"]], fun=function(x){ifelse(is.na(x), x, ifelse(x==0, x, x+10))})
SFAge2020_avoiddeforest2[] <- ifelse(SF2020_avoiddeforest2[]==0, 0, SFAge2020_avoiddeforest2[])
writeRaster(SFAge2020_avoiddeforest2, "rasters/PGM/input/SFAge2020_avoiddeforest2.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_avoiddeforest2.px <- focal(SFAge2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_avoiddeforest2.px)<-"SFAgepx"
SFAge2020_avoiddeforest2.px[is.nan(SFAge2020_avoiddeforest2.px)] <- 0
SFAge2020_avoiddeforest2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_avoiddeforest2.px[])
#saving
writeRaster(SFAge2020_avoiddeforest2.px, "rasters/PGM/2020_avoiddeforest2/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_avoiddeforest2.ls <- focal(SFAge2020_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_avoiddeforest2.ls)<-"SFAgels"
SFAge2020_avoiddeforest2.ls[is.nan(SFAge2020_avoiddeforest2.ls)] <- 0
SFAge2020_avoiddeforest2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_avoiddeforest2.ls[])
#saving
writeRaster(SFAge2020_avoiddeforest2.ls, "rasters/PGM/2020_avoiddeforest2/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_avoiddeforest2.px", "SFAge2020_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young_avoiddeforest2 <- SFAge2020_avoiddeforest2
SF2020young_avoiddeforest2[] <- ifelse(SF2020young_avoiddeforest2[]>2, 1, 0)
writeRaster(SF2020young_avoiddeforest2, "rasters/PGM/input/SF2020young_avoiddeforest2.tif", format="GTiff", overwrite=T)

TF2020_avoiddeforest2 <- sum(UPF2020_avoiddeforest2, DPF2020_avoiddeforest2)
TF2020_avoiddeforest2 <- sum(TF2020_avoiddeforest2, SF2020young_avoiddeforest2)
writeRaster(TF2020_avoiddeforest2, "rasters/PGM/input/TF2020_avoiddeforest2.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_avoiddeforest2.px <- focal(TF2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_avoiddeforest2.px)<-"TFpx"
TF2020_avoiddeforest2.px[is.nan(TF2020_avoiddeforest2.px)] <- 0
TF2020_avoiddeforest2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_avoiddeforest2.px[])
#saving
writeRaster(TF2020_avoiddeforest2.px, "rasters/PGM/2020_avoiddeforest2/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_avoiddeforest2.ls <- focal(TF2020_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_avoiddeforest2.ls)<-"TFls"
TF2020_avoiddeforest2.ls[is.nan(TF2020_avoiddeforest2.ls)] <- 0
TF2020_avoiddeforest2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_avoiddeforest2.ls[])
#saving
writeRaster(TF2020_avoiddeforest2.ls, "rasters/PGM/2020_avoiddeforest2/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddeforest2young", "TF2020_avoiddeforest2.px", "TF2020_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature_avoiddeforest2 <- SFAge2020_avoiddeforest2
SF2020mature_avoiddeforest2[] <- ifelse(SF2020mature_avoiddeforest2[]>5, 1, 0)
writeRaster(SF2020mature_avoiddeforest2, "rasters/PGM/input/SF2020mature_avoiddeforest2.tif", format="GTiff", overwrite=T)

MF2020_avoiddeforest2 <- sum(UPF2020_avoiddeforest2, DPF2020_avoiddeforest2)
MF2020_avoiddeforest2 <- sum(MF2020_avoiddeforest2, SF2020mature_avoiddeforest2)
writeRaster(MF2020_avoiddeforest2, "rasters/PGM/input/MF2020_avoiddeforest2.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_avoiddeforest2.px <- focal(MF2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_avoiddeforest2.px)<-"MFpx"
MF2020_avoiddeforest2.px[is.nan(MF2020_avoiddeforest2.px)] <- 0
MF2020_avoiddeforest2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_avoiddeforest2.px[])
#saving
writeRaster(MF2020_avoiddeforest2.px, "rasters/PGM/2020_avoiddeforest2/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_avoiddeforest2.ls <- focal(MF2020_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_avoiddeforest2.ls)<-"MFls"
MF2020_avoiddeforest2.ls[is.nan(MF2020_avoiddeforest2.ls)] <- 0
MF2020_avoiddeforest2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_avoiddeforest2.ls[])
#saving
writeRaster(MF2020_avoiddeforest2.ls, "rasters/PGM/2020_avoiddeforest2/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoiddeforest2mature", "MF2020_avoiddeforest2.px", "MF2020_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_avoiddeforest2 <- raster("rasters/PGM/input/MF2020_avoiddeforest2.tif")
inv.MF2020_avoiddeforest2 <- MF2020_avoiddeforest2
inv.MF2020_avoiddeforest2[inv.MF2020_avoiddeforest2==1]<-NA
#cheking
#inv.MF2020_avoiddeforest2
#plot(inv.MF2020_avoiddeforest2)

edge.dist.2020_avoiddeforest2 <- distance(inv.MF2020_avoiddeforest2, doEdge=T)
names(edge.dist.2020_avoiddeforest2)<-"edgedist"
edge.dist.2020_avoiddeforest2[is.nan(edge.dist.2020_avoiddeforest2)] <- 0
edge.dist.2020_avoiddeforest2[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge.dist.2020_avoiddeforest2[])
#saving
writeRaster(edge.dist.2020_avoiddeforest2, "rasters/PGM/2020_avoiddeforest2/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_avoiddeforest2")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_avoiddeforest2 <- edge.dist.2020_avoiddeforest2
edge2020_avoiddeforest2[] <- ifelse(edge2020_avoiddeforest2[] < 200, 0, ifelse(edge2020_avoiddeforest2[]>300, 0, 1))
writeRaster(edge2020_avoiddeforest2, "rasters/PGM/input/edge2020_avoiddeforest2.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_avoiddeforest2.px <- focal(edge2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_avoiddeforest2.px)<-"edgepx"
edge2020_avoiddeforest2.px[is.nan(edge2020_avoiddeforest2.px)] <- 0
edge2020_avoiddeforest2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_avoiddeforest2.px[])
#saving
writeRaster(edge2020_avoiddeforest2.px, "rasters/PGM/2020_avoiddeforest2/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_avoiddeforest2.ls <- focal(edge2020_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_avoiddeforest2.ls)<-"edgels"
edge2020_avoiddeforest2.ls[is.nan(edge2020_avoiddeforest2.ls)] <- 0
edge2020_avoiddeforest2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_avoiddeforest2.ls[])
#saving
writeRaster(edge2020_avoiddeforest2.ls, "rasters/PGM/2020_avoiddeforest2/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_avoiddeforest2.px", "edge2020_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Restoration without avoiding ======================================
### Undegraded primary forest
#UPF2020_restor_wo_avoid <- raster("rasters/PGM/input/LULC/UPF2020_restor_wo_avoid.tif")

### mean upf cover in local scale (90m)
UPF2020_restor_wo_avoid.px <- focal(UPF2020_restor_wo_avoid, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_restor_wo_avoid.px)<-"UPFpx"
UPF2020_restor_wo_avoid.px[is.nan(UPF2020_restor_wo_avoid.px)] <- 0
UPF2020_restor_wo_avoid.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_restor_wo_avoid.px[])
#saving
writeRaster(UPF2020_restor_wo_avoid.px, "rasters/PGM/2020_restor_wo_avoid/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_restor_wo_avoid.ls <- focal(UPF2020_restor_wo_avoid, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_restor_wo_avoid.ls)<-"UPFls"
UPF2020_restor_wo_avoid.ls[is.nan(UPF2020_restor_wo_avoid.ls)] <- 0
UPF2020_restor_wo_avoid.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_restor_wo_avoid.ls[])
#saving
writeRaster(UPF2020_restor_wo_avoid.ls, "rasters/PGM/2020_restor_wo_avoid/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_restor_wo_avoid.px", "UPF2020_restor_wo_avoid.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_restor_wo_avoid <- raster("rasters/PGM/input/LULC/uDPF2020_restor_wo_avoid.tif")
#RDPF2020_restor_wo_avoid <- raster("rasters/PGM/input/LULC/RDPF2020_restor_wo_avoid.tif")
DPF2020_restor_wo_avoid <- sum(uDPF2020_restor_wo_avoid, RDPF2020_restor_wo_avoid)

### mean dpf cover in local scale (90m)
DPF2020_restor_wo_avoid.px <- focal(DPF2020_restor_wo_avoid, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_restor_wo_avoid.px)<-"DPFpx"
DPF2020_restor_wo_avoid.px[is.nan(DPF2020_restor_wo_avoid.px)] <- 0
DPF2020_restor_wo_avoid.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_restor_wo_avoid.px[])
#saving
writeRaster(DPF2020_restor_wo_avoid.px, "rasters/PGM/2020_restor_wo_avoid/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_restor_wo_avoid.ls <- focal(DPF2020_restor_wo_avoid, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_restor_wo_avoid.ls)<-"DPFls"
DPF2020_restor_wo_avoid.ls[is.nan(DPF2020_restor_wo_avoid.ls)] <- 0
DPF2020_restor_wo_avoid.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_restor_wo_avoid.ls[])
#saving
writeRaster(DPF2020_restor_wo_avoid.ls, "rasters/PGM/2020_restor_wo_avoid/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_restor_wo_avoid.px", "DPF2020_restor_wo_avoid.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_restor_wo_avoid <- TSD2020
writeRaster(TSD2020_restor_wo_avoid, "rasters/PGM/input/TSD2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_restor_wo_avoid.px <- focal(TSD2020_restor_wo_avoid, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_restor_wo_avoid.px)<-"TSDpx"
TSD2020_restor_wo_avoid.px[is.nan(TSD2020_restor_wo_avoid.px)] <- 0
TSD2020_restor_wo_avoid.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_restor_wo_avoid.px[])
#saving
writeRaster(TSD2020_restor_wo_avoid.px, "rasters/PGM/2020_restor_wo_avoid/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_restor_wo_avoid.ls <- focal(TSD2020_restor_wo_avoid, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_restor_wo_avoid.ls)<-"TSDls"
TSD2020_restor_wo_avoid.ls[is.nan(TSD2020_restor_wo_avoid.ls)] <- 0
TSD2020_restor_wo_avoid.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_restor_wo_avoid.ls[])
#saving
writeRaster(TSD2020_restor_wo_avoid.ls, "rasters/PGM/2020_restor_wo_avoid/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_restor_wo_avoid.px", "TSD2020_restor_wo_avoid.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_restor_wo_avoid <- raster("rasters/PGM/input/LULC/uSF2020_restor_wo_avoid.tif")
#DSF2020_restor_wo_avoid <- raster("rasters/PGM/input/LULC/DSF2020_restor_wo_avoid.tif")
SF2020_restor_wo_avoid <- sum(uSF2020_restor_wo_avoid, DSF2020_restor_wo_avoid)

### mean sf cover in local scale (90m)
SF2020_restor_wo_avoid.px <- focal(SF2020_restor_wo_avoid, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_restor_wo_avoid.px)<-"SFpx"
SF2020_restor_wo_avoid.px[is.nan(SF2020_restor_wo_avoid.px)] <- 0
SF2020_restor_wo_avoid.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_restor_wo_avoid.px[])
#saving
writeRaster(SF2020_restor_wo_avoid.px, "rasters/PGM/2020_restor_wo_avoid/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_restor_wo_avoid.ls <- focal(SF2020_restor_wo_avoid, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_restor_wo_avoid.ls)<-"SFls"
SF2020_restor_wo_avoid.ls[is.nan(SF2020_restor_wo_avoid.ls)] <- 0
SF2020_restor_wo_avoid.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_restor_wo_avoid.ls[])
#saving
writeRaster(SF2020_restor_wo_avoid.ls, "rasters/PGM/2020_restor_wo_avoid/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_wo_avoid.px", "SF2020_restor_wo_avoid.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_restor_wo_avoid <- pgm.sfage[["pgm.sfage.2020real"]]
SFAge2020_restor_wo_avoid[] <- ifelse(SFAge2020_restor_wo_avoid[]==0 & SF2020_restor_wo_avoid[]==1, 10, SFAge2020_restor_wo_avoid[])
writeRaster(SFAge2020_restor_wo_avoid, "rasters/PGM/input/SFAge2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_restor_wo_avoid.px <- focal(SFAge2020_restor_wo_avoid, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_restor_wo_avoid.px)<-"SFAgepx"
SFAge2020_restor_wo_avoid.px[is.nan(SFAge2020_restor_wo_avoid.px)] <- 0
SFAge2020_restor_wo_avoid.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_restor_wo_avoid.px[])
#saving
writeRaster(SFAge2020_restor_wo_avoid.px, "rasters/PGM/2020_restor_wo_avoid/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_restor_wo_avoid.ls <- focal(SFAge2020_restor_wo_avoid, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_restor_wo_avoid.ls)<-"SFAgels"
SFAge2020_restor_wo_avoid.ls[is.nan(SFAge2020_restor_wo_avoid.ls)] <- 0
SFAge2020_restor_wo_avoid.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_restor_wo_avoid.ls[])
#saving
writeRaster(SFAge2020_restor_wo_avoid.ls, "rasters/PGM/2020_restor_wo_avoid/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_restor_wo_avoid.px", "SFAge2020_restor_wo_avoid.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020_restor_wo_avoidyoung <- SFAge2020_restor_wo_avoid
SF2020_restor_wo_avoidyoung[] <- ifelse(SF2020_restor_wo_avoidyoung[]>2, 1, 0)
writeRaster(SF2020_restor_wo_avoidyoung, "rasters/PGM/input/SF2020_restor_wo_avoidyoung.tif", format="GTiff", overwrite=T)

TF2020_restor_wo_avoid <- sum(UPF2020_restor_wo_avoid, DPF2020_restor_wo_avoid)
TF2020_restor_wo_avoid <- sum(TF2020_restor_wo_avoid, SF2020_restor_wo_avoidyoung)
writeRaster(TF2020_restor_wo_avoid, "rasters/PGM/input/TF2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_restor_wo_avoid.px <- focal(TF2020_restor_wo_avoid, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_restor_wo_avoid.px)<-"TFpx"
TF2020_restor_wo_avoid.px[is.nan(TF2020_restor_wo_avoid.px)] <- 0
TF2020_restor_wo_avoid.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_restor_wo_avoid.px[])
#saving
writeRaster(TF2020_restor_wo_avoid.px, "rasters/PGM/2020_restor_wo_avoid/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_restor_wo_avoid.ls <- focal(TF2020_restor_wo_avoid, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_restor_wo_avoid.ls)<-"TFls"
TF2020_restor_wo_avoid.ls[is.nan(TF2020_restor_wo_avoid.ls)] <- 0
TF2020_restor_wo_avoid.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_restor_wo_avoid.ls[])
#saving
writeRaster(TF2020_restor_wo_avoid.ls, "rasters/PGM/2020_restor_wo_avoid/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_wo_avoidyoung", "TF2020_restor_wo_avoid.px", "TF2020_restor_wo_avoid.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020_restor_wo_avoidmature <- SFAge2020_restor_wo_avoid
SF2020_restor_wo_avoidmature[] <- ifelse(SF2020_restor_wo_avoidmature[]>5, 1, 0)
writeRaster(SF2020_restor_wo_avoidmature, "rasters/PGM/input/SF2020_restor_wo_avoidmature.tif", format="GTiff", overwrite=T)

MF2020_restor_wo_avoid <- sum(UPF2020_restor_wo_avoid, DPF2020_restor_wo_avoid)
MF2020_restor_wo_avoid <- sum(MF2020_restor_wo_avoid, SF2020_restor_wo_avoidmature)
writeRaster(MF2020_restor_wo_avoid, "rasters/PGM/input/MF2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_restor_wo_avoid.px <- focal(MF2020_restor_wo_avoid, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_restor_wo_avoid.px)<-"MFpx"
MF2020_restor_wo_avoid.px[is.nan(MF2020_restor_wo_avoid.px)] <- 0
MF2020_restor_wo_avoid.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_restor_wo_avoid.px[])
#saving
writeRaster(MF2020_restor_wo_avoid.px, "rasters/PGM/2020_restor_wo_avoid/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_restor_wo_avoid.ls <- focal(MF2020_restor_wo_avoid, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_restor_wo_avoid.ls)<-"MFls"
MF2020_restor_wo_avoid.ls[is.nan(MF2020_restor_wo_avoid.ls)] <- 0
MF2020_restor_wo_avoid.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_restor_wo_avoid.ls[])
#saving
writeRaster(MF2020_restor_wo_avoid.ls, "rasters/PGM/2020_restor_wo_avoid/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_wo_avoidmature", "MF2020_restor_wo_avoid.px", "MF2020_restor_wo_avoid.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_restor_wo_avoid <- raster("rasters/PGM/input/MF2020_restor_wo_avoid.tif")
inv.MF2020_restor_wo_avoid <- MF2020_restor_wo_avoid
inv.MF2020_restor_wo_avoid[inv.MF2020_restor_wo_avoid==1]<-NA
#cheking
#inv.MF2020_restor_wo_avoid
#plot(inv.MF2020_restor_wo_avoid)

edge.dist.2020_restor_wo_avoid <- distance(inv.MF2020_restor_wo_avoid, doEdge=T)
names(edge.dist.2020_restor_wo_avoid)<-"edgedist"
edge.dist.2020_restor_wo_avoid[is.nan(edge.dist.2020_restor_wo_avoid)] <- 0
edge.dist.2020_restor_wo_avoid[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge.dist.2020_restor_wo_avoid[])
#saving
writeRaster(edge.dist.2020_restor_wo_avoid, "rasters/PGM/2020_restor_wo_avoid/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_restor_wo_avoid")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_restor_wo_avoid <- edge.dist.2020_restor_wo_avoid
edge2020_restor_wo_avoid[] <- ifelse(edge2020_restor_wo_avoid[] < 200, 0, ifelse(edge2020_restor_wo_avoid[]>300, 0, 1))
writeRaster(edge2020_restor_wo_avoid, "rasters/PGM/input/edge2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_restor_wo_avoid.px <- focal(edge2020_restor_wo_avoid, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_restor_wo_avoid.px)<-"edgepx"
edge2020_restor_wo_avoid.px[is.nan(edge2020_restor_wo_avoid.px)] <- 0
edge2020_restor_wo_avoid.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_restor_wo_avoid.px[])
#saving
writeRaster(edge2020_restor_wo_avoid.px, "rasters/PGM/2020_restor_wo_avoid/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_restor_wo_avoid.ls <- focal(edge2020_restor_wo_avoid, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_restor_wo_avoid.ls)<-"edgels"
edge2020_restor_wo_avoid.ls[is.nan(edge2020_restor_wo_avoid.ls)] <- 0
edge2020_restor_wo_avoid.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_restor_wo_avoid.ls[])
#saving
writeRaster(edge2020_restor_wo_avoid.ls, "rasters/PGM/2020_restor_wo_avoid/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_restor_wo_avoid.px", "edge2020_restor_wo_avoid.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Avoid both (all) ==================================================
### Undegraded primary forest
#UPF2020_avoidboth <- raster("rasters/PGM/input/LULC/UPF2020_avoidboth.tif")

### mean upf cover in local scale (90m)
UPF2020_avoidboth.px <- focal(UPF2020_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_avoidboth.px)<-"UPFpx"
UPF2020_avoidboth.px[is.nan(UPF2020_avoidboth.px)] <- 0
UPF2020_avoidboth.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_avoidboth.px[])
#saving
writeRaster(UPF2020_avoidboth.px, "rasters/PGM/2020_avoidboth/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_avoidboth.ls <- focal(UPF2020_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_avoidboth.ls)<-"UPFls"
UPF2020_avoidboth.ls[is.nan(UPF2020_avoidboth.ls)] <- 0
UPF2020_avoidboth.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_avoidboth.ls[])
#saving
writeRaster(UPF2020_avoidboth.ls, "rasters/PGM/2020_avoidboth/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_avoidboth.px", "UPF2020_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_avoidboth <- raster("rasters/PGM/input/LULC/uDPF2020_avoidboth.tif")
#RDPF2020_avoidboth <- raster("rasters/PGM/input/LULC/RDPF2020_avoidboth.tif")
DPF2020_avoidboth <- sum(uDPF2020_avoidboth, RDPF2020_avoidboth)

### mean dpf cover in local scale (90m)
DPF2020_avoidboth.px <- focal(DPF2020_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_avoidboth.px)<-"DPFpx"
DPF2020_avoidboth.px[is.nan(DPF2020_avoidboth.px)] <- 0
DPF2020_avoidboth.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_avoidboth.px[])
#saving
writeRaster(DPF2020_avoidboth.px, "rasters/PGM/2020_avoidboth/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_avoidboth.ls <- focal(DPF2020_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_avoidboth.ls)<-"DPFls"
DPF2020_avoidboth.ls[is.nan(DPF2020_avoidboth.ls)] <- 0
DPF2020_avoidboth.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_avoidboth.ls[])
#saving
writeRaster(DPF2020_avoidboth.ls, "rasters/PGM/2020_avoidboth/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_avoidboth.px", "DPF2020_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_avoidboth <- TSD2020_avoiddegrad
writeRaster(TSD2020_avoidboth, "rasters/PGM/input/TSD2020_avoidboth.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_avoidboth.px <- focal(TSD2020_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_avoidboth.px)<-"TSDpx"
TSD2020_avoidboth.px[is.nan(TSD2020_avoidboth.px)] <- 0
TSD2020_avoidboth.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_avoidboth.px[])
#saving
writeRaster(TSD2020_avoidboth.px, "rasters/PGM/2020_avoidboth/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_avoidboth.ls <- focal(TSD2020_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_avoidboth.ls)<-"TSDls"
TSD2020_avoidboth.ls[is.nan(TSD2020_avoidboth.ls)] <- 0
TSD2020_avoidboth.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_avoidboth.ls[])
#saving
writeRaster(TSD2020_avoidboth.ls, "rasters/PGM/2020_avoidboth/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_avoidboth.px", "TSD2020_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_avoidboth <- raster("rasters/PGM/input/LULC/uSF2020_avoidboth.tif")
#DSF2020_avoidboth <- raster("rasters/PGM/input/LULC/DSF2020_avoidboth.tif")
SF2020_avoidboth <- sum(uSF2020_avoidboth, DSF2020_avoidboth)

### mean sf cover in local scale (90m)
SF2020_avoidboth.px <- focal(SF2020_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_avoidboth.px)<-"SFpx"
SF2020_avoidboth.px[is.nan(SF2020_avoidboth.px)] <- 0
SF2020_avoidboth.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_avoidboth.px[])
#saving
writeRaster(SF2020_avoidboth.px, "rasters/PGM/2020_avoidboth/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_avoidboth.ls <- focal(SF2020_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_avoidboth.ls)<-"SFls"
SF2020_avoidboth.ls[is.nan(SF2020_avoidboth.ls)] <- 0
SF2020_avoidboth.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_avoidboth.ls[])
#saving
writeRaster(SF2020_avoidboth.ls, "rasters/PGM/2020_avoidboth/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoidboth.px", "SF2020_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_avoidboth <- SFAge2020_avoiddeforest
writeRaster(SFAge2020_avoidboth, "rasters/PGM/input/SFAge2020_avoidboth.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_avoidboth.px <- focal(SFAge2020_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_avoidboth.px)<-"SFAgepx"
SFAge2020_avoidboth.px[is.nan(SFAge2020_avoidboth.px)] <- 0
SFAge2020_avoidboth.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_avoidboth.px[])
#saving
writeRaster(SFAge2020_avoidboth.px, "rasters/PGM/2020_avoidboth/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_avoidboth.ls <- focal(SFAge2020_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_avoidboth.ls)<-"SFAgels"
SFAge2020_avoidboth.ls[is.nan(SFAge2020_avoidboth.ls)] <- 0
SFAge2020_avoidboth.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_avoidboth.ls[])
#saving
writeRaster(SFAge2020_avoidboth.ls, "rasters/PGM/2020_avoidboth/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_avoidboth.px", "SFAge2020_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young_avoidboth <- SFAge2020_avoidboth
SF2020young_avoidboth[] <- ifelse(SF2020young_avoidboth[]>2, 1, 0)
writeRaster(SF2020young_avoidboth, "rasters/PGM/input/SF2020young_avoidboth.tif", format="GTiff", overwrite=T)

TF2020_avoidboth <- sum(UPF2020_avoidboth, DPF2020_avoidboth)
TF2020_avoidboth <- sum(TF2020_avoidboth, SF2020young_avoidboth)
writeRaster(TF2020_avoidboth, "rasters/PGM/input/TF2020_avoidboth.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_avoidboth.px <- focal(TF2020_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_avoidboth.px)<-"TFpx"
TF2020_avoidboth.px[is.nan(TF2020_avoidboth.px)] <- 0
TF2020_avoidboth.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_avoidboth.px[])
#saving
writeRaster(TF2020_avoidboth.px, "rasters/PGM/2020_avoidboth/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_avoidboth.ls <- focal(TF2020_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_avoidboth.ls)<-"TFls"
TF2020_avoidboth.ls[is.nan(TF2020_avoidboth.ls)] <- 0
TF2020_avoidboth.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_avoidboth.ls[])
#saving
writeRaster(TF2020_avoidboth.ls, "rasters/PGM/2020_avoidboth/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoidbothyoung", "TF2020_avoidboth.px", "TF2020_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature_avoidboth <- SFAge2020_avoidboth
SF2020mature_avoidboth[] <- ifelse(SF2020mature_avoidboth[]>5, 1, 0)
writeRaster(SF2020mature_avoidboth, "rasters/PGM/input/SF2020mature_avoidboth.tif", format="GTiff", overwrite=T)

MF2020_avoidboth <- sum(UPF2020_avoidboth, DPF2020_avoidboth)
MF2020_avoidboth <- sum(MF2020_avoidboth, SF2020mature_avoidboth)
writeRaster(MF2020_avoidboth, "rasters/PGM/input/MF2020_avoidboth.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_avoidboth.px <- focal(MF2020_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_avoidboth.px)<-"MFpx"
MF2020_avoidboth.px[is.nan(MF2020_avoidboth.px)] <- 0
MF2020_avoidboth.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_avoidboth.px[])
#saving
writeRaster(MF2020_avoidboth.px, "rasters/PGM/2020_avoidboth/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_avoidboth.ls <- focal(MF2020_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_avoidboth.ls)<-"MFls"
MF2020_avoidboth.ls[is.nan(MF2020_avoidboth.ls)] <- 0
MF2020_avoidboth.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_avoidboth.ls[])
#saving
writeRaster(MF2020_avoidboth.ls, "rasters/PGM/2020_avoidboth/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoidbothmature", "MF2020_avoidboth.px", "MF2020_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_avoidboth <- raster("rasters/PGM/input/MF2020_avoidboth.tif")
inv.MF2020_avoidboth <- MF2020_avoidboth
inv.MF2020_avoidboth[inv.MF2020_avoidboth==1]<-NA
#cheking
#inv.MF2020_avoidboth
#plot(inv.MF2020_avoidboth)

edge.dist.2020_avoidboth <- distance(inv.MF2020_avoidboth, doEdge=T)
names(edge.dist.2020_avoidboth)<-"edgedist"
edge.dist.2020_avoidboth[is.nan(edge.dist.2020_avoidboth)] <- 0
edge.dist.2020_avoidboth[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge.dist.2020_avoidboth[])
#saving
writeRaster(edge.dist.2020_avoidboth, "rasters/PGM/2020_avoidboth/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_avoidboth")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_avoidboth <- edge.dist.2020_avoidboth
edge2020_avoidboth[] <- ifelse(edge2020_avoidboth[] < 200, 0, ifelse(edge2020_avoidboth[]>300, 0, 1))
writeRaster(edge2020_avoidboth, "rasters/PGM/input/edge2020_avoidboth.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_avoidboth.px <- focal(edge2020_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_avoidboth.px)<-"edgepx"
edge2020_avoidboth.px[is.nan(edge2020_avoidboth.px)] <- 0
edge2020_avoidboth.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_avoidboth.px[])
#saving
writeRaster(edge2020_avoidboth.px, "rasters/PGM/2020_avoidboth/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_avoidboth.ls <- focal(edge2020_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_avoidboth.ls)<-"edgels"
edge2020_avoidboth.ls[is.nan(edge2020_avoidboth.ls)] <- 0
edge2020_avoidboth.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_avoidboth.ls[])
#saving
writeRaster(edge2020_avoidboth.ls, "rasters/PGM/2020_avoidboth/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_avoidboth.px", "edge2020_avoidboth.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Avoid both (Primary forests only) =================================
### Undegraded primary forest
#UPF2020_avoidboth2 <- raster("rasters/PGM/input/LULC/UPF2020_avoidboth2.tif")

### mean upf cover in local scale (90m)
UPF2020_avoidboth2.px <- focal(UPF2020_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_avoidboth2.px)<-"UPFpx"
UPF2020_avoidboth2.px[is.nan(UPF2020_avoidboth2.px)] <- 0
UPF2020_avoidboth2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_avoidboth2.px[])
#saving
writeRaster(UPF2020_avoidboth2.px, "rasters/PGM/2020_avoidboth2/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_avoidboth2.ls <- focal(UPF2020_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_avoidboth2.ls)<-"UPFls"
UPF2020_avoidboth2.ls[is.nan(UPF2020_avoidboth2.ls)] <- 0
UPF2020_avoidboth2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_avoidboth2.ls[])
#saving
writeRaster(UPF2020_avoidboth2.ls, "rasters/PGM/2020_avoidboth2/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_avoidboth2.px", "UPF2020_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_avoidboth2 <- raster("rasters/PGM/input/LULC/uDPF2020_avoidboth2.tif")
#RDPF2020_avoidboth2 <- raster("rasters/PGM/input/LULC/RDPF2020_avoidboth2.tif")
DPF2020_avoidboth2 <- sum(uDPF2020_avoidboth2, RDPF2020_avoidboth2)

### mean dpf cover in local scale (90m)
DPF2020_avoidboth2.px <- focal(DPF2020_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_avoidboth2.px)<-"DPFpx"
DPF2020_avoidboth2.px[is.nan(DPF2020_avoidboth2.px)] <- 0
DPF2020_avoidboth2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_avoidboth2.px[])
#saving
writeRaster(DPF2020_avoidboth2.px, "rasters/PGM/2020_avoidboth2/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_avoidboth2.ls <- focal(DPF2020_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_avoidboth2.ls)<-"DPFls"
DPF2020_avoidboth2.ls[is.nan(DPF2020_avoidboth2.ls)] <- 0
DPF2020_avoidboth2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_avoidboth2.ls[])
#saving
writeRaster(DPF2020_avoidboth2.ls, "rasters/PGM/2020_avoidboth2/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_avoidboth2.px", "DPF2020_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_avoidboth2 <- TSD2020_avoiddegrad2
writeRaster(TSD2020_avoidboth2, "rasters/PGM/input/TSD2020_avoidboth2.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_avoidboth2.px <- focal(TSD2020_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_avoidboth2.px)<-"TSDpx"
TSD2020_avoidboth2.px[is.nan(TSD2020_avoidboth2.px)] <- 0
TSD2020_avoidboth2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_avoidboth2.px[])
#saving
writeRaster(TSD2020_avoidboth2.px, "rasters/PGM/2020_avoidboth2/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_avoidboth2.ls <- focal(TSD2020_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_avoidboth2.ls)<-"TSDls"
TSD2020_avoidboth2.ls[is.nan(TSD2020_avoidboth2.ls)] <- 0
TSD2020_avoidboth2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_avoidboth2.ls[])
#saving
writeRaster(TSD2020_avoidboth2.ls, "rasters/PGM/2020_avoidboth2/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_avoidboth2.px", "TSD2020_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_avoidboth2 <- raster("rasters/PGM/input/LULC/uSF2020_avoidboth2.tif")
#DSF2020_avoidboth2 <- raster("rasters/PGM/input/LULC/DSF2020_avoidboth2.tif")
SF2020_avoidboth2 <- sum(uSF2020_avoidboth2, DSF2020_avoidboth2)

### mean sf cover in local scale (90m)
SF2020_avoidboth2.px <- focal(SF2020_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_avoidboth2.px)<-"SFpx"
SF2020_avoidboth2.px[is.nan(SF2020_avoidboth2.px)] <- 0
SF2020_avoidboth2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_avoidboth2.px[])
#saving
writeRaster(SF2020_avoidboth2.px, "rasters/PGM/2020_avoidboth2/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_avoidboth2.ls <- focal(SF2020_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_avoidboth2.ls)<-"SFls"
SF2020_avoidboth2.ls[is.nan(SF2020_avoidboth2.ls)] <- 0
SF2020_avoidboth2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_avoidboth2.ls[])
#saving
writeRaster(SF2020_avoidboth2.ls, "rasters/PGM/2020_avoidboth2/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoidboth2.px", "SF2020_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_avoidboth2 <- SFAge2020_avoiddeforest2
writeRaster(SFAge2020_avoidboth2, "rasters/PGM/input/SFAge2020_avoidboth2.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_avoidboth2.px <- focal(SFAge2020_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_avoidboth2.px)<-"SFAgepx"
SFAge2020_avoidboth2.px[is.nan(SFAge2020_avoidboth2.px)] <- 0
SFAge2020_avoidboth2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_avoidboth2.px[])
#saving
writeRaster(SFAge2020_avoidboth2.px, "rasters/PGM/2020_avoidboth2/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_avoidboth2.ls <- focal(SFAge2020_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_avoidboth2.ls)<-"SFAgels"
SFAge2020_avoidboth2.ls[is.nan(SFAge2020_avoidboth2.ls)] <- 0
SFAge2020_avoidboth2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_avoidboth2.ls[])
#saving
writeRaster(SFAge2020_avoidboth2.ls, "rasters/PGM/2020_avoidboth2/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_avoidboth2.px", "SFAge2020_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young_avoidboth2 <- SFAge2020_avoidboth2
SF2020young_avoidboth2[] <- ifelse(SF2020young_avoidboth2[]>2, 1, 0)
writeRaster(SF2020young_avoidboth2, "rasters/PGM/input/SF2020young_avoidboth2.tif", format="GTiff", overwrite=T)

TF2020_avoidboth2 <- sum(UPF2020_avoidboth2, DPF2020_avoidboth2)
TF2020_avoidboth2 <- sum(TF2020_avoidboth2, SF2020young_avoidboth2)
writeRaster(TF2020_avoidboth2, "rasters/PGM/input/TF2020_avoidboth2.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_avoidboth2.px <- focal(TF2020_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_avoidboth2.px)<-"TFpx"
TF2020_avoidboth2.px[is.nan(TF2020_avoidboth2.px)] <- 0
TF2020_avoidboth2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_avoidboth2.px[])
#saving
writeRaster(TF2020_avoidboth2.px, "rasters/PGM/2020_avoidboth2/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_avoidboth2.ls <- focal(TF2020_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_avoidboth2.ls)<-"TFls"
TF2020_avoidboth2.ls[is.nan(TF2020_avoidboth2.ls)] <- 0
TF2020_avoidboth2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_avoidboth2.ls[])
#saving
writeRaster(TF2020_avoidboth2.ls, "rasters/PGM/2020_avoidboth2/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoidboth2young", "TF2020_avoidboth2.px", "TF2020_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature_avoidboth2 <- SFAge2020_avoidboth2
SF2020mature_avoidboth2[] <- ifelse(SF2020mature_avoidboth2[]>5, 1, 0)
writeRaster(SF2020mature_avoidboth2, "rasters/PGM/input/SF2020mature_avoidboth2.tif", format="GTiff", overwrite=T)

MF2020_avoidboth2 <- sum(UPF2020_avoidboth2, DPF2020_avoidboth2)
MF2020_avoidboth2 <- sum(MF2020_avoidboth2, SF2020mature_avoidboth2)
writeRaster(MF2020_avoidboth2, "rasters/PGM/input/MF2020_avoidboth2.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_avoidboth2.px <- focal(MF2020_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_avoidboth2.px)<-"MFpx"
MF2020_avoidboth2.px[is.nan(MF2020_avoidboth2.px)] <- 0
MF2020_avoidboth2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_avoidboth2.px[])
#saving
writeRaster(MF2020_avoidboth2.px, "rasters/PGM/2020_avoidboth2/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_avoidboth2.ls <- focal(MF2020_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_avoidboth2.ls)<-"MFls"
MF2020_avoidboth2.ls[is.nan(MF2020_avoidboth2.ls)] <- 0
MF2020_avoidboth2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_avoidboth2.ls[])
#saving
writeRaster(MF2020_avoidboth2.ls, "rasters/PGM/2020_avoidboth2/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_avoidboth2mature", "MF2020_avoidboth2.px", "MF2020_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_avoidboth2 <- raster("rasters/PGM/input/MF2020_avoidboth2.tif")
inv.MF2020_avoidboth2 <- MF2020_avoidboth2
inv.MF2020_avoidboth2[inv.MF2020_avoidboth2==1]<-NA
#cheking
#inv.MF2020_avoidboth2
#plot(inv.MF2020_avoidboth2)

edge.dist.2020_avoidboth2 <- distance(inv.MF2020_avoidboth2, doEdge=T)
names(edge.dist.2020_avoidboth2)<-"edgedist"
edge.dist.2020_avoidboth2[is.nan(edge.dist.2020_avoidboth2)] <- 0
edge.dist.2020_avoidboth2[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge.dist.2020_avoidboth2[])
#saving
writeRaster(edge.dist.2020_avoidboth2, "rasters/PGM/2020_avoidboth2/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_avoidboth2")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_avoidboth2 <- edge.dist.2020_avoidboth2
edge2020_avoidboth2[] <- ifelse(edge2020_avoidboth2[] < 200, 0, ifelse(edge2020_avoidboth2[]>300, 0, 1))
writeRaster(edge2020_avoidboth2, "rasters/PGM/input/edge2020_avoidboth2.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_avoidboth2.px <- focal(edge2020_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_avoidboth2.px)<-"edgepx"
edge2020_avoidboth2.px[is.nan(edge2020_avoidboth2.px)] <- 0
edge2020_avoidboth2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_avoidboth2.px[])
#saving
writeRaster(edge2020_avoidboth2.px, "rasters/PGM/2020_avoidboth2/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_avoidboth2.ls <- focal(edge2020_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_avoidboth2.ls)<-"edgels"
edge2020_avoidboth2.ls[is.nan(edge2020_avoidboth2.ls)] <- 0
edge2020_avoidboth2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_avoidboth2.ls[])
#saving
writeRaster(edge2020_avoidboth2.ls, "rasters/PGM/2020_avoidboth2/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_avoidboth2.px", "edge2020_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Restoration and avoid deforestation (all) =========================
### Undegraded primary forest
#UPF2020_restor_n_avoiddeforest <- raster("rasters/PGM/input/LULC/UPF2020_restor_n_avoiddeforest.tif")

### mean upf cover in local scale (90m)
UPF2020_restor_n_avoiddeforest.px <- focal(UPF2020_restor_n_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_restor_n_avoiddeforest.px)<-"UPFpx"
UPF2020_restor_n_avoiddeforest.px[is.nan(UPF2020_restor_n_avoiddeforest.px)] <- 0
UPF2020_restor_n_avoiddeforest.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_restor_n_avoiddeforest.px[])
#saving
writeRaster(UPF2020_restor_n_avoiddeforest.px, "rasters/PGM/2020_restor_n_avoiddeforest/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_restor_n_avoiddeforest.ls <- focal(UPF2020_restor_n_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_restor_n_avoiddeforest.ls)<-"UPFls"
UPF2020_restor_n_avoiddeforest.ls[is.nan(UPF2020_restor_n_avoiddeforest.ls)] <- 0
UPF2020_restor_n_avoiddeforest.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_restor_n_avoiddeforest.ls[])
#saving
writeRaster(UPF2020_restor_n_avoiddeforest.ls, "rasters/PGM/2020_restor_n_avoiddeforest/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_restor_n_avoiddeforest.px", "UPF2020_restor_n_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_restor_n_avoiddeforest <- raster("rasters/PGM/input/LULC/uDPF2020_restor_n_avoiddeforest.tif")
#RDPF2020_restor_n_avoiddeforest <- raster("rasters/PGM/input/LULC/RDPF2020_restor_n_avoiddeforest.tif")
DPF2020_restor_n_avoiddeforest <- sum(uDPF2020_restor_n_avoiddeforest, RDPF2020_restor_n_avoiddeforest)

### mean dpf cover in local scale (90m)
DPF2020_restor_n_avoiddeforest.px <- focal(DPF2020_restor_n_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_restor_n_avoiddeforest.px)<-"DPFpx"
DPF2020_restor_n_avoiddeforest.px[is.nan(DPF2020_restor_n_avoiddeforest.px)] <- 0
DPF2020_restor_n_avoiddeforest.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_restor_n_avoiddeforest.px[])
#saving
writeRaster(DPF2020_restor_n_avoiddeforest.px, "rasters/PGM/2020_restor_n_avoiddeforest/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_restor_n_avoiddeforest.ls <- focal(DPF2020_restor_n_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_restor_n_avoiddeforest.ls)<-"DPFls"
DPF2020_restor_n_avoiddeforest.ls[is.nan(DPF2020_restor_n_avoiddeforest.ls)] <- 0
DPF2020_restor_n_avoiddeforest.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_restor_n_avoiddeforest.ls[])
#saving
writeRaster(DPF2020_restor_n_avoiddeforest.ls, "rasters/PGM/2020_restor_n_avoiddeforest/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_restor_n_avoiddeforest.px", "DPF2020_restor_n_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_restor_n_avoiddeforest <- TSD2020_avoiddeforest
writeRaster(TSD2020_restor_n_avoiddeforest, "rasters/PGM/input/TSD2020_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_restor_n_avoiddeforest.px <- focal(TSD2020_restor_n_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_restor_n_avoiddeforest.px)<-"TSDpx"
TSD2020_restor_n_avoiddeforest.px[is.nan(TSD2020_restor_n_avoiddeforest.px)] <- 0
TSD2020_restor_n_avoiddeforest.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_restor_n_avoiddeforest.px[])
#saving
writeRaster(TSD2020_restor_n_avoiddeforest.px, "rasters/PGM/2020_restor_n_avoiddeforest/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_restor_n_avoiddeforest.ls <- focal(TSD2020_restor_n_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_restor_n_avoiddeforest.ls)<-"TSDls"
TSD2020_restor_n_avoiddeforest.ls[is.nan(TSD2020_restor_n_avoiddeforest.ls)] <- 0
TSD2020_restor_n_avoiddeforest.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_restor_n_avoiddeforest.ls[])
#saving
writeRaster(TSD2020_restor_n_avoiddeforest.ls, "rasters/PGM/2020_restor_n_avoiddeforest/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_restor_n_avoiddeforest.px", "TSD2020_restor_n_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_restor_n_avoiddeforest <- raster("rasters/PGM/input/LULC/uSF2020_restor_n_avoiddeforest.tif")
#DSF2020_restor_n_avoiddeforest <- raster("rasters/PGM/input/LULC/DSF2020_restor_n_avoiddeforest.tif")
SF2020_restor_n_avoiddeforest <- sum(uSF2020_restor_n_avoiddeforest, DSF2020_restor_n_avoiddeforest)

### mean sf cover in local scale (90m)
SF2020_restor_n_avoiddeforest.px <- focal(SF2020_restor_n_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_restor_n_avoiddeforest.px)<-"SFpx"
SF2020_restor_n_avoiddeforest.px[is.nan(SF2020_restor_n_avoiddeforest.px)] <- 0
SF2020_restor_n_avoiddeforest.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_restor_n_avoiddeforest.px[])
#saving
writeRaster(SF2020_restor_n_avoiddeforest.px, "rasters/PGM/2020_restor_n_avoiddeforest/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_restor_n_avoiddeforest.ls <- focal(SF2020_restor_n_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_restor_n_avoiddeforest.ls)<-"SFls"
SF2020_restor_n_avoiddeforest.ls[is.nan(SF2020_restor_n_avoiddeforest.ls)] <- 0
SF2020_restor_n_avoiddeforest.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_restor_n_avoiddeforest.ls[])
#saving
writeRaster(SF2020_restor_n_avoiddeforest.ls, "rasters/PGM/2020_restor_n_avoiddeforest/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoiddeforest.px", "SF2020_restor_n_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_restor_n_avoiddeforest <- SFAge2020_avoiddeforest
SFAge2020_restor_n_avoiddeforest[] <- ifelse(SFAge2020_restor_n_avoiddeforest[]==0 & SF2020_restor_n_avoiddeforest[]==1, 10, SFAge2020_restor_n_avoiddeforest[])
writeRaster(SFAge2020_restor_n_avoiddeforest, "rasters/PGM/input/SFAge2020_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_restor_n_avoiddeforest.px <- focal(SFAge2020_restor_n_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_restor_n_avoiddeforest.px)<-"SFAgepx"
SFAge2020_restor_n_avoiddeforest.px[is.nan(SFAge2020_restor_n_avoiddeforest.px)] <- 0
SFAge2020_restor_n_avoiddeforest.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_restor_n_avoiddeforest.px[])
#saving
writeRaster(SFAge2020_restor_n_avoiddeforest.px, "rasters/PGM/2020_restor_n_avoiddeforest/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_restor_n_avoiddeforest.ls <- focal(SFAge2020_restor_n_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_restor_n_avoiddeforest.ls)<-"SFAgels"
SFAge2020_restor_n_avoiddeforest.ls[is.nan(SFAge2020_restor_n_avoiddeforest.ls)] <- 0
SFAge2020_restor_n_avoiddeforest.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_restor_n_avoiddeforest.ls[])
#saving
writeRaster(SFAge2020_restor_n_avoiddeforest.ls, "rasters/PGM/2020_restor_n_avoiddeforest/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_restor_n_avoiddeforest.px", "SFAge2020_restor_n_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young_restor_n_avoiddeforest <- SFAge2020_restor_n_avoiddeforest
SF2020young_restor_n_avoiddeforest[] <- ifelse(SF2020young_restor_n_avoiddeforest[]>2, 1, 0)
writeRaster(SF2020young_restor_n_avoiddeforest, "rasters/PGM/input/SF2020young_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)

TF2020_restor_n_avoiddeforest <- sum(UPF2020_restor_n_avoiddeforest, DPF2020_restor_n_avoiddeforest)
TF2020_restor_n_avoiddeforest <- sum(TF2020_restor_n_avoiddeforest, SF2020young_restor_n_avoiddeforest)
writeRaster(TF2020_restor_n_avoiddeforest, "rasters/PGM/input/TF2020_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_restor_n_avoiddeforest.px <- focal(TF2020_restor_n_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_restor_n_avoiddeforest.px)<-"TFpx"
TF2020_restor_n_avoiddeforest.px[is.nan(TF2020_restor_n_avoiddeforest.px)] <- 0
TF2020_restor_n_avoiddeforest.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_restor_n_avoiddeforest.px[])
#saving
writeRaster(TF2020_restor_n_avoiddeforest.px, "rasters/PGM/2020_restor_n_avoiddeforest/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_restor_n_avoiddeforest.ls <- focal(TF2020_restor_n_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_restor_n_avoiddeforest.ls)<-"TFls"
TF2020_restor_n_avoiddeforest.ls[is.nan(TF2020_restor_n_avoiddeforest.ls)] <- 0
TF2020_restor_n_avoiddeforest.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_restor_n_avoiddeforest.ls[])
#saving
writeRaster(TF2020_restor_n_avoiddeforest.ls, "rasters/PGM/2020_restor_n_avoiddeforest/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoiddeforestyoung", "TF2020_restor_n_avoiddeforest.px", "TF2020_restor_n_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature_restor_n_avoiddeforest <- SFAge2020_restor_n_avoiddeforest
SF2020mature_restor_n_avoiddeforest[] <- ifelse(SF2020mature_restor_n_avoiddeforest[]>5, 1, 0)
writeRaster(SF2020mature_restor_n_avoiddeforest, "rasters/PGM/input/SF2020mature_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)

MF2020_restor_n_avoiddeforest <- sum(UPF2020_restor_n_avoiddeforest, DPF2020_restor_n_avoiddeforest)
MF2020_restor_n_avoiddeforest <- sum(MF2020_restor_n_avoiddeforest, SF2020mature_restor_n_avoiddeforest)
writeRaster(MF2020_restor_n_avoiddeforest, "rasters/PGM/input/MF2020_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_restor_n_avoiddeforest.px <- focal(MF2020_restor_n_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_restor_n_avoiddeforest.px)<-"MFpx"
MF2020_restor_n_avoiddeforest.px[is.nan(MF2020_restor_n_avoiddeforest.px)] <- 0
MF2020_restor_n_avoiddeforest.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_restor_n_avoiddeforest.px[])
#saving
writeRaster(MF2020_restor_n_avoiddeforest.px, "rasters/PGM/2020_restor_n_avoiddeforest/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_restor_n_avoiddeforest.ls <- focal(MF2020_restor_n_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_restor_n_avoiddeforest.ls)<-"MFls"
MF2020_restor_n_avoiddeforest.ls[is.nan(MF2020_restor_n_avoiddeforest.ls)] <- 0
MF2020_restor_n_avoiddeforest.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_restor_n_avoiddeforest.ls[])
#saving
writeRaster(MF2020_restor_n_avoiddeforest.ls, "rasters/PGM/2020_restor_n_avoiddeforest/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoiddeforestmature", "MF2020_restor_n_avoiddeforest.px", "MF2020_restor_n_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_restor_n_avoiddeforest <- raster("rasters/PGM/input/MF2020_restor_n_avoiddeforest.tif")
inv.MF2020_restor_n_avoiddeforest <- MF2020_restor_n_avoiddeforest
inv.MF2020_restor_n_avoiddeforest[inv.MF2020_restor_n_avoiddeforest==1]<-NA
#cheking
#inv.MF2020_restor_n_avoiddeforest
#plot(inv.MF2020_restor_n_avoiddeforest)

edge.dist.2020_restor_n_avoiddeforest <- distance(inv.MF2020_restor_n_avoiddeforest, doEdge=T)
names(edge.dist.2020_restor_n_avoiddeforest)<-"edgedist"
edge.dist.2020_restor_n_avoiddeforest[is.nan(edge.dist.2020_restor_n_avoiddeforest)] <- 0
edge.dist.2020_restor_n_avoiddeforest[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge.dist.2020_restor_n_avoiddeforest[])
#saving
writeRaster(edge.dist.2020_restor_n_avoiddeforest, "rasters/PGM/2020_restor_n_avoiddeforest/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_restor_n_avoiddeforest")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_restor_n_avoiddeforest <- edge.dist.2020_restor_n_avoiddeforest
edge2020_restor_n_avoiddeforest[] <- ifelse(edge2020_restor_n_avoiddeforest[] < 200, 0, ifelse(edge2020_restor_n_avoiddeforest[]>300, 0, 1))
writeRaster(edge2020_restor_n_avoiddeforest, "rasters/PGM/input/edge2020_restor_n_avoiddeforest.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_restor_n_avoiddeforest.px <- focal(edge2020_restor_n_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_restor_n_avoiddeforest.px)<-"edgepx"
edge2020_restor_n_avoiddeforest.px[is.nan(edge2020_restor_n_avoiddeforest.px)] <- 0
edge2020_restor_n_avoiddeforest.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_restor_n_avoiddeforest.px[])
#saving
writeRaster(edge2020_restor_n_avoiddeforest.px, "rasters/PGM/2020_restor_n_avoiddeforest/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_restor_n_avoiddeforest.ls <- focal(edge2020_restor_n_avoiddeforest, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_restor_n_avoiddeforest.ls)<-"edgels"
edge2020_restor_n_avoiddeforest.ls[is.nan(edge2020_restor_n_avoiddeforest.ls)] <- 0
edge2020_restor_n_avoiddeforest.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_restor_n_avoiddeforest.ls[])
#saving
writeRaster(edge2020_restor_n_avoiddeforest.ls, "rasters/PGM/2020_restor_n_avoiddeforest/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_restor_n_avoiddeforest.px", "edge2020_restor_n_avoiddeforest.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Restoration and avoid deforestation (Primary forests only) ========
### Undegraded primary forest
#UPF2020_restor_n_avoiddeforest2 <- raster("rasters/PGM/input/LULC/UPF2020_restor_n_avoiddeforest2.tif")

### mean upf cover in local scale (90m)
UPF2020_restor_n_avoiddeforest2.px <- focal(UPF2020_restor_n_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_restor_n_avoiddeforest2.px)<-"UPFpx"
UPF2020_restor_n_avoiddeforest2.px[is.nan(UPF2020_restor_n_avoiddeforest2.px)] <- 0
UPF2020_restor_n_avoiddeforest2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_restor_n_avoiddeforest2.px[])
#saving
writeRaster(UPF2020_restor_n_avoiddeforest2.px, "rasters/PGM/2020_restor_n_avoiddeforest2/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_restor_n_avoiddeforest2.ls <- focal(UPF2020_restor_n_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_restor_n_avoiddeforest2.ls)<-"UPFls"
UPF2020_restor_n_avoiddeforest2.ls[is.nan(UPF2020_restor_n_avoiddeforest2.ls)] <- 0
UPF2020_restor_n_avoiddeforest2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_restor_n_avoiddeforest2.ls[])
#saving
writeRaster(UPF2020_restor_n_avoiddeforest2.ls, "rasters/PGM/2020_restor_n_avoiddeforest2/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_restor_n_avoiddeforest2.px", "UPF2020_restor_n_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_restor_n_avoiddeforest2 <- raster("rasters/PGM/input/LULC/uDPF2020_restor_n_avoiddeforest2.tif")
#RDPF2020_restor_n_avoiddeforest2 <- raster("rasters/PGM/input/LULC/RDPF2020_restor_n_avoiddeforest2.tif")
DPF2020_restor_n_avoiddeforest2 <- sum(uDPF2020_restor_n_avoiddeforest2, RDPF2020_restor_n_avoiddeforest2)

### mean dpf cover in local scale (90m)
DPF2020_restor_n_avoiddeforest2.px <- focal(DPF2020_restor_n_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_restor_n_avoiddeforest2.px)<-"DPFpx"
DPF2020_restor_n_avoiddeforest2.px[is.nan(DPF2020_restor_n_avoiddeforest2.px)] <- 0
DPF2020_restor_n_avoiddeforest2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_restor_n_avoiddeforest2.px[])
#saving
writeRaster(DPF2020_restor_n_avoiddeforest2.px, "rasters/PGM/2020_restor_n_avoiddeforest2/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_restor_n_avoiddeforest2.ls <- focal(DPF2020_restor_n_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_restor_n_avoiddeforest2.ls)<-"DPFls"
DPF2020_restor_n_avoiddeforest2.ls[is.nan(DPF2020_restor_n_avoiddeforest2.ls)] <- 0
DPF2020_restor_n_avoiddeforest2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_restor_n_avoiddeforest2.ls[])
#saving
writeRaster(DPF2020_restor_n_avoiddeforest2.ls, "rasters/PGM/2020_restor_n_avoiddeforest2/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_restor_n_avoiddeforest2.px", "DPF2020_restor_n_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_restor_n_avoiddeforest2 <- TSD2020_avoiddeforest2
writeRaster(TSD2020_restor_n_avoiddeforest2, "rasters/PGM/input/TSD2020_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_restor_n_avoiddeforest2.px <- focal(TSD2020_restor_n_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_restor_n_avoiddeforest2.px)<-"TSDpx"
TSD2020_restor_n_avoiddeforest2.px[is.nan(TSD2020_restor_n_avoiddeforest2.px)] <- 0
TSD2020_restor_n_avoiddeforest2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_restor_n_avoiddeforest2.px[])
#saving
writeRaster(TSD2020_restor_n_avoiddeforest2.px, "rasters/PGM/2020_restor_n_avoiddeforest2/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_restor_n_avoiddeforest2.ls <- focal(TSD2020_restor_n_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_restor_n_avoiddeforest2.ls)<-"TSDls"
TSD2020_restor_n_avoiddeforest2.ls[is.nan(TSD2020_restor_n_avoiddeforest2.ls)] <- 0
TSD2020_restor_n_avoiddeforest2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_restor_n_avoiddeforest2.ls[])
#saving
writeRaster(TSD2020_restor_n_avoiddeforest2.ls, "rasters/PGM/2020_restor_n_avoiddeforest2/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_restor_n_avoiddeforest2.px", "TSD2020_restor_n_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_restor_n_avoiddeforest2 <- raster("rasters/PGM/input/LULC/uSF2020_restor_n_avoiddeforest2.tif")
#DSF2020_restor_n_avoiddeforest2 <- raster("rasters/PGM/input/LULC/DSF2020_restor_n_avoiddeforest2.tif")
SF2020_restor_n_avoiddeforest2 <- sum(uSF2020_restor_n_avoiddeforest2, DSF2020_restor_n_avoiddeforest2)

### mean sf cover in local scale (90m)
SF2020_restor_n_avoiddeforest2.px <- focal(SF2020_restor_n_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_restor_n_avoiddeforest2.px)<-"SFpx"
SF2020_restor_n_avoiddeforest2.px[is.nan(SF2020_restor_n_avoiddeforest2.px)] <- 0
SF2020_restor_n_avoiddeforest2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_restor_n_avoiddeforest2.px[])
#saving
writeRaster(SF2020_restor_n_avoiddeforest2.px, "rasters/PGM/2020_restor_n_avoiddeforest2/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_restor_n_avoiddeforest2.ls <- focal(SF2020_restor_n_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_restor_n_avoiddeforest2.ls)<-"SFls"
SF2020_restor_n_avoiddeforest2.ls[is.nan(SF2020_restor_n_avoiddeforest2.ls)] <- 0
SF2020_restor_n_avoiddeforest2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_restor_n_avoiddeforest2.ls[])
#saving
writeRaster(SF2020_restor_n_avoiddeforest2.ls, "rasters/PGM/2020_restor_n_avoiddeforest2/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoiddeforest2.px", "SF2020_restor_n_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_restor_n_avoiddeforest2 <- SFAge2020_avoiddeforest2
SFAge2020_restor_n_avoiddeforest2[] <- ifelse(SFAge2020_restor_n_avoiddeforest2[]==0 & SF2020_restor_n_avoiddeforest2[]==1, 10, SFAge2020_restor_n_avoiddeforest2[])
writeRaster(SFAge2020_restor_n_avoiddeforest2, "rasters/PGM/input/SFAge2020_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_restor_n_avoiddeforest2.px <- focal(SFAge2020_restor_n_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_restor_n_avoiddeforest2.px)<-"SFAgepx"
SFAge2020_restor_n_avoiddeforest2.px[is.nan(SFAge2020_restor_n_avoiddeforest2.px)] <- 0
SFAge2020_restor_n_avoiddeforest2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_restor_n_avoiddeforest2.px[])
#saving
writeRaster(SFAge2020_restor_n_avoiddeforest2.px, "rasters/PGM/2020_restor_n_avoiddeforest2/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_restor_n_avoiddeforest2.ls <- focal(SFAge2020_restor_n_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_restor_n_avoiddeforest2.ls)<-"SFAgels"
SFAge2020_restor_n_avoiddeforest2.ls[is.nan(SFAge2020_restor_n_avoiddeforest2.ls)] <- 0
SFAge2020_restor_n_avoiddeforest2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_restor_n_avoiddeforest2.ls[])
#saving
writeRaster(SFAge2020_restor_n_avoiddeforest2.ls, "rasters/PGM/2020_restor_n_avoiddeforest2/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_restor_n_avoiddeforest2.px", "SFAge2020_restor_n_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young_restor_n_avoiddeforest2 <- SFAge2020_restor_n_avoiddeforest2
SF2020young_restor_n_avoiddeforest2[] <- ifelse(SF2020young_restor_n_avoiddeforest2[]>2, 1, 0)
writeRaster(SF2020young_restor_n_avoiddeforest2, "rasters/PGM/input/SF2020young_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)

TF2020_restor_n_avoiddeforest2 <- sum(UPF2020_restor_n_avoiddeforest2, DPF2020_restor_n_avoiddeforest2)
TF2020_restor_n_avoiddeforest2 <- sum(TF2020_restor_n_avoiddeforest2, SF2020young_restor_n_avoiddeforest2)
writeRaster(TF2020_restor_n_avoiddeforest2, "rasters/PGM/input/TF2020_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_restor_n_avoiddeforest2.px <- focal(TF2020_restor_n_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_restor_n_avoiddeforest2.px)<-"TFpx"
TF2020_restor_n_avoiddeforest2.px[is.nan(TF2020_restor_n_avoiddeforest2.px)] <- 0
TF2020_restor_n_avoiddeforest2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_restor_n_avoiddeforest2.px[])
#saving
writeRaster(TF2020_restor_n_avoiddeforest2.px, "rasters/PGM/2020_restor_n_avoiddeforest2/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_restor_n_avoiddeforest2.ls <- focal(TF2020_restor_n_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_restor_n_avoiddeforest2.ls)<-"TFls"
TF2020_restor_n_avoiddeforest2.ls[is.nan(TF2020_restor_n_avoiddeforest2.ls)] <- 0
TF2020_restor_n_avoiddeforest2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_restor_n_avoiddeforest2.ls[])
#saving
writeRaster(TF2020_restor_n_avoiddeforest2.ls, "rasters/PGM/2020_restor_n_avoiddeforest2/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoiddeforest2young", "TF2020_restor_n_avoiddeforest2.px", "TF2020_restor_n_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature_restor_n_avoiddeforest2 <- SFAge2020_restor_n_avoiddeforest2
SF2020mature_restor_n_avoiddeforest2[] <- ifelse(SF2020mature_restor_n_avoiddeforest2[]>5, 1, 0)
writeRaster(SF2020mature_restor_n_avoiddeforest2, "rasters/PGM/input/SF2020mature_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)

MF2020_restor_n_avoiddeforest2 <- sum(UPF2020_restor_n_avoiddeforest2, DPF2020_restor_n_avoiddeforest2)
MF2020_restor_n_avoiddeforest2 <- sum(MF2020_restor_n_avoiddeforest2, SF2020mature_restor_n_avoiddeforest2)
writeRaster(MF2020_restor_n_avoiddeforest2, "rasters/PGM/input/MF2020_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_restor_n_avoiddeforest2.px <- focal(MF2020_restor_n_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_restor_n_avoiddeforest2.px)<-"MFpx"
MF2020_restor_n_avoiddeforest2.px[is.nan(MF2020_restor_n_avoiddeforest2.px)] <- 0
MF2020_restor_n_avoiddeforest2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_restor_n_avoiddeforest2.px[])
#saving
writeRaster(MF2020_restor_n_avoiddeforest2.px, "rasters/PGM/2020_restor_n_avoiddeforest2/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_restor_n_avoiddeforest2.ls <- focal(MF2020_restor_n_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_restor_n_avoiddeforest2.ls)<-"MFls"
MF2020_restor_n_avoiddeforest2.ls[is.nan(MF2020_restor_n_avoiddeforest2.ls)] <- 0
MF2020_restor_n_avoiddeforest2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_restor_n_avoiddeforest2.ls[])
#saving
writeRaster(MF2020_restor_n_avoiddeforest2.ls, "rasters/PGM/2020_restor_n_avoiddeforest2/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoiddeforest2mature", "MF2020_restor_n_avoiddeforest2.px", "MF2020_restor_n_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_restor_n_avoiddeforest2 <- raster("rasters/PGM/input/MF2020_restor_n_avoiddeforest2.tif")
inv.MF2020_restor_n_avoiddeforest2 <- MF2020_restor_n_avoiddeforest2
inv.MF2020_restor_n_avoiddeforest2[inv.MF2020_restor_n_avoiddeforest2==1]<-NA
#cheking
#inv.MF2020_restor_n_avoiddeforest2
#plot(inv.MF2020_restor_n_avoiddeforest2)

edge.dist.2020_restor_n_avoiddeforest2 <- distance(inv.MF2020_restor_n_avoiddeforest2, doEdge=T)
names(edge.dist.2020_restor_n_avoiddeforest2)<-"edgedist"
edge.dist.2020_restor_n_avoiddeforest2[is.nan(edge.dist.2020_restor_n_avoiddeforest2)] <- 0
edge.dist.2020_restor_n_avoiddeforest2[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge.dist.2020_restor_n_avoiddeforest2[])
#saving
writeRaster(edge.dist.2020_restor_n_avoiddeforest2, "rasters/PGM/2020_restor_n_avoiddeforest2/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_restor_n_avoiddeforest2")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_restor_n_avoiddeforest2 <- edge.dist.2020_restor_n_avoiddeforest2
edge2020_restor_n_avoiddeforest2[] <- ifelse(edge2020_restor_n_avoiddeforest2[] < 200, 0, ifelse(edge2020_restor_n_avoiddeforest2[]>300, 0, 1))
writeRaster(edge2020_restor_n_avoiddeforest2, "rasters/PGM/input/edge2020_restor_n_avoiddeforest2.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_restor_n_avoiddeforest2.px <- focal(edge2020_restor_n_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_restor_n_avoiddeforest2.px)<-"edgepx"
edge2020_restor_n_avoiddeforest2.px[is.nan(edge2020_restor_n_avoiddeforest2.px)] <- 0
edge2020_restor_n_avoiddeforest2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_restor_n_avoiddeforest2.px[])
#saving
writeRaster(edge2020_restor_n_avoiddeforest2.px, "rasters/PGM/2020_restor_n_avoiddeforest2/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_restor_n_avoiddeforest2.ls <- focal(edge2020_restor_n_avoiddeforest2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_restor_n_avoiddeforest2.ls)<-"edgels"
edge2020_restor_n_avoiddeforest2.ls[is.nan(edge2020_restor_n_avoiddeforest2.ls)] <- 0
edge2020_restor_n_avoiddeforest2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_restor_n_avoiddeforest2.ls[])
#saving
writeRaster(edge2020_restor_n_avoiddeforest2.ls, "rasters/PGM/2020_restor_n_avoiddeforest2/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_restor_n_avoiddeforest2.px", "edge2020_restor_n_avoiddeforest2.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Restoration and avoid both (all) ==================================
### Undegraded primary forest
#UPF2020_restor_n_avoidboth <- raster("rasters/PGM/input/LULC/UPF2020_restor_n_avoidboth.tif")

### mean upf cover in local scale (90m)
UPF2020_restor_n_avoidboth.px <- focal(UPF2020_restor_n_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_restor_n_avoidboth.px)<-"UPFpx"
UPF2020_restor_n_avoidboth.px[is.nan(UPF2020_restor_n_avoidboth.px)] <- 0
UPF2020_restor_n_avoidboth.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_restor_n_avoidboth.px[])
#saving
writeRaster(UPF2020_restor_n_avoidboth.px, "rasters/PGM/2020_restor_n_avoidboth/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_restor_n_avoidboth.ls <- focal(UPF2020_restor_n_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_restor_n_avoidboth.ls)<-"UPFls"
UPF2020_restor_n_avoidboth.ls[is.nan(UPF2020_restor_n_avoidboth.ls)] <- 0
UPF2020_restor_n_avoidboth.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_restor_n_avoidboth.ls[])
#saving
writeRaster(UPF2020_restor_n_avoidboth.ls, "rasters/PGM/2020_restor_n_avoidboth/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_restor_n_avoidboth.px", "UPF2020_restor_n_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_restor_n_avoidboth <- raster("rasters/PGM/input/LULC/uDPF2020_restor_n_avoidboth.tif")
#RDPF2020_restor_n_avoidboth <- raster("rasters/PGM/input/LULC/RDPF2020_restor_n_avoidboth.tif")
DPF2020_restor_n_avoidboth <- sum(uDPF2020_restor_n_avoidboth, RDPF2020_restor_n_avoidboth)

### mean dpf cover in local scale (90m)
DPF2020_restor_n_avoidboth.px <- focal(DPF2020_restor_n_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_restor_n_avoidboth.px)<-"DPFpx"
DPF2020_restor_n_avoidboth.px[is.nan(DPF2020_restor_n_avoidboth.px)] <- 0
DPF2020_restor_n_avoidboth.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_restor_n_avoidboth.px[])
#saving
writeRaster(DPF2020_restor_n_avoidboth.px, "rasters/PGM/2020_restor_n_avoidboth/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_restor_n_avoidboth.ls <- focal(DPF2020_restor_n_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_restor_n_avoidboth.ls)<-"DPFls"
DPF2020_restor_n_avoidboth.ls[is.nan(DPF2020_restor_n_avoidboth.ls)] <- 0
DPF2020_restor_n_avoidboth.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_restor_n_avoidboth.ls[])
#saving
writeRaster(DPF2020_restor_n_avoidboth.ls, "rasters/PGM/2020_restor_n_avoidboth/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_restor_n_avoidboth.px", "DPF2020_restor_n_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_restor_n_avoidboth <- TSD2020_avoidboth
writeRaster(TSD2020_restor_n_avoidboth, "rasters/PGM/input/TSD2020_restor_n_avoidboth.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_restor_n_avoidboth.px <- focal(TSD2020_restor_n_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_restor_n_avoidboth.px)<-"TSDpx"
TSD2020_restor_n_avoidboth.px[is.nan(TSD2020_restor_n_avoidboth.px)] <- 0
TSD2020_restor_n_avoidboth.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_restor_n_avoidboth.px[])
#saving
writeRaster(TSD2020_restor_n_avoidboth.px, "rasters/PGM/2020_restor_n_avoidboth/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_restor_n_avoidboth.ls <- focal(TSD2020_restor_n_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_restor_n_avoidboth.ls)<-"TSDls"
TSD2020_restor_n_avoidboth.ls[is.nan(TSD2020_restor_n_avoidboth.ls)] <- 0
TSD2020_restor_n_avoidboth.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_restor_n_avoidboth.ls[])
#saving
writeRaster(TSD2020_restor_n_avoidboth.ls, "rasters/PGM/2020_restor_n_avoidboth/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_restor_n_avoidboth.px", "TSD2020_restor_n_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_restor_n_avoidboth <- raster("rasters/PGM/input/LULC/uSF2020_restor_n_avoidboth.tif")
#DSF2020_restor_n_avoidboth <- raster("rasters/PGM/input/LULC/DSF2020_restor_n_avoidboth.tif")
SF2020_restor_n_avoidboth <- sum(uSF2020_restor_n_avoidboth, DSF2020_restor_n_avoidboth)

### mean sf cover in local scale (90m)
SF2020_restor_n_avoidboth.px <- focal(SF2020_restor_n_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_restor_n_avoidboth.px)<-"SFpx"
SF2020_restor_n_avoidboth.px[is.nan(SF2020_restor_n_avoidboth.px)] <- 0
SF2020_restor_n_avoidboth.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_restor_n_avoidboth.px[])
#saving
writeRaster(SF2020_restor_n_avoidboth.px, "rasters/PGM/2020_restor_n_avoidboth/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_restor_n_avoidboth.ls <- focal(SF2020_restor_n_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_restor_n_avoidboth.ls)<-"SFls"
SF2020_restor_n_avoidboth.ls[is.nan(SF2020_restor_n_avoidboth.ls)] <- 0
SF2020_restor_n_avoidboth.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_restor_n_avoidboth.ls[])
#saving
writeRaster(SF2020_restor_n_avoidboth.ls, "rasters/PGM/2020_restor_n_avoidboth/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoidboth.px", "SF2020_restor_n_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_restor_n_avoidboth <- SFAge2020_avoidboth
SFAge2020_restor_n_avoidboth[] <- ifelse(SFAge2020_restor_n_avoidboth[]==0 & SF2020_restor_n_avoidboth[]==1, 10, SFAge2020_restor_n_avoidboth[])
writeRaster(SFAge2020_restor_n_avoidboth, "rasters/PGM/input/SFAge2020_restor_n_avoidboth.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_restor_n_avoidboth.px <- focal(SFAge2020_restor_n_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_restor_n_avoidboth.px)<-"SFAgepx"
SFAge2020_restor_n_avoidboth.px[is.nan(SFAge2020_restor_n_avoidboth.px)] <- 0
SFAge2020_restor_n_avoidboth.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_restor_n_avoidboth.px[])
#saving
writeRaster(SFAge2020_restor_n_avoidboth.px, "rasters/PGM/2020_restor_n_avoidboth/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_restor_n_avoidboth.ls <- focal(SFAge2020_restor_n_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_restor_n_avoidboth.ls)<-"SFAgels"
SFAge2020_restor_n_avoidboth.ls[is.nan(SFAge2020_restor_n_avoidboth.ls)] <- 0
SFAge2020_restor_n_avoidboth.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_restor_n_avoidboth.ls[])
#saving
writeRaster(SFAge2020_restor_n_avoidboth.ls, "rasters/PGM/2020_restor_n_avoidboth/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_restor_n_avoidboth.px", "SFAge2020_restor_n_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young_restor_n_avoidboth <- SFAge2020_restor_n_avoidboth
SF2020young_restor_n_avoidboth[] <- ifelse(SF2020young_restor_n_avoidboth[]>2, 1, 0)
writeRaster(SF2020young_restor_n_avoidboth, "rasters/PGM/input/SF2020young_restor_n_avoidboth.tif", format="GTiff", overwrite=T)

TF2020_restor_n_avoidboth <- sum(UPF2020_restor_n_avoidboth, DPF2020_restor_n_avoidboth)
TF2020_restor_n_avoidboth <- sum(TF2020_restor_n_avoidboth, SF2020young_restor_n_avoidboth)
writeRaster(TF2020_restor_n_avoidboth, "rasters/PGM/input/TF2020_restor_n_avoidboth.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_restor_n_avoidboth.px <- focal(TF2020_restor_n_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_restor_n_avoidboth.px)<-"TFpx"
TF2020_restor_n_avoidboth.px[is.nan(TF2020_restor_n_avoidboth.px)] <- 0
TF2020_restor_n_avoidboth.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_restor_n_avoidboth.px[])
#saving
writeRaster(TF2020_restor_n_avoidboth.px, "rasters/PGM/2020_restor_n_avoidboth/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_restor_n_avoidboth.ls <- focal(TF2020_restor_n_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_restor_n_avoidboth.ls)<-"TFls"
TF2020_restor_n_avoidboth.ls[is.nan(TF2020_restor_n_avoidboth.ls)] <- 0
TF2020_restor_n_avoidboth.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_restor_n_avoidboth.ls[])
#saving
writeRaster(TF2020_restor_n_avoidboth.ls, "rasters/PGM/2020_restor_n_avoidboth/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoidbothyoung", "TF2020_restor_n_avoidboth.px", "TF2020_restor_n_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature_restor_n_avoidboth <- SFAge2020_restor_n_avoidboth
SF2020mature_restor_n_avoidboth[] <- ifelse(SF2020mature_restor_n_avoidboth[]>5, 1, 0)
writeRaster(SF2020mature_restor_n_avoidboth, "rasters/PGM/input/SF2020mature_restor_n_avoidboth.tif", format="GTiff", overwrite=T)

MF2020_restor_n_avoidboth <- sum(UPF2020_restor_n_avoidboth, DPF2020_restor_n_avoidboth)
MF2020_restor_n_avoidboth <- sum(MF2020_restor_n_avoidboth, SF2020mature_restor_n_avoidboth)
writeRaster(MF2020_restor_n_avoidboth, "rasters/PGM/input/MF2020_restor_n_avoidboth.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_restor_n_avoidboth.px <- focal(MF2020_restor_n_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_restor_n_avoidboth.px)<-"MFpx"
MF2020_restor_n_avoidboth.px[is.nan(MF2020_restor_n_avoidboth.px)] <- 0
MF2020_restor_n_avoidboth.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_restor_n_avoidboth.px[])
#saving
writeRaster(MF2020_restor_n_avoidboth.px, "rasters/PGM/2020_restor_n_avoidboth/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_restor_n_avoidboth.ls <- focal(MF2020_restor_n_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_restor_n_avoidboth.ls)<-"MFls"
MF2020_restor_n_avoidboth.ls[is.nan(MF2020_restor_n_avoidboth.ls)] <- 0
MF2020_restor_n_avoidboth.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_restor_n_avoidboth.ls[])
#saving
writeRaster(MF2020_restor_n_avoidboth.ls, "rasters/PGM/2020_restor_n_avoidboth/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoidbothmature", "MF2020_restor_n_avoidboth.px", "MF2020_restor_n_avoidboth.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_restor_n_avoidboth <- raster("rasters/PGM/input/MF2020_restor_n_avoidboth.tif")
inv.MF2020_restor_n_avoidboth <- MF2020_restor_n_avoidboth
inv.MF2020_restor_n_avoidboth[inv.MF2020_restor_n_avoidboth==1]<-NA
#cheking
#inv.MF2020_restor_n_avoidboth
#plot(inv.MF2020_restor_n_avoidboth)

edge.dist.2020_restor_n_avoidboth <- distance(inv.MF2020_restor_n_avoidboth, doEdge=T)
names(edge.dist.2020_restor_n_avoidboth)<-"edgedist"
edge.dist.2020_restor_n_avoidboth[is.nan(edge.dist.2020_restor_n_avoidboth)] <- 0
edge.dist.2020_restor_n_avoidboth[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge.dist.2020_restor_n_avoidboth[])
#saving
writeRaster(edge.dist.2020_restor_n_avoidboth, "rasters/PGM/2020_restor_n_avoidboth/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_restor_n_avoidboth")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_restor_n_avoidboth <- edge.dist.2020_restor_n_avoidboth
edge2020_restor_n_avoidboth[] <- ifelse(edge2020_restor_n_avoidboth[] < 200, 0, ifelse(edge2020_restor_n_avoidboth[]>300, 0, 1))
writeRaster(edge2020_restor_n_avoidboth, "rasters/PGM/input/edge2020_restor_n_avoidboth.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_restor_n_avoidboth.px <- focal(edge2020_restor_n_avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_restor_n_avoidboth.px)<-"edgepx"
edge2020_restor_n_avoidboth.px[is.nan(edge2020_restor_n_avoidboth.px)] <- 0
edge2020_restor_n_avoidboth.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_restor_n_avoidboth.px[])
#saving
writeRaster(edge2020_restor_n_avoidboth.px, "rasters/PGM/2020_restor_n_avoidboth/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_restor_n_avoidboth.ls <- focal(edge2020_restor_n_avoidboth, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_restor_n_avoidboth.ls)<-"edgels"
edge2020_restor_n_avoidboth.ls[is.nan(edge2020_restor_n_avoidboth.ls)] <- 0
edge2020_restor_n_avoidboth.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_restor_n_avoidboth.ls[])
#saving
writeRaster(edge2020_restor_n_avoidboth.ls, "rasters/PGM/2020_restor_n_avoidboth/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_restor_n_avoidboth.px", "edge2020_restor_n_avoidboth.ls")]) #keeping only raster stack
gc()



#
#



## Scenario: Restoration and avoid both (Primary forests only) =================
### Undegraded primary forest
#UPF2020_restor_n_avoidboth2 <- raster("rasters/PGM/input/LULC/UPF2020_restor_n_avoidboth2.tif")

### mean upf cover in local scale (90m)
UPF2020_restor_n_avoidboth2.px <- focal(UPF2020_restor_n_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(UPF2020_restor_n_avoidboth2.px)<-"UPFpx"
UPF2020_restor_n_avoidboth2.px[is.nan(UPF2020_restor_n_avoidboth2.px)] <- 0
UPF2020_restor_n_avoidboth2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_restor_n_avoidboth2.px[])
#saving
writeRaster(UPF2020_restor_n_avoidboth2.px, "rasters/PGM/2020_restor_n_avoidboth2/UPFpx.tif", format="GTiff", overwrite=T)

# mean upf cover in landscape scale (~1000m)
UPF2020_restor_n_avoidboth2.ls <- focal(UPF2020_restor_n_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(UPF2020_restor_n_avoidboth2.ls)<-"UPFls"
UPF2020_restor_n_avoidboth2.ls[is.nan(UPF2020_restor_n_avoidboth2.ls)] <- 0
UPF2020_restor_n_avoidboth2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, UPF2020_restor_n_avoidboth2.ls[])
#saving
writeRaster(UPF2020_restor_n_avoidboth2.ls, "rasters/PGM/2020_restor_n_avoidboth2/UPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("UPF2020_restor_n_avoidboth2.px", "UPF2020_restor_n_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




### Degraded primary forest
#uDPF2020_restor_n_avoidboth2 <- raster("rasters/PGM/input/LULC/uDPF2020_restor_n_avoidboth2.tif")
#RDPF2020_restor_n_avoidboth2 <- raster("rasters/PGM/input/LULC/RDPF2020_restor_n_avoidboth2.tif")
DPF2020_restor_n_avoidboth2 <- sum(uDPF2020_restor_n_avoidboth2, RDPF2020_restor_n_avoidboth2)

### mean dpf cover in local scale (90m)
DPF2020_restor_n_avoidboth2.px <- focal(DPF2020_restor_n_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(DPF2020_restor_n_avoidboth2.px)<-"DPFpx"
DPF2020_restor_n_avoidboth2.px[is.nan(DPF2020_restor_n_avoidboth2.px)] <- 0
DPF2020_restor_n_avoidboth2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_restor_n_avoidboth2.px[])
#saving
writeRaster(DPF2020_restor_n_avoidboth2.px, "rasters/PGM/2020_restor_n_avoidboth2/DPFpx.tif", format="GTiff", overwrite=T)

# mean dpf cover in landscape scale (~1000m)
DPF2020_restor_n_avoidboth2.ls <- focal(DPF2020_restor_n_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(DPF2020_restor_n_avoidboth2.ls)<-"DPFls"
DPF2020_restor_n_avoidboth2.ls[is.nan(DPF2020_restor_n_avoidboth2.ls)] <- 0
DPF2020_restor_n_avoidboth2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, DPF2020_restor_n_avoidboth2.ls[])
#saving
writeRaster(DPF2020_restor_n_avoidboth2.ls, "rasters/PGM/2020_restor_n_avoidboth2/DPFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("DPF2020_restor_n_avoidboth2.px", "DPF2020_restor_n_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




## Time since degradation
TSD2020_restor_n_avoidboth2 <- TSD2020_avoidboth2
writeRaster(TSD2020_restor_n_avoidboth2, "rasters/PGM/input/TSD2020_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)

### mean time in local scale (90m)
TSD2020_restor_n_avoidboth2.px <- focal(TSD2020_restor_n_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TSD2020_restor_n_avoidboth2.px)<-"TSDpx"
TSD2020_restor_n_avoidboth2.px[is.nan(TSD2020_restor_n_avoidboth2.px)] <- 0
TSD2020_restor_n_avoidboth2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_restor_n_avoidboth2.px[])
#saving
writeRaster(TSD2020_restor_n_avoidboth2.px, "rasters/PGM/2020_restor_n_avoidboth2/TSDpx.tif", format="GTiff", overwrite=T)

# mean time in landscape scale (~1000m)
TSD2020_restor_n_avoidboth2.ls <- focal(TSD2020_restor_n_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TSD2020_restor_n_avoidboth2.ls)<-"TSDls"
TSD2020_restor_n_avoidboth2.ls[is.nan(TSD2020_restor_n_avoidboth2.ls)] <- 0
TSD2020_restor_n_avoidboth2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TSD2020_restor_n_avoidboth2.ls[])
#saving
writeRaster(TSD2020_restor_n_avoidboth2.ls, "rasters/PGM/2020_restor_n_avoidboth2/TSDls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("TSD2020_restor_n_avoidboth2.px", "TSD2020_restor_n_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest
#uSF2020_restor_n_avoidboth2 <- raster("rasters/PGM/input/LULC/uSF2020_restor_n_avoidboth2.tif")
#DSF2020_restor_n_avoidboth2 <- raster("rasters/PGM/input/LULC/DSF2020_restor_n_avoidboth2.tif")
SF2020_restor_n_avoidboth2 <- sum(uSF2020_restor_n_avoidboth2, DSF2020_restor_n_avoidboth2)

### mean sf cover in local scale (90m)
SF2020_restor_n_avoidboth2.px <- focal(SF2020_restor_n_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SF2020_restor_n_avoidboth2.px)<-"SFpx"
SF2020_restor_n_avoidboth2.px[is.nan(SF2020_restor_n_avoidboth2.px)] <- 0
SF2020_restor_n_avoidboth2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_restor_n_avoidboth2.px[])
#saving
writeRaster(SF2020_restor_n_avoidboth2.px, "rasters/PGM/2020_restor_n_avoidboth2/SFpx.tif", format="GTiff", overwrite=T)

# mean sf cover in landscape scale (~1000m)
SF2020_restor_n_avoidboth2.ls <- focal(SF2020_restor_n_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SF2020_restor_n_avoidboth2.ls)<-"SFls"
SF2020_restor_n_avoidboth2.ls[is.nan(SF2020_restor_n_avoidboth2.ls)] <- 0
SF2020_restor_n_avoidboth2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SF2020_restor_n_avoidboth2.ls[])
#saving
writeRaster(SF2020_restor_n_avoidboth2.ls, "rasters/PGM/2020_restor_n_avoidboth2/SFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoidboth2.px", "SF2020_restor_n_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




## Secondary forest age
SFAge2020_restor_n_avoidboth2 <- SFAge2020_avoidboth2
SFAge2020_restor_n_avoidboth2[] <- ifelse(SFAge2020_restor_n_avoidboth2[]==0 & SF2020_restor_n_avoidboth2[]==1, 10, SFAge2020_restor_n_avoidboth2[])
writeRaster(SFAge2020_restor_n_avoidboth2, "rasters/PGM/input/SFAge2020_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)

### mean age in local scale (90m)
SFAge2020_restor_n_avoidboth2.px <- focal(SFAge2020_restor_n_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(SFAge2020_restor_n_avoidboth2.px)<-"SFAgepx"
SFAge2020_restor_n_avoidboth2.px[is.nan(SFAge2020_restor_n_avoidboth2.px)] <- 0
SFAge2020_restor_n_avoidboth2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_restor_n_avoidboth2.px[])
#saving
writeRaster(SFAge2020_restor_n_avoidboth2.px, "rasters/PGM/2020_restor_n_avoidboth2/SFAgepx.tif", format="GTiff", overwrite=T)

# mean age in landscape scale (~1000m)
SFAge2020_restor_n_avoidboth2.ls <- focal(SFAge2020_restor_n_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(SFAge2020_restor_n_avoidboth2.ls)<-"SFAgels"
SFAge2020_restor_n_avoidboth2.ls[is.nan(SFAge2020_restor_n_avoidboth2.ls)] <- 0
SFAge2020_restor_n_avoidboth2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, SFAge2020_restor_n_avoidboth2.ls[])
#saving
writeRaster(SFAge2020_restor_n_avoidboth2.ls, "rasters/PGM/2020_restor_n_avoidboth2/SFAgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SFAge2020_restor_n_avoidboth2.px", "SFAge2020_restor_n_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




## Total forest (UPF + DPF + SF>2yr)
SF2020young_restor_n_avoidboth2 <- SFAge2020_restor_n_avoidboth2
SF2020young_restor_n_avoidboth2[] <- ifelse(SF2020young_restor_n_avoidboth2[]>2, 1, 0)
writeRaster(SF2020young_restor_n_avoidboth2, "rasters/PGM/input/SF2020young_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)

TF2020_restor_n_avoidboth2 <- sum(UPF2020_restor_n_avoidboth2, DPF2020_restor_n_avoidboth2)
TF2020_restor_n_avoidboth2 <- sum(TF2020_restor_n_avoidboth2, SF2020young_restor_n_avoidboth2)
writeRaster(TF2020_restor_n_avoidboth2, "rasters/PGM/input/TF2020_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)

### mean tf cover in local scale (90m)
TF2020_restor_n_avoidboth2.px <- focal(TF2020_restor_n_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(TF2020_restor_n_avoidboth2.px)<-"TFpx"
TF2020_restor_n_avoidboth2.px[is.nan(TF2020_restor_n_avoidboth2.px)] <- 0
TF2020_restor_n_avoidboth2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_restor_n_avoidboth2.px[])
#saving
writeRaster(TF2020_restor_n_avoidboth2.px, "rasters/PGM/2020_restor_n_avoidboth2/TFpx.tif", format="GTiff", overwrite=T)

# mean tf cover in landscape scale (~1000m)
TF2020_restor_n_avoidboth2.ls <- focal(TF2020_restor_n_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(TF2020_restor_n_avoidboth2.ls)<-"TFls"
TF2020_restor_n_avoidboth2.ls[is.nan(TF2020_restor_n_avoidboth2.ls)] <- 0
TF2020_restor_n_avoidboth2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, TF2020_restor_n_avoidboth2.ls[])
#saving
writeRaster(TF2020_restor_n_avoidboth2.ls, "rasters/PGM/2020_restor_n_avoidboth2/TFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoidboth2young", "TF2020_restor_n_avoidboth2.px", "TF2020_restor_n_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




## Total mature forest (UPF + DPF + SF>5yr)
SF2020mature_restor_n_avoidboth2 <- SFAge2020_restor_n_avoidboth2
SF2020mature_restor_n_avoidboth2[] <- ifelse(SF2020mature_restor_n_avoidboth2[]>5, 1, 0)
writeRaster(SF2020mature_restor_n_avoidboth2, "rasters/PGM/input/SF2020mature_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)

MF2020_restor_n_avoidboth2 <- sum(UPF2020_restor_n_avoidboth2, DPF2020_restor_n_avoidboth2)
MF2020_restor_n_avoidboth2 <- sum(MF2020_restor_n_avoidboth2, SF2020mature_restor_n_avoidboth2)
writeRaster(MF2020_restor_n_avoidboth2, "rasters/PGM/input/MF2020_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)

### mean mf cover in local scale (90m)
MF2020_restor_n_avoidboth2.px <- focal(MF2020_restor_n_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(MF2020_restor_n_avoidboth2.px)<-"MFpx"
MF2020_restor_n_avoidboth2.px[is.nan(MF2020_restor_n_avoidboth2.px)] <- 0
MF2020_restor_n_avoidboth2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_restor_n_avoidboth2.px[])
#saving
writeRaster(MF2020_restor_n_avoidboth2.px, "rasters/PGM/2020_restor_n_avoidboth2/MFpx.tif", format="GTiff", overwrite=T)

# mean mf cover in landscape scale (~1000m)
MF2020_restor_n_avoidboth2.ls <- focal(MF2020_restor_n_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(MF2020_restor_n_avoidboth2.ls)<-"MFls"
MF2020_restor_n_avoidboth2.ls[is.nan(MF2020_restor_n_avoidboth2.ls)] <- 0
MF2020_restor_n_avoidboth2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, MF2020_restor_n_avoidboth2.ls[])
#saving
writeRaster(MF2020_restor_n_avoidboth2.ls, "rasters/PGM/2020_restor_n_avoidboth2/MFls.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("SF2020_restor_n_avoidboth2mature", "MF2020_restor_n_avoidboth2.px", "MF2020_restor_n_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




# [edgedist] distance to forest edge
#' this variable is the euclidean distance between mature forest 
#' to the nearest cell that is not MF
#MF2020_restor_n_avoidboth2 <- raster("rasters/PGM/input/MF2020_restor_n_avoidboth2.tif")
inv.MF2020_restor_n_avoidboth2 <- MF2020_restor_n_avoidboth2
inv.MF2020_restor_n_avoidboth2[inv.MF2020_restor_n_avoidboth2==1]<-NA
#cheking
#inv.MF2020_restor_n_avoidboth2
#plot(inv.MF2020_restor_n_avoidboth2)

edge.dist.2020_restor_n_avoidboth2 <- distance(inv.MF2020_restor_n_avoidboth2, doEdge=T)
names(edge.dist.2020_restor_n_avoidboth2)<-"edgedist"
edge.dist.2020_restor_n_avoidboth2[is.nan(edge.dist.2020_restor_n_avoidboth2)] <- 0
edge.dist.2020_restor_n_avoidboth2[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge.dist.2020_restor_n_avoidboth2[])
#saving
writeRaster(edge.dist.2020_restor_n_avoidboth2, "rasters/PGM/2020_restor_n_avoidboth2/edgedist.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("inv.MF2020_restor_n_avoidboth2")])
gc()



#
#




# Forest edge
# marking edge and core areas
edge2020_restor_n_avoidboth2 <- edge.dist.2020_restor_n_avoidboth2
edge2020_restor_n_avoidboth2[] <- ifelse(edge2020_restor_n_avoidboth2[] < 200, 0, ifelse(edge2020_restor_n_avoidboth2[]>300, 0, 1))
writeRaster(edge2020_restor_n_avoidboth2, "rasters/PGM/input/edge2020_restor_n_avoidboth2.tif", format="GTiff", overwrite=T)

### mean edge cover in local scale (90m)
edge2020_restor_n_avoidboth2.px <- focal(edge2020_restor_n_avoidboth2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
names(edge2020_restor_n_avoidboth2.px)<-"edgepx"
edge2020_restor_n_avoidboth2.px[is.nan(edge2020_restor_n_avoidboth2.px)] <- 0
edge2020_restor_n_avoidboth2.px[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_restor_n_avoidboth2.px[])
#saving
writeRaster(edge2020_restor_n_avoidboth2.px, "rasters/PGM/2020_restor_n_avoidboth2/edgepx.tif", format="GTiff", overwrite=T)

# mean edge cover in landscape scale (~1000m)
edge2020_restor_n_avoidboth2.ls <- focal(edge2020_restor_n_avoidboth2, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
names(edge2020_restor_n_avoidboth2.ls)<-"edgels"
edge2020_restor_n_avoidboth2.ls[is.nan(edge2020_restor_n_avoidboth2.ls)] <- 0
edge2020_restor_n_avoidboth2.ls[] <- ifelse(pgm.lulc[[1]][]==0, NA, edge2020_restor_n_avoidboth2.ls[])
#saving
writeRaster(edge2020_restor_n_avoidboth2.ls, "rasters/PGM/2020_restor_n_avoidboth2/edgels.tif", format="GTiff", overwrite=T)

rm(list=ls()[ls() %in% c("edge2020_restor_n_avoidboth2.px", "edge2020_restor_n_avoidboth2.ls")]) #keeping only raster stack
gc()



#
#




## import rural properties shapefiles and data from SISCAR
#' https://www.car.gov.br/publico/municipios/downloads
#' see the "restoration_candidate_pgm.R" on the script folder
#' for details on reducing overlap between properties and 
pgm.car <- readOGR(dsn = "rasters/PGM/raw", layer = "pgm_car_after_restoration_2010")
names(pgm.car@data) <- c("COD_IMOVEL","NUM_AREA","COD_ESTADO","NOM_MUNICI","NUM_MODULO","TIPO_IMOVE","SITUACAO","CONDICAO_I",
                         "num_area_ha","num_modulo_new","num_area_flag","FOREST_COVER_2007","FOREST_COVER_2007_PP","FOREST_COVER_2010","FOREST_COVER_2010_PP",
                         "NEED_INCREMENT","APP_FOREST_COVER_INCREMENT","APP_FOREST_COVER_INCREMENT_PP","NEED_INCREMENT_AFTER_APP",
                         "ARL_FOREST_COVER_INCREMENT","ARL_FOREST_COVER_INCREMENT_PP","NEED_INCREMENT_AFTER_ARL")
pgm.car <- spTransform(pgm.car, crs(std.proj))

#checking
#nrow(pgm.car)
#nrow(pgm.car@data[pgm.car@data$NEED_INCREMENT==1,])
#head(pgm.car@data[,"num_area_ha"])

#note:
#according to Brazilian Forest Code
#small properties have less than or equal to 4 fiscal modules
#medium/big properties heave more than 4 fiscal modules
#in PGM, 1 fiscal module = 55ha [or 550000m2]


## [propertysize]
# transforming the properties polygons into rasters
# cell value is the property size
pgm.car.raster <- rasterize(pgm.car, pgm.lulc[[1]], field = "num_area_ha", fun = "mean", background = 0)
pgm.car.raster[] <- ifelse(pgm.lulc[[1]][]==0, NA, pgm.car.raster[])
#checking
#st_crs(pgm.car.raster)==st_crs(pgm.shp)
#plot(pgm.car.raster)
#plot(pgm.car, add=T)

#saving
writeRaster(pgm.car.raster, "rasters/PGM/2010_real/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(pgm.car.raster, "rasters/PGM/2020_real/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(pgm.car.raster, "rasters/PGM/2020_avoiddeforest/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(pgm.car.raster, "rasters/PGM/2020_avoiddeforest2/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(pgm.car.raster, "rasters/PGM/2020_avoiddegrad/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(pgm.car.raster, "rasters/PGM/2020_avoiddegrad2/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(pgm.car.raster, "rasters/PGM/2020_avoidboth/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(pgm.car.raster, "rasters/PGM/2020_avoidboth2/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(pgm.car.raster, "rasters/PGM/2020_restor_wo_avoid/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(pgm.car.raster, "rasters/PGM/2020_restor_n_avoiddeforest/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(pgm.car.raster, "rasters/PGM/2020_restor_n_avoiddeforest2/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(pgm.car.raster, "rasters/PGM/2020_restor_n_avoidboth/propertysize.tif", format="GTiff", overwrite=T)
writeRaster(pgm.car.raster, "rasters/PGM/2020_restor_n_avoidboth2/propertysize.tif", format="GTiff", overwrite=T)



#
#



## Fixed costs =================================================================
## avoid degradation costs -- adding fire control costs
#' creating and maintaining (i.e., every four years on average) 
#' 6m wide fire breaks at the edge of forested areas

#calculating forest perimeter in each property
#adding variable for forest cover
pgm.car@data$FOREST_PERIMETER <- NA

j=nrow(pgm.car@data)
for (i in pgm.car$COD_IMOVEL) {
  
  rural.property <- pgm.car[pgm.car$COD_IMOVEL==i,]
  
  rural.property.edge <- crop(pgm.lulc.2010.forest.class, extent(rural.property))
  rural.property.edge <- mask(rural.property.edge, rural.property)
  
  if(all(is.na(values(rural.property.edge)))) next
  #convert raster to polygons
  rural.property.edge.shp <- as_Spatial(st_as_sf(st_as_stars(rural.property.edge),
                                                 as_points = FALSE, merge = TRUE))
  #cheking & adjupgments
  #st_crs(forest.cover.shp)==st_crs(pgm.shp)
  #gIsValid(forest.cover.shp)
  #FALSE here means that you'll need to run the buffer routine:
  #forest.cover.shp <- rgeos::gBuffer(forest.cover.shp, byid = TRUE, width = 0)
  
  rural.property.edge.shp <- gDifference(as(rural.property.edge.shp,"SpatialLines"),
                                         as(gUnaryUnion(rural.property.edge.shp),"SpatialLines"),
                                         byid=TRUE)
  
  if(is.null(rural.property.edge.shp)) next
  
  #estimating the perimeter
  pgm.car@data[pgm.car$COD_IMOVEL==i,"FOREST_PERIMETER"] <- gLength(rural.property.edge.shp)*111111
  
  
  j=j-1
  cat("\n>", j, "out of", nrow(pgm.car@data), "properties left<\n")
  
}


#' @description fire breaks could be cleared at rate of 33.333 meters per day, costing R$100 per day according to IPAM
#' source: https://ipam.org.br/wp-content/uploads/2009/05/te%CC%81cnicas_de_prevenc%CC%A7a%CC%83o_de_fogo_acidental_.pdf
#' cost of fire control was (P/33.33333) x 100, where P is the perimeter in meters of forested area in the property
pgm.car@data$cost <- (as.numeric((pgm.car@data$FOREST_PERIMETER/33.33333) * 100))
pgm.car@data$cost <- ifelse(is.na(pgm.car@data$cost), 0, pgm.car@data$cost)

#calculating the costs per area per year
pgm.car@data$final_cost <- (pgm.car@data$cost/pgm.car@data$num_area_ha)/4

#convert to raster
avoid.degrad.cost <- rasterize(pgm.car, pgm.lulc.2010.forest.class, field = "final_cost", fun = mean)
avoid.degrad.cost[] <- ifelse(pgm.lulc[[1]][]==0, NA, avoid.degrad.cost[])
#plot(avoid.degrad.cost)

#saving
writeRaster(avoid.degrad.cost, "models.output/costs/PGM_2010_real_base_firecontrol.tif", format="GTiff", overwrite=T)



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

restor.cost1 <- pgm.lulc[["pgm.lulc.2010real"]]

values(restor.cost1)[values(restor.cost1) %in% agro.class] <- 1
values(restor.cost1)[values(restor.cost1) > 1] <- 0
names(restor.cost1) <- "restoration.no.fences"


restor.cost1.deforest.dist <- pgm.deforest.dist.copy
values(restor.cost1.deforest.dist)[values(restor.cost1.deforest.dist) == 0] <- NA
values(restor.cost1.deforest.dist)[values(restor.cost1.deforest.dist) > 500] <- NA
values(restor.cost1.deforest.dist)[values(restor.cost1.deforest.dist) <= 500] <- 1

restor.cost1 <- mask(restor.cost1, restor.cost1.deforest.dist)
restor.cost1[is.na(restor.cost1)] <- 0
restor.cost1[restor.cost1==1] <- 189.13 * 0.09

#natural regeneration with fence
restor.cost2 <- pgm.lulc[["pgm.lulc.2010real"]]

restor.cost2[restor.cost2 == 15] <- 1
restor.cost2[restor.cost2 > 1] <- 0
names(restor.cost2) <- "restoration.fences"


restor.cost2 <- mask(restor.cost2, restor.cost1.deforest.dist)
restor.cost2[is.na(restor.cost2)] <- 0
restor.cost2[restor.cost2==1] <- 1331.55 * 0.09

#active restoration
restor.cost3 <- pgm.lulc[["pgm.lulc.2010real"]]

values(restor.cost3)[values(restor.cost3) %in% deforestation.class.list] <- 1
values(restor.cost3)[values(restor.cost3) > 1] <- 0
names(restor.cost3) <- "restoration.active"


restor.cost3.deforest.dist <- pgm.deforest.dist.copy
values(restor.cost3.deforest.dist)[values(restor.cost3.deforest.dist) <= 500] <- NA
values(restor.cost3.deforest.dist)[values(restor.cost3.deforest.dist) > 500] <- 1

restor.cost3 <- mask(restor.cost3, restor.cost3.deforest.dist)
restor.cost3[is.na(restor.cost3)] <- 0
restor.cost3[restor.cost3==1] <- 7899.71 * 0.09

#restoration cost layer
restor.cost.final <- sum(restor.cost1, restor.cost2, restor.cost3)

candidate.areas.final.mask <- candidate.areas.final
candidate.areas.final.mask[candidate.areas.final.mask==0] <- NA

restor.cost.final <- mask(restor.cost.final, candidate.areas.final.mask)

restor.cost.final[is.na(restor.cost.final)] <- 0
restor.cost.final[] <- ifelse(pgm.lulc[[1]][]==0, NA, restor.cost.final[])

writeRaster(restor.cost.final, paste0("models.output/costs/PGM_2010_real_base_restoration.tif"), format = "GTiff", overwrite = T)



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




# end =================================










































































#==============================| previous modeling approach


#### time since degradation 2010 data from RAS 
#### quantitative comparison of manual inspection of satellite images [150m resolution]
#### and field observations done by two observers (TG and SN)
#### see RAS environmental explanatory variable guideline document for details
#### 2010-2020 data from DETER
#### 
###pgm.degrad.2010 <- raster("rasters/PGM/raw/pgm-2010-deg_tsince0_150m.grd")
###
#### Conversion of rasters into same extent
###pgm.degrad.2010 <- projectRaster(pgm.degrad.2010, crs = std.proj, method='ngb')
###pgm.degrad.2010 <- resample(pgm.degrad.2010, pgm.lulc.2010.forest.class, method='ngb')
###
####checking
####st_crs(pgm.degrad.2010)==st_crs(pgm.shp)
####plot(pgm.degrad.2010)
####range(values(pgm.degrad.2010), na.rm=T)
###
###
####excluding non-forest areas
###pgm.degrad.2010 <- mask(pgm.degrad.2010, pgm.lulc.2010.forest.mask)
###pgm.degrad.2010[is.na(pgm.degrad.2010)] <- 0
###names(pgm.degrad.2010) <- c("pgm.degrad.2010real")
###
#### calculating time since degradation for 2020
###pgm.degrad.temp <- pgm.degrad.2010
###
####creating mask to detect repeated degradation
###pgm.repeateddegrad.mask <- pgm.degrad.temp
###pgm.repeateddegrad.mask <- mask(pgm.repeateddegrad.mask, pgm.lulc.2010.forest.mask)
###pgm.repeateddegrad.mask[pgm.repeateddegrad.mask[]<24] <- 1000
###pgm.repeateddegrad.mask[pgm.repeateddegrad.mask[]!=1000] <- NA
###
###
#### deter data between 2011 and 2015
###deter.2011.15 <- load_degrad(dataset = "degrad", raw_data = T, time_period = 2011:2015)
###
###for (year in 1:5) {   #1=2011; 5=2015
###  
###  #selecting year-by-year
###  deter.yearx <- deter.2011.15[[year]] 
###  deter.yearx <- sf:::as_Spatial(deter.yearx$geometry)
###  
###  #croping to study area
###  pgm.deter.yearx <- crop(deter.yearx, extent(pgm.lulc.2010.forest.mask))
###  
###  #converting to raster
###  pgm.deter.yearx <- rasterize(pgm.deter.yearx, pgm.lulc.2010.forest.mask, field=1000)
###  pgm.deter.yearx[is.na(pgm.deter.yearx)]<-0
###  
###  pgm.deter.yearx <- mask(pgm.deter.yearx, pgm.lulc.2010.forest.mask)
###  
###  #count sites with repeated degradation
###  pgm.repeateddegrad.mask <- sum(pgm.repeateddegrad.mask, pgm.deter.yearx, na.rm = T)
###  pgm.repeateddegrad.mask[pgm.repeateddegrad.mask[]<1000] <- NA
###  
###  
###  #counting time since degradation
###  pgm.degrad.temp <- pgm.degrad.temp+1
###  pgm.degrad.temp <- sum(pgm.degrad.temp, pgm.deter.yearx, na.rm = T)
###  pgm.degrad.temp <- mask(pgm.degrad.temp, pgm.lulc.2010.forest.mask)
###  pgm.degrad.temp[pgm.degrad.temp[]>1000] <- 0
###  
###  cat("\n> year", year, "done! <\n")
###}
###
###
####par(mfrow=c(1,3))
####plot(pgm.degrad.temp)
####plot(pgm.deter.yearx, col=c("white", "black"))
###
###
#### deter data between 2016 and 2020
###deter.2016.20 <- readOGR(dsn = "rasters/PGM/raw", layer = "deter_public")
###pgm.deter.2016.20 <- crop(deter.2016.20, extent(pgm.degrad.temp))
###degradation_cat <- c('CICATRIZ_DE_QUEIMADA', 'CS_DESORDENADO', 'CS_GEOMETRICO') #'DEGRADACAO'
###pgm.deter.2016.20 <- pgm.deter.2016.20[pgm.deter.2016.20$CLASSNAME %in% degradation_cat,]
###
###rm(deter.2016.20); gc()
###
###for (year in 2016:2020) {
###  
###  #selecting year-by-year
###  pgm.deter.yearx <- pgm.deter.2016.20[grep(year, pgm.deter.2016.20$VIEW_DATE),]
###  
###  #converting to raster
###  pgm.deter.yearx <- rasterize(pgm.deter.yearx, pgm.lulc.2010.forest.mask, field=1000)
###  pgm.deter.yearx[is.na(pgm.deter.yearx)]<-0
###  
###  pgm.deter.yearx <- mask(pgm.deter.yearx, pgm.lulc.2010.forest.mask)
###  
###  #count sites with repeated degradation
###  pgm.repeateddegrad.mask <- sum(pgm.repeateddegrad.mask, pgm.deter.yearx, na.rm = T)
###  pgm.repeateddegrad.mask[pgm.repeateddegrad.mask[]<1000] <- NA
###  
###  
###  #counting time since degradation
###  pgm.degrad.temp <- pgm.degrad.temp+1
###  pgm.degrad.temp <- sum(pgm.degrad.temp, pgm.deter.yearx, na.rm = T)
###  pgm.degrad.temp <- mask(pgm.degrad.temp, pgm.lulc.2010.forest.mask)
###  pgm.degrad.temp[pgm.degrad.temp[]>1000] <- 0
###  
###  cat("\n> year", year, "done! <\n")
###}
###
###
####creating mask to detect repeated degradation
###pgm.repeateddegrad.mask <- mask(pgm.repeateddegrad.mask, pgm.lulc.2020.forest.mask)
###pgm.repeateddegrad.mask[pgm.repeateddegrad.mask[]>1000] <- 1
###pgm.repeateddegrad.mask[pgm.repeateddegrad.mask[]!=1] <- NA
###
###
###pgm.degrad.2020 <- pgm.degrad.temp
###
####excluding non-forest areas
###pgm.degrad.2020 <- mask(pgm.degrad.2020, pgm.lulc.2020.forest.mask)
###pgm.degrad.2020[is.na(pgm.degrad.2020)] <- 0
###names(pgm.degrad.2020) <- "pgm.degrad.2020real"
###
####writeRaster(pgm.degrad.2020, "rasters/PGM/raw/pgm-2020-deg_tsince0.tif", format = "GTiff", overwrite = T)
####pgm.degrad.2020 <- raster("rasters/PGM/raw/pgm-2020-deg_tsince0.tif")
###
###
###
###pgm.degrad <- stack(pgm.degrad.2010, pgm.degrad.2020)
####checking
####st_crs(pgm.degrad)==st_crs(pgm.shp)
####plot(pgm.degrad)
####range(values(pgm.degrad[["pgm.degrad.2010real"]]), na.rm=T)
###
####non-degraded sites will be considered with 300 years following (BIB)
###pgm.degrad[["pgm.degrad.2010real"]][pgm.degrad[["pgm.degrad.2010real"]]>23] <- 300
###pgm.degrad[["pgm.degrad.2020real"]][pgm.degrad[["pgm.degrad.2020real"]]>33] <- 300
###
#### isolating degraded primary forest and degraded secondary forest class pixels
###pgm.degrad.2010.forest.class <- pgm.degrad[["pgm.degrad.2010real"]]
###pgm.degrad.2010.forest.class <- mask(pgm.degrad.2010.forest.class, pgm.lulc.2010.forest.mask)
###pgm.degrad.2010.forest.class[pgm.degrad.2010.forest.class>23]<-NA
###pgm.degrad.2010.forest.class[pgm.degrad.2010.forest.class<=23]<-1
###
###pgm.degrad.2010.mask <- pgm.degrad.2010.forest.class
###
###pgm.degrad.2010.forest.class[is.na(pgm.degrad.2010.forest.class)]<-0
###
###
###
###pgm.degrad.2020.forest.class <- pgm.degrad[["pgm.degrad.2020real"]]
###pgm.degrad.2020.forest.class <- mask(pgm.degrad.2020.forest.class, pgm.lulc.2020.forest.mask)
###pgm.degrad.2020.forest.class[pgm.degrad.2020.forest.class>33]<-NA
###pgm.degrad.2020.forest.class[pgm.degrad.2020.forest.class<=33]<-1
###
###pgm.degrad.2020.mask <- pgm.degrad.2020.forest.class
###
###pgm.degrad.2020.forest.class[is.na(pgm.degrad.2020.forest.class)]<-0
###
###
###rm(list=ls()[ls() %in% c("deter.2011.15", "pgm.deter.2016.20", "degradation_cat", "pgm.degrad.2010", "pgm.degrad.2020",
###                         "pgm.degrad.temp", "deter.yearx", "year")])
###gc()
###
###
###
####
####
#
#
#
#
#
#
########################################################################################################################.
########################################.
##### setting explanatory variables ####.
########################################.
#####         2010 Real             ####.
########################################.
#
## [SFage] secondary forest age -- pixel: 3x3 (60m); and landscape: 35x35 (1000m)
## this variable is the mean age of secondary forest
#
#SFage2010 <- pgm.sfage[["pgm.sfage.2010real"]]
##plot(SFage2010)
#
##saving
#writeRaster(SFage2010, "rasters/PGM/input/SFage2010_real.tif", format="GTiff", overwrite=T)
#
#
## mean sf age in pixel scale (60m)
#SFage2010.px <- focal(SFage2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##SFage2010.px
##anyNA(SFage2010.px[])
##plot(SFage2010.px)
#
#names(SFage2010.px)<-"SFagepx"
#SFage2010.px[is.nan(SFage2010.px)] <- 0
##SFage2010.px <- mask(SFage2010.px, pgm.shp)
#
##saving
#writeRaster(SFage2010.px, "rasters/PGM/2010_real/SFagepx.tif", format="GTiff", overwrite=T)
#
## mean sf age in landscape scale (1000m)
#SFage2010.ls <- focal(SFage2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
###cheking
##SFage2010.ls
##anyNA(SFage2010.ls[])
##plot(SFage2010.ls)
#
#names(SFage2010.ls)<-"SFagels"
#SFage2010.ls[is.nan(SFage2010.ls)] <- 0
##SFage2010.ls <- mask(SFage2010.ls, pgm.shp)
#
##saving
#writeRaster(SFage2010.ls, "rasters/PGM/2010_real/SFagels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("SFage2010.px", "SFage2010.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [SF] secondary forest
## this variable includes forest pixels in SFage raster
## which has less than 25 years for 2010 or less than 35 for 2020
#
#SF2010 <- pgm.sfage.2010.all.class
##plot(SF2010)
#
##saving
#writeRaster(SF2010, "rasters/PGM/input/SF2010_real.tif", format="GTiff", overwrite=T)
#
#
## mean sf cover in pixel scale (60m)
#SF2010.px <- focal(SF2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##SF2010.px
##anyNA(SF2010.px[])
##plot(SF2010.px)
#
#names(SF2010.px)<-"SFpx"
#SF2010.px[is.nan(SF2010.px)] <- 0
##SF2010.px <- mask(SF2010.px, pgm.shp)
#
##saving
#writeRaster(SF2010.px, "rasters/PGM/2010_real/SFpx.tif", format="GTiff", overwrite=T)
#
## mean sf cover in landscape scale (1000m)
#SF2010.ls <- focal(SF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
###cheking
##SF2010.ls
##anyNA(SF2010.ls[])
##plot(SF2010.ls)
#
#names(SF2010.ls)<-"SFls"
#SF2010.ls[is.nan(SF2010.ls)] <- 0
##SF2010.ls <- mask(SF2010.ls, pgm.shp)
#
##saving
#writeRaster(SF2010.ls, "rasters/PGM/2010_real/SFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("SF2010.px", "SF2010.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [TSD] time since degradation
## this variable is the mean time since a degradation event
#
#TSD2010 <- pgm.degrad[["pgm.degrad.2010real"]]
##plot(TSD2010)
#
##saving
#writeRaster(TSD2010, "rasters/PGM/input/TSD2010_real.tif", format="GTiff", overwrite=T)
#
#
## mean tsd cover in pixel scale (60m)
#TSD2010.px <- focal(TSD2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TSD2010.px
##anyNA(TSD2010.px[])
##plot(TSD2010.px)
#
#names(TSD2010.px)<-"TSDpx"
#TSD2010.px[is.nan(TSD2010.px)] <- 0
##TSD2010.px <- mask(TSD2010.px, pgm.shp)
#
##saving
#writeRaster(TSD2010.px, "rasters/PGM/2010_real/TSDpx.tif", format="GTiff", overwrite=T)
#
## mean tsd cover in landscape scale (1000m)
#TSD2010.ls <- focal(TSD2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
###cheking
##TSD2010.ls
##anyNA(TSD2010.ls[])
##plot(TSD2010.ls)
#
#names(TSD2010.ls)<-"TSDls"
#TSD2010.ls[is.nan(TSD2010.ls)] <- 0
##TSD2010.ls <- mask(TSD2010.ls, pgm.shp)
#
##saving
#writeRaster(TSD2010.ls, "rasters/PGM/2010_real/TSDls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TSD2010.px", "TSD2010.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [DPF] degraded primary forest
## this variable includes forest pixels in TSD raster
## which has less than 24 years for 2010 or less than 34 for 2020
#
#DPF2010 <- pgm.degrad.2010.forest.class
##plot(DPF2010)
#
##saving
#writeRaster(DPF2010, "rasters/PGM/input/DPF2010_real.tif", format="GTiff", overwrite=T)
#
#
## mean dpf cover in pixel scale (60m)
#DPF2010.px <- focal(DPF2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##DPF2010.px
##anyNA(DPF2010.px[])
##plot(DPF2010.px)
#
#names(DPF2010.px)<-"DPFpx"
#DPF2010.px[is.nan(DPF2010.px)] <- 0
##DPF2010.px <- mask(DPF2010.px, pgm.shp)
#
##saving
#writeRaster(DPF2010.px, "rasters/PGM/2010_real/DPFpx.tif", format="GTiff", overwrite=T)
#
## mean dpf cover in landscape scale (1000m)
#DPF2010.ls <- focal(DPF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
###cheking
##DPF2010.ls
##anyNA(DPF2010.ls[])
##plot(DPF2010.ls)
#
#names(DPF2010.ls)<-"DPFls"
#DPF2010.ls[is.nan(DPF2010.ls)] <- 0
##DPF2010.ls <- mask(DPF2010.ls, pgm.shp)
#
##saving
#writeRaster(DPF2010.ls, "rasters/PGM/2010_real/DPFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("DPF2010.px", "DPF2010.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [UPF] undisturbed primary forest
## this variable includes all forest pixels in LULC raster (value == 3)
## excluding degraded and secondary forest
#
#UPF2010<-pgm.lulc.2010.forest.class
#UPF2010<-mask(UPF2010, pgm.sfage.2010.mask, inverse=T)
#UPF2010<-mask(UPF2010, pgm.degrad.2010.mask, inverse=T)
###cheking
##unique(UPF2010[])
##plot(UPF2010)
#
##saving
#writeRaster(UPF2010, "rasters/PGM/input/UPF2010_real.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (60m)
#UPF2010.px <- focal(UPF2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##UPF2010.px
##anyNA(UPF2010.px[])
##plot(UPF2010.px)
#
#names(UPF2010.px)<-"UPFpx"
#UPF2010.px[is.nan(UPF2010.px)] <- 0
##UPF2010.px <- mask(UPF2010.px, pgm.shp)
#
##saving
#writeRaster(UPF2010.px, "rasters/PGM/2010_real/UPFpx.tif", format="GTiff", overwrite=T)
#
## mean upf cover in landscape scale (1000m)
#UPF2010.ls <- focal(UPF2010, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
###cheking
##UPF2010.ls
##anyNA(UPF2010.ls[])
##plot(UPF2010.ls)
#
#names(UPF2010.ls)<-"UPFls"
#UPF2010.ls[is.nan(UPF2010.ls)] <- 0
##UPF2010.ls <- mask(UPF2010.ls, pgm.shp)
#
##saving
#writeRaster(UPF2010.ls, "rasters/PGM/2010_real/UPFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("UPF2010.px", "UPF2010.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [F3] Forest type 3 or Total forest
## this variable includes forest pixels in LULC raster
## including degraded forest and secondary forest older than 2 years
#
#SF2010.young <- SFage2010
#SF2010.young[SF2010.young <= 2] <- 0
#SF2010.young[SF2010.young > 2] <- 1
#
##saving
#writeRaster(SF2010.young, "rasters/PGM/input/SF2010_real-young.tif", format="GTiff", overwrite=T)
#
#
#TF2010 <- sum(UPF2010, DPF2010, SF2010.young, na.rm = T)
#TF2010[TF2010>1] <- 1
###cheking
##TF2010
##plot(TF2010)
#
##saving
#writeRaster(TF2010, "rasters/PGM/input/TF2010_real.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#TF2010.px <- focal(TF2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TF2010.px
##anyNA(TF2010.px[])
##plot(TF2010.px)
#
#names(TF2010.px)<-"TFpx"
#TF2010.px[is.nan(TF2010.px)] <- 0
#TF2010.px <- mask(TF2010.px, pgm.shp)
#
##saving
#writeRaster(TF2010.px, "rasters/PGM/2010_real/TFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TF2010.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [F1] Forest type 1 or Mature forest
## this variable includes forest pixels in LULC raster
## including degraded forest and secondary forest older than 10 years
#
#SF2010.mature <- SFage2010
#SF2010.mature[SF2010.mature < 10] <- 0
#SF2010.mature[SF2010.mature >= 10] <- 1
#
##saving
#writeRaster(SF2010.mature, "rasters/PGM/input/SFAge2010_real-mature.tif", format="GTiff", overwrite=T)
#
#
#MF2010 <- sum(UPF2010, DPF2010, SF2010.mature, na.rm = T)
#MF2010[MF2010>1] <- 1
###cheking
##MF2010
##plot(MF2010)
#
##saving
#writeRaster(MF2010, "rasters/PGM/input/MF2010_real.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#MF2010.px <- focal(MF2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##MF2010.px
##anyNA(MF2010.px[])
##plot(MF2010.px)
#
#names(MF2010.px)<-"MFpx"
#MF2010.px[is.nan(MF2010.px)] <- 0
#MF2010.px <- mask(MF2010.px, pgm.shp)
#
##saving
#writeRaster(MF2010.px, "rasters/PGM/2010_real/MFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("MF2010.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [edgedist] distance to forest edge
## this variable is the euclidean distance between mature forest 
## to the nearest cell that is not MF
#
#inv.MF2010 <- MF2010
#inv.MF2010[inv.MF2010==1]<-NA
##cheking
##inv.MF2010
##plot(inv.MF2010)
#
#edge.dist.2010 <- distance(inv.MF2010)
###cheking
##edge.dist.2010
##anyNA(edge.dist.2010[])
##plot(edge.dist.2010)
#
#names(edge.dist.2010)<-"edgedist"
#edge.dist.2010[is.nan(edge.dist.2010)] <- 0
#edge.dist.2010 <- mask(edge.dist.2010, pgm.shp)
#
##saving
#writeRaster(edge.dist.2010, "rasters/PGM/2010_real/edgedist.tif", format="GTiff", overwrite=T)
##writeRaster(edge.dist.2010, "rasters/PGM/2020_avoidboth/edgedist.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("inv.MF2010")])
#gc()
#
#
#
##
#
#
## [edge] forest edge
## this variable is the mean of edge area based on mature forest
#
## marking edge and core areas
#edge2010 <- edge.dist.2010
#edge2010[edge2010>300] <- 0
#edge2010[edge2010!=0] <- 1
#
##saving
#writeRaster(edge2010, "rasters/PGM/input/edge2010_real.tif", format="GTiff", overwrite=T)
#
##alternative method
##edge2010 <- focal(MF2010, matrix(1,ncol=3,nrow=3), fun = mean, na.rm = T, pad = T)
##edge2010 <- edge2010 + MF2010
##edge2010[edge2010 == 2] <- 0                  # core area
##edge2010[edge2010 > 0 & edge2010 < 2] <- 1    # edge
#
##cheking
##edge2010
##unique(edge2010[])
##plot(edge2010)
#
## mean edge in pixel scale (200m)
#edge2010.px <- focal(edge2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##edge2010.px
##anyNA(edge2010.px[])
##plot(edge2010.px)
#
#names(edge2010.px)<-"edgepx"
#edge2010.px[is.nan(edge2010.px)] <- 0
#edge2010.px <- mask(edge2010.px, pgm.shp)
#
##saving
#writeRaster(edge2010.px, "rasters/PGM/2010_real/edgepx.tif", format="GTiff", overwrite=T)
##
#
## mean edge in landscape scale (1000m)
#edge2010.ls <- focal(edge2010, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##edge2010.ls
##anyNA(edge2010.ls[])
##plot(edge2010.ls)
#
#names(edge2010.ls)<-"edgels"
#edge2010.ls[is.nan(edge2010.ls)] <- 0
#edge2010.ls <- mask(edge2010.ls, pgm.shp)
#
##saving
#writeRaster(edge2010.px, "rasters/PGM/2010_real/edgels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("edge2010.px", "edge2010.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
#################################################################################################################################.
#####         2020 Real             ####.
########################################.
#
## [SFage] secondary forest age -- pixel: 3x3 (60m); and landscape: 35x35 (1000m)
## this variable is the mean age of secondary forest
#
#SFage2020 <- pgm.sfage[["pgm.sfage.2020real"]]
##plot(SFage2020)
#
##saving
#writeRaster(SFage2020, "rasters/PGM/input/SFage2020_real.tif", format="GTiff", overwrite=T)
#
#
## mean sf age in pixel scale (60m)
#SFage2020.px <- focal(SFage2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##SFage2020.px
##anyNA(SFage2020.px[])
##plot(SFage2020.px)
#
#names(SFage2020.px)<-"SFagepx"
#SFage2020.px[is.nan(SFage2020.px)] <- 0
##SFage2020.px <- mask(SFage2020.px, pgm.shp)
#
##saving
#writeRaster(SFage2020.px, "rasters/PGM/2020_real/SFagepx.tif", format="GTiff", overwrite=T)
#
## mean sf age in landscape scale (1000m)
#SFage2020.ls <- focal(SFage2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
###cheking
##SFage2020.ls
##anyNA(SFage2020.ls[])
##plot(SFage2020.ls)
#
#names(SFage2020.ls)<-"SFagels"
#SFage2020.ls[is.nan(SFage2020.ls)] <- 0
##SFage2020.ls <- mask(SFage2020.ls, pgm.shp)
#
##saving
#writeRaster(SFage2020.ls, "rasters/PGM/2020_real/SFagels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("SFage2020.px", "SFage2020.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [SF] secondary forest
## this variable includes forest pixels in SFage raster
## which has less than 25 years for 2020 or less than 35 for 2020
#
#SF2020 <- pgm.sfage.2020.all.class
##plot(SF2020)
#
##saving
#writeRaster(SF2020, "rasters/PGM/input/SF2020_real.tif", format="GTiff", overwrite=T)
#
#
## mean sf cover in pixel scale (60m)
#SF2020.px <- focal(SF2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##SF2020.px
##anyNA(SF2020.px[])
##plot(SF2020.px)
#
#names(SF2020.px)<-"SFpx"
#SF2020.px[is.nan(SF2020.px)] <- 0
##SF2020.px <- mask(SF2020.px, pgm.shp)
#
##saving
#writeRaster(SF2020.px, "rasters/PGM/2020_real/SFpx.tif", format="GTiff", overwrite=T)
#
## mean sf cover in landscape scale (1000m)
#SF2020.ls <- focal(SF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
###cheking
##SF2020.ls
##anyNA(SF2020.ls[])
##plot(SF2020.ls)
#
#names(SF2020.ls)<-"SFls"
#SF2020.ls[is.nan(SF2020.ls)] <- 0
##SF2020.ls <- mask(SF2020.ls, pgm.shp)
#
##saving
#writeRaster(SF2020.ls, "rasters/PGM/2020_real/SFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("SF2020.px", "SF2020.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [TSD] time since degradation
## this variable is the mean time since a degradation event
#
#TSD2020 <- pgm.degrad[["pgm.degrad.2020real"]]
##plot(TSD2020)
#
##saving
#writeRaster(TSD2020, "rasters/PGM/input/TSD2020_real.tif", format="GTiff", overwrite=T)
#
#
## mean tsd cover in pixel scale (60m)
#TSD2020.px <- focal(TSD2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TSD2020.px
##anyNA(TSD2020.px[])
##plot(TSD2020.px)
#
#names(TSD2020.px)<-"TSDpx"
#TSD2020.px[is.nan(TSD2020.px)] <- 0
##TSD2020.px <- mask(TSD2020.px, pgm.shp)
#
##saving
#writeRaster(TSD2020.px, "rasters/PGM/2020_real/TSDpx.tif", format="GTiff", overwrite=T)
#
## mean tsd cover in landscape scale (1000m)
#TSD2020.ls <- focal(TSD2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
###cheking
##TSD2020.ls
##anyNA(TSD2020.ls[])
##plot(TSD2020.ls)
#
#names(TSD2020.ls)<-"TSDls"
#TSD2020.ls[is.nan(TSD2020.ls)] <- 0
##TSD2020.ls <- mask(TSD2020.ls, pgm.shp)
#
##saving
#writeRaster(TSD2020.ls, "rasters/PGM/2020_real/TSDls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TSD2020.px", "TSD2020.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [DPF] degraded primary forest
## this variable includes forest pixels in TSD raster
## which has less than 24 years for 2020 or less than 34 for 2020
#
#DPF2020 <- pgm.degrad.2020.forest.class
##plot(DPF2020)
#
##saving
#writeRaster(DPF2020, "rasters/PGM/input/DPF2020_real.tif", format="GTiff", overwrite=T)
#
#
## mean dpf cover in pixel scale (60m)
#DPF2020.px <- focal(DPF2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##DPF2020.px
##anyNA(DPF2020.px[])
##plot(DPF2020.px)
#
#names(DPF2020.px)<-"DPFpx"
#DPF2020.px[is.nan(DPF2020.px)] <- 0
##DPF2020.px <- mask(DPF2020.px, pgm.shp)
#
##saving
#writeRaster(DPF2020.px, "rasters/PGM/2020_real/DPFpx.tif", format="GTiff", overwrite=T)
#
## mean dpf cover in landscape scale (1000m)
#DPF2020.ls <- focal(DPF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
###cheking
##DPF2020.ls
##anyNA(DPF2020.ls[])
##plot(DPF2020.ls)
#
#names(DPF2020.ls)<-"DPFls"
#DPF2020.ls[is.nan(DPF2020.ls)] <- 0
##DPF2020.ls <- mask(DPF2020.ls, pgm.shp)
#
##saving
#writeRaster(DPF2020.ls, "rasters/PGM/2020_real/DPFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("DPF2020.px", "DPF2020.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [UPF] undisturbed primary forest
## this variable includes all forest pixels in LULC raster (value == 3)
## excluding degraded and secondary forest
#
#UPF2020<-pgm.lulc.2020.forest.class
#UPF2020<-mask(UPF2020, pgm.sfage.2020.mask, inverse=T)
#UPF2020<-mask(UPF2020, pgm.degrad.2020.mask, inverse=T)
###cheking
##unique(UPF2020[])
##plot(UPF2020)
#
##saving
#writeRaster(UPF2020, "rasters/PGM/input/UPF2020_real.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (60m)
#UPF2020.px <- focal(UPF2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##UPF2020.px
##anyNA(UPF2020.px[])
##plot(UPF2020.px)
#
#names(UPF2020.px)<-"UPFpx"
#UPF2020.px[is.nan(UPF2020.px)] <- 0
##UPF2020.px <- mask(UPF2020.px, pgm.shp)
#
##saving
#writeRaster(UPF2020.px, "rasters/PGM/2020_real/UPFpx.tif", format="GTiff", overwrite=T)
#
## mean upf cover in landscape scale (1000m)
#UPF2020.ls <- focal(UPF2020, matrix(1,ncol=35,nrow=35), fun=mean, na.rm=T)
###cheking
##UPF2020.ls
##anyNA(UPF2020.ls[])
##plot(UPF2020.ls)
#
#names(UPF2020.ls)<-"UPFls"
#UPF2020.ls[is.nan(UPF2020.ls)] <- 0
##UPF2020.ls <- mask(UPF2020.ls, pgm.shp)
#
##saving
#writeRaster(UPF2020.ls, "rasters/PGM/2020_real/UPFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("UPF2020.px", "UPF2020.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [F3] Forest type 3 or Total forest
## this variable includes forest pixels in LULC raster
## including degraded forest and secondary forest older than 2 years
#
#SF2020.young <- SFage2020
#SF2020.young[SF2020.young <= 2] <- 0
#SF2020.young[SF2020.young > 2] <- 1
#
##saving
#writeRaster(SF2020.young, "rasters/PGM/input/SF2020_real-young.tif", format="GTiff", overwrite=T)
#
#
#TF2020 <- sum(UPF2020, DPF2020, SF2020.young, na.rm = T)
#TF2020[TF2020>1] <- 1
###cheking
##TF2020
##plot(TF2020)
#
##saving
#writeRaster(TF2020, "rasters/PGM/input/TF2020_real.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#TF2020.px <- focal(TF2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TF2020.px
##anyNA(TF2020.px[])
##plot(TF2020.px)
#
#names(TF2020.px)<-"TFpx"
#TF2020.px[is.nan(TF2020.px)] <- 0
#TF2020.px <- mask(TF2020.px, pgm.shp)
#
##saving
#writeRaster(TF2020.px, "rasters/PGM/2020_real/TFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TF2020.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [F1] Forest type 1 or Mature forest
## this variable includes forest pixels in LULC raster
## including degraded forest and secondary forest older than 10 years
#
#SF2020.mature <- SFage2020
#SF2020.mature[SF2020.mature < 10] <- 0
#SF2020.mature[SF2020.mature >= 10] <- 1
#
##saving
#writeRaster(SF2020.mature, "rasters/PGM/input/SFAge2020_real-mature.tif", format="GTiff", overwrite=T)
#
#
#MF2020 <- sum(UPF2020, DPF2020, SF2020.mature, na.rm = T)
#MF2020[MF2020>1] <- 1
###cheking
##MF2020
##plot(MF2020)
#
##saving
#writeRaster(MF2020, "rasters/PGM/input/MF2020_real.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#MF2020.px <- focal(MF2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##MF2020.px
##anyNA(MF2020.px[])
##plot(MF2020.px)
#
#names(MF2020.px)<-"MFpx"
#MF2020.px[is.nan(MF2020.px)] <- 0
#MF2020.px <- mask(MF2020.px, pgm.shp)
#
##saving
#writeRaster(MF2020.px, "rasters/PGM/2020_real/MFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("MF2020.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [edgedist] distance to forest edge
## this variable is the euclidean distance between mature forest 
## to the nearest cell that is not MF
#
#inv.MF2020 <- MF2020
#inv.MF2020[inv.MF2020==1]<-NA
##cheking
##inv.MF2020
##plot(inv.MF2020)
#
#edge.dist.2020 <- distance(inv.MF2020)
###cheking
##edge.dist.2020
##anyNA(edge.dist.2020[])
##plot(edge.dist.2020)
#
#names(edge.dist.2020)<-"edgedist"
#edge.dist.2020[is.nan(edge.dist.2020)] <- 0
#edge.dist.2020 <- mask(edge.dist.2020, pgm.shp)
#
##saving
#writeRaster(edge.dist.2020, "rasters/PGM/2020_real/edgedist.tif", format="GTiff", overwrite=T)
##writeRaster(edge.dist.2020, "rasters/PGM/2020_avoidboth/edgedist.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("inv.MF2020")])
#gc()
#
#
#
##
#
#
## [edge] forest edge
## this variable is the mean of edge area based on mature forest
#
## marking edge and core areas
#edge2020 <- edge.dist.2020
#edge2020[edge2020>300] <- 0
#edge2020[edge2020!=0] <- 1
#
##saving
#writeRaster(edge2020, "rasters/PGM/input/edge2020_real.tif", format="GTiff", overwrite=T)
#
##alternative method
##edge2020 <- focal(MF2020, matrix(1,ncol=3,nrow=3), fun = mean, na.rm = T, pad = T)
##edge2020 <- edge2020 + MF2020
##edge2020[edge2020 == 2] <- 0                  # core area
##edge2020[edge2020 > 0 & edge2020 < 2] <- 1    # edge
#
##cheking
##edge2020
##unique(edge2020[])
##plot(edge2020)
#
## mean edge in pixel scale (200m)
#edge2020.px <- focal(edge2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##edge2020.px
##anyNA(edge2020.px[])
##plot(edge2020.px)
#
#names(edge2020.px)<-"edgepx"
#edge2020.px[is.nan(edge2020.px)] <- 0
#edge2020.px <- mask(edge2020.px, pgm.shp)
#
##saving
#writeRaster(edge2020.px, "rasters/PGM/2020_real/edgepx.tif", format="GTiff", overwrite=T)
##
#
## mean edge in landscape scale (1000m)
#edge2020.ls <- focal(edge2020, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##edge2020.ls
##anyNA(edge2020.ls[])
##plot(edge2020.ls)
#
#names(edge2020.ls)<-"edgels"
#edge2020.ls[is.nan(edge2020.ls)] <- 0
#edge2020.ls <- mask(edge2020.ls, pgm.shp)
#
##saving
#writeRaster(edge2020.px, "rasters/PGM/2020_real/edgels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("edge2020.px", "edge2020.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
#################################################################################################################################.
#####         2020 Avoid degradation             ####.
#
## [SFage] secondary forest age -- pixel: 3x3 (200m); and landscape: 11x11 (1000m)
## this variable is the mean age of secondary forest
#
#SFage2020_avoiddegrad <- pgm.sfage[["pgm.sfage.2020real"]]
##plot(SFage2020_avoiddegrad)
#
##saving
##the same as 2020 Real
#
## mean sf age in pixel scale (200m)
#SFage2020_avoiddegrad.px <- focal(SFage2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##SFage2020_avoiddegrad.px
##anyNA(SFage2020_avoiddegrad.px[])
##plot(SFage2020_avoiddegrad.px)
#
#names(SFage2020_avoiddegrad.px)<-"SFagepx"
#SFage2020_avoiddegrad.px[is.nan(SFage2020_avoiddegrad.px)] <- 0
#SFage2020_avoiddegrad.px <- mask(SFage2020_avoiddegrad.px, pgm.shp)
#
##saving
#writeRaster(SFage2020_avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/SFagepx.tif", format="GTiff", overwrite=T)
#
## mean sf age in landscape scale (1000m)
#SFage2020_avoiddegrad.ls <- focal(SFage2020_avoiddegrad, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##SFage2020_avoiddegrad.ls
##anyNA(SFage2020_avoiddegrad.ls[])
##plot(SFage2020_avoiddegrad.ls)
#
#names(SFage2020_avoiddegrad.ls)<-"SFagels"
#SFage2020_avoiddegrad.ls[is.nan(SFage2020_avoiddegrad.ls)] <- 0
#SFage2020_avoiddegrad.ls <- mask(SFage2020_avoiddegrad.ls, pgm.shp)
#
##saving
#writeRaster(SFage2020_avoiddegrad.ls, "rasters/PGM/2020_avoiddegrad/SFagels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("SFage2020_avoiddegrad.px", "SFage2020_avoiddegrad.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [SF] secondary forest
## this variable includes forest pixels in SFage raster
## which has less than 25 years for 2020_avoiddegrad or less than 35 for 2020_avoiddegrad
#
#SF2020_avoiddegrad <- pgm.sfage[["pgm.sfage.2020real"]]
#SF2020_avoiddegrad[SF2020_avoiddegrad>34] <- 0
#SF2020_avoiddegrad[SF2020_avoiddegrad>0] <- 1
#SF2020_avoiddegrad[SF2020_avoiddegrad<1] <- 0
#
##plot(SF2020_avoiddegrad)
#
##saving
##the same as 2020 Real
#
#
## mean sf cover in pixel scale (200m)
#SF2020_avoiddegrad.px <- focal(SF2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##SF2020_avoiddegrad.px
##anyNA(SF2020_avoiddegrad.px[])
##plot(SF2020_avoiddegrad.px)
#
#names(SF2020_avoiddegrad.px)<-"SFpx"
#SF2020_avoiddegrad.px[is.nan(SF2020_avoiddegrad.px)] <- 0
#SF2020_avoiddegrad.px <- mask(SF2020_avoiddegrad.px, pgm.shp)
#
##saving
#writeRaster(SF2020_avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/SFpx.tif", format="GTiff", overwrite=T)
#
## mean sf cover in landscape scale (1000m)
#SF2020_avoiddegrad.ls <- focal(SF2020_avoiddegrad, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##SF2020_avoiddegrad.ls
##anyNA(SF2020_avoiddegrad.ls[])
##plot(SF2020_avoiddegrad.ls)
#
#names(SF2020_avoiddegrad.ls)<-"SFls"
#SF2020_avoiddegrad.ls[is.nan(SF2020_avoiddegrad.ls)] <- 0
#SF2020_avoiddegrad.ls <- mask(SF2020_avoiddegrad.ls, pgm.shp)
#
##saving
#writeRaster(SF2020_avoiddegrad.ls, "rasters/PGM/2020_avoiddegrad/SFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("SF2020_avoiddegrad.px", "SF2020_avoiddegrad.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [TSD] time since degradation
## this variable is the mean time since a degradation event
#
#TSD2020_avoiddegrad <- calc(TSD2010, fun=function(x){ifelse(is.na(x), x, ifelse(x==0, x, ifelse(x==300, x, x+10)))})
#TSD2020_avoiddegrad[get("SF2020_avoiddegrad")[] == 1] <- 0
##plot(TSD2020_avoiddegrad)
#
##saving
#writeRaster(TSD2020_avoiddegrad, "rasters/PGM/input/TSD2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#
#
## mean tsd cover in pixel scale (200m)
#TSD2020_avoiddegrad.px <- focal(TSD2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TSD2020_avoiddegrad.px
##anyNA(TSD2020_avoiddegrad.px[])
##plot(TSD2020_avoiddegrad.px)
#
#names(TSD2020_avoiddegrad.px)<-"TSDpx"
#TSD2020_avoiddegrad.px[is.nan(TSD2020_avoiddegrad.px)] <- 0
#TSD2020_avoiddegrad.px <- mask(TSD2020_avoiddegrad.px, pgm.shp)
#
##saving
#writeRaster(TSD2020_avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/TSDpx.tif", format="GTiff", overwrite=T)
#
## mean tsd cover in landscape scale (1000m)
#TSD2020_avoiddegrad.ls <- focal(TSD2020_avoiddegrad, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##TSD2020_avoiddegrad.ls
##anyNA(TSD2020_avoiddegrad.ls[])
##plot(TSD2020_avoiddegrad.ls)
#
#names(TSD2020_avoiddegrad.ls)<-"TSDls"
#TSD2020_avoiddegrad.ls[is.nan(TSD2020_avoiddegrad.ls)] <- 0
#TSD2020_avoiddegrad.ls <- mask(TSD2020_avoiddegrad.ls, pgm.shp)
#
##saving
#writeRaster(TSD2020_avoiddegrad.ls, "rasters/PGM/2020_avoiddegrad/TSDls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TSD2020_avoiddegrad.px", "TSD2020_avoiddegrad.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [DPF] degraded primary forest
## this variable includes forest pixels in TSD raster
## which has less than 24 years for 2020_avoiddegrad or less than 34 for 2020_avoiddegrad
#
#DPF2020_avoiddegrad <- TSD2020_avoiddegrad
#DPF2020_avoiddegrad[DPF2020_avoiddegrad>33]<-NA
#DPF2020_avoiddegrad[DPF2020_avoiddegrad==0]<-NA
#DPF2020_avoiddegrad[!is.na(DPF2020_avoiddegrad)] <- 1
#DPF2020_avoiddegrad<-sum(pgm.lulc.2020.forest.class, DPF2020_avoiddegrad, na.rm=T)
#DPF2020_avoiddegrad[DPF2020_avoiddegrad!=2]<-0
#DPF2020_avoiddegrad[DPF2020_avoiddegrad==2]<-1
##plot(DPF2020_avoiddegrad)
#
##saving
#writeRaster(DPF2020_avoiddegrad, "rasters/PGM/input/DPF2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#
#
## mean dpf cover in pixel scale (200m)
#DPF2020_avoiddegrad.px <- focal(DPF2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##DPF2020_avoiddegrad.px
##anyNA(DPF2020_avoiddegrad.px[])
##plot(DPF2020_avoiddegrad.px)
#
#names(DPF2020_avoiddegrad.px)<-"DPFpx"
#DPF2020_avoiddegrad.px[is.nan(DPF2020_avoiddegrad.px)] <- 0
#DPF2020_avoiddegrad.px <- mask(DPF2020_avoiddegrad.px, pgm.shp)
#
##saving
#writeRaster(DPF2020_avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/DPFpx.tif", format="GTiff", overwrite=T)
#
## mean dpf cover in landscape scale (1000m)
#DPF2020_avoiddegrad.ls <- focal(DPF2020_avoiddegrad, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##DPF2020_avoiddegrad.ls
##anyNA(DPF2020_avoiddegrad.ls[])
##plot(DPF2020_avoiddegrad.ls)
#
#names(DPF2020_avoiddegrad.ls)<-"DPFls"
#DPF2020_avoiddegrad.ls[is.nan(DPF2020_avoiddegrad.ls)] <- 0
#DPF2020_avoiddegrad.ls <- mask(DPF2020_avoiddegrad.ls, pgm.shp)
#
##saving
#writeRaster(DPF2020_avoiddegrad.ls, "rasters/PGM/2020_avoiddegrad/DPFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("pgm.degrad.2020.forest.class", "DPF2020_avoiddegrad.px", "DPF2020_avoiddegrad.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [UPF] undisturbed primary forest
## this variable includes all forest pixels in LULC raster (value == 3)
## excluding degraded and secondary forest
#
#UPF2020_avoiddegrad<-sum(pgm.lulc.2020.forest.class, SF2020_avoiddegrad, DPF2020_avoiddegrad, na.rm = T)
#UPF2020_avoiddegrad[UPF2020_avoiddegrad>1]<-0
###cheking
##unique(UPF2020_avoiddegrad[])
##plot(UPF2020_avoiddegrad)
#
##saving
#writeRaster(UPF2020_avoiddegrad, "rasters/PGM/input/UPF2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#UPF2020_avoiddegrad.px <- focal(UPF2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##UPF2020_avoiddegrad.px
##anyNA(UPF2020_avoiddegrad.px[])
##plot(UPF2020_avoiddegrad.px)
#
#names(UPF2020_avoiddegrad.px)<-"UPFpx"
#UPF2020_avoiddegrad.px[is.nan(UPF2020_avoiddegrad.px)] <- 0
#UPF2020_avoiddegrad.px <- mask(UPF2020_avoiddegrad.px, pgm.shp)
#
##saving
#writeRaster(UPF2020_avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/UPFpx.tif", format="GTiff", overwrite=T)
#
## mean upf cover in landscape scale (1000m)
#UPF2020_avoiddegrad.ls <- focal(UPF2020_avoiddegrad, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##UPF2020_avoiddegrad.ls
##anyNA(UPF2020_avoiddegrad.ls[])
##plot(UPF2020_avoiddegrad.ls)
#
#names(UPF2020_avoiddegrad.ls)<-"UPFls"
#UPF2020_avoiddegrad.ls[is.nan(UPF2020_avoiddegrad.ls)] <- 0
#UPF2020_avoiddegrad.ls <- mask(UPF2020_avoiddegrad.ls, pgm.shp)
#
##saving
#writeRaster(UPF2020_avoiddegrad.ls, "rasters/PGM/2020_avoiddegrad/UPFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("UPF2020_avoiddegrad.px", "UPF2020_avoiddegrad.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [F3] Forest type 3 or Total forest
## this variable includes forest pixels in LULC raster
## including degraded forest and secondary forest older than 2 years
#
#SF2020_avoiddegrad.young <- SFage2020_avoiddegrad
#SF2020_avoiddegrad.young[SF2020_avoiddegrad.young <= 2] <- 0
#SF2020_avoiddegrad.young[SF2020_avoiddegrad.young > 2] <- 1
#
##saving
##same as SF2020.young
#
#
#TF2020_avoiddegrad <- sum(UPF2020_avoiddegrad, DPF2020_avoiddegrad, SF2020_avoiddegrad.young, na.rm = T)
#TF2020_avoiddegrad[TF2020_avoiddegrad>1] <- 1
###cheking
##TF2020_avoiddegrad
##plot(TF2020_avoiddegrad)
#
##saving
#writeRaster(TF2020_avoiddegrad, "rasters/PGM/input/TF2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#TF2020_avoiddegrad.px <- focal(TF2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TF2020_avoiddegrad.px
##anyNA(TF2020_avoiddegrad.px[])
##plot(TF2020_avoiddegrad.px)
#
#names(TF2020_avoiddegrad.px)<-"TFpx"
#TF2020_avoiddegrad.px[is.nan(TF2020_avoiddegrad.px)] <- 0
#TF2020_avoiddegrad.px <- mask(TF2020_avoiddegrad.px, pgm.shp)
#
##saving
#writeRaster(TF2020_avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/TFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TF2020_avoiddegrad.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [F1] Forest type 1 or Mature forest
## this variable includes forest pixels in LULC raster
## including degraded forest and secondary forest older than 10 years
#
#SF2020_avoiddegrad.mature <- SFage2020_avoiddegrad
#SF2020_avoiddegrad.mature[SF2020_avoiddegrad.mature < 10] <- 0
#SF2020_avoiddegrad.mature[SF2020_avoiddegrad.mature >= 10] <- 1
#
##saving
##same as SF2020.mature
#
#
#MF2020_avoiddegrad <- sum(UPF2020_avoiddegrad, DPF2020_avoiddegrad, SF2020_avoiddegrad.mature, na.rm = T)
#MF2020_avoiddegrad[MF2020_avoiddegrad>1] <- 1
###cheking
##MF2020_avoiddegrad
##plot(MF2020_avoiddegrad)
#
##saving
#writeRaster(MF2020_avoiddegrad, "rasters/PGM/input/MF2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#MF2020_avoiddegrad.px <- focal(MF2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##MF2020_avoiddegrad.px
##anyNA(MF2020_avoiddegrad.px[])
##plot(MF2020_avoiddegrad.px)
#
#names(MF2020_avoiddegrad.px)<-"MFpx"
#MF2020_avoiddegrad.px[is.nan(MF2020_avoiddegrad.px)] <- 0
#MF2020_avoiddegrad.px <- mask(MF2020_avoiddegrad.px, pgm.shp)
#
##saving
#writeRaster(MF2020_avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/MFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("MF2020_avoiddegrad.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [edgedist] distance to forest edge
## this variable is the euclidean distance between mature forest 
## to the nearest cell that is not MF
#
#inv.MF2020_avoiddegrad <- MF2020_avoiddegrad
#inv.MF2020_avoiddegrad[inv.MF2020_avoiddegrad==1]<-NA
##cheking
##inv.MF2020_avoiddegrad
##plot(inv.MF2020_avoiddegrad)
#
#edge.dist.2020_avoiddegrad <- distance(inv.MF2020_avoiddegrad)
###cheking
##edge.dist.2020_avoiddegrad
##anyNA(edge.dist.2020_avoiddegrad[])
##plot(edge.dist.2020_avoiddegrad)
#
#names(edge.dist.2020_avoiddegrad)<-"edgedist"
#edge.dist.2020_avoiddegrad[is.nan(edge.dist.2020_avoiddegrad)] <- 0
#edge.dist.2020_avoiddegrad <- mask(edge.dist.2020_avoiddegrad, pgm.shp)
#
##saving
#writeRaster(edge.dist.2020_avoiddegrad, "rasters/PGM/2020_avoiddegrad/edgedist.tif", format="GTiff", overwrite=T)
##writeRaster(edge.dist.2020_avoiddegrad, "rasters/PGM/2020_avoiddegrad_avoidboth/edgedist.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("inv.MF2020_avoiddegrad")])
#gc()
#
#
#
##
#
#
## [edge] forest edge
## this variable is the mean of edge area based on mature forest
#
## marking edge and core areas
#edge2020_avoiddegrad <- edge.dist.2020_avoiddegrad
#edge2020_avoiddegrad[edge2020_avoiddegrad>300] <- 0
#edge2020_avoiddegrad[edge2020_avoiddegrad!=0] <- 1
#
##saving
#writeRaster(edge2020_avoiddegrad, "rasters/PGM/input/edge2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#
##alternative method
##edge2020_avoiddegrad <- focal(MF2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun = mean, na.rm = T, pad = T)
##edge2020_avoiddegrad <- edge2020_avoiddegrad + MF2020_avoiddegrad
##edge2020_avoiddegrad[edge2020_avoiddegrad == 2] <- 0                  # core area
##edge2020_avoiddegrad[edge2020_avoiddegrad > 0 & edge2020_avoiddegrad < 2] <- 1    # edge
#
##cheking
##edge2020_avoiddegrad
##unique(edge2020_avoiddegrad[])
##plot(edge2020_avoiddegrad)
#
## mean edge in pixel scale (200m)
#edge2020_avoiddegrad.px <- focal(edge2020_avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##edge2020_avoiddegrad.px
##anyNA(edge2020_avoiddegrad.px[])
##plot(edge2020_avoiddegrad.px)
#
#names(edge2020_avoiddegrad.px)<-"edgepx"
#edge2020_avoiddegrad.px[is.nan(edge2020_avoiddegrad.px)] <- 0
#edge2020_avoiddegrad.px <- mask(edge2020_avoiddegrad.px, pgm.shp)
#
##saving
#writeRaster(edge2020_avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/edgepx.tif", format="GTiff", overwrite=T)
##
#
## mean edge in landscape scale (1000m)
#edge2020_avoiddegrad.ls <- focal(edge2020_avoiddegrad, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##edge2020_avoiddegrad.ls
##anyNA(edge2020_avoiddegrad.ls[])
##plot(edge2020_avoiddegrad.ls)
#
#names(edge2020_avoiddegrad.ls)<-"edgels"
#edge2020_avoiddegrad.ls[is.nan(edge2020_avoiddegrad.ls)] <- 0
#edge2020_avoiddegrad.ls <- mask(edge2020_avoiddegrad.ls, pgm.shp)
#
##saving
#writeRaster(edge2020_avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/edgels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("edge2020_avoiddegrad.px", "edge2020_avoiddegrad.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
#################################################################################################################################.
#####         2020 Avoid deforestation             ####.
#
## [SFage] secondary forest age -- pixel: 3x3 (200m); and landscape: 11x11 (1000m)
## this variable is the mean age of secondary forest
#
#SFage2020_avoiddeforest <- calc(SFage2010, fun=function(x){ifelse(is.na(x), x, ifelse(x==0, x, x+10))})
##plot(SFage2020_avoiddeforest)
#
##saving
#writeRaster(SFage2020_avoiddeforest, "rasters/PGM/input/SFage2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#
## mean sf age in pixel scale (200m)
#SFage2020_avoiddeforest.px <- focal(SFage2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##SFage2020_avoiddeforest.px
##anyNA(SFage2020_avoiddeforest.px[])
##plot(SFage2020_avoiddeforest.px)
#
#names(SFage2020_avoiddeforest.px)<-"SFagepx"
#SFage2020_avoiddeforest.px[is.nan(SFage2020_avoiddeforest.px)] <- 0
#SFage2020_avoiddeforest.px <- mask(SFage2020_avoiddeforest.px, pgm.shp)
#
##saving
#writeRaster(SFage2020_avoiddeforest.px, "rasters/PGM/2020_avoiddeforest/SFagepx.tif", format="GTiff", overwrite=T)
#
## mean sf age in landscape scale (1000m)
#SFage2020_avoiddeforest.ls <- focal(SFage2020_avoiddeforest, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##SFage2020_avoiddeforest.ls
##anyNA(SFage2020_avoiddeforest.ls[])
##plot(SFage2020_avoiddeforest.ls)
#
#names(SFage2020_avoiddeforest.ls)<-"SFagels"
#SFage2020_avoiddeforest.ls[is.nan(SFage2020_avoiddeforest.ls)] <- 0
#SFage2020_avoiddeforest.ls <- mask(SFage2020_avoiddeforest.ls, pgm.shp)
#
##saving
#writeRaster(SFage2020_avoiddeforest.ls, "rasters/PGM/2020_avoiddeforest/SFagels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("SFage2020_avoiddeforest.px", "SFage2020_avoiddeforest.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [SF] secondary forest
## this variable includes forest pixels in SFage raster
## which has less than 25 years for 2020_avoiddeforest or less than 35 for 2020_avoiddeforest
#
#SF2020_avoiddeforest <- SFage2020_avoiddeforest
#SF2020_avoiddeforest[SF2020_avoiddeforest>34] <- 0
#SF2020_avoiddeforest[SF2020_avoiddeforest>0] <- 1
#SF2020_avoiddeforest[SF2020_avoiddeforest<1] <- 0
#
##plot(SF2020_avoiddeforest)
#
##saving
#writeRaster(SF2020_avoiddeforest, "rasters/PGM/input/SF2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#
#
## mean sf cover in pixel scale (200m)
#SF2020_avoiddeforest.px <- focal(SF2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##SF2020_avoiddeforest.px
##anyNA(SF2020_avoiddeforest.px[])
##plot(SF2020_avoiddeforest.px)
#
#names(SF2020_avoiddeforest.px)<-"SFpx"
#SF2020_avoiddeforest.px[is.nan(SF2020_avoiddeforest.px)] <- 0
#SF2020_avoiddeforest.px <- mask(SF2020_avoiddeforest.px, pgm.shp)
#
##saving
#writeRaster(SF2020_avoiddeforest.px, "rasters/PGM/2020_avoiddeforest/SFpx.tif", format="GTiff", overwrite=T)
#
## mean sf cover in landscape scale (1000m)
#SF2020_avoiddeforest.ls <- focal(SF2020_avoiddeforest, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##SF2020_avoiddeforest.ls
##anyNA(SF2020_avoiddeforest.ls[])
##plot(SF2020_avoiddeforest.ls)
#
#names(SF2020_avoiddeforest.ls)<-"SFls"
#SF2020_avoiddeforest.ls[is.nan(SF2020_avoiddeforest.ls)] <- 0
#SF2020_avoiddeforest.ls <- mask(SF2020_avoiddeforest.ls, pgm.shp)
#
##saving
#writeRaster(SF2020_avoiddeforest.ls, "rasters/PGM/2020_avoiddeforest/SFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("SF2020_avoiddeforest.px", "SF2020_avoiddeforest.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [TSD] time since degradation
## this variable is the mean time since a degradation event
#
#TSD2020_avoiddeforest <- pgm.degrad[["pgm.degrad.2020real"]]
#TSD2020_avoiddeforest[get("SF2020_avoiddeforest")[] == 1] <- 0
##plot(TSD2020_avoiddeforest)
#
##saving
#writeRaster(TSD2020_avoiddeforest, "rasters/PGM/input/TSD2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#
#
## mean tsd cover in pixel scale (200m)
#TSD2020_avoiddeforest.px <- focal(TSD2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TSD2020_avoiddeforest.px
##anyNA(TSD2020_avoiddeforest.px[])
##plot(TSD2020_avoiddeforest.px)
#
#names(TSD2020_avoiddeforest.px)<-"TSDpx"
#TSD2020_avoiddeforest.px[is.nan(TSD2020_avoiddeforest.px)] <- 0
#TSD2020_avoiddeforest.px <- mask(TSD2020_avoiddeforest.px, pgm.shp)
#
##saving
#writeRaster(TSD2020_avoiddeforest.px, "rasters/PGM/2020_avoiddeforest/TSDpx.tif", format="GTiff", overwrite=T)
#
## mean tsd cover in landscape scale (1000m)
#TSD2020_avoiddeforest.ls <- focal(TSD2020_avoiddeforest, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##TSD2020_avoiddeforest.ls
##anyNA(TSD2020_avoiddeforest.ls[])
##plot(TSD2020_avoiddeforest.ls)
#
#names(TSD2020_avoiddeforest.ls)<-"TSDls"
#TSD2020_avoiddeforest.ls[is.nan(TSD2020_avoiddeforest.ls)] <- 0
#TSD2020_avoiddeforest.ls <- mask(TSD2020_avoiddeforest.ls, pgm.shp)
#
##saving
#writeRaster(TSD2020_avoiddeforest.ls, "rasters/PGM/2020_avoiddeforest/TSDls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TSD2020_avoiddeforest.px", "TSD2020_avoiddeforest.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [DPF] degraded primary forest
## this variable includes forest pixels in TSD raster
## which has less than 24 years for 2020_avoiddeforest or less than 34 for 2020_avoiddeforest
#
#DPF2020_avoiddeforest <- TSD2020_avoiddeforest
#DPF2020_avoiddeforest[DPF2020_avoiddeforest>33]<-NA
#DPF2020_avoiddeforest[DPF2020_avoiddeforest==0]<-NA
#DPF2020_avoiddeforest[!is.na(DPF2020_avoiddeforest)] <- 1
#DPF2020_avoiddeforest<-sum(pgm.lulc.2020.forest.class, DPF2020_avoiddeforest, na.rm=T)
#DPF2020_avoiddeforest[DPF2020_avoiddeforest!=2]<-0
#DPF2020_avoiddeforest[DPF2020_avoiddeforest==2]<-1
##plot(DPF2020_avoiddeforest)
#
##saving
#writeRaster(DPF2020_avoiddeforest, "rasters/PGM/input/DPF2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#
#
## mean dpf cover in pixel scale (200m)
#DPF2020_avoiddeforest.px <- focal(DPF2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##DPF2020_avoiddeforest.px
##anyNA(DPF2020_avoiddeforest.px[])
##plot(DPF2020_avoiddeforest.px)
#
#names(DPF2020_avoiddeforest.px)<-"DPFpx"
#DPF2020_avoiddeforest.px[is.nan(DPF2020_avoiddeforest.px)] <- 0
#DPF2020_avoiddeforest.px <- mask(DPF2020_avoiddeforest.px, pgm.shp)
#
##saving
#writeRaster(DPF2020_avoiddeforest.px, "rasters/PGM/2020_avoiddeforest/DPFpx.tif", format="GTiff", overwrite=T)
#
## mean dpf cover in landscape scale (1000m)
#DPF2020_avoiddeforest.ls <- focal(DPF2020_avoiddeforest, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##DPF2020_avoiddeforest.ls
##anyNA(DPF2020_avoiddeforest.ls[])
##plot(DPF2020_avoiddeforest.ls)
#
#names(DPF2020_avoiddeforest.ls)<-"DPFls"
#DPF2020_avoiddeforest.ls[is.nan(DPF2020_avoiddeforest.ls)] <- 0
#DPF2020_avoiddeforest.ls <- mask(DPF2020_avoiddeforest.ls, pgm.shp)
#
##saving
#writeRaster(DPF2020_avoiddeforest.ls, "rasters/PGM/2020_avoiddeforest/DPFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("pgm.degrad.2020.forest.class", "DPF2020_avoiddeforest.px", "DPF2020_avoiddeforest.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [UPF] undisturbed primary forest
## this variable includes all forest pixels in LULC raster (value == 3)
## excluding degraded and secondary forest
#
#UPF2020_avoiddeforest<-sum(pgm.lulc.2010.forest.class, SF2020_avoiddeforest, DPF2020_avoiddeforest, na.rm = T)
#UPF2020_avoiddeforest[UPF2020_avoiddeforest>1]<-0
###cheking
##unique(UPF2020_avoiddeforest[])
##plot(UPF2020_avoiddeforest)
#
##saving
#writeRaster(UPF2020_avoiddeforest, "rasters/PGM/input/UPF2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#UPF2020_avoiddeforest.px <- focal(UPF2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##UPF2020_avoiddeforest.px
##anyNA(UPF2020_avoiddeforest.px[])
##plot(UPF2020_avoiddeforest.px)
#
#names(UPF2020_avoiddeforest.px)<-"UPFpx"
#UPF2020_avoiddeforest.px[is.nan(UPF2020_avoiddeforest.px)] <- 0
#UPF2020_avoiddeforest.px <- mask(UPF2020_avoiddeforest.px, pgm.shp)
#
##saving
#writeRaster(UPF2020_avoiddeforest.px, "rasters/PGM/2020_avoiddeforest/UPFpx.tif", format="GTiff", overwrite=T)
#
## mean upf cover in landscape scale (1000m)
#UPF2020_avoiddeforest.ls <- focal(UPF2020_avoiddeforest, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##UPF2020_avoiddeforest.ls
##anyNA(UPF2020_avoiddeforest.ls[])
##plot(UPF2020_avoiddeforest.ls)
#
#names(UPF2020_avoiddeforest.ls)<-"UPFls"
#UPF2020_avoiddeforest.ls[is.nan(UPF2020_avoiddeforest.ls)] <- 0
#UPF2020_avoiddeforest.ls <- mask(UPF2020_avoiddeforest.ls, pgm.shp)
#
##saving
#writeRaster(UPF2020_avoiddeforest.ls, "rasters/PGM/2020_avoiddeforest/UPFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("UPF2020_avoiddeforest.px", "UPF2020_avoiddeforest.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [F3] Forest type 3 or Total forest
## this variable includes forest pixels in LULC raster
## including degraded forest and secondary forest older than 2 years
#
#SF2020_avoiddeforest.young <- SFage2020_avoiddeforest
#SF2020_avoiddeforest.young[SF2020_avoiddeforest.young <= 2] <- 0
#SF2020_avoiddeforest.young[SF2020_avoiddeforest.young > 2] <- 1
#
##saving
#writeRaster(SF2020_avoiddeforest.young, "rasters/PGM/input/SF2020_avoiddeforest.young.tif", format="GTiff", overwrite=T)
#
#
#TF2020_avoiddeforest <- sum(UPF2020_avoiddeforest, DPF2020_avoiddeforest, SF2020_avoiddeforest.young, na.rm = T)
#TF2020_avoiddeforest[TF2020_avoiddeforest>1] <- 1
###cheking
##TF2020_avoiddeforest
##plot(TF2020_avoiddeforest)
#
##saving
#writeRaster(TF2020_avoiddeforest, "rasters/PGM/input/TF2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#TF2020_avoiddeforest.px <- focal(TF2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TF2020_avoiddeforest.px
##anyNA(TF2020_avoiddeforest.px[])
##plot(TF2020_avoiddeforest.px)
#
#names(TF2020_avoiddeforest.px)<-"TFpx"
#TF2020_avoiddeforest.px[is.nan(TF2020_avoiddeforest.px)] <- 0
#TF2020_avoiddeforest.px <- mask(TF2020_avoiddeforest.px, pgm.shp)
#
##saving
#writeRaster(TF2020_avoiddeforest.px, "rasters/PGM/2020_avoiddeforest/TFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TF2020_avoiddeforest.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [F1] Forest type 1 or Mature forest
## this variable includes forest pixels in LULC raster
## including degraded forest and secondary forest older than 10 years
#
#SF2020_avoiddeforest.mature <- SFage2020_avoiddeforest
#SF2020_avoiddeforest.mature[SF2020_avoiddeforest.mature < 10] <- 0
#SF2020_avoiddeforest.mature[SF2020_avoiddeforest.mature >= 10] <- 1
#
##saving
#writeRaster(SF2020_avoiddeforest.mature, "rasters/PGM/input/SF2020_avoiddeforest.mature.tif", format="GTiff", overwrite=T)
#
#
#MF2020_avoiddeforest <- sum(UPF2020_avoiddeforest, DPF2020_avoiddeforest, SF2020_avoiddeforest.mature, na.rm = T)
#MF2020_avoiddeforest[MF2020_avoiddeforest>1] <- 1
###cheking
##MF2020_avoiddeforest
##plot(MF2020_avoiddeforest)
#
##saving
#writeRaster(MF2020_avoiddeforest, "rasters/PGM/input/MF2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#MF2020_avoiddeforest.px <- focal(MF2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##MF2020_avoiddeforest.px
##anyNA(MF2020_avoiddeforest.px[])
##plot(MF2020_avoiddeforest.px)
#
#names(MF2020_avoiddeforest.px)<-"MFpx"
#MF2020_avoiddeforest.px[is.nan(MF2020_avoiddeforest.px)] <- 0
#MF2020_avoiddeforest.px <- mask(MF2020_avoiddeforest.px, pgm.shp)
#
##saving
#writeRaster(MF2020_avoiddeforest.px, "rasters/PGM/2020_avoiddeforest/MFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("MF2020_avoiddeforest.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [edgedist] distance to forest edge
## this variable is the euclidean distance between mature forest 
## to the nearest cell that is not MF
#
#inv.MF2020_avoiddeforest <- MF2020_avoiddeforest
#inv.MF2020_avoiddeforest[inv.MF2020_avoiddeforest==1]<-NA
##cheking
##inv.MF2020_avoiddeforest
##plot(inv.MF2020_avoiddeforest)
#
#edge.dist.2020_avoiddeforest <- distance(inv.MF2020_avoiddeforest)
###cheking
##edge.dist.2020_avoiddeforest
##anyNA(edge.dist.2020_avoiddeforest[])
##plot(edge.dist.2020_avoiddeforest)
#
#names(edge.dist.2020_avoiddeforest)<-"edgedist"
#edge.dist.2020_avoiddeforest[is.nan(edge.dist.2020_avoiddeforest)] <- 0
#edge.dist.2020_avoiddeforest <- mask(edge.dist.2020_avoiddeforest, pgm.shp)
#
##saving
#writeRaster(edge.dist.2020_avoiddeforest, "rasters/PGM/2020_avoiddeforest/edgedist.tif", format="GTiff", overwrite=T)
##writeRaster(edge.dist.2020_avoiddeforest, "rasters/PGM/2020_avoiddeforest_avoidboth/edgedist.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("inv.MF2020_avoiddeforest")])
#gc()
#
#
#
##
#
#
## [edge] forest edge
## this variable is the mean of edge area based on mature forest
#
## marking edge and core areas
#edge2020_avoiddeforest <- edge.dist.2020_avoiddeforest
#edge2020_avoiddeforest[edge2020_avoiddeforest>300] <- 0
#edge2020_avoiddeforest[edge2020_avoiddeforest!=0] <- 1
#
##saving
#writeRaster(edge2020_avoiddeforest, "rasters/PGM/input/edge2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#
##alternative method
##edge2020_avoiddeforest <- focal(MF2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun = mean, na.rm = T, pad = T)
##edge2020_avoiddeforest <- edge2020_avoiddeforest + MF2020_avoiddeforest
##edge2020_avoiddeforest[edge2020_avoiddeforest == 2] <- 0                  # core area
##edge2020_avoiddeforest[edge2020_avoiddeforest > 0 & edge2020_avoiddeforest < 2] <- 1    # edge
#
##cheking
##edge2020_avoiddeforest
##unique(edge2020_avoiddeforest[])
##plot(edge2020_avoiddeforest)
#
## mean edge in pixel scale (200m)
#edge2020_avoiddeforest.px <- focal(edge2020_avoiddeforest, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##edge2020_avoiddeforest.px
##anyNA(edge2020_avoiddeforest.px[])
##plot(edge2020_avoiddeforest.px)
#
#names(edge2020_avoiddeforest.px)<-"edgepx"
#edge2020_avoiddeforest.px[is.nan(edge2020_avoiddeforest.px)] <- 0
#edge2020_avoiddeforest.px <- mask(edge2020_avoiddeforest.px, pgm.shp)
#
##saving
#writeRaster(edge2020_avoiddeforest.px, "rasters/PGM/2020_avoiddeforest/edgepx.tif", format="GTiff", overwrite=T)
##
#
## mean edge in landscape scale (1000m)
#edge2020_avoiddeforest.ls <- focal(edge2020_avoiddeforest, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##edge2020_avoiddeforest.ls
##anyNA(edge2020_avoiddeforest.ls[])
##plot(edge2020_avoiddeforest.ls)
#
#names(edge2020_avoiddeforest.ls)<-"edgels"
#edge2020_avoiddeforest.ls[is.nan(edge2020_avoiddeforest.ls)] <- 0
#edge2020_avoiddeforest.ls <- mask(edge2020_avoiddeforest.ls, pgm.shp)
#
##saving
#writeRaster(edge2020_avoiddeforest.px, "rasters/PGM/2020_avoiddeforest/edgels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("edge2020_avoiddeforest.px", "edge2020_avoiddeforest.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
#################################################################################################################################.
#####         2020 Avoid deforestation             ####.
#####          primary forest only                 ####.
#
## [SFage] secondary forest age -- pixel: 3x3 (200m); and landscape: 11x11 (1000m)
## this variable is the mean age of secondary forest
#
#SFage2020_avoiddeforest2 <- SFage2020
##plot(SFage2020_avoiddeforest2)
#
##saving
##same as 2020 Real
#
## mean sf age in pixel scale (200m)
#SFage2020_avoiddeforest2.px <- focal(SFage2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##SFage2020_avoiddeforest2.px
##anyNA(SFage2020_avoiddeforest2.px[])
##plot(SFage2020_avoiddeforest2.px)
#
#names(SFage2020_avoiddeforest2.px)<-"SFagepx"
#SFage2020_avoiddeforest2.px[is.nan(SFage2020_avoiddeforest2.px)] <- 0
#SFage2020_avoiddeforest2.px <- mask(SFage2020_avoiddeforest2.px, pgm.shp)
#
##saving
#writeRaster(SFage2020_avoiddeforest2.px, "rasters/PGM/2020_avoiddeforest2/SFagepx.tif", format="GTiff", overwrite=T)
#
## mean sf age in landscape scale (1000m)
#SFage2020_avoiddeforest2.ls <- focal(SFage2020_avoiddeforest2, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##SFage2020_avoiddeforest2.ls
##anyNA(SFage2020_avoiddeforest2.ls[])
##plot(SFage2020_avoiddeforest2.ls)
#
#names(SFage2020_avoiddeforest2.ls)<-"SFagels"
#SFage2020_avoiddeforest2.ls[is.nan(SFage2020_avoiddeforest2.ls)] <- 0
#SFage2020_avoiddeforest2.ls <- mask(SFage2020_avoiddeforest2.ls, pgm.shp)
#
##saving
#writeRaster(SFage2020_avoiddeforest2.ls, "rasters/PGM/2020_avoiddeforest2/SFagels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("SFage2020_avoiddeforest2.px", "SFage2020_avoiddeforest2.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [SF] secondary forest
## this variable includes forest pixels in SFage raster
## which has less than 25 years for 2020_avoiddeforest2 or less than 35 for 2020_avoiddeforest2
#
#SF2020_avoiddeforest2 <- SFage2020_avoiddeforest2
#SF2020_avoiddeforest2[SF2020_avoiddeforest2>34] <- 0
#SF2020_avoiddeforest2[SF2020_avoiddeforest2>0] <- 1
#SF2020_avoiddeforest2[SF2020_avoiddeforest2<1] <- 0
#
##plot(SF2020_avoiddeforest2)
#
##saving
##same as 2020 Real
#
#
## mean sf cover in pixel scale (200m)
#SF2020_avoiddeforest2.px <- focal(SF2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##SF2020_avoiddeforest2.px
##anyNA(SF2020_avoiddeforest2.px[])
##plot(SF2020_avoiddeforest2.px)
#
#names(SF2020_avoiddeforest2.px)<-"SFpx"
#SF2020_avoiddeforest2.px[is.nan(SF2020_avoiddeforest2.px)] <- 0
#SF2020_avoiddeforest2.px <- mask(SF2020_avoiddeforest2.px, pgm.shp)
#
##saving
#writeRaster(SF2020_avoiddeforest2.px, "rasters/PGM/2020_avoiddeforest2/SFpx.tif", format="GTiff", overwrite=T)
#
## mean sf cover in landscape scale (1000m)
#SF2020_avoiddeforest2.ls <- focal(SF2020_avoiddeforest2, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##SF2020_avoiddeforest2.ls
##anyNA(SF2020_avoiddeforest2.ls[])
##plot(SF2020_avoiddeforest2.ls)
#
#names(SF2020_avoiddeforest2.ls)<-"SFls"
#SF2020_avoiddeforest2.ls[is.nan(SF2020_avoiddeforest2.ls)] <- 0
#SF2020_avoiddeforest2.ls <- mask(SF2020_avoiddeforest2.ls, pgm.shp)
#
##saving
#writeRaster(SF2020_avoiddeforest2.ls, "rasters/PGM/2020_avoiddeforest2/SFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("SF2020_avoiddeforest2.px", "SF2020_avoiddeforest2.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [TSD] time since degradation
## this variable is the mean time since a degradation event
#
#TSD2020_avoiddeforest2 <- pgm.degrad[["pgm.degrad.2020real"]]
#TSD2020_avoiddeforest2[get("SF2020_avoiddeforest2")[] == 1] <- 0
##plot(TSD2020_avoiddeforest2)
#
##saving
#writeRaster(TSD2020_avoiddeforest2, "rasters/PGM/input/TSD2020_avoiddeforest2.tif", format="GTiff", overwrite=T)
#
#
## mean tsd cover in pixel scale (200m)
#TSD2020_avoiddeforest2.px <- focal(TSD2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TSD2020_avoiddeforest2.px
##anyNA(TSD2020_avoiddeforest2.px[])
##plot(TSD2020_avoiddeforest2.px)
#
#names(TSD2020_avoiddeforest2.px)<-"TSDpx"
#TSD2020_avoiddeforest2.px[is.nan(TSD2020_avoiddeforest2.px)] <- 0
#TSD2020_avoiddeforest2.px <- mask(TSD2020_avoiddeforest2.px, pgm.shp)
#
##saving
#writeRaster(TSD2020_avoiddeforest2.px, "rasters/PGM/2020_avoiddeforest2/TSDpx.tif", format="GTiff", overwrite=T)
#
## mean tsd cover in landscape scale (1000m)
#TSD2020_avoiddeforest2.ls <- focal(TSD2020_avoiddeforest2, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##TSD2020_avoiddeforest2.ls
##anyNA(TSD2020_avoiddeforest2.ls[])
##plot(TSD2020_avoiddeforest2.ls)
#
#names(TSD2020_avoiddeforest2.ls)<-"TSDls"
#TSD2020_avoiddeforest2.ls[is.nan(TSD2020_avoiddeforest2.ls)] <- 0
#TSD2020_avoiddeforest2.ls <- mask(TSD2020_avoiddeforest2.ls, pgm.shp)
#
##saving
#writeRaster(TSD2020_avoiddeforest2.ls, "rasters/PGM/2020_avoiddeforest2/TSDls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TSD2020_avoiddeforest2.px", "TSD2020_avoiddeforest2.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [DPF] degraded primary forest
## this variable includes forest pixels in TSD raster
## which has less than 24 years for 2020_avoiddeforest2 or less than 34 for 2020_avoiddeforest2
#
#DPF2020_avoiddeforest2 <- TSD2020_avoiddeforest2
#DPF2020_avoiddeforest2[DPF2020_avoiddeforest2>33]<-NA
#DPF2020_avoiddeforest2[DPF2020_avoiddeforest2==0]<-NA
#DPF2020_avoiddeforest2[!is.na(DPF2020_avoiddeforest2)] <- 1
#DPF2020_avoiddeforest2<-sum(pgm.lulc.2020.forest.class, DPF2020_avoiddeforest2, na.rm=T)
#DPF2020_avoiddeforest2[DPF2020_avoiddeforest2!=2]<-0
#DPF2020_avoiddeforest2[DPF2020_avoiddeforest2==2]<-1
##plot(DPF2020_avoiddeforest2)
#
##saving
#writeRaster(DPF2020_avoiddeforest2, "rasters/PGM/input/DPF2020_avoiddeforest2.tif", format="GTiff", overwrite=T)
#
#
## mean dpf cover in pixel scale (200m)
#DPF2020_avoiddeforest2.px <- focal(DPF2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##DPF2020_avoiddeforest2.px
##anyNA(DPF2020_avoiddeforest2.px[])
##plot(DPF2020_avoiddeforest2.px)
#
#names(DPF2020_avoiddeforest2.px)<-"DPFpx"
#DPF2020_avoiddeforest2.px[is.nan(DPF2020_avoiddeforest2.px)] <- 0
#DPF2020_avoiddeforest2.px <- mask(DPF2020_avoiddeforest2.px, pgm.shp)
#
##saving
#writeRaster(DPF2020_avoiddeforest2.px, "rasters/PGM/2020_avoiddeforest2/DPFpx.tif", format="GTiff", overwrite=T)
#
## mean dpf cover in landscape scale (1000m)
#DPF2020_avoiddeforest2.ls <- focal(DPF2020_avoiddeforest2, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##DPF2020_avoiddeforest2.ls
##anyNA(DPF2020_avoiddeforest2.ls[])
##plot(DPF2020_avoiddeforest2.ls)
#
#names(DPF2020_avoiddeforest2.ls)<-"DPFls"
#DPF2020_avoiddeforest2.ls[is.nan(DPF2020_avoiddeforest2.ls)] <- 0
#DPF2020_avoiddeforest2.ls <- mask(DPF2020_avoiddeforest2.ls, pgm.shp)
#
##saving
#writeRaster(DPF2020_avoiddeforest2.ls, "rasters/PGM/2020_avoiddeforest2/DPFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("pgm.degrad.2020.forest.class", "DPF2020_avoiddeforest2.px", "DPF2020_avoiddeforest2.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [UPF] undisturbed primary forest
## this variable includes all forest pixels in LULC raster (value == 3)
## excluding degraded and secondary forest
#
#UPF2020_avoiddeforest2<-sum(pgm.lulc.2010.forest.class, SF2020_avoiddeforest2, DPF2020_avoiddeforest2, na.rm = T)
#UPF2020_avoiddeforest2[UPF2020_avoiddeforest2>1]<-0
###cheking
##unique(UPF2020_avoiddeforest2[])
##plot(UPF2020_avoiddeforest2)
#
##saving
#writeRaster(UPF2020_avoiddeforest2, "rasters/PGM/input/UPF2020_avoiddeforest2.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#UPF2020_avoiddeforest2.px <- focal(UPF2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##UPF2020_avoiddeforest2.px
##anyNA(UPF2020_avoiddeforest2.px[])
##plot(UPF2020_avoiddeforest2.px)
#
#names(UPF2020_avoiddeforest2.px)<-"UPFpx"
#UPF2020_avoiddeforest2.px[is.nan(UPF2020_avoiddeforest2.px)] <- 0
#UPF2020_avoiddeforest2.px <- mask(UPF2020_avoiddeforest2.px, pgm.shp)
#
##saving
#writeRaster(UPF2020_avoiddeforest2.px, "rasters/PGM/2020_avoiddeforest2/UPFpx.tif", format="GTiff", overwrite=T)
#
## mean upf cover in landscape scale (1000m)
#UPF2020_avoiddeforest2.ls <- focal(UPF2020_avoiddeforest2, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##UPF2020_avoiddeforest2.ls
##anyNA(UPF2020_avoiddeforest2.ls[])
##plot(UPF2020_avoiddeforest2.ls)
#
#names(UPF2020_avoiddeforest2.ls)<-"UPFls"
#UPF2020_avoiddeforest2.ls[is.nan(UPF2020_avoiddeforest2.ls)] <- 0
#UPF2020_avoiddeforest2.ls <- mask(UPF2020_avoiddeforest2.ls, pgm.shp)
#
##saving
#writeRaster(UPF2020_avoiddeforest2.ls, "rasters/PGM/2020_avoiddeforest2/UPFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("UPF2020_avoiddeforest2.px", "UPF2020_avoiddeforest2.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [F3] Forest type 3 or Total forest
## this variable includes forest pixels in LULC raster
## including degraded forest and secondary forest older than 2 years
#
#SF2020_avoiddeforest2.young <- SFage2020_avoiddeforest2
#SF2020_avoiddeforest2.young[SF2020_avoiddeforest2.young <= 2] <- 0
#SF2020_avoiddeforest2.young[SF2020_avoiddeforest2.young > 2] <- 1
#
##saving
##same as 2020 Real
#
#
#TF2020_avoiddeforest2 <- sum(UPF2020_avoiddeforest2, DPF2020_avoiddeforest2, SF2020_avoiddeforest2.young, na.rm = T)
#TF2020_avoiddeforest2[TF2020_avoiddeforest2>1] <- 1
###cheking
##TF2020_avoiddeforest2
##plot(TF2020_avoiddeforest2)
#
##saving
#writeRaster(TF2020_avoiddeforest2, "rasters/PGM/input/TF2020_avoiddeforest2.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#TF2020_avoiddeforest2.px <- focal(TF2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TF2020_avoiddeforest2.px
##anyNA(TF2020_avoiddeforest2.px[])
##plot(TF2020_avoiddeforest2.px)
#
#names(TF2020_avoiddeforest2.px)<-"TFpx"
#TF2020_avoiddeforest2.px[is.nan(TF2020_avoiddeforest2.px)] <- 0
#TF2020_avoiddeforest2.px <- mask(TF2020_avoiddeforest2.px, pgm.shp)
#
##saving
#writeRaster(TF2020_avoiddeforest2.px, "rasters/PGM/2020_avoiddeforest2/TFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TF2020_avoiddeforest2.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [F1] Forest type 1 or Mature forest
## this variable includes forest pixels in LULC raster
## including degraded forest and secondary forest older than 10 years
#
#SF2020_avoiddeforest2.mature <- SFage2020_avoiddeforest2
#SF2020_avoiddeforest2.mature[SF2020_avoiddeforest2.mature < 10] <- 0
#SF2020_avoiddeforest2.mature[SF2020_avoiddeforest2.mature >= 10] <- 1
#
##saving
##same as 2020 Real
#
#
#MF2020_avoiddeforest2 <- sum(UPF2020_avoiddeforest2, DPF2020_avoiddeforest2, SF2020_avoiddeforest2.mature, na.rm = T)
#MF2020_avoiddeforest2[MF2020_avoiddeforest2>1] <- 1
###cheking
##MF2020_avoiddeforest2
##plot(MF2020_avoiddeforest2)
#
##saving
#writeRaster(MF2020_avoiddeforest2, "rasters/PGM/input/MF2020_avoiddeforest2.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#MF2020_avoiddeforest2.px <- focal(MF2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##MF2020_avoiddeforest2.px
##anyNA(MF2020_avoiddeforest2.px[])
##plot(MF2020_avoiddeforest2.px)
#
#names(MF2020_avoiddeforest2.px)<-"MFpx"
#MF2020_avoiddeforest2.px[is.nan(MF2020_avoiddeforest2.px)] <- 0
#MF2020_avoiddeforest2.px <- mask(MF2020_avoiddeforest2.px, pgm.shp)
#
##saving
#writeRaster(MF2020_avoiddeforest2.px, "rasters/PGM/2020_avoiddeforest2/MFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("MF2020_avoiddeforest2.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## [edgedist] distance to forest edge
## this variable is the euclidean distance between mature forest 
## to the nearest cell that is not MF
#
#inv.MF2020_avoiddeforest2 <- MF2020_avoiddeforest2
#inv.MF2020_avoiddeforest2[inv.MF2020_avoiddeforest2==1]<-NA
##cheking
##inv.MF2020_avoiddeforest2
##plot(inv.MF2020_avoiddeforest2)
#
#edge.dist.2020_avoiddeforest2 <- distance(inv.MF2020_avoiddeforest2)
###cheking
##edge.dist.2020_avoiddeforest2
##anyNA(edge.dist.2020_avoiddeforest2[])
##plot(edge.dist.2020_avoiddeforest2)
#
#names(edge.dist.2020_avoiddeforest2)<-"edgedist"
#edge.dist.2020_avoiddeforest2[is.nan(edge.dist.2020_avoiddeforest2)] <- 0
#edge.dist.2020_avoiddeforest2 <- mask(edge.dist.2020_avoiddeforest2, pgm.shp)
#
##saving
#writeRaster(edge.dist.2020_avoiddeforest2, "rasters/PGM/2020_avoiddeforest2/edgedist.tif", format="GTiff", overwrite=T)
##writeRaster(edge.dist.2020_avoiddeforest2, "rasters/PGM/2020_avoiddeforest2_avoidboth/edgedist.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("inv.MF2020_avoiddeforest2")])
#gc()
#
#
#
##
#
#
## [edge] forest edge
## this variable is the mean of edge area based on mature forest
#
## marking edge and core areas
#edge2020_avoiddeforest2 <- edge.dist.2020_avoiddeforest2
#edge2020_avoiddeforest2[edge2020_avoiddeforest2>300] <- 0
#edge2020_avoiddeforest2[edge2020_avoiddeforest2!=0] <- 1
#
##saving
#writeRaster(edge2020_avoiddeforest2, "rasters/PGM/input/edge2020_avoiddeforest2.tif", format="GTiff", overwrite=T)
#
##alternative method
##edge2020_avoiddeforest2 <- focal(MF2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun = mean, na.rm = T, pad = T)
##edge2020_avoiddeforest2 <- edge2020_avoiddeforest2 + MF2020_avoiddeforest2
##edge2020_avoiddeforest2[edge2020_avoiddeforest2 == 2] <- 0                  # core area
##edge2020_avoiddeforest2[edge2020_avoiddeforest2 > 0 & edge2020_avoiddeforest2 < 2] <- 1    # edge
#
##cheking
##edge2020_avoiddeforest2
##unique(edge2020_avoiddeforest2[])
##plot(edge2020_avoiddeforest2)
#
## mean edge in pixel scale (200m)
#edge2020_avoiddeforest2.px <- focal(edge2020_avoiddeforest2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##edge2020_avoiddeforest2.px
##anyNA(edge2020_avoiddeforest2.px[])
##plot(edge2020_avoiddeforest2.px)
#
#names(edge2020_avoiddeforest2.px)<-"edgepx"
#edge2020_avoiddeforest2.px[is.nan(edge2020_avoiddeforest2.px)] <- 0
#edge2020_avoiddeforest2.px <- mask(edge2020_avoiddeforest2.px, pgm.shp)
#
##saving
#writeRaster(edge2020_avoiddeforest2.px, "rasters/PGM/2020_avoiddeforest2/edgepx.tif", format="GTiff", overwrite=T)
##
#
## mean edge in landscape scale (1000m)
#edge2020_avoiddeforest2.ls <- focal(edge2020_avoiddeforest2, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##edge2020_avoiddeforest2.ls
##anyNA(edge2020_avoiddeforest2.ls[])
##plot(edge2020_avoiddeforest2.ls)
#
#names(edge2020_avoiddeforest2.ls)<-"edgels"
#edge2020_avoiddeforest2.ls[is.nan(edge2020_avoiddeforest2.ls)] <- 0
#edge2020_avoiddeforest2.ls <- mask(edge2020_avoiddeforest2.ls, pgm.shp)
#
##saving
#writeRaster(edge2020_avoiddeforest2.px, "rasters/PGM/2020_avoiddeforest2/edgels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("edge2020_avoiddeforest2.px", "edge2020_avoiddeforest2.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
## [UPF] undisturbed primary forest -- pixel: 5x5 (200m); and landscape: 35x35 (1000m)
## this variable includes all forest pixels in LULC raster (value == 3)
## excluding those with age < 25 in 2010 SF raster or age < 35 in 2020 SF raster
## excluding pixels degraded
#
## scenario 2010
#UPF2010<-sum(pgm.lulc.2010.forest.class, pgm.sfage.2010.all.class, pgm.degrad.2010.forest.class, na.rm = T)
#UPF2010[UPF2010>1]<-0
###cheking
##unique(UPF2010[])
##plot(UPF2010)
#
##saving
#writeRaster(UPF2010, "rasters/PGM/input/UPF2010_real.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#UPF2010.px <- focal(UPF2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##UPF2010.px
##anyNA(UPF2010.px[])
##plot(UPF2010.px)
#
#names(UPF2010.px)<-"UPFpx"
#UPF2010.px[is.nan(UPF2010.px)] <- 0
#UPF2010.px <- mask(UPF2010.px, pgm.shp)
#
##saving
#writeRaster(UPF2010.px, "rasters/PGM/2010_real/UPFpx.tif", format="GTiff", overwrite=T)
#
## mean upf cover in landscape scale (1000m)
#UPF2010.ls <- focal(UPF2010, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##UPF2010.ls
##anyNA(UPF2010.ls[])
##plot(UPF2010.ls)
#
#names(UPF2010.ls)<-"UPFls"
#UPF2010.ls[is.nan(UPF2010.ls)] <- 0
#UPF2010.ls <- mask(UPF2010.ls, pgm.shp)
#
##saving
#writeRaster(UPF2010.ls, "rasters/PGM/2010_real/UPFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("UPF2010.px", "UPF2010.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario avoid degradation
#UPF.avoiddegrad<-sum(pgm.lulc.2020.forest.class, pgm.sfage.2020.all.class, pgm.degrad.2010.forest.class, na.rm = T)
#UPF.avoiddegrad[UPF.avoiddegrad>1]<-0
###cheking
##unique(UPF.avoiddegrad[])
##plot(UPF.avoiddegrad)
#
##saving
#writeRaster(UPF.avoiddegrad, "rasters/PGM/input/UPF2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#UPF.avoiddegrad.px <- focal(UPF.avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##UPF.avoiddegrad.px
##anyNA(UPF.avoiddegrad.px[])
##plot(UPF.avoiddegrad.px)
#
#names(UPF.avoiddegrad.px)<-"UPFpx"
#UPF.avoiddegrad.px[is.nan(UPF.avoiddegrad.px)] <- 0
#UPF.avoiddegrad.px <- mask(UPF.avoiddegrad.px, pgm.shp)
#
##saving
#writeRaster(UPF.avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/UPFpx.tif", format="GTiff", overwrite=T)
#
## mean upf cover in landscape scale (1000m)
#UPF.avoiddegrad.ls <- focal(UPF.avoiddegrad, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##UPF.avoiddegrad.ls
##anyNA(UPF.avoiddegrad.ls[])
##plot(UPF.avoiddegrad.ls)
#
#names(UPF.avoiddegrad.ls)<-"UPFls"
#UPF.avoiddegrad.ls[is.nan(UPF.avoiddegrad.ls)] <- 0
#UPF.avoiddegrad.ls <- mask(UPF.avoiddegrad.ls, pgm.shp)
#
##saving
#writeRaster(UPF.avoiddegrad.ls, "rasters/PGM/2020_avoiddegrad/UPFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("UPF.avoiddegrad.px", "UPF.avoiddegrad.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario avoid deforestation
#UPF.avoiddefor<-sum(pgm.lulc.2020.forest.class, pgm.sfage.2010.all.class, pgm.degrad.2020.forest.class, na.rm = T)
#UPF.avoiddefor[UPF.avoiddefor>1]<-0
###cheking
##unique(UPF.avoiddefor[])
##plot(UPF.avoiddefor)
#
##saving
#writeRaster(UPF.avoiddefor, "rasters/PGM/input/UPF2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#UPF.avoiddefor.px <- focal(UPF.avoiddefor, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##UPF.avoiddefor.px
##anyNA(UPF.avoiddefor.px[])
##plot(UPF.avoiddefor.px)
#
#names(UPF.avoiddefor.px)<-"UPFpx"
#UPF.avoiddefor.px[is.nan(UPF.avoiddefor.px)] <- 0
#UPF.avoiddefor.px <- mask(UPF.avoiddefor.px, pgm.shp)
#
##saving
#writeRaster(UPF.avoiddefor.px, "rasters/PGM/2020_avoiddeforest/UPFpx.tif", format="GTiff", overwrite=T)
#writeRaster(UPF.avoiddefor.px, "rasters/PGM/2020_restor_n_avoid_deforest/UPFpx.tif", format="GTiff", overwrite=T)
#
## mean upf cover in landscape scale (1000m)
#UPF.avoiddefor.ls <- focal(UPF.avoiddefor, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##UPF.avoiddefor.ls
##anyNA(UPF.avoiddefor.ls[])
##plot(UPF.avoiddefor.ls)
#
#names(UPF.avoiddefor.ls)<-"UPFls"
#UPF.avoiddefor.ls[is.nan(UPF.avoiddefor.ls)] <- 0
#UPF.avoiddefor.ls <- mask(UPF.avoiddefor.ls, pgm.shp)
#
##saving
#writeRaster(UPF.avoiddefor.ls, "rasters/PGM/2020_avoiddeforest/UPFls.tif", format="GTiff", overwrite=T)
#writeRaster(UPF.avoiddefor.ls, "rasters/PGM/2020_restor_n_avoid_deforest/UPFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("UPF.avoiddefor.px", "UPF.avoiddefor.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario avoid deforestation upf only
#UPF.avoiddefor2<-sum(pgm.lulc.2010.forest.class, pgm.sfage.2020.all.class, pgm.degrad.2020.forest.class, na.rm = T)
#UPF.avoiddefor2[UPF.avoiddefor2>1]<-0
###cheking
##unique(UPF.avoiddefor2[])
##plot(UPF.avoiddefor2)
#
##saving
#writeRaster(UPF.avoiddefor2, "rasters/PGM/input/UPF2020_avoiddeforest2.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#UPF.avoiddefor2.px <- focal(UPF.avoiddefor2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##UPF.avoiddefor2.px
##anyNA(UPF.avoiddefor2.px[])
##plot(UPF.avoiddefor2.px)
#
#names(UPF.avoiddefor2.px)<-"UPFpx"
#UPF.avoiddefor2.px[is.nan(UPF.avoiddefor2.px)] <- 0
#UPF.avoiddefor2.px <- mask(UPF.avoiddefor2.px, pgm.shp)
#
##saving
#writeRaster(UPF.avoiddefor2.px, "rasters/PGM/2020_avoiddeforest2/UPFpx.tif", format="GTiff", overwrite=T)
#
## mean upf cover in landscape scale (1000m)
#UPF.avoiddefor2.ls <- focal(UPF.avoiddefor2, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##UPF.avoiddefor2.ls
##anyNA(UPF.avoiddefor2.ls[])
##plot(UPF.avoiddefor2.ls)
#
#names(UPF.avoiddefor2.ls)<-"UPFls"
#UPF.avoiddefor2.ls[is.nan(UPF.avoiddefor2.ls)] <- 0
#UPF.avoiddefor2.ls <- mask(UPF.avoiddefor2.ls, pgm.shp)
#
##saving
#writeRaster(UPF.avoiddefor2.ls, "rasters/PGM/2020_avoiddeforest2/UPFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("UPF.avoiddefor2.px", "UPF.avoiddefor2.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario avoid both
#UPF.avoidboth<-sum(UPF.avoiddegrad, UPF.avoiddefor, na.rm = T)
#UPF.avoidboth[UPF.avoidboth>1]<-1
###cheking
##unique(UPF.avoidboth[])
##plot(UPF.avoidboth)
#
##saving
#writeRaster(UPF.avoidboth, "rasters/PGM/input/UPF2020_avoidboth.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#UPF.avoidboth.px <- focal(UPF.avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##UPF.avoidboth.px
##anyNA(UPF.avoidboth.px[])
##plot(UPF.avoidboth.px)
#
#names(UPF.avoidboth.px)<-"UPFpx"
#UPF.avoidboth.px[is.nan(UPF.avoidboth.px)] <- 0
#UPF.avoidboth.px <- mask(UPF.avoidboth.px, pgm.shp)
#
##saving
#writeRaster(UPF.avoidboth.px, "rasters/PGM/2020_avoidboth/UPFpx.tif", format="GTiff", overwrite=T)
#writeRaster(UPF.avoidboth.px, "rasters/PGM/2020_restor_n_avoid_both/UPFpx.tif", format="GTiff", overwrite=T)
#
## mean upf cover in landscape scale (1000m)
#UPF.avoidboth.ls <- focal(UPF.avoidboth, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##UPF.avoidboth.ls
##anyNA(UPF.avoidboth.ls[])
##plot(UPF.avoidboth.ls)
#
#names(UPF.avoidboth.ls)<-"UPFls"
#UPF.avoidboth.ls[is.nan(UPF.avoidboth.ls)] <- 0
#UPF.avoidboth.ls <- mask(UPF.avoidboth.ls, pgm.shp)
#
##saving
#writeRaster(UPF.avoidboth.ls, "rasters/PGM/2020_avoidboth/UPFls.tif", format="GTiff", overwrite=T)
#writeRaster(UPF.avoidboth.ls, "rasters/PGM/2020_restor_n_avoid_both/UPFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("UPF.avoidboth.px", "UPF.avoidboth.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario 2020
#UPF2020<-sum(pgm.lulc.2020.forest.class, pgm.sfage.2020.all.class, pgm.degrad.2020.forest.class, na.rm = T)
#UPF2020[UPF2020>1]<-0
###cheking
##UPF2020
##plot(UPF2020)
#
##saving
#writeRaster(UPF2020, "rasters/PGM/input/UPF2020_real.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#UPF2020.px <- focal(UPF2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##UPF2020.px
##anyNA(UPF2020.px[])
##plot(UPF2020.px)
#
#names(UPF2020.px)<-"UPFpx"
#UPF2020.px[is.nan(UPF2020.px)] <- 0
#UPF2020.px <- mask(UPF2020.px, pgm.shp)
#
##saving
#writeRaster(UPF2020.px, "rasters/PGM/2020_real/UPFpx.tif", format="GTiff", overwrite=T)
#writeRaster(UPF2020.px, "rasters/PGM/2020_restor_wo_avoid/UPFpx.tif", format="GTiff", overwrite=T)
#
## mean upf cover in landscape scale (1000m)
#UPF2020.ls <- focal(UPF2020, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##UPF2020.ls
##anyNA(UPF2020.ls[])
##plot(UPF2020.ls)
#
#names(UPF2020.ls)<-"UPFls"
#UPF2020.ls[is.nan(UPF2020.ls)] <- 0
#UPF2020.ls <- mask(UPF2020.ls, pgm.shp)
#
##saving
#writeRaster(UPF2020.ls, "rasters/PGM/2020_real/UPFls.tif", format="GTiff", overwrite=T)
#writeRaster(UPF2020.ls, "rasters/PGM/2020_restor_wo_avoid/UPFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("UPF2020.px", "UPF2020.ls")]) #keeping only raster stack
#gc()
#
#
##
#
#
########################################################################################################################.
#
## [DPF] degraded primary forest -- pixel: 5x5 (200m); and landscape: 35x35 (1000m)
## this variable includes forest pixels in LULC raster (value == 3)
## which overlaps with pixels with fire (burned) and/or pixels degraded (burned and logged / logged)
#
## scenario 2010
#DPF2010 <- pgm.degrad.2010.forest.class
##plot(DPF2010)
#
##saving
#writeRaster(DPF2010, "rasters/PGM/input/DPF2010_real.tif", format="GTiff", overwrite=T)
#
#
## mean dpf cover in pixel scale (200m)
#DPF2010.px <- focal(DPF2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##DPF2010.px
##anyNA(DPF2010.px[])
##plot(DPF2010.px)
#
#names(DPF2010.px)<-"DPFpx"
#DPF2010.px[is.nan(DPF2010.px)] <- 0
#DPF2010.px <- mask(DPF2010.px, pgm.shp)
#
##saving
#writeRaster(DPF2010.px, "rasters/PGM/2010_real/DPFpx.tif", format="GTiff", overwrite=T)
#writeRaster(DPF2010.px, "rasters/PGM/2020_avoiddegrad/DPFpx.tif", format="GTiff", overwrite=T)
#writeRaster(DPF2010.px, "rasters/PGM/2020_avoidboth/DPFpx.tif", format="GTiff", overwrite=T)
#writeRaster(DPF2010.px, "rasters/PGM/2020_restor_n_avoid_both/DPFpx.tif", format="GTiff", overwrite=T)
#
## mean dpf cover in landscape scale (1000m)
#DPF2010.ls <- focal(DPF2010, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##DPF2010.ls
##anyNA(DPF2010.ls[])
##plot(DPF2010.ls)
#
#names(DPF2010.ls)<-"DPFls"
#DPF2010.ls[is.nan(DPF2010.ls)] <- 0
#DPF2010.ls <- mask(DPF2010.ls, pgm.shp)
#
##saving
#writeRaster(DPF2010.ls, "rasters/PGM/2010_real/DPFls.tif", format="GTiff", overwrite=T)
#writeRaster(DPF2010.ls, "rasters/PGM/2020_avoiddegrad/DPFls.tif", format="GTiff", overwrite=T)
#writeRaster(DPF2010.ls, "rasters/PGM/2020_avoidboth/DPFls.tif", format="GTiff", overwrite=T)
#writeRaster(DPF2010.ls, "rasters/PGM/2020_restor_n_avoid_both/DPFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("pgm.degrad.2010.forest.class", "DPF2010.px", "DPF2010.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario 2020
#DPF2020 <- pgm.degrad.2020.forest.class
##plot(DPF2020)
#
##saving
#writeRaster(DPF2020, "rasters/PGM/input/DPF2020_real.tif", format="GTiff", overwrite=T)
#
#
## mean dpf cover in pixel scale (200m)
#DPF2020.px <- focal(DPF2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##DPF2020.px
##anyNA(DPF2020.px[])
##plot(DPF2020.px)
#
#names(DPF2020.px)<-"DPFpx"
#DPF2020.px[is.nan(DPF2020.px)] <- 0
#DPF2020.px <- mask(DPF2020.px, pgm.shp)
#
##saving
#writeRaster(DPF2020.px, "rasters/PGM/2020_real/DPFpx.tif", format="GTiff", overwrite=T)
#writeRaster(DPF2020.px, "rasters/PGM/2020_avoiddeforest/DPFpx.tif", format="GTiff", overwrite=T)
#writeRaster(DPF2020.px, "rasters/PGM/2020_avoiddeforest2/DPFpx.tif", format="GTiff", overwrite=T)
#writeRaster(DPF2020.px, "rasters/PGM/2020_restor_wo_avoid/DPFpx.tif", format="GTiff", overwrite=T)
#writeRaster(DPF2020.px, "rasters/PGM/2020_restor_n_avoid_deforest/DPFpx.tif", format="GTiff", overwrite=T)
#
## mean dpf cover in landscape scale (1000m)
#DPF2020.ls <- focal(DPF2020, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##DPF2020.ls
##anyNA(DPF2020.ls[])
##plot(DPF2020.ls)
#
#names(DPF2020.ls)<-"DPFls"
#DPF2020.ls[is.nan(DPF2020.ls)] <- 0
#DPF2020.ls <- mask(DPF2020.ls, pgm.shp)
#
##saving
#writeRaster(DPF2020.ls, "rasters/PGM/2020_real/DPFls.tif", format="GTiff", overwrite=T)
#writeRaster(DPF2020.ls, "rasters/PGM/2020_avoiddeforest/DPFls.tif", format="GTiff", overwrite=T)
#writeRaster(DPF2020.ls, "rasters/PGM/2020_avoiddeforest2/DPFls.tif", format="GTiff", overwrite=T)
#writeRaster(DPF2020.ls, "rasters/PGM/2020_restor_wo_avoid/DPFls.tif", format="GTiff", overwrite=T)
#writeRaster(DPF2020.ls, "rasters/PGM/2020_restor_n_avoid_deforest/DPFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("pgm.degrad.2020.forest.class", "DPF2020.px", "DPF2020.ls")]) #keeping only raster stack
#gc()
#
#
##
#
#
########################################################################################################################.
#
## [TSD] time since degradation -- pixel: 5x5 (200m); and landscape: 35x35 (1000m)
## this variable is the mean time since a degradation event
#
## scenario 2010
#TSD2010 <- pgm.degrad[["pgm.degrad.2010real"]]
##plot(TSD2010)
#
##saving
#writeRaster(TSD2010, "rasters/PGM/input/TSD2010_real.tif", format="GTiff", overwrite=T)
#
#
## mean tsd cover in pixel scale (200m)
#TSD2010.px <- focal(TSD2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TSD2010.px
##anyNA(TSD2010.px[])
##plot(TSD2010.px)
#
#names(TSD2010.px)<-"TSDpx"
#TSD2010.px[is.nan(TSD2010.px)] <- 0
#TSD2010.px <- mask(TSD2010.px, pgm.shp)
#
##saving
#writeRaster(TSD2010.px, "rasters/PGM/2010_real/TSDpx.tif", format="GTiff", overwrite=T)
#
## mean tsd cover in landscape scale (1000m)
#TSD2010.ls <- focal(TSD2010, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##TSD2010.ls
##anyNA(TSD2010.ls[])
##plot(TSD2010.ls)
#
#names(TSD2010.ls)<-"TSDls"
#TSD2010.ls[is.nan(TSD2010.ls)] <- 0
#TSD2010.ls <- mask(TSD2010.ls, pgm.shp)
#
##saving
#writeRaster(TSD2010.ls, "rasters/PGM/2010_real/TSDls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TSD2010.px", "TSD2010.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario 2020 -- avoid degradation
#TSD2010.recovery10 <- calc(TSD2010, fun=function(x){ifelse(is.na(x), x, ifelse(x==0, x, ifelse(x==300, x, x+10)))})
##plot(TSD2010.recovery10)
#
##saving
#writeRaster(TSD2010.recovery10, "rasters/PGM/input/TSD2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#
#
## mean dpf cover in pixel scale (200m)
#TSD2010.recovery10.px <- focal(TSD2010.recovery10, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TSD2010.recovery10.px
##anyNA(TSD2010.recovery10.px[])
##plot(TSD2010.recovery10.px)
#
#names(TSD2010.recovery10.px)<-"TSDpx"
#TSD2010.recovery10.px[is.nan(TSD2010.recovery10.px)] <- 0
#TSD2010.recovery10.px <- mask(TSD2010.recovery10.px, pgm.shp)
#
##saving
#writeRaster(TSD2010.recovery10.px, "rasters/PGM/2020_avoiddegrad/TSDpx.tif", format="GTiff", overwrite=T)
#writeRaster(TSD2010.recovery10.px, "rasters/PGM/2020_avoidboth/TSDpx.tif", format="GTiff", overwrite=T)
#writeRaster(TSD2010.recovery10.px, "rasters/PGM/2020_restor_n_avoid_both/TSDpx.tif", format="GTiff", overwrite=T)
#
## mean dpf cover in landscape scale (1000m)
#TSD2010.recovery10.ls <- focal(TSD2010.recovery10, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##TSD2010.recovery10.ls
##anyNA(TSD2010.recovery10.ls[])
##plot(TSD2010.recovery10.ls)
#
#names(TSD2010.recovery10.ls)<-"TSDls"
#TSD2010.recovery10.ls[is.nan(TSD2010.recovery10.ls)] <- 0
#TSD2010.recovery10.ls <- mask(TSD2010.recovery10.ls, pgm.shp)
#
##saving
#writeRaster(TSD2010.recovery10.ls, "rasters/PGM/2020_avoiddegrad/TSDls.tif", format="GTiff", overwrite=T)
#writeRaster(TSD2010.recovery10.ls, "rasters/PGM/2020_avoidboth/TSDls.tif", format="GTiff", overwrite=T)
#writeRaster(TSD2010.recovery10.ls, "rasters/PGM/2020_restor_n_avoid_both/TSDls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TSD2010.recovery10.px", "TSD2010.recovery10.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario 2020
#TSD2020 <- pgm.degrad[["pgm.degrad.2020real"]]
##plot(TSD2020)
#
##saving
#writeRaster(TSD2020, "rasters/PGM/input/TSD2020_real.tif", format="GTiff", overwrite=T)
#
#
## mean tsd cover in pixel scale (200m)
#TSD2020.px <- focal(TSD2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TSD2020.px
##anyNA(TSD2020.px[])
##plot(TSD2020.px)
#
#names(TSD2020.px)<-"TSDpx"
#TSD2020.px[is.nan(TSD2020.px)] <- 0
#TSD2020.px <- mask(TSD2020.px, pgm.shp)
#
##saving
#writeRaster(TSD2020.px, "rasters/PGM/2020_real/TSDpx.tif", format="GTiff", overwrite=T)
#writeRaster(TSD2020.px, "rasters/PGM/2020_avoiddeforest/TSDpx.tif", format="GTiff", overwrite=T)
#writeRaster(TSD2020.px, "rasters/PGM/2020_avoiddeforest2/TSDpx.tif", format="GTiff", overwrite=T)
#writeRaster(TSD2020.px, "rasters/PGM/2020_restor_wo_avoid/TSDpx.tif", format="GTiff", overwrite=T)
#writeRaster(TSD2020.px, "rasters/PGM/2020_restor_n_avoid_deforest/TSDpx.tif", format="GTiff", overwrite=T)
#
## mean tsd cover in landscape scale (1000m)
#TSD2020.ls <- focal(TSD2020, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##TSD2020.ls
##anyNA(TSD2020.ls[])
##plot(TSD2020.ls)
#
#names(TSD2020.ls)<-"TSDls"
#TSD2020.ls[is.nan(TSD2020.ls)] <- 0
#TSD2020.ls <- mask(TSD2020.ls, pgm.shp)
#
##saving
#writeRaster(TSD2020.ls, "rasters/PGM/2020_real/TSDls.tif", format="GTiff", overwrite=T)
#writeRaster(TSD2020.ls, "rasters/PGM/2020_avoiddeforest/TSDls.tif", format="GTiff", overwrite=T)
#writeRaster(TSD2020.ls, "rasters/PGM/2020_avoiddeforest2/TSDls.tif", format="GTiff", overwrite=T)
#writeRaster(TSD2020.ls, "rasters/PGM/2020_restor_wo_avoid/TSDls.tif", format="GTiff", overwrite=T)
#writeRaster(TSD2020.ls, "rasters/PGM/2020_restor_n_avoid_deforest/TSDls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TSD2020.px", "TSD2020.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
########################################################################################################################.
#
## [SF] secondary forest -- pixel: 5x5 (200m); and landscape: 35x35 (1000m)
## this variable includes forest pixels in SFage raster
## which has less than 25 years for 2010 or less than 35 for 2020
#
## scenario 2010
#SF2010 <- pgm.sfage.2010.all.class
##plot(SF2010)
#
##saving
#writeRaster(SF2010, "rasters/PGM/input/SF2010_real.tif", format="GTiff", overwrite=T)
#
#
## mean sf cover in pixel scale (200m)
#SF2010.px <- focal(SF2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##SF2010.px
##anyNA(SF2010.px[])
##plot(SF2010.px)
#
#names(SF2010.px)<-"SFpx"
#SF2010.px[is.nan(SF2010.px)] <- 0
#SF2010.px <- mask(SF2010.px, pgm.shp)
#
##saving
#writeRaster(SF2010.px, "rasters/PGM/2010_real/SFpx.tif", format="GTiff", overwrite=T)
#writeRaster(SF2010.px, "rasters/PGM/2020_avoiddeforest/SFpx.tif", format="GTiff", overwrite=T)
#writeRaster(SF2010.px, "rasters/PGM/2020_avoidboth/SFpx.tif", format="GTiff", overwrite=T)
#
## mean sf cover in landscape scale (1000m)
#SF2010.ls <- focal(SF2010, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##SF2010.ls
##anyNA(SF2010.ls[])
##plot(SF2010.ls)
#
#names(SF2010.ls)<-"SFls"
#SF2010.ls[is.nan(SF2010.ls)] <- 0
#SF2010.ls <- mask(SF2010.ls, pgm.shp)
#
##saving
#writeRaster(SF2010.ls, "rasters/PGM/2010_real/SFls.tif", format="GTiff", overwrite=T)
#writeRaster(SF2010.ls, "rasters/PGM/2020_avoiddeforest/SFls.tif", format="GTiff", overwrite=T)
#writeRaster(SF2010.ls, "rasters/PGM/2020_avoidboth/SFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("SF2010.px", "SF2010.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario 2020
#SF2020 <- pgm.sfage.2020.all.class
##plot(SF2020)
#
##saving
#writeRaster(SF2020, "rasters/PGM/input/SF2020_real.tif", format="GTiff", overwrite=T)
#
#
## mean sf cover in pixel scale (200m)
#SF2020.px <- focal(SF2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##SF2020.px
##anyNA(SF2020.px[])
##plot(SF2020.px)
#
#names(SF2020.px)<-"SFpx"
#SF2020.px[is.nan(SF2020.px)] <- 0
#SF2020.px <- mask(SF2020.px, pgm.shp)
#
##saving
#writeRaster(SF2020.px, "rasters/PGM/2020_real/SFpx.tif", format="GTiff", overwrite=T)
#writeRaster(SF2020.px, "rasters/PGM/2020_avoiddegrad/SFpx.tif", format="GTiff", overwrite=T)
#writeRaster(SF2020.px, "rasters/PGM/2020_avoiddeforest2/SFpx.tif", format="GTiff", overwrite=T)
#
#
## mean sf cover in landscape scale (1000m)
#SF2020.ls <- focal(SF2020, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##SF2020.ls
##anyNA(SF2020.ls[])
##plot(SF2020.ls)
#
#names(SF2020.ls)<-"SFls"
#SF2020.ls[is.nan(SF2020.ls)] <- 0
#SF2020.ls <- mask(SF2020.ls, pgm.shp)
#
##saving
#writeRaster(SF2020.ls, "rasters/PGM/2020_real/SFls.tif", format="GTiff", overwrite=T)
#writeRaster(SF2020.ls, "rasters/PGM/2020_avoiddegrad/SFls.tif", format="GTiff", overwrite=T)
#writeRaster(SF2020.ls, "rasters/PGM/2020_avoiddeforest2/SFls.tif", format="GTiff", overwrite=T)
#
#
#rm(list=ls()[ls() %in% c("SF2020.px", "SF2020.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario -- restoration AND avoid
#SF2010.restore10 <- sum(SF2010, candidate.areas.final, na.rm = T)
#SF2010.restore10[SF2010.restore10>1]<-1
##plot(SF2010.restore10)
#
##saving
#writeRaster(SF2010.restore10, "rasters/PGM/input/SF2020_restor_n_avoid.tif", format="GTiff", overwrite=T)
#
#
## mean sf cover in pixel scale (200m)
#SF2010.restore10.px <- focal(SF2010.restore10, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##SF2010.restore10.px
##anyNA(SF2010.restore10.px[])
##plot(SF2010.restore10.px)
#
#names(SF2010.restore10.px)<-"SFpx"
#SF2010.restore10.px[is.nan(SF2010.restore10.px)] <- 0
#SF2010.restore10.px <- mask(SF2010.restore10.px, pgm.shp)
#
##saving
#writeRaster(SF2010.restore10.px, "rasters/PGM/2020_restor_n_avoid_deforest/SFpx.tif", format="GTiff", overwrite=T)
#writeRaster(SF2010.restore10.px, "rasters/PGM/2020_restor_n_avoid_both/SFpx.tif", format="GTiff", overwrite=T)
#
## mean sf cover in landscape scale (1000m)
#SF2010.restore10.ls <- focal(SF2010.restore10, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##SF2010.restore10.ls
##anyNA(SF2010.restore10.ls[])
##plot(SF2010.restore10.ls)
#
#names(SF2010.restore10.ls)<-"SFls"
#SF2010.restore10.ls[is.nan(SF2010.restore10.ls)] <- 0
#SF2010.restore10.ls <- mask(SF2010.restore10.ls, pgm.shp)
#
##saving
#writeRaster(SF2010.restore10.ls, "rasters/PGM/2020_restor_n_avoid_deforest/SFls.tif", format="GTiff", overwrite=T)
#writeRaster(SF2010.restore10.ls, "rasters/PGM/2020_restor_n_avoid_both/SFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("SF2010.restore10.px", "SF2010.restore10.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario -- restoration WITHOUT avoid
#SF2020.restore10 <- sum(SF2020, candidate.areas.final, na.rm = T)
#SF2020.restore10[SF2020.restore10>1]<-1
##plot(SF2020.restore10)
#
##saving
#writeRaster(SF2020.restore10, "rasters/PGM/input/SF2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)
#
#
## mean sf cover in pixel scale (200m)
#SF2020.restore10.px <- focal(SF2020.restore10, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##SF2020.restore10.px
##anyNA(SF2020.restore10.px[])
##plot(SF2020.restore10.px)
#
#names(SF2020.restore10.px)<-"SFpx"
#SF2020.restore10.px[is.nan(SF2020.restore10.px)] <- 0
#SF2020.restore10.px <- mask(SF2020.restore10.px, pgm.shp)
#
##saving
#writeRaster(SF2020.restore10.px, "rasters/PGM/2020_restor_wo_avoid/SFpx.tif", format="GTiff", overwrite=T)
#
## mean sf cover in landscape scale (1000m)
#SF2020.restore10.ls <- focal(SF2020.restore10, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##SF2020.restore10.ls
##anyNA(SF2020.restore10.ls[])
##plot(SF2020.restore10.ls)
#
#names(SF2020.restore10.ls)<-"SFls"
#SF2020.restore10.ls[is.nan(SF2020.restore10.ls)] <- 0
#SF2020.restore10.ls <- mask(SF2020.restore10.ls, pgm.shp)
#
##saving
#writeRaster(SF2020.restore10.ls, "rasters/PGM/2020_restor_wo_avoid/SFls.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("SF2020.restore10.px", "SF2020.restore10.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
########################################################################################################################.
#
## [SFage] secondary forest age -- pixel: 5x5 (200m); and landscape: 35x35 (1000m)
## this variable is the mean age of secondary forest
#
## scenario 2010
#SFage2010 <- pgm.sfage[["pgm.sfage.2010real"]]
##plot(SFage2010)
#
##saving
#writeRaster(SFage2010, "rasters/PGM/input/SFage2010_real.tif", format="GTiff", overwrite=T)
#
#
## mean sf age in pixel scale (200m)
#SFage2010.px <- focal(SFage2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##SFage2010.px
##anyNA(SFage2010.px[])
##plot(SFage2010.px)
#
#names(SFage2010.px)<-"SFagepx"
#SFage2010.px[is.nan(SFage2010.px)] <- 0
#SFage2010.px <- mask(SFage2010.px, pgm.shp)
#
##saving
#writeRaster(SFage2010.px, "rasters/PGM/2010_real/SFagepx.tif", format="GTiff", overwrite=T)
#
## mean sf age in landscape scale (1000m)
#SFage2010.ls <- focal(SFage2010, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##SFage2010.ls
##anyNA(SFage2010.ls[])
##plot(SFage2010.ls)
#
#names(SFage2010.ls)<-"SFagels"
#SFage2010.ls[is.nan(SFage2010.ls)] <- 0
#SFage2010.ls <- mask(SFage2010.ls, pgm.shp)
#
##saving
#writeRaster(SFage2010.ls, "rasters/PGM/2010_real/SFagels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("SFage2010.px", "SFage2010.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario 2020 -- avoid deforestation
#SFage2010.recovery10 <- calc(SFage2010, fun=function(x){ifelse(is.na(x), x, ifelse(x==0, x, x+10))})
##plot(SFage2010.recovery10)
#
##saving
#writeRaster(SFage2010.recovery10, "rasters/PGM/input/SFage2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#
#
## mean dpf cover in pixel scale (200m)
#SFage2010.recovery10.px <- focal(SFage2010.recovery10, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##SFage2010.recovery10.px
##anyNA(SFage2010.recovery10.px[])
##plot(SFage2010.recovery10.px)
#
#names(SFage2010.recovery10.px)<-"SFagepx"
#SFage2010.recovery10.px[is.nan(SFage2010.recovery10.px)] <- 0
#SFage2010.recovery10.px <- mask(SFage2010.recovery10.px, pgm.shp)
#
##saving
#writeRaster(SFage2010.recovery10.px, "rasters/PGM/2020_avoiddeforest/SFagepx.tif", format="GTiff", overwrite=T)
#writeRaster(SFage2010.recovery10.px, "rasters/PGM/2020_avoidboth/SFagepx.tif", format="GTiff", overwrite=T)
#
## mean dpf cover in landscape scale (1000m)
#SFage2010.recovery10.ls <- focal(SFage2010.recovery10, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##SFage2010.recovery10.ls
##anyNA(SFage2010.recovery10.ls[])
##plot(SFage2010.recovery10.ls)
#
#names(SFage2010.recovery10.ls)<-"SFagels"
#SFage2010.recovery10.ls[is.nan(SFage2010.recovery10.ls)] <- 0
#SFage2010.recovery10.ls <- mask(SFage2010.recovery10.ls, pgm.shp)
#
##saving
#writeRaster(SFage2010.recovery10.ls, "rasters/PGM/2020_avoiddeforest/SFagels.tif", format="GTiff", overwrite=T)
#writeRaster(SFage2010.recovery10.ls, "rasters/PGM/2020_avoidboth/SFagels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("SFage2010.recovery10.px", "SFage2010.recovery10.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario 2020
#SFage2020 <- pgm.sfage[["pgm.sfage.2020real"]]
##plot(SFage2010)
#
##saving
#writeRaster(SFage2020, "rasters/PGM/input/SFage2020_real.tif", format="GTiff", overwrite=T)
#
#
## mean sf age in pixel scale (200m)
#SFage2020.px <- focal(SFage2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##SFage2020.px
##anyNA(SFage2020.px[])
##plot(SFage2020.px)
#
#names(SFage2020.px)<-"SFagepx"
#SFage2020.px[is.nan(SFage2020.px)] <- 0
#SFage2020.px <- mask(SFage2020.px, pgm.shp)
#
##saving
#writeRaster(SFage2020.px, "rasters/PGM/2020_real/SFagepx.tif", format="GTiff", overwrite=T)
#writeRaster(SFage2020.px, "rasters/PGM/2020_avoiddegrad/SFagepx.tif", format="GTiff", overwrite=T)
#writeRaster(SFage2020.px, "rasters/PGM/2020_avoiddeforest2/SFagepx.tif", format="GTiff", overwrite=T)
#
#
## mean sf cover in landscape scale (1000m)
#SFage2020.ls <- focal(SFage2020, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##SFage2020.ls
##anyNA(SFage2020.ls[])
##plot(SFage2020.ls)
#
#names(SFage2020.ls)<-"SFagels"
#SFage2020.ls[is.nan(SFage2020.ls)] <- 0
#SFage2020.ls <- mask(SFage2020.ls, pgm.shp)
#
##saving
#writeRaster(SFage2020.ls, "rasters/PGM/2020_real/SFagels.tif", format="GTiff", overwrite=T)
#writeRaster(SFage2020.ls, "rasters/PGM/2020_avoiddegrad/SFagels.tif", format="GTiff", overwrite=T)
#writeRaster(SFage2020.ls, "rasters/PGM/2020_avoiddeforest2/SFagels.tif", format="GTiff", overwrite=T)
#
#
#rm(list=ls()[ls() %in% c("SFage2020.px", "SFage2020.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario -- restoration AND avoid
#candidate.areas.final.age <- calc(candidate.areas.final, fun=function(x){ifelse(x==1, x+9, x)})
##plot(candidate.areas.final.age)
#
#SFAge2010.restore10 <- sum(SFage2010.recovery10, candidate.areas.final.age, na.rm = T)
#values(SFAge2010.restore10)[values(SFAge2010.restore10)>=35]<-35
##plot(SFAge2010.restore10)
#
##saving
#writeRaster(SFAge2010.restore10, "rasters/PGM/input/SFage2020_restor_n_avoid.tif", format="GTiff", overwrite=T)
#
#
## mean sf cover in pixel scale (200m)
#SFAge2010.restore10.px <- focal(SFAge2010.restore10, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##SFAge2010.restore10.px
##anyNA(SFAge2010.restore10.px[])
##plot(SFAge2010.restore10.px)
#
#names(SFAge2010.restore10.px)<-"SFagepx"
#SFAge2010.restore10.px[is.nan(SFAge2010.restore10.px)] <- 0
#SFAge2010.restore10.px <- mask(SFAge2010.restore10.px, pgm.shp)
#
##saving
#writeRaster(SFAge2010.restore10.px, "rasters/PGM/2020_restor_n_avoid_deforest/SFagepx.tif", format="GTiff", overwrite=T)
#writeRaster(SFAge2010.restore10.px, "rasters/PGM/2020_restor_n_avoid_both/SFagepx.tif", format="GTiff", overwrite=T)
#
## mean sf cover in landscape scale (1000m)
#SFAge2010.restore10.ls <- focal(SFAge2010.restore10, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##SFAge2010.restore10.ls
##anyNA(SFAge2010.restore10.ls[])
##plot(SFAge2010.restore10.ls)
#
#names(SFAge2010.restore10.ls)<-"SFagels"
#SFAge2010.restore10.ls[is.nan(SFAge2010.restore10.ls)] <- 0
#SFAge2010.restore10.ls <- mask(SFAge2010.restore10.ls, pgm.shp)
#
##saving
#writeRaster(SFAge2010.restore10.ls, "rasters/PGM/2020_restor_n_avoid_deforest/SFagels.tif", format="GTiff", overwrite=T)
#writeRaster(SFAge2010.restore10.ls, "rasters/PGM/2020_restor_n_avoid_both/SFagels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("SFAge2010.restore10.px", "SFAge2010.restore10.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario -- restoration WITHOUT avoid
#SFAge2020.restore10 <- sum(SFage2020, candidate.areas.final.age, na.rm = T)
#values(SFAge2010.restore10)[values(SFAge2010.restore10)>=35]<-35
##plot(SFAge2020.restore10)
#
##saving
#writeRaster(SFAge2020.restore10, "rasters/PGM/input/SFage2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)
#
#
## mean sf cover in pixel scale (200m)
#SFAge2020.restore10.px <- focal(SFAge2020.restore10, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##SFAge2020.restore10.px
##anyNA(SFAge2020.restore10.px[])
##plot(SFAge2020.restore10.px)
#
#names(SFAge2020.restore10.px)<-"SFagepx"
#SFAge2020.restore10.px[is.nan(SFAge2020.restore10.px)] <- 0
#SFAge2020.restore10.px <- mask(SFAge2020.restore10.px, pgm.shp)
#
##saving
#writeRaster(SFAge2020.restore10.px, "rasters/PGM/2020_restor_wo_avoid/SFagepx.tif", format="GTiff", overwrite=T)
#
## mean sf cover in landscape scale (1000m)
#SFAge2020.restore10.ls <- focal(SFAge2020.restore10, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##SFAge2020.restore10.ls
##anyNA(SFAge2020.restore10.ls[])
##plot(SFAge2020.restore10.ls)
#
#names(SFAge2020.restore10.ls)<-"SFagels"
#SFAge2020.restore10.ls[is.nan(SFAge2020.restore10.ls)] <- 0
#SFAge2020.restore10.ls <- mask(SFAge2020.restore10.ls, pgm.shp)
#
##saving
#writeRaster(SFAge2020.restore10.ls, "rasters/PGM/2020_restor_wo_avoid/SFagels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("SFAge2020.restore10.px", "SFAge2020.restore10.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
########################################################################################################################.
#
## [F3] Forest type 3 or Total forest -- pixel: 5x5 (200m)
## this variable includes forest pixels in LULC raster
## including degraded forest and secondary forest older than 2 years
#
## scenario 2010
#SF2010.young <- SFage2010
#SF2010.young[SF2010.young <= 2] <- 0
#SF2010.young[SF2010.young > 2] <- 1
#
##saving
#writeRaster(SF2010.young, "rasters/PGM/input/SF2010_real-young.tif", format="GTiff", overwrite=T)
#
#
#TF2010 <- sum(UPF2010, DPF2010, SF2010.young, na.rm = T)
#TF2010[TF2010>1] <- 1
###cheking
##TF2010
##plot(TF2010)
#
##saving
#writeRaster(TF2010, "rasters/PGM/input/TF2010_real.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#TF2010.px <- focal(TF2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TF2010.px
##anyNA(TF2010.px[])
##plot(TF2010.px)
#
#names(TF2010.px)<-"TFpx"
#TF2010.px[is.nan(TF2010.px)] <- 0
#TF2010.px <- mask(TF2010.px, pgm.shp)
#
##saving
#writeRaster(TF2010.px, "rasters/PGM/2010_real/TFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TF2010.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario 2020
#SF2020.young <- SFage2020
#SF2020.young[SF2020.young <= 2] <- 0
#SF2020.young[SF2020.young > 2] <- 1
#
##saving
#writeRaster(SF2020.young, "rasters/PGM/input/SF2020_real-young.tif", format="GTiff", overwrite=T)
#
#
#TF2020 <- sum(UPF2020, DPF2020, SF2020.young, na.rm = T)
#TF2020[TF2020>1] <- 1
###cheking
##TF2020
##plot(TF2020)
#
##saving
#writeRaster(TF2020, "rasters/PGM/input/TF2020_real.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#TF2020.px <- focal(TF2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TF2020.px
##anyNA(TF2020.px[])
##plot(TF2020.px)
#
#names(TF2020.px)<-"TFpx"
#TF2020.px[is.nan(TF2020.px)] <- 0
#TF2020.px <- mask(TF2020.px, pgm.shp)
#
##saving
#writeRaster(TF2020.px, "rasters/PGM/2020_real/TFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TF2020.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario avoid degradation
#TF.avoiddegrad <- sum(UPF.avoiddegrad, DPF2010, SF2020.young, na.rm = T)
#TF.avoiddegrad[TF.avoiddegrad>1] <- 1
###cheking
##TF.avoiddegrad
##plot(TF.avoiddegrad)
#
##saving
#writeRaster(TF.avoiddegrad, "rasters/PGM/input/TF2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#TF.avoiddegrad.px <- focal(TF.avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TF.avoiddegrad.px
##anyNA(TF.avoiddegrad.px[])
##plot(TF.avoiddegrad.px)
#
#names(TF.avoiddegrad.px)<-"TFpx"
#TF.avoiddegrad.px[is.nan(TF.avoiddegrad.px)] <- 0
#TF.avoiddegrad.px <- mask(TF.avoiddegrad.px, pgm.shp)
#
##saving
#writeRaster(TF.avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/TFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TF.avoiddegrad.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario avoid deforestation
#TF.avoiddefor <- sum(UPF.avoiddefor, DPF2020, SF2010.young, na.rm = T)
#TF.avoiddefor[TF.avoiddefor>1] <- 1
###cheking
##TF.avoiddefor
##plot(TF.avoiddefor)
#
##saving
#writeRaster(TF.avoiddefor, "rasters/PGM/input/TF2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#TF.avoiddefor.px <- focal(TF.avoiddefor, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TF.avoiddefor.px
##anyNA(TF.avoiddefor.px[])
##plot(TF.avoiddefor.px)
#
#names(TF.avoiddefor.px)<-"TFpx"
#TF.avoiddefor.px[is.nan(TF.avoiddefor.px)] <- 0
#TF.avoiddefor.px <- mask(TF.avoiddefor.px, pgm.shp)
#
##saving
#writeRaster(TF.avoiddefor.px, "rasters/PGM/2020_avoiddeforest/TFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TF.avoiddefor.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario avoid deforestation upf only
#TF.avoiddefor2 <- sum(UPF.avoiddefor2, DPF2020, SF2020.young, na.rm = T)
#TF.avoiddefor2[TF.avoiddefor2>1] <- 1
###cheking
##TF.avoiddefor2
##plot(TF.avoiddefor2)
#
##saving
#writeRaster(TF.avoiddefor2, "rasters/PGM/input/TF2020_avoiddeforest2.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#TF.avoiddefor2.px <- focal(TF.avoiddefor2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TF.avoiddefor2.px
##anyNA(TF.avoiddefor2.px[])
##plot(TF.avoiddefor2.px)
#
#names(TF.avoiddefor2.px)<-"TFpx"
#TF.avoiddefor2.px[is.nan(TF.avoiddefor2.px)] <- 0
#TF.avoiddefor2.px <- mask(TF.avoiddefor2.px, pgm.shp)
#
##saving
#writeRaster(TF.avoiddefor2.px, "rasters/PGM/2020_avoiddeforest2/TFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TF.avoiddefor.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario avoid both
#TF.avoidboth <- sum(TF.avoiddegrad, TF.avoiddefor, na.rm = T)
#TF.avoidboth[TF.avoidboth>1] <- 1
###cheking
##TF.avoidboth
##plot(TF.avoidboth)
#
##saving
#writeRaster(TF.avoidboth, "rasters/PGM/input/TF2020_avoidboth.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#TF.avoidboth.px <- focal(TF.avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TF.avoidboth.px
##anyNA(TF.avoidboth.px[])
##plot(TF.avoidboth.px)
#
#names(TF.avoidboth.px)<-"TFpx"
#TF.avoidboth.px[is.nan(TF.avoidboth.px)] <- 0
#TF.avoidboth.px <- mask(TF.avoidboth.px, pgm.shp)
#
##saving
#writeRaster(TF.avoidboth.px, "rasters/PGM/2020_avoidboth/TFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TF.avoidboth.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario restoration WITHOUT avoid
#SFAge2020.restore10.young <- SFAge2020.restore10
#SFAge2020.restore10.young[SFAge2020.restore10.young <= 2] <- 0
#SFAge2020.restore10.young[SFAge2020.restore10.young > 2] <- 1
#
##saving
#writeRaster(SFAge2020.restore10.young, "rasters/PGM/input/SFAge2020_restor_wo_avoid-young.tif", format="GTiff", overwrite=T)
#
#
#TF.restore10.a <- sum(UPF2020, DPF2020, SFAge2020.restore10.young, na.rm = T)
#TF.restore10.a[TF.restore10.a>1] <- 1
###cheking
##TF.restore10.a
##plot(TF.restore10.a)
#
##saving
#writeRaster(TF.restore10.a, "rasters/PGM/input/TF2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#TF.restore10.a.px <- focal(TF.restore10.a, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TF.restore10.a.px
##anyNA(TF.restore10.a.px[])
##plot(TF.restore10.a.px)
#
#names(TF.restore10.a.px)<-"TFpx"
#TF.restore10.a.px[is.nan(TF.restore10.a.px)] <- 0
#TF.restore10.a.px <- mask(TF.restore10.a.px, pgm.shp)
#
##saving
#writeRaster(TF.restore10.a.px, "rasters/PGM/2020_restor_wo_avoid/TFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TF.restore10.a.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario restoration AND avoid deforestation
#SFAge2010.restore10.young <- SFAge2010.restore10
#SFAge2010.restore10.young[SFAge2010.restore10.young <= 2] <- 0
#SFAge2010.restore10.young[SFAge2010.restore10.young > 2] <- 1
#
##saving
#writeRaster(SFAge2010.restore10.young, "rasters/PGM/input/SFAge2020_restor_n_avoid_deforest-young.tif", format="GTiff", overwrite=T)
#
#
#TF.restore10.b <- sum(UPF.avoiddefor, DPF2020, SFAge2010.restore10.young, na.rm = T)
#TF.restore10.b[TF.restore10.b>1] <- 1
###cheking
##TF.restore10.b
##plot(TF.restore10.b)
#
##saving
#writeRaster(TF.restore10.b, "rasters/PGM/input/TF2020_restor_n_avoid_deforest.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#TF.restore10.b.px <- focal(TF.restore10.b, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TF.restore10.b.px
##anyNA(TF.restore10.b.px[])
##plot(TF.restore10.b.px)
#
#names(TF.restore10.b.px)<-"TFpx"
#TF.restore10.b.px[is.nan(TF.restore10.b.px)] <- 0
#TF.restore10.b.px <- mask(TF.restore10.b.px, pgm.shp)
#
##saving
#writeRaster(TF.restore10.b.px, "rasters/PGM/2020_restor_n_avoid_deforest/TFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TF.restore10.b.pxv")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario restoration AND avoid both
#TF.restore10.c <- sum(UPF.avoidboth, DPF2010, SFAge2010.restore10.young, na.rm = T)
#TF.restore10.c[TF.restore10.c>1] <- 1
###cheking
##TF.restore10.c
##plot(TF.restore10.c)
#
##saving
#writeRaster(TF.restore10.c, "rasters/PGM/input/TF2020_restor_n_avoid_both.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#TF.restore10.c.px <- focal(TF.restore10.c, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##TF.restore10.c.px
##anyNA(TF.restore10.c.px[])
##plot(TF.restore10.c.px)
#
#names(TF.restore10.c.px)<-"TFpx"
#TF.restore10.c.px[is.nan(TF.restore10.c.px)] <- 0
#TF.restore10.c.px <- mask(TF.restore10.c.px, pgm.shp)
#
##saving
#writeRaster(TF.restore10.c.px, "rasters/PGM/2020_restor_n_avoid_both/TFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("TF.restore10.c.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
########################################################################################################################.
#
## [F1] Forest type 1 or Mature forest -- pixel: 5x5 (200m)
## this variable includes forest pixels in LULC raster
## including degraded forest and secondary forest older than 10 years
#
## scenario 2010
#SF2010.mature <- SFage2010
#SF2010.mature[SF2010.mature < 10] <- 0
#SF2010.mature[SF2010.mature >= 10] <- 1
#
##saving
#writeRaster(SF2010.mature, "rasters/PGM/input/SFAge2010_real-mature.tif", format="GTiff", overwrite=T)
#
#
#MF2010 <- sum(UPF2010, DPF2010, SF2010.mature, na.rm = T)
#MF2010[MF2010>1] <- 1
###cheking
##MF2010
##plot(MF2010)
#
##saving
#writeRaster(MF2010, "rasters/PGM/input/MF2010_real.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#MF2010.px <- focal(MF2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##MF2010.px
##anyNA(MF2010.px[])
##plot(MF2010.px)
#
#names(MF2010.px)<-"MFpx"
#MF2010.px[is.nan(MF2010.px)] <- 0
#MF2010.px <- mask(MF2010.px, pgm.shp)
#
##saving
#writeRaster(MF2010.px, "rasters/PGM/2010_real/MFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("MF2010.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario 2020
#SF2020.mature <- SFage2020
#SF2020.mature[SF2020.mature < 10] <- 0
#SF2020.mature[SF2020.mature >= 10] <- 1
#
##saving
#writeRaster(SF2020.mature, "rasters/PGM/input/SFAge2020_real-mature.tif", format="GTiff", overwrite=T)
#
#
#MF2020 <- sum(UPF2020, DPF2020, SF2020.mature, na.rm = T)
#MF2020[MF2020>1] <- 1
###cheking
##MF2020
##plot(MF2020)
#
##saving
#writeRaster(MF2020, "rasters/PGM/input/MF2020_real.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#MF2020.px <- focal(MF2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##MF2020.px
##anyNA(MF2020.px[])
##plot(MF2020.px)
#
#names(MF2020.px)<-"MFpx"
#MF2020.px[is.nan(MF2020.px)] <- 0
#MF2020.px <- mask(MF2020.px, pgm.shp)
#
##saving
#writeRaster(MF2020.px, "rasters/PGM/2020_real/MFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("MF2020.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario avoid degradation
#MF.avoiddegrad <- sum(UPF.avoiddegrad, DPF2010, SF2020.mature, na.rm = T)
#MF.avoiddegrad[MF.avoiddegrad>1] <- 1
###cheking
##MF.avoiddegrad
##plot(MF.avoiddegrad)
#
##saving
#writeRaster(MF.avoiddegrad, "rasters/PGM/input/MF2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#MF.avoiddegrad.px <- focal(MF.avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##MF.avoiddegrad.px
##anyNA(MF.avoiddegrad.px[])
##plot(MF.avoiddegrad.px)
#
#names(MF.avoiddegrad.px)<-"MFpx"
#MF.avoiddegrad.px[is.nan(MF.avoiddegrad.px)] <- 0
#MF.avoiddegrad.px <- mask(MF.avoiddegrad.px, pgm.shp)
#
##saving
#writeRaster(MF.avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/MFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("MF.avoiddegrad.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario avoid deforestation
#MF.avoiddefor <- sum(UPF.avoiddefor, DPF2020, SF2010.mature, na.rm = T)
#MF.avoiddefor[MF.avoiddefor>1] <- 1
###cheking
##MF.avoiddefor
##plot(MF.avoiddefor)
#
##saving
#writeRaster(MF.avoiddefor, "rasters/PGM/input/MF2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#MF.avoiddefor.px <- focal(MF.avoiddefor, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##MF.avoiddefor.px
##anyNA(MF.avoiddefor.px[])
##plot(MF.avoiddefor.px)
#
#names(MF.avoiddefor.px)<-"MFpx"
#MF.avoiddefor.px[is.nan(MF.avoiddefor.px)] <- 0
#MF.avoiddefor.px <- mask(MF.avoiddefor.px, pgm.shp)
#
##saving
#writeRaster(MF.avoiddefor.px, "rasters/PGM/2020_avoiddeforest/MFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("MF.avoiddefor.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario avoid deforestation upf only
#MF.avoiddefor2 <- sum(UPF.avoiddefor2, DPF2020, SF2020.mature, na.rm = T)
#MF.avoiddefor2[MF.avoiddefor2>1] <- 1
###cheking
##MF.avoiddefor2
##plot(MF.avoiddefor2)
#
##saving
#writeRaster(MF.avoiddefor2, "rasters/PGM/input/MF2020_avoiddeforest2.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#MF.avoiddefor2.px <- focal(MF.avoiddefor2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##MF.avoiddefor2.px
##anyNA(MF.avoiddefor2.px[])
##plot(MF.avoiddefor2.px)
#
#names(MF.avoiddefor2.px)<-"MFpx"
#MF.avoiddefor2.px[is.nan(MF.avoiddefor2.px)] <- 0
#MF.avoiddefor2.px <- mask(MF.avoiddefor2.px, pgm.shp)
#
##saving
#writeRaster(MF.avoiddefor2.px, "rasters/PGM/2020_avoiddeforest2/MFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("MF.avoiddefor2.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario avoid both
#MF.avoidboth <- sum(MF.avoiddegrad, MF.avoiddefor, na.rm = T)
#MF.avoidboth[MF.avoidboth>1] <- 1
###cheking
##MF.avoidboth
##plot(MF.avoidboth)
#
##saving
#writeRaster(MF.avoidboth, "rasters/PGM/input/MF2020_avoidboth.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#MF.avoidboth.px <- focal(MF.avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##MF.avoidboth.px
##anyNA(MF.avoidboth.px[])
##plot(MF.avoidboth.px)
#
#names(MF.avoidboth.px)<-"MFpx"
#MF.avoidboth.px[is.nan(MF.avoidboth.px)] <- 0
#MF.avoidboth.px <- mask(MF.avoidboth.px, pgm.shp)
#
##saving
#writeRaster(MF.avoidboth.px, "rasters/PGM/2020_avoidboth/MFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("MF.avoidboth.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario restoration WITHOUT avoid
#SFAge2020.restore10.mature <- SFAge2020.restore10
#SFAge2020.restore10.mature[SFAge2020.restore10.mature < 10] <- 0
#SFAge2020.restore10.mature[SFAge2020.restore10.mature >= 10] <- 1
#
##saving
#writeRaster(SFAge2020.restore10.mature, "rasters/PGM/input/SFAge2020_restor_wo_avoid-mature.tif", format="GTiff", overwrite=T)
#
#
#MF.restore10.a <- sum(UPF2020, DPF2020, SFAge2020.restore10.mature, na.rm = T)
#MF.restore10.a[MF.restore10.a>1] <- 1
###cheking
##MF.restore10.a
##plot(MF.restore10.a)
#
##saving
#writeRaster(MF.restore10.a, "rasters/PGM/input/MF2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#MF.restore10.a.px <- focal(MF.restore10.a, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##MF.restore10.a.px
##anyNA(MF.restore10.a.px[])
##plot(MF.restore10.a.px)
#
#names(MF.restore10.a.px)<-"MFpx"
#MF.restore10.a.px[is.nan(MF.restore10.a.px)] <- 0
#MF.restore10.a.px <- mask(MF.restore10.a.px, pgm.shp)
#
##saving
#writeRaster(MF.restore10.a.px, "rasters/PGM/2020_restor_wo_avoid/MFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("MF.restore10.a.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario restoration AND avoid deforestation
#SFAge2010.restore10.mature <- SFAge2010.restore10
#SFAge2010.restore10.mature[SFAge2010.restore10.mature < 10] <- 0
#SFAge2010.restore10.mature[SFAge2010.restore10.mature >= 10] <- 1
#
##saving
#writeRaster(SFAge2010.restore10.mature, "rasters/PGM/input/SFAge2020_restor_n_avoid_deforest-mature.tif", format="GTiff", overwrite=T)
#
#
#MF.restore10.b <- sum(UPF.avoiddefor, DPF2020, SFAge2010.restore10.mature, na.rm = T)
#MF.restore10.b[MF.restore10.b>1] <- 1
###cheking
##MF.restore10.b
##plot(MF.restore10.b)
#
##saving
#writeRaster(MF.restore10.b, "rasters/PGM/input/MF2020_restor_n_avoid_deforest.tif", format="GTiff", overwrite=T)
#
#
##mean upf cover in pixel scale (200m)
#MF.restore10.b.px <- focal(MF.restore10.b, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##MF.restore10.b.px
##anyNA(MF.restore10.b.px[])
##plot(MF.restore10.b.px)
#
#names(MF.restore10.b.px)<-"MFpx"
#MF.restore10.b.px[is.nan(MF.restore10.b.px)] <- 0
#MF.restore10.b.px <- mask(MF.restore10.b.px, pgm.shp)
#
##saving
#writeRaster(MF.restore10.b.px, "rasters/PGM/2020_restor_n_avoid_deforest/MFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("MF.restore10.b.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario restoration AND avoid both
#MF.restore10.c <- sum(UPF.avoidboth, DPF2010, SFAge2010.restore10.mature, na.rm = T)
#MF.restore10.c[MF.restore10.c>1] <- 1
###cheking
##MF.restore10.c
##plot(MF.restore10.c)
#
##saving
#writeRaster(MF.restore10.c, "rasters/PGM/input/MF2020_restor_n_avoid_both.tif", format="GTiff", overwrite=T)
#
#
## mean upf cover in pixel scale (200m)
#MF.restore10.c.px <- focal(MF.restore10.c, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##MF.restore10.c.px
##anyNA(MF.restore10.c.px[])
##plot(MF.restore10.c.px)
#
#names(MF.restore10.c.px)<-"MFpx"
#MF.restore10.c.px[is.nan(MF.restore10.c.px)] <- 0
#MF.restore10.c.px <- mask(MF.restore10.c.px, pgm.shp)
#
##saving
#writeRaster(MF.restore10.c.px, "rasters/PGM/2020_restor_n_avoid_both/MFpx.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("MF.restore10.c.px")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
########################################################################################################################.#######################################################################################################################.
#
## [edgedist] distance to forest edge
## this variable is the euclidean distance between mature forest 
## to the nearest cell that is not MF
#
## scenario 2010
#inv.MF2010 <- MF2010
#inv.MF2010[inv.MF2010==1]<-NA
##cheking
##inv.MF2010
##plot(inv.MF2010)
#
#edge.dist.2010 <- distance(inv.MF2010)
###cheking
##edge.dist.2010
##anyNA(edge.dist.2010[])
##plot(edge.dist.2010)
#
#names(edge.dist.2010)<-"edgedist"
#edge.dist.2010[is.nan(edge.dist.2010)] <- 0
#edge.dist.2010 <- mask(edge.dist.2010, pgm.shp)
#
##saving
#writeRaster(edge.dist.2010, "rasters/PGM/2010_real/edgedist.tif", format="GTiff", overwrite=T)
##writeRaster(edge.dist.2010, "rasters/PGM/2020_avoidboth/edgedist.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("inv.MF2010")])
#gc()
#
#
#
##
#
#
## scenario 2020
#inv.MF2020 <- MF2020
#inv.MF2020[inv.MF2020==1]<-NA
##cheking
##inv.MF2020
##plot(inv.MF2020)
#
#edge.dist.2020 <- distance(inv.MF2020)
###cheking
##edge.dist.2020
##anyNA(edge.dist.2020[])
##plot(edge.dist.2020)
#
#names(edge.dist.2020)<-"edgedist"
#edge.dist.2020[is.nan(edge.dist.2020)] <- 0
#edge.dist.2020 <- mask(edge.dist.2020, pgm.shp)
#
##saving
#writeRaster(edge.dist.2020, "rasters/PGM/2020_real/edgedist.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("inv.MF2020")])
#gc()
#
#
#
##
#
#
## scenario avoid degradation
#inv.MF.avoiddegrad <- MF.avoiddegrad
#inv.MF.avoiddegrad[inv.MF.avoiddegrad==1]<-NA
##cheking
##inv.MF.avoiddegrad
##plot(inv.MF.avoiddegrad)
#
#edge.dist.avoiddegrad <- distance(inv.MF2020)
###cheking
##edge.dist.avoiddegrad
##anyNA(edge.dist.avoiddegrad[])
##plot(edge.dist.avoiddegrad)
#
#names(edge.dist.avoiddegrad)<-"edgedist"
#edge.dist.avoiddegrad[is.nan(edge.dist.avoiddegrad)] <- 0
#edge.dist.avoiddegrad <- mask(edge.dist.avoiddegrad, pgm.shp)
#
##saving
#writeRaster(edge.dist.avoiddegrad, "rasters/PGM/2020_avoiddegrad/edgedist.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("inv.MF.avoiddegrad")])
#gc()
#
#
#
##
#
#
## scenario avoid deforestation
#inv.MF.avoiddefor <- MF.avoiddefor
#inv.MF.avoiddefor[inv.MF.avoiddefor==1]<-NA
##cheking
##inv.MF.avoiddefor
##plot(inv.MF.avoiddefor)
#
#edge.dist.avoiddefor <- distance(inv.MF.avoiddefor)
###cheking
##edge.dist.avoiddefor
##anyNA(edge.dist.avoiddefor[])
##plot(edge.dist.avoiddefor)
#
#names(edge.dist.avoiddefor)<-"edgedist"
#edge.dist.avoiddefor[is.nan(edge.dist.avoiddefor)] <- 0
#edge.dist.avoiddefor <- mask(edge.dist.avoiddefor, pgm.shp)
#
##saving
#writeRaster(edge.dist.avoiddefor, "rasters/PGM/2020_avoiddeforest/edgedist.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("inv.MF.avoiddefor")])
#gc()
#
#
#
##
#
#
## scenario avoid deforestation upf only
#inv.MF.avoiddefor2 <- MF.avoiddefor2
#inv.MF.avoiddefor2[inv.MF.avoiddefor2==1]<-NA
##cheking
##inv.MF.avoiddefor2
##plot(inv.MF.avoiddefor2)
#
#edge.dist.avoiddefor2 <- distance(inv.MF.avoiddefor2)
###cheking
##edge.dist.avoiddefor2
##anyNA(edge.dist.avoiddefor2[])
##plot(edge.dist.avoiddefor2)
#
#names(edge.dist.avoiddefor2)<-"edgedist"
#edge.dist.avoiddefor2[is.nan(edge.dist.avoiddefor2)] <- 0
#edge.dist.avoiddefor2 <- mask(edge.dist.avoiddefor2, pgm.shp)
#
##saving
#writeRaster(edge.dist.avoiddefor, "rasters/PGM/2020_avoiddeforest2/edgedist.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("inv.MF.avoiddefor")])
#gc()
#
#
#
##
#
#
## scenario avoid both
#inv.MF.avoidboth <- MF.avoidboth
#inv.MF.avoidboth[inv.MF.avoidboth==1]<-NA
##cheking
##inv.MF.avoidboth
##plot(inv.MF.avoidboth)
#
#edge.dist.avoidboth <- distance(inv.MF.avoidboth)
###cheking
##edge.dist.avoidboth
##anyNA(edge.dist.avoidboth[])
##plot(edge.dist.avoidboth)
#
#names(edge.dist.avoidboth)<-"edgedist"
#edge.dist.avoidboth[is.nan(edge.dist.avoidboth)] <- 0
#edge.dist.avoidboth <- mask(edge.dist.avoidboth, pgm.shp)
#
##saving
#writeRaster(edge.dist.avoidboth, "rasters/PGM/2020_avoidboth/edgedist.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("inv.MF.avoidboth")])
#gc()
#
#
#
##
#
#
## scenario restoration WITHOUT avoid
#inv.MF.restore10.a <- MF.restore10.a
#inv.MF.restore10.a[inv.MF.restore10.a==1]<-NA
##cheking
##inv.MF.avoidboth
##plot(inv.MF.avoidboth)
#
#edge.dist.restore10.a <- distance(inv.MF.restore10.a)
###cheking
##edge.dist.restore10.a
##anyNA(edge.dist.restore10.a[])
##plot(edge.dist.restore10.a)
#
#names(edge.dist.restore10.a)<-"edgedist"
#edge.dist.restore10.a[is.nan(edge.dist.restore10.a)] <- 0
#edge.dist.restore10.a <- mask(edge.dist.restore10.a, pgm.shp)
#
##saving
#writeRaster(edge.dist.restore10.a, "rasters/PGM/2020_restor_wo_avoid/edgedist.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("inv.MF.restore10.a")])
#gc()
#
#
#
##
#
#
## scenario restoration AND avoid deforest
#inv.MF.restore10.b <- MF.restore10.b
#inv.MF.restore10.b[inv.MF.restore10.b==1]<-NA
##cheking
##inv.MF.avoidboth
##plot(inv.MF.avoidboth)
#
#edge.dist.restore10.b <- distance(inv.MF.restore10.b)
###cheking
##edge.dist.restore10.b
##anyNA(edge.dist.restore10.b[])
##plot(edge.dist.restore10.b)
#
#names(edge.dist.restore10.b)<-"edgedist"
#edge.dist.restore10.b[is.nan(edge.dist.restore10.b)] <- 0
#edge.dist.restore10.b <- mask(edge.dist.restore10.b, pgm.shp)
#
##saving
#writeRaster(edge.dist.restore10.b, "rasters/PGM/2020_restor_n_avoid_deforest/edgedist.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("inv.MF.restore10.b")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario restoration AND avoid both
#inv.MF.restore10.c <- MF.restore10.c
#inv.MF.restore10.c[inv.MF.restore10.c==1]<-NA
##cheking
##inv.MF.avoidboth
##plot(inv.MF.avoidboth)
#
#edge.dist.restore10.c <- distance(inv.MF.restore10.c)
###cheking
##edge.dist.restore10.b
##anyNA(edge.dist.restore10.b[])
##plot(edge.dist.restore10.b)
#
#names(edge.dist.restore10.c)<-"edgedist"
#edge.dist.restore10.c[is.nan(edge.dist.restore10.c)] <- 0
#edge.dist.restore10.c <- mask(edge.dist.restore10.c, pgm.shp)
#
##saving
#writeRaster(edge.dist.restore10.c, "rasters/PGM/2020_restor_n_avoid_both/edgedist.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("inv.MF.restore10.c")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
########################################################################################################################.#######################################################################################################################.
#
## [edge] forest edge -- pixel: 5x5 (200m); and landscape: 35x35 (1000m)
## this variable is the mean of edge area based on mature forest
#
## scenario 2010
## marking edge and core areas
#edge2010 <- edge.dist.2010
#edge2010[edge2010>300] <- 0
#edge2010[edge2010!=0] <- 1
#
##saving
#writeRaster(edge2010, "rasters/PGM/input/edge2010_real.tif", format="GTiff", overwrite=T)
#
#
##alternative method
##edge2010 <- focal(MF2010, matrix(1,ncol=3,nrow=3), fun = mean, na.rm = T, pad = T)
##edge2010 <- edge2010 + MF2010
##edge2010[edge2010 == 2] <- 0                  # core area
##edge2010[edge2010 > 0 & edge2010 < 2] <- 1    # edge
#
##cheking
##edge2010
##unique(edge2010[])
##plot(edge2010)
#
## mean edge in pixel scale (200m)
#edge2010.px <- focal(edge2010, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##edge2010.px
##anyNA(edge2010.px[])
##plot(edge2010.px)
#
#names(edge2010.px)<-"edgepx"
#edge2010.px[is.nan(edge2010.px)] <- 0
#edge2010.px <- mask(edge2010.px, pgm.shp)
#
##saving
#writeRaster(edge2010.px, "rasters/PGM/2010_real/edgepx.tif", format="GTiff", overwrite=T)
##
#
## mean sf cover in landscape scale (1000m)
#edge2010.ls <- focal(edge2010, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##edge2010.ls
##anyNA(edge2010.ls[])
##plot(edge2010.ls)
#
#names(edge2010.ls)<-"edgels"
#edge2010.ls[is.nan(edge2010.ls)] <- 0
#edge2010.ls <- mask(edge2010.ls, pgm.shp)
#
##saving
#writeRaster(edge2010.px, "rasters/PGM/2010_real/edgels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("edge2010.px", "edge2010.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario 2020
## marking edge and core areas
#edge2020 <- edge.dist.2020
#edge2020[edge2020>300] <- 0
#edge2020[edge2020!=0] <- 1
#
##saving
#writeRaster(edge2020, "rasters/PGM/input/edge2020_real.tif", format="GTiff", overwrite=T)
#
#
##cheking
##edge2020
##unique(edge2020[])
##plot(edge2020)
#
## mean edge in pixel scale (200m)
#edge2020.px <- focal(edge2020, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##edge2020.px
##anyNA(edge2020.px[])
##plot(edge2020.px)
#
#names(edge2020.px)<-"edgepx"
#edge2020.px[is.nan(edge2020.px)] <- 0
#edge2020.px <- mask(edge2020.px, pgm.shp)
#
##saving
#writeRaster(edge2020.px, "rasters/PGM/2020_real/edgepx.tif", format="GTiff", overwrite=T)
##
#
## mean sf cover in landscape scale (1000m)
#edge2020.ls <- focal(edge2020, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##edge2020.ls
##anyNA(edge2020.ls[])
##plot(edge2020.ls)
#
#names(edge2020.ls)<-"edgels"
#edge2020.ls[is.nan(edge2020.ls)] <- 0
#edge2020.ls <- mask(edge2020.ls, pgm.shp)
#
##saving
#writeRaster(edge2020.ls, "rasters/PGM/2020_real/edgels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("edge2020.px", "edge2020.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario avoid degradation
## marking edge and core areas
#edge.avoiddegrad <- edge.dist.avoiddegrad
#edge.avoiddegrad[edge.avoiddegrad>300] <- 0
#edge.avoiddegrad[edge.avoiddegrad!=0] <- 1
#
##saving
#writeRaster(edge.avoiddegrad, "rasters/PGM/input/edge2020_avoiddegrad.tif", format="GTiff", overwrite=T)
#
#
##cheking
##edge.avoiddegrad
##unique(edge.avoiddegrad[])
##plot(edge.avoiddegrad)
#
## mean edge in pixel scale (200m)
#edge.avoiddegrad.px <- focal(edge.avoiddegrad, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##edge.avoiddegrad.px
##anyNA(edge.avoiddegrad.px[])
##plot(edge.avoiddegrad.px)
#
#names(edge.avoiddegrad.px)<-"edgepx"
#edge.avoiddegrad.px[is.nan(edge.avoiddegrad.px)] <- 0
#edge.avoiddegrad.px <- mask(edge.avoiddegrad.px, pgm.shp)
#
##saving
#writeRaster(edge.avoiddegrad.px, "rasters/PGM/2020_avoiddegrad/edgepx.tif", format="GTiff", overwrite=T)
##
#
## mean sf cover in landscape scale (1000m)
#edge.avoiddegrad.ls <- focal(edge.avoiddegrad, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##edge.avoiddegrad.ls
##anyNA(edge.avoiddegrad.ls[])
##plot(edge.avoiddegrad.ls)
#
#names(edge.avoiddegrad.ls)<-"edgels"
#edge.avoiddegrad.ls[is.nan(edge.avoiddegrad.ls)] <- 0
#edge.avoiddegrad.ls <- mask(edge.avoiddegrad.ls, pgm.shp)
#
##saving
#writeRaster(edge.avoiddegrad.ls, "rasters/PGM/2020_avoiddegrad/edgels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("edge.avoiddegrad.px", "edge.avoiddegrad.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario avoid deforestation
## marking edge and core areas
#edge.avoiddefor <- edge.dist.avoiddefor
#edge.avoiddefor[edge.avoiddefor>300] <- 0
#edge.avoiddefor[edge.avoiddefor!=0] <- 1
#
##saving
#writeRaster(edge.avoiddefor, "rasters/PGM/input/edge2020_avoiddeforest.tif", format="GTiff", overwrite=T)
#
#
##cheking
##edge.avoiddefor
##unique(edge.avoiddefor[])
##plot(edge.avoiddefor)
#
## mean edge in pixel scale (200m)
#edge.avoiddefor.px <- focal(edge.avoiddefor, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##edge.avoiddefor.px
##anyNA(edge.avoiddefor.px[])
##plot(edge.avoiddefor.px)
#
#names(edge.avoiddefor.px)<-"edgepx"
#edge.avoiddefor.px[is.nan(edge.avoiddefor.px)] <- 0
#edge.avoiddefor.px <- mask(edge.avoiddefor.px, pgm.shp)
#
##saving
#writeRaster(edge.avoiddefor.px, "rasters/PGM/2020_avoiddeforest/edgepx.tif", format="GTiff", overwrite=T)
##
#
## mean sf cover in landscape scale (1000m)
#edge.avoiddefor.ls <- focal(edge.avoiddefor, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##edge.avoiddefor.ls
##anyNA(edge.avoiddefor.ls[])
##plot(edge.avoiddefor.ls)
#
#names(edge.avoiddefor.ls)<-"edgels"
#edge.avoiddefor.ls[is.nan(edge.avoiddefor.ls)] <- 0
#edge.avoiddefor.ls <- mask(edge.avoiddefor.ls, pgm.shp)
#
##saving
#writeRaster(edge.avoiddefor.ls, "rasters/PGM/2020_avoiddeforest/edgels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("edge.avoiddefor.px", "edge.avoiddefor.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario avoid deforestation upf only
## marking edge and core areas
#edge.avoiddefor2 <- edge.dist.avoiddefor2
#edge.avoiddefor2[edge.avoiddefor2>300] <- 0
#edge.avoiddefor2[edge.avoiddefor2!=0] <- 1
#
##saving
#writeRaster(edge.avoiddefor2, "rasters/PGM/input/edge2020_avoiddeforest2.tif", format="GTiff", overwrite=T)
#
#
##cheking
##edge.avoiddefor2
##unique(edge.avoiddefor2[])
##plot(edge.avoiddefor2)
#
## mean edge in pixel scale (200m)
#edge.avoiddefor2.px <- focal(edge.avoiddefor2, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##edge.avoiddefor2.px
##anyNA(edge.avoiddefor2.px[])
##plot(edge.avoiddefor2.px)
#
#names(edge.avoiddefor2.px)<-"edgepx"
#edge.avoiddefor2.px[is.nan(edge.avoiddefor2.px)] <- 0
#edge.avoiddefor2.px <- mask(edge.avoiddefor2.px, pgm.shp)
#
##saving
#writeRaster(edge.avoiddefor2.px, "rasters/PGM/2020_avoiddeforest2/edgepx.tif", format="GTiff", overwrite=T)
##
#
## mean sf cover in landscape scale (1000m)
#edge.avoiddefor2.ls <- focal(edge.avoiddefor2, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##edge.avoiddefor2.ls
##anyNA(edge.avoiddefor2.ls[])
##plot(edge.avoiddefor2.ls)
#
#names(edge.avoiddefor2.ls)<-"edgels"
#edge.avoiddefor2.ls[is.nan(edge.avoiddefor2.ls)] <- 0
#edge.avoiddefor2.ls <- mask(edge.avoiddefor2.ls, pgm.shp)
#
##saving
#writeRaster(edge.avoiddefor2.ls, "rasters/PGM/2020_avoiddeforest2/edgels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("edge.avoiddefor2.px", "edge.avoiddefor2.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario avoid both
## marking edge and core areas
#edge.avoidboth <- edge.dist.avoidboth
#edge.avoidboth[edge.avoidboth>300] <- 0
#edge.avoidboth[edge.avoidboth!=0] <- 1
#
##saving
#writeRaster(edge.avoidboth, "rasters/PGM/input/edge2020_avoidboth.tif", format="GTiff", overwrite=T)
#
#
##cheking
##edge.avoidboth
##unique(edge.avoidboth[])
##plot(edge.avoidboth)
#
## mean edge in pixel scale (200m)
#edge.avoidboth.px <- focal(edge.avoidboth, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##edge.avoidboth.px
##anyNA(edge.avoidboth.px[])
##plot(edge.avoidboth.px)
#
#names(edge.avoidboth.px)<-"edgepx"
#edge.avoidboth.px[is.nan(edge.avoidboth.px)] <- 0
#edge.avoidboth.px <- mask(edge.avoidboth.px, pgm.shp)
#
##saving
#writeRaster(edge.avoidboth.px, "rasters/PGM/2020_avoidboth/edgepx.tif", format="GTiff", overwrite=T)
##
#
## mean sf cover in landscape scale (1000m)
#edge.avoidboth.ls <- focal(edge.avoidboth, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##edge.avoidboth.ls
##anyNA(edge.avoidboth.ls[])
##plot(edge.avoidboth.ls)
#
#names(edge.avoidboth.ls)<-"edgels"
#edge.avoidboth.ls[is.nan(edge.avoidboth.ls)] <- 0
#edge.avoidboth.ls <- mask(edge.avoidboth.ls, pgm.shp)
#
##saving
#writeRaster(edge.avoidboth.ls, "rasters/PGM/2020_avoidboth/edgels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("edge.avoidboth.px", "edge.avoidboth.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario restoration WITHOUT avoid
## marking edge and core areas
#edge.restore10.a <- edge.dist.edge.restore10.a
#edge.restore10.a[edge.restore10.a>300] <- 0
#edge.restore10.a[edge.restore10.a!=0] <- 1
#
##saving
#writeRaster(edge.restore10.a, "rasters/PGM/input/edge2020_restor_wo_avoid.tif", format="GTiff", overwrite=T)
#
#
##cheking
##edge.restore10.a
##unique(edge.restore10.a[])
##plot(edge.restore10.a)
#
## mean edge in pixel scale (200m)
#edge.restore10.a.px <- focal(edge.restore10.a, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##edge.restore10.a.px
##anyNA(edge.restore10.a.px[])
##plot(edge.restore10.a.px)
#
#names(edge.restore10.a.px)<-"edgepx"
#edge.restore10.a.px[is.nan(edge.restore10.a.px)] <- 0
#edge.restore10.a.px <- mask(edge.restore10.a.px, pgm.shp)
#
##saving
#writeRaster(edge.restore10.a.px, "rasters/PGM/2020_restor_wo_avoid/edgepx.tif", format="GTiff", overwrite=T)
##
#
## mean sf cover in landscape scale (1000m)
#edge.restore10.a.ls <- focal(edge.restore10.a, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##edge.restore10.a.ls
##anyNA(edge.restore10.a.ls[])
##plot(edge.restore10.a.ls)
#
#names(edge.restore10.a.ls)<-"edgels"
#edge.restore10.a.ls[is.nan(edge.restore10.a.ls)] <- 0
#edge.restore10.a.ls <- mask(edge.restore10.a.ls, pgm.shp)
#
##saving
#writeRaster(edge.restore10.a.ls, "rasters/PGM/2020_restor_wo_avoid/edgels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("edge.restore10.a.px", "edge.restore10.a.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario restoration AND avoid deforest
## marking edge and core areas
#edge.restore10.b <- edge.dist.edge.restore10.b
#edge.restore10.b[edge.restore10.b>300] <- 0
#edge.restore10.b[edge.restore10.b!=0] <- 1
#
##saving
#writeRaster(edge.restore10.b, "rasters/PGM/input/edge2020_restor_n_avoid_deforest.tif", format="GTiff", overwrite=T)
#
#
##cheking
##edge.restore10.b
##unique(edge.restore10.b[])
##plot(edge.restore10.b)
#
## mean edge in pixel scale (200m)
#edge.restore10.b.px <- focal(edge.restore10.b, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##edge.restore10.b.px
##anyNA(edge.restore10.b.px[])
##plot(edge.restore10.b.px)
#
#names(edge.restore10.b.px)<-"edgepx"
#edge.restore10.b.px[is.nan(edge.restore10.b.px)] <- 0
#edge.restore10.b.px <- mask(edge.restore10.b.px, pgm.shp)
#
##saving
#writeRaster(edge.restore10.b.px, "rasters/PGM/2020_restor_n_avoid_deforest/edgepx.tif", format="GTiff", overwrite=T)
##
#
## mean sf cover in landscape scale (1000m)
#edge.restore10.b.ls <- focal(edge.restore10.b, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##edge.restore10.b.ls
##anyNA(edge.restore10.b.ls[])
##plot(edge.restore10.b.ls)
#
#names(edge.restore10.b.ls)<-"edgels"
#edge.restore10.b.ls[is.nan(edge.restore10.b.ls)] <- 0
#edge.restore10.b.ls <- mask(edge.restore10.b.ls, pgm.shp)
#
##saving
#writeRaster(edge.restore10.b.ls, "rasters/PGM/2020_restor_n_avoid_deforest/edgels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("edge.restore10.b.px", "edge.restore10.b.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
## scenario restoration AND avoid both
## marking edge and core areas
#edge.restore10.c <- edge.dist.edge.restore10.c
#edge.restore10.c[edge.restore10.c>300] <- 0
#edge.restore10.c[edge.restore10.c!=0] <- 1
#
##saving
#writeRaster(edge.restore10.c, "rasters/PGM/input/edge2020_restor_n_avoid_both.tif", format="GTiff", overwrite=T)
#
#
##cheking
##edge.restore10.c
##unique(edge.restore10.c[])
##plot(edge.restore10.c)
#
## mean edge in pixel scale (200m)
#edge.restore10.c.px <- focal(edge.restore10.c, matrix(1,ncol=3,nrow=3), fun=mean, na.rm=T)
###cheking
##edge.restore10.c.px
##anyNA(edge.restore10.c.px[])
##plot(edge.restore10.c.px)
#
#names(edge.restore10.c.px)<-"edgepx"
#edge.restore10.c.px[is.nan(edge.restore10.c.px)] <- 0
#edge.restore10.c.px <- mask(edge.restore10.c.px, pgm.shp)
#
##saving
#writeRaster(edge.restore10.c.px, "rasters/PGM/2020_restor_n_avoid_both/edgepx.tif", format="GTiff", overwrite=T)
##
#
## mean sf cover in landscape scale (1000m)
#edge.restore10.c.ls <- focal(edge.restore10.c, matrix(1,ncol=11,nrow=11), fun=mean, na.rm=T)
###cheking
##edge.restore10.c.ls
##anyNA(edge.restore10.c.ls[])
##plot(edge.restore10.c.ls)
#
#names(edge.restore10.c.ls)<-"edgels"
#edge.restore10.c.ls[is.nan(edge.restore10.c.ls)] <- 0
#edge.restore10.c.ls <- mask(edge.restore10.c.ls, pgm.shp)
#
##saving
#writeRaster(edge.restore10.c.ls, "rasters/PGM/2020_restor_n_avoid_both/edgels.tif", format="GTiff", overwrite=T)
#
#rm(list=ls()[ls() %in% c("edge.restore10.c.px", "edge.restore10.c.ls")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
########################################################################################################################.#######################################################################################################################.
#
## [meantemp] annual average temperature from nasa earth observation
#
## download and save global data
##urls <- c("https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755469&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-01"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755470&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-02"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755471&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-03"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755472&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-04"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755473&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-05"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755474&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-06"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755475&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-07"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755476&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-08"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755477&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-09"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755478&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-10"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755479&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-11"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1755480&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-12"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1784090&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-01"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1785058&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-02"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1785890&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-03"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1786979&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-04"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1794500&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-05"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1795357&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-06"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1796358&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-07"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1797155&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-08"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1799174&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-09"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1799941&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-10"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1800680&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-11"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1801650&cs=rgb&format=TIFF&width=3600&height=1800") #"2020-12"
##
##dir.create("rasters/PGM/input/climate")
#
#
## scenario 2010
#temp.list <- list.files("rasters/PGM/raw/climate", "LSTD", full.names = T, recursive = T)
#
#temp2010.list <- grep("2010", temp.list, value = T)
#temp2010 <- stack(temp2010.list)
##plot(temp2010)
#
#meantemp2010 <- mean(temp2010, na.rm=T)
#pgm.meantemp2010 <- crop(meantemp2010, extent(pgm.lulc.2010.forest.class))
##plot(pgm.meantemp2010)
##plot(pgm.shp, add=T)
#
#pgm.meantemp2010 <- resample(pgm.meantemp2010, pgm.lulc.2010.forest.class, method='bilinear')
#pgm.meantemp2010 <- mask(pgm.meantemp2010, pgm.shp)
##plot(pgm.meantemp2010)
#
##saving
#writeRaster(pgm.meantemp2010, "rasters/PGM/2010_real/meantemps.tif", format="GTiff", overwrite=T)
##
#
## scenario 2020
#temp2020.list <- grep("2020", temp.list, value = T)
#temp2020 <- stack(temp2020.list)
##plot(temp2020)
#
#meantemp2020 <- mean(temp2020, na.rm=T)
#pgm.meantemp2020 <- crop(meantemp2020, extent(pgm.lulc.2020.forest.class))
##plot(pgm.meantemp2020)
##plot(pgm.shp, add=T)
#
#pgm.meantemp2020 <- resample(pgm.meantemp2020, pgm.lulc.2020.forest.class, method='bilinear')
#pgm.meantemp2020 <- mask(pgm.meantemp2020, pgm.shp)
##plot(pgm.meantemp2020)
#
##saving
#writeRaster(pgm.meantemp2020, "rasters/PGM/2020_real/meantemps.tif", format="GTiff", overwrite=T)
#writeRaster(pgm.meantemp2020, "rasters/PGM/2020_avoiddeforest/meantemps.tif", format="GTiff", overwrite=T)
#writeRaster(pgm.meantemp2020, "rasters/PGM/2020_avoiddeforest2/meantemps.tif", format="GTiff", overwrite=T)
#writeRaster(pgm.meantemp2020, "rasters/PGM/2020_avoiddegrad/meantemps.tif", format="GTiff", overwrite=T)
#writeRaster(pgm.meantemp2020, "rasters/PGM/2020_avoidboth/meantemps.tif", format="GTiff", overwrite=T)
#writeRaster(pgm.meantemp2020, "rasters/PGM/2020_restor_wo_avoid/meantemps.tif", format="GTiff", overwrite=T)
#writeRaster(pgm.meantemp2020, "rasters/PGM/2020_restor_n_avoid_deforest/meantemps.tif", format="GTiff", overwrite=T)
#writeRaster(pgm.meantemp2020, "rasters/PGM/2020_restor_n_avoid_both/meantemps.tif", format="GTiff", overwrite=T)
##
#
#rm(list=ls()[ls() %in% c("temp.list", "temp2010.list", "meantemp2010", "temp2020.list", "meantemp2020")]) #keeping only raster stack
#gc()
#
#
#
##
#
#
########################################################################################################################.#######################################################################################################################.
#
## [meanprecip] annual average precipitation from nasa earth observation
#
## download and save global data
##urls <- c("https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843747&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-01"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843749&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-02"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843751&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-03"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843755&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-04"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843735&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-05"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843745&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-06"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843753&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-07"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843759&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-08"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843761&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-09"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843763&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-10"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843765&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-11"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843769&cs=rgb&format=TIFF&width=3600&height=1800", #"2010-12"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843983&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-01"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843985&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-02"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843987&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-03"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843989&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-04"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843991&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-05"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843993&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-06"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843995&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-07"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843997&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-08"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1843999&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-09"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1844001&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-10"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1844003&cs=rgb&format=TIFF&width=3600&height=1800", #"2020-11"
##          "https://neo.gsfc.nasa.gov/servlet/RenderData?si=1844005&cs=rgb&format=TIFF&width=3600&height=1800") #"2020-12"
##
#
#
## scenario 2010
#precip.list <- list.files("rasters/PGM/raw/climate", "GPM", full.names = T, recursive = T)
#
#precip2010.list <- grep("2010", precip.list, value = T)
#precip2010 <- stack(precip2010.list)
##plot(precip2010)
#rm(precip2010.list)
#
#meanprecip2010 <- mean(precip2010, na.rm=T)
#pgm.meanprecip2010 <- crop(meanprecip2010, extent(pgm.lulc.2010.forest.class))
##plot(pgm.meanprecip2010)
##plot(pgm.shp, add=T)
#
#pgm.meanprecip2010 <- resample(pgm.meanprecip2010, pgm.lulc.2010.forest.class, method='bilinear')
#pgm.meanprecip2010 <- mask(pgm.meanprecip2010, pgm.shp)
##plot(pgm.meanprecip2010)
#
##saving
#writeRaster(pgm.meanprecip2010, "rasters/PGM/2010_real/meanprecips.tif", format="GTiff", overwrite=T)
##
#
## scenario 2020
#precip2020.list <- grep("2020", precip.list, value = T)
#precip2020 <- stack(precip2020.list)
##plot(precip2020)
#
#
#meanprecip2020 <- mean(precip2020, na.rm=T)
#pgm.meanprecip2020 <- crop(meanprecip2020, extent(pgm.lulc.2020.forest.class))
##plot(pgm.meanprecip2020)
##plot(pgm.shp, add=T)
#
#pgm.meanprecip2020 <- resample(pgm.meanprecip2020, pgm.lulc.2020.forest.class, method='bilinear')
#pgm.meanprecip2020 <- mask(pgm.meanprecip2020, pgm.shp)
##plot(pgm.meanprecip2020)
#
##saving
#writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_real/meanprecips.tif", format="GTiff", overwrite=T)
#writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_avoiddeforest/meanprecips.tif", format="GTiff", overwrite=T)
#writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_avoiddeforest2/meanprecips.tif", format="GTiff", overwrite=T)
#writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_avoiddegrad/meanprecips.tif", format="GTiff", overwrite=T)
#writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_avoidboth/meanprecips.tif", format="GTiff", overwrite=T)
#writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_restor_wo_avoid/meanprecips.tif", format="GTiff", overwrite=T)
#writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_restor_n_avoid_deforest/meanprecips.tif", format="GTiff", overwrite=T)
#writeRaster(pgm.meanprecip2020, "rasters/PGM/2020_restor_n_avoid_both/meanprecips.tif", format="GTiff", overwrite=T)
##
#
#rm(list=ls()[ls() %in% c("precip.list", "precip2010.list", "meanprecip2010", "precip2020.list", "meanprecip2020")])
#gc()
#
#
#
##
#
#
########################################################################################################################.#######################################################################################################################.
## [elevation]
#elevation <- raster("rasters/PGM/raw/elevation_pgm.tif")
#elevation <- projectRaster(elevation, crs = std.proj)
#elevation <- resample(elevation, pgm.lulc.2010.forest.class, method='bilinear')
#values(elevation)[is.na(values(elevation))] = 0
#elevation <- mask(elevation, pgm.shp)
#
##saving
#writeRaster(elevation, "rasters/PGM/2010_real/elevation.tif", format="GTiff", overwrite=T)
#writeRaster(elevation, "rasters/PGM/2020_real/elevation.tif", format="GTiff", overwrite=T)
#writeRaster(elevation, "rasters/PGM/2020_avoiddeforest/elevation.tif", format="GTiff", overwrite=T)
#writeRaster(elevation, "rasters/PGM/2020_avoiddeforest2/elevation.tif", format="GTiff", overwrite=T)
#writeRaster(elevation, "rasters/PGM/2020_avoiddegrad/elevation.tif", format="GTiff", overwrite=T)
#writeRaster(elevation, "rasters/PGM/2020_avoidboth/elevation.tif", format="GTiff", overwrite=T)
#writeRaster(elevation, "rasters/PGM/2020_restor_wo_avoid/elevation.tif", format="GTiff", overwrite=T)
#writeRaster(elevation, "rasters/PGM/2020_restor_n_avoid_deforest/elevation.tif", format="GTiff", overwrite=T)
#writeRaster(elevation, "rasters/PGM/2020_restor_n_avoid_both/elevation.tif", format="GTiff", overwrite=T)
#
##
#
## [distriver] distance to rivers
#pgm.river.raster <- rasterize(pgm.river, pgm.lulc, field = "buffer", background = 0)
##checking
##st_crs(pgm.river.raster)==st_crs(pgm.shp)
#
#inv.pgm.river <- pgm.river.raster
#inv.pgm.river[inv.pgm.river==0] <- NA
##cheking
##inv.pgm.river
##plot(inv.pgm.river)
#
#dist.river <- distance(inv.pgm.river)
#dist.river <- mask(dist.river, pgm.shp)
###cheking
##dist.river
##anyNA(dist.river[])
##plot(dist.river)
#
##saving
#writeRaster(dist.river, "rasters/PGM/2010_real/distriver.tif", format="GTiff", overwrite=T)
#writeRaster(dist.river, "rasters/PGM/2020_real/distriver.tif", format="GTiff", overwrite=T)
#writeRaster(dist.river, "rasters/PGM/2020_avoiddeforest/distriver.tif", format="GTiff", overwrite=T)
#writeRaster(dist.river, "rasters/PGM/2020_avoiddeforest2/distriver.tif", format="GTiff", overwrite=T)
#writeRaster(dist.river, "rasters/PGM/2020_avoiddegrad/distriver.tif", format="GTiff", overwrite=T)
#writeRaster(dist.river, "rasters/PGM/2020_avoidboth/distriver.tif", format="GTiff", overwrite=T)
#writeRaster(dist.river, "rasters/PGM/2020_restor_wo_avoid/distriver.tif", format="GTiff", overwrite=T)
#writeRaster(dist.river, "rasters/PGM/2020_restor_n_avoid_deforest/distriver.tif", format="GTiff", overwrite=T)
#writeRaster(dist.river, "rasters/PGM/2020_restor_n_avoid_both/distriver.tif", format="GTiff", overwrite=T)
#
#
## [distroad] distance to roads
#dist.road <- raster("rasters/PGM/raw/dist_road_pgm.tif")
#dist.road <- projectRaster(dist.road, crs = std.proj)
#dist.road <- resample(dist.road, pgm.lulc.2010.forest.class, method='bilinear')
#values(dist.road)[is.na(values(dist.road))] <- 0
#dist.road <- mask(dist.road, pgm.shp)
#
##saving
#writeRaster(dist.road, "rasters/PGM/2010_real/distroad.tif", format="GTiff", overwrite=T)
#writeRaster(dist.road, "rasters/PGM/2020_real/distroad.tif", format="GTiff", overwrite=T)
#writeRaster(dist.road, "rasters/PGM/2020_avoiddeforest/distroad.tif", format="GTiff", overwrite=T)
#writeRaster(dist.road, "rasters/PGM/2020_avoiddeforest2/distroad.tif", format="GTiff", overwrite=T)
#writeRaster(dist.road, "rasters/PGM/2020_avoiddegrad/distroad.tif", format="GTiff", overwrite=T)
#writeRaster(dist.road, "rasters/PGM/2020_avoidboth/distroad.tif", format="GTiff", overwrite=T)
#writeRaster(dist.road, "rasters/PGM/2020_restor_wo_avoid/distroad.tif", format="GTiff", overwrite=T)
#writeRaster(dist.road, "rasters/PGM/2020_restor_n_avoid_deforest/distroad.tif", format="GTiff", overwrite=T)
#writeRaster(dist.road, "rasters/PGM/2020_restor_n_avoid_both/distroad.tif", format="GTiff", overwrite=T)
#
##
#
## [distmarket] distance to municipality nucleus
#
#pgm.munic.nucleus <- data.frame(ID = "pgm", long = -47.35311, lat = -3.00249)
#
#pgm.munic.nucleus.coord <- SpatialPointsDataFrame(coords = pgm.munic.nucleus[,c("long","lat")], 
#                                                  data = pgm.munic.nucleus, 
#                                                  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
#
#pgm.munic.nucleus.coord <- spTransform(pgm.munic.nucleus.coord, crs(std.proj))
#
#distmarket <- distanceFromPoints(object = pgm.lulc, xy = pgm.munic.nucleus.coord)
#
#distmarket <- mask(distmarket, pgm.shp)
#
##saving
#writeRaster(distmarket, "rasters/PGM/2010_real/distmarket.tif", format="GTiff", overwrite=T)
#writeRaster(distmarket, "rasters/PGM/2020_real/distmarket.tif", format="GTiff", overwrite=T)
#writeRaster(distmarket, "rasters/PGM/2020_avoiddeforest/distmarket.tif", format="GTiff", overwrite=T)
#writeRaster(distmarket, "rasters/PGM/2020_avoiddeforest2/distmarket.tif", format="GTiff", overwrite=T)
#writeRaster(distmarket, "rasters/PGM/2020_avoiddegrad/distmarket.tif", format="GTiff", overwrite=T)
#writeRaster(distmarket, "rasters/PGM/2020_avoidboth/distmarket.tif", format="GTiff", overwrite=T)
#writeRaster(distmarket, "rasters/PGM/2020_restor_wo_avoid/distmarket.tif", format="GTiff", overwrite=T)
#writeRaster(distmarket, "rasters/PGM/2020_restor_n_avoid_deforest/distmarket.tif", format="GTiff", overwrite=T)
#writeRaster(distmarket, "rasters/PGM/2020_restor_n_avoid_both/distmarket.tif", format="GTiff", overwrite=T)
#
#
#
##
#
#
########################################################################################################################.
#
## Creating forest mask for conservation action costs
## each land cover category is assigned a value
#dir.create("rasters/PGM/all_forest_mask", recursive = T)
#
#
## scenario 2010
#UPF2010.mask <- UPF2010
#
#DPF2010.mask <- DPF2010
#
#SF2010.mask <- SF2010
#
#TF2010.mask <- sum(UPF2010.mask, DPF2010.mask, SF2010.mask, na.rm = T)
###cheking
##sort(unique(values(TF2010.mask)))
###handle overlaps
#TF2010.mask[TF2010.mask>1] <- 1
##plot(TF2010.mask)
#writeRaster(TF2010.mask, "rasters/PGM/all_forest_mask/PGM_2010_real.tif", format = "GTiff", overwrite = T)
#
#
##
#
#
## scenario 2020
#UPF2020.mask <- UPF2020
#
#DPF2020.mask <- DPF2020
#
#SF2020.mask <- SF2020
#
#TF2020.mask <- sum(UPF2020.mask, DPF2020.mask, SF2020.mask, na.rm = T)
###cheking
##sort(unique(values(TF2020.mask)))
###handle overlaps
#TF2020.mask[TF2020.mask>1] <- 1 #sf over upf
##plot(TF2020.mask)
#writeRaster(TF2020.mask, "rasters/PGM/all_forest_mask/PGM_2020_real.tif", format = "GTiff", overwrite = T)
#
#
##
#
#
## scenario avoid degradation
#UPF.avoiddegrad.mask <- UPF.avoiddegrad
#
#TF.avoiddegrad.mask <- sum(UPF.avoiddegrad.mask, DPF2010.mask, SF2020.mask, na.rm = T)
###cheking
##sort(unique(values(TF.avoiddegrad.mask)))
###handle overlaps
#TF.avoiddegrad.mask[TF.avoiddegrad.mask>1] <- 1
##plot(TF.avoiddegrad.mask)
#writeRaster(TF.avoiddegrad.mask, "rasters/PGM/all_forest_mask/PGM_2020_avoiddegrad.tif", format = "GTiff", overwrite = T)
#
#
##
#
#
## scenario avoid deforestation
#UPF.avoiddefor.mask <- UPF.avoiddefor
#
#TF.avoiddefor.mask <- sum(UPF.avoiddefor, DPF2020.mask, SF2010.mask, na.rm = T)
###cheking
##sort(unique(values(TF.avoiddefor.mask)))
###handle overlaps
#TF.avoiddefor.mask[TF.avoiddefor.mask>1] <- 1
##plot(TF.avoiddefor.mask)
#writeRaster(TF.avoiddefor.mask, "rasters/PGM/all_forest_mask/PGM_2020_avoiddeforest.tif", format = "GTiff", overwrite = T)
#
#
##
#
#
## scenario avoid deforestation upf only
#UPF.avoiddefor2.mask <- UPF.avoiddefor2
#
#TF.avoiddefor2.mask <- sum(UPF.avoiddefor2.mask, DPF2020.mask, SF2020.mask, na.rm = T)
###cheking
##sort(unique(values(TF.avoiddefor2.mask)))
###handle overlaps
#TF.avoiddefor2.mask[TF.avoiddefor2.mask>1] <- 1
##plot(TF.avoiddefor2.mask)
#writeRaster(TF.avoiddefor2.mask, "rasters/PGM/all_forest_mask/PGM_2020_avoiddeforest2.tif", format = "GTiff", overwrite = T)
#
#
##
#
#
## scenario avoid both
#UPF.avoidboth.mask <- UPF.avoidboth
#
#TF.avoidboth.mask <- sum(UPF.avoidboth.mask, DPF2010.mask, SF2010.mask, na.rm = T)
###cheking
##sort(unique(values(TF.avoidboth.mask)))
###handle overlaps
#TF.avoidboth.mask[TF.avoidboth.mask>1] <- 1
##plot(TF.avoidboth.mask)
#writeRaster(TF.avoidboth.mask, "rasters/PGM/all_forest_mask/PGM_2020_avoidboth.tif", format = "GTiff", overwrite = T)
#
#
##
#
#
## scenario restoration without avoid
#SF2020.restore10.mask <- SF2020.restore10
#
#TF.restore10.a <- sum(UPF2020.mask, DPF2020.mask, SF2020.restore10.mask, na.rm = T)
###cheking
##sort(unique(values(TF.restore10.a)))
#TF.restore10.a[TF.restore10.a>1] <- 1
##plot(TF.restore10.a)
#writeRaster(TF.restore10.a, "rasters/PGM/all_forest_mask/PGM_2020_restor_wo_avoid.tif", format = "GTiff", overwrite = T)
#
##
#
#
## scenario restoration and avoid deforestation
#SF2010.restore10.mask <- SF2010.restore10
#
#TF.restore10.b <- sum(UPF.avoiddefor.mask, DPF2020.mask, SF2010.restore10.mask, na.rm = T)
###cheking
##sort(unique(values(TF.restore10.b)))
#TF.restore10.b[TF.restore10.b>1] <- 1
##plot(TF.restore10.b)
#writeRaster(TF.restore10.b, "rasters/PGM/all_forest_mask/PGM_2020_restor_n_avoid_deforest.tif", format = "GTiff", overwrite = T)
#
##
#
#
## scenario restoration and avoid both
#TF.restore10.c <- sum(UPF.avoidboth.mask, DPF2010.mask, SF2010.restore10.mask, na.rm = T)
###cheking
##sort(unique(values(TF.restore10.c)))
#TF.restore10.c[TF.restore10.c>1] <- 1
##plot(TF.restore10.c)
#writeRaster(TF.restore10.c, "rasters/PGM/all_forest_mask/PGM_2020_restor_n_avoid_both.tif", format = "GTiff", overwrite = T)
#
##
##optional
##save.image("~/conserv_opportunities_jamesthomson/github_repo/pgm_environment.RData")
#rm(list=ls())
#
##
#
#
########################################################################################################################.
## creating land-use land-cover maps for each scenario
##dir.create("rasters/PGM/lulc", recursive = T)
#
#
##PGM 2010
##undegradded primary forest == 3
#UPF2010 <- raster("rasters/PGM/input/UPF2010_real.tif")
#UPF2010[UPF2010==1]<-3
#
##degradded primary forest == 30
#DPF2010 <- raster("rasters/PGM/input/DPF2010_real.tif")
#DPF2010[DPF2010==1]<-30
#
##secondary forest == 12
#SF2010 <- raster("rasters/PGM/input/SF2010_real.tif")
#SF2010[SF2010==1]<-12
#
#LULC2010 <- UPF2010 + DPF2010 + SF2010
##LULC2010[LULC2010==15]<-12
##plot(LULC2010)
##table(values(LULC2010))
#
#
##PGM 2020
##undegradded primary forest == 3
#UPF2020 <- raster("rasters/PGM/input/UPF2020_real.tif")
#UPF2020[UPF2020==1]<-3
#
##degradded primary forest == 30
#DPF2020 <- raster("rasters/PGM/input/DPF2020_real.tif")
#DPF2020[DPF2020==1]<-30
#
##secondary forest == 12
#SF2020 <- raster("rasters/PGM/input/SF2020_real.tif")
#SF2020[SF2020==1]<-12
#
#LULC2020 <- UPF2020 + DPF2020 + SF2020
##LULC2020[LULC2020==15]<-12
##plot(LULC2020)
##table(values(LULC2020))
#
#
##PGM 2020_avoiddegrad
##undegradded primary forest == 3
#UPF2020_avoiddegrad <- raster("rasters/PGM/input/UPF2020_avoiddegrad.tif")
#UPF2020_avoiddegrad[UPF2020_avoiddegrad==1]<-3
#
#LULC2020_avoiddegrad <- UPF2020_avoiddegrad + DPF2010 + SF2020
##LULC2020_avoiddegrad[LULC2020_avoiddegrad==33]<-30
##LULC2020_avoiddegrad[LULC2020_avoiddegrad==15]<-12
##LULC2020_avoiddegrad[LULC2020_avoiddegrad==42]<-12
##plot(LULC2020_avoiddegrad)
##table(values(LULC2020_avoiddegrad))
#
#
##PGM 2020_avoiddeforest
##undegradded primary forest == 3
#UPF2020_avoiddeforest <- raster("rasters/PGM/input/UPF2020_avoiddeforest.tif")
#UPF2020_avoiddeforest[UPF2020_avoiddeforest==1]<-3
#
#LULC2020_avoiddeforest <- UPF2020_avoiddeforest + DPF2020 + SF2010
##LULC2020_avoiddeforest[LULC2020_avoiddeforest==15]<-12
##LULC2020_avoiddeforest[LULC2020_avoiddeforest==42]<-12
##plot(LULC2020_avoiddeforest)
##table(values(LULC2020_avoiddeforest))
#
#
##PGM 2020_avoiddeforest upf only
##undegradded primary forest == 3
#UPF2020_avoiddeforest2 <- raster("rasters/PGM/input/UPF2020_avoiddeforest2.tif")
#UPF2020_avoiddeforest2[UPF2020_avoiddeforest2==1]<-3
#
#LULC2020_avoiddeforest2 <- UPF2020_avoiddeforest2 + DPF2020 + SF2010
##LULC2020_avoiddeforest2[LULC2020_avoiddeforest2==33]<-30
##LULC2020_avoiddeforest2[LULC2020_avoiddeforest2==45]<-30
##LULC2020_avoiddeforest2[LULC2020_avoiddeforest2==15]<-12
##LULC2020_avoiddeforest2[LULC2020_avoiddeforest2==42]<-12
##plot(LULC2020_avoiddeforest2)
##table(values(LULC2020_avoiddeforest2))
#
#
##PGM 2020_avoidboth
##undegradded primary forest == 3
#UPF2020_avoidboth <- raster("rasters/PGM/input/UPF2020_avoidboth.tif")
#UPF2020_avoidboth[UPF2020_avoidboth==1]<-3
#
#LULC2020_avoidboth <- UPF2020_avoidboth + DPF2010 + SF2010
##LULC2020_avoidboth[LULC2020_avoidboth==33]<-30
##LULC2020_avoidboth[LULC2020_avoidboth==45]<-30
##LULC2020_avoidboth[LULC2020_avoidboth==15]<-12
##LULC2020_avoidboth[LULC2020_avoidboth==42]<-12
##plot(LULC2020_avoidboth)
##table(values(LULC2020_avoidboth))
#
#
##PGM 2020_restor_wo_avoid
##secondary forest == 12
#SF2020.restore10 <- raster("rasters/PGM/input/SF2020_restor_wo_avoid.tif")
#SF2020.restore10[SF2020.restore10==1]<-12
#
#LULC2020_restor_wo_avoid <- UPF2020 + DPF2020 + SF2020.restore10
##LULC2020_restor_wo_avoid[LULC2020_restor_wo_avoid==15]<-12
##LULC2020_restor_wo_avoid[LULC2020_restor_wo_avoid==42]<-12
##plot(LULC2020_restor_wo_avoid)
##table(values(LULC2020_restor_wo_avoid))
#
#
##PGM 2020_restor_n_avoid_deforest
##secondary forest == 12
#SF2010.restore10 <- raster("rasters/PGM/input/SF2020_restor_n_avoid.tif")
#SF2010.restore10[SF2010.restore10==1]<-12
#
#LULC2020_restor_n_avoid_deforest <- UPF2020_avoiddeforest + DPF2020 + SF2010.restore10
##LULC2020_restor_n_avoid_deforest[LULC2020_restor_n_avoid_deforest==15]<-12
##LULC2020_restor_n_avoid_deforest[LULC2020_restor_n_avoid_deforest==42]<-12
##plot(LULC2020_restor_n_avoid_deforest)
##table(values(LULC2020_restor_n_avoid_deforest))
#
#
##PGM 2020_restor_n_avoid_both
#LULC2020_restor_n_avoid_both <- UPF2020_avoidboth + DPF2010 + SF2010.restore10
##LULC2020_restor_n_avoid_both[LULC2020_restor_n_avoid_both==33]<-30
##LULC2020_restor_n_avoid_both[LULC2020_restor_n_avoid_both==45]<-30
##LULC2020_restor_n_avoid_both[LULC2020_restor_n_avoid_both==15]<-12
##LULC2020_restor_n_avoid_both[LULC2020_restor_n_avoid_both==42]<-12
##plot(LULC2020_restor_n_avoid_both)
##table(values(LULC2020_restor_n_avoid_both))
#
#
## [transition-plot] sankey diagram
#LULC <- stack(LULC2010, LULC2020, LULC2020_avoiddegrad, LULC2020_avoiddeforest, LULC2020_avoiddeforest2, LULC2020_avoidboth, 
#              LULC2020_restor_wo_avoid, LULC2020_restor_n_avoid_deforest, LULC2020_restor_n_avoid_both)
#LULC <- mask(LULC, pgm.shp)
#names(LULC) <- c("2010_real", "2020_real", "2020_avoiddegrad", "2020_avoiddeforest", "2020_avoiddeforest2", "2020_avoidboth", 
#             "2020_restor_wo_avoid", "2020_restor_n_avoid_deforest", "2020_restor_n_avoid_both")
##plot(LULC, nr=2)
#LULC.df <- as.data.frame(LULC, xy = TRUE)
#
#
#pgm.lulc.df <- LULC.df %>% dplyr::select(x, y, X2010_real, X2020_real) %>% 
#  pivot_longer(
#    X2010_real:X2020_real,
#    names_to = "ID",
#    values_to = "Class"
#  ) %>% 
#  mutate(
#    Scenario = factor(case_when(str_detect(ID, "X2010_real")~ "2010 Real",
#                                str_detect(ID, "X2020_real")~ "2020 Real"),
#                      levels = c("2010 Real",
#                                 "2020 Real")),
#    Class = factor(Class, 
#                   levels = c(3, 30, 12, 0), 
#                   labels = c("Undegraded Primary Forest (UPF)", "Degraded Forest (DPF)", "Secondary Forest (SF)", "Non-Forest (Def)"))
#  )
#
#
#real.tp <- data.frame(
#  Period1 = factor(LULC[["X2010_real"]][], levels = c(3, 30, 12, 0), labels = c("UPF", "DPF", "SF", "Def")),
#  Period2 = factor(LULC[["X2020_real"]][], levels = c(3, 30, 12, 0), labels = c("UPF", "DPF", "SF", "Def"))
#)
#
#
#
#library(ggalluvial)
##> Loading required package: ggplot2
#library(ggfittext)
#
#
#addline_format <- function(x,...){
#  gsub('\\s','\n',x)
#}
#
#gridExtra::grid.arrange(
#  
#  pgm.lulc.df %>% filter(ID == "X2010_real") %>% 
#    ggplot() +
#    geom_raster(aes(x = x, y = y, fill = Class)) +
#    scale_fill_manual(values = c("#263A29", "#65451F", "#83764F", "#F97B22"), na.value = NA) +
#    labs(title = "2010 Real", x = "", y = "Latitude") +
#    theme_minimal() +
#    theme(legend.position = c(.2,.8)),
#  
#  real.tp %>% drop_na() %>% sample_n(size = 50000, replace = T) %>%
#    ggplot(aes(axis1 = Period1, axis2 = Period2)) +
#    geom_flow(aes(fill = Period1), width = .15, curve_type = "quintic") +
#    geom_stratum(width = .15) +
#    scale_x_discrete(limits = c("Period1", "Period2"), 
#                     breaks=c("Period1", "Period2"), 
#                     labels=addline_format(c("2010 Real", "2020 Real")),
#                     expand = c(.05, .05)) +
#    scale_fill_manual(values = c("#263A29", "#65451F", "#83764F", "#F97B22")) +
#    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
#    theme_minimal()+
#    theme(axis.text.y= element_blank(), legend.position = "none"),
#  
#  pgm.lulc.df %>% filter(ID == "X2020_real") %>% 
#    ggplot() +
#    geom_raster(aes(x = x, y = y, fill = Class)) +
#    scale_fill_manual(values = c("#263A29", "#65451F", "#83764F", "#F97B22"), na.value = NA) +
#    labs(title = "2020 Real", x = "Longitude", y = "Latitude") +
#    theme_minimal() +
#    theme(legend.position = "none"),
#  
#  # box plot and scatter plot
#  ncol = 2, nrow = 2, 
#  layout_matrix = rbind(c(1,2), c(3,2))
#  
#)
#
#
#
#avoiddegrad.tp <- data.frame(
#  Period1 = factor(LULC[["X2010_real"]][], levels = c(3, 30, 12, 0), labels = c("UPF", "DPF", "SF", "Def")),
#  Period2 = factor(LULC[["X2020_avoiddegrad"]][], levels = c(3, 30, 12, 0), labels = c("UPF", "DPF", "SF", "Def"))
#)
#
#avoiddeforest.tp <- data.frame(
#  Period1 = factor(LULC[["X2010_real"]][], levels = c(3, 30, 12, 0), labels = c("UPF", "DPF", "SF", "Def")),
#  Period2 = factor(LULC[["X2020_avoiddeforest"]][], levels = c(3, 30, 12, 0), labels = c("UPF", "DPF", "SF", "Def"))
#)
#
#avoiddeforest2.tp <- data.frame(
#  Period1 = factor(LULC[["X2010_real"]][], levels = c(3, 30, 12, 0), labels = c("UPF", "DPF", "SF", "Def")),
#  Period2 = factor(LULC[["X2020_avoiddeforest2"]][], levels = c(3, 30, 12, 0), labels = c("UPF", "DPF", "SF", "Def"))
#)
#
#avoidboth.tp <- data.frame(
#  Period1 = factor(LULC[["X2010_real"]][], levels = c(3, 30, 12, 0), labels = c("UPF", "DPF", "SF", "Def")),
#  Period2 = factor(LULC[["X2020_avoidboth"]][], levels = c(3, 30, 12, 0), labels = c("UPF", "DPF", "SF", "Def"))
#)
#
#restor_wo_avoid.tp <- data.frame(
#  Period1 = factor(LULC[["X2010_real"]][], levels = c(3, 30, 12, 0), labels = c("UPF", "DPF", "SF", "Def")),
#  Period2 = factor(LULC[["X2020_restor_wo_avoid"]][], levels = c(3, 30, 12, 0), labels = c("UPF", "DPF", "SF", "Def"))
#)
#
#restor_n_avoid_deforest.tp <- data.frame(
#  Period1 = factor(LULC[["X2010_real"]][], levels = c(3, 30, 12, 0), labels = c("UPF", "DPF", "SF", "Def")),
#  Period2 = factor(LULC[["X2020_restor_n_avoid_deforest"]][], levels = c(3, 30, 12, 0), labels = c("UPF", "DPF", "SF", "Def"))
#)
#
#restor_n_avoid_both.tp <- data.frame(
#  Period1 = factor(LULC[["X2010_real"]][], levels = c(3, 30, 12, 0), labels = c("UPF", "DPF", "SF", "Def")),
#  Period2 = factor(LULC[["X2020_restor_n_avoid_both"]][], levels = c(3, 30, 12, 0), labels = c("UPF", "DPF", "SF", "Def"))
#)
#
#
#
#
#
## Create the alluvial plot
#
#addline_format <- function(x,...){
#  gsub('\\s','\n',x)
#}
#
#cowplot::plot_grid(
#  
#  real.tp %>% drop_na() %>% sample_n(size = 50000, replace = T) %>%
#  ggplot(aes(axis1 = Period1, axis2 = Period2)) +
#  geom_flow(aes(fill = Period1), width = .15, curve_type = "quintic") +
#  geom_stratum(width = .15) +
#  scale_x_discrete(limits = c("Period1", "Period2"), 
#                   breaks=c("Period1", "Period2"), 
#                   labels=addline_format(c("2010 Real", "2020 Real")),
#                   expand = c(.05, .05)) +
#  scale_fill_manual(values = c("#263A29", "#65451F", "#83764F", "#F97B22")) +
#  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
#  theme_minimal()+
#  theme(axis.text.y= element_blank(), legend.position = "none"),
#  
#  avoiddegrad.tp %>% drop_na() %>% sample_n(size = 50000, replace = T) %>%
#    ggplot(aes(axis1 = Period1, axis2 = Period2)) +
#    geom_flow(aes(fill = Period1), width = .15, curve_type = "quintic") +
#    geom_stratum(width = .15) +
#    scale_x_discrete(limits = c("Period1", "Period2"), 
#                     breaks=c("Period1", "Period2"), 
#                     labels=addline_format(c("2010 Real", "Avoid degradation")),
#                     expand = c(.05, .05)) +
#    scale_fill_manual(values = c("#263A29", "#65451F", "#83764F", "#F97B22")) +
#    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
#    theme_minimal()+
#    theme(axis.text.y= element_blank(), legend.position = "none"),
#  
#  avoiddeforest.tp %>% drop_na() %>% sample_n(size = 50000, replace = T) %>%
#    ggplot(aes(axis1 = Period1, axis2 = Period2)) +
#    geom_flow(aes(fill = Period1), width = .15, curve_type = "quintic") +
#    geom_stratum(width = .15) +
#    scale_x_discrete(limits = c("Period1", "Period2"), 
#                     breaks=c("Period1", "Period2"), 
#                     labels=addline_format(c("2010 Real", "Avoid deforestation")),
#                     expand = c(.05, .05)) +
#    scale_fill_manual(values = c("#263A29", "#65451F", "#83764F", "#F97B22")) +
#    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
#    theme_minimal()+
#    theme(axis.text.y= element_blank(), legend.position = "none"),
#  
#  avoiddeforest2.tp %>% drop_na() %>% sample_n(size = 50000, replace = T) %>%
#    ggplot(aes(axis1 = Period1, axis2 = Period2)) +
#    geom_flow(aes(fill = Period1), width = .15, curve_type = "quintic") +
#    geom_stratum(width = .15) +
#    scale_x_discrete(limits = c("Period1", "Period2"), 
#                     breaks=c("Period1", "Period2"), 
#                     labels=addline_format(c("2010 Real", "Avoid deforestation primary forest")),
#                     expand = c(.05, .05)) +
#    scale_fill_manual(values = c("#263A29", "#65451F", "#83764F", "#F97B22")) +
#    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
#    theme_minimal()+
#    theme(axis.text.y= element_blank(), legend.position = "none"),
#  
#  avoidboth.tp %>% drop_na() %>% sample_n(size = 50000, replace = T) %>%
#    ggplot(aes(axis1 = Period1, axis2 = Period2)) +
#    geom_flow(aes(fill = Period1), width = .15, curve_type = "quintic") +
#    geom_stratum(width = .15) +
#    scale_x_discrete(limits = c("Period1", "Period2"), 
#                     breaks=c("Period1", "Period2"), 
#                     labels=addline_format(c("2010 Real", "Avoid both")),
#                     expand = c(.05, .05)) +
#    scale_fill_manual(values = c("#263A29", "#65451F", "#83764F", "#F97B22")) +
#    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
#    theme_minimal()+
#    theme(axis.text.y= element_blank(), legend.position = "none"),
#  
#  restor_wo_avoid.tp %>% drop_na() %>% sample_n(size = 50000, replace = T) %>%
#    ggplot(aes(axis1 = Period1, axis2 = Period2)) +
#    geom_flow(aes(fill = Period1), width = .15, curve_type = "quintic") +
#    geom_stratum(width = .15) +
#    scale_x_discrete(limits = c("Period1", "Period2"), 
#                     breaks=c("Period1", "Period2"), 
#                     labels=addline_format(c("2010 Real", "Restoration without avoid")),
#                     expand = c(.05, .05)) +
#    scale_fill_manual(values = c("#263A29", "#65451F", "#83764F", "#F97B22")) +
#    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
#    theme_minimal()+
#    theme(axis.text.y= element_blank(), legend.position = "none"),
#  
#  restor_n_avoid_deforest.tp %>% drop_na() %>% sample_n(size = 50000, replace = T) %>%
#    ggplot(aes(axis1 = Period1, axis2 = Period2)) +
#    geom_flow(aes(fill = Period1), width = .15, curve_type = "quintic") +
#    geom_stratum(width = .15) +
#    scale_x_discrete(limits = c("Period1", "Period2"), 
#                     breaks=c("Period1", "Period2"), 
#                     labels=addline_format(c("2010 Real", "Restoration and avoid deforestation")),
#                     expand = c(.05, .05)) +
#    scale_fill_manual(values = c("#263A29", "#65451F", "#83764F", "#F97B22")) +
#    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
#    theme_minimal()+
#    theme(axis.text.y= element_blank(), legend.position = "none"),
#  
#  restor_n_avoid_both.tp %>% drop_na() %>% sample_n(size = 50000, replace = T) %>%
#    ggplot(aes(axis1 = Period1, axis2 = Period2)) +
#    geom_flow(aes(fill = Period1), width = .15, curve_type = "quintic") +
#    geom_stratum(width = .15) +
#    scale_x_discrete(limits = c("Period1", "Period2"), 
#                     breaks=c("Period1", "Period2"), 
#                     labels=addline_format(c("2010 Real", "Restoration and avoid both")),
#                     expand = c(.05, .05)) +
#    scale_fill_manual(values = c("#263A29", "#65451F", "#83764F", "#F97B22")) +
#    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
#    theme_minimal()+
#    theme(axis.text.y= element_blank(), legend.position = "none"),
#  
#  ncol = 4, align = "hv")
#  
#
#
#
#
#
#
#
#
#
#
#
########################################################################################################################.#######################################################################################################################.
#
###### detecting multicollinearity between exploratory variables ####.
##env.explanatory.var.list <- list.files("rasters/PGM/2010_real", pattern = ".tif", full.names = T, recursive = T)
##
##env.explanatory.var <- stack(env.explanatory.var.list)
##names(env.explanatory.var) <- unlist(strsplit(env.explanatory.var.list, "/|.tif"))[seq(4,80,4)]
####cheking
###env.explanatory.var
###plot(env.explanatory.var[[1:10]], nc=2)
###plot(env.explanatory.var[[11:20]], nc=2)
##
### visual inspection of aggregation using removeCollinearity() function from virtualspecies package
##correlated.var <- removeCollinearity(env.explanatory.var, multicollinearity.cutoff = 0.7, sample.points = T, nb.points = 999999, method = "pearson", plot = T)
####cheking
###correlated.var
##
### variation inflation factor
##inflated.var <- vifcor(env.explanatory.var, th = 0.7, maxobservations = 999999)
####cheking
###inflated.var@results
###inflated.var@excluded
##
##sel.var.df <- data.frame(rbind(cbind(VAR=inflated.var@results$Variables, VIF=inflated.var@results$VIF),
##                               cbind(VAR=inflated.var@excluded, VIF=NA)))
##
##write.csv(sel.var.df, paste0("rasters/PGM/selected_environmental_explanatory_variables_byVIF.csv", sep=""), row.names = F)
##
####cheking
###env.explanatory.var <- env.explanatory.var[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
###plot(env.explanatory.var)
##
##
###
##
##rm(list=ls()) #keeping only raster stack
##gc()
##
##
##
###


