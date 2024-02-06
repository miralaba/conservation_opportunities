

# importing raw rasters ===================================|
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


# creating time since degradation based on mapbiomas fire, degrad and deter ====================|
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





# import fire time series from mapbiomas =============================|
pgm.raw.fire.list <- list.files("C:/Users/miral/Dropbox/MAPBIOMAS-EXPORT/fogo-time-series/PGM",
                                pattern = "mapbiomas", full.names = T, recursive = T)


pgm.raw.fire <- raster::stack(pgm.raw.fire.list)


pgm.degrad.list <- list.files("C:/Users/miral/Dropbox/MAPBIOMAS-EXPORT/fogo-time-series/PGM",
                                pattern = "degrad", full.names = T, recursive = T)


pgm.degrad <- raster::stack(pgm.degrad.list)


pgm.deter.list <- list.files("C:/Users/miral/Dropbox/MAPBIOMAS-EXPORT/fogo-time-series/PGM",
                              pattern = "deter", full.names = T, recursive = T)


pgm.deter <- raster::stack(pgm.deter.list)







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


##2007-2010 (combining fire and degradation)
for (i in 23:26) {
  
  #selecting year-by-year
    #fire
    fire.yearx <- pgm.raw.fire[[i]]
    
    #degradation
    inpe.yearx1 <- pgm.degrad[[i-22]]
    inpe.yearx2 <- pgm.degrad[[i-21]]
    inpe.yearx3 <- pgm.degrad[[i-20]]
    
    inpe.yearx <- sum(inpe.yearx1, inpe.yearx2)
    inpe.yearx <- sum(inpe.yearx, inpe.yearx3)
    inpe.yearx[inpe.yearx!=3] <- 0
    inpe.yearx[inpe.yearx==3] <- 1
    
    degrad.yearx <- sum(inpe.yearx, fire.yearx)
    degrad.yearx[degrad.yearx>1] <- 1
    
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

pgm.2010.time.since.degrad <- pgm.time.since.degrad
pgm.2010.time.since.degrad <-  mask(pgm.2010.time.since.degrad, pgm.lulc.2010.forest.mask)
writeRaster(pgm.2010.time.since.degrad, "rasters/PGM/raw/pgm-2010-deg_tsince0.tif", format = "GTiff", overwrite = T)


pgm.1985.2010.freq.degrad <- pgm.freq.degrad
pgm.1985.2010.freq.degrad <-  mask(pgm.1985.2010.freq.degrad, pgm.lulc.2010.forest.mask)
writeRaster(pgm.1985.2010.freq.degrad, "rasters/PGM/raw/pgm-firefreq-1985_2010.tif", format = "GTiff", overwrite = T)

##2011-2015
for (i in 27:31) {
  
  #selecting year-by-year
  #fire
  fire.yearx <- pgm.raw.fire[[i]]
  
  #degradation
  inpe.yearx1 <- pgm.degrad[[i-22]]
  
  if(i!=31){
    inpe.yearx2 <- pgm.degrad[[i-21]]
    inpe.yearx3 <- pgm.degrad[[i-20]]
    
    inpe.yearx <- sum(inpe.yearx1, inpe.yearx2)
    inpe.yearx <- sum(inpe.yearx, inpe.yearx3)
    inpe.yearx[inpe.yearx!=3] <- 0
    inpe.yearx[inpe.yearx==3] <- 1
  }
  else
    {
      inpe.yearx2a <- pgm.degrad[[i-21]]
      inpe.yearx2b <- pgm.deter[[i-30]]
      inpe.yearx3 <- pgm.deter[[i-29]]
    
      inpe.yearx <- sum(inpe.yearx1, inpe.yearx2a)
      inpe.yearx <- sum(inpe.yearx, inpe.yearx2b)
      inpe.yearx <- sum(inpe.yearx, inpe.yearx3)
      inpe.yearx[inpe.yearx!=4] <- 0
      inpe.yearx[inpe.yearx==4] <- 1
  }
  
  
  degrad.yearx <- sum(inpe.yearx, fire.yearx)
  degrad.yearx[degrad.yearx>1] <- 1
  
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

##2016-2020
for (i in 32:36) {
  
  #selecting year-by-year
  #fire
  fire.yearx <- pgm.raw.fire[[i]]
  
  #degradation
  inpe.yearx1 <- pgm.deter[[i-31]]
  
  if(i<35){
    inpe.yearx2 <- pgm.deter[[i-30]]
    inpe.yearx3 <- pgm.deter[[i-29]]
    
    inpe.yearx <- sum(inpe.yearx1, inpe.yearx2)
    inpe.yearx <- sum(inpe.yearx, inpe.yearx3)
    inpe.yearx[inpe.yearx!=3] <- 0
    inpe.yearx[inpe.yearx==3] <- 1
  }
  if(i==35){
    inpe.yearx2 <- pgm.deter[[i-32]]
    inpe.yearx3 <- pgm.deter[[i-30]]
    
    inpe.yearx <- sum(inpe.yearx1, inpe.yearx2)
    inpe.yearx <- sum(inpe.yearx, inpe.yearx3)
    inpe.yearx[inpe.yearx!=3] <- 0
    inpe.yearx[inpe.yearx==3] <- 1
  }
  if(i==36){
    inpe.yearx2 <- pgm.deter[[i-32]]
    inpe.yearx3 <- pgm.deter[[i-33]]
    
    inpe.yearx <- sum(inpe.yearx1, inpe.yearx2)
    inpe.yearx <- sum(inpe.yearx, inpe.yearx3)
    inpe.yearx[inpe.yearx!=3] <- 0
    inpe.yearx[inpe.yearx==3] <- 1
  }
  
  
  degrad.yearx <- sum(inpe.yearx, fire.yearx)
  degrad.yearx[degrad.yearx>1] <- 1
  
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

pgm.2020.time.since.degrad <- pgm.time.since.degrad
pgm.2020.time.since.degrad <-  mask(pgm.2020.time.since.degrad, pgm.lulc.2020.forest.mask)
writeRaster(pgm.2020.time.since.degrad, "rasters/PGM/raw/pgm-2020-deg_tsince0.tif", format = "GTiff", overwrite = T)


pgm.1985.2020.freq.degrad <- pgm.freq.degrad
pgm.1985.2020.freq.degrad <-  mask(pgm.1985.2020.freq.degrad, pgm.lulc.2020.forest.mask)
writeRaster(pgm.1985.2020.freq.degrad, "rasters/PGM/raw/pgm-firefreq-1985_2020.tif", format = "GTiff", overwrite = T)

rm(list=ls()[ls() %in% c("pgm.raw.fire.list", "pgm.raw.fire", "pgm.degrad.list", "pgm.degrad",
                         "pgm.deter.list", "pgm.deter", "degrad.yearx", "fire.yearx", "inpe.yearx", "inpe.yearx1", 
                         "inpe.yearx2", "inpe.yearx2a", "inpe.yearx2b", "inpe.yearx3", "i")])
gc()
#
#
















#use this part of the script only after load objects in "layer_buide_[...].R"

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







































