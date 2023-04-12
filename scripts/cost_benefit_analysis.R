
#' @title Cost-effectiveness of conservation actions in Amazon
#' @description 

#### pre-setting ####
memory.limit(1000000)

#### loading required packages ####
library(tidyverse)
library(raster)
library(rgdal)
library(rgeos)
library(sf)
library(scales)
library(reshape)
library(stringr)


# importing species list and methods and evaluation metrics data frame
forestdep.spplist <- read.csv("data/species_summary_final.csv")

#### adding variable on conservation value ####

# bird conservation value is inverse occurrence area size
# scaled from 0 [the biggest] to 1 [the smallest]
# based on birdlife_v2017b shapefiles
bird.distribution.shp <- readOGR(dsn="D:/GIS/shapes/bird_distr_BirdLifeInternational2016", layer="Birds_of_Amazonia")
head(bird.distribution.shp@data)
# adding variable to match with species list dataframe
bird.distribution.shp@data$Binomial <- gsub(" ", "", bird.distribution.shp@data$SCINAME)

# adding shape area from birdlife_v2017b shapefiles
forestdep.spplist <- forestdep.spplist %>% left_join(bird.distribution.shp@data[,c("Binomial", "Shape_Area")])

# scaling
forestdep.spplist$Shape_Area_scaled <- rescale(forestdep.spplist$Shape_Area, to = c(.999,.001))

#

rm(bird.distribution.shp)


# tree conservation value is wood density
# scaled from 0 [the biggest] to 1 [the smallest]
# based on the World Checklist of Vascular Plants
## trees
pgm.treedata <- read.csv("data/input/Flora.composition.and.biomass_PGM_Erika_23.01.2013.csv")
stm.treedata <- read.csv("data/input/Flora.composition.and.biomass_STM_Erika_23.01.2013.csv")
# merging tree data
treedata <- rbind(pgm.treedata, stm.treedata)
rm(pgm.treedata); rm(stm.treedata)
# new var with spp binomial
treedata$Binomial <- paste(treedata$Genera, treedata$Species, sep="")
treedata <- treedata %>% dplyr::select(Genera.SP.SPP, Binomial) %>% 
  distinct(.keep_all = T) %>% 
  filter(Binomial %in% forestdep.spplist$Binomial) %>% 
  mutate(Genera.SP.SPP = str_trim(str_replace_all(Genera.SP.SPP, "_", " ")))

treedata$Shape_Area <- NA


library(rWCVP)

j=nrow(treedata)
for (sp in treedata$Genera.SP.SPP) {
  
  test <- try(wcvp_distribution(treedata[treedata$Genera.SP.SPP==sp,"Genera.SP.SPP"], taxon_rank="species"), silent=TRUE)
  distribution <- if(class(test) %in% 'try-error') { next } 
  else {wcvp_distribution(treedata[treedata$Genera.SP.SPP==sp,"Genera.SP.SPP"], taxon_rank="species")}
  
  treedata[treedata$Genera.SP.SPP==sp,"Shape_Area"] <- sum(st_area(distribution))
  
  j=j-1
  cat("\n> done:", sp, "now", j, "species left <\n")
}

# scaling
treedata <- treedata %>% 
  mutate(Shape_Area_scaled = rescale(Shape_Area, to = c(.999,.001))) %>% 
  dplyr::select(Binomial:Shape_Area_scaled) %>% 
  distinct(Binomial, .keep_all = T)

# merging data
forestdep.spplist <- forestdep.spplist %>% rows_update(treedata, by = c("Binomial"))


# saving
write.csv(forestdep.spplist, "data/species_summary_final.csv")

#### biodiversity benefit ####









