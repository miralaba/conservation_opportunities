
#' @title Cost-effectiveness of conservation actions in Amazon
#' @description 

#### loading required packages ####
library(tidyverse)
library(raster)
library(rgdal)
library(rgeos)
library(sf)
library(scales)
library(reshape)
library(spThin)
library(usdm)
library(virtualspecies)

#### detecting multicollinearity between exploratory variables ####
pgm.env.explanatory.var.list <- list.files("rasters/PGM/2010_real", pattern = ".tif", full.names = T, recursive = T)
pgm.env.explanatory.var <- stack(pgm.env.explanatory.var.list)
#names(pgm.env.explanatory.var) <- unlist(strsplit(pgm.env.explanatory.var.list, "/|.tif"))[seq(4,88,4)]
##cheking
#pgm.env.explanatory.var
#plot(pgm.env.explanatory.var[[1:10]], nc=2)
#plot(pgm.env.explanatory.var[[11:20]], nc=2)



stm.env.explanatory.var.list <- list.files("rasters/STM/2010_real", pattern = ".tif", full.names = T, recursive = T)
stm.env.explanatory.var <- stack(stm.env.explanatory.var.list)
#names(stm.env.explanatory.var) <- unlist(strsplit(stm.env.explanatory.var.list, "/|.tif"))[seq(4,88,4)]
##cheking
#stm.env.explanatory.var
#plot(stm.env.explanatory.var[[1:8]], nc=2)
#plot(stm.env.explanatory.var[[9:16]], nc=2)
#plot(stm.env.explanatory.var[[17:22]], nc=2)


# merging pgm and stm 2010
env.explanatory.var <- stack()

for (var in names(pgm.env.explanatory.var)) {
  
  cat("\n>working on", var, "now<\n")
  
  var.x <- merge(pgm.env.explanatory.var[[var]], stm.env.explanatory.var[[var]])
  names(var.x)<-var
  env.explanatory.var <- addLayer(env.explanatory.var, var.x)
  
  
}

##cheking
#env.explanatory.var
#plot(env.explanatory.var[[1:10]], nc=2)
#plot(env.explanatory.var[[11:20]], nc=2)



#


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

sel.var.df[sel.var.df$VAR=="TSDls", "VIF"] <- 999
sel.var.df[sel.var.df$VAR=="UPFpx", "VIF"] <- 999
sel.var.df$Type <- "environment"
sel.var.df[sel.var.df$VAR %in% c("distmarket", "propertysize"), "Type"] <- "economic"

write.csv(sel.var.df, "rasters/selected_environmental_explanatory_variables_byVIF.csv", row.names = F)

##cheking
#env.explanatory.var <- env.explanatory.var[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
#plot(env.explanatory.var)


#

rm(list=ls()[ls() %in% c("env.explanatory.var.low", "correlated.var", "inflated.var")]) #keeping only raster stack
gc()



#


#######################################################################################################################

#### import data ####
## transect data
transectdata <- read.csv("data/raw/RAS_transects_environment_all.csv")
#head(transectdata)
#str(transectdata)
#summary(transectdata)

# excluding Varzea transects PGM
exclude <- c("100_1", "100_4","100_7","81_12","423_2") #transect code
transectdata <- transectdata[!transectdata$Transectcode %in% exclude,]


## birds
birddata <- read.csv("data/raw/Bird_data_Standard_NM_26022013_Final.csv")
#head(birddata)
#str(birddata)
#summary(birddata)

# rename columns to match with transect data
names(birddata)[1] <- "Region"
names(birddata)[3] <- "Catchment"
names(birddata)[5] <- "Transectcode"

# new var to mark data frame as birds
birddata$Group <- "Birds"

## trees
pgm.treedata <- read.csv("data/raw/Flora.composition.and.biomass_PGM_Erika_23.01.2013.csv")
stm.treedata <- read.csv("data/raw/Flora.composition.and.biomass_STM_Erika_23.01.2013.csv")
# include region code
pgm.treedata$Region <- "PGM"
stm.treedata$Region <- "STM"
# merging tree data
treedata <- rbind(pgm.treedata, stm.treedata)
#head(treedata)
#str(treedata)
#summary(treedata)

# new var to mark data frame as birds
treedata$Group <- "Trees"
# new var with code for transect by catchment
treedata$Transectcode <- paste(treedata$Catchment, treedata$Transect, sep="_")
# new var with spp binomial
treedata$Binomial <- paste(treedata$Genera, treedata$Species, sep="")


## merging spp data with binomial, location and lulc code
transectdata <- transectdata[,c(1:6,10)]
birddata <- birddata[,c(23,1,3,4,5,15)]
treedata <- treedata[,c(36,35,1,2,37,38)]
sppdata <- rbind(birddata, treedata)
sppdata <- merge(sppdata, transectdata)
#head(treedata)
#str(treedata)
#summary(treedata)


# new var with longlat
pgm.utm <- data.frame(x=sppdata[sppdata$Region=="PGM", "UTM_X"], y=sppdata[sppdata$Region=="PGM", "UTM_Y"]) 
coordinates(pgm.utm) <- ~x+y 
#class(pgm.utm)
proj4string(pgm.utm) <- CRS("+proj=utm +zone=23 +south +datum=WGS84 +units=m +ellps=WGS84") 
pgm.longlat <- spTransform(pgm.utm, CRS("+proj=longlat +datum=WGS84"))

stm.utm <- data.frame(x=sppdata[sppdata$Region=="STM", "UTM_X"], y=sppdata[sppdata$Region=="STM", "UTM_Y"]) 
coordinates(stm.utm) <- ~x+y 
#class(stm.utm)
proj4string(stm.utm) <- CRS("+proj=utm +zone=21 +south +datum=WGS84 +units=m +ellps=WGS84") 
stm.longlat <- spTransform(stm.utm, CRS("+proj=longlat +datum=WGS84"))

longlat <- rbind(pgm.longlat, stm.longlat)

sppdata$Longitude <- longlat@coords[,1]
sppdata$Latitude <- longlat@coords[,2]
#plot(env.explanatory.var[["UPFls"]])
#points(longlat)


# removing withespace and blank cells
sppdata$Binomial <- iconv(sppdata$Binomial, from = "ISO-8859-1", to = "UTF-8")
sppdata$Binomial <- gsub("\\s+","",sppdata$Binomial)
sppdata <- sppdata[!(is.na(sppdata$Binomial) | sppdata$Binomial==""), ]

# removing unidentifed species
#sort(unique(sppdata[grep(".sp$", sppdata$Binomial),"Binomial"]))
sppdata <- sppdata[-grep(".sp$", sppdata$Binomial),]
sppdata <- sppdata[-grep(".sp1$", sppdata$Binomial),]

#
## filter 1: selecting forest dependent species using undisturbed primary forest layer
UFPls.core <- env.explanatory.var[["UPFls"]]

UFPls.core[UFPls.core>=.8]<-1
UFPls.core[UFPls.core<.8]<-0
## checking
#plot(UFPls.core)
#points(SpatialPoints(sppdata[,c("Longitude", "Latitude")]))


sppdata$forestdep <- extract(UFPls.core, SpatialPoints(sppdata[, c("Longitude", "Latitude")]))

forestdep.spplist <- as.data.frame(sppdata %>% group_by(Binomial, forestdep) %>% 
                                       summarise(Group=first(Group), n=n()) %>% 
                                       mutate(freq=n/sum(n)) %>% 
                                       filter(forestdep==1 & freq >= .25) %>% 
                                       ungroup())

forestdep.sppdata <- sppdata[sppdata$Binomial %in% forestdep.spplist$Binomial,]

#
# filter 2: excluding records to avoid spatial autocorrelation
sppdata.final <- forestdep.sppdata
sppdata.final <- sppdata.final[-(1:nrow(sppdata.final)),]
for (i in forestdep.spplist$Binomial) {
  
  # filter one species at a time
  spp_x <- forestdep.sppdata[forestdep.sppdata$Binomial==i,]
  
  # record spatial thinning
  thinned_occur <- thin(loc.data = spp_x, 
                        lat.col = "Latitude", long.col = "Longitude", 
                        spec.col = "Binomial", 
                                        # the distance (km) that records are separated by -- 
                        thin.par = 1.5, # equals the distance between transects 
                                        # (see study transect distribution session in material e methods)
                        reps = 30, 
                        locs.thinned.list.return = TRUE, 
                        write.files = FALSE, 
                        write.log.file = FALSE)
  
  occur <- thinned_occur[[sample.int(30,1)]]
  occur$Binomial <- i
  
  forestdep.spplist[forestdep.spplist$Binomial==i,"Nrec"] <- nrow(occur)
  sppdata.final <- rbind(sppdata.final, occur)
  
  
  cat("\n>occurences for species", i, "thined<\n")
  
}


sppdata.final <- left_join(sppdata.final, forestdep.sppdata)
sppdata.final <- sppdata.final[!duplicated(sppdata.final),]
sppdata.final <- sppdata.final[,c(4:8,3,9,10,1,2,11)]
## checking
#nrow(as.data.frame(table(sppdata.final[sppdata.final$Group=="Trees","Binomial"])))
#nrow(as.data.frame(table(sppdata.final[sppdata.final$Group=="Birds","Binomial"])))

# add new var for ensemble model scores
forestdep.spplist<-forestdep.spplist[,-c(2,4,5)]
forestdep.spplist$job <- c(rep(1:12, each=floor(length(unique(sppdata.final$Binomial))/12)),12,12,12,12,12)
forestdep.spplist$job <- ifelse(forestdep.spplist$Nrec<10, 13, forestdep.spplist$job)
forestdep.spplist$Done <- "FALSE"

rm(list= ls()[!(ls() %in% c("env.explanatory.var", "sel.var.df", "forestdep.spplist", "sppdata.final"))])
gc()
#



write.csv(forestdep.spplist, "data/species_summary.csv", row.names = F)
write.csv(sppdata.final, "data/presence_records.csv", row.names = F)


