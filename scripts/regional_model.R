
#' @title Cost-effectiveness of conservation actions in Amazon
#' @description 

# loading required packages ====================================================
library(tidyverse)
library(raster)
library(rgdal)
library(rgeos)
library(sf)
library(scales)
library(reshape)
library(spThin)
library(usdm)
library(biomod2)
library(fitdistrplus)
library(dismo)
library(caret)
library(randomForest)
library(gridExtra)
library(grDevices)

# loading input data ===========================================================
## standard projection
std.proj <- "+proj=longlat +datum=WGS84 +units=m +no_defs"

## shapefile paragominas
pgm.shp <- readOGR(dsn = "shapes", layer = "Paragominas_Mask_R3")
proj4string(pgm.shp) <- CRS(std.proj)
pgm.shp <- spTransform(pgm.shp, crs(std.proj))

## shapefile santarem
stm.shp <- readOGR(dsn = "shapes", layer = "Santarem")
stm.shp <- spTransform(stm.shp, crs(std.proj))


#explanatory variables
sel.var.df <- read.csv("rasters/selected_environmental_explanatory_variables_byVIF.csv")

pgm.env.explanatory.var.list <- list.files("rasters2/PGM/2010_real", pattern = ".tif", full.names = T, recursive = T)

pgm.env.explanatory.var <- stack(pgm.env.explanatory.var.list)
#names(pgm.env.explanatory.var) <- unlist(strsplit(pgm.env.explanatory.var.list, "/|.tif"))[seq(4,88,4)]

pgm.env.explanatory.var <- pgm.env.explanatory.var[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
##cheking
#pgm.env.explanatory.var
#plot(pgm.env.explanatory.var[[1:10]], nc=2)
#plot(pgm.env.explanatory.var[[11:20]], nc=2)



stm.env.explanatory.var.list <- list.files("rasters2/STM/2010_real", pattern = ".tif", full.names = T, recursive = T)

stm.env.explanatory.var <- stack(stm.env.explanatory.var.list)
#names(stm.env.explanatory.var) <- unlist(strsplit(stm.env.explanatory.var.list, "/|.tif"))[seq(4,88,4)]

stm.env.explanatory.var <- stm.env.explanatory.var[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
##cheking
#stm.env.explanatory.var
#plot(stm.env.explanatory.var[[1:10]], nc=2)
#plot(stm.env.explanatory.var[[11:20]], nc=2)



# merging pgm and stm 2010
env.explanatory.var <- stack()

for (var in names(pgm.env.explanatory.var)) {
  
  cat("\n>working on", var, "now<\n")
  
  var.x <- merge(pgm.env.explanatory.var[[var]], stm.env.explanatory.var[[var]])
  names(var.x)<-var
  env.explanatory.var <- addLayer(env.explanatory.var, var.x)
  
  
}


# creating directories =========================================================
dir.create("models.output/biodiversity.maps/PGM", recursive = T)
dir.create("models.output/biodiversity.maps/PGM/2010_real", recursive = T)
dir.create("models.output/biodiversity.maps/PGM/2020_real", recursive = T)
dir.create("models.output/biodiversity.maps/PGM/2020_avoiddeforest", recursive = T)
dir.create("models.output/biodiversity.maps/PGM/2020_avoiddeforest2", recursive = T)
dir.create("models.output/biodiversity.maps/PGM/2020_avoiddegrad", recursive = T)
dir.create("models.output/biodiversity.maps/PGM/2020_avoiddegrad2", recursive = T)
dir.create("models.output/biodiversity.maps/PGM/2020_avoidboth", recursive = T)
dir.create("models.output/biodiversity.maps/PGM/2020_avoidboth2", recursive = T)
dir.create("models.output/biodiversity.maps/PGM/2020_restor_wo_avoid", recursive = T)
dir.create("models.output/biodiversity.maps/PGM/2020_restor_n_avoiddeforest", recursive = T)
dir.create("models.output/biodiversity.maps/PGM/2020_restor_n_avoiddeforest2", recursive = T)
dir.create("models.output/biodiversity.maps/PGM/2020_restor_n_avoidboth", recursive = T)
dir.create("models.output/biodiversity.maps/PGM/2020_restor_n_avoidboth2", recursive = T)

dir.create("models.output/biodiversity.maps/STM", recursive = T)
dir.create("models.output/biodiversity.maps/STM/2010_real", recursive = T)
dir.create("models.output/biodiversity.maps/STM/2020_real", recursive = T)
dir.create("models.output/biodiversity.maps/STM/2020_avoiddeforest", recursive = T)
dir.create("models.output/biodiversity.maps/STM/2020_avoiddeforest2", recursive = T)
dir.create("models.output/biodiversity.maps/STM/2020_avoiddegrad", recursive = T)
dir.create("models.output/biodiversity.maps/STM/2020_avoiddegrad2", recursive = T)
dir.create("models.output/biodiversity.maps/STM/2020_avoidboth", recursive = T)
dir.create("models.output/biodiversity.maps/STM/2020_avoidboth2", recursive = T)
dir.create("models.output/biodiversity.maps/STM/2020_restor_wo_avoid", recursive = T)
dir.create("models.output/biodiversity.maps/STM/2020_restor_n_avoiddeforest", recursive = T)
dir.create("models.output/biodiversity.maps/STM/2020_restor_n_avoiddeforest2", recursive = T)
dir.create("models.output/biodiversity.maps/STM/2020_restor_n_avoidboth", recursive = T)
dir.create("models.output/biodiversity.maps/STM/2020_restor_n_avoidboth2", recursive = T)

#
#




# biodiversity benefit: for species with more than 5 records ===================
#response variables
forestdep.spplist <- read.csv("data/species_summary.csv")
#forestdep.spplist <- read.csv("models.output/evaluation/updated_species_summary_v1.csv")

sppdata.final <- read.csv("data/presence_records.csv")

# environment for pseudo-absence data 
# (i.e., areas with less than 20% of undegraded primary forest in landscape)
UPFls.invcore <- merge(mask(raster(pgm.env.explanatory.var.list[23]), pgm.shp), 
                       mask(raster(stm.env.explanatory.var.list[23]), stm.shp))
UPFls.invcore[UPFls.invcore>.2]<-NA
UPFls.invcore[UPFls.invcore<=.2]<-1


#dataframe to store model parameters, quality assessment, and variable importance
namevector<-c("mtry", "ntree", "TP", "FP", "TN", "FN", "ROC", names(env.explanatory.var))
forestdep.spplist[,namevector] <- as.numeric(NA)



#model fitting and cross-validated predictions for Random Forest
#setting the algorithm to use in caret by defining a list that contains 
#a number of custom named elements that the caret package looks for:
customRF <- list(type = "Classification", library = "randomForest", loop = NULL)

customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))

customRF$grid <- function(x, y, len = NULL, search = "grid") {}

customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}

customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)

customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")

customRF$sort <- function(x) x[order(x[,1]),]

customRF$levels <- function(x) x$classes


rm(list=ls()[!ls() %in% c("std.proj", "pgm.shp", "stm.shp", "sel.var.df", "env.explanatory.var", 
                         "forestdep.spplist", "sppdata.final", "UPFls.invcore", "customRF")])
gc()



#' species distribution modelling
#i <- as.character(forestdep.spplist$Binomial[269])
for (i in forestdep.spplist$Binomial) {
  
  #spliting the job
  if(forestdep.spplist[forestdep.spplist$Binomial==i, "job"] != 1) next
  
  #checking if model is done
  if(forestdep.spplist[forestdep.spplist$Binomial==i, "Done"] == TRUE) next
  
  cat("\n> starting", i, "<\n")
  #starting modeling procedure
  ##data preparation
  occur.p <- as.data.frame(sppdata.final %>% 
                           filter(Binomial==i) %>% 
                           dplyr::select(Binomial, Region, Catchment, Longitude, Latitude) %>% 
                           mutate(response=1))
  
  occur <- occur.p
  
        #pseudo-absence dataset
        #cells with less than 20% of undegraded primary forest in landscape
        #for each catchment -- 10km buffer
        for (c in unique(occur.p$Catchment)) {
        ptx <- SpatialPoints(occur.p[occur.p$Catchment==c,c("Longitude", "Latitude")])
        buffer <- gBuffer(ptx, width=0.09009009)
        pa.area <- mask(UPFls.invcore, buffer)
        pa.data <- as.data.frame(xyFromCell(pa.area, cell = which(pa.area[]==1), spatial = F))
        names(pa.data) <- c("Longitude", "Latitude")
        if(nrow(pa.data)==0) next
        pa.data$response <- 0
        set.seed(999)
        occur.a <- as.data.frame(pa.data %>% sample_n(as.numeric(length(ptx)*5), replace = T) %>% 
                                   mutate(Binomial=i,
                                          Catchment=c,
                                          Region=as.character(unique(occur[occur$Catchment==c, "Region"]))))
        occur <- rbind(occur, occur.a)
        }
        cat("\n> pseudo-absence for ", i, "created <\n")
  
  #extracting values from environmental explanatory variables
  env.var <- extract(env.explanatory.var, SpatialPoints(occur[,c("Longitude", "Latitude")]))
  occur <- cbind(occur, env.var)
  occur$Catchment <- as.factor(occur$Catchment)
  occur$response <- as.factor(occur$response)
  levels(occur$response) <- c("Abs","Pres")
  
  #train model
  fold<-groupKFold(occur$Catchment, k = length(unique(occur$Catchment)))
  
  ctrl <- trainControl(method="repeatedcv", number=9, repeats=30,
                       index = fold,
                       classProbs=T, summaryFunction = twoClassSummary,
                       savePredictions = T)  
  
  tunegrid <- expand.grid(.mtry=c(3:14), .ntree=seq(100,2000, 50))
  
  set.seed(999)
  mod.rf.fn1 <- train(response ~ distriver + distroad + DPFpx + edgedist + edgels + edgepx + elevation +
                        meanprecips + meantemps + MFls + SFAgepx + SFls + TSDls + UPFpx, data = occur,
                      method=customRF, importance=T, metric="ROC", tuneGrid=tunegrid, trControl=ctrl) #, ntree=ntree
    
  
  forestdep.spplist[forestdep.spplist$Binomial==i,"mtry"] <- mod.rf.fn1$bestTune[,"mtry"]
  forestdep.spplist[forestdep.spplist$Binomial==i,"ntree"] <- mod.rf.fn1$bestTune[,"ntree"]
  forestdep.spplist[forestdep.spplist$Binomial==i,"TP"] <- mod.rf.fn1$finalModel$confusion[2,2]
  forestdep.spplist[forestdep.spplist$Binomial==i,"FP"] <- mod.rf.fn1$finalModel$confusion[2,1]
  forestdep.spplist[forestdep.spplist$Binomial==i,"TN"] <- mod.rf.fn1$finalModel$confusion[1,1]
  forestdep.spplist[forestdep.spplist$Binomial==i,"FN"] <- mod.rf.fn1$finalModel$confusion[1,2]
  forestdep.spplist[forestdep.spplist$Binomial==i,"ROC"] <- mod.rf.fn1$results[mod.rf.fn1$results$mtry==forestdep.spplist[forestdep.spplist$Binomial==i,"mtry"] &
                                                                                 mod.rf.fn1$results$ntree==forestdep.spplist[forestdep.spplist$Binomial==i,"ntree"], "ROC"]
  
  forestdep.spplist <- forestdep.spplist %>% rows_update(
                              mod.rf.fn1$finalModel$importance %>% 
                                               as_tibble(rownames="VAR") %>% 
                                               dplyr::select(VAR, MeanDecreaseAccuracy) %>% 
                                               arrange(desc(MeanDecreaseAccuracy)) %>% 
                                               mutate(Binomial=i,
                                                      MeanDecreaseAccuracy=rescale(MeanDecreaseAccuracy,
                                                                                   to=c(0,1))) %>% 
                                               pivot_wider(names_from=VAR, values_from = MeanDecreaseAccuracy))
                    
  cat("\n> RF model framework for ", i, " completed; starting projections <\n")
  #projections
  
  if (any(occur$Region=="PGM")) {
    #PGM 2010 real
    pgm.2010real.raster.list <- list.files("rasters/PGM/2010_real/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2010real <- stack(pgm.2010real.raster.list)
    pgm.2010real.rf <- pgm.2010real[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_pgm.2010real <- predict(pgm.2010real.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_pgm.2010real, paste0("models.output/biodiversity.maps/PGM/2010_real/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #PGM 2020 real
    pgm.2020real.raster.list <- list.files("rasters/PGM/2020_real/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020real <- stack(pgm.2020real.raster.list)
    pgm.2020real.rf <- pgm.2020real[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_pgm.2020real <- predict(pgm.2020real.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_pgm.2020real, paste0("models.output/biodiversity.maps/PGM/2020_real/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #PGM 2020 avoiddeforest (all)
    pgm.2020_avoiddeforest.raster.list <- list.files("rasters/PGM/2020_avoiddeforest/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020_avoiddeforest <- stack(pgm.2020_avoiddeforest.raster.list)
    pgm.2020_avoiddeforest.rf <- pgm.2020_avoiddeforest[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_pgm.2020_avoiddeforest <- predict(pgm.2020_avoiddeforest.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_pgm.2020_avoiddeforest, paste0("models.output/biodiversity.maps/PGM/2020_avoiddeforest/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #PGM 2020 avoiddeforest (primary forest only)
    pgm.2020_avoiddeforest2.raster.list <- list.files("rasters/PGM/2020_avoiddeforest2/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020_avoiddeforest2 <- stack(pgm.2020_avoiddeforest2.raster.list)
    pgm.2020_avoiddeforest2.rf <- pgm.2020_avoiddeforest2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_pgm.2020_avoiddeforest2 <- predict(pgm.2020_avoiddeforest2.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_pgm.2020_avoiddeforest2, paste0("models.output/biodiversity.maps/PGM/2020_avoiddeforest2/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #PGM 2020 avoiddegrad (all)
    pgm.2020_avoiddegrad.raster.list <- list.files("rasters/PGM/2020_avoiddegrad/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020_avoiddegrad <- stack(pgm.2020_avoiddegrad.raster.list)
    pgm.2020_avoiddegrad.rf <- pgm.2020_avoiddegrad[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_pgm.2020_avoiddegrad <- predict(pgm.2020_avoiddegrad.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_pgm.2020_avoiddegrad, paste0("models.output/biodiversity.maps/PGM/2020_avoiddegrad/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #PGM 2020 avoiddegrad (primary forest only)
    pgm.2020_avoiddegrad2.raster.list <- list.files("rasters/PGM/2020_avoiddegrad2/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020_avoiddegrad2 <- stack(pgm.2020_avoiddegrad2.raster.list)
    pgm.2020_avoiddegrad2.rf <- pgm.2020_avoiddegrad2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_pgm.2020_avoiddegrad2 <- predict(pgm.2020_avoiddegrad2.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_pgm.2020_avoiddegrad2, paste0("models.output/biodiversity.maps/PGM/2020_avoiddegrad2/", i, ".tif"), format = "GTiff", overwrite = T)
    #writeRaster(proj_pgm.2020_avoiddegrad2.conbywm, paste0("models.output/biodiversity.maps/PGM/2020_avoiddegrad2/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #PGM 2020 restor_wo_avoid
    pgm.2020_restor_wo_avoid.raster.list <- list.files("rasters/PGM/2020_restor_wo_avoid/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020_restor_wo_avoid <- stack(pgm.2020_restor_wo_avoid.raster.list)
    pgm.2020_restor_wo_avoid.rf <- pgm.2020_restor_wo_avoid[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_pgm.2020_restor_wo_avoid <- predict(pgm.2020_restor_wo_avoid.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_pgm.2020_restor_wo_avoid, paste0("models.output/biodiversity.maps/PGM/2020_restor_wo_avoid/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #PGM 2020 avoidboth (all)
    pgm.2020_avoidboth.raster.list <- list.files("rasters/PGM/2020_avoidboth/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020_avoidboth <- stack(pgm.2020_avoidboth.raster.list)
    pgm.2020_avoidboth.rf <- pgm.2020_avoidboth[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_pgm.2020_avoidboth <- predict(pgm.2020_avoidboth.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_pgm.2020_avoidboth, paste0("models.output/biodiversity.maps/PGM/2020_avoidboth/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #PGM 2020 avoidboth (primary forest only)
    pgm.2020_avoidboth2.raster.list <- list.files("rasters/PGM/2020_avoidboth2/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020_avoidboth2 <- stack(pgm.2020_avoidboth2.raster.list)
    pgm.2020_avoidboth2.rf <- pgm.2020_avoidboth2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_pgm.2020_avoidboth2 <- predict(pgm.2020_avoidboth2.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_pgm.2020_avoidboth2, paste0("models.output/biodiversity.maps/PGM/2020_avoidboth2/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #PGM 2020 restor_n_avoiddeforest (all)
    pgm.2020_restor_n_avoiddeforest.raster.list <- list.files("rasters/PGM/2020_restor_n_avoiddeforest/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020_restor_n_avoiddeforest <- stack(pgm.2020_restor_n_avoiddeforest.raster.list)
    pgm.2020_restor_n_avoiddeforest.rf <- pgm.2020_restor_n_avoiddeforest[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_pgm.2020_restor_n_avoiddeforest <- predict(pgm.2020_restor_n_avoiddeforest.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_pgm.2020_restor_n_avoiddeforest, paste0("models.output/biodiversity.maps/PGM/2020_restor_n_avoiddeforest/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #PGM 2020 restor_n_avoiddeforest (primary forest only)
    pgm.2020_restor_n_avoiddeforest2.raster.list <- list.files("rasters/PGM/2020_restor_n_avoiddeforest2/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020_restor_n_avoiddeforest2 <- stack(pgm.2020_restor_n_avoiddeforest2.raster.list)
    pgm.2020_restor_n_avoiddeforest2.rf <- pgm.2020_restor_n_avoiddeforest2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_pgm.2020_restor_n_avoiddeforest2 <- predict(pgm.2020_restor_n_avoiddeforest2.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_pgm.2020_restor_n_avoiddeforest2, paste0("models.output/biodiversity.maps/PGM/2020_restor_n_avoiddeforest2/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #PGM 2020 restor_n_avoidboth (all)
    pgm.2020_restor_n_avoidboth.raster.list <- list.files("rasters/PGM/2020_restor_n_avoidboth/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020_restor_n_avoidboth <- stack(pgm.2020_restor_n_avoidboth.raster.list)
    pgm.2020_restor_n_avoidboth.rf <- pgm.2020_restor_n_avoidboth[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_pgm.2020_restor_n_avoidboth <- predict(pgm.2020_restor_n_avoidboth.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_pgm.2020_restor_n_avoidboth, paste0("models.output/biodiversity.maps/PGM/2020_restor_n_avoidboth/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #PGM 2020 restor_n_avoidboth (primary forest only)
    pgm.2020_restor_n_avoidboth2.raster.list <- list.files("rasters/PGM/2020_restor_n_avoidboth2/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020_restor_n_avoidboth2 <- stack(pgm.2020_restor_n_avoidboth2.raster.list)
    pgm.2020_restor_n_avoidboth2.rf <- pgm.2020_restor_n_avoidboth2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_pgm.2020_restor_n_avoidboth2 <- predict(pgm.2020_restor_n_avoidboth2.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_pgm.2020_restor_n_avoidboth2, paste0("models.output/biodiversity.maps/PGM/2020_restor_n_avoidboth2/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    #
    
  }
  
  
  if (any(occur$Region=="STM")) {
    #STM 2010 real
    stm.2010real.raster.list <- list.files("rasters/STM/2010_real/", pattern = ".tif", full.names = T, recursive = T)
    stm.2010real <- stack(stm.2010real.raster.list)
    stm.2010real.rf <- stm.2010real[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_stm.2010real <- predict(stm.2010real.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_stm.2010real, paste0("models.output/biodiversity.maps/STM/2010_real/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #STM 2020 real
    stm.2020real.raster.list <- list.files("rasters/STM/2020_real/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020real <- stack(stm.2020real.raster.list)
    stm.2020real.rf <- stm.2020real[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_stm.2020real <- predict(stm.2020real.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_stm.2020real, paste0("models.output/biodiversity.maps/STM/2020_real/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #STM 2020 avoiddeforest (all)
    stm.2020_avoiddeforest.raster.list <- list.files("rasters/STM/2020_avoiddeforest/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020_avoiddeforest <- stack(stm.2020_avoiddeforest.raster.list)
    stm.2020_avoiddeforest.rf <- stm.2020_avoiddeforest[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_stm.2020_avoiddeforest <- predict(stm.2020_avoiddeforest.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_stm.2020_avoiddeforest, paste0("models.output/biodiversity.maps/STM/2020_avoiddeforest/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #STM 2020 avoiddeforest (primary forest only)
    stm.2020_avoiddeforest2.raster.list <- list.files("rasters/STM/2020_avoiddeforest2/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020_avoiddeforest2 <- stack(stm.2020_avoiddeforest2.raster.list)
    stm.2020_avoiddeforest2.rf <- stm.2020_avoiddeforest2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_stm.2020_avoiddeforest2 <- predict(stm.2020_avoiddeforest2.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_stm.2020_avoiddeforest2, paste0("models.output/biodiversity.maps/STM/2020_avoiddeforest2/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #STM 2020 avoiddegrad (all)
    stm.2020_avoiddegrad.raster.list <- list.files("rasters/STM/2020_avoiddegrad/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020_avoiddegrad <- stack(stm.2020_avoiddegrad.raster.list)
    stm.2020_avoiddegrad.rf <- stm.2020_avoiddegrad[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_stm.2020_avoiddegrad <- predict(stm.2020_avoiddegrad.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_stm.2020_avoiddegrad, paste0("models.output/biodiversity.maps/STM/2020_avoiddegrad/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #STM 2020 avoiddegrad (primary forest only)
    stm.2020_avoiddegrad2.raster.list <- list.files("rasters/STM/2020_avoiddegrad2/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020_avoiddegrad2 <- stack(stm.2020_avoiddegrad2.raster.list)
    stm.2020_avoiddegrad2.rf <- stm.2020_avoiddegrad2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_stm.2020_avoiddegrad2 <- predict(stm.2020_avoiddegrad2.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_stm.2020_avoiddegrad2, paste0("models.output/biodiversity.maps/STM/2020_avoiddegrad2/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #STM 2020 restor_wo_avoid
    stm.2020_restor_wo_avoid.raster.list <- list.files("rasters/STM/2020_restor_wo_avoid/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020_restor_wo_avoid <- stack(stm.2020_restor_wo_avoid.raster.list)
    stm.2020_restor_wo_avoid.rf <- stm.2020_restor_wo_avoid[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_stm.2020_restor_wo_avoid <- predict(stm.2020_restor_wo_avoid.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_stm.2020_restor_wo_avoid, paste0("models.output/biodiversity.maps/STM/2020_restor_wo_avoid/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #STM 2020 avoidboth (all)
    stm.2020_avoidboth.raster.list <- list.files("rasters/STM/2020_avoidboth/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020_avoidboth <- stack(stm.2020_avoidboth.raster.list)
    stm.2020_avoidboth.rf <- stm.2020_avoidboth[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_stm.2020_avoidboth <- predict(stm.2020_avoidboth.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_stm.2020_avoidboth, paste0("models.output/biodiversity.maps/STM/2020_avoidboth/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #STM 2020 avoidboth (primary forest only)
    stm.2020_avoidboth2.raster.list <- list.files("rasters/STM/2020_avoidboth2/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020_avoidboth2 <- stack(stm.2020_avoidboth2.raster.list)
    stm.2020_avoidboth2.rf <- stm.2020_avoidboth2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_stm.2020_avoidboth2 <- predict(stm.2020_avoidboth2.rf, mod.rf.fn1, type="prob", index=2)
    
    
    
    writeRaster(mod.rf.proj_stm.2020_avoidboth2, paste0("models.output/biodiversity.maps/STM/2020_avoidboth2/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #STM 2020 restor_n_avoiddeforest (all)
    stm.2020_restor_n_avoiddeforest.raster.list <- list.files("rasters/STM/2020_restor_n_avoiddeforest/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020_restor_n_avoiddeforest <- stack(stm.2020_restor_n_avoiddeforest.raster.list)
    stm.2020_restor_n_avoiddeforest.rf <- stm.2020_restor_n_avoiddeforest[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_stm.2020_restor_n_avoiddeforest <- predict(stm.2020_restor_n_avoiddeforest.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_stm.2020_restor_n_avoiddeforest, paste0("models.output/biodiversity.maps/STM/2020_restor_n_avoiddeforest/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #STM 2020 restor_n_avoiddeforest (primary forest only)
    stm.2020_restor_n_avoiddeforest2.raster.list <- list.files("rasters/STM/2020_restor_n_avoiddeforest2/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020_restor_n_avoiddeforest2 <- stack(stm.2020_restor_n_avoiddeforest2.raster.list)
    stm.2020_restor_n_avoiddeforest2.rf <- stm.2020_restor_n_avoiddeforest2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_stm.2020_restor_n_avoiddeforest2 <- predict(stm.2020_restor_n_avoiddeforest2.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_stm.2020_restor_n_avoiddeforest2, paste0("models.output/biodiversity.maps/STM/2020_restor_n_avoiddeforest2/", i, ".tif"), format = "GTiff", overwrite = T)
    #writeRaster(proj_stm.2020_restor_n_avoiddeforest2.conbywm, paste0("models.output/biodiversity.maps/STM/2020_restor_n_avoiddeforest2/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #STM 2020 restor_n_avoidboth (all)
    stm.2020_restor_n_avoidboth.raster.list <- list.files("rasters/STM/2020_restor_n_avoidboth/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020_restor_n_avoidboth <- stack(stm.2020_restor_n_avoidboth.raster.list)
    stm.2020_restor_n_avoidboth.rf <- stm.2020_restor_n_avoidboth[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_stm.2020_restor_n_avoidboth <- predict(stm.2020_restor_n_avoidboth.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_stm.2020_restor_n_avoidboth, paste0("models.output/biodiversity.maps/STM/2020_restor_n_avoidboth/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    
    
    
    #STM 2020 restor_n_avoidboth (primary forest only)
    stm.2020_restor_n_avoidboth2.raster.list <- list.files("rasters/STM/2020_restor_n_avoidboth2/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020_restor_n_avoidboth2 <- stack(stm.2020_restor_n_avoidboth2.raster.list)
    stm.2020_restor_n_avoidboth2.rf <- stm.2020_restor_n_avoidboth2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    
    set.seed(999)
    mod.rf.proj_stm.2020_restor_n_avoidboth2 <- predict(stm.2020_restor_n_avoidboth2.rf, mod.rf.fn1, type="prob", index=2)
    
    
    writeRaster(mod.rf.proj_stm.2020_restor_n_avoidboth2, paste0("models.output/biodiversity.maps/STM/2020_restor_n_avoidboth2/", i, ".tif"), format = "GTiff", overwrite = T)
    
    #
    #
    
  }
  
  forestdep.spplist[forestdep.spplist$Binomial==i,"Done"] <- "TRUE"
  write.csv(forestdep.spplist, "models.output/biodiversity.maps/preliminar_species_summary.csv", row.names = F)
  
  cat("\n> ", i, "concluded <\n")
  rm(list=ls()[!ls() %in% c("std.proj", "pgm.shp", "stm.shp", "sel.var.df", "env.explanatory.var", 
                            "forestdep.spplist", "sppdata.final", "UPFls.invcore", "customRF")])
  gc()
}
  
  
  

#
#



#  biodiversity benefit: for species with less than 5 records ==================
#' % of UPFls

#forestdep.spplist <- read.csv("data/updated_species_summary_edby_visual_inspection.csv")
#sppdata.final <- read.csv("data/presence_records.csv")

#i <- as.character(forestdep.spplist$Binomial[2])
for (i in forestdep.spplist$Binomial) {
  
  #spliting the job
  if(forestdep.spplist[forestdep.spplist$Binomial==i, "job"] != 13) next
  
  #checking if model is done
  if(forestdep.spplist[forestdep.spplist$Binomial==i, "Done"] != FALSE) next
  
  #starting modeling procedure
  occur <- sppdata.final[sppdata.final$Binomial==i,]
  
  #  
  
if (any(occur$Region=="PGM")) {
  
  pgm.occur <- occur[occur$Region=="PGM",]
  pgm.occur <- SpatialPoints(occur[,c("Longitude", "Latitude")])
  pres.bkg <- gBuffer(pgm.occur, width = 0.09)
  
  # 2010 real
  pgm.2010real <- raster("rasters/PGM/2010_real/UPFls.tif")
  pgm.pres.2010real <- mask(pgm.2010real, pres.bkg, updatevalue=0)
  
  writeRaster(pgm.pres.2010real, paste0("models.output/biodiversity.maps/PGM/2010_real/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 real
  pgm.2020real <- raster("rasters/PGM/2020_real/UPFls.tif")
  pgm.pres.2020real <- mask(pgm.2020real, pres.bkg, updatevalue=0)
  
  writeRaster(pgm.pres.2020real, paste0("models.output/biodiversity.maps/PGM/2010_real/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 avoiddeforest (all)
  pgm.2020_avoiddeforest <- raster("rasters/PGM/2020_avoiddeforest/UPFls.tif")
  pgm.pres.2020_avoiddeforest <- mask(pgm.2020_avoiddeforest, pres.bkg, updatevalue=0)
  
  writeRaster(pgm.pres.2020_avoiddeforest, paste0("models.output/biodiversity.maps/PGM/2020_avoiddeforest/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 avoiddeforest (primary forest only)
  pgm.2020_avoiddeforest2 <- raster("rasters/PGM/2020_avoiddeforest2/UPFls.tif")
  pgm.pres.2020_avoiddeforest2 <- mask(pgm.2020_avoiddeforest2, pres.bkg, updatevalue=0)
  
  writeRaster(pgm.pres.2020_avoiddeforest2, paste0("models.output/biodiversity.maps/PGM/2020_avoiddeforest2/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 avoiddegrad (all)
  pgm.2020_avoiddegrad <- raster("rasters/PGM/2020_avoiddegrad/UPFls.tif")
  pgm.pres.2020_avoiddegrad <- mask(pgm.2020_avoiddegrad, pres.bkg, updatevalue=0)
  
  writeRaster(pgm.pres.2020_avoiddegrad, paste0("models.output/biodiversity.maps/PGM/2020_avoiddegrad/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 avoiddegrad (primary forest only)
  pgm.2020_avoiddegrad2 <- raster("rasters/PGM/2020_avoiddegrad2/UPFls.tif")
  pgm.pres.2020_avoiddegrad2 <- mask(pgm.2020_avoiddegrad2, pres.bkg, updatevalue=0)
  
  writeRaster(pgm.pres.2020_avoiddegrad2, paste0("models.output/biodiversity.maps/PGM/2020_avoiddegrad2/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 restor_wo_avoid
  pgm.2020_restor_wo_avoid <- raster("rasters/PGM/2020_restor_wo_avoid/UPFls.tif")
  pgm.pres.2020_restor_wo_avoid <- mask(pgm.2020_restor_wo_avoid, pres.bkg, updatevalue=0)
  
  writeRaster(pgm.pres.2020_restor_wo_avoid, paste0("models.output/biodiversity.maps/PGM/2020_restor_wo_avoid/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 avoidboth (all)
  pgm.2020_avoidboth <- raster("rasters/PGM/2020_avoidboth/UPFls.tif")
  pgm.pres.2020_avoidboth <- mask(pgm.2020_avoidboth, pres.bkg, updatevalue=0)
  
  writeRaster(pgm.pres.2020_avoidboth, paste0("models.output/biodiversity.maps/PGM/2020_avoidboth/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 avoidboth (primary forest only)
  pgm.2020_avoidboth2 <- raster("rasters/PGM/2020_avoidboth2/UPFls.tif")
  pgm.pres.2020_avoidboth2 <- mask(pgm.2020_avoidboth2, pres.bkg, updatevalue=0)
  
  writeRaster(pgm.pres.2020_avoidboth2, paste0("models.output/biodiversity.maps/PGM/2020_avoidboth2/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 restor_n_avoiddeforest (all)
  pgm.2020_restor_n_avoiddeforest <- raster("rasters/PGM/2020_restor_n_avoiddeforest/UPFls.tif")
  pgm.pres.2020_restor_n_avoiddeforest <- mask(pgm.2020_restor_n_avoiddeforest, pres.bkg, updatevalue=0)
  
  writeRaster(pgm.pres.2020_restor_n_avoiddeforest, paste0("models.output/biodiversity.maps/PGM/2020_restor_n_avoiddeforest/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 restor_n_avoiddeforest (primary forest only)
  pgm.2020_restor_n_avoiddeforest2 <- raster("rasters/PGM/2020_restor_n_avoiddeforest2/UPFls.tif")
  pgm.pres.2020_restor_n_avoiddeforest2 <- mask(pgm.2020_restor_n_avoiddeforest2, pres.bkg, updatevalue=0)
  
  writeRaster(pgm.pres.2020_restor_n_avoiddeforest2, paste0("models.output/biodiversity.maps/PGM/2020_restor_n_avoiddeforest2/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 restor_n_avoidboth (all)
  pgm.2020_restor_n_avoidboth <- raster("rasters/PGM/2020_restor_n_avoidboth/UPFls.tif")
  pgm.pres.2020_restor_n_avoidboth <- mask(pgm.2020_restor_n_avoidboth, pres.bkg, updatevalue=0)
  
  writeRaster(pgm.pres.2020_restor_n_avoidboth, paste0("models.output/biodiversity.maps/PGM/2020_restor_n_avoidboth/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 restor_n_avoidboth (primary forest only)
  pgm.2020_restor_n_avoidboth2 <- raster("rasters/PGM/2020_restor_n_avoidboth2/UPFls.tif")
  pgm.pres.2020_restor_n_avoidboth2 <- mask(pgm.2020_restor_n_avoidboth2, pres.bkg, updatevalue=0)
  
  writeRaster(pgm.pres.2020_restor_n_avoidboth2, paste0("models.output/biodiversity.maps/PGM/2020_restor_n_avoidboth2/", i, ".tif"), format = "GTiff", overwrite = T)

  
  rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "occur"))])
  gc()
  #
  
}
  
  
  #  

if (any(occur$Region=="STM")) {
  
  stm.occur <- occur[occur$Region=="STM",]
  stm.occur <- SpatialPoints(occur[,c("Longitude", "Latitude")])
  pres.bkg <- gBuffer(stm.occur, width = 0.09)
  
  # 2010 real
  stm.2010real <- raster("rasters/STM/2010_real/UPFls.tif")
  stm.pres.2010real <- mask(stm.2010real, pres.bkg, updatevalue=0)
  
  writeRaster(stm.pres.2010real, paste0("models.output/biodiversity.maps/STM/2010_real/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 real
  stm.2020real <- raster("rasters/STM/2020_real/UPFls.tif")
  stm.pres.2020real <- mask(stm.2020real, pres.bkg, updatevalue=0)
  
  writeRaster(stm.pres.2020real, paste0("models.output/biodiversity.maps/STM/2010_real/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 avoiddeforest (all)
  stm.2020_avoiddeforest <- raster("rasters/STM/2020_avoiddeforest/UPFls.tif")
  stm.pres.2020_avoiddeforest <- mask(stm.2020_avoiddeforest, pres.bkg, updatevalue=0)
  
  writeRaster(stm.pres.2020_avoiddeforest, paste0("models.output/biodiversity.maps/STM/2020_avoiddeforest/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 avoiddeforest (primary forest only)
  stm.2020_avoiddeforest2 <- raster("rasters/STM/2020_avoiddeforest2/UPFls.tif")
  stm.pres.2020_avoiddeforest2 <- mask(stm.2020_avoiddeforest2, pres.bkg, updatevalue=0)
  
  writeRaster(stm.pres.2020_avoiddeforest2, paste0("models.output/biodiversity.maps/STM/2020_avoiddeforest2/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 avoiddegrad (all)
  stm.2020_avoiddegrad <- raster("rasters/STM/2020_avoiddegrad/UPFls.tif")
  stm.pres.2020_avoiddegrad <- mask(stm.2020_avoiddegrad, pres.bkg, updatevalue=0)
  
  writeRaster(stm.pres.2020_avoiddegrad, paste0("models.output/biodiversity.maps/STM/2020_avoiddegrad/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 avoiddegrad (primary forest only)
  stm.2020_avoiddegrad2 <- raster("rasters/STM/2020_avoiddegrad2/UPFls.tif")
  stm.pres.2020_avoiddegrad2 <- mask(stm.2020_avoiddegrad2, pres.bkg, updatevalue=0)
  
  writeRaster(stm.pres.2020_avoiddegrad2, paste0("models.output/biodiversity.maps/STM/2020_avoiddegrad2/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 restor_wo_avoid
  stm.2020_restor_wo_avoid <- raster("rasters/STM/2020_restor_wo_avoid/UPFls.tif")
  stm.pres.2020_restor_wo_avoid <- mask(stm.2020_restor_wo_avoid, pres.bkg, updatevalue=0)
  
  writeRaster(stm.pres.2020_restor_wo_avoid, paste0("models.output/biodiversity.maps/STM/2020_restor_wo_avoid/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 avoidboth (all)
  stm.2020_avoidboth <- raster("rasters/STM/2020_avoidboth/UPFls.tif")
  stm.pres.2020_avoidboth <- mask(stm.2020_avoidboth, pres.bkg, updatevalue=0)
  
  writeRaster(stm.pres.2020_avoidboth, paste0("models.output/biodiversity.maps/STM/2020_avoidboth/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 avoidboth (primary forest only)
  stm.2020_avoidboth2 <- raster("rasters/STM/2020_avoidboth2/UPFls.tif")
  stm.pres.2020_avoidboth2 <- mask(stm.2020_avoidboth2, pres.bkg, updatevalue=0)
  
  writeRaster(stm.pres.2020_avoidboth2, paste0("models.output/biodiversity.maps/STM/2020_avoidboth2/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 restor_n_avoiddeforest (all)
  stm.2020_restor_n_avoiddeforest <- raster("rasters/STM/2020_restor_n_avoiddeforest/UPFls.tif")
  stm.pres.2020_restor_n_avoiddeforest <- mask(stm.2020_restor_n_avoiddeforest, pres.bkg, updatevalue=0)
  
  writeRaster(stm.pres.2020_restor_n_avoiddeforest, paste0("models.output/biodiversity.maps/STM/2020_restor_n_avoiddeforest/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 restor_n_avoiddeforest (primary forest only)
  stm.2020_restor_n_avoiddeforest2 <- raster("rasters/STM/2020_restor_n_avoiddeforest2/UPFls.tif")
  stm.pres.2020_restor_n_avoiddeforest2 <- mask(stm.2020_restor_n_avoiddeforest2, pres.bkg, updatevalue=0)
  
  writeRaster(stm.pres.2020_restor_n_avoiddeforest2, paste0("models.output/biodiversity.maps/STM/2020_restor_n_avoiddeforest2/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 restor_n_avoidboth (all)
  stm.2020_restor_n_avoidboth <- raster("rasters/STM/2020_restor_n_avoidboth/UPFls.tif")
  stm.pres.2020_restor_n_avoidboth <- mask(stm.2020_restor_n_avoidboth, pres.bkg, updatevalue=0)
  
  writeRaster(stm.pres.2020_restor_n_avoidboth, paste0("models.output/biodiversity.maps/STM/2020_restor_n_avoidboth/", i, ".tif"), format = "GTiff", overwrite = T)
  
  # 2020 restor_n_avoidboth (primary forest only)
  stm.2020_restor_n_avoidboth2 <- raster("rasters/STM/2020_restor_n_avoidboth2/UPFls.tif")
  stm.pres.2020_restor_n_avoidboth2 <- mask(stm.2020_restor_n_avoidboth2, pres.bkg, updatevalue=0)
  
  writeRaster(stm.pres.2020_restor_n_avoidboth2, paste0("models.output/biodiversity.maps/STM/2020_restor_n_avoidboth2/", i, ".tif"), format = "GTiff", overwrite = T)
  
  
  rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "occur"))])
  gc()
  #
  
}
  
  forestdep.spplist[forestdep.spplist$Binomial==i,"Done"] <- "TRUE"
  write.csv(forestdep.spplist, "models.output/biodiversity.maps/preliminar_species_summary.csv", row.names = F)
  
  
  cat("\n> ", i, "concluded <\n")
}


# importing species list and methods and evaluation metrics data frame
forestdep.spplist <- read.csv("data/species_summary.csv")


#  biodiversity benefit:  adding conservation value ============================
#' bird conservation value is inverse occurrence area size
#' scaled from 0 [the biggest] to 1 [the smallest]
#' based on birdlife_v2017b shapefiles
bird.distribution.shp <- readOGR(dsn="~/GIS/bird_distr_BirdLifeInternational2016", layer="Birds_of_Amazonia")
head(bird.distribution.shp@data)
# adding variable to match with species list dataframe
bird.distribution.shp@data$Binomial <- gsub(" ", "", bird.distribution.shp@data$SCINAME)

# adding shape area from birdlife_v2017b shapefiles
forestdep.spplist <- forestdep.spplist %>% left_join(bird.distribution.shp@data[,c("Binomial", "Shape_Area")]) %>%
                        group_by(Binomial) %>% mutate(Shape_Area = sum(Shape_Area)) %>% distinct(Binomial, .keep_all = T)

# scaling
forestdep.spplist$Shape_Area_scaled <- rescale(forestdep.spplist$Shape_Area, to = c(.999,.001))

#

rm(bird.distribution.shp)


#' tree conservation value is inverse occurrence area size
#' scaled from 0 [the biggest] to 1 [the smallest]
#' based on birdlife_v2017b shapefiles
pgm.treedata <- read.csv("data/raw/Flora.composition.and.biomass_PGM_Erika_23.01.2013.csv")
stm.treedata <- read.csv("data/raw/Flora.composition.and.biomass_STM_Erika_23.01.2013.csv")
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

rm(treedata)


# saving
write.csv(forestdep.spplist, "data/species_summary.csv")
#forestdep.spplist <- read.csv("data/species_summary.csv")



# syntethic map: biodiversity benefit ==========================================
dir.create("models.output/biodiversity.benefits", recursive = T)


scenarios <- c("PGM/2010_real/", "PGM/2020_real/", "PGM/2020_avoiddegrad/", "PGM/2020_avoiddegrad2/", 
               "PGM/2020_avoiddeforest/", "PGM/2020_avoiddeforest2/", "PGM/2020_avoidboth/", "PGM/2020_avoidboth2/",
               "PGM/2020_restor_wo_avoid/", "PGM/2020_restor_n_avoiddeforest/", "PGM/2020_restor_n_avoiddeforest2/",
               "PGM/2020_restor_n_avoidboth/", "PGM/2020_restor_n_avoidboth2",
               "STM/2010_real/", "STM/2020_real/", "STM/2020_avoiddegrad/", "STM/2020_avoiddegrad2/", 
               "STM/2020_avoiddeforest/", "STM/2020_avoiddeforest2/", "STM/2020_avoidboth/", "STM/2020_avoidboth2/",
               "STM/2020_restor_wo_avoid/", "STM/2020_restor_n_avoiddeforest/", "STM/2020_restor_n_avoiddeforest2/",
               "STM/2020_restor_n_avoidboth/", "STM/2020_restor_n_avoidboth2")


for (s in scenarios) {
  
  cat("\n>>>>>>>>>>>>> STARTING SCENARIO", s, " <<<<<<<<<<<<<<\n")
  maps.list <- list.files(paste0("models.output/biodiversity.maps/", s), pattern = ".tif", full.names = T, recursive = T)
  maps.list <- grep(paste(forestdep.spplist$Binomial, collapse = "|"), maps.list, value = T)
  
  #biodiversity.maps <- stack(maps.list)
  
  #conservation.value <- forestdep.spplist[grep(paste(names(biodiversity.maps), collapse = "|"), forestdep.spplist$Binomial), "Shape_Area_scaled"]
  
  i=0
  biodiversity.benefit <- NULL
  for (m in forestdep.spplist$Binomial) {
    
    # Load the raster
    spx <- try(raster(grep(m, maps.list, value = T)), silent = T)
    
    # Multiply the raster by the corresponding conservation value
    spx <- if(class(spx) %in% 'try-error') { next } 
    else {spx * forestdep.spplist[forestdep.spplist$Binomial==m,"Shape_Area_scaled"]}
    
    # Add the weighted raster to the biodiversity benefit raster
    if (is.null(biodiversity.benefit)) {
      biodiversity.benefit <- spx
    } else {
      biodiversity.benefit <- biodiversity.benefit + spx
    }
    
    
    # Calculate progress percentage
    i=i+1
    progress <- (i / length(maps.list)) * 100
    
    # Check if progress is a multiple of 10
    if (round(progress %% 10, 0) == 0) {
      cat("Processing scenario", s, ":", progress, "% complete\n")
    }
  }
  
  rm(spx); gc()
  
  writeRaster(biodiversity.benefit, paste0("D:/projetos/lancaster/conservation_opportunities/models.output/biodiversity.benefits/", gsub("/", "_", s), "biodiversity_benefit.tif"), format="GTiff", overwrite=T)
  
  rm(list=ls()[!ls() %in% c("forestdep.spplist", "biodiversity.benefit", "s")]) #keeping only raster stack
  gc()
  cat("\n>>>>>>>>>>>>> SCENARIO", s, "DONE <<<<<<<<<<<<<<\n")
  
  
}

#
#





# carbon benefit ===============================================================

dir.create("models.output/carbon.benefits2", recursive = T)

## transect data
transectdata <- read.csv("data/raw/RAS_transects_environment_all.csv")
#head(transectdata)
#str(transectdata)
#summary(transectdata)

#excluding Varzea transects PGM
exclude <- c("100_1", "100_4","100_7","81_12","423_2") #transect code
transectdata <- transectdata[!transectdata$Transectcode %in% exclude,]

#adding variable considering 50% of AGB as carbon stock
carbon <- transectdata %>% dplyr::select(Region:UTM_Y, LU_FT_Code, AGB_Trees10) %>% 
  mutate(carbon_stock = as.numeric(AGB_Trees10)/2) %>% drop_na()

carbon$non_zero <- ifelse(carbon$carbon_stock > 0, 1, 0)

# new var with longlat
stm.utm <- data.frame(x=carbon[carbon$Region=="STM", "UTM_X"], y=carbon[carbon$Region=="STM", "UTM_Y"]) 
coordinates(stm.utm) <- ~x+y 
#class(stm.utm)
proj4string(stm.utm) <- CRS("+proj=utm +zone=21 +south +datum=WGS84 +units=m +ellps=WGS84") 
stm.longlat <- spTransform(stm.utm, CRS("+proj=longlat +datum=WGS84"))

pgm.utm <- data.frame(x=carbon[carbon$Region=="PGM", "UTM_X"], y=carbon[carbon$Region=="PGM", "UTM_Y"]) 
coordinates(pgm.utm) <- ~x+y 
#class(pgm.utm)
proj4string(pgm.utm) <- CRS("+proj=utm +zone=23 +south +datum=WGS84 +units=m +ellps=WGS84") 
pgm.longlat <- spTransform(pgm.utm, CRS("+proj=longlat +datum=WGS84"))

longlat <- rbind(stm.longlat, pgm.longlat)

carbon$Longitude <- longlat@coords[,1]
carbon$Latitude <- longlat@coords[,2]
#plot(env.explanatory.var[["UPFpx"]])
#points(longlat)

#extracting values from environmental explanatory variables
env.var <- extract(env.explanatory.var, SpatialPoints(carbon[,c("Longitude", "Latitude")]))
carbon <- cbind(carbon, env.var)


##scaling predictors -- z-scores
#carbon <- carbon %>% mutate(distriver_z = ((distriver - mean(distriver))/sd(distriver)),
#                            distroad_z = ((distroad - mean(distroad))/sd(distroad)),
#                            DPFpx_z = ((DPFpx - mean(DPFpx))/sd(DPFpx)),
#                            edgedist_z = ((edgedist - mean(edgedist))/sd(edgedist)),
#                            edgels_z = ((edgels - mean(edgels))/sd(edgels)),
#                            elevation_z = ((elevation - mean(elevation))/sd(elevation)),
#                            meanprecips_z = ((meanprecips - mean(meanprecips))/sd(meanprecips)),
#                            meantemps_z = ((meantemps - mean(meantemps))/sd(meantemps)),
#                            SFagels_z = ((SFagels - mean(SFagels))/sd(SFagels)),
#                            SFpx_z = ((SFpx - mean(SFpx))/sd(SFpx)),
#                            TFpx_z = ((TFpx - mean(TFpx))/sd(TFpx)),
#                            UPFls_z = ((UPFls - mean(UPFls))/sd(UPFls)),
#                            TSDls_z = ((TSDls - mean(TSDls))/sd(TSDls)))
#



write.csv(carbon, "data/carbon.csv", row.names = F)
#carbon <- read.csv("data/carbon.csv")

#rm(list= ls()[!(ls() %in% c("sel.var.df", "carbon", "env.explanatory.var"))])
gc()


# distribution family
carbon.norm <- fitdist(carbon[carbon$non_zero==1,"carbon_stock"], "norm", method="mle")
summary(carbon.norm)

carbon.gamma <- fitdist(carbon[carbon$non_zero==1,"carbon_stock"], "gamma", method="mle")
summary(carbon.gamma)

carbon.expo <- fitdist(carbon[carbon$non_zero==1,"carbon_stock"], "exp", method="mle")
summary(carbon.expo)

par(mfrow=c(2,2))
plot.legend <- c("normal", "gamma", "expo")
denscomp(list(carbon.norm, carbon.gamma, carbon.expo), legendtext = plot.legend)
cdfcomp(list(carbon.norm, carbon.gamma, carbon.expo), legendtext = plot.legend)
qqcomp(list(carbon.norm, carbon.gamma, carbon.expo), legendtext = plot.legend)
ppcomp(list(carbon.norm, carbon.gamma, carbon.expo), legendtext = plot.legend)


#plot(Fn <- ecdf(carbon[carbon$non_zero==1,"carbon_stock"]))
#
#curve(pnorm(x, mean(carbon[carbon$non_zero==1,"carbon_stock"]), sd(carbon[carbon$non_zero==1,"carbon_stock"])), 
#      from = 0, to = 350, add = TRUE, col='red', lwd = 2)
#
#curve(pgamma(x, shape=1.08401844, rate=0.01162663), 
#      from = 0, to = 350, add = TRUE, col='blue', lwd = 2)




#model fiting and cross-validated predictions for Random Forest
#setting the algorithm to use in caret by defining a list that contains 
#a number of custom named elements that the caret package looks for:
customRF <- list(type = "Regression", library = "randomForest", loop = NULL)

customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))

customRF$grid <- function(x, y, len = NULL, search = "grid") {}

customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}

customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)

customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")

customRF$sort <- function(x) x[order(x[,1]),]

customRF$levels <- function(x) x$classes


#train model
#custom resampling method in caret: 'leave-pair-out-cross-validation'
#catchments are used as training set, and tested against the 2 left out
num_samps <- length(unique(carbon$Catchment))
holdout <- combn(num_samps, 2)
in_training <- apply(holdout, 2, 
                     function(x, all) all[!(all %in% x)],
                     all = 1:num_samps)

## need a more effcient way of doing this:
index <- vector(mode = "list", length = ncol(in_training))
for(i in 1:ncol(in_training)) index[[i]] <- in_training[,i]
## cosmetic:
names(index) <- caret:::prettySeq(seq(along = index))

ctrl <- trainControl(method = "cv", ## this will be ignored since we supply index
                     index = index,
                     savePredictions = T)  
#control <- trainControl(method="repeatedcv", number=10, repeats=3)
tunegrid <- expand.grid(.mtry=c(3:13), .ntree=seq(100,2000, 50))

set.seed(999)
mod.rf.fn2 <- train(carbon_stock ~ distriver + distroad + DPFpx + edgedist + edgels + edgepx + elevation +
                     meanprecips + meantemps + MFls + SFAgepx + TSDls + UPFpx, data = subset(carbon, non_zero == 1),
                method=customRF, metric="RMSE", tuneGrid=tunegrid, trControl=ctrl)

print(mod.rf.fn2)
plot(mod.rf.fn2)


# Fig S7
layout(matrix(c(1,2,3,4,1,5,6,7), 2, 4, byrow = TRUE))
varImpPlot(mod.rf.fn2$finalModel)
partialPlot(mod.rf.fn2$finalModel, carbon, TSDls)
partialPlot(mod.rf.fn2$finalModel, carbon, edgedist)
partialPlot(mod.rf.fn2$finalModel, carbon, UPFpx)
partialPlot(mod.rf.fn2$finalModel, carbon, MFls)
partialPlot(mod.rf.fn2$finalModel, carbon, SFls)
partialPlot(mod.rf.fn2$finalModel, carbon, distroad)


#rm(list= ls()[!(ls() %in% c("sel.var.df", "carbon", "env.explanatory.var", "mod.mlr.fn2", "mod.rf.fn2"))])
#gc()




#PGM 2010 real
pgm.2010real.raster.list <- list.files("rasters/PGM/2010_real/", pattern = ".tif", full.names = T, recursive = T)
pgm.2010real <- stack(pgm.2010real.raster.list)
pgm.2010real.rf <- pgm.2010real[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_pgm.2010real <- predict(pgm.2010real.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_pgm.2010real, paste0("models.output/carbon.benefits/PGM_2010_real_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#PGM 2020 real
pgm.2020real.raster.list <- list.files("rasters/PGM/2020_real/", pattern = ".tif", full.names = T, recursive = T)
pgm.2020real <- stack(pgm.2020real.raster.list)
pgm.2020real.rf <- pgm.2020real[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_pgm.2020real <- predict(pgm.2020real.rf, mod.rf.fn2, type="raw")
  
  
writeRaster(mod.rf.proj_pgm.2020real, paste0("models.output/carbon.benefits/PGM_2020_real_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#PGM 2020 avoiddeforest (all)
pgm.2020_avoiddeforest.raster.list <- list.files("rasters/PGM/2020_avoiddeforest/", pattern = ".tif", full.names = T, recursive = T)
pgm.2020_avoiddeforest <- stack(pgm.2020_avoiddeforest.raster.list)
pgm.2020_avoiddeforest.rf <- pgm.2020_avoiddeforest[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_pgm.2020_avoiddeforest <- predict(pgm.2020_avoiddeforest.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_pgm.2020_avoiddeforest, paste0("models.output/carbon.benefits/PGM_2020_avoiddeforest_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#PGM 2020 avoiddeforest (primary forest only)
pgm.2020_avoiddeforest2.raster.list <- list.files("rasters/PGM/2020_avoiddeforest2/", pattern = ".tif", full.names = T, recursive = T)
pgm.2020_avoiddeforest2 <- stack(pgm.2020_avoiddeforest2.raster.list)
pgm.2020_avoiddeforest2.rf <- pgm.2020_avoiddeforest2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_pgm.2020_avoiddeforest2 <- predict(pgm.2020_avoiddeforest2.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_pgm.2020_avoiddeforest2, paste0("models.output/carbon.benefits/PGM_2020_avoiddeforest2_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#PGM 2020 avoiddegrad (all)
pgm.2020_avoiddegrad.raster.list <- list.files("rasters/PGM/2020_avoiddegrad/", pattern = ".tif", full.names = T, recursive = T)
pgm.2020_avoiddegrad <- stack(pgm.2020_avoiddegrad.raster.list)
pgm.2020_avoiddegrad.rf <- pgm.2020_avoiddegrad[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_pgm.2020_avoiddegrad <- predict(pgm.2020_avoiddegrad.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_pgm.2020_avoiddegrad, paste0("models.output/carbon.benefits/PGM_2020_avoiddegrad_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#PGM 2020 avoiddegrad (primary forest only)
pgm.2020_avoiddegrad2.raster.list <- list.files("rasters/PGM/2020_avoiddegrad2/", pattern = ".tif", full.names = T, recursive = T)
pgm.2020_avoiddegrad2 <- stack(pgm.2020_avoiddegrad2.raster.list)
pgm.2020_avoiddegrad2.rf <- pgm.2020_avoiddegrad2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_pgm.2020_avoiddegrad2 <- predict(pgm.2020_avoiddegrad2.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_pgm.2020_avoiddegrad2, paste0("models.output/carbon.benefits/PGM_2020_avoiddegrad2_carbon_benefit.tif"), format = "GTiff", overwrite = T)
#writeRaster(proj_pgm.2020_avoiddegrad2.conbywm, paste0("models.output/carbon.benefits/PGM_2020_avoiddegrad2_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#PGM 2020 restor_wo_avoid
pgm.2020_restor_wo_avoid.raster.list <- list.files("rasters/PGM/2020_restor_wo_avoid/", pattern = ".tif", full.names = T, recursive = T)
pgm.2020_restor_wo_avoid <- stack(pgm.2020_restor_wo_avoid.raster.list)
pgm.2020_restor_wo_avoid.rf <- pgm.2020_restor_wo_avoid[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_pgm.2020_restor_wo_avoid <- predict(pgm.2020_restor_wo_avoid.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_pgm.2020_restor_wo_avoid, paste0("models.output/carbon.benefits/PGM_2020_restor_wo_avoid_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#PGM 2020 avoidboth (all)
pgm.2020_avoidboth.raster.list <- list.files("rasters/PGM/2020_avoidboth/", pattern = ".tif", full.names = T, recursive = T)
pgm.2020_avoidboth <- stack(pgm.2020_avoidboth.raster.list)
pgm.2020_avoidboth.rf <- pgm.2020_avoidboth[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_pgm.2020_avoidboth <- predict(pgm.2020_avoidboth.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_pgm.2020_avoidboth, paste0("models.output/carbon.benefits/PGM_2020_avoidboth_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#PGM 2020 avoidboth (primary forest only)
pgm.2020_avoidboth2.raster.list <- list.files("rasters/PGM/2020_avoidboth2/", pattern = ".tif", full.names = T, recursive = T)
pgm.2020_avoidboth2 <- stack(pgm.2020_avoidboth2.raster.list)
pgm.2020_avoidboth2.rf <- pgm.2020_avoidboth2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_pgm.2020_avoidboth2 <- predict(pgm.2020_avoidboth2.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_pgm.2020_avoidboth2, paste0("models.output/carbon.benefits/PGM_2020_avoidboth2_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#PGM 2020 restor_n_avoiddeforest (all)
pgm.2020_restor_n_avoiddeforest.raster.list <- list.files("rasters/PGM/2020_restor_n_avoiddeforest/", pattern = ".tif", full.names = T, recursive = T)
pgm.2020_restor_n_avoiddeforest <- stack(pgm.2020_restor_n_avoiddeforest.raster.list)
pgm.2020_restor_n_avoiddeforest.rf <- pgm.2020_restor_n_avoiddeforest[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_pgm.2020_restor_n_avoiddeforest <- predict(pgm.2020_restor_n_avoiddeforest.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_pgm.2020_restor_n_avoiddeforest, paste0("models.output/carbon.benefits/PGM_2020_restor_n_avoiddeforest_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#PGM 2020 restor_n_avoiddeforest (primary forest only)
pgm.2020_restor_n_avoiddeforest2.raster.list <- list.files("rasters/PGM/2020_restor_n_avoiddeforest2/", pattern = ".tif", full.names = T, recursive = T)
pgm.2020_restor_n_avoiddeforest2 <- stack(pgm.2020_restor_n_avoiddeforest2.raster.list)
pgm.2020_restor_n_avoiddeforest2.rf <- pgm.2020_restor_n_avoiddeforest2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_pgm.2020_restor_n_avoiddeforest2 <- predict(pgm.2020_restor_n_avoiddeforest2.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_pgm.2020_restor_n_avoiddeforest2, paste0("models.output/carbon.benefits/PGM_2020_restor_n_avoiddeforest2_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#PGM 2020 restor_n_avoidboth (all)
pgm.2020_restor_n_avoidboth.raster.list <- list.files("rasters/PGM/2020_restor_n_avoidboth/", pattern = ".tif", full.names = T, recursive = T)
pgm.2020_restor_n_avoidboth <- stack(pgm.2020_restor_n_avoidboth.raster.list)
pgm.2020_restor_n_avoidboth.rf <- pgm.2020_restor_n_avoidboth[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_pgm.2020_restor_n_avoidboth <- predict(pgm.2020_restor_n_avoidboth.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_pgm.2020_restor_n_avoidboth, paste0("models.output/carbon.benefits/PGM_2020_restor_n_avoidboth_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#PGM 2020 restor_n_avoidboth (primary forest only)
pgm.2020_restor_n_avoidboth2.raster.list <- list.files("rasters/PGM/2020_restor_n_avoidboth2/", pattern = ".tif", full.names = T, recursive = T)
pgm.2020_restor_n_avoidboth2 <- stack(pgm.2020_restor_n_avoidboth2.raster.list)
pgm.2020_restor_n_avoidboth2.rf <- pgm.2020_restor_n_avoidboth2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_pgm.2020_restor_n_avoidboth2 <- predict(pgm.2020_restor_n_avoidboth2.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_pgm.2020_restor_n_avoidboth2, paste0("models.output/carbon.benefits/PGM_2020_restor_n_avoidboth2_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#
#



#STM 2010 real
stm.2010real.raster.list <- list.files("rasters/STM/2010_real/", pattern = ".tif", full.names = T, recursive = T)
stm.2010real <- stack(stm.2010real.raster.list)
stm.2010real.rf <- stm.2010real[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_stm.2010real <- predict(stm.2010real.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_stm.2010real, paste0("models.output/carbon.benefits/STM_2010_real_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#STM 2020 real
stm.2020real.raster.list <- list.files("rasters/STM/2020_real/", pattern = ".tif", full.names = T, recursive = T)
stm.2020real <- stack(stm.2020real.raster.list)
stm.2020real.rf <- stm.2020real[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_stm.2020real <- predict(stm.2020real.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_stm.2020real, paste0("models.output/carbon.benefits/STM_2020_real_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#STM 2020 avoiddeforest (all)
stm.2020_avoiddeforest.raster.list <- list.files("rasters/STM/2020_avoiddeforest/", pattern = ".tif", full.names = T, recursive = T)
stm.2020_avoiddeforest <- stack(stm.2020_avoiddeforest.raster.list)
stm.2020_avoiddeforest.rf <- stm.2020_avoiddeforest[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_stm.2020_avoiddeforest <- predict(stm.2020_avoiddeforest.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_stm.2020_avoiddeforest, paste0("models.output/carbon.benefits/STM_2020_avoiddeforest_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#STM 2020 avoiddeforest (primary forest only)
stm.2020_avoiddeforest2.raster.list <- list.files("rasters/STM/2020_avoiddeforest2/", pattern = ".tif", full.names = T, recursive = T)
stm.2020_avoiddeforest2 <- stack(stm.2020_avoiddeforest2.raster.list)
stm.2020_avoiddeforest2.rf <- stm.2020_avoiddeforest2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_stm.2020_avoiddeforest2 <- predict(stm.2020_avoiddeforest2.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_stm.2020_avoiddeforest2, paste0("models.output/carbon.benefits/STM_2020_avoiddeforest2_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#STM 2020 avoiddegrad (all)
stm.2020_avoiddegrad.raster.list <- list.files("rasters/STM/2020_avoiddegrad/", pattern = ".tif", full.names = T, recursive = T)
stm.2020_avoiddegrad <- stack(stm.2020_avoiddegrad.raster.list)
stm.2020_avoiddegrad.rf <- stm.2020_avoiddegrad[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_stm.2020_avoiddegrad <- predict(stm.2020_avoiddegrad.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_stm.2020_avoiddegrad, paste0("models.output/carbon.benefits/STM_2020_avoiddegrad_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#STM 2020 avoiddegrad (primary forest only)
stm.2020_avoiddegrad2.raster.list <- list.files("rasters/STM/2020_avoiddegrad2/", pattern = ".tif", full.names = T, recursive = T)
stm.2020_avoiddegrad2 <- stack(stm.2020_avoiddegrad2.raster.list)
stm.2020_avoiddegrad2.rf <- stm.2020_avoiddegrad2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_stm.2020_avoiddegrad2 <- predict(stm.2020_avoiddegrad2.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_stm.2020_avoiddegrad2, paste0("models.output/carbon.benefits/STM_2020_avoiddegrad2_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#STM 2020 restor_wo_avoid
stm.2020_restor_wo_avoid.raster.list <- list.files("rasters/STM/2020_restor_wo_avoid/", pattern = ".tif", full.names = T, recursive = T)
stm.2020_restor_wo_avoid <- stack(stm.2020_restor_wo_avoid.raster.list)
stm.2020_restor_wo_avoid.rf <- stm.2020_restor_wo_avoid[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_stm.2020_restor_wo_avoid <- predict(stm.2020_restor_wo_avoid.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_stm.2020_restor_wo_avoid, paste0("models.output/carbon.benefits/STM_2020_restor_wo_avoid_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#STM 2020 avoidboth (all)
stm.2020_avoidboth.raster.list <- list.files("rasters/STM/2020_avoidboth/", pattern = ".tif", full.names = T, recursive = T)
stm.2020_avoidboth <- stack(stm.2020_avoidboth.raster.list)
stm.2020_avoidboth.rf <- stm.2020_avoidboth[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_stm.2020_avoidboth <- predict(stm.2020_avoidboth.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_stm.2020_avoidboth, paste0("models.output/carbon.benefits/STM_2020_avoidboth_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#STM 2020 avoidboth (primary forest only)
stm.2020_avoidboth2.raster.list <- list.files("rasters/STM/2020_avoidboth2/", pattern = ".tif", full.names = T, recursive = T)
stm.2020_avoidboth2 <- stack(stm.2020_avoidboth2.raster.list)
stm.2020_avoidboth2.rf <- stm.2020_avoidboth2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_stm.2020_avoidboth2 <- predict(stm.2020_avoidboth2.rf, mod.rf.fn2, type="raw")



writeRaster(mod.rf.proj_stm.2020_avoidboth2, paste0("models.output/carbon.benefits/STM_2020_avoidboth2_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#STM 2020 restor_n_avoiddeforest (all)
stm.2020_restor_n_avoiddeforest.raster.list <- list.files("rasters/STM/2020_restor_n_avoiddeforest/", pattern = ".tif", full.names = T, recursive = T)
stm.2020_restor_n_avoiddeforest <- stack(stm.2020_restor_n_avoiddeforest.raster.list)
stm.2020_restor_n_avoiddeforest.rf <- stm.2020_restor_n_avoiddeforest[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_stm.2020_restor_n_avoiddeforest <- predict(stm.2020_restor_n_avoiddeforest.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_stm.2020_restor_n_avoiddeforest, paste0("models.output/carbon.benefits/STM_2020_restor_n_avoiddeforest_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#STM 2020 restor_n_avoiddeforest (primary forest only)
stm.2020_restor_n_avoiddeforest2.raster.list <- list.files("rasters/STM/2020_restor_n_avoiddeforest2/", pattern = ".tif", full.names = T, recursive = T)
stm.2020_restor_n_avoiddeforest2 <- stack(stm.2020_restor_n_avoiddeforest2.raster.list)
stm.2020_restor_n_avoiddeforest2.rf <- stm.2020_restor_n_avoiddeforest2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_stm.2020_restor_n_avoiddeforest2 <- predict(stm.2020_restor_n_avoiddeforest2.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_stm.2020_restor_n_avoiddeforest2, paste0("models.output/carbon.benefits/STM_2020_restor_n_avoiddeforest2_carbon_benefit.tif"), format = "GTiff", overwrite = T)
#writeRaster(proj_stm.2020_restor_n_avoiddeforest2.conbywm, paste0("models.output/carbon.benefits/STM_2020_restor_n_avoiddeforest2_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#STM 2020 restor_n_avoidboth (all)
stm.2020_restor_n_avoidboth.raster.list <- list.files("rasters/STM/2020_restor_n_avoidboth/", pattern = ".tif", full.names = T, recursive = T)
stm.2020_restor_n_avoidboth <- stack(stm.2020_restor_n_avoidboth.raster.list)
stm.2020_restor_n_avoidboth.rf <- stm.2020_restor_n_avoidboth[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_stm.2020_restor_n_avoidboth <- predict(stm.2020_restor_n_avoidboth.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_stm.2020_restor_n_avoidboth, paste0("models.output/carbon.benefits/STM_2020_restor_n_avoidboth_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#



#STM 2020 restor_n_avoidboth (primary forest only)
stm.2020_restor_n_avoidboth2.raster.list <- list.files("rasters/STM/2020_restor_n_avoidboth2/", pattern = ".tif", full.names = T, recursive = T)
stm.2020_restor_n_avoidboth2 <- stack(stm.2020_restor_n_avoidboth2.raster.list)
stm.2020_restor_n_avoidboth2.rf <- stm.2020_restor_n_avoidboth2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]

set.seed(999)
mod.rf.proj_stm.2020_restor_n_avoidboth2 <- predict(stm.2020_restor_n_avoidboth2.rf, mod.rf.fn2, type="raw")


writeRaster(mod.rf.proj_stm.2020_restor_n_avoidboth2, paste0("models.output/carbon.benefits/STM_2020_restor_n_avoidboth2_carbon_benefit.tif"), format = "GTiff", overwrite = T)

#
#




# opportunity cost =============================================================

## property data
property <- read.csv("data/raw/areas.csv")
#head(property)
#str(property)
#summary(property)

#data on property area are from 2006 and 2009
#some properties have no information on their area in 2009 (n=25)
#in these cases we are using the information from 2006 (n=20)
property[which(is.na(property$totalarea_2009)),"totalarea_2009"] <- property[which(is.na(property$totalarea_2009)),"totalarea_2006"]
#we are excluding properties with no information on their total area (n=5)
property <- property[-which(is.na(property$totalarea_2009)),c("cd_propriedade", "totalarea_2009")]
names(property) <- c("id", "propertysize")

## income/profit data
revenue <- read.csv("data/raw/revenue.csv")
#head(revenue)
#str(revenue)
#summary(revenue)
names(revenue)[1] <- "id"
revenue$catchment <- as.factor(revenue$catchment)

#merging and calculating profit/ha
property <- property %>% left_join(revenue) %>% mutate(profit_ha = profit/propertysize)
#table(property$catchment)

#extracting the explanatory variables by catchment
expl.var <- carbon[,c(2,13:23,25:28)]
names(expl.var)[1] <- names(property)[4]

#merging
property <- property %>% left_join(expl.var, multiple = "all", relationship = "many-to-many") %>% group_by(id) %>% 
  sample_n(1, replace = T) %>% ungroup()



##scaling predictors -- z-scores
#property <- property %>% mutate(dominant = factor(dominant, levels = c("Cattlebeef","Cattlemilk","Animalother","Soy","Annual","Perennial","Mixed","Other")),
#                                distmarket_z = ((distmarket - mean(distmarket))/sd(distmarket)),
#                                property_z = ((property - mean(property))/sd(property)),
#                                distriver_z = ((distriver - mean(distriver))/sd(distriver)),
#                                distroad_z = ((distroad - mean(distroad))/sd(distroad)),
#                                DPFpx_z = ((DPFpx - mean(DPFpx))/sd(DPFpx)),
#                                edgedist_z = ((edgedist - mean(edgedist))/sd(edgedist)),
#                                edgels_z = ((edgels - mean(edgels))/sd(edgels)),
#                                elevation_z = ((elevation - mean(elevation))/sd(elevation)),
#                                meanprecips_z = ((meanprecips - mean(meanprecips))/sd(meanprecips)),
#                                meantemps_z = ((meantemps - mean(meantemps))/sd(meantemps)),
#                                SFagels_z = ((SFagels - mean(SFagels))/sd(SFagels)),
#                                SFpx_z = ((SFpx - mean(SFpx))/sd(SFpx)),
#                                TFpx_z = ((TFpx - mean(TFpx))/sd(TFpx)),
#                                UPFls_z = ((UPFls - mean(UPFls))/sd(UPFls)),
#                                TSDls_z = ((TSDls - mean(TSDls))/sd(TSDls)))

write.csv(property, "data/property.csv", row.names = F)
#property <- read.csv("data/property.csv")

rm(list= ls()[!(ls() %in% c("sel.var.df", "carbon", "property", "env.explanatory.var"))])
gc()


#keeping only properties with positive profit (n=53)
property.total <- property
property <- property %>% dplyr::filter(profit_ha > 0 & !is.nan(profit_ha) & !is.infinite(profit_ha))

#excluding outliers -- three properties with area lower than 6ha but profits greater than $10k/ha:
#property %>% dplyr::filter(propertysize < 10 & profit_ha > 10000)
#id == 419 (PGM), 461 and 322 (STM)
property <- property %>% dplyr::filter(!id %in% c(419, 461, 322))

property <- as.data.frame(property)

# distribution family
property.norm <- fitdist(property$profit_ha, "norm", method="mle")
summary(property.norm)

property.lnorm <- fitdist(property$profit_ha, "lnorm", method="mle")
summary(property.lnorm)

property.gamma <- fitdist(property$profit_ha, "gamma", method="mme")
summary(property.gamma)

property.weibull <- fitdist(property$profit_ha, "weibull", method="mle")
summary(property.weibull)


par(mfrow=c(2,2))
plot.legend <- c("lognormal", "gamma", "weibull")
denscomp(list(property.lnorm, property.gamma, property.weibull), legendtext = plot.legend)
cdfcomp(list(property.lnorm, property.gamma, property.weibull), legendtext = plot.legend)
qqcomp(list(property.lnorm, property.gamma, property.weibull), legendtext = plot.legend)
ppcomp(list(property.lnorm, property.gamma, property.weibull), legendtext = plot.legend)


# log of profit
property <- property %>% mutate(profit_halog = log(profit_ha))
property$dominant <- as.factor(property$dominant)


par(mfrow=c(1,1))
plot(Fn <- ecdf(property$profit_halog))
curve(pnorm(x, mean(property$profit_halog), sd(property$profit_halog)), from = 0, to = 10, add = TRUE, col='darkgreen', lwd = 2)
#
#curve(pgamma(x, shape=0.1680618477, rate=0.0002775644), from = 0, to = 17000, add = TRUE, col='blue', lwd = 2)
#
#curve(pweibull(x, shape=0.6446108, scale=405.1428443), from = 0, to = 17000, add = TRUE, col='red', lwd = 2)


#model fiting and cross-validated predictions for Random Forest
#setting the algorithm to use in caret by defining a list that contains a number of custom named elements that the caret package looks for
customRF <- list(type = "Regression", library = "randomForest", loop = NULL)

customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))

customRF$grid <- function(x, y, len = NULL, search = "grid") {}

customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}

customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)

customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")

customRF$sort <- function(x) x[order(x[,1]),]

customRF$levels <- function(x) x$classes


#train model
#custom resampling method in caret: 'leave-pair-out-cross-validation'
#catchments are used as training set, and tested against the 2 left out
num_samps <- length(unique(property$catchment))
holdout <- combn(num_samps, 2)
in_training <- apply(holdout, 2, 
                     function(x, all) all[!(all %in% x)],
                     all = 1:num_samps)

## need a more effcient way of doing this:
index <- vector(mode = "list", length = ncol(in_training))
for(i in 1:ncol(in_training)) index[[i]] <- in_training[,i]
## cosmetic:
names(index) <- caret:::prettySeq(seq(along = index))

ctrl <- trainControl(method = "cv", ## this will be ignored since we supply index
                     index = index,
                     savePredictions = T)  
#control <- trainControl(method="repeatedcv", number=10, repeats=3)
tunegrid <- expand.grid(.mtry=c(3:14), .ntree=seq(100,2000, 50))

set.seed(999)
mod.rf.fn3 <- train(profit_halog ~ propertysize + distmarket + distriver + distroad + DPFpx + 
                     edgedist + edgepx + edgels + elevation + meanprecips +
                     meantemps + MFls + SFAgepx + SFls + TSDls + UPFpx, data = property,
                   method=customRF, metric="RMSE", tuneGrid=tunegrid, trControl=ctrl)
print(mod.rf.fn3)
plot(mod.rf.fn3)


layout(matrix(c(1,2,3,4,1,5,6,7), 2, 4, byrow = TRUE))
varImpPlot(mod.rf.fn3$finalModel)
partialPlot(mod.rf.fn3$finalModel, property, propertysize)
partialPlot(mod.rf.fn3$finalModel, property, distmarket)
partialPlot(mod.rf.fn3$finalModel, property, meantemps)
partialPlot(mod.rf.fn3$finalModel, property, elevation)
partialPlot(mod.rf.fn3$finalModel, property, meanprecips)
partialPlot(mod.rf.fn3$finalModel, property, distroad)




#predctions
#PGM
pgm.2010real.raster.list <- list.files("rasters/PGM/2010_real/", pattern = ".tif", full.names = T, recursive = T)
pgm.2010real <- stack(pgm.2010real.raster.list)
pgm.2010real.rf <- pgm.2010real[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]


set.seed(999)
mod.rf.proj_pgm.2010real <- exp(predict(pgm.2010real.rf, mod.rf.fn3, type="raw"))


writeRaster(mod.rf.proj_pgm.2010real, "models.output/costs/PGM_2010_real_base_opportunity.tif", format = "GTiff", overwrite = T)

#



#STM
stm.2010real.raster.list <- list.files("rasters/STM/2010_real/", pattern = ".tif", full.names = T, recursive = T)
stm.2010real <- stack(stm.2010real.raster.list)
stm.2010real.rf <- stm.2010real[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]


set.seed(999)
mod.rf.proj_stm.2010real <- exp(predict(stm.2010real.rf, mod.rf.fn3, type="raw"))


writeRaster(mod.rf.proj_stm.2010real, "models.output/costs/STM_2010_real_base_opportunity.tif", format = "GTiff", overwrite = T)

#

#par(mfrow=c(1,2))
#plot(mod.rf.proj_pgm.2010real, col = terrain.colors(length(seq(0, 1500, by = 125)), rev = T), breaks= seq(0, 1500, by = 125))
#plot(mod.rf.proj_stm.2010real, col = terrain.colors(length(seq(0, 1500, by = 125)), rev = T), breaks= seq(0, 1500, by = 125))


#
#



# harvestable tree volume value and operational costs ==========================

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

# new var with code for transect by catchment
treedata$Transectcode <- paste(treedata$Catchment, treedata$Transect, sep="_")
# new var with spp binomial
treedata$Binomial <- paste(treedata$Genera, treedata$Species, sep="")

# excluding Varzea transects PGM
exclude <- c("100_1", "100_4","100_7","81_12","423_2") #transect code
treedata <- treedata[!treedata$Transectcode %in% exclude,]

# Pará commercial tree price data
# source: http://www.sefa.pa.gov.br/legislacao/interna/portaria/ps2015_00005.pdf
raw.timber.prices <- read.csv("data/raw/raw_timber_values_pa.csv")

# transect tree volume and price (trees with DBH>=40)
# source: Vidal et al. 2016 - http://dx.doi.org/10.1016/j.foreco.2016.06.003
commercial.treedata <- treedata %>% filter(DBH>=40) %>% 
                            left_join(raw.timber.prices) %>% 
                               mutate(Volume = -6.86 + (1.994*log(DBH)),
                                      Tree_value = Volume * Interno..R..)

# cumulative sum of the most valuable trees up to a maximum total harvestable volume of 30.1m3 over 35 years
# assuming forest productivity of 0.86 m3/ha/year
# source: BRASIL, RES. CONAMA No 406, DE 02 DE FEVEREIRO DE 2009
commercial.treedata <- commercial.treedata %>% 
                             group_by(Transectcode) %>% 
                                    arrange(desc(Tree_value)) %>% 
                                    mutate(y = cumsum(Volume)) %>% 
                                    dplyr::filter(y<=30.1 & Tree_value > 0)

# costs of planning and operational of the exploration
# the IMAZON estimated the general cost of planning and the operational costs in Paragominas-PA, 2013
# ~ US$26.49/m3 -- US$1 == R$2.5 -- R$66.23/m3
#source: https://imazon.org.br/custos-e-beneficios-do-manejo-florestal-para-a-producao-de-madeira-na-amazonia-oriental-n-10/
commercial.treedata <- commercial.treedata %>% mutate(costs = Volume * 66.23,
                                                      Tree_profit = Tree_value - costs)


# sum of the harvestable tree volume value by transect
transect.harvest.value <- commercial.treedata %>% group_by(Transectcode) %>% 
                                  summarise(profit = sum(Tree_profit, na.rm = T),
                                            Region = unique(Region),
                                            Catchment = unique(Catchment),
                                            Transect = unique(Transect)) %>% 
                                  mutate(profit_ha = profit/0.25,
                                         profit_ha_year = profit_ha/25) %>% 
                                  ungroup() %>% 
                                  dplyr::filter(profit > 0)

#head(transect.harvest.value)
#str(transect.harvest.value)
#summary(transect.harvest.value)

#extracting the explanatory variables by catchment
expl.var <- carbon[,c(2:4,13:28)]
expl.var$Transectcode <- gsub("ExtraPGM", "ExtraPGM_", expl.var$Transectcode)

#merging
transect.harvest.value <- transect.harvest.value %>% left_join(expl.var, multiple = "all", relationship = "many-to-many") %>% 
  group_by(Transectcode) %>% 
  sample_n(1, replace = T) %>% 
  ungroup()



##scaling predictors -- z-scores
#transect.harvest.value <- transect.harvest.value %>% mutate(distmarket_z = ((distmarket - mean(distmarket))/sd(distmarket)),
#                                                            distriver_z = ((distriver - mean(distriver))/sd(distriver)),
#                                                            distroad_z = ((distroad - mean(distroad))/sd(distroad)),
#                                                            DPFpx_z = ((DPFpx - mean(DPFpx))/sd(DPFpx)),
#                                                            edgedist_z = ((edgedist - mean(edgedist))/sd(edgedist)),
#                                                            edgels_z = ((edgels - mean(edgels))/sd(edgels)),
#                                                            elevation_z = ((elevation - mean(elevation))/sd(elevation)),
#                                                            meanprecips_z = ((meanprecips - mean(meanprecips))/sd(meanprecips)),
#                                                            meantemps_z = ((meantemps - mean(meantemps))/sd(meantemps)),
#                                                            SFagels_z = ((SFagels - mean(SFagels))/sd(SFagels)),
#                                                            SFpx_z = ((SFpx - mean(SFpx))/sd(SFpx)),
#                                                            TFpx_z = ((TFpx - mean(TFpx))/sd(TFpx)),
#                                                            UPFls_z = ((UPFls - mean(UPFls))/sd(UPFls)),
#                                                            TSDls_z = ((TSDls - mean(TSDls))/sd(TSDls)))

write.csv(transect.harvest.value, "data/transect_harvest_value.csv", row.names = F)
#transect.harvest.value <- read.csv("data/transect_harvest_value.csv")

rm(list= ls()[!(ls() %in% c("sel.var.df", "carbon", "property", "transect.harvest.value", "env.explanatory.var"))])
gc()


# distribution family
harvest.value.norm <- fitdist(transect.harvest.value$profit_ha_year, "norm", method="mle")
summary(harvest.value.norm)

harvest.value.lnorm <- fitdist(transect.harvest.value$profit_ha_year, "lnorm", method="mle")
summary(harvest.value.lnorm)

harvest.value.gamma <- fitdist(transect.harvest.value$profit_ha_year, "gamma", method="mme")
summary(harvest.value.gamma)

harvest.value.expo <- fitdist(transect.harvest.value$profit_ha_year, "exp", method="mme")
summary(harvest.value.expo)


par(mfrow=c(2,2))
plot.legend <- c("lognormal", "gamma", "exponential")
denscomp(list(harvest.value.lnorm, harvest.value.gamma, harvest.value.expo), legendtext = plot.legend)
cdfcomp(list(harvest.value.lnorm, harvest.value.gamma, harvest.value.expo), legendtext = plot.legend)
qqcomp(list(harvest.value.lnorm, harvest.value.gamma, harvest.value.expo), legendtext = plot.legend)
ppcomp(list(harvest.value.lnorm, harvest.value.gamma, harvest.value.expo), legendtext = plot.legend)


#plot(Fn <- ecdf(transect.harvest.value$profit_ha_year))
#
#curve(plnorm(x, mean(log(transect.harvest.value$profit_ha_year)), sd(log(transect.harvest.value$profit_ha_year))), from = 317, to = 13500, add = TRUE, col='darkgreen', lwd = 2)
#
#curve(pgamma(x, shape=1.5917114890, rate=0.0005607584), from = 317, to = 13500, add = TRUE, col='blue', lwd = 2)
#
#curve(pexp(x, rate=0.000352299), from = 317, to = 13500, add = TRUE, col='red', lwd = 2)

transect.harvest.value <- as.data.frame(transect.harvest.value)

#model fiting and cross-validated predictions for Random Forest
#setting the algorithm to use in caret by defining a list that contains a number of custom named elements that the caret package looks for
customRF <- list(type = "Regression", library = "randomForest", loop = NULL)

customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))

customRF$grid <- function(x, y, len = NULL, search = "grid") {}

customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}

customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)

customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")

customRF$sort <- function(x) x[order(x[,1]),]

customRF$levels <- function(x) x$classes


#train model
#custom resampling method in caret: 'leave-pair-out-cross-validation'
#catchments are used as training set, and tested against the 2 left out
num_samps <- length(unique(transect.harvest.value$Catchment))
holdout <- combn(num_samps, 2)
in_training <- apply(holdout, 2, 
                     function(x, all) all[!(all %in% x)],
                     all = 1:num_samps)

## need a more effcient way of doing this:
index <- vector(mode = "list", length = ncol(in_training))
for(i in 1:ncol(in_training)) index[[i]] <- in_training[,i]
## cosmetic:
names(index) <- caret:::prettySeq(seq(along = index))

ctrl <- trainControl(method = "cv", ## this will be ignored since we supply index
                     index = index,
                     savePredictions = T)  
#control <- trainControl(method="repeatedcv", number=10, repeats=3)
tunegrid <- expand.grid(.mtry=c(3:14), .ntree=seq(100,2000, 50))

set.seed(999)
mod.rf.fn4 <- train(profit_ha_year ~ propertysize + distmarket + distriver + distroad + DPFpx + 
                      edgedist + edgepx + edgels + elevation + meanprecips +
                      meantemps + MFls + SFAgepx + SFls + TSDls + UPFpx, data = transect.harvest.value,
                    method=customRF, metric="RMSE", tuneGrid=tunegrid, trControl=ctrl)
print(mod.rf.fn4)
plot(mod.rf.fn4)


layout(matrix(c(1,2,3,4,1,5,6,7), 2, 4, byrow = TRUE))
#plot(mod.rf.fn4)
varImpPlot(mod.rf.fn4$finalModel)
partialPlot(mod.rf.fn4$finalModel, transect.harvest.value, meanprecips)
partialPlot(mod.rf.fn4$finalModel, transect.harvest.value, meantemps)
partialPlot(mod.rf.fn4$finalModel, transect.harvest.value, edgedist)
partialPlot(mod.rf.fn4$finalModel, transect.harvest.value, distmarket)
partialPlot(mod.rf.fn4$finalModel, transect.harvest.value, distroad)
partialPlot(mod.rf.fn4$finalModel, transect.harvest.value, TSDls)




#predctions
#PGM
pgm.2010real.raster.list <- list.files("rasters/PGM/2010_real/", pattern = ".tif", full.names = T, recursive = T)
pgm.2010real <- stack(pgm.2010real.raster.list)
pgm.2010real.rf <- pgm.2010real[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]


set.seed(999)
mod.rf.proj_pgm.2010real <- predict(pgm.2010real.rf, mod.rf.fn4, type="raw")


writeRaster(mod.rf.proj_pgm.2010real, "models.output/costs/PGM_2010_real_base_haverst.tif", format = "GTiff", overwrite = T)

#



#STM
stm.2010real.raster.list <- list.files("rasters/STM/2010_real/", pattern = ".tif", full.names = T, recursive = T)
stm.2010real <- stack(stm.2010real.raster.list)
stm.2010real.rf <- stm.2010real[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]


set.seed(999)
mod.rf.proj_stm.2010real <- predict(stm.2010real.rf, mod.rf.fn4, type="raw")


writeRaster(mod.rf.proj_stm.2010real, "models.output/costs/STM_2010_real_base_haverst.tif", format = "GTiff", overwrite = T)

#

#par(mfrow=c(1,2))
#plot(mod.rf.proj_pgm.2010real, col = terrain.colors(length(seq(0, 200, by = 20)), rev = T), breaks= seq(0, 200, by = 20))
#plot(mod.rf.proj_stm.2010real, col = terrain.colors(length(seq(0, 200, by = 20)), rev = T), breaks= seq(0, 200, by = 20))


rm(list= ls())
gc()







# end =================================










































































#==============================| previous modeling approach


#' species distribution modelling


#i <- as.character(forestdep.spplist$Binomial[1])
for (i in forestdep.spplist$Binomial) {
  
  #spliting the job
  if(forestdep.spplist[forestdep.spplist$Binomial==i, "job"] != 1) next
  
  #checking if model is done
  if(forestdep.spplist[forestdep.spplist$Binomial==i, "Done"] == TRUE) next
  
  #starting modeling procedure
  occur <- sppdata.final[sppdata.final$Binomial==i,]
  
  ##data preparation
  myRespName <- as.character(i)
  myResp <- c(rep.int(1, times = as.numeric(nrow(occur))), rep.int(NA, times = as.numeric(nrow(occur)*30)))
  myRespXY <- rbind(occur[,c("Longitude", "Latitude")], sample_n(pa.data, as.numeric(nrow(occur)*30), replace = T))
  myExpl <- env.explanatory.var[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
  myPAtable <- data.frame(PA1 = ifelse(myResp == 1, TRUE, FALSE),
                          PA2 = ifelse(myResp == 1, TRUE, FALSE),
                          PA3 = ifelse(myResp == 1, TRUE, FALSE))
  
  for (m in 1:ncol(myPAtable)) myPAtable[sample(which(is.na(myPAtable[, m])), as.numeric(nrow(occur)*10)), m] = TRUE
  
  myPAtable <- as.data.frame(!is.na(myPAtable))
  
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = myExpl,
                                       resp.xy = myRespXY,
                                       resp.name = myRespName,
                                       PA.strategy = "user.defined",
                                       PA.user.table = myPAtable,
                                       dir.name = "models.output")
  
  ## checking
  #myBiomodData
  #nrow(myBiomodData@coord)
  #head(myBiomodData@coord, xx)
  #plot(env.explanatory.var[[1]])
  #points(myBiomodData@coord[c(1:xx),], col="red")
  #points(myBiomodData@coord[c(xx:nrow(myBiomodData@coord)),])
  
  # algorithms control
  myBiomodOption <- BIOMOD_ModelingOptions(RF = list(ntree=1000, nodesize=10),
                                           MAXENT.Phillips = list(memory_allocated = 4096)) #64000
  #
  # model fit and assessment parameters
  myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                      modeling.id = as.character(i),
                                      models = c('RF', 'MAXENT.Phillips'), #'GLM', 
                                      bm.options = myBiomodOption,
                                      nb.rep = 10,
                                      data.split.perc = 80,
                                      do.full.models = F,
                                      metric.eval = c('TSS', 'ROC'),
                                      var.import = 5,
                                      seed.val = 9999)
  
  # saving
  capture.output(myBiomodData, file = paste0("models.output/evaluation/", i,"_BiomodData.txt"))
  capture.output(myBiomodModelOut, file = paste0("models.output/evaluation/", i,"_BiomodModelOut.txt"))
  capture.output(get_evaluations(myBiomodModelOut), file = paste0("models.output/evaluation/", i,"_eval_BiomodModelOut.txt"))
  capture.output(get_variables_importance(myBiomodModelOut), file = paste0("models.output/evaluation/", i,"_var_importance.txt"))
  
  # adding var mean variable importance to selected exploratory variables data frame
  #varimport <- get_variables_importance(myBiomodModelOut)
  #sel.var.df[,i] <- NA
  #sel.var.df[!is.na(sel.var.df$VIF), i] <- rowMeans(varimport, na.rm=T)
  #write.csv(sel.var.df, "rasters/selected_environmental_explanatory_variables_byVIF.csv", row.names = F)
  
  # graphic of variables importance
  pdf(paste("models.output/evaluation/", i ,"_variables_importance.pdf", sep="" ))
  gVarImp <- bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c("algo", "run", "expl.var"))
  dev.off()
  
  # graphic of model scores
  pdf(paste("models.output/evaluation/", i ,"_models_scores.pdf", sep="" ))
  gEval <- bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'run'))
  dev.off()
  
  #' if the model is ready and need only project the maps
  #' start from here
  
  #sel.var.df <- read.csv(paste0("selected_exploratory_var/", i, "_VIF.csv"), header=T)
  #myBiomodModelOut <- load(paste0(i,"/",i,".",i,".models.out"))
  #myBiomodModelOut <- get(myBiomodModelOut)
  
  
  #updating evaluation threshold
  eval_threshold <- 0.7
  
  ev<-get_evaluations(myBiomodModelOut, as.data.frame=T)
  
  #selected_GLM <- subset (ev$Testing.data, grepl("GLM", ev$Model.name) & grepl("TSS", ev$Eval.metric)) 
  selected_RF <- subset (ev$Testing.data, grepl("RF", ev$Model.name) & grepl("TSS", ev$Eval.metric))
  selected_MAX <- subset (ev$Testing.data, grepl("MAXENT.Phillips", ev$Model.name) & grepl("TSS", ev$Eval.metric)) 
  
  new.eval_threshold <- min(c(max(selected_RF[!is.na(selected_RF)]),
                              max(selected_MAX[!is.na(selected_MAX)])))
  
  if (new.eval_threshold < 0.7) {eval_threshold = floor(new.eval_threshold*10)/10}
  
  
  all.models <- get_built_models(myBiomodModelOut)
  chosen.models.from.ev <- ev[ev$Eval.metric=="TSS" & ev$Testing.data>eval_threshold,]
  chosen.models.from.ev$Model.name2 <- paste(chosen.models.from.ev$Dataset, chosen.models.from.ev$Run, chosen.models.from.ev$Algo, sep="_")
  chosen.models.from.all.models <- grep(paste(chosen.models.from.ev$Model.name2,collapse="|"), all.models, value=TRUE)
  
  if(nrow(chosen.models.from.ev)<=5) next
  #nrow(chosen.models.from.ev)<5
  
  #
  #ensemble model scores to species list
  #EMeval <- as.data.frame(get_evaluations(myBiomodEM_all))
  EMeval <- ev[ev$Model.name %in% chosen.models.from.ev$Model.name,]
  forestdep.spplist[forestdep.spplist$Binomial==i,"ROC"] <- paste(round(mean(EMeval[EMeval$Eval.metric=="ROC","Testing.data"]),1), "+/-", round(sd(EMeval[EMeval$Eval.metric=="ROC","Testing.data"]),2))
  forestdep.spplist[forestdep.spplist$Binomial==i,"TSS"] <- paste(round(mean(EMeval[EMeval$Eval.metric=="TSS","Testing.data"]),1), "+/-", round(sd(EMeval[EMeval$Eval.metric=="TSS","Testing.data"]),2))
  forestdep.spplist[forestdep.spplist$Binomial==i,"Sensitivity"] <- paste(round(mean(EMeval[EMeval$Eval.metric=="TSS","Sensitivity"]),2), "+/-", round(sd(EMeval[EMeval$Eval.metric=="TSS","Sensitivity"]),2))
  forestdep.spplist[forestdep.spplist$Binomial==i,"Specificity"] <- paste(round(mean(EMeval[EMeval$Eval.metric=="TSS","Specificity"]),2), "+/-", round(sd(EMeval[EMeval$Eval.metric=="TSS","Specificity"]),2))
  forestdep.spplist[forestdep.spplist$Binomial==i,"Cutoff"] <- paste(round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2), "+/-", round((sd(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2))
  
  write.csv(forestdep.spplist, "models.output/evaluation/updated_species_summary_v1.csv")
  
  
  ##predctions
  if (any(occur$Region=="PGM")) {
    
    # 2020 real
    pgm.2020real.raster.list <- list.files("rasters/PGM/2020_real/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020real <- stack(pgm.2020real.raster.list)
    names(pgm.2020real) <- unlist(strsplit(pgm.2020real.raster.list, "/|.tif"))[seq(4,88,4)]
    rm("pgm.2020real.raster.list")
    
    pgm.2020real <- pgm.2020real[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    ## checking
    #names(pgm.2020real)
    #plot(pgm.2020real[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020real <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for real scenario in pgm")
      
      BIOMOD_LoadModels(bm.out = myBiomodModelOut, full.name = m, as = "mod")
      
      temp_workdir = NULL
      if (length(grep("MAXENT.Phillips$", m)) == 1) {
        temp_workdir = mod@model_output_dir
      }
      
      indivdual_proj_2020real <- predict(mod, pgm.2020real, temp_workdir = temp_workdir)
      names(indivdual_proj_2020real) <- m
      proj_2020real <- addLayer(proj_2020real, indivdual_proj_2020real)
      
    }
    
    
    
    # building a consensus map by mean weight
    proj_2020real.conbywm <- weighted.mean(proj_2020real, EMeval[EMeval$Eval.metric=="TSS","Testing.data"], na.rm=T)
    
    writeRaster(proj_2020real.conbywm, paste0("models.output/maps/PGM/2020_real/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
    
    #cutoff.th <- round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2)
    #
    #proj_2020real.conbywm.bin <- proj_2020real.conbywm
    #proj_2020real.conbywm.bin[proj_2020real.conbywm.bin>=cutoff.th]<-1
    #proj_2020real.conbywm.bin[proj_2020real.conbywm.bin<cutoff.th]<-0
    #
    #writeRaster(proj_2020real.conbywm.bin, paste0("models.output/maps/PGM/2020_real/", i, "_merged_algo_merged_dataset_merged_run_TSSBin.tif"), format = "GTiff", overwrite = T)
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "pa.data", "occur", "sel.var.df", "env.explanatory.var",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 avoid degradation
    pgm.2020avoiddegrad.raster.list <- list.files("rasters/PGM/2020_avoiddegrad/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020avoiddegrad <- stack(pgm.2020avoiddegrad.raster.list)
    names(pgm.2020avoiddegrad) <- unlist(strsplit(pgm.2020avoiddegrad.raster.list, "/|.tif"))[seq(4,88,4)]
    rm("pgm.2020avoiddegrad.raster.list")
    
    pgm.2020avoiddegrad <- pgm.2020avoiddegrad[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    ## checking
    #names(pgm.2020avoiddegrad)
    #plot(pgm.2020avoiddegrad[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020avoiddegrad <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for avoid degradation scenario in pgm")
      
      BIOMOD_LoadModels(bm.out = myBiomodModelOut, full.name = m, as = "mod")
      
      temp_workdir = NULL
      if (length(grep("MAXENT.Phillips$", m)) == 1) {
        temp_workdir = mod@model_output_dir
      }
      
      indivdual_proj_2020avoiddegrad <- predict(mod, pgm.2020avoiddegrad, temp_workdir = temp_workdir)
      names(indivdual_proj_2020avoiddegrad) <- m
      proj_2020avoiddegrad <- addLayer(proj_2020avoiddegrad, indivdual_proj_2020avoiddegrad)
      
    }
    
    
    
    # building a consensus map by mean weight
    proj_2020avoiddegrad.conbywm <- weighted.mean(proj_2020avoiddegrad, EMeval[EMeval$Eval.metric=="TSS","Testing.data"], na.rm=T)
    
    writeRaster(proj_2020avoiddegrad.conbywm, paste0("models.output/maps/PGM/2020_avoiddegrad/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
    
    #cutoff.th <- round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2)
    #
    #proj_2020avoiddegrad.conbywm.bin <- proj_2020avoiddegrad.conbywm
    #proj_2020avoiddegrad.conbywm.bin[proj_2020avoiddegrad.conbywm.bin>=cutoff.th]<-1
    #proj_2020avoiddegrad.conbywm.bin[proj_2020avoiddegrad.conbywm.bin<cutoff.th]<-0
    #
    #writeRaster(proj_2020avoiddegrad.conbywm.bin, paste0("models.output/maps/PGM/2020_avoiddegrad/", i, "_merged_algo_merged_dataset_merged_run_TSSBin.tif"), format = "GTiff", overwrite = T)
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "pa.data", "occur", "sel.var.df", "env.explanatory.var",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 avoid deforestation
    pgm.2020avoiddeforest.raster.list <- list.files("rasters/PGM/2020_avoiddeforest/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020avoiddeforest <- stack(pgm.2020avoiddeforest.raster.list)
    names(pgm.2020avoiddeforest) <- unlist(strsplit(pgm.2020avoiddeforest.raster.list, "/|.tif"))[seq(4,88,4)]
    rm("pgm.2020avoiddeforest.raster.list")
    
    pgm.2020avoiddeforest <- pgm.2020avoiddeforest[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    ## checking
    #names(pgm.2020avoiddeforest)
    #plot(pgm.2020avoiddeforest[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020avoiddeforest <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for avoid deforestation scenario in pgm")
      
      BIOMOD_LoadModels(bm.out = myBiomodModelOut, full.name = m, as = "mod")
      
      temp_workdir = NULL
      if (length(grep("MAXENT.Phillips$", m)) == 1) {
        temp_workdir = mod@model_output_dir
      }
      
      indivdual_proj_2020avoiddeforest <- predict(mod, pgm.2020avoiddeforest, temp_workdir = temp_workdir)
      names(indivdual_proj_2020avoiddeforest) <- m
      proj_2020avoiddeforest <- addLayer(proj_2020avoiddeforest, indivdual_proj_2020avoiddeforest)
      
    }
    
    
    
    # building a consensus map by mean weight
    proj_2020avoiddeforest.conbywm <- weighted.mean(proj_2020avoiddeforest, EMeval[EMeval$Eval.metric=="TSS","Testing.data"], na.rm=T)
    
    writeRaster(proj_2020avoiddeforest.conbywm, paste0("models.output/maps/PGM/2020_avoiddeforest/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
    
    #cutoff.th <- round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2)
    #
    #proj_2020avoiddeforest.conbywm.bin <- proj_2020avoiddeforest.conbywm
    #proj_2020avoiddeforest.conbywm.bin[proj_2020avoiddeforest.conbywm.bin>=cutoff.th]<-1
    #proj_2020avoiddeforest.conbywm.bin[proj_2020avoiddeforest.conbywm.bin<cutoff.th]<-0
    #
    #writeRaster(proj_2020avoiddeforest.conbywm.bin, paste0("models.output/maps/PGM/2020_avoiddeforest/", i, "_merged_algo_merged_dataset_merged_run_TSSBin.tif"), format = "GTiff", overwrite = T)
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "pa.data", "occur", "sel.var.df", "env.explanatory.var",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 avoid deforestation upf only
    pgm.2020avoiddeforest2.raster.list <- list.files("rasters/PGM/2020_avoiddeforest2/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020avoiddeforest2 <- stack(pgm.2020avoiddeforest2.raster.list)
    names(pgm.2020avoiddeforest2) <- unlist(strsplit(pgm.2020avoiddeforest2.raster.list, "/|.tif"))[seq(4,88,4)]
    rm("pgm.2020avoiddeforest2.raster.list")
    
    pgm.2020avoiddeforest2 <- pgm.2020avoiddeforest2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    ## checking
    #names(pgm.2020avoiddeforest2)
    #plot(pgm.2020avoiddeforest2[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020avoiddeforest2 <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for avoid deforestation upf only scenario in pgm")
      
      BIOMOD_LoadModels(bm.out = myBiomodModelOut, full.name = m, as = "mod")
      
      temp_workdir = NULL
      if (length(grep("MAXENT.Phillips$", m)) == 1) {
        temp_workdir = mod@model_output_dir
      }
      
      indivdual_proj_2020avoiddeforest2 <- predict(mod, pgm.2020avoiddeforest2, temp_workdir = temp_workdir)
      names(indivdual_proj_2020avoiddeforest2) <- m
      proj_2020avoiddeforest2 <- addLayer(proj_2020avoiddeforest2, indivdual_proj_2020avoiddeforest2)
      
    }
    
    
    
    # building a consensus map by mean weight
    proj_2020avoiddeforest2.conbywm <- weighted.mean(proj_2020avoiddeforest2, EMeval[EMeval$Eval.metric=="TSS","Testing.data"], na.rm=T)
    
    writeRaster(proj_2020avoiddeforest2.conbywm, paste0("models.output/maps/PGM/2020_avoiddeforest2/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
    
    #cutoff.th <- round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2)
    #
    #proj_2020avoiddeforest2.conbywm.bin <- proj_2020avoiddeforest2.conbywm
    #proj_2020avoiddeforest2.conbywm.bin[proj_2020avoiddeforest2.conbywm.bin>=cutoff.th]<-1
    #proj_2020avoiddeforest2.conbywm.bin[proj_2020avoiddeforest2.conbywm.bin<cutoff.th]<-0
    #
    #writeRaster(proj_2020avoiddeforest2.conbywm.bin, paste0("models.output/maps/PGM/2020_avoiddeforest2/", i, "_merged_algo_merged_dataset_merged_run_TSSBin.tif"), format = "GTiff", overwrite = T)
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "pa.data", "occur", "sel.var.df", "env.explanatory.var",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 avoid both
    pgm.2020avoidboth.raster.list <- list.files("rasters/PGM/2020_avoidboth/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020avoidboth <- stack(pgm.2020avoidboth.raster.list)
    names(pgm.2020avoidboth) <- unlist(strsplit(pgm.2020avoidboth.raster.list, "/|.tif"))[seq(4,88,4)]
    rm("pgm.2020avoidboth.raster.list")
    
    pgm.2020avoidboth <- pgm.2020avoidboth[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    ## checking
    #names(pgm.2020avoidboth)
    #plot(pgm.2020avoidboth[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020avoidboth <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for avoid both scenario in pgm")
      
      BIOMOD_LoadModels(bm.out = myBiomodModelOut, full.name = m, as = "mod")
      
      temp_workdir = NULL
      if (length(grep("MAXENT.Phillips$", m)) == 1) {
        temp_workdir = mod@model_output_dir
      }
      
      indivdual_proj_2020avoidboth <- predict(mod, pgm.2020avoidboth, temp_workdir = temp_workdir)
      names(indivdual_proj_2020avoidboth) <- m
      proj_2020avoidboth <- addLayer(proj_2020avoidboth, indivdual_proj_2020avoidboth)
      
    }
    
    
    
    # building a consensus map by mean weight
    proj_2020avoidboth.conbywm <- weighted.mean(proj_2020avoidboth, EMeval[EMeval$Eval.metric=="TSS","Testing.data"], na.rm=T)
    
    writeRaster(proj_2020avoidboth.conbywm, paste0("models.output/maps/PGM/2020_avoidboth/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
    
    #cutoff.th <- round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2)
    #
    #proj_2020avoidboth.conbywm.bin <- proj_2020avoidboth.conbywm
    #proj_2020avoidboth.conbywm.bin[proj_2020avoidboth.conbywm.bin>=cutoff.th]<-1
    #proj_2020avoidboth.conbywm.bin[proj_2020avoidboth.conbywm.bin<cutoff.th]<-0
    #
    #writeRaster(proj_2020avoidboth.conbywm.bin, paste0("models.output/maps/PGM/2020_avoidboth/", i, "_merged_algo_merged_dataset_merged_run_TSSBin.tif"), format = "GTiff", overwrite = T)
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "pa.data", "occur", "sel.var.df", "env.explanatory.var",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 restoration without avoid
    pgm.2020restor_wo_avoid.raster.list <- list.files("rasters/PGM/2020_restor_wo_avoid/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020restor_wo_avoid <- stack(pgm.2020restor_wo_avoid.raster.list)
    names(pgm.2020restor_wo_avoid) <- unlist(strsplit(pgm.2020restor_wo_avoid.raster.list, "/|.tif"))[seq(4,88,4)]
    rm("pgm.2020restor_wo_avoid.raster.list")
    
    pgm.2020restor_wo_avoid <- pgm.2020restor_wo_avoid[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    ## checking
    #names(pgm.2020restor_wo_avoid)
    #plot(pgm.2020restor_wo_avoid[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020restor_wo_avoid <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for restoration without avoid scenario in pgm")
      
      BIOMOD_LoadModels(bm.out = myBiomodModelOut, full.name = m, as = "mod")
      
      temp_workdir = NULL
      if (length(grep("MAXENT.Phillips$", m)) == 1) {
        temp_workdir = mod@model_output_dir
      }
      
      indivdual_proj_2020restor_wo_avoid <- predict(mod, pgm.2020restor_wo_avoid, temp_workdir = temp_workdir)
      names(indivdual_proj_2020restor_wo_avoid) <- m
      proj_2020restor_wo_avoid <- addLayer(proj_2020restor_wo_avoid, indivdual_proj_2020restor_wo_avoid)
      
    }
    
    
    
    # building a consensus map by mean weight
    proj_2020restor_wo_avoid.conbywm <- weighted.mean(proj_2020restor_wo_avoid, EMeval[EMeval$Eval.metric=="TSS","Testing.data"], na.rm=T)
    
    writeRaster(proj_2020restor_wo_avoid.conbywm, paste0("models.output/maps/PGM/2020_restor_wo_avoid/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
    
    #cutoff.th <- round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2)
    #
    #proj_2020restor_wo_avoid.conbywm.bin <- proj_2020restor_wo_avoid.conbywm
    #proj_2020restor_wo_avoid.conbywm.bin[proj_2020restor_wo_avoid.conbywm.bin>=cutoff.th]<-1
    #proj_2020restor_wo_avoid.conbywm.bin[proj_2020restor_wo_avoid.conbywm.bin<cutoff.th]<-0
    #
    #writeRaster(proj_2020restor_wo_avoid.conbywm.bin, paste0("models.output/maps/PGM/2020_restor_wo_avoid/", i, "_merged_algo_merged_dataset_merged_run_TSSBin.tif"), format = "GTiff", overwrite = T)
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "pa.data", "occur", "sel.var.df", "env.explanatory.var",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 restoration and avoid deforestation
    pgm.2020restor_n_avoid_deforest.raster.list <- list.files("rasters/PGM/2020_restor_n_avoid_deforest/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020restor_n_avoid_deforest <- stack(pgm.2020restor_n_avoid_deforest.raster.list)
    names(pgm.2020restor_n_avoid_deforest) <- unlist(strsplit(pgm.2020restor_n_avoid_deforest.raster.list, "/|.tif"))[seq(4,88,4)]
    rm("pgm.2020restor_n_avoid_deforest.raster.list")
    
    pgm.2020restor_n_avoid_deforest <- pgm.2020restor_n_avoid_deforest[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    ## checking
    #names(pgm.2020restor_n_avoid_deforest)
    #plot(pgm.2020restor_n_avoid_deforest[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020restor_n_avoid_deforest <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "restoration and avoid scenario in pgm")
      
      BIOMOD_LoadModels(bm.out = myBiomodModelOut, full.name = m, as = "mod")
      
      temp_workdir = NULL
      if (length(grep("MAXENT.Phillips$", m)) == 1) {
        temp_workdir = mod@model_output_dir
      }
      
      indivdual_proj_2020restor_n_avoid_deforest <- predict(mod, pgm.2020restor_n_avoid_deforest, temp_workdir = temp_workdir)
      names(indivdual_proj_2020restor_n_avoid_deforest) <- m
      proj_2020restor_n_avoid_deforest <- addLayer(proj_2020restor_n_avoid_deforest, indivdual_proj_2020restor_n_avoid_deforest)
      
    }
    
    
    
    # building a consensus map by mean weight
    proj_2020restor_n_avoid_deforest.conbywm <- weighted.mean(proj_2020restor_n_avoid_deforest, EMeval[EMeval$Eval.metric=="TSS","Testing.data"], na.rm=T)
    
    writeRaster(proj_2020restor_n_avoid_deforest.conbywm, paste0("models.output/maps/PGM/2020_restor_n_avoid_deforest/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
    
    #cutoff.th <- round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2)
    #
    #proj_2020restor_n_avoid_deforest.conbywm.bin <- proj_2020restor_n_avoid_deforest.conbywm
    #proj_2020restor_n_avoid_deforest.conbywm.bin[proj_2020restor_n_avoid_deforest.conbywm.bin>=cutoff.th]<-1
    #proj_2020restor_n_avoid_deforest.conbywm.bin[proj_2020restor_n_avoid_deforest.conbywm.bin<cutoff.th]<-0
    #
    #writeRaster(proj_2020restor_n_avoid_deforest.conbywm.bin, paste0("models.output/maps/PGM/2020_restor_n_avoid_deforest/", i, "_merged_algo_merged_dataset_merged_run_TSSBin.tif"), format = "GTiff", overwrite = T)
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "pa.data", "occur", "sel.var.df", "env.explanatory.var",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 restoration and avoid both
    pgm.2020restor_n_avoid_both.raster.list <- list.files("rasters/PGM/2020_restor_n_avoid_both/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020restor_n_avoid_both <- stack(pgm.2020restor_n_avoid_both.raster.list)
    names(pgm.2020restor_n_avoid_both) <- unlist(strsplit(pgm.2020restor_n_avoid_both.raster.list, "/|.tif"))[seq(4,88,4)]
    rm("pgm.2020restor_n_avoid_both.raster.list")
    
    pgm.2020restor_n_avoid_both <- pgm.2020restor_n_avoid_both[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    ## checking
    #names(pgm.2020restor_n_avoid_both)
    #plot(pgm.2020restor_n_avoid_both[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020restor_n_avoid_both <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "restoration and avoid scenario in pgm")
      
      BIOMOD_LoadModels(bm.out = myBiomodModelOut, full.name = m, as = "mod")
      
      temp_workdir = NULL
      if (length(grep("MAXENT.Phillips$", m)) == 1) {
        temp_workdir = mod@model_output_dir
      }
      
      indivdual_proj_2020restor_n_avoid_both <- predict(mod, pgm.2020restor_n_avoid_both, temp_workdir = temp_workdir)
      names(indivdual_proj_2020restor_n_avoid_both) <- m
      proj_2020restor_n_avoid_both <- addLayer(proj_2020restor_n_avoid_both, indivdual_proj_2020restor_n_avoid_both)
      
    }
    
    
    
    # building a consensus map by mean weight
    proj_2020restor_n_avoid_both.conbywm <- weighted.mean(proj_2020restor_n_avoid_both, EMeval[EMeval$Eval.metric=="TSS","Testing.data"], na.rm=T)
    
    writeRaster(proj_2020restor_n_avoid_both.conbywm, paste0("models.output/maps/PGM/2020_restor_n_avoid_both/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
    
    #cutoff.th <- round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2)
    #
    #proj_2020restor_n_avoid_both.conbywm.bin <- proj_2020restor_n_avoid_both.conbywm
    #proj_2020restor_n_avoid_both.conbywm.bin[proj_2020restor_n_avoid_both.conbywm.bin>=cutoff.th]<-1
    #proj_2020restor_n_avoid_both.conbywm.bin[proj_2020restor_n_avoid_both.conbywm.bin<cutoff.th]<-0
    #
    #writeRaster(proj_2020restor_n_avoid_both.conbywm.bin, paste0("models.output/maps/PGM/2020_restor_n_avoid_both/", i, "_merged_algo_merged_dataset_merged_run_TSSBin.tif"), format = "GTiff", overwrite = T)
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "pa.data", "occur", "sel.var.df", "env.explanatory.var",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
  }
  
  
  #  
  
  if (any(occur$Region=="STM")) {
    
    # 2020 real
    stm.2020real.raster.list <- list.files("rasters/STM/2020_real/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020real <- stack(stm.2020real.raster.list)
    names(stm.2020real) <- unlist(strsplit(stm.2020real.raster.list, "/|.tif"))[seq(4,88,4)]
    rm("stm.2020real.raster.list")
    
    stm.2020real <- stm.2020real[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    ## checking
    #names(stm.2020real)
    #plot(stm.2020real[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020real <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for real scenario in stm")
      
      BIOMOD_LoadModels(bm.out = myBiomodModelOut, full.name = m, as = "mod")
      
      temp_workdir = NULL
      if (length(grep("MAXENT.Phillips$", m)) == 1) {
        temp_workdir = mod@model_output_dir
      }
      
      indivdual_proj_2020real <- predict(mod, stm.2020real, temp_workdir = temp_workdir)
      names(indivdual_proj_2020real) <- m
      proj_2020real <- addLayer(proj_2020real, indivdual_proj_2020real)
      
    }
    
    
    
    # building a consensus map by mean weight
    proj_2020real.conbywm <- weighted.mean(proj_2020real, EMeval[EMeval$Eval.metric=="TSS","Testing.data"], na.rm=T)
    
    writeRaster(proj_2020real.conbywm, paste0("models.output/maps/STM/2020_real/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
    
    #cutoff.th <- round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2)
    #
    #proj_2020real.conbywm.bin <- proj_2020real.conbywm
    #proj_2020real.conbywm.bin[proj_2020real.conbywm.bin>=cutoff.th]<-1
    #proj_2020real.conbywm.bin[proj_2020real.conbywm.bin<cutoff.th]<-0
    #
    #writeRaster(proj_2020real.conbywm.bin, paste0("models.output/maps/STM/2020_real/", i, "_merged_algo_merged_dataset_merged_run_TSSBin.tif"), format = "GTiff", overwrite = T)
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "pa.data", "occur", "sel.var.df", "env.explanatory.var",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 avoid degradation
    stm.2020avoiddegrad.raster.list <- list.files("rasters/STM/2020_avoiddegrad/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020avoiddegrad <- stack(stm.2020avoiddegrad.raster.list)
    names(stm.2020avoiddegrad) <- unlist(strsplit(stm.2020avoiddegrad.raster.list, "/|.tif"))[seq(4,88,4)]
    rm("stm.2020avoiddegrad.raster.list")
    
    stm.2020avoiddegrad <- stm.2020avoiddegrad[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    ## checking
    #names(stm.2020avoiddegrad)
    #plot(stm.2020avoiddegrad[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020avoiddegrad <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for avoid degradation scenario in stm")
      
      BIOMOD_LoadModels(bm.out = myBiomodModelOut, full.name = m, as = "mod")
      
      temp_workdir = NULL
      if (length(grep("MAXENT.Phillips$", m)) == 1) {
        temp_workdir = mod@model_output_dir
      }
      
      indivdual_proj_2020avoiddegrad <- predict(mod, stm.2020avoiddegrad, temp_workdir = temp_workdir)
      names(indivdual_proj_2020avoiddegrad) <- m
      proj_2020avoiddegrad <- addLayer(proj_2020avoiddegrad, indivdual_proj_2020avoiddegrad)
      
    }
    
    
    
    # building a consensus map by mean weight
    proj_2020avoiddegrad.conbywm <- weighted.mean(proj_2020avoiddegrad, EMeval[EMeval$Eval.metric=="TSS","Testing.data"], na.rm=T)
    
    writeRaster(proj_2020avoiddegrad.conbywm, paste0("models.output/maps/STM/2020_avoiddegrad/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
    
    #cutoff.th <- round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2)
    #
    #proj_2020avoiddegrad.conbywm.bin <- proj_2020avoiddegrad.conbywm
    #proj_2020avoiddegrad.conbywm.bin[proj_2020avoiddegrad.conbywm.bin>=cutoff.th]<-1
    #proj_2020avoiddegrad.conbywm.bin[proj_2020avoiddegrad.conbywm.bin<cutoff.th]<-0
    #
    #writeRaster(proj_2020avoiddegrad.conbywm.bin, paste0("models.output/maps/STM/2020_avoiddegrad/", i, "_merged_algo_merged_dataset_merged_run_TSSBin.tif"), format = "GTiff", overwrite = T)
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "pa.data", "occur", "sel.var.df", "env.explanatory.var",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 avoid deforestation
    stm.2020avoiddeforest.raster.list <- list.files("rasters/STM/2020_avoiddeforest/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020avoiddeforest <- stack(stm.2020avoiddeforest.raster.list)
    names(stm.2020avoiddeforest) <- unlist(strsplit(stm.2020avoiddeforest.raster.list, "/|.tif"))[seq(4,88,4)]
    rm("stm.2020avoiddeforest.raster.list")
    
    stm.2020avoiddeforest <- stm.2020avoiddeforest[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    ## checking
    #names(stm.2020avoiddeforest)
    #plot(stm.2020avoiddeforest[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020avoiddeforest <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for avoid deforestation scenario in stm")
      
      BIOMOD_LoadModels(bm.out = myBiomodModelOut, full.name = m, as = "mod")
      
      temp_workdir = NULL
      if (length(grep("MAXENT.Phillips$", m)) == 1) {
        temp_workdir = mod@model_output_dir
      }
      
      indivdual_proj_2020avoiddeforest <- predict(mod, stm.2020avoiddeforest, temp_workdir = temp_workdir)
      names(indivdual_proj_2020avoiddeforest) <- m
      proj_2020avoiddeforest <- addLayer(proj_2020avoiddeforest, indivdual_proj_2020avoiddeforest)
      
    }
    
    
    
    # building a consensus map by mean weight
    proj_2020avoiddeforest.conbywm <- weighted.mean(proj_2020avoiddeforest, EMeval[EMeval$Eval.metric=="TSS","Testing.data"], na.rm=T)
    
    writeRaster(proj_2020avoiddeforest.conbywm, paste0("models.output/maps/STM/2020_avoiddeforest/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
    
    #cutoff.th <- round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2)
    #
    #proj_2020avoiddeforest.conbywm.bin <- proj_2020avoiddeforest.conbywm
    #proj_2020avoiddeforest.conbywm.bin[proj_2020avoiddeforest.conbywm.bin>=cutoff.th]<-1
    #proj_2020avoiddeforest.conbywm.bin[proj_2020avoiddeforest.conbywm.bin<cutoff.th]<-0
    #
    #writeRaster(proj_2020avoiddeforest.conbywm.bin, paste0("models.output/maps/STM/2020_avoiddeforest/", i, "_merged_algo_merged_dataset_merged_run_TSSBin.tif"), format = "GTiff", overwrite = T)
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "pa.data", "occur", "sel.var.df", "env.explanatory.var",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 avoid deforestation upf only
    stm.2020avoiddeforest2.raster.list <- list.files("rasters/STM/2020_avoiddeforest2/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020avoiddeforest2 <- stack(stm.2020avoiddeforest2.raster.list)
    names(stm.2020avoiddeforest2) <- unlist(strsplit(stm.2020avoiddeforest2.raster.list, "/|.tif"))[seq(4,88,4)]
    rm("stm.2020avoiddeforest2.raster.list")
    
    stm.2020avoiddeforest2 <- stm.2020avoiddeforest2[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    ## checking
    #names(stm.2020avoiddeforest2)
    #plot(stm.2020avoiddeforest2[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020avoiddeforest2 <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for avoid deforestation upf only scenario in stm")
      
      BIOMOD_LoadModels(bm.out = myBiomodModelOut, full.name = m, as = "mod")
      
      temp_workdir = NULL
      if (length(grep("MAXENT.Phillips$", m)) == 1) {
        temp_workdir = mod@model_output_dir
      }
      
      indivdual_proj_2020avoiddeforest2 <- predict(mod, stm.2020avoiddeforest2, temp_workdir = temp_workdir)
      names(indivdual_proj_2020avoiddeforest2) <- m
      proj_2020avoiddeforest2 <- addLayer(proj_2020avoiddeforest2, indivdual_proj_2020avoiddeforest2)
      
    }
    
    
    
    # building a consensus map by mean weight
    proj_2020avoiddeforest2.conbywm <- weighted.mean(proj_2020avoiddeforest2, EMeval[EMeval$Eval.metric=="TSS","Testing.data"], na.rm=T)
    
    writeRaster(proj_2020avoiddeforest2.conbywm, paste0("models.output/maps/STM/2020_avoiddeforest2/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
    
    #cutoff.th <- round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2)
    #
    #proj_2020avoiddeforest2.conbywm.bin <- proj_2020avoiddeforest2.conbywm
    #proj_2020avoiddeforest2.conbywm.bin[proj_2020avoiddeforest2.conbywm.bin>=cutoff.th]<-1
    #proj_2020avoiddeforest2.conbywm.bin[proj_2020avoiddeforest2.conbywm.bin<cutoff.th]<-0
    #
    #writeRaster(proj_2020avoiddeforest2.conbywm.bin, paste0("models.output/maps/STM/2020_avoiddeforest2/", i, "_merged_algo_merged_dataset_merged_run_TSSBin.tif"), format = "GTiff", overwrite = T)
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "pa.data", "occur", "sel.var.df", "env.explanatory.var",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 avoid both
    stm.2020avoidboth.raster.list <- list.files("rasters/STM/2020_avoidboth/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020avoidboth <- stack(stm.2020avoidboth.raster.list)
    names(stm.2020avoidboth) <- unlist(strsplit(stm.2020avoidboth.raster.list, "/|.tif"))[seq(4,88,4)]
    rm("stm.2020avoidboth.raster.list")
    
    stm.2020avoidboth <- stm.2020avoidboth[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    ## checking
    #names(stm.2020avoidboth)
    #plot(stm.2020avoidboth[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020avoidboth <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for avoid both scenario in stm")
      
      BIOMOD_LoadModels(bm.out = myBiomodModelOut, full.name = m, as = "mod")
      
      temp_workdir = NULL
      if (length(grep("MAXENT.Phillips$", m)) == 1) {
        temp_workdir = mod@model_output_dir
      }
      
      indivdual_proj_2020avoidboth <- predict(mod, stm.2020avoidboth, temp_workdir = temp_workdir)
      names(indivdual_proj_2020avoidboth) <- m
      proj_2020avoidboth <- addLayer(proj_2020avoidboth, indivdual_proj_2020avoidboth)
      
    }
    
    
    
    # building a consensus map by mean weight
    proj_2020avoidboth.conbywm <- weighted.mean(proj_2020avoidboth, EMeval[EMeval$Eval.metric=="TSS","Testing.data"], na.rm=T)
    
    writeRaster(proj_2020avoidboth.conbywm, paste0("models.output/maps/STM/2020_avoidboth/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
    
    #cutoff.th <- round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2)
    #
    #proj_2020avoidboth.conbywm.bin <- proj_2020avoidboth.conbywm
    #proj_2020avoidboth.conbywm.bin[proj_2020avoidboth.conbywm.bin>=cutoff.th]<-1
    #proj_2020avoidboth.conbywm.bin[proj_2020avoidboth.conbywm.bin<cutoff.th]<-0
    #
    #writeRaster(proj_2020avoidboth.conbywm.bin, paste0("models.output/maps/STM/2020_avoidboth/", i, "_merged_algo_merged_dataset_merged_run_TSSBin.tif"), format = "GTiff", overwrite = T)
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "pa.data", "occur", "sel.var.df", "env.explanatory.var",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 restoration without avoid
    stm.2020restor_wo_avoid.raster.list <- list.files("rasters/STM/2020_restor_wo_avoid/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020restor_wo_avoid <- stack(stm.2020restor_wo_avoid.raster.list)
    names(stm.2020restor_wo_avoid) <- unlist(strsplit(stm.2020restor_wo_avoid.raster.list, "/|.tif"))[seq(4,88,4)]
    rm("stm.2020restor_wo_avoid.raster.list")
    
    stm.2020restor_wo_avoid <- stm.2020restor_wo_avoid[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    ## checking
    #names(stm.2020restor_wo_avoid)
    #plot(stm.2020restor_wo_avoid[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020restor_wo_avoid <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for restoration without avoid scenario in stm")
      
      BIOMOD_LoadModels(bm.out = myBiomodModelOut, full.name = m, as = "mod")
      
      temp_workdir = NULL
      if (length(grep("MAXENT.Phillips$", m)) == 1) {
        temp_workdir = mod@model_output_dir
      }
      
      indivdual_proj_2020restor_wo_avoid <- predict(mod, stm.2020restor_wo_avoid, temp_workdir = temp_workdir)
      names(indivdual_proj_2020restor_wo_avoid) <- m
      proj_2020restor_wo_avoid <- addLayer(proj_2020restor_wo_avoid, indivdual_proj_2020restor_wo_avoid)
      
    }
    
    
    
    # building a consensus map by mean weight
    proj_2020restor_wo_avoid.conbywm <- weighted.mean(proj_2020restor_wo_avoid, EMeval[EMeval$Eval.metric=="TSS","Testing.data"], na.rm=T)
    
    writeRaster(proj_2020restor_wo_avoid.conbywm, paste0("models.output/maps/STM/2020_restor_wo_avoid/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
    
    #cutoff.th <- round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2)
    #
    #proj_2020restor_wo_avoid.conbywm.bin <- proj_2020restor_wo_avoid.conbywm
    #proj_2020restor_wo_avoid.conbywm.bin[proj_2020restor_wo_avoid.conbywm.bin>=cutoff.th]<-1
    #proj_2020restor_wo_avoid.conbywm.bin[proj_2020restor_wo_avoid.conbywm.bin<cutoff.th]<-0
    #
    #writeRaster(proj_2020restor_wo_avoid.conbywm.bin, paste0("models.output/maps/STM/2020_restor_wo_avoid/", i, "_merged_algo_merged_dataset_merged_run_TSSBin.tif"), format = "GTiff", overwrite = T)
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "pa.data", "occur", "sel.var.df", "env.explanatory.var",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 restoration and avoid deforestation
    stm.2020restor_n_avoid_deforest.raster.list <- list.files("rasters/STM/2020_restor_n_avoid_deforest/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020restor_n_avoid_deforest <- stack(stm.2020restor_n_avoid_deforest.raster.list)
    names(stm.2020restor_n_avoid_deforest) <- unlist(strsplit(stm.2020restor_n_avoid_deforest.raster.list, "/|.tif"))[seq(4,88,4)]
    rm("stm.2020restor_n_avoid_deforest.raster.list")
    
    stm.2020restor_n_avoid_deforest <- stm.2020restor_n_avoid_deforest[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    ## checking
    #names(stm.2020restor_n_avoid_deforest)
    #plot(stm.2020restor_n_avoid_deforest[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020restor_n_avoid_deforest <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "restoration and avoid scenario in stm")
      
      BIOMOD_LoadModels(bm.out = myBiomodModelOut, full.name = m, as = "mod")
      
      temp_workdir = NULL
      if (length(grep("MAXENT.Phillips$", m)) == 1) {
        temp_workdir = mod@model_output_dir
      }
      
      indivdual_proj_2020restor_n_avoid_deforest <- predict(mod, stm.2020restor_n_avoid_deforest, temp_workdir = temp_workdir)
      names(indivdual_proj_2020restor_n_avoid_deforest) <- m
      proj_2020restor_n_avoid_deforest <- addLayer(proj_2020restor_n_avoid_deforest, indivdual_proj_2020restor_n_avoid_deforest)
      
    }
    
    
    
    # building a consensus map by mean weight
    proj_2020restor_n_avoid_deforest.conbywm <- weighted.mean(proj_2020restor_n_avoid_deforest, EMeval[EMeval$Eval.metric=="TSS","Testing.data"], na.rm=T)
    
    writeRaster(proj_2020restor_n_avoid_deforest.conbywm, paste0("models.output/maps/STM/2020_restor_n_avoid_deforest/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
    
    #cutoff.th <- round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2)
    #
    #proj_2020restor_n_avoid_deforest.conbywm.bin <- proj_2020restor_n_avoid_deforest.conbywm
    #proj_2020restor_n_avoid_deforest.conbywm.bin[proj_2020restor_n_avoid_deforest.conbywm.bin>=cutoff.th]<-1
    #proj_2020restor_n_avoid_deforest.conbywm.bin[proj_2020restor_n_avoid_deforest.conbywm.bin<cutoff.th]<-0
    #
    #writeRaster(proj_2020restor_n_avoid_deforest.conbywm.bin, paste0("models.output/maps/STM/2020_restor_n_avoid_deforest/", i, "_merged_algo_merged_dataset_merged_run_TSSBin.tif"), format = "GTiff", overwrite = T)
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "pa.data", "occur", "sel.var.df", "env.explanatory.var",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 restoration and avoid both
    stm.2020restor_n_avoid_both.raster.list <- list.files("rasters/STM/2020_restor_n_avoid_both/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020restor_n_avoid_both <- stack(stm.2020restor_n_avoid_both.raster.list)
    names(stm.2020restor_n_avoid_both) <- unlist(strsplit(stm.2020restor_n_avoid_both.raster.list, "/|.tif"))[seq(4,88,4)]
    rm("stm.2020restor_n_avoid_both.raster.list")
    
    stm.2020restor_n_avoid_both <- stm.2020restor_n_avoid_both[[sel.var.df[!is.na(sel.var.df$VIF) & sel.var.df$Type=="environment","VAR"]]]
    ## checking
    #names(stm.2020restor_n_avoid_both)
    #plot(stm.2020restor_n_avoid_both[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020restor_n_avoid_both <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "restoration and avoid scenario in stm")
      
      BIOMOD_LoadModels(bm.out = myBiomodModelOut, full.name = m, as = "mod")
      
      temp_workdir = NULL
      if (length(grep("MAXENT.Phillips$", m)) == 1) {
        temp_workdir = mod@model_output_dir
      }
      
      indivdual_proj_2020restor_n_avoid_both <- predict(mod, stm.2020restor_n_avoid_both, temp_workdir = temp_workdir)
      names(indivdual_proj_2020restor_n_avoid_both) <- m
      proj_2020restor_n_avoid_both <- addLayer(proj_2020restor_n_avoid_both, indivdual_proj_2020restor_n_avoid_both)
      
    }
    
    
    
    # building a consensus map by mean weight
    proj_2020restor_n_avoid_both.conbywm <- weighted.mean(proj_2020restor_n_avoid_both, EMeval[EMeval$Eval.metric=="TSS","Testing.data"], na.rm=T)
    
    writeRaster(proj_2020restor_n_avoid_both.conbywm, paste0("models.output/maps/STM/2020_restor_n_avoid_both/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
    
    #cutoff.th <- round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2)
    #
    #proj_2020restor_n_avoid_both.conbywm.bin <- proj_2020restor_n_avoid_both.conbywm
    #proj_2020restor_n_avoid_both.conbywm.bin[proj_2020restor_n_avoid_both.conbywm.bin>=cutoff.th]<-1
    #proj_2020restor_n_avoid_both.conbywm.bin[proj_2020restor_n_avoid_both.conbywm.bin<cutoff.th]<-0
    #
    #writeRaster(proj_2020restor_n_avoid_both.conbywm.bin, paste0("models.output/maps/STM/2020_restor_n_avoid_both/", i, "_merged_algo_merged_dataset_merged_run_TSSBin.tif"), format = "GTiff", overwrite = T)
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "pa.data", "occur", "sel.var.df", "env.explanatory.var",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
  }  
  
  
  #
  
  forestdep.spplist[forestdep.spplist$Binomial==i,"Done"] <- TRUE
  write.csv(forestdep.spplist, "models.output/evaluation/updated_species_summary_v1.csv")
  
  rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "pa.data", "occur", "sel.var.df", "env.explanatory.var"))])
  
  
  #do.call(file.remove, list(list.files(paste0("models.output/", i), full.names = T, recursive = T)))
  unlink(paste0("models.output/", i, "/*"), recursive = T, force = T)
  
}

#
#
