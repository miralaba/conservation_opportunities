
#' @title Cost-effectiveness of conservation actions in Amazon
#' @description 

#### setting working directory ####
#setwd("")
#options(expression=200000)
#memory.limit(1000000)

##### loading required packages ####
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
library(gridExtra)
library(grDevices)
library(doParallel)

#### presetting ####



#### import data ####
## transect data
transectdata <- read.csv("data/RAS_transects_environment_all.csv")
# keeping only PGM
pgm.transectdata <- transectdata[transectdata$Region=="PGM",]

## birds
birddata <- read.csv("data/Bird_data_Standard_NM_26022013_Final.csv")
# rename columns to match with transect data
names(birddata)[1] <- "Region"
names(birddata)[3] <- "Catchment"
names(birddata)[5] <- "Transectcode"
# new var to mark data frame as birds
birddata$Group <- "Birds"
# keeping only PGM
pgm.birddata <- birddata[birddata$Region=="PGM",]

## trees
pgm.treedata <- read.csv("data/Flora.composition.and.biomass_PGM_Erika_23.01.2013.csv")
# include region code
pgm.treedata$Region <- "PGM"
# new var to mark data frame as birds
pgm.treedata$Group <- "Trees"
# new var with code for transect by catchment
pgm.treedata$Transectcode <- paste(pgm.treedata$Catchment,pgm.treedata$Transect,sep="_")
# new var with spp binomial
pgm.treedata$Binomial <- paste(pgm.treedata$Genera,pgm.treedata$Species,sep="")


## merging spp data with binomial, location and lulc code
pgm.transectdata <- pgm.transectdata[,c(1:6,10)]
pgm.birddata <- pgm.birddata[,c(23,1,3,4,5,15)]
pgm.treedata <- pgm.treedata[,c(36,35,1,2,37,38)]
pgm.sppdata <- rbind(pgm.birddata,pgm.treedata)
pgm.sppdata <- merge(pgm.sppdata, pgm.transectdata)
## checking
#str(pgm.sppdata)
#View(pgm.sppdata)

# excluding Varzea transects PGM
exclude <- c("100_1", "100_4","100_7","81_12","423_2") #transect code
pgm.sppdata <- pgm.sppdata[!pgm.sppdata$Transectcode %in% exclude,]

# new var with longlat
utm1 <- data.frame(x=pgm.sppdata$UTM_X, y=pgm.sppdata$UTM_Y) 
coordinates(utm1) <- ~x+y 
#class(utm1)
proj4string(utm1) <- CRS("+proj=utm +zone=23 +south +datum=WGS84 +units=m +ellps=WGS84") 
utm2 <- spTransform(utm1, CRS("+proj=longlat +datum=WGS84"))

pgm.sppdata$Longitude <- utm2@coords[,1]
pgm.sppdata$Latitude <- utm2@coords[,2]

#
## filter 1: selecting forest dependent species using undisturbed primary forest layer
pgm.2010.UFPpx <- raster("rasters/PGM/2010/UPFpx.tif")

pgm.2010.UFPpx.90more <- pgm.2010.UFPpx
pgm.2010.UFPpx.90more[pgm.2010.UFPpx.90more>=.9]<-1
pgm.2010.UFPpx.90more[pgm.2010.UFPpx.90more<.9]<-0
## checking
#plot(pgm.2010.UFPpx.90more)
#points(SpatialPoints(pgm.sppdata[,c("Longitude", "Latitude")]))
pgm.sppdata$forestdep <- extract(pgm.2010.UFPpx.90more, SpatialPoints(pgm.sppdata[,c("Longitude", "Latitude")]))
pgm.forestdep.unique.spplist <- unique(pgm.sppdata[pgm.sppdata$forestdep==1,"Binomial"])
pgm.forestdep.sppdata <- pgm.sppdata[pgm.sppdata$Binomial %in% pgm.forestdep.unique.spplist,]

pgm.forestdep.spplist <- as.data.frame(table(pgm.forestdep.sppdata$Binomial))


#
## filter 2: excluding species with low records
pgm.forestdep.spplist <- pgm.forestdep.spplist[pgm.forestdep.spplist$Freq>=10,]
names(pgm.forestdep.spplist) <- c("Binomial", "Nrec")
pgm.forestdep.spplist <- merge(pgm.forestdep.spplist, pgm.sppdata[,c(6, 5)])
pgm.forestdep.spplist <- pgm.forestdep.spplist[!duplicated(pgm.forestdep.spplist),]

pgm.sppdata.final <- pgm.sppdata[pgm.sppdata$Binomial %in% pgm.forestdep.spplist$Binomial,]
pgm.sppdata.final <- pgm.sppdata.final[,c(1:4,9,5:8,10,11)]
## checking
#nrow(as.data.frame(table(pgm.sppdata.final[pgm.sppdata.final$Group=="Trees","Binomial"])))
#nrow(as.data.frame(table(pgm.sppdata.final[pgm.sppdata.final$Group=="Birds","Binomial"])))

write.csv(pgm.forestdep.spplist, "forest_dependent_species_list11.csv")
#

rm(list= ls()[!(ls() %in% c("pgm.forestdep.spplist", "pgm.sppdata.final"))])
gc()
#

# add new var for ensemble model scores
pgm.forestdep.spplist$ROC<-NA
pgm.forestdep.spplist$TSS<-NA
pgm.forestdep.spplist$Sensitivity<-NA
pgm.forestdep.spplist$Specificity<-NA
pgm.forestdep.spplist$Cutoff<-NA

#cl <- makeCluster(16)
#doParallel::registerDoParallel(cl)
#foreach(i=iter(pgm.forestdep.spplist$Binomial),.packages=c("raster","biomod2","spThin","usdm"),.errorhandling="stop")%dopar%{  

#i <- as.character(pgm.forestdep.spplist$Binomial[1])
for (i in pgm.forestdep.spplist$Binomial[1:5]) {
  # filter 3: excluding records to avoid spatial autocorrelation
  # filter one species at a time
  spp_x <- pgm.sppdata.final[pgm.sppdata.final$Binomial==i,]
  
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
  
  pgm.forestdep.spplist[pgm.forestdep.spplist$Binomial==i,"Nrec"] <- nrow(occur)
  write.csv(pgm.forestdep.spplist, "forest_dependent_species_list11.csv")
  
  #if(nrow(occur)<10) next
  nrow(occur)<10
  
  # importing exploratory variables
  pgm.2010.raster.list <- list.files("rasters/PGM/2010/", pattern = ".tif", full.names = T, recursive = T)
  pgm.2010 <- stack(pgm.2010.raster.list)
  rm("pgm.2010.raster.list")
  ## checking
  #names(pgm.2010)
  #plot(pgm.2010[[17]])
  #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
  
  # removing colinearity between exploratory variables
  sel.var <- vifcor(pgm.2010, th=0.7, maxobservations=10000)
  sel.var.df <- data.frame(rbind(cbind(VAR=sel.var@results$Variables, VIF=sel.var@results$VIF),
                                 cbind(VAR=sel.var@excluded, VIF=NA)))
  
  write.csv(sel.var.df, paste0("selected_exploratory_var/", i,"_VIF.csv", sep=""), row.names = F)
  
  pgm.2010 <- pgm.2010[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
  
  rm(list= ls()[!(ls() %in% c("pgm.forestdep.spplist", "pgm.sppdata.final", "i", "spp_x", "occur", "sel.var.df", "pgm.2010"))])
  gc()
  
  ## input
  set.seed(3711)
  myBiomodData <- BIOMOD_FormatingData(resp.var = rep.int(1, times = nrow(occur)),
                                       expl.var = pgm.2010,
                                       resp.xy = occur[,c("Longitude", "Latitude")],
                                       resp.name = i,
                                       PA.nb.rep = 3,
                                       PA.nb.absences = as.numeric(nrow(occur)*10),
                                       PA.strategy = "disk",
                                       PA.dist.min = 5000)
  
  ## checking
  #myBiomodData
  #nrow(myBiomodData@coord)
  #head(myBiomodData@coord, xx)
  #plot(pgm.2010[[1]])
  #points(myBiomodData@coord[c(1:xx),], col="red")
  #points(myBiomodData@coord[c(xx:nrow(myBiomodData@coord)),])
  
  # algorithms control
  myBiomodOption <- BIOMOD_ModelingOptions(#GLM = list(control = glm.control(epsilon = 1e-08, maxit = 1000)),
                                           RF = list(ntree=1000, nodesize=10),
                                           MAXENT.Phillips = list(memory_allocated = 4096)) #12288
  #
  # model fit parameters
  set.seed(3711)
  myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                      modeling.id = i,
                                      models = c('RF', 'MAXENT.Phillips'), #'GLM', 
                                      bm.options = myBiomodOption,
                                      nb.rep = 10,
                                      data.split.perc = 80,
                                      do.full.models = F,
                                      metric.eval = c('TSS', 'ROC'),
                                      var.import = 5,
                                      nb.cpu = 3)
  
  # saving
  capture.output(myBiomodData, file = paste0(i ,"/", i,"_BiomodData.txt"))
  capture.output(myBiomodModelOut, file = paste0(i, "/", i,"_BiomodModelOut.txt"))
  capture.output(get_evaluations(myBiomodModelOut), file = paste0(i, "/", i,"_eval_BiomodModelOut.txt"))
  capture.output(get_variables_importance(myBiomodModelOut), file = paste0(i, "/", i,"_var_importance.txt"))
  
  # adding var mean variable importance to selected exploratory variables data frame
  varimport <- get_variables_importance(myBiomodModelOut)
  sel.var.df$varimport <- NA
  sel.var.df[!is.na(sel.var.df$VIF), "varimport"] <- rowMeans(varimport, na.rm=T)
  write.csv(sel.var.df, paste0("selected_exploratory_var/", i,"_VIF.csv", sep=""), row.names = F)
  
  # graphic of model scores
  pdf(paste(i, "/", i ,"_models_scores.pdf", sep="" ))
  gEval <- bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'run'))
  dev.off()
  
  ####################################
  ## caso o modelo ja esteja pronto ##
  ## e precises apenas projetar     ##
  ## os mapas, comecar daqui        ##
  ####################################
  #sel.var.df <- read.csv(paste0("selected_exploratory_var/", i, "_VIF.csv"), header=T)
  #myBiomodModelOut <- load(paste0(i,"/",i,".",i,".models.out"))
  #myBiomodModelOut <- get(myBiomodModelOut)
  
  
  # updating evaluation threshold
  eval_threshold <- 0.7
  
  ev<-get_evaluations(myBiomodModelOut, as.data.frame=T)
  
  #selected_GLM <- subset (ev$Testing.data, grepl("GLM", ev$Model.name) & grepl("TSS", ev$Eval.metric)) 
  selected_RF <- subset (ev$Testing.data, grepl("RF", ev$Model.name) & grepl("TSS", ev$Eval.metric))
  selected_MAX <- subset (ev$Testing.data, grepl("MAXENT.Phillips", ev$Model.name) & grepl("TSS", ev$Eval.metric)) 
  
  new.eval_threshold <- min(c(#max(selected_GLM[!is.na(selected_GLM)]),
                              max(selected_RF[!is.na(selected_RF)]),
                              max(selected_MAX[!is.na(selected_MAX)])))
  
  if (new.eval_threshold < 0.7) {eval_threshold = floor(new.eval_threshold*10)/10}
  
  
  all.models <- get_built_models(myBiomodModelOut)
  chosen.models.from.ev <- ev[ev$Eval.metric=="TSS" & ev$Testing.data>eval_threshold,]
  chosen.models.from.ev$Model.name2 <- paste(chosen.models.from.ev$Dataset, chosen.models.from.ev$Run, chosen.models.from.ev$Algo, sep="_")
  chosen.models.from.all.models <- grep(paste(chosen.models.from.ev$Model.name2,collapse="|"), all.models, value=TRUE)

  #if(nrow(chosen.models.from.ev)<5) next
  nrow(chosen.models.from.ev)<5
  
  ## ensemble: buiding a consensus between data sets (pseudo-absence), algorithms and runs
  #myBiomodEM_algo <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
  #                                          models.chosen = 'all',
  #                                          em.by = 'algo',
  #                                          metric.select = c('TSS'),
  #                                          metric.select.thresh = eval_threshold,
  #                                          metric.eval = c('TSS', 'ROC'),
  #                                          var.import = 0,
  #                                          prob.mean = F,
  #                                          prob.cv = F,
  #                                          prob.ci = F,
  #                                          prob.ci.alpha = 0.05,
  #                                          prob.median = F,
  #                                          committee.averaging = F,
  #                                          prob.mean.weight = T,
  #                                          prob.mean.weight.decay = 'proportional',
  #                                          nb.cpu = 3,
  #                                          seed.val = 712)
  #
  ## saving
  #capture.output(myBiomodEM_all, file = paste0(i, "/", i,"_EM_all.txt"))
  #capture.output(get_evaluations(myBiomodEM_all), file = paste0(i, "/", i,"_eval_EM_all.txt"))
  #
  # ensemble model scores to species list
  #EMeval <- as.data.frame(get_evaluations(myBiomodEM_all))
  EMeval <- ev[ev$Model.name %in% chosen.models.from.ev$Model.name,]
  pgm.forestdep.spplist[pgm.forestdep.spplist$Binomial==i,"ROC"] <- paste(round(mean(EMeval[EMeval$Eval.metric=="ROC","Testing.data"]),1), "+/-", round(sd(EMeval[EMeval$Eval.metric=="ROC","Testing.data"]),2))
  pgm.forestdep.spplist[pgm.forestdep.spplist$Binomial==i,"TSS"] <- paste(round(mean(EMeval[EMeval$Eval.metric=="TSS","Testing.data"]),1), "+/-", round(sd(EMeval[EMeval$Eval.metric=="TSS","Testing.data"]),2))
  pgm.forestdep.spplist[pgm.forestdep.spplist$Binomial==i,"Sensitivity"] <- paste(round(mean(EMeval[EMeval$Eval.metric=="TSS","Sensitivity"]),2), "+/-", round(sd(EMeval[EMeval$Eval.metric=="TSS","Sensitivity"]),2))
  pgm.forestdep.spplist[pgm.forestdep.spplist$Binomial==i,"Specificity"] <- paste(round(mean(EMeval[EMeval$Eval.metric=="TSS","Specificity"]),2), "+/-", round(sd(EMeval[EMeval$Eval.metric=="TSS","Specificity"]),2))
  pgm.forestdep.spplist[pgm.forestdep.spplist$Binomial==i,"Cutoff"] <- paste(round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2), "+/-", round((sd(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2))
  
  write.csv(pgm.forestdep.spplist, "forest_dependent_species_list11.csv")
  
  ## projecting individual maps to 2010
  #myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
  #                                  proj.name = '2010',
  #                                  new.env = pgm.2010,
  #                                  models.chosen = 'all',
  #                                  metric.binary = 'TSS',
  #                                  compress = F,
  #                                  build.clamping.mask = F,
  #                                  output.format = '.img',
  #                                  do.stack = F)
  #                                  #nb.cpu = 3,
  #                                  #seed.val = 712)
  
  #dir.create(paste0(i ,"/proj_2010"))
  
  proj_2010 <- stack()
  for(m in chosen.models.from.all.models){
    
    cat("\n\t> Projecting", m, "...")
    
    BIOMOD_LoadModels(bm.out = myBiomodModelOut, full.name = m, as = "mod")
    
    temp_workdir = NULL
    if (length(grep("MAXENT.Phillips$", m)) == 1) {
      temp_workdir = mod@model_output_dir
    }
    
    indivdual_proj_2010 <- predict(mod, pgm.2010, temp_workdir = temp_workdir)
    names(indivdual_proj_2010) <- m
    proj_2010 <- addLayer(proj_2010, indivdual_proj_2010)
    
  }
  
  
  
  ## projecting ensemble
  #myBiomodProj_EM_all <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM_all,
  #                                                  bm.proj = myBiomodProj,
  #                                                  #new.env = NULL,
  #                                                  #xy.new.env = NULL,
  #                                                  models.chosen = 'all',
  #                                                  proj.name = '2010_consenso',
  #                                                  metric.binary = 'TSS',
  #                                                  metric.filter = NULL,
  #                                                  compress = NULL,
  #                                                  #output.format = '.img',
  #                                                  total.consensus = T,
  #                                                  nb.cpu = 1)
  
  
  # building a consensus map by mean weight
  dir.create(paste0(i ,"/proj_2010_consenso"))
  
  proj_2010.conbywm <- weighted.mean(proj_2010, EMeval[EMeval$Eval.metric=="TSS","Testing.data"], na.rm=T)
  
  #writeRaster(proj_2010.conbywm, paste0(i,"/proj_2010_consenso/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
  
  cutoff.th <- round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2)
  
  proj_2010.conbywm.bin <- proj_2010.conbywm
  proj_2010.conbywm.bin[proj_2010.conbywm.bin>=cutoff.th]<-1
  proj_2010.conbywm.bin[proj_2010.conbywm.bin<cutoff.th]<-0
  
  writeRaster(proj_2010.conbywm.bin, paste0(i,"/proj_2010_consenso/", i, "_merged_algo_merged_dataset_merged_run_TSSBin.tif"), format = "GTiff", overwrite = T)
  
  
  rm(list= ls()[!(ls() %in% c("pgm.forestdep.spplist", "pgm.sppdata.final", "i", "spp_x", "occur", "sel.var.df", 
                              "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                              "chosen.models.from.ev", "chosen.models.from.all.models", "proj_2010"))])
  gc()

  # uncertainties 
  # pixel global mean
  Mtotal <- mean(proj_2010)
  # pixel mean by algorithm
  #Mglm <- mean(proj_2010[[grep("GLM", names(proj_2010), value = T)]])
  Mrf <- mean(proj_2010[[grep("RF", names(proj_2010), value = T)]])
  Mmaxent <- mean(proj_2010[[grep("MAXENT.Phillips", names(proj_2010), value = T)]])
  
  # total square sum
  SST <- sum((proj_2010-Mtotal)^2)
  
  # square sum of algorithm
  SSMET <- #sum((proj_2010[[grep("GLM", names(proj_2010), value = T)]]-Mglm)^2)+
    sum((proj_2010[[grep("RF", names(proj_2010), value = T)]]-Mrf)^2)+
    sum((proj_2010[[grep("MAXENT.Phillips", names(proj_2010), value = T)]]-Mmaxent)^2)
  
  # variance proportion between algorithms
  between.met <- SSMET/SST
  
  # square sum and variance proportion within glm
  #SSglm <-sum((proj_2010[[grep("GLM", names(proj_2010), value = T)]]-Mglm)^2)
  #within.glm <- SSglm/SST
  
  # square sum and variance proportion within rf
  SSrf <-sum((proj_2010[[grep("RF", names(proj_2010), value = T)]]-Mrf)^2)
  within.rf <- SSrf/SST
  # square sum and variance proportion within maxent
  SSmaxent <- sum((proj_2010[[grep("MAXENT.Phillips", names(proj_2010), value = T)]]-Mmaxent)^2)
  within.maxent <- SSmaxent/SST
  
  uncert.part<-stack(between.met, within.rf, within.maxent) #within.glm, 
  names(uncert.part)<-c("between_algo", "within_rf", "within_maxent") #"within_glm", 
  #plot(uncert.part)
  
  writeRaster(uncert.part, paste0(i,"/proj_2010_consenso/uncertainty2010.tif"), format="GTiff")
  
  rm(list= ls()[(ls() %in% c("proj_2010", "Mtotal", "Mglm", "Mrf", "Mmaxent", "SST", "SSMET", "SSglm", "SSrf", "SSmaxent",
                             "between.met", "within.glm", "within.rf", "within.maxent", "uncert.part"))])
  gc()
  
  # 2020
  pgm.2020.raster.list <- list.files("rasters/PGM/2020/", pattern = ".tif", full.names = T, recursive = T)
  pgm.2020 <- stack(pgm.2020.raster.list)
  rm("pgm.2020.raster.list")
  ## checking
  #names(pgm.2020)
  #plot(pgm.2020[[17]])
  #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
  
  pgm.2020 <- pgm.2020[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
  
  ## projecting individual maps to 2020
  #myBiomodProj_F1 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
  #                                     new.env = pgm.2020,
  #                                     proj.name = '2020',
  #                                     xy.new.env = NULL,
  #                                     selected.models = "all",
  #                                     binary.meth = 'TSS',
  #                                     compress = F,
  #                                     build.clamping.mask = F,
  #                                     do.stack=F,
  #                                     output.format = '.img')
  
  #dir.create(paste0(i ,"/proj_2020"))
  
  proj_2020 <- stack()
  for(m in chosen.models.from.all.models){
    
    cat("\n\t> Projecting", m, "...")
    
    BIOMOD_LoadModels(bm.out = myBiomodModelOut, full.name = m, as = "mod")
    
    temp_workdir = NULL
    if (length(grep("MAXENT.Phillips$", m)) == 1) {
      temp_workdir = mod@model_output_dir
    }
    
    indivdual_proj_2020 <- predict(mod, pgm.2020, temp_workdir = temp_workdir)
    names(indivdual_proj_2020) <- m
    proj_2020 <- addLayer(proj_2020, indivdual_proj_2020)
    
  }
  
  
  
  ## projecting ensemble
  #myBiomodProj_EM_F1_all <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_all,
  #                                                     projection.output = myBiomodProj_F1,
  #                                                     new.env = NULL,
  #                                                     xy.new.env = NULL,
  #                                                     selected.models = 'all',
  #                                                     proj.name = '2020_consenso',
  #                                                     binary.meth = 'TSS',
  #                                                     filtered.meth = NULL,
  #                                                     compress = NULL,
  #                                                     output.format = '.img',
  #                                                     total.consensus = TRUE)
  
  
  # building a consensus map by mean weight
  dir.create(paste0(i ,"/proj_2020_consenso"))
  
  proj_2020.conbywm <- weighted.mean(proj_2020, EMeval[EMeval$Eval.metric=="TSS","Testing.data"], na.rm=T)
  
  #writeRaster(proj_2020.conbywm, paste0(i,"/proj_2020_consenso/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
  
  cutoff.th <- round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2)
  
  proj_2020.conbywm.bin <- proj_2020.conbywm
  proj_2020.conbywm.bin[proj_2020.conbywm.bin>=cutoff.th]<-1
  proj_2020.conbywm.bin[proj_2020.conbywm.bin<cutoff.th]<-0
  
  writeRaster(proj_2020.conbywm.bin, paste0(i,"/proj_2020_consenso/", i, "_merged_algo_merged_dataset_merged_run_TSSBin.tif"), format = "GTiff", overwrite = T)
  
  
  rm(list= ls()[!(ls() %in% c("pgm.forestdep.spplist", "pgm.sppdata.final", "i", "spp_x", "occur", "sel.var.df",
                              "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                              "chosen.models.from.ev", "chosen.models.from.all.models", "proj_2020"))])
  gc()
  
  # uncertainties 
  # pixel global mean
  Mtotal <- mean(proj_2020)
  # pixel mean by algorithm
  #Mglm <- mean(proj_2020[[grep("GLM", names(proj_2020), value = T)]])
  Mrf <- mean(proj_2020[[grep("RF", names(proj_2020), value = T)]])
  Mmaxent <- mean(proj_2020[[grep("MAXENT.Phillips", names(proj_2020), value = T)]])
  
  # total square sum
  SST <- sum((proj_2020-Mtotal)^2)
  
  # square sum of algorithm
  SSMET <- #sum((proj_2020[[grep("GLM", names(proj_2020), value = T)]]-Mglm)^2)+
    sum((proj_2020[[grep("RF", names(proj_2020), value = T)]]-Mrf)^2)+
    sum((proj_2020[[grep("MAXENT.Phillips", names(proj_2020), value = T)]]-Mmaxent)^2)
  
  # variance proportion between algorithms
  between.met <- SSMET/SST
  
  # square sum and variance proportion within glm
  #SSglm <-sum((proj_2020[[grep("GLM", names(proj_2020), value = T)]]-Mglm)^2)
  #within.glm <- SSglm/SST
  
  # square sum and variance proportion within rf
  SSrf <-sum((proj_2020[[grep("RF", names(proj_2020), value = T)]]-Mrf)^2)
  within.rf <- SSrf/SST
  # square sum and variance proportion within maxent
  SSmaxent <- sum((proj_2020[[grep("MAXENT.Phillips", names(proj_2020), value = T)]]-Mmaxent)^2)
  within.maxent <- SSmaxent/SST
  
  uncert.part<-stack(between.met, within.rf, within.maxent) #within.glm, 
  names(uncert.part)<-c("between_algo", "within_rf", "within_maxent") #"within_glm", 
  #plot(uncert.part)
  
  writeRaster(uncert.part, paste0(as.character(i),"/proj_2020_consenso/uncertainty2020.grd"), format="raster")
  
  
  rm(list= ls()[!(ls() %in% c("pgm.forestdep.spplist", "pgm.sppdata.final", "i"))])
  gc()
  
  #

do.call(file.remove, list(list.files(paste0(i ,"/.BIOMOD_DATA"), full.names = TRUE, recursive = TRUE)))
do.call(file.remove, list(list.files(paste0(i ,"/models"), full.names = TRUE, recursive = TRUE)))

  
}
#}
#stopCluster()
#gc()
   







   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
##########################################################
##########################################################
## filter 1: selectinh forest dependent species using NMDS
# keeping only forest sites
pgm.forest.sppdata <- pgm.sppdata[!pgm.sppdata$LU_FT_Code %in% c("PA", "MA", "AP", "REF", "SHA"),]
#View(pgm.forest.sppdata)

# arrange [melt] data into sampling-by-species matrix
pgm.transect.species.matrix <- cast(pgm.forest.sppdata, Transectcode~Binomial, value = "Latitude")
rownames(pgm.transect.species.matrix) <- pgm.transect.species.matrix[,1]
pgm.transect.species.matrix <- pgm.transect.species.matrix[,2:ncol(pgm.transect.species.matrix)]
#View(pgm.transect.species.matrix)

# two-dimensional non-metric multidimensional scaling (MDS)
pgm.NMDS <- metaMDS(pgm.transect.species.matrix, k=2, trymax = 999, autotransform = F)
## checking
#goodness(pgm.NMDS)
#stressplot(pgm.NMDS)

# defining LULC codes
treat <- pgm.forest.sppdata[!duplicated(pgm.forest.sppdata$Transectcode),c(4,8)]
#rownames(treat) <- treat[,1]
#treat <- treat[,2:ncol(treat)]

# NMDS values to data frame::sites
pgm.NMDS.sites.df <- as.data.frame(scores(pgm.NMDS, display = "sites"))
pgm.NMDS.sites.df$Transectcode <- rownames(pgm.transect.species.matrix)
pgm.NMDS.sites.df$LU_FT_Code <- treat[,2]
#head(pgm.NMDS.sites.df)

# NMDS values to data frame::species
pgm.NMDS.species.df <- as.data.frame(scores(pgm.NMDS, display = "species"))
pgm.NMDS.species.df$Binomial <- colnames(pgm.transect.species.matrix)
pgm.NMDS.species.df$spID <- 1:nrow(pgm.NMDS.species.df)
#head(pgm.NMDS.species.df)

# NMDS values to data frame::convex hull of undisturbed forests &
# selecting species inside convex hull
X11()
plot(pgm.NMDS.species.df[, c("NMDS1","NMDS2")])
pgm.NMDS.sites.hull <- with(treat, ordihull(pgm.NMDS, LU_FT_Code, scaling = "symmetric", label = TRUE))
selectedPoints <- gatepoints::fhs(pgm.NMDS.species.df[, c("NMDS1","NMDS2")])
pgm.forest.sppdata.list <- pgm.NMDS.species.df[selectedPoints,]

pgm.forest.sppdata <- pgm.sppdata[pgm.sppdata$Binomial %in% pgm.forest.sppdata.list$Binomial,]
#
   
   