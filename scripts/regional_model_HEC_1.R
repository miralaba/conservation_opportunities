
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
library(spThin)
library(usdm)
library(biomod2)
library(gridExtra)
library(grDevices)

#### loading input data ####
#explanatory variables
sel.var.df <- read.csv("rasters/selected_environmental_explanatory_variables_byVIF.csv")

pgm.env.explanatory.var.list <- list.files("rasters/PGM/2010_real", pattern = ".tif", full.names = T, recursive = T)

pgm.env.explanatory.var <- stack(pgm.env.explanatory.var.list)
names(pgm.env.explanatory.var) <- unlist(strsplit(pgm.env.explanatory.var.list, "/|.tif"))[seq(4,80,4)]

pgm.env.explanatory.var <- pgm.env.explanatory.var[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
pgm.env.explanatory.var<-aggregate(pgm.env.explanatory.var,fact=35,fun=mean,na.rm=T); pgm.env.explanatory.var<-stack(pgm.env.explanatory.var)
##cheking
#pgm.env.explanatory.var
#plot(pgm.env.explanatory.var[[1:10]], nc=2)
#plot(pgm.env.explanatory.var[[11:20]], nc=2)



stm.env.explanatory.var.list <- list.files("rasters/STM/2010_real", pattern = ".tif", full.names = T, recursive = T)

stm.env.explanatory.var <- stack(stm.env.explanatory.var.list)
names(stm.env.explanatory.var) <- unlist(strsplit(stm.env.explanatory.var.list, "/|.tif"))[seq(4,80,4)]

stm.env.explanatory.var <- stm.env.explanatory.var[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
stm.env.explanatory.var<-aggregate(stm.env.explanatory.var,fact=35,fun=mean,na.rm=T); stm.env.explanatory.var<-stack(stm.env.explanatory.var)
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


#response variables
forestdep.spplist <- read.csv("data/species_summary.csv")

sppdata.final <- read.csv("data/presence_records.csv")

#pseudo-absence full dataset
#all cells with less than 20% of undegraded primary forest in landscape
UFPls.invcore <- env.explanatory.var[["UPFls"]]
UFPls.invcore[UFPls.invcore>.2]<-NA
UFPls.invcore[UFPls.invcore<=.2]<-1

pa.data <- as.data.frame(xyFromCell(UFPls.invcore, cell = which(UFPls.invcore[]==1), spatial = F))
names(pa.data) <- c("Longitude", "Latitude")

#### folders ####
dir.create("models.output/maps/PGM", recursive = T)
dir.create("models.output/maps/PGM/2020_real", recursive = T)
dir.create("models.output/maps/PGM/2020_avoiddeforest", recursive = T)
dir.create("models.output/maps/PGM/2020_avoiddegrad", recursive = T)
dir.create("models.output/maps/PGM/2020_avoidboth", recursive = T)
dir.create("models.output/maps/PGM/2020_restor_wo_avoid", recursive = T)
dir.create("models.output/maps/PGM/2020_restor_n_avoid", recursive = T)

dir.create("models.output/maps/STM", recursive = T)
dir.create("models.output/maps/STM/2020_real", recursive = T)
dir.create("models.output/maps/STM/2020_avoiddeforest", recursive = T)
dir.create("models.output/maps/STM/2020_avoiddegrad", recursive = T)
dir.create("models.output/maps/STM/2020_avoidboth", recursive = T)
dir.create("models.output/maps/STM/2020_restor_wo_avoid", recursive = T)
dir.create("models.output/maps/STM/2020_restor_n_avoid", recursive = T)

dir.create("models.output/evaluation", recursive = T)

#### species distribution modelling ####
#i <- as.character(forestdep.spplist$Binomial[1])
for (i in forestdep.spplist$Binomial) {
  
  if(forestdep.spplist[forestdep.spplist$Binomial==i, "Done"]==TRUE) next
  
  occur <- sppdata.final[sppdata.final$Binomial==i,]
  
  if(nrow(occur)<5) next
  #nrow(occur)<10
  
  myRespName <- as.character(i)
  myResp <- c(rep.int(1, times = as.numeric(nrow(occur))), rep.int(NA, times = as.numeric(nrow(occur)*30)))
  myRespXY <- rbind(occur[,c("Longitude", "Latitude")], sample_n(pa.data, as.numeric(nrow(occur)*30), replace = T))
  myExpl <- env.explanatory.var
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
  # model fit parameters
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
  
  # graphic of model scores
  pdf(paste("models.output/evaluation/", i ,"_models_scores.pdf", sep="" ))
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
  # ensemble model scores to species list
  #EMeval <- as.data.frame(get_evaluations(myBiomodEM_all))
  EMeval <- ev[ev$Model.name %in% chosen.models.from.ev$Model.name,]
  forestdep.spplist[forestdep.spplist$Binomial==i,"ROC"] <- paste(round(mean(EMeval[EMeval$Eval.metric=="ROC","Testing.data"]),1), "+/-", round(sd(EMeval[EMeval$Eval.metric=="ROC","Testing.data"]),2))
  forestdep.spplist[forestdep.spplist$Binomial==i,"TSS"] <- paste(round(mean(EMeval[EMeval$Eval.metric=="TSS","Testing.data"]),1), "+/-", round(sd(EMeval[EMeval$Eval.metric=="TSS","Testing.data"]),2))
  forestdep.spplist[forestdep.spplist$Binomial==i,"Sensitivity"] <- paste(round(mean(EMeval[EMeval$Eval.metric=="TSS","Sensitivity"]),2), "+/-", round(sd(EMeval[EMeval$Eval.metric=="TSS","Sensitivity"]),2))
  forestdep.spplist[forestdep.spplist$Binomial==i,"Specificity"] <- paste(round(mean(EMeval[EMeval$Eval.metric=="TSS","Specificity"]),2), "+/-", round(sd(EMeval[EMeval$Eval.metric=="TSS","Specificity"]),2))
  forestdep.spplist[forestdep.spplist$Binomial==i,"Cutoff"] <- paste(round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2), "+/-", round((sd(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2))
  
  write.csv(forestdep.spplist, "models.output/evaluation/species_summary_n.csv")
  
  if (any(occur$Region=="PGM")) {
    
    # 2020 real
    pgm.2020real.raster.list <- list.files("rasters/PGM/2020_real/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020real <- stack(pgm.2020real.raster.list)
    names(pgm.2020real) <- unlist(strsplit(pgm.2020real.raster.list, "/|.tif"))[seq(4,80,4)]
    rm("pgm.2020real.raster.list")
    
    pgm.2020real <- pgm.2020real[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
    pgm.2020real<-aggregate(pgm.2020real,fact=35,fun=mean,na.rm=T); pgm.2020real<-stack(pgm.2020real)
    ## checking
    #names(pgm.2020real)
    #plot(pgm.2020real[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020real <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for real scenario")
      
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
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "occur", "sel.var.df", "env.explanatory.var",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 avoid degradation
    pgm.2020avoiddegrad.raster.list <- list.files("rasters/PGM/2020_avoiddegrad/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020avoiddegrad <- stack(pgm.2020avoiddegrad.raster.list)
    names(pgm.2020avoiddegrad) <- unlist(strsplit(pgm.2020avoiddegrad.raster.list, "/|.tif"))[seq(4,80,4)]
    rm("pgm.2020avoiddegrad.raster.list")
    
    pgm.2020avoiddegrad <- pgm.2020avoiddegrad[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
    pgm.2020avoiddegrad<-aggregate(pgm.2020avoiddegrad,fact=35,fun=mean,na.rm=T); pgm.2020avoiddegrad<-stack(pgm.2020avoiddegrad)
    ## checking
    #names(pgm.2020avoiddegrad)
    #plot(pgm.2020avoiddegrad[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020avoiddegrad <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for avoid degradation scenario")
      
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
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "occur", "sel.var.df",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 avoid deforestation
    pgm.2020avoiddeforest.raster.list <- list.files("rasters/PGM/2020_avoiddeforest/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020avoiddeforest <- stack(pgm.2020avoiddeforest.raster.list)
    names(pgm.2020avoiddeforest) <- unlist(strsplit(pgm.2020avoiddeforest.raster.list, "/|.tif"))[seq(4,80,4)]
    rm("pgm.2020avoiddeforest.raster.list")
    
    pgm.2020avoiddeforest <- pgm.2020avoiddeforest[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
    pgm.2020avoiddeforest<-aggregate(pgm.2020avoiddeforest,fact=35,fun=mean,na.rm=T); pgm.2020avoiddeforest<-stack(pgm.2020avoiddeforest)
    ## checking
    #names(pgm.2020avoiddeforest)
    #plot(pgm.2020avoiddeforest[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020avoiddeforest <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for avoid deforestation scenario")
      
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
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "occur", "sel.var.df",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 avoid both
    pgm.2020avoidboth.raster.list <- list.files("rasters/PGM/2020_avoidboth/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020avoidboth <- stack(pgm.2020avoidboth.raster.list)
    names(pgm.2020avoidboth) <- unlist(strsplit(pgm.2020avoidboth.raster.list, "/|.tif"))[seq(4,80,4)]
    rm("pgm.2020avoidboth.raster.list")
    
    pgm.2020avoidboth <- pgm.2020avoidboth[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
    pgm.2020avoidboth<-aggregate(pgm.2020avoidboth,fact=35,fun=mean,na.rm=T); pgm.2020avoidboth<-stack(pgm.2020avoidboth)
    ## checking
    #names(pgm.2020avoidboth)
    #plot(pgm.2020avoidboth[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020avoidboth <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for avoid both scenario")
      
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
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "occur", "sel.var.df",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 restoration without avoid
    pgm.2020restor_wo_avoid.raster.list <- list.files("rasters/PGM/2020_restor_wo_avoid/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020restor_wo_avoid <- stack(pgm.2020restor_wo_avoid.raster.list)
    names(pgm.2020restor_wo_avoid) <- unlist(strsplit(pgm.2020restor_wo_avoid.raster.list, "/|.tif"))[seq(4,80,4)]
    rm("pgm.2020restor_wo_avoid.raster.list")
    
    pgm.2020restor_wo_avoid <- pgm.2020restor_wo_avoid[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
    pgm.2020restor_wo_avoid<-aggregate(pgm.2020restor_wo_avoid,fact=35,fun=mean,na.rm=T); pgm.2020restor_wo_avoid<-stack(pgm.2020restor_wo_avoid)
    ## checking
    #names(pgm.2020restor_wo_avoid)
    #plot(pgm.2020restor_wo_avoid[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020restor_wo_avoid <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for restoration without avoid scenario")
      
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
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "occur", "sel.var.df",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 restoration and avoid
    pgm.2020restor_n_avoid.raster.list <- list.files("rasters/PGM/2020_restor_n_avoid/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020restor_n_avoid <- stack(pgm.2020restor_n_avoid.raster.list)
    names(pgm.2020restor_n_avoid) <- unlist(strsplit(pgm.2020restor_n_avoid.raster.list, "/|.tif"))[seq(4,80,4)]
    rm("pgm.2020restor_n_avoid.raster.list")
    
    pgm.2020restor_n_avoid <- pgm.2020restor_n_avoid[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
    pgm.2020restor_n_avoid<-aggregate(pgm.2020restor_n_avoid,fact=35,fun=mean,na.rm=T); pgm.2020restor_n_avoid<-stack(pgm.2020restor_n_avoid)
    ## checking
    #names(pgm.2020restor_n_avoid)
    #plot(pgm.2020restor_n_avoid[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020restor_n_avoid <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "restoration and avoid scenario")
      
      BIOMOD_LoadModels(bm.out = myBiomodModelOut, full.name = m, as = "mod")
      
      temp_workdir = NULL
      if (length(grep("MAXENT.Phillips$", m)) == 1) {
        temp_workdir = mod@model_output_dir
      }
      
      indivdual_proj_2020restor_n_avoid <- predict(mod, pgm.2020restor_n_avoid, temp_workdir = temp_workdir)
      names(indivdual_proj_2020restor_n_avoid) <- m
      proj_2020restor_n_avoid <- addLayer(proj_2020restor_n_avoid, indivdual_proj_2020restor_n_avoid)
      
    }
    
    
    
    # building a consensus map by mean weight
    proj_2020restor_n_avoid.conbywm <- weighted.mean(proj_2020restor_n_avoid, EMeval[EMeval$Eval.metric=="TSS","Testing.data"], na.rm=T)
    
    writeRaster(proj_2020restor_n_avoid.conbywm, paste0("models.output/maps/PGM/2020_restor_n_avoid/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
    
    #cutoff.th <- round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2)
    #
    #proj_2020restor_n_avoid.conbywm.bin <- proj_2020restor_n_avoid.conbywm
    #proj_2020restor_n_avoid.conbywm.bin[proj_2020restor_n_avoid.conbywm.bin>=cutoff.th]<-1
    #proj_2020restor_n_avoid.conbywm.bin[proj_2020restor_n_avoid.conbywm.bin<cutoff.th]<-0
    #
    #writeRaster(proj_2020restor_n_avoid.conbywm.bin, paste0("models.output/maps/PGM/2020_restor_n_avoid/", i, "_merged_algo_merged_dataset_merged_run_TSSBin.tif"), format = "GTiff", overwrite = T)
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "occur", "sel.var.df",
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
    names(stm.2020real) <- unlist(strsplit(stm.2020real.raster.list, "/|.tif"))[seq(4,80,4)]
    rm("stm.2020real.raster.list")
    
    stm.2020real <- stm.2020real[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
    stm.2020real<-aggregate(stm.2020real,fact=35,fun=mean,na.rm=T); stm.2020real<-stack(stm.2020real)
    ## checking
    #names(stm.2020real)
    #plot(stm.2020real[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020real <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for real scenario")
      
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
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "occur", "sel.var.df",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 avoid degradation
    stm.2020avoiddegrad.raster.list <- list.files("rasters/STM/2020_avoiddegrad/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020avoiddegrad <- stack(stm.2020avoiddegrad.raster.list)
    names(stm.2020avoiddegrad) <- unlist(strsplit(stm.2020avoiddegrad.raster.list, "/|.tif"))[seq(4,80,4)]
    rm("stm.2020avoiddegrad.raster.list")
    
    stm.2020avoiddegrad <- stm.2020avoiddegrad[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
    stm.2020avoiddegrad<-aggregate(stm.2020avoiddegrad,fact=35,fun=mean,na.rm=T); stm.2020avoiddegrad<-stack(stm.2020avoiddegrad)
    ## checking
    #names(stm.2020avoiddegrad)
    #plot(stm.2020avoiddegrad[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020avoiddegrad <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for avoid degradation scenario")
      
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
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "occur", "sel.var.df",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 avoid deforestation
    stm.2020avoiddeforest.raster.list <- list.files("rasters/STM/2020_avoiddeforest/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020avoiddeforest <- stack(stm.2020avoiddeforest.raster.list)
    names(stm.2020avoiddeforest) <- unlist(strsplit(stm.2020avoiddeforest.raster.list, "/|.tif"))[seq(4,80,4)]
    rm("stm.2020avoiddeforest.raster.list")
    
    stm.2020avoiddeforest <- stm.2020avoiddeforest[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
    stm.2020avoiddeforest<-aggregate(stm.2020avoiddeforest,fact=35,fun=mean,na.rm=T); stm.2020avoiddeforest<-stack(stm.2020avoiddeforest)
    ## checking
    #names(stm.2020avoiddeforest)
    #plot(stm.2020avoiddeforest[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020avoiddeforest <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for avoid deforestation scenario")
      
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
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "occur", "sel.var.df",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 avoid both
    stm.2020avoidboth.raster.list <- list.files("rasters/STM/2020_avoidboth/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020avoidboth <- stack(stm.2020avoidboth.raster.list)
    names(stm.2020avoidboth) <- unlist(strsplit(stm.2020avoidboth.raster.list, "/|.tif"))[seq(4,80,4)]
    rm("stm.2020avoidboth.raster.list")
    
    stm.2020avoidboth <- stm.2020avoidboth[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
    stm.2020avoidboth<-aggregate(stm.2020avoidboth,fact=35,fun=mean,na.rm=T); stm.2020avoidboth<-stack(stm.2020avoidboth)
    ## checking
    #names(stm.2020avoidboth)
    #plot(stm.2020avoidboth[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020avoidboth <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for avoid both scenario")
      
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
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "occur", "sel.var.df",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 restoration without avoid
    stm.2020restor_wo_avoid.raster.list <- list.files("rasters/STM/2020_restor_wo_avoid/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020restor_wo_avoid <- stack(stm.2020restor_wo_avoid.raster.list)
    names(stm.2020restor_wo_avoid) <- unlist(strsplit(stm.2020restor_wo_avoid.raster.list, "/|.tif"))[seq(4,80,4)]
    rm("stm.2020restor_wo_avoid.raster.list")
    
    stm.2020restor_wo_avoid <- stm.2020restor_wo_avoid[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
    stm.2020restor_wo_avoid<-aggregate(stm.2020restor_wo_avoid,fact=35,fun=mean,na.rm=T); stm.2020restor_wo_avoid<-stack(stm.2020restor_wo_avoid)
    ## checking
    #names(stm.2020restor_wo_avoid)
    #plot(stm.2020restor_wo_avoid[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020restor_wo_avoid <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "for restoration without avoid scenario")
      
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
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "occur", "sel.var.df",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
    
    
    # 2020 restoration and avoid
    stm.2020restor_n_avoid.raster.list <- list.files("rasters/STM/2020_restor_n_avoid/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020restor_n_avoid <- stack(stm.2020restor_n_avoid.raster.list)
    names(stm.2020restor_n_avoid) <- unlist(strsplit(stm.2020restor_n_avoid.raster.list, "/|.tif"))[seq(4,80,4)]
    rm("stm.2020restor_n_avoid.raster.list")
    
    stm.2020restor_n_avoid <- stm.2020restor_n_avoid[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
    stm.2020restor_n_avoid<-aggregate(stm.2020restor_n_avoid,fact=35,fun=mean,na.rm=T); stm.2020restor_n_avoid<-stack(stm.2020restor_n_avoid)
    ## checking
    #names(stm.2020restor_n_avoid)
    #plot(stm.2020restor_n_avoid[[12]])
    #points(SpatialPoints(occur[,c("Longitude", "Latitude")]))
    
    proj_2020restor_n_avoid <- stack()
    for(m in chosen.models.from.all.models){
      
      cat("\n\t> Projecting", m, "restoration and avoid scenario")
      
      BIOMOD_LoadModels(bm.out = myBiomodModelOut, full.name = m, as = "mod")
      
      temp_workdir = NULL
      if (length(grep("MAXENT.Phillips$", m)) == 1) {
        temp_workdir = mod@model_output_dir
      }
      
      indivdual_proj_2020restor_n_avoid <- predict(mod, stm.2020restor_n_avoid, temp_workdir = temp_workdir)
      names(indivdual_proj_2020restor_n_avoid) <- m
      proj_2020restor_n_avoid <- addLayer(proj_2020restor_n_avoid, indivdual_proj_2020restor_n_avoid)
      
    }
    
    
    
    # building a consensus map by mean weight
    proj_2020restor_n_avoid.conbywm <- weighted.mean(proj_2020restor_n_avoid, EMeval[EMeval$Eval.metric=="TSS","Testing.data"], na.rm=T)
    
    writeRaster(proj_2020restor_n_avoid.conbywm, paste0("models.output/maps/STM/2020_restor_n_avoid/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
    
    #cutoff.th <- round((mean(EMeval[EMeval$Eval.metric=="TSS","Cutoff"]))/1000,2)
    #
    #proj_2020restor_n_avoid.conbywm.bin <- proj_2020restor_n_avoid.conbywm
    #proj_2020restor_n_avoid.conbywm.bin[proj_2020restor_n_avoid.conbywm.bin>=cutoff.th]<-1
    #proj_2020restor_n_avoid.conbywm.bin[proj_2020restor_n_avoid.conbywm.bin<cutoff.th]<-0
    #
    #writeRaster(proj_2020restor_n_avoid.conbywm.bin, paste0("models.output/maps/STM/2020_restor_n_avoid/", i, "_merged_algo_merged_dataset_merged_run_TSSBin.tif"), format = "GTiff", overwrite = T)
    
    
    rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "occur", "sel.var.df",
                                "myBiomodData", "myBiomodOption", "myBiomodModelOut", "ev", "eval_threshold", "EMeval",
                                "chosen.models.from.ev", "chosen.models.from.all.models"))])
    gc()
    #
    
  }  
  

#

  forestdep.spplist[forestdep.spplist$Binomial==i,"Done"] <- TRUE
  write.csv(forestdep.spplist, "models.output/evaluation/updated_species_summary_vn.csv")
  
  rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i"))])


  #do.call(file.remove, list(list.files(paste0("models.output/", i), full.names = T, recursive = T)))
  unlink(paste0("models.output/", i, "/*"), recursive = T, force = T)
  
}







