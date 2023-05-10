
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
names(pgm.env.explanatory.var) <- unlist(strsplit(pgm.env.explanatory.var.list, "/|.tif"))[seq(4,84,4)]

pgm.env.explanatory.var <- pgm.env.explanatory.var[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
##cheking
#pgm.env.explanatory.var
#plot(pgm.env.explanatory.var[[1:10]], nc=2)
#plot(pgm.env.explanatory.var[[11:20]], nc=2)



stm.env.explanatory.var.list <- list.files("rasters/STM/2010_real", pattern = ".tif", full.names = T, recursive = T)

stm.env.explanatory.var <- stack(stm.env.explanatory.var.list)
names(stm.env.explanatory.var) <- unlist(strsplit(stm.env.explanatory.var.list, "/|.tif"))[seq(4,84,4)]

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


#response variables
forestdep.spplist <- read.csv("data/species_summary.csv")
#forestdep.spplist <- read.csv("models.output/evaluation/updated_species_summary_v1.csv")

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
dir.create("models.output/maps/PGM/2020_restor_n_avoid_deforest", recursive = T)
dir.create("models.output/maps/PGM/2020_restor_n_avoid_both", recursive = T)

dir.create("models.output/maps/STM", recursive = T)
dir.create("models.output/maps/STM/2020_real", recursive = T)
dir.create("models.output/maps/STM/2020_avoiddeforest", recursive = T)
dir.create("models.output/maps/STM/2020_avoiddegrad", recursive = T)
dir.create("models.output/maps/STM/2020_avoidboth", recursive = T)
dir.create("models.output/maps/STM/2020_restor_wo_avoid", recursive = T)
dir.create("models.output/maps/STM/2020_restor_n_avoid_deforest", recursive = T)
dir.create("models.output/maps/STM/2020_restor_n_avoid_both", recursive = T)

dir.create("models.output/evaluation", recursive = T)








##############################################
####        biodiversity benefit          ####
##############################################
#### for species with more than 5 records ####
####    species distribution modelling    ####
##############################################

#i <- as.character(forestdep.spplist$Binomial[1])
for (i in forestdep.spplist$Binomial) {
  
  #spliting the job
  if(forestdep.spplist[forestdep.spplist$Binomial==i, "job"] != 2) next
  
  #checking if model is done
  if(forestdep.spplist[forestdep.spplist$Binomial==i, "Done"] == TRUE) next
  
  #starting modeling procedure
  occur <- sppdata.final[sppdata.final$Binomial==i,]
  
  if(nrow(occur)<5) next
  #nrow(occur)<10
  
  ##data preparation
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
  
  # graphic of model scores
  pdf(paste("models.output/evaluation/", i ,"_models_scores.pdf", sep="" ))
  gEval <- bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'run'))
  dev.off()
  
  ########################################
  ##       if the model is ready        ##
  ##   and need only project the maps   ##
  ##        start from here             ##
  ########################################
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
    names(pgm.2020real) <- unlist(strsplit(pgm.2020real.raster.list, "/|.tif"))[seq(4,84,4)]
    rm("pgm.2020real.raster.list")
    
    pgm.2020real <- pgm.2020real[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
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
    names(pgm.2020avoiddegrad) <- unlist(strsplit(pgm.2020avoiddegrad.raster.list, "/|.tif"))[seq(4,84,4)]
    rm("pgm.2020avoiddegrad.raster.list")
    
    pgm.2020avoiddegrad <- pgm.2020avoiddegrad[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
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
    names(pgm.2020avoiddeforest) <- unlist(strsplit(pgm.2020avoiddeforest.raster.list, "/|.tif"))[seq(4,84,4)]
    rm("pgm.2020avoiddeforest.raster.list")
    
    pgm.2020avoiddeforest <- pgm.2020avoiddeforest[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
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
    
    
    
    # 2020 avoid both
    pgm.2020avoidboth.raster.list <- list.files("rasters/PGM/2020_avoidboth/", pattern = ".tif", full.names = T, recursive = T)
    pgm.2020avoidboth <- stack(pgm.2020avoidboth.raster.list)
    names(pgm.2020avoidboth) <- unlist(strsplit(pgm.2020avoidboth.raster.list, "/|.tif"))[seq(4,84,4)]
    rm("pgm.2020avoidboth.raster.list")
    
    pgm.2020avoidboth <- pgm.2020avoidboth[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
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
    names(pgm.2020restor_wo_avoid) <- unlist(strsplit(pgm.2020restor_wo_avoid.raster.list, "/|.tif"))[seq(4,84,4)]
    rm("pgm.2020restor_wo_avoid.raster.list")
    
    pgm.2020restor_wo_avoid <- pgm.2020restor_wo_avoid[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
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
    names(pgm.2020restor_n_avoid_deforest) <- unlist(strsplit(pgm.2020restor_n_avoid_deforest.raster.list, "/|.tif"))[seq(4,84,4)]
    rm("pgm.2020restor_n_avoid_deforest.raster.list")
    
    pgm.2020restor_n_avoid_deforest <- pgm.2020restor_n_avoid_deforest[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
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
    names(pgm.2020restor_n_avoid_both) <- unlist(strsplit(pgm.2020restor_n_avoid_both.raster.list, "/|.tif"))[seq(4,84,4)]
    rm("pgm.2020restor_n_avoid_both.raster.list")
    
    pgm.2020restor_n_avoid_both <- pgm.2020restor_n_avoid_both[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
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
    names(stm.2020real) <- unlist(strsplit(stm.2020real.raster.list, "/|.tif"))[seq(4,84,4)]
    rm("stm.2020real.raster.list")
    
    stm.2020real <- stm.2020real[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
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
    names(stm.2020avoiddegrad) <- unlist(strsplit(stm.2020avoiddegrad.raster.list, "/|.tif"))[seq(4,84,4)]
    rm("stm.2020avoiddegrad.raster.list")
    
    stm.2020avoiddegrad <- stm.2020avoiddegrad[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
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
    names(stm.2020avoiddeforest) <- unlist(strsplit(stm.2020avoiddeforest.raster.list, "/|.tif"))[seq(4,84,4)]
    rm("stm.2020avoiddeforest.raster.list")
    
    stm.2020avoiddeforest <- stm.2020avoiddeforest[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
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
    
    
    
    # 2020 avoid both
    stm.2020avoidboth.raster.list <- list.files("rasters/STM/2020_avoidboth/", pattern = ".tif", full.names = T, recursive = T)
    stm.2020avoidboth <- stack(stm.2020avoidboth.raster.list)
    names(stm.2020avoidboth) <- unlist(strsplit(stm.2020avoidboth.raster.list, "/|.tif"))[seq(4,84,4)]
    rm("stm.2020avoidboth.raster.list")
    
    stm.2020avoidboth <- stm.2020avoidboth[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
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
    names(stm.2020restor_wo_avoid) <- unlist(strsplit(stm.2020restor_wo_avoid.raster.list, "/|.tif"))[seq(4,84,4)]
    rm("stm.2020restor_wo_avoid.raster.list")
    
    stm.2020restor_wo_avoid <- stm.2020restor_wo_avoid[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
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
    names(stm.2020restor_n_avoid_deforest) <- unlist(strsplit(stm.2020restor_n_avoid_deforest.raster.list, "/|.tif"))[seq(4,84,4)]
    rm("stm.2020restor_n_avoid_deforest.raster.list")
    
    stm.2020restor_n_avoid_deforest <- stm.2020restor_n_avoid_deforest[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
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
    names(stm.2020restor_n_avoid_both) <- unlist(strsplit(stm.2020restor_n_avoid_both.raster.list, "/|.tif"))[seq(4,84,4)]
    rm("stm.2020restor_n_avoid_both.raster.list")
    
    stm.2020restor_n_avoid_both <- stm.2020restor_n_avoid_both[[sel.var.df[!is.na(sel.var.df$VIF),"VAR"]]]
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






##############################################
####        biodiversity benefit          ####
##############################################
#### for species with less than 5 records ####
####              % of UPFls              ####
##############################################

forestdep.spplist <- read.csv("data/updated_species_summary_edby_visual_inspection.csv")

sppdata.final <- read.csv("data/presence_records.csv")


j=nrow(forestdep.spplist[forestdep.spplist$Done==F,])
for (i in forestdep.spplist$Binomial) {
  
  #checking if model is done
  if(forestdep.spplist[forestdep.spplist$Binomial==i, "Done"] != FALSE) next
  
  #starting modeling procedure
  occur <- sppdata.final[sppdata.final$Binomial==i,]
  
  #  
  
if (any(occur$Region=="PGM")) {
  
  pgm.occur <- occur[occur$Region=="PGM",]
  pgm.occur <- SpatialPoints(occur[,c("Longitude", "Latitude")])
  pres.bkg <- gBuffer(pgm.occur, width = 0.09)
  
  # 2020 real
  pgm.2020real <- raster("rasters/PGM/2020_real/UPFls.tif")
  pgm.pres.2020real <- mask(pgm.2020real, pres.bkg, updatevalue=0)
  
  writeRaster(pgm.pres.2020real, paste0("models.output/maps/PGM/2020_real/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
  
  
  # 2020 avoid degradation
  pgm.2020avoiddegrad <- raster("rasters/PGM/2020_avoiddegrad/UPFls.tif")
  pgm.2020avoiddegrad <- mask(pgm.2020avoiddegrad, pres.bkg, updatevalue=0)
  
  writeRaster(pgm.2020avoiddegrad, paste0("models.output/maps/PGM/2020_avoiddegrad/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
  
  
  # 2020 avoid deforestation
  pgm.2020avoiddeforest <- raster("rasters/PGM/2020_avoiddeforest/UPFls.tif")
  pgm.2020avoiddeforest <- mask(pgm.2020avoiddeforest, pres.bkg, updatevalue=0)
  
  writeRaster(pgm.2020avoiddeforest, paste0("models.output/maps/PGM/2020_avoiddeforest/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
  
  
  # 2020 avoid both
  pgm.2020avoidboth <- raster("rasters/PGM/2020_avoidboth/UPFls.tif")
  pgm.2020avoidboth <- mask(pgm.2020avoidboth, pres.bkg, updatevalue=0)
  
  writeRaster(pgm.2020avoidboth, paste0("models.output/maps/PGM/2020_avoidboth/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
  
  
  # 2020 restoration without avoid
  pgm.2020restor_wo_avoid <- raster("rasters/PGM/2020_restor_wo_avoid/UPFls.tif")
  pgm.2020restor_wo_avoid <- mask(pgm.2020restor_wo_avoid, pres.bkg, updatevalue=0)
  
  writeRaster(pgm.2020restor_wo_avoid, paste0("models.output/maps/PGM/2020_restor_wo_avoid/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
  
  
  # 2020 restoration and avoid restoration
  pgm.2020restor_n_avoid_deforest <- raster("rasters/PGM/2020_restor_n_avoid_deforest/UPFls.tif")
  pgm.2020restor_n_avoid_deforest <- mask(pgm.2020restor_n_avoid_deforest, pres.bkg, updatevalue=0)
  
  writeRaster(pgm.2020restor_n_avoid_deforest, paste0("models.output/maps/PGM/2020_restor_n_avoid_deforest/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
  
  
  # 2020 restoration and avoid both
  pgm.2020restor_n_avoid_both <- raster("rasters/PGM/2020_restor_n_avoid_both/UPFls.tif")
  pgm.2020restor_n_avoid_both <- mask(pgm.2020restor_n_avoid_both, pres.bkg, updatevalue=0)
  
  writeRaster(pgm.2020restor_n_avoid_both, paste0("models.output/maps/PGM/2020_restor_n_avoid_both/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
  
  
  rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "j", "occur"))])
  gc()
  #
  
}
  
  
  #  

if (any(occur$Region=="STM")) {
  
  stm.occur <- occur[occur$Region=="STM",]
  stm.occur <- SpatialPoints(occur[,c("Longitude", "Latitude")])
  pres.bkg <- gBuffer(stm.occur, width = 0.09)
  
  # 2020 real
  stm.2020real <- raster("rasters/STM/2020_real/UPFls.tif")
  stm.pres.2020real <- mask(stm.2020real, pres.bkg, updatevalue=0)
  
  writeRaster(stm.pres.2020real, paste0("models.output/maps/STM/2020_real/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
  
  
  # 2020 avoid degradation
  stm.2020avoiddegrad <- raster("rasters/STM/2020_avoiddegrad/UPFls.tif")
  stm.2020avoiddegrad <- mask(stm.2020avoiddegrad, pres.bkg, updatevalue=0)
  
  writeRaster(stm.2020avoiddegrad, paste0("models.output/maps/STM/2020_avoiddegrad/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
  
  
  # 2020 avoid deforestation
  stm.2020avoiddeforest <- raster("rasters/STM/2020_avoiddeforest/UPFls.tif")
  stm.2020avoiddeforest <- mask(stm.2020avoiddeforest, pres.bkg, updatevalue=0)
  
  writeRaster(stm.2020avoiddeforest, paste0("models.output/maps/STM/2020_avoiddeforest/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
  
  
  # 2020 avoid both
  stm.2020avoidboth <- raster("rasters/STM/2020_avoidboth/UPFls.tif")
  stm.2020avoidboth <- mask(stm.2020avoidboth, pres.bkg, updatevalue=0)
  
  writeRaster(stm.2020avoidboth, paste0("models.output/maps/STM/2020_avoidboth/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
  
  
  # 2020 restoration without avoid
  stm.2020restor_wo_avoid <- raster("rasters/STM/2020_restor_wo_avoid/UPFls.tif")
  stm.2020restor_wo_avoid <- mask(stm.2020restor_wo_avoid, pres.bkg, updatevalue=0)
  
  writeRaster(stm.2020restor_wo_avoid, paste0("models.output/maps/STM/2020_restor_wo_avoid/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
  
  
  # 2020 restoration and avoid restoration
  stm.2020restor_n_avoid_deforest <- raster("rasters/STM/2020_restor_n_avoid_deforest/UPFls.tif")
  stm.2020restor_n_avoid_deforest <- mask(stm.2020restor_n_avoid_deforest, pres.bkg, updatevalue=0)
  
  writeRaster(stm.2020restor_n_avoid_deforest, paste0("models.output/maps/STM/2020_restor_n_avoid_deforest/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
  
  
  # 2020 restoration and avoid both
  stm.2020restor_n_avoid_both <- raster("rasters/STM/2020_restor_n_avoid_both/UPFls.tif")
  stm.2020restor_n_avoid_both <- mask(stm.2020restor_n_avoid_both, pres.bkg, updatevalue=0)
  
  writeRaster(stm.2020restor_n_avoid_both, paste0("models.output/maps/STM/2020_restor_n_avoid_both/", i, "_merged_algo_merged_dataset_merged_run.tif"), format = "GTiff", overwrite = T)
  
  
  rm(list= ls()[!(ls() %in% c("forestdep.spplist", "sppdata.final", "i", "j", "occur"))])
  gc()
  #
  
}
  
  forestdep.spplist[forestdep.spplist$Binomial==i,"Done"] <- TRUE
  write.csv(forestdep.spplist, "models.output/evaluation/updated_species_summary_edby_visual_inspection_v2.csv.csv")
  
  
  j=j-1
  cat("\n> species", i, "done; now", j, "species to go <\n")
}


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
forestdep.spplist <- forestdep.spplist %>% left_join(bird.distribution.shp@data[,c("Binomial", "Shape_Area")]) %>%
                        group_by(Binomial) %>% mutate(Shape_Area = sum(Shape_Area)) %>% distinct(Binomial, .keep_all = T)

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

rm(treedata)


# saving
write.csv(forestdep.spplist, "data/species_summary_final.csv")
#forestdep.spplist <- read.csv("data/species_summary_final.csv")



#### biodiversity benefit ####
forestdep.spplist.total <- forestdep.spplist
forestdep.spplist <- forestdep.spplist %>% filter(Method != "Excluded" & !is.na(Shape_Area))

scenarios <- c("PGM/2020_real/", "PGM/2020_avoiddegrad/", "PGM/2020_avoiddeforest/", "PGM/2020_avoidboth/",
               "PGM/2020_restor_wo_avoid/", "PGM/2020_restor_n_avoid_deforest/", "PGM/2020_restor_n_avoid_both/",
               "STM/2020_real/", "STM/2020_avoiddegrad/", "STM/2020_avoiddeforest/", "STM/2020_avoidboth/",
               "STM/2020_restor_wo_avoid/", "STM/2020_restor_n_avoid_deforest/", "STM/2020_restor_n_avoid_both/")


for (s in scenarios) {
  
  maps.list <- list.files(paste0("models.output/maps/", s), pattern = ".tif", full.names = T, recursive = T)
  maps.list <- grep(paste(forestdep.spplist$Binomial, collapse = "|"), maps.list, value = T)
  
  biodiversity.maps <- stack(maps.list)
  names(biodiversity.maps) <- gsub("_merged_algo_merged_dataset_merged_run.tif", "", 
                                   gsub(paste0("models.output/maps/", s), "", maps.list))
  
  conservation.value <- forestdep.spplist[grep(paste(names(biodiversity.maps), collapse = "|"), forestdep.spplist$Binomial), "Shape_Area_scaled"]
  
  
  biodiversity.benefit.a <- stack()
  for (m in 1:length(conservation.value)) {
    biodiversity.benefit.x <- biodiversity.maps[[m]] * conservation.value[[m]]
    biodiversity.benefit.a <- addLayer(biodiversity.benefit.a, biodiversity.benefit.x)
    cat("\n> add conservation value to", names(biodiversity.maps[[m]]), "in scenario", s,  "<\n")
  }
  rm(biodiversity.benefit.x); gc()
  
  biodiversity.benefit <- sum(biodiversity.benefit.a, na.rm = T)
  
  writeRaster(biodiversity.benefit, paste0("models.output/biodiversity.benefits/", gsub("/", "_", s), "biodiversity_benefit.tif"), format="GTiff", overwrite=T)
  
  rm(list=ls()[!ls() %in% c("forestdep.spplist.total", "forestdep.spplist", "biodiversity.benefit", "s")]) #keeping only raster stack
  gc()
  cat("\n>>> SCENARIO", s, "DONE <<<\n")
  
  
}










##############################################
####            carbon benefit            ####
##############################################

## transect data
transectdata <- read.csv("data/input/RAS_transects_environment_all.csv")
#head(transectdata)
#str(transectdata)
#summary(transectdata)

#excluding Varzea transects PGM
exclude <- c("100_1", "100_4","100_7","81_12","423_2") #transect code
transectdata <- transectdata[!transectdata$Transectcode %in% exclude,]

#adding variable considering 50% of AGB as carbon stock
carbon <- transectdata %>% dplyr::select(Region:UTM_Y, LU_FT_Code, AGB_Trees10) %>% 
  mutate(carbon_stock = as.numeric(AGB_Trees10)/2)

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
#plot(env.explanatory.var[["UPFls"]])
#points(longlat)

#extracting values from environmental explanatory variables
env.var <- extract(env.explanatory.var, SpatialPoints(carbon[,c("Longitude", "Latitude")]))
carbon <- cbind(carbon, env.var)


#scaling predictors -- z-scores
carbon <- carbon %>% mutate(distriver_z = ((distriver - mean(distriver))/sd(distriver)),
                            distroad_z = ((distroad - mean(distroad))/sd(distroad)),
                            DPFpx_z = ((DPFpx - mean(DPFpx))/sd(DPFpx)),
                            edgedist_z = ((edgedist - mean(edgedist))/sd(edgedist)),
                            edgels_z = ((edgels - mean(edgels))/sd(edgels)),
                            elevation_z = ((elevation - mean(elevation))/sd(elevation)),
                            meanprecips_z = ((meanprecips - mean(meanprecips))/sd(meanprecips)),
                            meantemps_z = ((meantemps - mean(meantemps))/sd(meantemps)),
                            SFagels_z = ((SFagels - mean(SFagels))/sd(SFagels)),
                            SFpx_z = ((SFpx - mean(SFpx))/sd(SFpx)),
                            TFpx_z = ((TFpx - mean(TFpx))/sd(TFpx)),
                            UPFls_z = ((UPFls - mean(UPFls))/sd(UPFls)),
                            TSDls_z = ((TSDls - mean(TSDls))/sd(TSDls)))




carbon <- carbon[-which(is.na(carbon$carbon_stock)),]
write.csv(carbon, "data/carbon.csv")
#carbon <- read.csv("data/carbon.csv")

rm(list= ls()[!(ls() %in% c("sel.var.df", "carbon", "env.explanatory.var"))])
gc()


#loading libraries
library(fitdistrplus)
library(glmmTMB)
library(DHARMa)
library(lme4)
library(lmerTest)
library(mgcv)
library(randomForest)
library(MuMIn)
library(AICcmodavg)
library(lattice)
library(visreg)
library(mecofun)
library(dismo)
library(caret)




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




#null model
mod.null <- lm(carbon_stock ~ 1, data = subset(carbon, non_zero==1))
summary(mod.null)




#k-fold
set.seed(999)
group <- kfold(carbon, k=5, by = carbon$Catchment)
unique(group)

#model fiting and cross-validated predictions for GLM:
eval.mlr <- NULL
models.mlr.list <- list()

for (i in 1:5) {
  train <- carbon[group != i,c(9,10,26:38)]
  test <- carbon[group == i,c(9,10,26:38)]
  crosspred.mlr <- glm(carbon_stock ~ distriver_z + distroad_z + DPFpx_z + edgedist_z + edgels_z + elevation_z + meanprecips_z +
                         meantemps_z + SFagels_z + SFpx_z + TFpx_z + UPFls_z + TSDls_z, data = subset(train, non_zero == 1))
  model.sel <- step(crosspred.mlr, test="LRT")
  mod.mlr.fn <- glm(model.sel$formula, data = subset(train, non_zero == 1))
  models.mlr.list[[i]] <- mod.mlr.fn
  
  predict.mlr <- predict(mod.mlr.fn, subset(test, non_zero==1))
  #predict.mlr <- predict(crosspred.mlr, subset(test, non_zero==1))
  eval.mlr[i] <- caret::RMSE(predict.mlr, test[test$non_zero==1, "carbon_stock"])
  
}

model.sel(models.mlr.list)
eval.mlr

#best AICc [==2071.1]; not the best R2 [==39.75]
mod.mlr.fn1 <- glm(carbon_stock ~ SFpx_z + TFpx_z + TSDls_z, data = subset(carbon, non_zero==1))
summary(mod.mlr.fn1)
r.squaredLR(mod.mlr.fn1)
plot(mod.mlr.fn1)

#second best AICc [==2143.3]; best R2 [==46.64]
mod.mlr.fn2 <- glm(carbon_stock ~ edgedist_z + SFpx_z + TFpx_z + UPFls_z, data = subset(carbon, non_zero==1))
summary(mod.mlr.fn2)
r.squaredLR(mod.mlr.fn2)
plot(mod.mlr.fn2)

#adding TSDls_z to mod.mlr.fn2
mod.mlr.fn3 <- glm(carbon_stock ~ edgedist_z + SFpx_z + TFpx_z + UPFls_z + TSDls_z, data = subset(carbon, non_zero==1))
summary(mod.mlr.fn3)
r.squaredLR(mod.mlr.fn3)
plot(mod.mlr.fn3)

#adding TSDls_z to mod.mlr.fn2
mod.mlr.fn4 <- glm(carbon_stock ~ edgedist + SFpx + TFpx + UPFls + TSDls, data = subset(carbon, non_zero==1))
summary(mod.mlr.fn4)
r.squaredLR(mod.mlr.fn4)
plot(mod.mlr.fn4)


#visreg(mod.mlr.fn4, xvar = "edgedist_z", data = subset(carbon, non_zero==1))




#model fiting and cross-validated predictions for GLMM -- catchments as random effects
eval.glmm <- NULL
models.glmm.list <- list()

for (i in 1:5) {
  train <- carbon[group != i,c(2,9,10,26:38)]
  test <- carbon[group == i,c(2,9,10,26:38)]
  crosspred.glmm <- lmer(carbon_stock ~ distriver_z + distroad_z + DPFpx_z + edgedist_z + edgels_z + elevation_z + meanprecips_z +
                           meantemps_z + SFagels_z + SFpx_z + TFpx_z + UPFls_z + TSDls_z + (1|Catchment),
                         data = subset(train, non_zero == 1))
  step_res <- step(crosspred.glmm)
  models.glmm.list[[i]] <- get_model(step_res)
  
  predict.glmm <- predict(get_model(step_res), subset(test, non_zero==1))
  eval.glmm[i] <- caret::RMSE(predict.glmm, test[test$non_zero==1, "carbon_stock"])
  
}

model.sel(models.glmm.list)
eval.glmm
#obs. results were the same as glm; there is no influence of catchment as random variable




#model fiting for GAM
mod.gam.full <- gam(carbon_stock ~ s(distriver_z, bs = "cs") + s(distroad_z, bs = "cs") + s(DPFpx_z, bs = "cs") +
                      s(edgedist_z, bs = "cs") + s(edgels_z, bs = "cs") + s(elevation_z, bs = "cs") + s(meanprecips_z, bs = "cs") +
                      s(meantemps_z, bs = "cs") + s(SFagels_z, bs = "cs") + s(SFpx_z, bs = "cs") + s(TFpx_z, bs = "cs") +
                      s(UPFls_z, bs = "cs") + s(TSDls_z, bs = "cs"), method="REML", select = T,
                    data = subset(carbon, non_zero==1))
summary(mod.gam.full)
plot(mod.gam.full, pages=1, residuals=TRUE)
gam.check(mod.gam.full)

rsd <- residuals(mod.gam.full)
gam(rsd ~ s(UPFls_z, k=100, bs="cs"), gamma=1.4, data = subset(carbon, non_zero==1))


mod.gam1 <- gam(carbon_stock ~ s(edgedist_z, k = 6, bs = "cs") + s(SFpx_z, k = 6, bs = "cs") + 
                  s(TFpx_z, k = 6, bs = "cs") + s(UPFls_z, k = 6, bs = "cs") + s(TSDls_z, k = 6, bs = "cs"),
                method="REML", select = T, data = subset(carbon, non_zero==1))
summary(mod.gam1)
plot(mod.gam1, pages=1, residuals=TRUE)
gam.check(mod.gam1)


mod.gam2 <- gam(carbon_stock ~ s(edgedist_z, k = 6, fx = T, bs = "tp") + s(SFpx_z, k = 6, fx = T, bs = "tp") + 
                  s(TFpx_z, k = 6, fx = T, bs = "tp") + s(UPFls_z, k = 6, fx = T, bs = "tp") + s(TSDls_z, k = 6, fx = T, bs = "tp"),
                method="REML", select = T, data = subset(carbon, non_zero==1))
summary(mod.gam2)
plot(mod.gam2, pages=1, residuals=TRUE)
gam.check(mod.gam2)


mod.gam3 <- gam(carbon_stock ~ edgedist_z + SFpx_z +  s(TFpx_z, k = 6, bs = "cs") + s(UPFls_z, k = 7, bs = "cs") +
                  s(TSDls_z, k = 6, bs = "cs"), method="REML", select = T, data = subset(carbon, non_zero==1))
summary(mod.gam3)
plot(mod.gam3, pages=1, residuals=TRUE)
gam.check(mod.gam3)


mod.gam4 <- gam(carbon_stock ~ edgedist + SFpx +  s(TFpx, k = 6, bs = "cs") + s(UPFls, k = 7, bs = "cs") +
                  s(TSDls, k = 6, bs = "cs"), method="REML", select = T, data = subset(carbon, non_zero==1))
summary(mod.gam4)
plot(mod.gam4, pages=1, residuals=TRUE)
gam.check(mod.gam4)


model.sel(mod.gam1, mod.gam2, mod.gam3, mod.gam4, mod.mlr.fn2, mod.mlr.fn3, mod.mlr.fn4)

#obs. results were no batter than glm; the smooth parameter didn't result in better fit




#model fiting and cross-validated predictions for RF
eval.rf <- NULL
models.rf.list <- list()

for (i in 1:5) {
  train <- carbon[group != i,c(9,10,26:38)]
  test <- carbon[group == i,c(9,10,26:38)]
  models.rf.list[[i]] <- randomForest(x = train[train$non_zero==1, 3:15], y = train[train$non_zero==1,"carbon_stock"],
                                      xtest = test[test$non_zero==1, 3:15], ytest = test[test$non_zero==1,"carbon_stock"],
                                      ntree = 1000, nodesize = 10, importance =T, nPerm = 5)
  
}

par(mfrow=c(3,2))
plot(models.rf.list[[1]])
plot(models.rf.list[[2]])
plot(models.rf.list[[3]])
plot(models.rf.list[[4]])
plot(models.rf.list[[5]])

par(mfrow=c(3,2))
varImpPlot(models.rf.list[[1]], type = 1)
varImpPlot(models.rf.list[[2]], type = 1)
varImpPlot(models.rf.list[[3]], type = 1)
varImpPlot(models.rf.list[[4]], type = 1)
varImpPlot(models.rf.list[[5]], type = 1)

mod.rf.fn <- randomForest(x = carbon[carbon$non_zero==1, c(29,35:38)], y = carbon[carbon$non_zero==1,"carbon_stock"],
                          ntree=600, nodesize=10, importance =T, nPerm = 5)

plot(mod.rf.fn)
varImpPlot(mod.rf.fn, type = 1)

mod.rf.fn2 <- randomForest(x = carbon[carbon$non_zero==1, c(16,22:25)], y = carbon[carbon$non_zero==1,"carbon_stock"],
                           ntree=600, nodesize=10, importance =T, nPerm = 5)

plot(mod.rf.fn2)
varImpPlot(mod.rf.fn2, type = 1)


rm(list= ls()[!(ls() %in% c("sel.var.df", "carbon", "env.explanatory.var", "mod.mlr.fn4", "mod.rf.fn2"))])
gc()





#predctions
dir.create("models.output/carbon.benefits", recursive = T)
sel.var.df <- c("edgedist", "SFpx", "TFpx", "UPFls", "TSDls")

#PGM 2020 real
pgm.2020real.raster.list <- list.files("rasters/PGM/2020_real/", pattern = ".tif", full.names = T, recursive = T)
pgm.2020real <- stack(pgm.2020real.raster.list)
names(pgm.2020real) <- unlist(strsplit(pgm.2020real.raster.list, "/|.tif"))[seq(4,84,4)]
pgm.2020real <- pgm.2020real[[sel.var.df]]

set.seed(999)
mod.mlr.proj_pgm.2020real <- predict(pgm.2020real, mod.mlr.fn4)
set.seed(999)
mod.rf.proj_pgm.2020real <- predict(pgm.2020real, mod.rf.fn2)
  
  
#building a consensus map by mean weight
proj_pgm.2020real.conbywm <- weighted.mean(stack(mod.mlr.proj_pgm.2020real, mod.rf.proj_pgm.2020real), 
                                       c(RMSE(mod.mlr.fn4$fitted.values, mod.mlr.fn4$y),
                                         RMSE(mod.rf.fn2$predicted, mod.rf.fn2$y)), 
                                       na.rm=T)
  
writeRaster(proj_pgm.2020real.conbywm, paste0("models.output/carbon.benefits/PGM_2020_real_carbon_benefit.tif"), format = "GTiff", overwrite = T)
  
rm(list= ls()[!(ls() %in% c("sel.var.df", "carbon", "env.explanatory.var", "mod.mlr.fn4", "mod.rf.fn2",
                            "proj_pgm.2020real.conbywm"))])
gc()
#



#PGM 2020 avoid degradation
pgm.2020avoiddegrad.raster.list <- list.files("rasters/PGM/2020_avoiddegrad/", pattern = ".tif", full.names = T, recursive = T)
pgm.2020avoiddegrad <- stack(pgm.2020avoiddegrad.raster.list)
names(pgm.2020avoiddegrad) <- unlist(strsplit(pgm.2020avoiddegrad.raster.list, "/|.tif"))[seq(4,84,4)]
pgm.2020avoiddegrad <- pgm.2020avoiddegrad[[sel.var.df]]

set.seed(999)
mod.mlr.proj_pgm.2020avoiddegrad <- predict(pgm.2020avoiddegrad, mod.mlr.fn4)
set.seed(999)
mod.rf.proj_pgm.2020avoiddegrad <- predict(pgm.2020avoiddegrad, mod.rf.fn2)


#building a consensus map by mean weight
proj_pgm.2020avoiddegrad.conbywm <- weighted.mean(stack(mod.mlr.proj_pgm.2020avoiddegrad, mod.rf.proj_pgm.2020avoiddegrad), 
                                           c(RMSE(mod.mlr.fn4$fitted.values, mod.mlr.fn4$y),
                                             RMSE(mod.rf.fn2$predicted, mod.rf.fn2$y)), 
                                           na.rm=T)

writeRaster(proj_pgm.2020avoiddegrad.conbywm, paste0("models.output/carbon.benefits/PGM_2020_avoiddegrad_carbon_benefit.tif"), format = "GTiff", overwrite = T)

rm(list= ls()[!(ls() %in% c("sel.var.df", "carbon", "env.explanatory.var", "mod.mlr.fn4", "mod.rf.fn2",
                            "proj_pgm.2020real.conbywm", "proj_pgm.2020avoiddegrad.conbywm"))])
gc()
#



#PGM 2020 avoid deforestation
pgm.2020avoiddeforest.raster.list <- list.files("rasters/PGM/2020_avoiddeforest/", pattern = ".tif", full.names = T, recursive = T)
pgm.2020avoiddeforest <- stack(pgm.2020avoiddeforest.raster.list)
names(pgm.2020avoiddeforest) <- unlist(strsplit(pgm.2020avoiddeforest.raster.list, "/|.tif"))[seq(4,84,4)]
pgm.2020avoiddeforest <- pgm.2020avoiddeforest[[sel.var.df]]

set.seed(999)
mod.mlr.proj_pgm.2020avoiddeforest <- predict(pgm.2020avoiddeforest, mod.mlr.fn4)
set.seed(999)
mod.rf.proj_pgm.2020avoiddeforest <- predict(pgm.2020avoiddeforest, mod.rf.fn2)


#building a consensus map by mean weight
proj_pgm.2020avoiddeforest.conbywm <- weighted.mean(stack(mod.mlr.proj_pgm.2020avoiddeforest, mod.rf.proj_pgm.2020avoiddeforest), 
                                                  c(RMSE(mod.mlr.fn4$fitted.values, mod.mlr.fn4$y),
                                                    RMSE(mod.rf.fn2$predicted, mod.rf.fn2$y)), 
                                                  na.rm=T)

writeRaster(proj_pgm.2020avoiddeforest.conbywm, paste0("models.output/carbon.benefits/PGM_2020_avoiddeforest_carbon_benefit.tif"), format = "GTiff", overwrite = T)

rm(list= ls()[!(ls() %in% c("sel.var.df", "carbon", "env.explanatory.var", "mod.mlr.fn4", "mod.rf.fn2",
                            "proj_pgm.2020real.conbywm", "proj_pgm.2020avoiddegrad.conbywm", "proj_pgm.2020avoiddeforest.conbywm"))])
gc()
#



#PGM 2020 avoid both
pgm.2020avoidboth.raster.list <- list.files("rasters/PGM/2020_avoidboth/", pattern = ".tif", full.names = T, recursive = T)
pgm.2020avoidboth <- stack(pgm.2020avoidboth.raster.list)
names(pgm.2020avoidboth) <- unlist(strsplit(pgm.2020avoidboth.raster.list, "/|.tif"))[seq(4,84,4)]
pgm.2020avoidboth <- pgm.2020avoidboth[[sel.var.df]]

set.seed(999)
mod.mlr.proj_pgm.2020avoidboth <- predict(pgm.2020avoidboth, mod.mlr.fn4)
set.seed(999)
mod.rf.proj_pgm.2020avoidboth <- predict(pgm.2020avoidboth, mod.rf.fn2)


#building a consensus map by mean weight
proj_pgm.2020avoidboth.conbywm <- weighted.mean(stack(mod.mlr.proj_pgm.2020avoidboth, mod.rf.proj_pgm.2020avoidboth), 
                                                    c(RMSE(mod.mlr.fn4$fitted.values, mod.mlr.fn4$y),
                                                      RMSE(mod.rf.fn2$predicted, mod.rf.fn2$y)), 
                                                    na.rm=T)

writeRaster(proj_pgm.2020avoidboth.conbywm, paste0("models.output/carbon.benefits/PGM_2020_avoidboth_carbon_benefit.tif"), format = "GTiff", overwrite = T)

rm(list= ls()[!(ls() %in% c("sel.var.df", "carbon", "env.explanatory.var", "mod.mlr.fn4", "mod.rf.fn2",
                            "proj_pgm.2020real.conbywm", "proj_pgm.2020avoiddegrad.conbywm", "proj_pgm.2020avoiddeforest.conbywm",
                            "proj_pgm.2020avoidboth.conbywm"))])
gc()
#



#PGM 2020 restoration without avoid
pgm.2020restor_wo_avoid.raster.list <- list.files("rasters/PGM/2020_restor_wo_avoid/", pattern = ".tif", full.names = T, recursive = T)
pgm.2020restor_wo_avoid <- stack(pgm.2020restor_wo_avoid.raster.list)
names(pgm.2020restor_wo_avoid) <- unlist(strsplit(pgm.2020restor_wo_avoid.raster.list, "/|.tif"))[seq(4,84,4)]
pgm.2020restor_wo_avoid <- pgm.2020restor_wo_avoid[[sel.var.df]]

set.seed(999)
mod.mlr.proj_pgm.2020restor_wo_avoid <- predict(pgm.2020restor_wo_avoid, mod.mlr.fn4)
set.seed(999)
mod.rf.proj_pgm.2020restor_wo_avoid <- predict(pgm.2020restor_wo_avoid, mod.rf.fn2)


#building a consensus map by mean weight
proj_pgm.2020restor_wo_avoid.conbywm <- weighted.mean(stack(mod.mlr.proj_pgm.2020restor_wo_avoid, mod.rf.proj_pgm.2020restor_wo_avoid), 
                                                c(RMSE(mod.mlr.fn4$fitted.values, mod.mlr.fn4$y),
                                                  RMSE(mod.rf.fn2$predicted, mod.rf.fn2$y)), 
                                                na.rm=T)

writeRaster(proj_pgm.2020restor_wo_avoid.conbywm, paste0("models.output/carbon.benefits/PGM_2020_restor_wo_avoid_carbon_benefit.tif"), format = "GTiff", overwrite = T)

rm(list= ls()[!(ls() %in% c("sel.var.df", "carbon", "env.explanatory.var", "mod.mlr.fn4", "mod.rf.fn2",
                            "proj_pgm.2020real.conbywm", "proj_pgm.2020avoiddegrad.conbywm", "proj_pgm.2020avoiddeforest.conbywm",
                            "proj_pgm.2020avoidboth.conbywm", "proj_pgm.2020restor_wo_avoid.conbywm"))])
gc()
#



#PGM 2020 restoration and avoid deforest
pgm.2020restor_n_avoiddeforest.raster.list <- list.files("rasters/PGM/2020_restor_n_avoid_deforest/", pattern = ".tif", full.names = T, recursive = T)
pgm.2020restor_n_avoiddeforest <- stack(pgm.2020restor_n_avoiddeforest.raster.list)
names(pgm.2020restor_n_avoiddeforest) <- unlist(strsplit(pgm.2020restor_n_avoiddeforest.raster.list, "/|.tif"))[seq(4,84,4)]
pgm.2020restor_n_avoiddeforest <- pgm.2020restor_n_avoiddeforest[[sel.var.df]]

set.seed(999)
mod.mlr.proj_pgm.2020restor_n_avoiddeforest <- predict(pgm.2020restor_n_avoiddeforest, mod.mlr.fn4)
set.seed(999)
mod.rf.proj_pgm.2020restor_n_avoiddeforest <- predict(pgm.2020restor_n_avoiddeforest, mod.rf.fn2)


#building a consensus map by mean weight
proj_pgm.2020restor_n_avoiddeforest.conbywm <- weighted.mean(stack(mod.mlr.proj_pgm.2020restor_n_avoiddeforest, mod.rf.proj_pgm.2020restor_n_avoiddeforest), 
                                                      c(RMSE(mod.mlr.fn4$fitted.values, mod.mlr.fn4$y),
                                                        RMSE(mod.rf.fn2$predicted, mod.rf.fn2$y)), 
                                                      na.rm=T)

writeRaster(proj_pgm.2020restor_n_avoiddeforest.conbywm, paste0("models.output/carbon.benefits/PGM_2020_restor_n_avoid_deforest_carbon_benefit.tif"), format = "GTiff", overwrite = T)

rm(list= ls()[!(ls() %in% c("sel.var.df", "carbon", "env.explanatory.var", "mod.mlr.fn4", "mod.rf.fn2",
                            "proj_pgm.2020real.conbywm", "proj_pgm.2020avoiddegrad.conbywm", "proj_pgm.2020avoiddeforest.conbywm",
                            "proj_pgm.2020avoidboth.conbywm", "proj_pgm.2020restor_wo_avoid.conbywm", "proj_pgm.2020restor_n_avoiddeforest.conbywm"))])
gc()
#



#PGM 2020 restoration and avoid both
pgm.2020restor_n_avoidboth.raster.list <- list.files("rasters/PGM/2020_restor_n_avoid_both/", pattern = ".tif", full.names = T, recursive = T)
pgm.2020restor_n_avoidboth <- stack(pgm.2020restor_n_avoidboth.raster.list)
names(pgm.2020restor_n_avoidboth) <- unlist(strsplit(pgm.2020restor_n_avoidboth.raster.list, "/|.tif"))[seq(4,84,4)]
pgm.2020restor_n_avoidboth <- pgm.2020restor_n_avoidboth[[sel.var.df]]

set.seed(999)
mod.mlr.proj_pgm.2020restor_n_avoidboth <- predict(pgm.2020restor_n_avoidboth, mod.mlr.fn4)
set.seed(999)
mod.rf.proj_pgm.2020restor_n_avoidboth <- predict(pgm.2020restor_n_avoidboth, mod.rf.fn2)


#building a consensus map by mean weight
proj_pgm.2020restor_n_avoidboth.conbywm <- weighted.mean(stack(mod.mlr.proj_pgm.2020restor_n_avoidboth, mod.rf.proj_pgm.2020restor_n_avoidboth), 
                                                             c(RMSE(mod.mlr.fn4$fitted.values, mod.mlr.fn4$y),
                                                               RMSE(mod.rf.fn2$predicted, mod.rf.fn2$y)), 
                                                             na.rm=T)

writeRaster(proj_pgm.2020restor_n_avoidboth.conbywm, paste0("models.output/carbon.benefits/PGM_2020_restor_n_avoid_both_carbon_benefit.tif"), format = "GTiff", overwrite = T)

rm(list= ls()[!(ls() %in% c("sel.var.df", "carbon", "env.explanatory.var", "mod.mlr.fn4", "mod.rf.fn2",
                            "proj_pgm.2020real.conbywm", "proj_pgm.2020avoiddegrad.conbywm", "proj_pgm.2020avoiddeforest.conbywm",
                            "proj_pgm.2020avoidboth.conbywm", "proj_pgm.2020restor_wo_avoid.conbywm", "proj_pgm.2020restor_n_avoiddeforest.conbywm",
                            "proj_pgm.2020restor_n_avoidboth.conbywm"))])
gc()
#



#STM 2020 real
stm.2020real.raster.list <- list.files("rasters/STM/2020_real/", pattern = ".tif", full.names = T, recursive = T)
stm.2020real <- stack(stm.2020real.raster.list)
names(stm.2020real) <- unlist(strsplit(stm.2020real.raster.list, "/|.tif"))[seq(4,84,4)]
stm.2020real <- stm.2020real[[sel.var.df]]

set.seed(999)
mod.mlr.proj_stm.2020real <- predict(stm.2020real, mod.mlr.fn4)
set.seed(999)
mod.rf.proj_stm.2020real <- predict(stm.2020real, mod.rf.fn2)


#building a consensus map by mean weight
proj_stm.2020real.conbywm <- weighted.mean(stack(mod.mlr.proj_stm.2020real, mod.rf.proj_stm.2020real), 
                                           c(RMSE(mod.mlr.fn4$fitted.values, mod.mlr.fn4$y),
                                             RMSE(mod.rf.fn2$predicted, mod.rf.fn2$y)), 
                                           na.rm=T)

writeRaster(proj_stm.2020real.conbywm, paste0("models.output/carbon.benefits/STM_2020_real_carbon_benefit.tif"), format = "GTiff", overwrite = T)

rm(list= ls()[!(ls() %in% c("sel.var.df", "carbon", "env.explanatory.var", "mod.mlr.fn4", "mod.rf.fn2",
                            "proj_pgm.2020real.conbywm", "proj_pgm.2020avoiddegrad.conbywm", "proj_pgm.2020avoiddeforest.conbywm",
                            "proj_pgm.2020avoidboth.conbywm", "proj_pgm.2020restor_wo_avoid.conbywm", "proj_pgm.2020restor_n_avoiddeforest.conbywm",
                            "proj_pgm.2020restor_n_avoidboth.conbywm", "proj_stm.2020real.conbywm"))])
gc()
#



#STM 2020 avoid degradation
stm.2020avoiddegrad.raster.list <- list.files("rasters/STM/2020_avoiddegrad/", pattern = ".tif", full.names = T, recursive = T)
stm.2020avoiddegrad <- stack(stm.2020avoiddegrad.raster.list)
names(stm.2020avoiddegrad) <- unlist(strsplit(stm.2020avoiddegrad.raster.list, "/|.tif"))[seq(4,84,4)]
stm.2020avoiddegrad <- stm.2020avoiddegrad[[sel.var.df]]

set.seed(999)
mod.mlr.proj_stm.2020avoiddegrad <- predict(stm.2020avoiddegrad, mod.mlr.fn4)
set.seed(999)
mod.rf.proj_stm.2020avoiddegrad <- predict(stm.2020avoiddegrad, mod.rf.fn2)


#building a consensus map by mean weight
proj_stm.2020avoiddegrad.conbywm <- weighted.mean(stack(mod.mlr.proj_stm.2020avoiddegrad, mod.rf.proj_stm.2020avoiddegrad), 
                                                  c(RMSE(mod.mlr.fn4$fitted.values, mod.mlr.fn4$y),
                                                    RMSE(mod.rf.fn2$predicted, mod.rf.fn2$y)), 
                                                  na.rm=T)

writeRaster(proj_stm.2020avoiddegrad.conbywm, paste0("models.output/carbon.benefits/STM_2020_avoiddegrad_carbon_benefit.tif"), format = "GTiff", overwrite = T)

rm(list= ls()[!(ls() %in% c("sel.var.df", "carbon", "env.explanatory.var", "mod.mlr.fn4", "mod.rf.fn2",
                            "proj_pgm.2020real.conbywm", "proj_pgm.2020avoiddegrad.conbywm", "proj_pgm.2020avoiddeforest.conbywm",
                            "proj_pgm.2020avoidboth.conbywm", "proj_pgm.2020restor_wo_avoid.conbywm", "proj_pgm.2020restor_n_avoiddeforest.conbywm",
                            "proj_pgm.2020restor_n_avoidboth.conbywm", "proj_stm.2020real.conbywm", "proj_stm.2020avoiddegrad.conbywm"))])
gc()
#



#STM 2020 avoid deforestation
stm.2020avoiddeforest.raster.list <- list.files("rasters/STM/2020_avoiddeforest/", pattern = ".tif", full.names = T, recursive = T)
stm.2020avoiddeforest <- stack(stm.2020avoiddeforest.raster.list)
names(stm.2020avoiddeforest) <- unlist(strsplit(stm.2020avoiddeforest.raster.list, "/|.tif"))[seq(4,84,4)]
stm.2020avoiddeforest <- stm.2020avoiddeforest[[sel.var.df]]

set.seed(999)
mod.mlr.proj_stm.2020avoiddeforest <- predict(stm.2020avoiddeforest, mod.mlr.fn4)
set.seed(999)
mod.rf.proj_stm.2020avoiddeforest <- predict(stm.2020avoiddeforest, mod.rf.fn2)


#building a consensus map by mean weight
proj_stm.2020avoiddeforest.conbywm <- weighted.mean(stack(mod.mlr.proj_stm.2020avoiddeforest, mod.rf.proj_stm.2020avoiddeforest), 
                                                    c(RMSE(mod.mlr.fn4$fitted.values, mod.mlr.fn4$y),
                                                      RMSE(mod.rf.fn2$predicted, mod.rf.fn2$y)), 
                                                    na.rm=T)

writeRaster(proj_stm.2020avoiddeforest.conbywm, paste0("models.output/carbon.benefits/STM_2020_avoiddeforest_carbon_benefit.tif"), format = "GTiff", overwrite = T)

rm(list= ls()[!(ls() %in% c("sel.var.df", "carbon", "env.explanatory.var", "mod.mlr.fn4", "mod.rf.fn2",
                            "proj_pgm.2020real.conbywm", "proj_pgm.2020avoiddegrad.conbywm", "proj_pgm.2020avoiddeforest.conbywm",
                            "proj_pgm.2020avoidboth.conbywm", "proj_pgm.2020restor_wo_avoid.conbywm", "proj_pgm.2020restor_n_avoiddeforest.conbywm",
                            "proj_pgm.2020restor_n_avoidboth.conbywm", "proj_stm.2020real.conbywm", "proj_stm.2020avoiddegrad.conbywm",
                            "proj_stm.2020avoiddeforest.conbywm"))])
gc()
#



#STM 2020 avoid both
stm.2020avoidboth.raster.list <- list.files("rasters/STM/2020_avoidboth/", pattern = ".tif", full.names = T, recursive = T)
stm.2020avoidboth <- stack(stm.2020avoidboth.raster.list)
names(stm.2020avoidboth) <- unlist(strsplit(stm.2020avoidboth.raster.list, "/|.tif"))[seq(4,84,4)]
stm.2020avoidboth <- stm.2020avoidboth[[sel.var.df]]

set.seed(999)
mod.mlr.proj_stm.2020avoidboth <- predict(stm.2020avoidboth, mod.mlr.fn4)
set.seed(999)
mod.rf.proj_stm.2020avoidboth <- predict(stm.2020avoidboth, mod.rf.fn2)


#building a consensus map by mean weight
proj_stm.2020avoidboth.conbywm <- weighted.mean(stack(mod.mlr.proj_stm.2020avoidboth, mod.rf.proj_stm.2020avoidboth), 
                                                c(RMSE(mod.mlr.fn4$fitted.values, mod.mlr.fn4$y),
                                                  RMSE(mod.rf.fn2$predicted, mod.rf.fn2$y)), 
                                                na.rm=T)

writeRaster(proj_stm.2020avoidboth.conbywm, paste0("models.output/carbon.benefits/STM_2020_avoidboth_carbon_benefit.tif"), format = "GTiff", overwrite = T)

rm(list= ls()[!(ls() %in% c("sel.var.df", "carbon", "env.explanatory.var", "mod.mlr.fn4", "mod.rf.fn2",
                            "proj_pgm.2020real.conbywm", "proj_pgm.2020avoiddegrad.conbywm", "proj_pgm.2020avoiddeforest.conbywm",
                            "proj_pgm.2020avoidboth.conbywm", "proj_pgm.2020restor_wo_avoid.conbywm", "proj_pgm.2020restor_n_avoiddeforest.conbywm",
                            "proj_pgm.2020restor_n_avoidboth.conbywm", "proj_stm.2020real.conbywm", "proj_stm.2020avoiddegrad.conbywm",
                            "proj_stm.2020avoiddeforest.conbywm", "proj_stm.2020avoidboth.conbywm"))])
gc()
#



#STM 2020 restoration without avoid
stm.2020restor_wo_avoid.raster.list <- list.files("rasters/STM/2020_restor_wo_avoid/", pattern = ".tif", full.names = T, recursive = T)
stm.2020restor_wo_avoid <- stack(stm.2020restor_wo_avoid.raster.list)
names(stm.2020restor_wo_avoid) <- unlist(strsplit(stm.2020restor_wo_avoid.raster.list, "/|.tif"))[seq(4,84,4)]
stm.2020restor_wo_avoid <- stm.2020restor_wo_avoid[[sel.var.df]]

set.seed(999)
mod.mlr.proj_stm.2020restor_wo_avoid <- predict(stm.2020restor_wo_avoid, mod.mlr.fn4)
set.seed(999)
mod.rf.proj_stm.2020restor_wo_avoid <- predict(stm.2020restor_wo_avoid, mod.rf.fn2)


#building a consensus map by mean weight
proj_stm.2020restor_wo_avoid.conbywm <- weighted.mean(stack(mod.mlr.proj_stm.2020restor_wo_avoid, mod.rf.proj_stm.2020restor_wo_avoid), 
                                                      c(RMSE(mod.mlr.fn4$fitted.values, mod.mlr.fn4$y),
                                                        RMSE(mod.rf.fn2$predicted, mod.rf.fn2$y)), 
                                                      na.rm=T)

writeRaster(proj_stm.2020restor_wo_avoid.conbywm, paste0("models.output/carbon.benefits/STM_2020_restor_wo_avoid_carbon_benefit.tif"), format = "GTiff", overwrite = T)

rm(list= ls()[!(ls() %in% c("sel.var.df", "carbon", "env.explanatory.var", "mod.mlr.fn4", "mod.rf.fn2",
                            "proj_pgm.2020real.conbywm", "proj_pgm.2020avoiddegrad.conbywm", "proj_pgm.2020avoiddeforest.conbywm",
                            "proj_pgm.2020avoidboth.conbywm", "proj_pgm.2020restor_wo_avoid.conbywm", "proj_pgm.2020restor_n_avoiddeforest.conbywm",
                            "proj_pgm.2020restor_n_avoidboth.conbywm", "proj_stm.2020real.conbywm", "proj_stm.2020avoiddegrad.conbywm",
                            "proj_stm.2020avoiddeforest.conbywm", "proj_stm.2020avoidboth.conbywm", "proj_stm.2020restor_wo_avoid.conbywm"))])
gc()
#



#STM 2020 restoration and avoid deforest
stm.2020restor_n_avoiddeforest.raster.list <- list.files("rasters/STM/2020_restor_n_avoid_deforest/", pattern = ".tif", full.names = T, recursive = T)
stm.2020restor_n_avoiddeforest <- stack(stm.2020restor_n_avoiddeforest.raster.list)
names(stm.2020restor_n_avoiddeforest) <- unlist(strsplit(stm.2020restor_n_avoiddeforest.raster.list, "/|.tif"))[seq(4,84,4)]
stm.2020restor_n_avoiddeforest <- stm.2020restor_n_avoiddeforest[[sel.var.df]]

set.seed(999)
mod.mlr.proj_stm.2020restor_n_avoiddeforest <- predict(stm.2020restor_n_avoiddeforest, mod.mlr.fn4)
set.seed(999)
mod.rf.proj_stm.2020restor_n_avoiddeforest <- predict(stm.2020restor_n_avoiddeforest, mod.rf.fn2)


#building a consensus map by mean weight
proj_stm.2020restor_n_avoiddeforest.conbywm <- weighted.mean(stack(mod.mlr.proj_stm.2020restor_n_avoiddeforest, mod.rf.proj_stm.2020restor_n_avoiddeforest), 
                                                             c(RMSE(mod.mlr.fn4$fitted.values, mod.mlr.fn4$y),
                                                               RMSE(mod.rf.fn2$predicted, mod.rf.fn2$y)), 
                                                             na.rm=T)

writeRaster(proj_stm.2020restor_n_avoiddeforest.conbywm, paste0("models.output/carbon.benefits/STM_2020_restor_n_avoid_deforest_carbon_benefit.tif"), format = "GTiff", overwrite = T)

rm(list= ls()[!(ls() %in% c("sel.var.df", "carbon", "env.explanatory.var", "mod.mlr.fn4", "mod.rf.fn2",
                            "proj_pgm.2020real.conbywm", "proj_pgm.2020avoiddegrad.conbywm", "proj_pgm.2020avoiddeforest.conbywm",
                            "proj_pgm.2020avoidboth.conbywm", "proj_pgm.2020restor_wo_avoid.conbywm", "proj_pgm.2020restor_n_avoiddeforest.conbywm",
                            "proj_pgm.2020restor_n_avoidboth.conbywm", "proj_stm.2020real.conbywm", "proj_stm.2020avoiddegrad.conbywm",
                            "proj_stm.2020avoiddeforest.conbywm", "proj_stm.2020avoidboth.conbywm", "proj_stm.2020restor_wo_avoid.conbywm",
                            "proj_stm.2020restor_n_avoiddeforest.conbywm"))])
gc()
#



#STM 2020 restoration and avoid both
stm.2020restor_n_avoidboth.raster.list <- list.files("rasters/STM/2020_restor_n_avoid_both/", pattern = ".tif", full.names = T, recursive = T)
stm.2020restor_n_avoidboth <- stack(stm.2020restor_n_avoidboth.raster.list)
names(stm.2020restor_n_avoidboth) <- unlist(strsplit(stm.2020restor_n_avoidboth.raster.list, "/|.tif"))[seq(4,84,4)]
stm.2020restor_n_avoidboth <- stm.2020restor_n_avoidboth[[sel.var.df]]

set.seed(999)
mod.mlr.proj_stm.2020restor_n_avoidboth <- predict(stm.2020restor_n_avoidboth, mod.mlr.fn4)
set.seed(999)
mod.rf.proj_stm.2020restor_n_avoidboth <- predict(stm.2020restor_n_avoidboth, mod.rf.fn2)


#building a consensus map by mean weight
proj_stm.2020restor_n_avoidboth.conbywm <- weighted.mean(stack(mod.mlr.proj_stm.2020restor_n_avoidboth, mod.rf.proj_stm.2020restor_n_avoidboth), 
                                                         c(RMSE(mod.mlr.fn4$fitted.values, mod.mlr.fn4$y),
                                                           RMSE(mod.rf.fn2$predicted, mod.rf.fn2$y)), 
                                                         na.rm=T)

writeRaster(proj_stm.2020restor_n_avoidboth.conbywm, paste0("models.output/carbon.benefits/STM_2020_restor_n_avoid_both_carbon_benefit.tif"), format = "GTiff", overwrite = T)

rm(list= ls()[!(ls() %in% c("sel.var.df", "carbon", "env.explanatory.var", "mod.mlr.fn4", "mod.rf.fn2",
                            "proj_pgm.2020real.conbywm", "proj_pgm.2020avoiddegrad.conbywm", "proj_pgm.2020avoiddeforest.conbywm",
                            "proj_pgm.2020avoidboth.conbywm", "proj_pgm.2020restor_wo_avoid.conbywm", "proj_pgm.2020restor_n_avoiddeforest.conbywm",
                            "proj_pgm.2020restor_n_avoidboth.conbywm", "proj_stm.2020real.conbywm", "proj_stm.2020avoiddegrad.conbywm",
                            "proj_stm.2020avoiddeforest.conbywm", "proj_stm.2020avoidboth.conbywm", "proj_stm.2020restor_wo_avoid.conbywm",
                            "proj_stm.2020restor_n_avoiddeforest.conbywm", "proj_stm.2020restor_n_avoidboth.conbywm"))])
gc()
#










##############################################
####         opportunity cost             ####
##############################################

## property data
property <- read.csv("data/input/areas.csv")
#head(property)
#str(property)
#summary(property)
property[which(is.na(property$totalarea_2009)),"totalarea_2009"] <- property[which(is.na(property$totalarea_2009)),"totalarea_2006"]
property <- property[-which(is.na(property$totalarea_2009)),c("cd_propriedade", "totalarea_2009")]
names(property) <- c("id", "property")

## cost data
revenue <- read.csv("data/input/revenue.csv")
#head(revenue)
#str(revenue)
#summary(revenue)
names(revenue)[1] <- "id"
revenue$catchment <- as.factor(revenue$catchment)

#merging
property <- property %>% left_join(revenue) %>% mutate(profit_ha = profit/property)
table(property$catchment)

#extracting the explanatory variables by catchment
expl.var <- carbon[,c(2,13:25)]
names(expl.var)[1] <- names(property)[4]

#merging
property <- property %>% left_join(expl.var) %>% distinct(id, .keep_all = T)





#scaling predictors -- z-scores
property <- property %>% mutate(dominant = factor(dominant, levels = c("Cattlebeef","Cattlemilk","Animalother","Soy","Annual","Perennial","Mixed","Other")),
                                property_z = ((property - mean(property))/sd(property)),
                                distriver_z = ((distriver - mean(distriver))/sd(distriver)),
                                distroad_z = ((distroad - mean(distroad))/sd(distroad)),
                                DPFpx_z = ((DPFpx - mean(DPFpx))/sd(DPFpx)),
                                edgedist_z = ((edgedist - mean(edgedist))/sd(edgedist)),
                                edgels_z = ((edgels - mean(edgels))/sd(edgels)),
                                elevation_z = ((elevation - mean(elevation))/sd(elevation)),
                                meanprecips_z = ((meanprecips - mean(meanprecips))/sd(meanprecips)),
                                meantemps_z = ((meantemps - mean(meantemps))/sd(meantemps)),
                                SFagels_z = ((SFagels - mean(SFagels))/sd(SFagels)),
                                SFpx_z = ((SFpx - mean(SFpx))/sd(SFpx)),
                                TFpx_z = ((TFpx - mean(TFpx))/sd(TFpx)),
                                UPFls_z = ((UPFls - mean(UPFls))/sd(UPFls)),
                                TSDls_z = ((TSDls - mean(TSDls))/sd(TSDls)))

write.csv(property, "data/property.csv", row.names = F)
#property <- read.csv("data/property.csv")

rm(list= ls()[!(ls() %in% c("property"))])
gc()

#keeping only properties with positive profit
property.total <- property
property <- property %>% dplyr::filter(profit_ha > 0 & !is.nan(profit_ha) & !is.infinite(profit_ha))

#excluding outliers -- two properties with area lower than 2ha but profits greater than $20k/ha:
#id == 419 (PGM), 461 and 322 (STM)
property <- property %>% dplyr::filter(!id %in% c(419, 461, 322))


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
plot.legend <- c("gamma", "lognormal", "gamma")
denscomp(list(property.gamma, property.lnorm, property.gamma), legendtext = plot.legend)
cdfcomp(list(property.gamma, property.lnorm, property.gamma), legendtext = plot.legend)
qqcomp(list(property.gamma, property.lnorm, property.gamma), legendtext = plot.legend)
ppcomp(list(property.gamma, property.lnorm, property.gamma), legendtext = plot.legend)


#plot(Fn <- ecdf(property$profit_ha))
#
#curve(plnorm(x, mean(log(property$profit)), sd(log(property$profit))), from = 128, to = 9000000, add = TRUE, col='darkgreen', lwd = 2)
#
#curve(pgamma(x, shape=6.044639e-02, rate=4.413056e-07), from = 128, to = 9000000, add = TRUE, col='blue', lwd = 2)
#
#curve(pweibull(x, shape=4.160007e-01, scale=2.931072e+04), from = 128, to = 9000000, add = TRUE, col='red', lwd = 2)


# log of profit
property <- property %>% mutate(profit_halog = log(profit_ha))

#null model
mod.null <- lm(profit_halog ~ 1, data = property)
summary(mod.null)




#k-fold
set.seed(999)
group <- kfold(property, k=5, by = property$catchment)
unique(group)

#model fiting and cross-validated predictions for GLM:
eval.mlr <- NULL
models.mlr.list <- list()

for (i in 1:5) {
  train <- property[group != i,]
  test <- property[group == i,]
  crosspred.mlr <- glm(profit_halog ~ dominant + property_z + distriver_z + distroad_z + DPFpx_z + edgedist_z + edgels_z + elevation_z +
                         meanprecips_z + meantemps_z + SFagels_z + SFpx_z + TFpx_z + UPFls_z + TSDls_z,
                       data = train)
  model.sel <- step(crosspred.mlr, test="LRT")
  mod.mlr.fn <- glm(model.sel$formula, data = train)
  models.mlr.list[[i]] <- mod.mlr.fn
  
  predict.mlr <- predict(mod.mlr.fn, test)
  #predict.mlr <- predict(crosspred.mlr, subset(test, non_zero==1))
  eval.mlr[i] <- caret::RMSE(predict.mlr, test$profit_halog)
  
}

model.sel(models.mlr.list)
eval.mlr

#best AICc [==1164.0]; not the best RMSE [==1.83]
mod.mlr.fn1 <- glm(profit_halog ~ property_z + distroad_z + DPFpx_z + edgedist_z + edgels_z + elevation_z +
                     meantemps_z + SFagels_z + UPFls_z + TSDls_z, data = property)
summary(mod.mlr.fn1)
r.squaredLR(mod.mlr.fn1)
plot(mod.mlr.fn1)

#considering non-transformed var
mod.mlr.fn2 <- glm(profit_halog ~ property + distroad + DPFpx + edgedist + edgels + elevation +
                     meantemps + SFagels + UPFls + TSDls, data = property)
summary(mod.mlr.fn2)
r.squaredLR(mod.mlr.fn2)
plot(mod.mlr.fn2)


#visreg(mod.mlr.fn4, xvar = "edgedist_z", data = property)




#model fiting and cross-validated predictions for GLMM -- catchments as random effects
eval.glmm <- NULL
models.glmm.list <- list()

for (i in 1:5) {
  train <- property[group != i,]
  test <- property[group == i,]
  crosspred.glmm <- lmer(profit_halog ~ distriver_z + distroad_z + DPFpx_z + edgedist_z + edgels_z + elevation_z + meanprecips_z +
                           meantemps_z + SFagels_z + SFpx_z + TFpx_z + UPFls_z + TSDls_z + (property_z|dominant),
                         data = train)
  step_res <- step(crosspred.glmm)
  models.glmm.list[[i]] <- get_model(step_res)
  
  predict.glmm <- predict(get_model(step_res), test)
  eval.glmm[i] <- caret::RMSE(predict.glmm, test$profit_halog)
  
}

model.sel(models.glmm.list)
eval.glmm

#best AICc [==1200.4]; not the best RMSE [==1.84]
mod.glmm.fn1 <- lmer(profit_halog ~ distroad_z + edgedist_z + edgels_z + elevation_z + 
                       meantemps_z + SFagels_z + UPFls_z + TSDls_z + (property_z|dominant),
                     data = property)
summary(mod.glmm.fn1)
r.squaredLR(mod.glmm.fn1)
plot(mod.glmm.fn1)

#considering non-transformed var
mod.glmm.fn2 <- lmer(profit_halog ~ distroad + edgedist + edgels + elevation + 
                       meantemps + SFagels + UPFls + TSDls + (property|dominant),
                     data = property)
summary(mod.glmm.fn2)
r.squaredLR(mod.glmm.fn2)
plot(mod.glmm.fn2)


#visreg(mod.mlr.fn4, xvar = "edgedist_z", data = property)




#model fiting for GAM
mod.gam.full <- gam(profit_halog ~ s(property_z, k = 6, bs = "cs") + s(distriver_z, k = 6, bs = "cs") + s(distroad_z, k = 6, bs = "cs") +
                      s(DPFpx_z, k = 6, bs = "cs") + s(edgedist_z, k = 6, bs = "cs") + s(edgels_z, k = 6, bs = "cs") + s(elevation_z, k = 6, bs = "cs") + 
                      s(meanprecips_z, k = 6, bs = "cs") + s(meantemps_z, k = 6, bs = "cs") + s(SFagels_z, k = 6, bs = "cs") + s(SFpx_z, k = 6, bs = "cs") + 
                      s(TFpx_z, k = 6, bs = "cs") + s(UPFls_z, k = 6, bs = "cs") + s(TSDls_z, k = 6, bs = "cs"), method="REML", select = T,
                    data = property)
summary(mod.gam.full)
plot(mod.gam.full, pages=1, residuals=TRUE)
gam.check(mod.gam.full)


mod.gam1 <- gam(profit_halog ~ s(property_z, k = 6, bs = "cs") + s(distriver_z, k = 6, bs = "cs") +
                  s(edgedist_z, k = 6, bs = "cs") + s(edgels_z, k = 6, bs = "cs") + s(elevation_z, k = 6, bs = "cs") + 
                  s(meantemps_z, k = 6, bs = "cs"), method="REML", select = T, data = property)
summary(mod.gam1)
plot(mod.gam1, pages=1, residuals=TRUE)
gam.check(mod.gam1)


mod.gam2 <- gam(profit_halog ~ s(property_z, k = 6, fx = T, bs = "tp") + s(distriver_z, k = 6, fx = T, bs = "tp") +
                  s(edgedist_z, k = 6, fx = T, bs = "tp") + s(edgels_z, k = 6, fx = T, bs = "tp") + 
                  s(elevation_z, k = 6, fx = T, bs = "tp") + s(meantemps_z, k = 6, fx = T, bs = "tp"),
                method="REML", select = T, data = property)
summary(mod.gam2)
plot(mod.gam2, pages=1, residuals=TRUE)
gam.check(mod.gam2)


mod.gam3 <- gam(profit_halog ~ s(property, k = 6, bs = "cs") + s(distriver, k = 6, bs = "cs") +
                  s(edgedist, k = 6, bs = "cs") + s(edgels, k = 6, bs = "cs") + s(elevation, k = 6, bs = "cs") + 
                  s(meantemps, k = 6, bs = "cs"), method="REML", select = T, data = property)
summary(mod.gam3)
plot(mod.gam3, pages=1, residuals=TRUE)
gam.check(mod.gam3)

mod.gam4 <- gam(profit_halog ~ s(property, k = 6, bs = "cs") + s(distroad, k = 6, bs = "cs") +
                  s(edgedist, k = 6, bs = "cs") + s(elevation, k = 6, bs = "cs") + s(meanprecips, k = 6, bs = "cs") +
                  s(meantemps, k = 6, bs = "cs"), method="REML", select = T, data = property)
summary(mod.gam4)
plot(mod.gam4, pages=1, residuals=TRUE)
gam.check(mod.gam4)


model.sel(mod.gam1, mod.gam2, mod.gam3, mod.gam4, mod.glmm.fn1, mod.glmm.fn2, mod.mlr.fn1, mod.mlr.fn2) #



#model fiting and cross-validated predictions for RF
eval.rf <- NULL
models.rf.list <- list()

for (i in 1:5) {
  train <- property[group != i,]
  test <- property[group == i,]
  models.rf.list[[i]] <- randomForest(x = train[,c(3,21:34)], y = train$profit_halog,
                                      xtest = test[,c(3,21:34)], ytest = test$profit_halog,
                                      ntree = 1000, nodesize = 10, importance =T, nPerm = 5)
  
}

par(mfrow=c(3,2))
plot(models.rf.list[[1]])
plot(models.rf.list[[2]])
plot(models.rf.list[[3]])
plot(models.rf.list[[4]])
plot(models.rf.list[[5]])

par(mfrow=c(3,2))
varImpPlot(models.rf.list[[1]], type = 1)
varImpPlot(models.rf.list[[2]], type = 1)
varImpPlot(models.rf.list[[3]], type = 1)
varImpPlot(models.rf.list[[4]], type = 1)
varImpPlot(models.rf.list[[5]], type = 1)

mod.rf.fn1 <- randomForest(x = property[,c(21,29,22,25,26,27)], y = property$profit_halog,
                          ntree=600, nodesize=10, importance =T, nPerm = 5)

par(mfrow=c(1,2))
plot(mod.rf.fn1)
varImpPlot(mod.rf.fn1, type = 1)

mod.rf.fn2 <- randomForest(x = property[, c(2,15,8,12,11,13)], y = property$profit_halog,
                           ntree=600, nodesize=10, importance =T, nPerm = 5)

plot(mod.rf.fn2)
varImpPlot(mod.rf.fn2, type = 1)


mod.rf.fn3 <- randomForest(x = property[,c(2,15,9,14,13,11)], y = property$profit_halog,
                          ntree=600, nodesize=10, importance =T, nPerm = 5)

plot(mod.rf.fn3)
varImpPlot(mod.rf.fn3, type = 1)


rm(list= ls()[!(ls() %in% c("property", "mod.gam3", "mod.rf.fn2"))])
gc()





#predctions
#PGM
dir.create("rasters/PGM/costs", recursive = T)



pgm.2010real.raster.list <- list.files("rasters/PGM/2010_real/", pattern = ".tif", full.names = T, recursive = T)
pgm.2010real <- stack(pgm.2010real.raster.list)
names(pgm.2010real) <- unlist(strsplit(pgm.2010real.raster.list, "/|.tif"))[seq(4,84,4)]

pgm.2010real <- pgm.2010real[[c("property", "distriver", "edgedist", "edgels", "elevation", "meantemps")]]
#pgm.2010real <- pgm.2010real[[c("property", "distriver", "DPFpx", "elevation", "meanprecips", "meantemps", "SFagels")]]


set.seed(999)
mod.gam.proj_pgm.2010real <- exp(predict(pgm.2010real, mod.gam3))
set.seed(999)
mod.rf.proj_pgm.2010real <- exp(predict(pgm.2010real, mod.rf.fn2))


#building a consensus map by mean weight
proj_pgm.2010real.conbywm <- weighted.mean(stack(mod.gam.proj_pgm.2010real, mod.rf.proj_pgm.2010real), 
                                           c(RMSE(mod.gam3$fitted.values, mod.gam3$y),
                                             RMSE(mod.rf.fn2$predicted, mod.rf.fn2$y)), 
                                           na.rm=T)

writeRaster(proj_pgm.2010real.conbywm, paste0("rasters/PGM/costs/PGM_2010_real_base_opportunity_cost.tif"), format = "GTiff", overwrite = T)

rm(list= ls()[!(ls() %in% c("property", "mod.gam3", "mod.rf.fn2", "proj_pgm.2010real.conbywm"))])
gc()
#



#STM
dir.create("rasters/STM/costs", recursive = T)



stm.2010real.raster.list <- list.files("rasters/STM/2010_real/", pattern = ".tif", full.names = T, recursive = T)
stm.2010real <- stack(stm.2010real.raster.list)
names(stm.2010real) <- unlist(strsplit(stm.2010real.raster.list, "/|.tif"))[seq(4,84,4)]

stm.2010real <- stm.2010real[[c("property", "DPFpx", "edgedist", "elevation", "meantemps", "SFagels", "TSDls")]]
#stm.2010real <- stm.2010real[[c("property", "distriver", "DPFpx", "elevation", "meanprecips", "meantemps", "SFagels")]]


set.seed(999)
mod.gam.proj_stm.2010real <- exp(predict(stm.2010real, mod.gam3))
set.seed(999)
mod.rf.proj_stm.2010real <- exp(predict(stm.2010real, mod.rf.fn2))


#building a consensus map by mean weight
proj_stm.2010real.conbywm <- weighted.mean(stack(mod.gam.proj_stm.2010real, mod.rf.proj_stm.2010real), 
                                           c(RMSE(mod.gam3$fitted.values, mod.gam3$y),
                                             RMSE(mod.rf.fn2$predicted, mod.rf.fn2$y)), 
                                           na.rm=T)

writeRaster(proj_stm.2010real.conbywm, paste0("rasters/PGM/costs/STM_2010_real_base_opportunity_cost.tif"), format = "GTiff", overwrite = T)

rm(list= ls()[!(ls() %in% c("property", "mod.gam3", "mod.rf.fn2", "proj_pgm.2010real.conbywm", "proj_stm.2010real.conbywm"))])
gc()
#




# scenario avoid degradation (mean timber value + fire control)
writeRaster(TF.avoiddegrad.mask, "rasters/PGM/all_forest_mask/PGM_2020_avoiddegrad.tif", format = "GTiff", overwrite = T)


pgm.2020avoiddegrad <- sum(TF.avoiddegrad.mask, TF2020.mask, na.rm = T)
##cheking
#sort(unique(values(pgm.2020avoiddegrad)))
pgm.2020avoiddegrad[pgm.2020avoiddegrad>1] <- 0
pgm.2020avoiddegrad[pgm.2020avoiddegrad==1] <- 0
#plot(pgm.2020avoiddegrad)








































































