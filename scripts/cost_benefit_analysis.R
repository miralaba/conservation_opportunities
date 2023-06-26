
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
library(ggridges)
#library(ggblend)


#filter the cells with conservation action for each scenario
# [the difference in total forests between scenario and 2020 Real]

pgm.allforest.list <- list.files("rasters/PGM/all_forest_mask/", pattern = ".tif", full.names = T, recursive = T)

pgm.allforest <- stack(pgm.allforest.list)
names(pgm.allforest) <- unlist(strsplit(pgm.allforest.list, "/|.tif"))[seq(4,32,4)]


pgm.conservationaction.avoiddegrad <- pgm.allforest[["PGM_2020_avoiddegrad"]] - pgm.allforest[["PGM_2020_real"]]
pgm.conservationaction.avoiddegrad[pgm.conservationaction.avoiddegrad!=1]<-NA

pgm.conservationaction.avoiddeforest <- pgm.allforest[["PGM_2020_avoiddeforest"]] - pgm.allforest[["PGM_2020_real"]]
pgm.conservationaction.avoiddeforest[pgm.conservationaction.avoiddeforest!=1]<-NA

pgm.conservationaction.avoidboth <- pgm.allforest[["PGM_2020_avoidboth"]] - pgm.allforest[["PGM_2020_real"]]
pgm.conservationaction.avoidboth[pgm.conservationaction.avoidboth!=1]<-NA

pgm.conservationaction.restor_wo_avoid <- pgm.allforest[["PGM_2020_restor_wo_avoid"]] - pgm.allforest[["PGM_2020_real"]]
pgm.conservationaction.restor_wo_avoid[pgm.conservationaction.restor_wo_avoid!=1]<-NA

pgm.conservationaction.restor_n_avoid_deforest <- pgm.allforest[["PGM_2020_restor_n_avoid_deforest"]] - pgm.allforest[["PGM_2020_real"]]
pgm.conservationaction.restor_n_avoid_deforest[pgm.conservationaction.restor_n_avoid_deforest!=1]<-NA

pgm.conservationaction.restor_n_avoid_both <- pgm.allforest[["PGM_2020_restor_n_avoid_both"]] - pgm.allforest[["PGM_2020_real"]]
pgm.conservationaction.restor_n_avoid_both[pgm.conservationaction.restor_n_avoid_both!=1]<-NA


stm.allforest.list <- list.files("rasters/STM/all_forest_mask/", pattern = ".tif", full.names = T, recursive = T)

stm.allforest <- stack(stm.allforest.list)
names(stm.allforest) <- unlist(strsplit(stm.allforest.list, "/|.tif"))[seq(4,32,4)]


stm.conservationaction.avoiddegrad <- stm.allforest[["STM_2020_avoiddegrad"]] - stm.allforest[["STM_2020_real"]]
stm.conservationaction.avoiddegrad[stm.conservationaction.avoiddegrad!=1]<-NA

stm.conservationaction.avoiddeforest <- stm.allforest[["STM_2020_avoiddeforest"]] - stm.allforest[["STM_2020_real"]]
stm.conservationaction.avoiddeforest[stm.conservationaction.avoiddeforest!=1]<-NA

stm.conservationaction.avoidboth <- stm.allforest[["STM_2020_avoidboth"]] - stm.allforest[["STM_2020_real"]]
stm.conservationaction.avoidboth[stm.conservationaction.avoidboth!=1]<-NA

stm.conservationaction.restor_wo_avoid <- stm.allforest[["STM_2020_restor_wo_avoid"]] - stm.allforest[["STM_2020_real"]]
stm.conservationaction.restor_wo_avoid[stm.conservationaction.restor_wo_avoid!=1]<-NA

stm.conservationaction.restor_n_avoid_deforest <- stm.allforest[["STM_2020_restor_n_avoid_deforest"]] - stm.allforest[["STM_2020_real"]]
stm.conservationaction.restor_n_avoid_deforest[stm.conservationaction.restor_n_avoid_deforest!=1]<-NA

stm.conservationaction.restor_n_avoid_both <- stm.allforest[["STM_2020_restor_n_avoid_both"]] - stm.allforest[["STM_2020_real"]]
stm.conservationaction.restor_n_avoid_both[stm.conservationaction.restor_n_avoid_both!=1]<-NA


#########################################
#### comparing biodiversity benefits ####
####      between scenarios and      ####
####         real projections        ####
####           (landscape)           ####
#########################################
biodiversity.benefit.list <- list.files("models.output/biodiversity.benefits/", pattern = ".tif", full.names = T, recursive = T)

pgm.biodiversity.benefit <- stack(grep("PGM", biodiversity.benefit.list, value = T))
names(pgm.biodiversity.benefit) <- unlist(strsplit(biodiversity.benefit.list, "/|.tif"))[seq(3,21,3)]
pgm.biodiversity.benefit.df <- as.data.frame(pgm.biodiversity.benefit, xy = TRUE)
pgm.biodiversity.benefit.df <- pgm.biodiversity.benefit.df %>% 
                                    pivot_longer(
                                      PGM_2020_avoidboth_biodiversity_benefit:PGM_2020_restor_wo_avoid_biodiversity_benefit,
                                      names_to = "ID",
                                      values_to = "Biodiversity"
                                      ) %>% 
                                    mutate(
                                      across(c('x', 'y'), round, 6),
                                      Region = "PGM",
                                      Scenario = factor(case_when(str_detect(ID, "avoiddegrad")~ "Avoid degradation",
                                                                  str_detect(ID, "avoiddeforest")~ "Avoid deforestation",
                                                                  str_detect(ID, "avoidboth")~ "Avoid both",
                                                                  str_detect(ID, "restor_wo_avoid")~ "Restoration without avoid",
                                                                  str_detect(ID, "restor_n_avoid_deforest")~ "Restoration and avoid deforestation",
                                                                  str_detect(ID, "restor_n_avoid_both")~ "Restoration and avoid both",
                                                                  str_detect(ID, "real")~ "Real"),
                                                  levels = c("Avoid degradation",
                                                             "Avoid deforestation",
                                                             "Avoid both",
                                                             "Restoration without avoid",
                                                             "Restoration and avoid deforestation",
                                                             "Restoration and avoid both",
                                                             "Real")),
                                      Biodiversity = na_if(Biodiversity, 0)
                                      )


stm.biodiversity.benefit <- stack(grep("STM", biodiversity.benefit.list, value = T))
names(stm.biodiversity.benefit) <- unlist(strsplit(biodiversity.benefit.list, "/|.tif"))[seq(24,42,3)]
stm.biodiversity.benefit.df <- as.data.frame(stm.biodiversity.benefit, xy = TRUE)
stm.biodiversity.benefit.df <- stm.biodiversity.benefit.df %>% 
                                    pivot_longer(
                                      STM_2020_avoidboth_biodiversity_benefit:STM_2020_restor_wo_avoid_biodiversity_benefit,
                                      names_to = "ID",
                                      values_to = "Biodiversity"
                                      ) %>% 
                                    mutate(
                                      across(c('x', 'y'), round, 6),
                                      Region = "STM",
                                      Scenario = factor(case_when(str_detect(ID, "avoiddegrad")~ "Avoid degradation",
                                                                  str_detect(ID, "avoiddeforest")~ "Avoid deforestation",
                                                                  str_detect(ID, "avoidboth")~ "Avoid both",
                                                                  str_detect(ID, "restor_wo_avoid")~ "Restoration without avoid",
                                                                  str_detect(ID, "restor_n_avoid_deforest")~ "Restoration and avoid deforestation",
                                                                  str_detect(ID, "restor_n_avoid_both")~ "Restoration and avoid both",
                                                                  str_detect(ID, "real")~ "Real"),
                                                        levels = c("Avoid degradation",
                                                                   "Avoid deforestation",
                                                                   "Avoid both",
                                                                   "Restoration without avoid",
                                                                   "Restoration and avoid deforestation",
                                                                   "Restoration and avoid both",
                                                                   "Real")),
                                      Biodiversity = na_if(Biodiversity, 0)
                                      )


biodiversity.benefit <- rbind(pgm.biodiversity.benefit.df, stm.biodiversity.benefit.df)

biodiversity.benefit %>% group_by(Region, Scenario) %>% summarise(mean.benefit = mean(Biodiversity, na.rm=T))

cowplot::plot_grid(

  biodiversity.benefit %>% filter(ID == "PGM_2020_avoiddegrad_biodiversity_benefit") %>% 
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Avoid degradation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(ID == "PGM_2020_avoiddegrad_biodiversity_benefit" | ID == "PGM_2020_real_biodiversity_benefit") %>% 
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot(aes(y=Biodiversity, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=138), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=174), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(ID == "PGM_2020_restor_wo_avoid_biodiversity_benefit") %>%  
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Restoration without avoid", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(ID == "PGM_2020_restor_wo_avoid_biodiversity_benefit" | ID == "PGM_2020_real_biodiversity_benefit") %>% 
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot(aes(y=Biodiversity, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=138), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=143), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(ID == "PGM_2020_avoiddeforest_biodiversity_benefit") %>%  
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Avoid deforestation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(ID == "PGM_2020_avoiddeforest_biodiversity_benefit" | ID == "PGM_2020_real_biodiversity_benefit") %>% 
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot(aes(y=Biodiversity, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=138), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=163), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(ID == "PGM_2020_restor_n_avoid_deforest_biodiversity_benefit") %>%  
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Restoration and avoid deforestation", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(ID == "PGM_2020_restor_n_avoid_deforest_biodiversity_benefit" | ID == "PGM_2020_real_biodiversity_benefit") %>% 
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot(aes(y=Biodiversity, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=138), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=168), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(ID == "PGM_2020_avoidboth_biodiversity_benefit") %>%  
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Avoid both", x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(ID == "PGM_2020_avoidboth_biodiversity_benefit" | ID == "PGM_2020_real_biodiversity_benefit") %>% 
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot(aes(y=Biodiversity, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=138), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=174), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(ID == "PGM_2020_restor_n_avoid_both_biodiversity_benefit") %>%  
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Restoration and avoid both", x = "Longitude", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(ID == "PGM_2020_restor_n_avoid_both_biodiversity_benefit" | ID == "PGM_2020_real_biodiversity_benefit") %>% 
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot(aes(y=Biodiversity, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725"), labels = c("Scenario", "Real")) +
    scale_fill_manual(values = c("#440154", "#fde725"), labels = c("Scenario", "Real")) +
    geom_hline(aes(yintercept=138), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=174), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "Biodiversity benefit values") +
    theme_minimal() +
    theme(legend.title = element_blank(), legend.position = c(.8,.8)),

ncol = 4, align = "hv")




cowplot::plot_grid(
  
  biodiversity.benefit %>% filter(ID == "STM_2020_avoiddegrad_biodiversity_benefit") %>%  
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Avoid degradation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(ID == "STM_2020_avoiddegrad_biodiversity_benefit" | ID == "STM_2020_real_biodiversity_benefit") %>% 
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot(aes(y=Biodiversity, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=179), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=206), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(ID == "STM_2020_restor_wo_avoid_biodiversity_benefit") %>%  
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Restoration without avoid", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(ID == "STM_2020_restor_wo_avoid_biodiversity_benefit" | ID == "STM_2020_real_biodiversity_benefit") %>% 
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot(aes(y=Biodiversity, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=179), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=188), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(ID == "STM_2020_avoiddeforest_biodiversity_benefit") %>%  
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Avoid deforestation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Region == "STM" & ID == "STM_2020_avoiddeforest_biodiversity_benefit" | ID == "STM_2020_real_biodiversity_benefit") %>% 
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot(aes(y=Biodiversity, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=179), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=200), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(ID == "STM_2020_restor_n_avoid_deforest_biodiversity_benefit") %>%  
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Restoration and avoid deforestation", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(ID == "STM_2020_restor_n_avoid_deforest_biodiversity_benefit" | ID == "STM_2020_real_biodiversity_benefit") %>% 
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot(aes(y=Biodiversity, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=179), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=209), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(ID == "STM_2020_avoidboth_biodiversity_benefit") %>%  
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Avoid both", x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(ID == "STM_2020_avoidboth_biodiversity_benefit" | ID == "STM_2020_real_biodiversity_benefit") %>% 
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot(aes(y=Biodiversity, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=179), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=216), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(ID == "STM_2020_restor_n_avoid_both_biodiversity_benefit") %>%  
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Restoration and avoid both", x = "Longitude", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(ID == "STM_2020_restor_n_avoid_both_biodiversity_benefit" | ID == "STM_2020_real_biodiversity_benefit") %>% 
    #mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot(aes(y=Biodiversity, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725"), labels = c("Scenario", "Real")) +
    scale_fill_manual(values = c("#440154", "#fde725"), labels = c("Scenario", "Real")) +
    geom_hline(aes(yintercept=179), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=214), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "Biodiversity benefit values") +
    theme_minimal() +
    theme(legend.title = element_blank(), legend.position = c(.8,.8)),
  
  ncol = 4, rel_widths = c(1, 1.8), align = "hv")


#rm(list= ls()[!(ls() %in% c("pgm.biodiversity.benefit", "stm.biodiversity.benefit", "biodiversity.benefit"))])
#gc()


#########################################
#### comparing biodiversity benefits ####
####      between scenarios and      ####
####         real projections        ####
####  (conservation action areas)    ####
#########################################
##biodiversity
pgm.conservact.avoiddegrad.biodbenefitmask <- mask(pgm.biodiversity.benefit[["PGM_2020_avoiddegrad_biodiversity_benefit"]], pgm.conservationaction.avoiddegrad)
pgm.conservact.avoiddeforest.biodbenefitmask <- mask(pgm.biodiversity.benefit[["PGM_2020_avoiddeforest_biodiversity_benefit"]], pgm.conservationaction.avoiddeforest)
pgm.conservact.avoidboth.biodbenefitmask <- mask(pgm.biodiversity.benefit[["PGM_2020_avoidboth_biodiversity_benefit"]], pgm.conservationaction.avoidboth)
pgm.conservact.restor_wo_avoid.biodbenefitmask <- mask(pgm.biodiversity.benefit[["PGM_2020_restor_wo_avoid_biodiversity_benefit"]], pgm.conservationaction.restor_wo_avoid)
pgm.conservact.restor_n_avoid_deforest.biodbenefitmask <- mask(pgm.biodiversity.benefit[["PGM_2020_restor_n_avoid_deforest_biodiversity_benefit"]], pgm.conservationaction.restor_n_avoid_deforest)
pgm.conservact.restor_n_avoid_both.biodbenefitmask <- mask(pgm.biodiversity.benefit[["PGM_2020_restor_n_avoid_both_biodiversity_benefit"]], pgm.conservationaction.restor_n_avoid_both)

pgm.conservact.biodbenefitmask <- stack(pgm.conservact.avoiddegrad.biodbenefitmask,
                                        pgm.conservact.avoiddeforest.biodbenefitmask,
                                        pgm.conservact.avoidboth.biodbenefitmask,
                                        pgm.conservact.restor_wo_avoid.biodbenefitmask,
                                        pgm.conservact.restor_n_avoid_deforest.biodbenefitmask,
                                        pgm.conservact.restor_n_avoid_both.biodbenefitmask)

names(pgm.conservact.biodbenefitmask) <- c("pgm_avoiddegrad_biodbenefitmask",
                                           "pgm_avoiddeforest_biodbenefitmask",
                                           "pgm_avoidboth_biodbenefitmask",
                                           "pgm_restor_wo_avoid_biodbenefitmask",
                                           "pgm_restor_n_avoid_deforest_biodbenefitmask",
                                           "pgm_restor_n_avoid_both_biodbenefitmask")

pgm.conservact.biodbenefitmask.df <- as.data.frame(pgm.conservact.biodbenefitmask, xy = TRUE)
pgm.conservact.biodbenefitmask.df <- pgm.conservact.biodbenefitmask.df %>% 
  pivot_longer(
    pgm_avoiddegrad_biodbenefitmask:pgm_restor_n_avoid_both_biodbenefitmask,
    names_to = "ID",
    values_to = "Biodiversity"
  ) %>% 
  mutate(
    across(c('x', 'y'), round, 6),
    Region = "PGM",
    Scenario = factor(case_when(str_detect(ID, "avoiddegrad")~ "Avoid degradation",
                                str_detect(ID, "avoiddeforest")~ "Avoid deforestation",
                                str_detect(ID, "avoidboth")~ "Avoid both",
                                str_detect(ID, "restor_wo_avoid")~ "Restoration without avoid",
                                str_detect(ID, "restor_n_avoid_deforest")~ "Restoration and avoid deforestation",
                                str_detect(ID, "restor_n_avoid_both")~ "Restoration and avoid both",
                                str_detect(ID, "real")~ "Real"),
                      levels = c("Avoid degradation",
                                 "Avoid deforestation",
                                 "Avoid both",
                                 "Restoration without avoid",
                                 "Restoration and avoid deforestation",
                                 "Restoration and avoid both",
                                 "Real")),
    Biodiversity = na_if(Biodiversity, 0)
  )


stm.conservact.avoiddegrad.biodbenefitmask <- mask(stm.biodiversity.benefit[["STM_2020_avoiddegrad_biodiversity_benefit"]], stm.conservationaction.avoiddegrad)
stm.conservact.avoiddeforest.biodbenefitmask <- mask(stm.biodiversity.benefit[["STM_2020_avoiddeforest_biodiversity_benefit"]], stm.conservationaction.avoiddeforest)
stm.conservact.avoidboth.biodbenefitmask <- mask(stm.biodiversity.benefit[["STM_2020_avoidboth_biodiversity_benefit"]], stm.conservationaction.avoidboth)
stm.conservact.restor_wo_avoid.biodbenefitmask <- mask(stm.biodiversity.benefit[["STM_2020_restor_wo_avoid_biodiversity_benefit"]], stm.conservationaction.restor_wo_avoid)
stm.conservact.restor_n_avoid_deforest.biodbenefitmask <- mask(stm.biodiversity.benefit[["STM_2020_restor_n_avoid_deforest_biodiversity_benefit"]], stm.conservationaction.restor_n_avoid_deforest)
stm.conservact.restor_n_avoid_both.biodbenefitmask <- mask(stm.biodiversity.benefit[["STM_2020_restor_n_avoid_both_biodiversity_benefit"]], stm.conservationaction.restor_n_avoid_both)

stm.conservact.biodbenefitmask <- stack(stm.conservact.avoiddegrad.biodbenefitmask,
                                        stm.conservact.avoiddeforest.biodbenefitmask,
                                        stm.conservact.avoidboth.biodbenefitmask,
                                        stm.conservact.restor_wo_avoid.biodbenefitmask,
                                        stm.conservact.restor_n_avoid_deforest.biodbenefitmask,
                                        stm.conservact.restor_n_avoid_both.biodbenefitmask)

names(stm.conservact.biodbenefitmask) <- c("stm_avoiddegrad_biodbenefitmask",
                                           "stm_avoiddeforest_biodbenefitmask",
                                           "stm_avoidboth_biodbenefitmask",
                                           "stm_restor_wo_avoid_biodbenefitmask",
                                           "stm_restor_n_avoid_deforest_biodbenefitmask",
                                           "stm_restor_n_avoid_both_biodbenefitmask")

stm.conservact.biodbenefitmask.df <- as.data.frame(stm.conservact.biodbenefitmask, xy = TRUE)
stm.conservact.biodbenefitmask.df <- stm.conservact.biodbenefitmask.df %>% 
  pivot_longer(
    stm_avoiddegrad_biodbenefitmask:stm_restor_n_avoid_both_biodbenefitmask,
    names_to = "ID",
    values_to = "Biodiversity"
  ) %>% 
  mutate(
    across(c('x', 'y'), round, 6),
    Region = "STM",
    Scenario = factor(case_when(str_detect(ID, "avoiddegrad")~ "Avoid degradation",
                                str_detect(ID, "avoiddeforest")~ "Avoid deforestation",
                                str_detect(ID, "avoidboth")~ "Avoid both",
                                str_detect(ID, "restor_wo_avoid")~ "Restoration without avoid",
                                str_detect(ID, "restor_n_avoid_deforest")~ "Restoration and avoid deforestation",
                                str_detect(ID, "restor_n_avoid_both")~ "Restoration and avoid both",
                                str_detect(ID, "real")~ "Real"),
                      levels = c("Avoid degradation",
                                 "Avoid deforestation",
                                 "Avoid both",
                                 "Restoration without avoid",
                                 "Restoration and avoid deforestation",
                                 "Restoration and avoid both",
                                 "Real")),
    Biodiversity = na_if(Biodiversity, 0)
  )


conservact.biodbenefitmask.df <- rbind(pgm.conservact.biodbenefitmask.df, stm.conservact.biodbenefitmask.df)


conservact.biodbenefitmask.df %>% 
  ggplot(aes(y=Biodiversity, fill = Scenario, colour = Scenario)) +
  geom_violin(aes(x=Scenario), show.legend = F) +
  coord_flip() +
  facet_wrap(~Region)


#########################################
####    comparing carbon benefits    ####
####      between scenarios and      ####
####         real projections        ####
####           (landscape)           ####
#########################################
carbon.benefit.list <- list.files("models.output/carbon.benefits/", pattern = ".tif", full.names = T, recursive = T)

pgm.carbon.benefit <- stack(grep("PGM", carbon.benefit.list, value = T))
names(pgm.carbon.benefit) <- unlist(strsplit(carbon.benefit.list, "/|.tif"))[seq(3,21,3)]
pgm.carbon.benefit.df <- as.data.frame(pgm.carbon.benefit, xy = TRUE)
pgm.carbon.benefit.df <- pgm.carbon.benefit.df %>% 
                            pivot_longer(
                              PGM_2020_avoidboth_carbon_benefit:PGM_2020_restor_wo_avoid_carbon_benefit,
                              names_to = "ID",
                              values_to = "Carbon"
                              ) %>% 
                            mutate(
                              across(c('x', 'y'), round, 6),
                              Region = "PGM",
                              Scenario = factor(case_when(str_detect(ID, "avoiddegrad")~ "Avoid degradation",
                                                          str_detect(ID, "avoiddeforest")~ "Avoid deforestation",
                                                          str_detect(ID, "avoidboth")~ "Avoid both",
                                                          str_detect(ID, "restor_wo_avoid")~ "Restoration without avoid",
                                                          str_detect(ID, "restor_n_avoid_deforest")~ "Restoration and avoid deforestation",
                                                          str_detect(ID, "restor_n_avoid_both")~ "Restoration and avoid both",
                                                          str_detect(ID, "real")~ "Real"),
                                                levels = c("Avoid degradation",
                                                           "Avoid deforestation",
                                                           "Avoid both",
                                                           "Restoration without avoid",
                                                           "Restoration and avoid deforestation",
                                                           "Restoration and avoid both",
                                                           "Real")),
                              Carbon = na_if(Carbon, 0)
                              )


stm.carbon.benefit <- stack(grep("STM", carbon.benefit.list, value = T))
names(stm.carbon.benefit) <- unlist(strsplit(carbon.benefit.list, "/|.tif"))[seq(24,42,3)]
stm.carbon.benefit.df <- as.data.frame(stm.carbon.benefit, xy = TRUE)
stm.carbon.benefit.df <- stm.carbon.benefit.df %>% 
                               pivot_longer(
                                 STM_2020_avoidboth_carbon_benefit:STM_2020_restor_wo_avoid_carbon_benefit,
                                 names_to = "ID",
                                 values_to = "Carbon"
                                 ) %>% 
                               mutate(
                                 across(c('x', 'y'), round, 6),
                                 Region = "STM",
                                 Scenario = factor(case_when(str_detect(ID, "avoiddegrad")~ "Avoid degradation",
                                                             str_detect(ID, "avoiddeforest")~ "Avoid deforestation",
                                                             str_detect(ID, "avoidboth")~ "Avoid both",
                                                             str_detect(ID, "restor_wo_avoid")~ "Restoration without avoid",
                                                             str_detect(ID, "restor_n_avoid_deforest")~ "Restoration and avoid deforestation",
                                                             str_detect(ID, "restor_n_avoid_both")~ "Restoration and avoid both",
                                                             str_detect(ID, "real")~ "Real"),
                                                   levels = c("Avoid degradation",
                                                              "Avoid deforestation",
                                                              "Avoid both",
                                                              "Restoration without avoid",
                                                              "Restoration and avoid deforestation",
                                                              "Restoration and avoid both",
                                                              "Real")),
                                 Carbon = na_if(Carbon, 0)
                               )


carbon.benefit <- rbind(pgm.carbon.benefit.df, stm.carbon.benefit.df)

carbon.benefit %>% group_by(Region, Scenario) %>% summarise(mean.benefit = mean(Carbon, na.rm=T))

cowplot::plot_grid(
  
  carbon.benefit %>% filter(ID == "PGM_2020_avoiddegrad_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid degradation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(ID == "PGM_2020_avoiddegrad_carbon_benefit" | ID == "PGM_2020_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
  ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=79.9), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=86.9), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(ID == "PGM_2020_restor_wo_avoid_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>% 
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration without avoid", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(ID == "PGM_2020_restor_wo_avoid_carbon_benefit" | ID == "PGM_2020_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
  ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=79.9), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=81.9), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(ID == "PGM_2020_avoiddeforest_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>% 
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid deforestation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_avoiddeforest_carbon_benefit" | ID == "PGM_2020_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
  ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=79.9), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=86.5), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(ID == "PGM_2020_restor_n_avoid_deforest_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>% 
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration and avoid deforestation", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(ID == "PGM_2020_restor_n_avoid_deforest_carbon_benefit" | ID == "PGM_2020_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
  ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=79.9), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=87.7), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(ID == "PGM_2020_avoidboth_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid both", x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(ID == "PGM_2020_avoidboth_carbon_benefit" | ID == "PGM_2020_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
  ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=79.9), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=91.6), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "Carbon benefit values") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(ID == "PGM_2020_restor_n_avoid_both_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>% 
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration and avoid both", x = "Longitude", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(ID == "PGM_2020_restor_n_avoid_both_carbon_benefit" | ID == "PGM_2020_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
  ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725"), labels = c("ID", "Real")) +
    scale_fill_manual(values = c("#374f00", "#fde725"), labels = c("ID", "Real")) +
    geom_hline(aes(yintercept=79.9), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=85.0), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "Carbon benefit values") +
    theme_minimal() +
    theme(legend.position = c(.8,.8)),
  
  ncol = 4, align = "hv")




cowplot::plot_grid(
  
  carbon.benefit %>% filter(ID == "STM_2020_avoiddegrad_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>% 
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid degradation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(ID == "STM_2020_avoiddegrad_carbon_benefit" | ID == "STM_2020_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
  ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=76.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=80.3), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(ID == "STM_2020_restor_wo_avoid_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>% 
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration without avoid", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(ID == "STM_2020_restor_wo_avoid_carbon_benefit" | ID == "STM_2020_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
  ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=76.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=78.8), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(ID == "STM_2020_avoiddeforest_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>% 
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid deforestation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(ID == "STM_2020_avoiddeforest_carbon_benefit" | ID == "STM_2020_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
  ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=76.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=81.7), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(ID == "STM_2020_restor_n_avoid_deforest_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>% 
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration and avoid deforestation", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(ID == "STM_2020_restor_n_avoid_deforest_carbon_benefit" | ID == "STM_2020_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
  ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=76.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=83.9), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(ID == "STM_2020_avoidboth_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>% 
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid both", x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(ID == "STM_2020_avoidboth_carbon_benefit" | ID == "STM_2020_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
  ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=76.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=85.7), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "Carbon benefit values") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(ID == "STM_2020_restor_n_avoid_both_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>% 
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration and avoid both", x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(ID == "STM_2020_restor_n_avoid_both_carbon_benefit" | ID == "STM_2020_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
  ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725"), labels = c("ID", "Real")) +
    scale_fill_manual(values = c("#374f00", "#fde725"), labels = c("ID", "Real")) +
    geom_hline(aes(yintercept=76.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=81.4), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "Carbon benefit values") +
    theme_minimal() +
    theme(legend.position = c(.8,.8)),
  
  ncol = 4, rel_widths = c(1, 1.8), align = "hv")


#rm(list= ls()[!(ls() %in% c("pgm.biodiversity.benefit", "stm.biodiversity.benefit", "biodiversity.benefit", 
#                            "pgm.carbon.benefit", "stm.carbon.benefit", "carbon.benefit"))])
#gc()


#########################################
####    comparing carbon benefits    ####
####      between scenarios and      ####
####         real projections        ####
####  (conservation action areas)    ####
#########################################
##carbon
pgm.conservact.avoiddegrad.carbbenefitmask <- mask(pgm.carbon.benefit[["PGM_2020_avoiddegrad_carbon_benefit"]], pgm.conservationaction.avoiddegrad)
pgm.conservact.avoiddeforest.carbbenefitmask <- mask(pgm.carbon.benefit[["PGM_2020_avoiddeforest_carbon_benefit"]], pgm.conservationaction.avoiddeforest)
pgm.conservact.avoidboth.carbbenefitmask <- mask(pgm.carbon.benefit[["PGM_2020_avoidboth_carbon_benefit"]], pgm.conservationaction.avoidboth)
pgm.conservact.restor_wo_avoid.carbbenefitmask <- mask(pgm.carbon.benefit[["PGM_2020_restor_wo_avoid_carbon_benefit"]], pgm.conservationaction.restor_wo_avoid)
pgm.conservact.restor_n_avoid_deforest.carbbenefitmask <- mask(pgm.carbon.benefit[["PGM_2020_restor_n_avoid_deforest_carbon_benefit"]], pgm.conservationaction.restor_n_avoid_deforest)
pgm.conservact.restor_n_avoid_both.carbbenefitmask <- mask(pgm.carbon.benefit[["PGM_2020_restor_n_avoid_both_carbon_benefit"]], pgm.conservationaction.restor_n_avoid_both)

pgm.conservact.carbbenefitmask <- stack(pgm.conservact.avoiddegrad.carbbenefitmask,
                                        pgm.conservact.avoiddeforest.carbbenefitmask,
                                        pgm.conservact.avoidboth.carbbenefitmask,
                                        pgm.conservact.restor_wo_avoid.carbbenefitmask,
                                        pgm.conservact.restor_n_avoid_deforest.carbbenefitmask,
                                        pgm.conservact.restor_n_avoid_both.carbbenefitmask)

names(pgm.conservact.carbbenefitmask) <- c("pgm_avoiddegrad_carbbenefitmask",
                                           "pgm_avoiddeforest_carbbenefitmask",
                                           "pgm_avoidboth_carbbenefitmask",
                                           "pgm_restor_wo_avoid_carbbenefitmask",
                                           "pgm_restor_n_avoid_deforest_carbbenefitmask",
                                           "pgm_restor_n_avoid_both_carbbenefitmask")

pgm.conservact.carbbenefitmask.df <- as.data.frame(pgm.conservact.carbbenefitmask, xy = TRUE)
pgm.conservact.carbbenefitmask.df <- pgm.conservact.carbbenefitmask.df %>% 
  pivot_longer(
    pgm_avoiddegrad_carbbenefitmask:pgm_restor_n_avoid_both_carbbenefitmask,
    names_to = "ID",
    values_to = "Carbon"
  ) %>% 
  mutate(
    across(c('x', 'y'), round, 6),
    Region = "PGM",
    Scenario = factor(case_when(str_detect(ID, "avoiddegrad")~ "Avoid degradation",
                                str_detect(ID, "avoiddeforest")~ "Avoid deforestation",
                                str_detect(ID, "avoidboth")~ "Avoid both",
                                str_detect(ID, "restor_wo_avoid")~ "Restoration without avoid",
                                str_detect(ID, "restor_n_avoid_deforest")~ "Restoration and avoid deforestation",
                                str_detect(ID, "restor_n_avoid_both")~ "Restoration and avoid both",
                                str_detect(ID, "real")~ "Real"),
                      levels = c("Avoid degradation",
                                 "Avoid deforestation",
                                 "Avoid both",
                                 "Restoration without avoid",
                                 "Restoration and avoid deforestation",
                                 "Restoration and avoid both",
                                 "Real")),
    Carbon = na_if(Carbon, 0)
  )


stm.conservact.avoiddegrad.carbbenefitmask <- mask(stm.carbon.benefit[["STM_2020_avoiddegrad_carbon_benefit"]], stm.conservationaction.avoiddegrad)
stm.conservact.avoiddeforest.carbbenefitmask <- mask(stm.carbon.benefit[["STM_2020_avoiddeforest_carbon_benefit"]], stm.conservationaction.avoiddeforest)
stm.conservact.avoidboth.carbbenefitmask <- mask(stm.carbon.benefit[["STM_2020_avoidboth_carbon_benefit"]], stm.conservationaction.avoidboth)
stm.conservact.restor_wo_avoid.carbbenefitmask <- mask(stm.carbon.benefit[["STM_2020_restor_wo_avoid_carbon_benefit"]], stm.conservationaction.restor_wo_avoid)
stm.conservact.restor_n_avoid_deforest.carbbenefitmask <- mask(stm.carbon.benefit[["STM_2020_restor_n_avoid_deforest_carbon_benefit"]], stm.conservationaction.restor_n_avoid_deforest)
stm.conservact.restor_n_avoid_both.carbbenefitmask <- mask(stm.carbon.benefit[["STM_2020_restor_n_avoid_both_carbon_benefit"]], stm.conservationaction.restor_n_avoid_both)

stm.conservact.carbbenefitmask <- stack(stm.conservact.avoiddegrad.carbbenefitmask,
                                        stm.conservact.avoiddeforest.carbbenefitmask,
                                        stm.conservact.avoidboth.carbbenefitmask,
                                        stm.conservact.restor_wo_avoid.carbbenefitmask,
                                        stm.conservact.restor_n_avoid_deforest.carbbenefitmask,
                                        stm.conservact.restor_n_avoid_both.carbbenefitmask)

names(stm.conservact.carbbenefitmask) <- c("stm_avoiddegrad_carbbenefitmask",
                                           "stm_avoiddeforest_carbbenefitmask",
                                           "stm_avoidboth_carbbenefitmask",
                                           "stm_restor_wo_avoid_carbbenefitmask",
                                           "stm_restor_n_avoid_deforest_carbbenefitmask",
                                           "stm_restor_n_avoid_both_carbbenefitmask")

stm.conservact.carbbenefitmask.df <- as.data.frame(stm.conservact.carbbenefitmask, xy = TRUE)
stm.conservact.carbbenefitmask.df <- stm.conservact.carbbenefitmask.df %>% 
  pivot_longer(
    stm_avoiddegrad_carbbenefitmask:stm_restor_n_avoid_both_carbbenefitmask,
    names_to = "ID",
    values_to = "Carbon"
  ) %>% 
  mutate(
    across(c('x', 'y'), round, 6),
    Region = "STM",
    Scenario = factor(case_when(str_detect(ID, "avoiddegrad")~ "Avoid degradation",
                                str_detect(ID, "avoiddeforest")~ "Avoid deforestation",
                                str_detect(ID, "avoidboth")~ "Avoid both",
                                str_detect(ID, "restor_wo_avoid")~ "Restoration without avoid",
                                str_detect(ID, "restor_n_avoid_deforest")~ "Restoration and avoid deforestation",
                                str_detect(ID, "restor_n_avoid_both")~ "Restoration and avoid both",
                                str_detect(ID, "real")~ "Real"),
                      levels = c("Avoid degradation",
                                 "Avoid deforestation",
                                 "Avoid both",
                                 "Restoration without avoid",
                                 "Restoration and avoid deforestation",
                                 "Restoration and avoid both",
                                 "Real")),
    Carbon = na_if(Carbon, 0)
  )


conservact.carbbenefitmask.df <- rbind(pgm.conservact.carbbenefitmask.df, stm.conservact.carbbenefitmask.df)


conservact.carbbenefitmask.df %>% 
  ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
  geom_violin(aes(x=Scenario), show.legend = F) +
  coord_flip() +
  facet_wrap(~Region)


##########################################
#### analyzing the benefit-cost ratio ####
##########################################
dir.create("models.output/biodiversity.cost.benefits", recursive = T)
dir.create("models.output/carbon.cost.benefits", recursive = T)

##PGM
#import costs
pgm.opportunity.cost <- raster("models.output/opportunity.costs/PGM_2010_real_base_opportunity_cost.tif")
pgm.opportunity.cost[is.na(pgm.opportunity.cost)]<-0
pgm.firecontrol.cost <- raster("models.output/opportunity.costs/PGM_2010_real_base_firecontrol.tif")
pgm.passiverestor.cost <- raster("models.output/opportunity.costs/PGM_2010_real_base_passiverestoration.tif")


#setting scenario costs
pgm.avoiddegrad.cost <- pgm.opportunity.cost + pgm.firecontrol.cost
pgm.avoiddeforest.cost <- pgm.opportunity.cost
pgm.avoidboth.cost <- pgm.opportunity.cost + pgm.firecontrol.cost
pgm.restor_wo_avoid.cost <- pgm.opportunity.cost + pgm.passiverestor.cost
pgm.restor_n_avoid_deforest.cost <- pgm.opportunity.cost + pgm.passiverestor.cost
pgm.restor_n_avoid_both.cost <- pgm.opportunity.cost + pgm.passiverestor.cost + pgm.firecontrol.cost


#converting in dataframe
pgm.costs.list <- c(pgm.avoiddegrad.cost, pgm.avoiddeforest.cost, pgm.avoidboth.cost, pgm.restor_wo_avoid.cost,
                    pgm.restor_n_avoid_deforest.cost, pgm.restor_n_avoid_both.cost)

pgm.costs <- stack(pgm.costs.list)
names(pgm.costs) <- c("pgm_avoiddegrad_cost", "pgm_avoiddeforest_cost", "pgm_avoidboth_cost", "pgm_restor_wo_avoid_cost",
                    "pgm_restor_n_avoid_deforest.cost", "pgm_restor_n_avoid_both_cost")

pgm.costs.df <- as.data.frame(pgm.costs, xy = TRUE)

pgm.costs.df <- pgm.costs.df %>% 
                    pivot_longer(
                      pgm_avoiddegrad_cost:pgm_restor_n_avoid_both_cost,
                      names_to = "ID",
                      values_to = "Costs"
                      ) %>% 
                    mutate(
                      across(c('x', 'y'), round, 6),
                      Region = "PGM",
                      Scenario = factor(case_when(str_detect(ID, "avoiddegrad")~ "Avoid degradation",
                                                  str_detect(ID, "avoiddeforest")~ "Avoid deforestation",
                                                  str_detect(ID, "avoidboth")~ "Avoid both",
                                                  str_detect(ID, "restor_wo_avoid")~ "Restoration without avoid",
                                                  str_detect(ID, "restor_n_avoid_deforest")~ "Restoration and avoid deforestation",
                                                  str_detect(ID, "restor_n_avoid_both")~ "Restoration and avoid both"),
                                        levels = c("Avoid degradation",
                                                   "Avoid deforestation",
                                                   "Avoid both",
                                                   "Restoration without avoid",
                                                   "Restoration and avoid deforestation",
                                                   "Restoration and avoid both")),
                      Costs = na_if(Costs, 0)
                      )


##STM
#import costs
stm.opportunity.cost <- raster("models.output/opportunity.costs/STM_2010_real_base_opportunity_cost.tif")
stm.opportunity.cost[is.na(stm.opportunity.cost)]<-0
stm.firecontrol.cost <- raster("models.output/opportunity.costs/STM_2010_real_base_firecontrol.tif")
stm.passiverestor.cost <- raster("models.output/opportunity.costs/STM_2010_real_base_passiverestoration.tif")


#setting scenario costs
stm.avoiddegrad.cost <- stm.opportunity.cost + stm.firecontrol.cost
stm.avoiddeforest.cost <- stm.opportunity.cost
stm.avoidboth.cost <- stm.opportunity.cost + stm.firecontrol.cost
stm.restor_wo_avoid.cost <- stm.opportunity.cost + stm.passiverestor.cost
stm.restor_n_avoid_deforest.cost <- stm.opportunity.cost + stm.passiverestor.cost
stm.restor_n_avoid_both.cost <- stm.opportunity.cost + stm.passiverestor.cost + stm.firecontrol.cost


#converting in dataframe
stm.costs.list <- c(stm.avoiddegrad.cost, stm.avoiddeforest.cost, stm.avoidboth.cost, stm.restor_wo_avoid.cost,
                    stm.restor_n_avoid_deforest.cost, stm.restor_n_avoid_both.cost)

stm.costs <- stack(stm.costs.list)
names(stm.costs) <- c("stm_avoiddegrad_cost", "stm_avoiddeforest_cost", "stm_avoidboth_cost", "stm_restor_wo_avoid_cost",
                      "stm_restor_n_avoid_deforest_cost", "stm_restor_n_avoid_both_cost")

stm.costs.df <- as.data.frame(stm.costs, xy = TRUE)

stm.costs.df <- stm.costs.df %>% 
                    pivot_longer(
                      stm_avoiddegrad_cost:stm_restor_n_avoid_both_cost,
                      names_to = "ID",
                      values_to = "Costs"
                      ) %>% 
                    mutate(
                      across(c('x', 'y'), round, 6),
                      Region = "STM",
                      Scenario = factor(case_when(str_detect(ID, "avoiddegrad")~ "Avoid degradation",
                                                  str_detect(ID, "avoiddeforest")~ "Avoid deforestation",
                                                  str_detect(ID, "avoidboth")~ "Avoid both",
                                                  str_detect(ID, "restor_wo_avoid")~ "Restoration without avoid",
                                                  str_detect(ID, "restor_n_avoid_deforest")~ "Restoration and avoid deforestation",
                                                  str_detect(ID, "restor_n_avoid_both")~ "Restoration and avoid both"),
                                        levels = c("Avoid degradation",
                                                   "Avoid deforestation",
                                                   "Avoid both",
                                                   "Restoration without avoid",
                                                   "Restoration and avoid deforestation",
                                                   "Restoration and avoid both")),
                      Costs = na_if(Costs, 0)
                      )


costs <- rbind(pgm.costs.df, stm.costs.df)

costs %>% group_by(Region, Scenario) %>% summarise(mean.costs = mean(Costs, na.rm=T))





benefit.cost.ratio.df <- costs %>% dplyr::select(-ID) %>% 
  left_join(
    (biodiversity.benefit %>% dplyr::select(-ID) %>% filter(Scenario!="Real")),
    by = c("x", "y", "Region", "Scenario")
    )

benefit.cost.ratio.df <- benefit.cost.ratio.df %>% 
  left_join(
    (carbon.benefit %>% dplyr::select(-ID)%>% filter(Scenario!="Real")),
    by = c("x", "y", "Region", "Scenario")
  )

#rm(list= ls()[!(ls() %in% c("pgm.biodiversity.benefit", "stm.biodiversity.benefit", "biodiversity.benefit", 
#                            "pgm.carbon.benefit", "stm.carbon.benefit", "carbon.benefit", 
#                            "pgm.costs", "stm.costs", "costs", "benefit.cost.ratio.df"))])
#gc()
#


benefit.cost.ratio.df <- benefit.cost.ratio.df %>% 
  mutate(Biodiversity.BCR = ifelse(is.na(Costs) | is.na(Biodiversity), NA, Biodiversity/Costs),
         Carbon.BCR = ifelse(is.na(Costs) | is.na(Carbon), NA, Carbon/Costs))

#benefit.cost.ratio.df <- read.csv("data/benefit_cost_ratio.csv")





benefit.cost.ratio.overview <- benefit.cost.ratio.df %>% group_by(Region, Scenario) %>% 
                                    summarise(
                                      Biodiversity_mean = mean(Biodiversity, na.rm=T),
                                      Carbon_mean = mean(Carbon, na.rm=T),
                                      Costs_mean = mean(Costs, na.rm=T),
                                      Biodiversity_sd = sd(Biodiversity, na.rm=T),
                                      Carbon_sd = sd(Carbon, na.rm=T),
                                      Costs_sd = sd(Costs, na.rm=T),
                                      )


library(ggrepel)

empty_theme <- theme(                              
  plot.background = element_blank(), 
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  panel.border = element_blank(), 
  panel.background = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  axis.text.y = element_text(angle = 90)
)


cowplot::plot_grid(

ggplot(data = benefit.cost.ratio.overview, 
       aes(x = Costs_mean, y = Biodiversity_mean, 
           fill = Scenario, 
           shape = Scenario,
           label = Scenario)) + 
  geom_point(aes(size = 8)) +
  #geom_errorbar(aes(ymin = Biodiversity_mean-Biodiversity_sd, ymax = Biodiversity_mean+Biodiversity_sd)) + 
  #geom_errorbarh(aes(xmin = Costs_mean-Costs_sd, xmax = Costs_mean+Costs_sd)) +
  geom_hline(data=filter(benefit.cost.ratio.df, Region=="PGM"), 
             aes(yintercept=163.0597), linetype="dashed") + 
  geom_hline(data=filter(benefit.cost.ratio.df, Region=="STM"), 
             aes(yintercept=201.396), linetype="dashed") +
  geom_vline(data=filter(benefit.cost.ratio.df, Region=="PGM"), 
             aes(xintercept=225.7966), linetype="dashed") + 
  geom_vline(data=filter(benefit.cost.ratio.df, Region=="STM"), 
             aes(xintercept=329.7966), linetype="dashed") +
  scale_x_continuous(breaks = c(160,400), 
                     labels=c("160" = "Low", "400" = "High")) +
  scale_y_continuous(breaks = c(160,210), 
                     labels=c("160" = "Low", "210" = "High")) +
  scale_shape_manual(values = c(21, 21, 21, 24, 24, 24)) +
  scale_fill_manual(values = c("#41644A", "#41644A", "#41644A", 
                               "#D3756B", "#D3756B", "#D3756B")) +
  geom_label_repel(size = 5,
                   fill = NA,
                   label.size = NA,
                   min.segment.length = 0,
                   segment.size = .5,
                   segment.color=NA,
                   box.padding = unit(.25, "lines"),
                   nudge_y = 1.0E-6) +
  facet_grid(~Region)+
  labs(title = "Biodiversity Benefit x Cost",
       x = "Costs",
       y = "Benefits") +
  empty_theme +
  theme(legend.position = "none"),



ggplot(data = benefit.cost.ratio.overview, 
       aes(x = Costs_mean, y = Carbon_mean, 
           fill = Scenario, 
           shape = Scenario,
           label = Scenario)) + 
  geom_point(aes(size = 8)) +
  #geom_errorbar(aes(ymin = Carbon_mean-Carbon_sd, ymax = Carbon_mean+Carbon_sd)) + 
  #geom_errorbarh(aes(xmin = Costs_mean-Costs_sd, xmax = Costs_mean+Costs_sd)) +
  geom_hline(data=filter(benefit.cost.ratio.df, Region=="PGM"), 
             aes(yintercept=85.04721), linetype="dashed") + 
  geom_hline(data=filter(benefit.cost.ratio.df, Region=="STM"), 
             aes(yintercept=80.79686), linetype="dashed") +
  geom_vline(data=filter(benefit.cost.ratio.df, Region=="PGM"), 
             aes(xintercept=225.7966), linetype="dashed") + 
  geom_vline(data=filter(benefit.cost.ratio.df, Region=="STM"), 
             aes(xintercept=329.7966), linetype="dashed") +
  scale_x_continuous(breaks = c(160,400), 
                     labels=c("160" = "Low", "400" = "High")) +
  scale_y_continuous(breaks = c(81,90), 
                     labels=c("81" = "Low", "90" = "High")) +
  scale_shape_manual(values = c(21, 21, 21, 24, 24, 24)) +
  scale_fill_manual(values = c("#41644A", "#41644A", "#41644A", 
                               "#D3756B", "#D3756B", "#D3756B")) +
  geom_label_repel(size = 5,
                   fill = NA,
                   label.size = NA,
                   min.segment.length = 0,
                   segment.size = .5,
                   segment.color=NA,
                   box.padding = unit(.25, "lines"),
                   nudge_y = 1.0E-6) +
  facet_grid(~Region)+
  labs(title = "Carbon Benefit x Cost",
       x = "Costs",
       y = "Benefits") +
  empty_theme +
  theme(legend.position = "none"),

ncol = 1, align = "hv")
  




#next steps




#step 3 - mask with costs

##biodiversity
pgm.conservact.avoiddegrad.biodcostmask <- mask(pgm.costs[["pgm_avoiddegrad_cost"]], pgm.conservact.avoiddegrad.biodbenefitmask)
pgm.conservact.avoiddeforest.biodcostmask <- mask(pgm.costs[["pgm_avoiddeforest_cost"]], pgm.conservact.avoiddeforest.biodbenefitmask)
pgm.conservact.avoidboth.biodcostmask <- mask(pgm.costs[["pgm_avoidboth_cost"]], pgm.conservact.avoidboth.biodbenefitmask)
pgm.conservact.restor_wo_avoid.biodcostmask <- mask(pgm.costs[["pgm_restor_wo_avoid_cost"]], pgm.conservact.restor_wo_avoid.biodbenefitmask)
pgm.conservact.restor_n_avoid_deforest.biodcostmask <- mask(pgm.costs[["pgm_restor_n_avoid_deforest.cost"]], pgm.conservact.restor_n_avoid_deforest.biodbenefitmask)
pgm.conservact.restor_n_avoid_both.biodcostmask <- mask(pgm.costs[["pgm_restor_n_avoid_both_cost"]], pgm.conservact.restor_n_avoid_both.biodbenefitmask)

pgm.conservact.biodcostmask <- stack(pgm.conservact.avoiddegrad.biodcostmask,
                                     pgm.conservact.avoiddeforest.biodcostmask,
                                     pgm.conservact.avoidboth.biodcostmask,
                                     pgm.conservact.restor_wo_avoid.biodcostmask,
                                     pgm.conservact.restor_n_avoid_deforest.biodcostmask,
                                     pgm.conservact.restor_n_avoid_both.biodcostmask)

names(pgm.conservact.carbbenefitmask) <- c("pgm_avoiddegrad_biodcostmask",
                                           "pgm_avoiddeforest_biodcostmask",
                                           "pgm_avoidboth_biodcostmask",
                                           "pgm_restor_wo_avoid_biodcostmask",
                                           "pgm_restor_n_avoid_deforest_biodcostmask",
                                           "pgm_restor_n_avoid_both_biodcostmask")




stm.conservact.avoiddegrad.biodcostmask <- mask(stm.costs[["stm_avoiddegrad_cost"]], stm.conservact.avoiddegrad.biodbenefitmask)
stm.conservact.avoiddeforest.biodcostmask <- mask(stm.costs[["stm_avoiddeforest_cost"]], stm.conservact.avoiddeforest.biodbenefitmask)
stm.conservact.avoidboth.biodcostmask <- mask(stm.costs[["stm_avoidboth_cost"]], stm.conservact.avoidboth.biodbenefitmask)
stm.conservact.restor_wo_avoid.biodcostmask <- mask(stm.costs[["stm_restor_wo_avoid_cost"]], stm.conservact.restor_wo_avoid.biodbenefitmask)
stm.conservact.restor_n_avoid_deforest.biodcostmask <- mask(stm.costs[["stm_restor_n_avoid_deforest_cost"]], stm.conservact.restor_n_avoid_deforest.biodbenefitmask)
stm.conservact.restor_n_avoid_both.biodcostmask <- mask(stm.costs[["stm_restor_n_avoid_both_cost"]], stm.conservact.restor_n_avoid_both.biodbenefitmask)

stm.conservact.biodcostmask <- stack(stm.conservact.avoiddegrad.biodcostmask,
                                     stm.conservact.avoiddeforest.biodcostmask,
                                     stm.conservact.avoidboth.biodcostmask,
                                     stm.conservact.restor_wo_avoid.biodcostmask,
                                     stm.conservact.restor_n_avoid_deforest.biodcostmask,
                                     stm.conservact.restor_n_avoid_both.biodcostmask)

names(stm.conservact.carbbenefitmask) <- c("stm_avoiddegrad_biodcostmask",
                                           "stm_avoiddeforest_biodcostmask",
                                           "stm_avoidboth_biodcostmask",
                                           "stm_restor_wo_avoid_biodcostmask",
                                           "stm_restor_n_avoid_deforest_biodcostmask",
                                           "stm_restor_n_avoid_both_biodcostmask")



##carbon
pgm.conservact.avoiddegrad.carbcostmask <- mask(pgm.costs[["pgm_avoiddegrad_cost"]], pgm.conservact.avoiddegrad.carbbenefitmask)
pgm.conservact.avoiddeforest.carbcostmask <- mask(pgm.costs[["pgm_avoiddeforest_cost"]], pgm.conservact.avoiddeforest.carbbenefitmask)
pgm.conservact.avoidboth.carbcostmask <- mask(pgm.costs[["pgm_avoidboth_cost"]], pgm.conservact.avoidboth.carbbenefitmask)
pgm.conservact.restor_wo_avoid.carbcostmask <- mask(pgm.costs[["pgm_restor_wo_avoid_cost"]], pgm.conservact.restor_wo_avoid.carbbenefitmask)
pgm.conservact.restor_n_avoid_deforest.carbcostmask <- mask(pgm.costs[["pgm_restor_n_avoid_deforest.cost"]], pgm.conservact.restor_n_avoid_deforest.carbbenefitmask)
pgm.conservact.restor_n_avoid_both.carbcostmask <- mask(pgm.costs[["pgm_restor_n_avoid_both_cost"]], pgm.conservact.restor_n_avoid_both.carbbenefitmask)

pgm.conservact.carbcostmask <- stack(pgm.conservact.avoiddegrad.carbcostmask,
                                     pgm.conservact.avoiddeforest.carbcostmask,
                                     pgm.conservact.avoidboth.carbcostmask,
                                     pgm.conservact.restor_wo_avoid.carbcostmask,
                                     pgm.conservact.restor_n_avoid_deforest.carbcostmask,
                                     pgm.conservact.restor_n_avoid_both.carbcostmask)

names(pgm.conservact.carbbenefitmask) <- c("pgm_avoiddegrad_carbcostmask",
                                           "pgm_avoiddeforest_carbcostmask",
                                           "pgm_avoidboth_carbcostmask",
                                           "pgm_restor_wo_avoid_carbcostmask",
                                           "pgm_restor_n_avoid_deforest_carbcostmask",
                                           "pgm_restor_n_avoid_both_carbcostmask")



stm.conservact.avoiddegrad.carbcostmask <- mask(stm.costs[["stm_avoiddegrad_cost"]], stm.conservact.avoiddegrad.carbbenefitmask)
stm.conservact.avoiddeforest.carbcostmask <- mask(stm.costs[["stm_avoiddeforest_cost"]], stm.conservact.avoiddeforest.carbbenefitmask)
stm.conservact.avoidboth.carbcostmask <- mask(stm.costs[["stm_avoidboth_cost"]], stm.conservact.avoidboth.carbbenefitmask)
stm.conservact.restor_wo_avoid.carbcostmask <- mask(stm.costs[["stm_restor_wo_avoid_cost"]], stm.conservact.restor_wo_avoid.carbbenefitmask)
stm.conservact.restor_n_avoid_deforest.carbcostmask <- mask(stm.costs[["stm_restor_n_avoid_deforest_cost"]], stm.conservact.restor_n_avoid_deforest.carbbenefitmask)
stm.conservact.restor_n_avoid_both.carbcostmask <- mask(stm.costs[["stm_restor_n_avoid_both_cost"]], stm.conservact.restor_n_avoid_both.carbbenefitmask)

stm.conservact.carbcostmask <- stack(stm.conservact.avoiddegrad.carbcostmask,
                                     stm.conservact.avoiddeforest.carbcostmask,
                                     stm.conservact.avoidboth.carbcostmask,
                                     stm.conservact.restor_wo_avoid.carbcostmask,
                                     stm.conservact.restor_n_avoid_deforest.carbcostmask,
                                     stm.conservact.restor_n_avoid_both.carbcostmask)

names(stm.conservact.carbbenefitmask) <- c("stm_avoiddegrad_carbcostmask",
                                           "stm_avoiddeforest_carbcostmask",
                                           "stm_avoidboth_carbcostmask",
                                           "stm_restor_wo_avoid_carbcostmask",
                                           "stm_restor_n_avoid_deforest_carbcostmask",
                                           "stm_restor_n_avoid_both_carbcostmask")




#step 4 - choose cells ranked by costs until a budget constraint
#and calculate the proportion of land of each scenario in each simulation
#[simulations!]

constraint.sim <- data.frame(Scenario=NA, N_cells=NA, Regiao=NA, Constraint=NA)

##biodiversity
r <- pgm.conservact.biodcostmask[[c(1,2,4)]]
plot(r)
values(r)[1,]

r.ord <- calc(r, fun=function(x, na.rm) x[order(x, decreasing = F)])
plot(r.ord)
values(r.ord)[1,]

r.scen <- calc(r, function(x, na.rm) order(x, decreasing = F))
r.scen <- mask(r.scen, r.ord[[1]])
plot(r.scen)
values(r.scen)[1,]
table(values(r.scen[[1]]))

ss <- sort(as.vector(r.ord[[1]]), decreasing = F)

for (i in seq(0,100000000,2000000)) {
  
  s_lim <- ss[cumsum(ss) <= i]
  
  constraint <- r.ord[[1]] 
  constraint[!constraint %in% s_lim ] <- NA
  #cellStats(constraint, sum)
  #plot(constraint)
  
  constraint_scen <- mask(r.scen[[1]], constraint)
  #plot(constraint_scen)
  constraint_df <- as.data.frame(table(values(constraint_scen)))
  constraint_df$Regiao <- "PGM"
  constraint_df$Constraint <- i
  colnames(constraint_df) <- c("Scenario", "N_cells", "Regiao", "Constraint")
  
  constraint.sim <- rbind(constraint.sim, constraint_df)
  
}

constraint.sim <- constraint.sim[-1,]



r2 <- stm.conservact.biodcostmask[[c(1,2,4)]]
plot(r2)
values(r2)[1,]

r2.ord <- calc(r2, fun=function(x, na.rm) x[order(x, decreasing = F)])
plot(r2.ord)
values(r2.ord)[1,]

r2.scen <- calc(r2, function(x, na.rm) order(x, decreasing = F))
r2.scen <- mask(r2.scen, r2.ord[[1]])
plot(r2.scen)
values(r2.scen)[1,]
table(values(r2.scen[[1]]))

ss2 <- sort(as.vector(r2.ord[[1]]), decreasing = F)

for (i in seq(0,100000000,2000000)) {
  
  s2_lim <- ss2[cumsum(ss2) <= i]
  
  constraint <- r2.ord[[1]] 
  constraint[!constraint %in% s2_lim ] <- NA
  #cellStats(constraint, sum)
  #plot(constraint)
  
  constraint_scen <- mask(r2.scen[[1]], constraint)
  #plot(constraint_scen)
  constraint_df <- as.data.frame(table(values(constraint_scen)))
  constraint_df$Regiao <- "STM"
  constraint_df$Constraint <- i
  colnames(constraint_df) <- c("Scenario", "N_cells", "Regiao", "Constraint")
  
  constraint.sim <- rbind(constraint.sim, constraint_df)
  
}



constraint.sim$Scenario <- ifelse(constraint.sim$Scenario == 1, "Degradation", 
                                  ifelse(constraint.sim$Scenario == 2, "Deforestation", "Restoration"))



constraint.sim %>% mutate(Scenario = factor(Scenario, levels = c("Degradation", 
                                                                 "Deforestation", 
                                                                 "Restoration"))) %>% 
  ggplot(aes(x=Constraint, y=N_cells, color=Scenario))+
  geom_smooth(se=F) +
  facet_wrap(~Regiao) +
  labs(x="Budget constraint", y="N_cells")+
  theme_bw() +
  theme(plot.title = element_text(size = 14, family = "sans", face = "bold"),
        text = element_text(size = 12, family = "sans"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11))



#















  
  


























