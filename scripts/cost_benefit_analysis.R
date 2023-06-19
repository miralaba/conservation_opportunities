
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

#########################################
#### comparing biodiversity benefits ####
####      between scenarios and      ####
####         real projections        ####
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
                                      across(c('x', 'y'), round, 5),
                                      Region = "PGM",
                                      Scenario = case_when(str_detect(ID, "avoiddegrad")~ "Avoid degradation",
                                                           str_detect(ID, "avoiddeforest")~ "Avoid deforestation",
                                                           str_detect(ID, "avoidboth")~ "Avoid both",
                                                           str_detect(ID, "restor_wo_avoid")~ "Restoration without avoid",
                                                           str_detect(ID, "restor_n_avoid_deforest")~ "Restoration and avoid deforestation",
                                                           str_detect(ID, "restor_n_avoid_both")~ "Restoration and avoid both",
                                                           str_detect(ID, "real")~ "Real"),
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
                                      across(c('x', 'y'), round, 5),
                                      Region = "STM",
                                      Scenario = case_when(str_detect(ID, "avoiddegrad")~ "Avoid degradation",
                                                           str_detect(ID, "avoiddeforest")~ "Avoid deforestation",
                                                           str_detect(ID, "avoidboth")~ "Avoid both",
                                                           str_detect(ID, "restor_wo_avoid")~ "Restoration without avoid",
                                                           str_detect(ID, "restor_n_avoid_deforest")~ "Restoration and avoid deforestation",
                                                           str_detect(ID, "restor_n_avoid_both")~ "Restoration and avoid both",
                                                           str_detect(ID, "real")~ "Real"),
                                      Biodiversity = na_if(Biodiversity, 0)
                                      )


biodiversity.benefit <- rbind(pgm.biodiversity.benefit.df, stm.biodiversity.benefit.df)

biodiversity.benefit %>% group_by(Region, Scenario) %>% summarise(mean.benefit = mean(Biodiversity, na.rm=T))

cowplot::plot_grid(

  biodiversity.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_avoiddegrad_biodiversity_benefit") %>% 
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Avoid degradation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(ID == "PGM_2020_avoiddegrad_biodiversity_benefit" | ID == "PGM_2020_real_biodiversity_benefit") %>% mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
  ggplot(aes(y=Biodiversity, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=138), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=174), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_restor_wo_avoid_biodiversity_benefit") %>%  
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Restoration without avoid", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_restor_wo_avoid_biodiversity_benefit" | ID == "PGM_2020_real_biodiversity_benefit") %>% mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
    ggplot(aes(y=Biodiversity, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=138), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=143), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_avoiddeforest_biodiversity_benefit") %>%  
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Avoid deforestation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_avoiddeforest_biodiversity_benefit" | ID == "PGM_2020_real_biodiversity_benefit") %>% mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
    ggplot(aes(y=Biodiversity, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=138), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=163), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_restor_n_avoid_deforest_biodiversity_benefit") %>%  
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Restoration and avoid deforestation", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_restor_n_avoid_deforest_biodiversity_benefit" | ID == "PGM_2020_real_biodiversity_benefit") %>% mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
    ggplot(aes(y=Biodiversity, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=138), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=168), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_avoidboth_biodiversity_benefit") %>%  
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Avoid both", x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_avoidboth_biodiversity_benefit" | ID == "PGM_2020_real_biodiversity_benefit") %>% mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
    ggplot(aes(y=Biodiversity, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=138), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=187), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_restor_n_avoid_both_biodiversity_benefit") %>%  
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Restoration and avoid both", x = "Longitude", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_restor_n_avoid_both_biodiversity_benefit" | ID == "PGM_2020_real_biodiversity_benefit") %>% mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
    ggplot(aes(y=Biodiversity, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725"), labels = c("ID", "Real")) +
    scale_fill_manual(values = c("#440154", "#fde725"), labels = c("ID", "Real")) +
    geom_hline(aes(yintercept=138), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=174), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "Biodiversity benefit values") +
    theme_minimal() +
    theme(legend.title = element_blank(), legend.position = c(.8,.8)),

ncol = 4, align = "hv")




cowplot::plot_grid(
  
  biodiversity.benefit %>% filter(Region == "STM" & ID == "STM_2020_avoiddegrad_biodiversity_benefit") %>%  
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Avoid degradation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(ID == "STM_2020_avoiddegrad_biodiversity_benefit" | ID == "STM_2020_real_biodiversity_benefit") %>% mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
    ggplot(aes(y=Biodiversity, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=179), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=206), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(Region == "STM" & ID == "STM_2020_restor_wo_avoid_biodiversity_benefit") %>%  
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Restoration without avoid", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Region == "STM" & ID == "STM_2020_restor_wo_avoid_biodiversity_benefit" | ID == "STM_2020_real_biodiversity_benefit") %>% mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
    ggplot(aes(y=Biodiversity, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=179), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=188), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(Region == "STM" & ID == "STM_2020_avoiddeforest_biodiversity_benefit") %>%  
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Avoid deforestation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Region == "STM" & ID == "STM_2020_avoiddeforest_biodiversity_benefit" | ID == "STM_2020_real_biodiversity_benefit") %>% mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
    ggplot(aes(y=Biodiversity, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=179), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=200), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(Region == "STM" & ID == "STM_2020_restor_n_avoid_deforest_biodiversity_benefit") %>%  
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Restoration and avoid deforestation", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Region == "STM" & ID == "STM_2020_restor_n_avoid_deforest_biodiversity_benefit" | ID == "STM_2020_real_biodiversity_benefit") %>% mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
    ggplot(aes(y=Biodiversity, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=179), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=209), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(Region == "STM" & ID == "STM_2020_avoidboth_biodiversity_benefit") %>%  
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Avoid both", x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Region == "STM" & ID == "STM_2020_avoidboth_biodiversity_benefit" | ID == "STM_2020_real_biodiversity_benefit") %>% mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
    ggplot(aes(y=Biodiversity, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=179), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=222), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(Region == "STM" & ID == "STM_2020_restor_n_avoid_both_biodiversity_benefit") %>%  
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Biodiversity)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Restoration and avoid both", x = "Longitude", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Region == "STM" & ID == "STM_2020_restor_n_avoid_both_biodiversity_benefit" | ID == "STM_2020_real_biodiversity_benefit") %>% mutate(Biodiversity = na_if(Biodiversity, 0)) %>%
    ggplot(aes(y=Biodiversity, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725"), labels = c("ID", "Real")) +
    scale_fill_manual(values = c("#440154", "#fde725"), labels = c("ID", "Real")) +
    geom_hline(aes(yintercept=179), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=214), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "Biodiversity benefit values") +
    theme_minimal() +
    theme(legend.title = element_blank(), legend.position = c(.8,.8)),
  
  ncol = 4, rel_widths = c(1, 1.8), align = "hv")


#rm(list= ls()[!(ls() %in% c("biodiversity.benefit"))])
#gc()



#########################################
####    comparing carbon benefits    ####
####      between scenarios and      ####
####         real projections        ####
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
                              across(c('x', 'y'), round, 5),
                              Region = "PGM",
                              Scenario = case_when(str_detect(ID, "avoiddegrad")~ "Avoid degradation",
                                                   str_detect(ID, "avoiddeforest")~ "Avoid deforestation",
                                                   str_detect(ID, "avoidboth")~ "Avoid both",
                                                   str_detect(ID, "restor_wo_avoid")~ "Restoration without avoid",
                                                   str_detect(ID, "restor_n_avoid_deforest")~ "Restoration and avoid deforestation",
                                                   str_detect(ID, "restor_n_avoid_both")~ "Restoration and avoid both",
                                                   str_detect(ID, "real")~ "Real"),
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
                                 across(c('x', 'y'), round, 5),
                                 Region = "STM",
                                 Scenario = case_when(str_detect(ID, "avoiddegrad")~ "Avoid degradation",
                                                      str_detect(ID, "avoiddeforest")~ "Avoid deforestation",
                                                      str_detect(ID, "avoidboth")~ "Avoid both",
                                                      str_detect(ID, "restor_wo_avoid")~ "Restoration without avoid",
                                                      str_detect(ID, "restor_n_avoid_deforest")~ "Restoration and avoid deforestation",
                                                      str_detect(ID, "restor_n_avoid_both")~ "Restoration and avoid both",
                                                      str_detect(ID, "real")~ "Real"),
                                 Carbon = na_if(Carbon, 0)
                               )


carbon.benefit <- rbind(pgm.carbon.benefit.df, stm.carbon.benefit.df)

carbon.benefit %>% group_by(Region, Scenario) %>% summarise(mean.benefit = mean(Carbon, na.rm=T))

cowplot::plot_grid(
  
  carbon.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_avoiddegrad_carbon_benefit") %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid degradation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(ID == "PGM_2020_avoiddegrad_carbon_benefit" | ID == "PGM_2020_real_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=79.9), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=86.9), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_restor_wo_avoid_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration without avoid", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_restor_wo_avoid_carbon_benefit" | ID == "PGM_2020_real_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=79.9), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=81.9), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_avoiddeforest_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid deforestation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_avoiddeforest_carbon_benefit" | ID == "PGM_2020_real_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=79.9), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=86.5), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_restor_n_avoid_deforest_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration and avoid deforestation", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_restor_n_avoid_deforest_carbon_benefit" | ID == "PGM_2020_real_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=79.9), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=87.7), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_avoidboth_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid both", x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_avoidboth_carbon_benefit" | ID == "PGM_2020_real_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=79.9), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=92.3), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "Carbon benefit values") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_restor_n_avoid_both_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration and avoid both", x = "Longitude", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_restor_n_avoid_both_carbon_benefit" | ID == "PGM_2020_real_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = ID, colour = ID)) +
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
  
  carbon.benefit %>% filter(Region == "STM" & ID == "STM_2020_avoiddegrad_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid degradation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(ID == "STM_2020_avoiddegrad_carbon_benefit" | ID == "STM_2020_real_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=76.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=80.3), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "STM" & ID == "STM_2020_restor_wo_avoid_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration without avoid", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Region == "STM" & ID == "STM_2020_restor_wo_avoid_carbon_benefit" | ID == "STM_2020_real_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=76.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=78.8), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "STM" & ID == "STM_2020_avoiddeforest_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid deforestation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Region == "STM" & ID == "STM_2020_avoiddeforest_carbon_benefit" | ID == "STM_2020_real_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=76.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=81.7), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "STM" & ID == "STM_2020_restor_n_avoid_deforest_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration and avoid deforestation", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Region == "STM" & ID == "STM_2020_restor_n_avoid_deforest_carbon_benefit" | ID == "STM_2020_real_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=76.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=83.9), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "STM" & ID == "STM_2020_avoidboth_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid both", x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Region == "STM" & ID == "STM_2020_avoidboth_carbon_benefit" | ID == "STM_2020_real_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=76.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=88.6), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "Carbon benefit values") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "STM" & ID == "STM_2020_restor_n_avoid_both_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration and avoid both", x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Region == "STM" & ID == "STM_2020_restor_n_avoid_both_carbon_benefit" | ID == "STM_2020_real_carbon_benefit") %>% mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = ID, colour = ID)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725"), labels = c("ID", "Real")) +
    scale_fill_manual(values = c("#374f00", "#fde725"), labels = c("ID", "Real")) +
    geom_hline(aes(yintercept=76.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=81.4), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "Carbon benefit values") +
    theme_minimal() +
    theme(legend.position = c(.8,.8)),
  
  ncol = 4, rel_widths = c(1, 1.8), align = "hv")


#rm(list= ls()[!(ls() %in% c("biodiversity.benefit", "carbon.benefit"))])
#gc()



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
                      across(c('x', 'y'), round, 5),
                      Region = "PGM",
                      Scenario = case_when(str_detect(ID, "avoiddegrad")~ "Avoid degradation",
                                           str_detect(ID, "avoiddeforest")~ "Avoid deforestation",
                                           str_detect(ID, "avoidboth")~ "Avoid both",
                                           str_detect(ID, "restor_wo_avoid")~ "Restoration without avoid",
                                           str_detect(ID, "restor_n_avoid_deforest")~ "Restoration and avoid deforestation",
                                           str_detect(ID, "restor_n_avoid_both")~ "Restoration and avoid both"),
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
                      across(c('x', 'y'), round, 5),
                      Region = "STM",
                      Scenario = case_when(str_detect(ID, "avoiddegrad")~ "Avoid degradation",
                                           str_detect(ID, "avoiddeforest")~ "Avoid deforestation",
                                           str_detect(ID, "avoidboth")~ "Avoid both",
                                           str_detect(ID, "restor_wo_avoid")~ "Restoration without avoid",
                                           str_detect(ID, "restor_n_avoid_deforest")~ "Restoration and avoid deforestation",
                                           str_detect(ID, "restor_n_avoid_both")~ "Restoration and avoid both"),
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

#rm(list= ls()[!(ls() %in% c("biodiversity.benefit", "carbon.benefit", "costs", "benefit.cost.ratio.df"))])
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
  geom_point(aes(size = 5)) +
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
  geom_label_repel(size = 3,
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
  geom_point(aes(size = 5)) +
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
  geom_label_repel(size = 3,
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
#step 1 - filter the cells with conservation action for each scenario
# [the difference in total forests between scenario and 2020 Real]


#step 2 - mask and select cells with higher conservation values
# [top 20% ???]


#step 3 - mask with costs


#step 4 - choose cells ranked by costs until a budget constraint
#and calculate the proportion of land of each scenario in each simulation
#[simulations!]


#















  
  


























