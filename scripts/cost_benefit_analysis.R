
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



biodiversity.benefit.list <- list.files("models.output/biodiversity.benefits/", pattern = ".tif", full.names = T, recursive = T)

pgm.biodiversity.benefit <- stack(grep("PGM", biodiversity.benefit.list, value = T))
names(pgm.biodiversity.benefit) <- unlist(strsplit(biodiversity.benefit.list, "/|.tif"))[seq(3,21,3)]
pgm.biodiversity.benefit.df <- as.data.frame(pgm.biodiversity.benefit, xy = TRUE)
pgm.biodiversity.benefit.df <- pgm.biodiversity.benefit.df %>% 
                                    pivot_longer(PGM_2020_avoidboth_biodiversity_benefit:PGM_2020_restor_wo_avoid_biodiversity_benefit,
                                                 names_to = "Scenario",
                                                 values_to = "Benefit") %>% 
                                    mutate(Region = "PGM")


stm.biodiversity.benefit <- stack(grep("STM", biodiversity.benefit.list, value = T))
names(stm.biodiversity.benefit) <- unlist(strsplit(biodiversity.benefit.list, "/|.tif"))[seq(24,42,3)]
stm.biodiversity.benefit.df <- as.data.frame(stm.biodiversity.benefit, xy = TRUE)
stm.biodiversity.benefit.df <- stm.biodiversity.benefit.df %>% 
                                    pivot_longer(STM_2020_avoidboth_biodiversity_benefit:STM_2020_restor_wo_avoid_biodiversity_benefit,
                                                 names_to = "Scenario",
                                                 values_to = "Benefit") %>% 
                                    mutate(Region = "STM")


biodiversity.benefit <- rbind(pgm.biodiversity.benefit.df, stm.biodiversity.benefit.df)

biodiversity.benefit %>% group_by(Region, Scenario) %>% summarise(mean.benefit = mean(Benefit, na.rm=T))

cowplot::plot_grid(

  biodiversity.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_avoiddegrad_biodiversity_benefit") %>% na_if(0) %>% 
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Avoid degradation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Scenario == "PGM_2020_avoiddegrad_biodiversity_benefit" | Scenario == "PGM_2020_real_biodiversity_benefit") %>% na_if(0) %>%
  ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=60.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=76.3), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_restor_wo_avoid_biodiversity_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Restoration without avoid", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_restor_wo_avoid_biodiversity_benefit" | Scenario == "PGM_2020_real_biodiversity_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=60.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=62.6), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_avoiddeforest_biodiversity_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Avoid deforestation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_avoiddeforest_biodiversity_benefit" | Scenario == "PGM_2020_real_biodiversity_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=60.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=71.4), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_restor_n_avoid_deforest_biodiversity_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Restoration and avoid deforestation", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_restor_n_avoid_deforest_biodiversity_benefit" | Scenario == "PGM_2020_real_biodiversity_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=60.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=73.5), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_avoidboth_biodiversity_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Avoid both", x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_avoidboth_biodiversity_benefit" | Scenario == "PGM_2020_real_biodiversity_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=60.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=78.7), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_restor_n_avoid_both_biodiversity_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Restoration and avoid both", x = "Longitude", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_restor_n_avoid_both_biodiversity_benefit" | Scenario == "PGM_2020_real_biodiversity_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725"), labels = c("Scenario", "Real")) +
    scale_fill_manual(values = c("#440154", "#fde725"), labels = c("Scenario", "Real")) +
    geom_hline(aes(yintercept=60.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=76.2), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "Biodiversity benefit values") +
    theme_minimal() +
    theme(legend.title = element_blank(), legend.position = c(.8,.8)),

ncol = 4, align = "hv")




cowplot::plot_grid(
  
  biodiversity.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_avoiddegrad_biodiversity_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Avoid degradation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Scenario == "STM_2020_avoiddegrad_biodiversity_benefit" | Scenario == "STM_2020_real_biodiversity_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=103), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=118), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_restor_wo_avoid_biodiversity_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Restoration without avoid", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_restor_wo_avoid_biodiversity_benefit" | Scenario == "STM_2020_real_biodiversity_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=103), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=108), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_avoiddeforest_biodiversity_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Avoid deforestation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_avoiddeforest_biodiversity_benefit" | Scenario == "STM_2020_real_biodiversity_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=103), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=115), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_restor_n_avoid_deforest_biodiversity_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Restoration and avoid deforestation", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_restor_n_avoid_deforest_biodiversity_benefit" | Scenario == "STM_2020_real_biodiversity_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=103), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=120), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_avoidboth_biodiversity_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Avoid both", x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_avoidboth_biodiversity_benefit" | Scenario == "STM_2020_real_biodiversity_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725")) +
    scale_fill_manual(values = c("#440154", "#fde725")) +
    geom_hline(aes(yintercept=103), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=120), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  biodiversity.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_restor_n_avoid_both_biodiversity_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_viridis_b(direction=-1, na.value = "white") +
    labs(title = "Restoration and avoid both", x = "Longitude", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  biodiversity.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_restor_n_avoid_both_biodiversity_benefit" | Scenario == "STM_2020_real_biodiversity_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#440154", "#fde725"), labels = c("Scenario", "Real")) +
    scale_fill_manual(values = c("#440154", "#fde725"), labels = c("Scenario", "Real")) +
    geom_hline(aes(yintercept=103), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=123), linewidth=1, linetype='dashed', color = "#440154") +
    labs(title = "", x = "", y = "Biodiversity benefit values") +
    theme_minimal() +
    theme(legend.title = element_blank(), legend.position = c(.8,.8)),
  
  ncol = 4, rel_widths = c(1, 2), align = "hv")



carbon.benefit.list <- list.files("models.output/carbon.benefits/", pattern = ".tif", full.names = T, recursive = T)

pgm.carbon.benefit <- stack(grep("PGM", carbon.benefit.list, value = T))
names(pgm.carbon.benefit) <- unlist(strsplit(carbon.benefit.list, "/|.tif"))[seq(3,21,3)]
pgm.carbon.benefit.df <- as.data.frame(pgm.carbon.benefit, xy = TRUE)
pgm.carbon.benefit.df <- pgm.carbon.benefit.df %>% 
  pivot_longer(PGM_2020_avoidboth_carbon_benefit:PGM_2020_restor_wo_avoid_carbon_benefit,
               names_to = "Scenario",
               values_to = "Benefit") %>% 
  mutate(Region = "PGM")


stm.carbon.benefit <- stack(grep("STM", carbon.benefit.list, value = T))
names(stm.carbon.benefit) <- unlist(strsplit(carbon.benefit.list, "/|.tif"))[seq(24,42,3)]
stm.carbon.benefit.df <- as.data.frame(stm.carbon.benefit, xy = TRUE)
stm.carbon.benefit.df <- stm.carbon.benefit.df %>% 
  pivot_longer(STM_2020_avoidboth_carbon_benefit:STM_2020_restor_wo_avoid_carbon_benefit,
               names_to = "Scenario",
               values_to = "Benefit") %>% 
  mutate(Region = "STM")


carbon.benefit <- rbind(pgm.carbon.benefit.df, stm.carbon.benefit.df)

carbon.benefit %>% group_by(Region, Scenario) %>% summarise(mean.benefit = mean(Benefit, na.rm=T))

cowplot::plot_grid(
  
  carbon.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_avoiddegrad_carbon_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid degradation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Scenario == "PGM_2020_avoiddegrad_carbon_benefit" | Scenario == "PGM_2020_real_carbon_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=79.9), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=86.9), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_restor_wo_avoid_carbon_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration without avoid", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_restor_wo_avoid_carbon_benefit" | Scenario == "PGM_2020_real_carbon_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=79.9), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=81.9), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_avoiddeforest_carbon_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid deforestation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_avoiddeforest_carbon_benefit" | Scenario == "PGM_2020_real_carbon_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=79.9), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=86.5), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_restor_n_avoid_deforest_carbon_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration and avoid deforestation", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_restor_n_avoid_deforest_carbon_benefit" | Scenario == "PGM_2020_real_carbon_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=79.9), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=87.7), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_avoidboth_carbon_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid both", x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_avoidboth_carbon_benefit" | Scenario == "PGM_2020_real_carbon_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=79.9), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=92.3), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "Biodiversity benefit values") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_restor_n_avoid_both_carbon_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration and avoid both", x = "Longitude", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Region == "PGM" & Scenario == "PGM_2020_restor_n_avoid_both_carbon_benefit" | Scenario == "PGM_2020_real_carbon_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725"), labels = c("Scenario", "Real")) +
    scale_fill_manual(values = c("#374f00", "#fde725"), labels = c("Scenario", "Real")) +
    geom_hline(aes(yintercept=79.9), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=85.0), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "Biodiversity benefit values") +
    theme_minimal() +
    theme(legend.position = c(.8,.8)),
  
  ncol = 4, align = "hv")




cowplot::plot_grid(
  
  carbon.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_avoiddegrad_carbon_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid degradation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Scenario == "STM_2020_avoiddegrad_carbon_benefit" | Scenario == "STM_2020_real_carbon_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=76.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=80.3), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_restor_wo_avoid_carbon_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration without avoid", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_restor_wo_avoid_carbon_benefit" | Scenario == "STM_2020_real_carbon_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=76.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=78.8), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_avoiddeforest_carbon_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid deforestation", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_avoiddeforest_carbon_benefit" | Scenario == "STM_2020_real_carbon_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=76.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=81.7), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_restor_n_avoid_deforest_carbon_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration and avoid deforestation", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_restor_n_avoid_deforest_carbon_benefit" | Scenario == "STM_2020_real_carbon_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=76.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=83.9), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_avoidboth_carbon_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid both", x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_avoidboth_carbon_benefit" | Scenario == "STM_2020_real_carbon_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725")) +
    scale_fill_manual(values = c("#374f00", "#fde725")) +
    geom_hline(aes(yintercept=76.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=88.6), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "Biodiversity benefit values") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_restor_n_avoid_both_carbon_benefit") %>% na_if(0) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Benefit)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration and avoid both", x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  carbon.benefit %>% filter(Region == "STM" & Scenario == "STM_2020_restor_n_avoid_both_carbon_benefit" | Scenario == "STM_2020_real_carbon_benefit") %>% na_if(0) %>%
    ggplot(aes(y=Benefit, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#374f00", "#fde725"), labels = c("Scenario", "Real")) +
    scale_fill_manual(values = c("#374f00", "#fde725"), labels = c("Scenario", "Real")) +
    geom_hline(aes(yintercept=76.3), linewidth=1, color = "#fde725") +
    geom_hline(aes(yintercept=81.4), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "", x = "", y = "Biodiversity benefit values") +
    theme_minimal() +
    theme(legend.position = c(.8,.8)),
  
  ncol = 4, rel_widths = c(1, 2), align = "hv")



pgm.opportunity.cost <- raster("models.output/biodiversity.benefits/PGM_cost.TIF")










