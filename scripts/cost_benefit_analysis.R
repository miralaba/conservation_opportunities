
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
library(stringr)
library(ggridges)
library(ggalluvial)
library(ggfittext)
#library(ggblend)

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}


# importing raw rasters ========================================================
## standard projection
std.proj <- "+proj=longlat +datum=WGS84 +units=m +no_defs"

## shapefile paragominas
pgm.shp <- readOGR(dsn = "shapes", layer = "Paragominas_Mask_R3")
proj4string(pgm.shp) <- CRS(std.proj)
pgm.shp <- spTransform(pgm.shp, crs(std.proj))

## shapefile santarem
stm.shp <- readOGR(dsn = "shapes", layer = "Santarem")
stm.shp <- spTransform(stm.shp, crs(std.proj))


# loading areas that change ====================================================
pgm.area.change.list <- list.files("rasters/PGM/input/LULC/", pattern = ".tif", full.names = T, recursive = T)
pgm.area.change.list <- grep("area_change", pgm.area.change.list, value = T)
pgm.area.change.list <- grep("bin", pgm.area.change.list, value = T, invert = T)

pgm.area.change <- stack(pgm.area.change.list)
pgm.area.change <- mask(pgm.area.change, pgm.shp)
#names(pgm.area.change)
#reordering :: real, avoid deforestation, avoid degradation, restoration
pgm.area.change <- pgm.area.change[[c(1,8,4,6,13,2,11,9,5,7,3,12,10)]]

rm(pgm.area.change.list)

#converting to dataframe
pgm.area.change.df.principals <- as.data.frame(pgm.area.change[[3:5]], xy = TRUE)
pgm.area.change.df.principals <- pgm.area.change.df.principals %>% 
  pivot_longer(
    area_change_2020_avoiddeforest:area_change_2020_restor_wo_avoid,
    names_to = "Scenario",
    values_to = "Cat"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "PGM",
    Scenario = factor(Scenario,
                      levels = c("area_change_2020_avoiddeforest",
                                 "area_change_2020_avoiddegrad",
                                 "area_change_2020_restor_wo_avoid"),
                      labels = c("Avoid Deforestation", "Avoid Degradation",
                                 "Restoration without avoid")),
    #Cat = na_if(Cat, 0),
    Cat = factor(Cat,
                 levels = c(1,10,25,100,125),
                 labels = c("UPF", "DPF", "RDPF", "SF", "DSF"))
  ) %>% ungroup() %>% drop_na()
#

stm.area.change.list <- list.files("rasters/STM/input/LULC/", pattern = ".tif", full.names = T, recursive = T)
stm.area.change.list <- grep("area_change", stm.area.change.list, value = T)
stm.area.change.list <- grep("bin", stm.area.change.list, value = T, invert = T)

stm.area.change <- stack(stm.area.change.list)
stm.area.change <- mask(stm.area.change, stm.shp)
stm.area.change <- stm.area.change[[c(1,8,4,6,13,2,11,9,5,7,3,12,10)]]

rm(stm.area.change.list)

stm.area.change.df.principals <- as.data.frame(stm.area.change[[3:5]], xy = TRUE)
stm.area.change.df.principals <- stm.area.change.df.principals %>% 
  pivot_longer(
    area_change_2020_avoiddeforest:area_change_2020_restor_wo_avoid,
    names_to = "Scenario",
    values_to = "Cat"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "STM",
    Scenario = factor(Scenario,
                      levels = c("area_change_2020_avoiddeforest",
                                 "area_change_2020_avoiddegrad",
                                 "area_change_2020_restor_wo_avoid"),
                      labels = c("Avoid Deforestation", "Avoid Degradation",
                                 "Restoration without avoid")),
    #Cat = na_if(Cat, 0),
    Cat = factor(Cat,
                 levels = c(1,10,25,100,125),
                 labels = c("UPF", "DPF", "RDPF", "SF", "DSF"))
  ) %>% ungroup() %>% drop_na()



# comparing biodiversity benefits ====================================================
biodiversity.benefit.list <- list.files("models.output/biodiversity.benefits/", pattern = ".tif", full.names = T, recursive = T)

pgm.biodiversity.benefit.list <- grep("PGM", biodiversity.benefit.list, value = T)
#pgm.biodiversity.benefit2.list <- grep("2_", pgm.biodiversity.benefit.list, value = T)
#pgm.biodiversity.benefit.list <- grep("2_", pgm.biodiversity.benefit.list, value = T, invert = T)

pgm.biodiversity.benefit.total <- stack(pgm.biodiversity.benefit.list)
pgm.biodiversity.benefit.total <- pgm.biodiversity.benefit.total[[c(1,8,4,6,13,2,11,9,5,7,3,12,10)]]
#plot(pgm.biodiversity.benefit.total[[2:13]], nr=3, col = colorRampPalette(c("#FCDAB7", "#1E5F74", "#133B5C", "#1D2D50"))(length(seq(0, 520, by = 10))), breaks= seq(0, 520, by = 10)) ## res = 1673 x 881


rm(pgm.biodiversity.benefit.list)

pgm.biodiversity.benefit <- stack()

for (i in c(3:13)) {
  
  layer.x <- pgm.biodiversity.benefit.total[[1]]
  layer.x[which(!is.na(layer.x[]))]<-NA
  layer.x[] <- pgm.biodiversity.benefit.total[[i]][] - pgm.biodiversity.benefit.total[[2]][]
  layer.x <- mask(layer.x, pgm.area.change[[i]])
  names(layer.x) <- names(pgm.biodiversity.benefit.total[[i]])
  
  pgm.biodiversity.benefit <- addLayer(pgm.biodiversity.benefit, layer.x)
  
  cat("\n>working on layer", i, "now<\n")
}

pgm.biodiversity.benefit <- mask(pgm.biodiversity.benefit, pgm.shp)
#names(pgm.biodiversity.benefit) <- unlist(strsplit(biodiversity.benefit.list, "/|.tif"))[seq(3,21,3)]


pgm.biodiversity.benefit.df.principals <- as.data.frame(pgm.biodiversity.benefit[[1:3]], xy = TRUE)
pgm.biodiversity.benefit.df.principals <- pgm.biodiversity.benefit.df.principals %>% 
  pivot_longer(
    PGM_2020_avoiddeforest_biodiversity_benefit:PGM_2020_restor_wo_avoid_biodiversity_benefit,
    names_to = "Scenario",
    values_to = "BBenefit"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "PGM",
    Scenario = factor(Scenario,
                      levels = c("PGM_2020_avoiddeforest_biodiversity_benefit",
                                 "PGM_2020_avoiddegrad_biodiversity_benefit",
                                 "PGM_2020_restor_wo_avoid_biodiversity_benefit"),
                      labels = c("Avoid Deforestation", "Avoid Degradation",
                                 "Restoration without avoid")),
    #Benefit = na_if(Benefit, 0)
  ) %>% ungroup() %>% drop_na() %>% 
  filter(BBenefit>=0) %>% 
  mutate(rescaled.BBenefit = rescale(BBenefit, to = c(0, 1))
  ) %>% 
  left_join(pgm.area.change.df.principals)
#head(pgm.biodiversity.benefit.df.principals)


#



stm.biodiversity.benefit.list <- grep("STM", biodiversity.benefit.list, value = T)
#stm.biodiversity.benefit2.list <- grep("2_", stm.biodiversity.benefit.list, value = T)
#stm.biodiversity.benefit.list <- grep("2_", stm.biodiversity.benefit.list, value = T, invert = T)

stm.biodiversity.benefit.total <- stack(stm.biodiversity.benefit.list)
stm.biodiversity.benefit.total <- stm.biodiversity.benefit.total[[c(1,8,4,6,13,2,11,9,5,7,3,12,10)]]
#plot(stm.biodiversity.benefit.total[[2:13]], nr=3, col = colorRampPalette(c("#FCDAB7", "#1E5F74", "#133B5C", "#1D2D50"))(length(seq(0, 580, by = 10))), breaks= seq(0, 580, by = 10)) ## res = 1673 x 881

rm(stm.biodiversity.benefit.list); rm(biodiversity.benefit.list)

stm.biodiversity.benefit <- stack()

for (i in c(3:13)) {
  
  layer.x <- stm.biodiversity.benefit.total[[1]]
  layer.x[which(!is.na(layer.x[]))]<-NA
  layer.x[] <- stm.biodiversity.benefit.total[[i]][] - stm.biodiversity.benefit.total[[2]][]
  layer.x <- mask(layer.x, stm.area.change[[i]])
  names(layer.x) <- names(stm.biodiversity.benefit.total[[i]])
  
  stm.biodiversity.benefit <- addLayer(stm.biodiversity.benefit, layer.x)
  
  cat("\n>working on layer", i, "now<\n")
}

stm.biodiversity.benefit <- mask(stm.biodiversity.benefit, stm.shp)
#names(stm.biodiversity.benefit) <- unlist(strsplit(biodiversity.benefit.list, "/|.tif"))[seq(3,21,3)]


stm.biodiversity.benefit.df.principals <- as.data.frame(stm.biodiversity.benefit[[1:3]], xy = TRUE)
stm.biodiversity.benefit.df.principals <- stm.biodiversity.benefit.df.principals %>% 
  pivot_longer(
    STM_2020_avoiddeforest_biodiversity_benefit:STM_2020_restor_wo_avoid_biodiversity_benefit,
    names_to = "Scenario",
    values_to = "BBenefit"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "STM",
    Scenario = factor(Scenario,
                      levels = c("STM_2020_avoiddeforest_biodiversity_benefit",
                                 "STM_2020_avoiddegrad_biodiversity_benefit",
                                 "STM_2020_restor_wo_avoid_biodiversity_benefit"),
                      labels = c("Avoid Deforestation", "Avoid Degradation",
                                 "Restoration without avoid")),
    #Benefit = na_if(Benefit, 0)
  ) %>% ungroup() %>% drop_na() %>% 
  filter(BBenefit>=0) %>% 
  mutate(rescaled.BBenefit = rescale(BBenefit, to = c(0, 1))
  ) %>% 
  left_join(stm.area.change.df.principals)
#head(stm.biodiversity.benefit.df.principals)


#


biodiversity.benefit.principals <- rbind(pgm.biodiversity.benefit.df.principals, stm.biodiversity.benefit.df.principals)
biodiversity.benefit.principals <- biodiversity.benefit.principals %>% mutate(Region = factor(Region, levels = c("PGM", "STM")))




##violin plot + jitter [dots colored by LULC category]
library(ggforce)
g1 <- biodiversity.benefit.principals %>% #filter(BBenefit>=0) %>%
  ggplot(aes(x=Scenario, y=rescaled.BBenefit)) +
  geom_violin(scale="width") + 
  geom_sina(data=biodiversity.benefit.principals %>% sample_n(100000), #%>% filter(BBenefit>=0),
            aes(color=Cat, group=Scenario), alpha=.1, scale="width")+
  scale_x_discrete(labels=addline_format(levels(biodiversity.benefit.principals$Scenario)),
                   expand = c(.05, .05)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("#294B29", "#50623A", "#76453B", "#789461", "#B19470")) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  labs(title = "", x = "", y = "Biodiversity benefit") +
  theme_minimal()+
  theme(text = element_text(size = 16, family = "sans"),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(face="bold"),
        axis.text.x= element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")




#
#



# comparing carbon benefits ====================================================
carbon.benefit.list <- list.files("models.output/carbon.benefits/", pattern = ".tif", full.names = T, recursive = T)

pgm.carbon.benefit.list <- grep("PGM", carbon.benefit.list, value = T)
#pgm.carbon.benefit2.list <- grep("2_", pgm.carbon.benefit.list, value = T)
#pgm.carbon.benefit.list <- grep("2_", pgm.carbon.benefit.list, value = T, invert = T)

pgm.carbon.benefit.total <- stack(pgm.carbon.benefit.list)
pgm.carbon.benefit.total <- pgm.carbon.benefit.total[[c(1,8,4,6,13,2,11,9,5,7,3,12,10)]]
#plot(pgm.carbon.benefit.total[[2:13]], nr=3, col = colorRampPalette(c("#F3F6F4", "#8fce00", "#374f00"))(length(seq(0, 225, by = 25))), breaks= seq(0, 225, by = 25)) ## res = 1673 x 881


rm(pgm.carbon.benefit.list)

pgm.carbon.benefit <- stack()

for (i in c(3:13)) {
  
  layer.x <- pgm.carbon.benefit.total[[1]]
  layer.x[which(!is.na(layer.x[]))]<-NA
  layer.x[] <- pgm.carbon.benefit.total[[i]][] - pgm.carbon.benefit.total[[2]][]
  layer.x <- mask(layer.x, pgm.area.change[[i]])
  names(layer.x) <- names(pgm.carbon.benefit.total[[i]])
  
  pgm.carbon.benefit <- addLayer(pgm.carbon.benefit, layer.x)
  
  cat("\n>working on layer", i, "now<\n")
}

pgm.carbon.benefit <- mask(pgm.carbon.benefit, pgm.shp)
#names(pgm.carbon.benefit) <- unlist(strsplit(carbon.benefit.list, "/|.tif"))[seq(3,21,3)]


pgm.carbon.benefit.df.principals <- as.data.frame(pgm.carbon.benefit[[1:3]], xy = TRUE)
pgm.carbon.benefit.df.principals <- pgm.carbon.benefit.df.principals %>% 
  pivot_longer(
    PGM_2020_avoiddeforest_carbon_benefit:PGM_2020_restor_wo_avoid_carbon_benefit,
    names_to = "Scenario",
    values_to = "CBenefit"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "PGM",
    Scenario = factor(Scenario,
                      levels = c("PGM_2020_avoiddeforest_carbon_benefit",
                                 "PGM_2020_avoiddegrad_carbon_benefit",
                                 "PGM_2020_restor_wo_avoid_carbon_benefit"),
                      labels = c("Avoid Deforestation", "Avoid Degradation",
                                 "Restoration without avoid")),
    #Benefit = na_if(Benefit, 0)
  ) %>% ungroup() %>% drop_na() %>% 
  filter(CBenefit>=0) %>% 
  mutate(rescaled.CBenefit = rescale(CBenefit, to = c(0, 1))
  ) %>% 
  left_join(pgm.area.change.df.principals)
#head(pgm.carbon.benefit.df.principals)


#



stm.carbon.benefit.list <- grep("STM", carbon.benefit.list, value = T)
#stm.carbon.benefit2.list <- grep("2_", stm.carbon.benefit.list, value = T)
#stm.carbon.benefit.list <- grep("2_", stm.carbon.benefit.list, value = T, invert = T)

stm.carbon.benefit.total <- stack(stm.carbon.benefit.list)
stm.carbon.benefit.total <- stm.carbon.benefit.total[[c(1,8,4,6,13,2,11,9,5,7,3,12,10)]]
#plot(stm.carbon.benefit.total[[2:13]], nr=3, col = colorRampPalette(c("#F3F6F4", "#8fce00", "#374f00"))(length(seq(0, 225, by = 25))), breaks= seq(0, 225, by = 25)) ## res = 1673 x 881

rm(stm.carbon.benefit.list); rm(carbon.benefit.list)

stm.carbon.benefit <- stack()

for (i in c(3:13)) {
  
  layer.x <- stm.carbon.benefit.total[[1]]
  layer.x[which(!is.na(layer.x[]))]<-NA
  layer.x[] <- stm.carbon.benefit.total[[i]][] - stm.carbon.benefit.total[[2]][]
  layer.x <- mask(layer.x, stm.area.change[[i]])
  names(layer.x) <- names(stm.carbon.benefit.total[[i]])
  
  stm.carbon.benefit <- addLayer(stm.carbon.benefit, layer.x)
  
  cat("\n>working on layer", i, "now<\n")
}

stm.carbon.benefit <- mask(stm.carbon.benefit, stm.shp)
#names(stm.carbon.benefit) <- unlist(strsplit(carbon.benefit.list, "/|.tif"))[seq(3,21,3)]


stm.carbon.benefit.df.principals <- as.data.frame(stm.carbon.benefit[[1:3]], xy = TRUE)
stm.carbon.benefit.df.principals <- stm.carbon.benefit.df.principals %>% 
  pivot_longer(
    STM_2020_avoiddeforest_carbon_benefit:STM_2020_restor_wo_avoid_carbon_benefit,
    names_to = "Scenario",
    values_to = "CBenefit"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "STM",
    Scenario = factor(Scenario,
                      levels = c("STM_2020_avoiddeforest_carbon_benefit",
                                 "STM_2020_avoiddegrad_carbon_benefit",
                                 "STM_2020_restor_wo_avoid_carbon_benefit"),
                      labels = c("Avoid Deforestation", "Avoid Degradation",
                                 "Restoration without avoid")),
    #Benefit = na_if(Benefit, 0)
  ) %>% ungroup() %>% drop_na() %>% 
  filter(CBenefit>=0) %>% 
  mutate(rescaled.CBenefit = rescale(CBenefit, to = c(0, 1))
  ) %>% 
  left_join(stm.area.change.df.principals)
#head(stm.carbon.benefit.df.principals)


#


carbon.benefit.principals <- rbind(pgm.carbon.benefit.df.principals, stm.carbon.benefit.df.principals)
carbon.benefit.principals <- carbon.benefit.principals %>% mutate(Region = factor(Region, levels = c("PGM", "STM")))




##violin plot + jitter [dots colored by LULC category]
ylim.prim <- c(0, 1)
ylim.sec <- c(0, 120)
b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]

g2 <- carbon.benefit.principals %>% #filter(Benefit>=0) %>%
  ggplot(aes(x=Scenario, y=rescaled.CBenefit)) +
  geom_violin(scale="width") + 
  geom_sina(data=carbon.benefit.principals %>% sample_n(100000), #%>% filter(CBenefit>=0),
            aes(color=Cat, group=Scenario), alpha=.1, scale="width")+
  scale_x_discrete(labels=addline_format(levels(carbon.benefit.principals$Scenario)),
                   expand = c(.05, .05)) +
  scale_y_continuous("Carbon benefit", limits = c(0, 1),
                     sec.axis = sec_axis(~ (. - a)/b, name = "MgC/ha")) +
  scale_color_manual(values = c("#294B29", "#50623A", "#76453B", "#789461", "#B19470")) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  labs(title = "", x = "") +
  theme_minimal()+
  theme(text = element_text(size = 16, family = "sans"),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(face="bold"),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")




#
#



# comparing costs ==============================================================
costs.list <- list.files("models.output/costs/", pattern = ".tif", full.names = T, recursive = T)

pgm.costs.list <- grep("PGM", costs.list, value = T)
#pgm.costs2.list <- grep("2_", pgm.costs.list, value = T)
#pgm.costs.list <- grep("2_", pgm.costs.list, value = T, invert = T)

pgm.costs.total <- stack(pgm.costs.list)
pgm.costs.total[["PGM_2010_real_base_firecontrol"]][] <- ifelse(is.na(pgm.costs.total[["PGM_2010_real_base_firecontrol"]][]), 0, pgm.costs.total[["PGM_2010_real_base_firecontrol"]][])
pgm.costs.total <- mask(pgm.costs.total, pgm.shp)
#names(pgm.costs.total) <- unlist(strsplit(costs.list, "/|.tif"))[seq(3,21,3)]
#plot(pgm.costs.total, nr=2, col = terrain.colors(length(seq(0, 800, by = 50)), rev = T), breaks= seq(0, 800, by = 50)) ## res = 1673 x 881
#

rm(pgm.costs.list)

pgm.max.opportunity.harvest.cost <- calc(pgm.costs.total[[c(2,3)]], function(x){max(x)})



pgm.avoiddeforest.cost.total <- pgm.max.opportunity.harvest.cost
pgm.avoiddeforest.cost <- mask(pgm.avoiddeforest.cost.total, pgm.area.change[["area_change_2020_avoiddeforest"]])
names(pgm.avoiddeforest.cost) <- "PGM_2020_avoiddeforest_costs"
pgm.avoiddeforest2.cost <- mask(pgm.avoiddeforest.cost.total, pgm.area.change[["area_change_2020_avoiddeforest2"]])
names(pgm.avoiddeforest2.cost) <- "PGM_2020_avoiddeforest2_costs"

pgm.avoiddegrad.cost.total <- pgm.costs.total[["PGM_2010_real_base_firecontrol"]] + pgm.costs.total[["PGM_2010_real_base_haverst"]]
pgm.avoiddegrad.cost <- mask(pgm.avoiddegrad.cost.total, pgm.area.change[["area_change_2020_avoiddegrad"]])
names(pgm.avoiddegrad.cost) <- "PGM_2020_avoiddegrad_costs"
pgm.avoiddegrad2.cost <- mask(pgm.avoiddegrad.cost.total, pgm.area.change[["area_change_2020_avoiddegrad2"]])
names(pgm.avoiddegrad2.cost) <- "PGM_2020_avoiddegrad2_costs"

pgm.restor_wo_avoid.cost.total <- pgm.costs.total[["PGM_2010_real_base_opportunity"]] + pgm.costs.total[["PGM_2010_real_base_restoration"]]
pgm.restor_wo_avoid.cost <- mask(pgm.restor_wo_avoid.cost.total, pgm.area.change[["area_change_2020_restor_wo_avoid"]])
names(pgm.restor_wo_avoid.cost) <- "PGM_2020_restor_wo_avoid_costs"

pgm.avoidboth.cost.total <- pgm.costs.total[["PGM_2010_real_base_firecontrol"]] + pgm.max.opportunity.harvest.cost
pgm.avoidboth.cost <- mask(pgm.avoidboth.cost.total, pgm.area.change[["area_change_2020_avoidboth"]])
names(pgm.avoidboth.cost) <- "PGM_2020_avoidboth_costs"
pgm.avoidboth2.cost <- mask(pgm.avoidboth.cost.total, pgm.area.change[["area_change_2020_avoidboth2"]])
names(pgm.avoidboth2.cost) <- "PGM_2020_avoidboth2_costs"

pgm.restor_n_avoid_deforest.cost.total <- pgm.costs.total[["PGM_2010_real_base_restoration"]] + pgm.max.opportunity.harvest.cost
pgm.restor_n_avoid_deforest.cost <- mask(pgm.restor_n_avoid_deforest.cost.total, pgm.area.change[["area_change_2020_restor_n_avoiddeforest"]])
names(pgm.restor_n_avoid_deforest.cost) <- "PGM_2020_restor_n_avoiddeforest_costs"
pgm.restor_n_avoid_deforest2.cost <- mask(pgm.restor_n_avoid_deforest.cost.total, pgm.area.change[["area_change_2020_restor_n_avoiddeforest2"]])
names(pgm.restor_n_avoid_deforest2.cost) <- "PGM_2020_restor_n_avoiddeforest2_costs"

pgm.restor_n_avoid_both.cost.total <- pgm.costs.total[["PGM_2010_real_base_firecontrol"]] + pgm.costs.total[["PGM_2010_real_base_restoration"]] + pgm.max.opportunity.harvest.cost
pgm.restor_n_avoid_both.cost <- mask(pgm.restor_n_avoid_both.cost.total, pgm.area.change[["area_change_2020_restor_n_avoidboth"]])
names(pgm.restor_n_avoid_both.cost) <- "PGM_2020_restor_n_avoidboth_costs"
pgm.restor_n_avoid_both2.cost <- mask(pgm.restor_n_avoid_both.cost.total, pgm.area.change[["area_change_2020_restor_n_avoidboth2"]])
names(pgm.restor_n_avoid_both2.cost) <- "PGM_2020_restor_n_avoidboth2_costs"


pgm.costs <- stack(pgm.avoiddeforest.cost, pgm.avoiddegrad.cost,
                   pgm.restor_wo_avoid.cost, pgm.avoidboth.cost,
                   pgm.restor_n_avoid_deforest.cost, pgm.restor_n_avoid_both.cost, 
                   pgm.avoiddeforest2.cost, pgm.avoiddegrad2.cost, 
                   pgm.avoidboth2.cost, pgm.restor_n_avoid_deforest2.cost,
                   pgm.restor_n_avoid_both2.cost)


rm(list=ls()[ls() %in% c("pgm.max.opportunity.harvest.cost", "pgm.avoiddeforest.cost.total", "pgm.avoiddeforest.cost", "pgm.avoiddeforest2.cost",
                         "pgm.avoiddegrad.cost.total", "pgm.avoiddegrad.cost", "pgm.avoiddegrad2.cost",
                         "pgm.restor_wo_avoid.cost.total","pgm.restor_wo_avoid.cost",
                         "pgm.avoidboth.cost.total", "pgm.avoidboth.cost", "pgm.avoidboth2.cost",
                         "pgm.restor_n_avoid_deforest.cost.total", "pgm.restor_n_avoid_deforest.cost", "pgm.restor_n_avoid_deforest2.cost",
                         "pgm.restor_n_avoid_both.cost.total", "pgm.restor_n_avoid_both.cost", "pgm.restor_n_avoid_both2.cost")])
gc()


pgm.costs.df.principals <- as.data.frame(pgm.costs[[1:3]], xy = TRUE)
pgm.costs.df.principals <- pgm.costs.df.principals %>% 
  pivot_longer(
    PGM_2020_avoiddeforest_costs:PGM_2020_restor_wo_avoid_costs,
    names_to = "Scenario",
    values_to = "Costs"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "PGM",
    Scenario = factor(Scenario,
                levels = c("PGM_2020_avoiddeforest_costs",
                           "PGM_2020_avoiddegrad_costs",
                           "PGM_2020_restor_wo_avoid_costs"),
                labels = c("Avoid Deforestation", "Avoid Degradation",
                           "Restoration without avoid"))
    #Costs = na_if(Costs, 0)
  ) %>% ungroup() %>% drop_na() %>% 
  left_join(pgm.area.change.df.principals)
#head(pgm.costs.df.principals)


#



stm.costs.list <- grep("STM", costs.list, value = T)
#stm.costs2.list <- grep("2_", stm.costs.list, value = T)
#stm.costs.list <- grep("2_", stm.costs.list, value = T, invert = T)

stm.costs.total <- stack(stm.costs.list)
stm.costs.total[["STM_2010_real_base_firecontrol"]][] <- ifelse(is.na(stm.costs.total[["STM_2010_real_base_firecontrol"]][]), 0, stm.costs.total[["STM_2010_real_base_firecontrol"]][])
stm.costs.total <- mask(stm.costs.total, stm.shp)
#names(stm.costs.total) <- unlist(strsplit(costs.list, "/|.tif"))[seq(3,21,3)]
#plot(stm.costs.total, nr=2, col = terrain.colors(length(seq(0, 800, by = 50)), rev = T), breaks= seq(0, 800, by = 50)) ## res = 1673 x 881
#

rm(stm.costs.list); rm(costs.list)

stm.max.opportunity.harvest.cost <- calc(stm.costs.total[[c(2,3)]], function(x){max(x)})



stm.avoiddeforest.cost.total <- stm.max.opportunity.harvest.cost
stm.avoiddeforest.cost <- mask(stm.avoiddeforest.cost.total, stm.area.change[["area_change_2020_avoiddeforest"]])
names(stm.avoiddeforest.cost) <- "STM_2020_avoiddeforest_costs"
stm.avoiddeforest2.cost <- mask(stm.avoiddeforest.cost.total, stm.area.change[["area_change_2020_avoiddeforest2"]])
names(stm.avoiddeforest2.cost) <- "STM_2020_avoiddeforest2_costs"

stm.avoiddegrad.cost.total <- stm.costs.total[["STM_2010_real_base_firecontrol"]] + stm.costs.total[["STM_2010_real_base_haverst"]]
stm.avoiddegrad.cost <- mask(stm.avoiddegrad.cost.total, stm.area.change[["area_change_2020_avoiddegrad"]])
names(stm.avoiddegrad.cost) <- "STM_2020_avoiddegrad_costs"
stm.avoiddegrad2.cost <- mask(stm.avoiddegrad.cost.total, stm.area.change[["area_change_2020_avoiddegrad2"]])
names(stm.avoiddegrad2.cost) <- "STM_2020_avoiddegrad2_costs"

stm.restor_wo_avoid.cost.total <- stm.costs.total[["STM_2010_real_base_opportunity"]] + stm.costs.total[["STM_2010_real_base_restoration"]]
stm.restor_wo_avoid.cost <- mask(stm.restor_wo_avoid.cost.total, stm.area.change[["area_change_2020_restor_wo_avoid"]])
names(stm.restor_wo_avoid.cost) <- "STM_2020_restor_wo_avoid_costs"

stm.avoidboth.cost.total <- stm.costs.total[["STM_2010_real_base_firecontrol"]] + stm.max.opportunity.harvest.cost
stm.avoidboth.cost <- mask(stm.avoidboth.cost.total, stm.area.change[["area_change_2020_avoidboth"]])
names(stm.avoidboth.cost) <- "STM_2020_avoidboth_costs"
stm.avoidboth2.cost <- mask(stm.avoidboth.cost.total, stm.area.change[["area_change_2020_avoidboth2"]])
names(stm.avoidboth2.cost) <- "STM_2020_avoidboth2_costs"

stm.restor_n_avoid_deforest.cost.total <- stm.costs.total[["STM_2010_real_base_restoration"]] + stm.max.opportunity.harvest.cost
stm.restor_n_avoid_deforest.cost <- mask(stm.restor_n_avoid_deforest.cost.total, stm.area.change[["area_change_2020_restor_n_avoiddeforest"]])
names(stm.restor_n_avoid_deforest.cost) <- "STM_2020_restor_n_avoiddeforest_costs"
stm.restor_n_avoid_deforest2.cost <- mask(stm.restor_n_avoid_deforest.cost.total, stm.area.change[["area_change_2020_restor_n_avoiddeforest2"]])
names(stm.restor_n_avoid_deforest2.cost) <- "STM_2020_restor_n_avoiddeforest2_costs"

stm.restor_n_avoid_both.cost.total <- stm.costs.total[["STM_2010_real_base_firecontrol"]] + stm.costs.total[["STM_2010_real_base_restoration"]] + stm.max.opportunity.harvest.cost
stm.restor_n_avoid_both.cost <- mask(stm.restor_n_avoid_both.cost.total, stm.area.change[["area_change_2020_restor_n_avoidboth"]])
names(stm.restor_n_avoid_both.cost) <- "STM_2020_restor_n_avoidboth_costs"
stm.restor_n_avoid_both2.cost <- mask(stm.restor_n_avoid_both.cost.total, stm.area.change[["area_change_2020_restor_n_avoidboth2"]])
names(stm.restor_n_avoid_both2.cost) <- "STM_2020_restor_n_avoidboth2_costs"


stm.costs <- stack(stm.avoiddeforest.cost, stm.avoiddegrad.cost,
                   stm.restor_wo_avoid.cost, stm.avoidboth.cost,
                   stm.restor_n_avoid_deforest.cost, stm.restor_n_avoid_both.cost, 
                   stm.avoiddeforest2.cost, stm.avoiddegrad2.cost, 
                   stm.avoidboth2.cost, stm.restor_n_avoid_deforest2.cost,
                   stm.restor_n_avoid_both2.cost)


rm(list=ls()[ls() %in% c("stm.max.opportunity.harvest.cost", "stm.avoiddeforest.cost.total", "stm.avoiddeforest.cost", "stm.avoiddeforest2.cost",
                         "stm.avoiddegrad.cost.total", "stm.avoiddegrad.cost", "stm.avoiddegrad2.cost",
                         "stm.restor_wo_avoid.cost.total","stm.restor_wo_avoid.cost",
                         "stm.avoidboth.cost.total", "stm.avoidboth.cost", "stm.avoidboth2.cost",
                         "stm.restor_n_avoid_deforest.cost.total", "stm.restor_n_avoid_deforest.cost", "stm.restor_n_avoid_deforest2.cost",
                         "stm.restor_n_avoid_both.cost.total", "stm.restor_n_avoid_both.cost", "stm.restor_n_avoid_both2.cost")])
gc()


stm.costs.df.principals <- as.data.frame(stm.costs[[1:3]], xy = TRUE)
stm.costs.df.principals <- stm.costs.df.principals %>% 
  pivot_longer(
    STM_2020_avoiddeforest_costs:STM_2020_restor_wo_avoid_costs,
    names_to = "Scenario",
    values_to = "Costs"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "STM",
    Scenario = factor(Scenario,
                levels = c("STM_2020_avoiddeforest_costs",
                           "STM_2020_avoiddegrad_costs",
                           "STM_2020_restor_wo_avoid_costs"),
                labels = c("Avoid Deforestation", "Avoid Degradation",
                           "Restoration without avoid"))
    #Costs = na_if(Costs, 0)
  ) %>% ungroup() %>% drop_na() %>% 
  left_join(stm.area.change.df.principals)
#head(stm.costs.df.principals)


#


costs.principals <- rbind(pgm.costs.df.principals, stm.costs.df.principals)
costs.principals <- costs.principals %>% 
  mutate(Region = factor(Region, levels = c("PGM", "STM"))) %>% 
  left_join(biodiversity.benefit.principals) %>% 
  mutate(BBenefit.year = BBenefit/10,
         B.CBr = ifelse(BBenefit.year == 0, NA, Costs/BBenefit.year)) %>% 
  left_join(carbon.benefit.principals) %>% 
  mutate(CBenefit.year = CBenefit/10,
         C.CBr = ifelse(CBenefit.year == 0, NA, Costs/CBenefit.year))




##boxplot
g3 <- costs.principals %>% #filter(Benefit>=0) %>% 
  ggplot(aes(x=Scenario, y=B.CBr)) +
  geom_boxplot(varwidth = T, outlier.shape = NA, fill="#603a62", color="#A4A4A4", show.legend = F) + 
  #geom_hline(yintercept=1, linetype='dashed', linewidth=1)+
  scale_x_discrete(labels=addline_format(levels(costs.principals$Scenario)),
                   expand = c(.05, .05)) +
  scale_y_continuous(limits = c(0, 400)) +
  labs(title = "", x = "", y = "Cost-Benefit ratio (R$ / Biodiversity)") +
  theme_minimal()+
  theme(text = element_text(size = 16, family = "sans"),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "bottom")



g4 <- costs.principals %>% #filter(Benefit>=0) %>% 
  ggplot(aes(x=Scenario, y=C.CBr)) +
  geom_boxplot(varwidth = T, outlier.shape = NA, fill="#603a62", color="#A4A4A4", show.legend = F) + 
  #geom_hline(yintercept=1, linetype='dashed', linewidth=1)+
  scale_x_discrete(labels=addline_format(levels(costs.principals$Scenario)),
                   expand = c(.05, .05)) +
  scale_y_continuous(limits = c(0, 400)) +
  labs(title = "", x = "", y = "Cost-Benefit ratio (R$ / MgC)") +
  theme_minimal()+
  theme(text = element_text(size = 16, family = "sans"),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "bottom")




library(ggpubr)
ggarrange(g1 + theme(plot.margin = margin(1, 20, 1, 1, unit = "pt")), 
          g2, 
          g3 + theme(plot.margin = margin(1, 45, 1, 1, unit = "pt")), 
          g4 + theme(plot.margin = margin(1, 70, 1, 1, unit = "pt")), 
          ncol = 2, nrow = 2, #align = "hv", 
          labels = c("A", "B", "C", "D"), common.legend = T, legend = "top")
#


# best strategy under budget constraints =======================================
cost.biodiv <- costs.principals %>% 
  dplyr::select(Scenario, Costs, Region, BBenefit, rescaled.BBenefit, B.CBr) %>% 
  filter(BBenefit>=0) %>% 
  group_by(Region) %>% 
  arrange(B.CBr) %>% 
  mutate(CumBudget = cumsum(Costs)) %>% 
  ungroup()


result_df <- data.frame(Region = factor(c("PGM", "STM")),
                        Budget = c(0,0),
                        Total_area = c(0,0),
                        Deforestation_area = c(0,0), 
                        Degradation_area = c(0,0), 
                        Restoration_area = c(0,0))


for (b in seq(0,565,.05)) {
  
  suppressMessages(
    result_df <- cost.biodiv %>% 
    group_by(Region, Scenario) %>% 
    filter(CumBudget<=b*1000000) %>% 
    mutate(Scenario = factor(case_when(str_detect(Scenario, "Deforestation")~ "Deforestation_area",
                                       str_detect(Scenario, "Degradation")~ "Degradation_area",
                                       str_detect(Scenario, "Restoration")~ "Restoration_area"),
                       levels = c("Deforestation_area",
                                  "Degradation_area",
                                  "Restoration_area"))) %>% 
    summarise(ncell = n()) %>% 
    pivot_wider(names_from = Scenario, values_from = ncell) %>% 
    ungroup() %>% 
    mutate(Budget = b*1000000,
           Total_area = rowSums(select_if(., is.numeric), na.rm=T)) %>% 
    full_join(result_df)
    )
  
  #cat("\n\t> R$", b*1000000, "budget constraint -- Done!")
  
}




cost.carb <- costs.principals %>% 
  dplyr::select(Scenario, Costs, Region, CBenefit, C.CBr) %>% 
  filter(CBenefit>=0) %>% 
  group_by(Region) %>% 
  arrange(C.CBr) %>% 
  mutate(CumBudget = cumsum(Costs)) %>% 
  ungroup()



for (b in seq(0,565,.05)) {
  
  suppressMessages(
    result_df <- cost.carb %>% 
      group_by(Region, Scenario) %>% 
      filter(CumBudget<=b*1000000) %>% 
      mutate(Scenario = factor(case_when(str_detect(Scenario, "Deforestation")~ "Deforestation_area",
                                         str_detect(Scenario, "Degradation")~ "Degradation_area",
                                         str_detect(Scenario, "Restoration")~ "Restoration_area"),
                               levels = c("Deforestation_area",
                                          "Degradation_area",
                                          "Restoration_area"))) %>% 
      summarise(ncell = n()) %>% 
      pivot_wider(names_from = Scenario, values_from = ncell) %>% 
      ungroup() %>% 
      mutate(Budget = b*1000000,
             Total_area = rowSums(select_if(., is.numeric), na.rm=T)) %>% 
      full_join(result_df)
  )
  
  #cat("\n\t> R$", b*1000000, "budget constraint -- Done!")
  
}


result_df$benefit <- NA
result_df[1:22600,"benefit"] <- "Carbon"
result_df[22601:nrow(result_df),"benefit"] <- "Biodiversity"

pgm.bio.def.area <- result_df %>% filter(Region=="PGM" & benefit=="Biodiversity") %>% first() %>% dplyr::select(Deforestation_area) %>% pull()
pgm.carb.def.area <- result_df %>% filter(Region=="PGM" & benefit=="Carbon") %>% first() %>% dplyr::select(Deforestation_area) %>% pull()
pgm.bio.deg.area <- result_df %>% filter(Region=="PGM" & benefit=="Biodiversity") %>% first() %>% dplyr::select(Degradation_area) %>% pull()
pgm.carb.deg.area <- result_df %>% filter(Region=="PGM" & benefit=="Carbon") %>% first() %>% dplyr::select(Degradation_area) %>% pull()
pgm.bio.rest.area <- result_df %>% filter(Region=="PGM" & benefit=="Biodiversity") %>% first() %>% dplyr::select(Restoration_area) %>% pull()
pgm.carb.rest.area <- result_df %>% filter(Region=="PGM" & benefit=="Carbon") %>% first() %>% dplyr::select(Restoration_area) %>% pull()


stm.bio.def.area <- result_df %>% filter(Region=="STM" & benefit=="Biodiversity") %>% first() %>% dplyr::select(Deforestation_area) %>% pull()
stm.carb.def.area <- result_df %>% filter(Region=="STM" & benefit=="Carbon") %>% first() %>% dplyr::select(Deforestation_area) %>% pull()
stm.bio.deg.area <- result_df %>% filter(Region=="STM" & benefit=="Biodiversity") %>% first() %>% dplyr::select(Degradation_area) %>% pull()
stm.carb.deg.area <- result_df %>% filter(Region=="STM" & benefit=="Carbon") %>% first() %>% dplyr::select(Degradation_area) %>% pull()
stm.bio.rest.area <- result_df %>% filter(Region=="STM" & benefit=="Biodiversity") %>% first() %>% dplyr::select(Restoration_area) %>% pull()
stm.carb.rest.area <- result_df %>% filter(Region=="STM" & benefit=="Carbon") %>% first() %>% dplyr::select(Restoration_area) %>% pull()


result_df <- result_df %>% 
  mutate(
    Deforestation_area_pp = ifelse(Region=="PGM" & benefit=="Biodiversity", Deforestation_area/pgm.bio.def.area,
                                   ifelse(Region=="PGM" & benefit=="Carbon", Deforestation_area/pgm.carb.def.area,
                                          ifelse(Region=="STM" & benefit=="Biodiversity", Deforestation_area/stm.bio.def.area,
                                                 Deforestation_area/stm.carb.def.area))),
    Degradation_area_pp = ifelse(Region=="PGM" & benefit=="Biodiversity", Degradation_area/pgm.bio.deg.area,
                                   ifelse(Region=="PGM" & benefit=="Carbon", Degradation_area/pgm.carb.deg.area,
                                          ifelse(Region=="STM" & benefit=="Biodiversity", Degradation_area/stm.bio.deg.area,
                                                 Degradation_area/stm.carb.deg.area))),
    Restoration_area_pp = ifelse(Region=="PGM" & benefit=="Biodiversity", Restoration_area/pgm.bio.rest.area,
                                 ifelse(Region=="PGM" & benefit=="Carbon", Restoration_area/pgm.carb.rest.area,
                                        ifelse(Region=="STM" & benefit=="Biodiversity", Restoration_area/stm.bio.rest.area,
                                               Restoration_area/stm.carb.rest.area)))
  ) %>% 
  group_by(Region, benefit) %>% 
  arrange(Budget) %>% 
  mutate(
    Total = Total_area - lag(Total_area, default = first(Total_area)),
    Degradation = Degradation_area - lag(Degradation_area, default = first(Degradation_area)),
    Deforestation = Deforestation_area - lag(Deforestation_area, default = first(Deforestation_area)),
    Restoration = Restoration_area - lag(Restoration_area, default = first(Restoration_area)),
    Proportion_Degradation = Degradation/Total,
    Proportion_Deforestation = Deforestation/Total,
    Proportion_Restoration = Restoration/Total
  )




##create a plot
g5 <- result_df %>%  
  ggplot(aes(x = Budget)) +
  #geom_line(aes(y = Proportion_Degradation, color = "Avoid Degradation"), linewidth = 1) +
  geom_smooth(aes(y = Proportion_Degradation, color = "Avoid Degradation"), linewidth = 2, method = "gam") +
  #geom_line(aes(y = Proportion_Deforestation, color = "Avoid Deforestation"), linewidth = 1) +
  geom_smooth(aes(y = Proportion_Deforestation, color = "Avoid Deforestation"), linewidth = 2, method = "gam") +
  #geom_line(aes(y = Proportion_Restoration, color = "Passive Restoration"), linewidth = 1) +
  geom_smooth(aes(y = Proportion_Restoration, color = "Restoration"), linewidth = 2, method = "gam") +
  scale_color_manual(values = c("Avoid Deforestation" = "#294B29", "Restoration" = "#789461", "Avoid Degradation" = "#76453B")) +
  facet_wrap(~benefit, ncol=1) +
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  scale_x_continuous(limits = c(50000, 500000000), 
                     labels = as.character(paste0(c(0.05, seq(100, 500, 100)),"M")), 
                     breaks = c(50000, seq(100000000, 500000000, 100000000))) +
  labs(title = "", x = "Brazilian Reais (Millions)", y = "Proportion of Total Area") +
  theme_minimal()+
  theme(text = element_text(size = 16, family = "sans"),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "bottom"
  )



g6 <- result_df %>%  
  ggplot(aes(x = Budget)) +
  geom_line(aes(y = Degradation_area_pp, color = "Avoid Degradation", linetype=Region), linewidth = 1) +
  #geom_smooth(aes(y = Proportion_Degradation, color = "Avoid Degradation"), linewidth = 2, method = "gam") +
  geom_line(aes(y = Deforestation_area_pp, color = "Avoid Deforestation", linetype=Region), linewidth = 1) +
  #geom_smooth(aes(y = Proportion_Deforestation, color = "Avoid Deforestation"), linewidth = 2, method = "gam") +
  geom_line(aes(y = Restoration_area_pp, color = "Restoration", linetype=Region), linewidth = 1) +
  #geom_smooth(aes(y = Proportion_Restoration, color = "Restoration"), linewidth = 2, method = "gam") +
  labs(x = "Budget (Brazilian Reais)", y = "Proportion of Area") +
  scale_color_manual(values = c("Avoid Deforestation" = "#294B29", "Restoration" = "#789461", "Avoid Degradation" = "#76453B")) +
  facet_wrap(~benefit, ncol=1) +
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  scale_x_continuous(limits = c(50000, 500000000), 
                     labels = as.character(paste0(c(0.05, seq(100, 500, 100)),"M")), 
                     breaks = c(50000, seq(100000000, 500000000, 100000000))) +
  labs(title = "", x = "Brazilian Reais (Millions)", y = "Proportion of Area by Action") +
  theme_minimal()+
  theme(text = element_text(size = 16, family = "sans"),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "bottom"
  )


ggarrange(g5, g6, ncol = 2, common.legend = T, legend = "bottom")


rm(list=ls()[ls() %in% c("pgm.area.change.df.principals", "stm.area.change.df.principals",
                         "pgm.carbon.benefit.df.principals", "stm.carbon.benefit.df.principals", "carbon.benefit.principals",
                         "pgm.costs.df.principals", "stm.costs.df.principals", "costs.principals")])
gc()
#
#




# comparing benefits for combinations and for primary forests only =============

pgm.area.change.df.all <- as.data.frame(pgm.area.change[[3:8]], xy = TRUE)
pgm.area.change.df.all <- pgm.area.change.df.all %>% 
  pivot_longer(
    area_change_2020_avoiddeforest:area_change_2020_restor_n_avoidboth,
    names_to = "Scenario",
    values_to = "Cat"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "PGM",
    Scenario = factor(Scenario,
                      levels = c("area_change_2020_avoiddeforest",
                                 "area_change_2020_avoiddegrad",
                                 "area_change_2020_restor_wo_avoid",
                                 "area_change_2020_avoidboth",
                                 "area_change_2020_restor_n_avoiddeforest",
                                 "area_change_2020_restor_n_avoidboth"),
                      labels = c("Avoid Deforestation", "Avoid Degradation",
                                 "Restoration without avoid", "Avoid both",
                                 "Restoration and avoid deforestation",
                                 "Restoration and avoid both")),
    #Cat = na_if(Cat, 0),
    Cat = factor(Cat,
                 levels = c(1,10,25,100,125),
                 labels = c("UPF", "DPF", "RDPF", "SF", "DSF"))
  ) %>% ungroup() %>% drop_na()



stm.area.change.df.all <- as.data.frame(stm.area.change[[3:8]], xy = TRUE)
stm.area.change.df.all <- stm.area.change.df.all %>% 
  pivot_longer(
    area_change_2020_avoiddeforest:area_change_2020_restor_n_avoidboth,
    names_to = "Scenario",
    values_to = "Cat"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "STM",
    Scenario = factor(Scenario,
                      levels = c("area_change_2020_avoiddeforest",
                                 "area_change_2020_avoiddegrad",
                                 "area_change_2020_restor_wo_avoid",
                                 "area_change_2020_avoidboth",
                                 "area_change_2020_restor_n_avoiddeforest",
                                 "area_change_2020_restor_n_avoidboth"),
                      labels = c("Avoid Deforestation", "Avoid Degradation",
                                 "Restoration without avoid", "Avoid both",
                                 "Restoration and avoid deforestation",
                                 "Restoration and avoid both")),
    #Cat = na_if(Cat, 0),
    Cat = factor(Cat,
                 levels = c(1,10,25,100,125),
                 labels = c("UPF", "DPF", "RDPF", "SF", "DSF"))
  ) %>% ungroup() %>% drop_na()



pgm.biodiversity.benefit.df.all <- as.data.frame(pgm.biodiversity.benefit[[1:6]], xy = TRUE)
pgm.biodiversity.benefit.df.all <- pgm.biodiversity.benefit.df.all %>% 
  pivot_longer(
    PGM_2020_avoiddeforest_biodiversity_benefit:PGM_2020_restor_n_avoidboth_biodiversity_benefit,
    names_to = "Scenario",
    values_to = "Benefit"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "PGM",
    Scenario = factor(Scenario,
                      levels = c("PGM_2020_avoiddeforest_biodiversity_benefit",
                                 "PGM_2020_avoiddegrad_biodiversity_benefit",
                                 "PGM_2020_restor_wo_avoid_biodiversity_benefit",
                                 "PGM_2020_avoidboth_biodiversity_benefit",
                                 "PGM_2020_restor_n_avoiddeforest_biodiversity_benefit",
                                 "PGM_2020_restor_n_avoidboth_biodiversity_benefit"),
                      labels = c("Avoid Deforestation", "Avoid Degradation",
                                 "Restoration without avoid", "Avoid both",
                                 "Restoration and avoid deforestation",
                                 "Restoration and avoid both")),
    #Benefit = na_if(Benefit, 0)
  ) %>% ungroup() %>% drop_na() %>% 
  left_join(pgm.area.change.df.all)
#head(pgm.biodiversity.benefit.df.all)


stm.biodiversity.benefit.df.all <- as.data.frame(stm.biodiversity.benefit[[1:6]], xy = TRUE)
stm.biodiversity.benefit.df.all <- stm.biodiversity.benefit.df.all %>% 
  pivot_longer(
    STM_2020_avoiddeforest_biodiversity_benefit:STM_2020_restor_n_avoidboth_biodiversity_benefit,
    names_to = "Scenario",
    values_to = "Benefit"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "STM",
    Scenario = factor(Scenario,
                      levels = c("STM_2020_avoiddeforest_biodiversity_benefit",
                                 "STM_2020_avoiddegrad_biodiversity_benefit",
                                 "STM_2020_restor_wo_avoid_biodiversity_benefit",
                                 "STM_2020_avoidboth_biodiversity_benefit",
                                 "STM_2020_restor_n_avoiddeforest_biodiversity_benefit",
                                 "STM_2020_restor_n_avoidboth_biodiversity_benefit"),
                      labels = c("Avoid Deforestation", "Avoid Degradation",
                                 "Restoration without avoid", "Avoid both",
                                 "Restoration and avoid deforestation",
                                 "Restoration and avoid both")),
    #Benefit = na_if(Benefit, 0)
  ) %>% ungroup() %>% drop_na() %>% 
  left_join(stm.area.change.df.all)
#head(stm.biodiversity.benefit.df.all)


#


biodiversity.benefit.all <- rbind(pgm.biodiversity.benefit.df.all, stm.biodiversity.benefit.df.all)
biodiversity.benefit.all <- biodiversity.benefit.all %>% mutate(Region = factor(Region, levels = c("PGM", "STM")))




##violin plot + jitter [dots colored by LULC category]
biodiversity.benefit.all %>% filter(Benefit>=0) %>%
  ggplot(aes(x=Scenario, y=Benefit)) +
  geom_violin(scale="width") + 
  geom_sina(data=biodiversity.benefit.all %>% sample_n(100000) %>% filter(Benefit>=0),
            aes(color=Cat, group=Scenario), alpha=.1, scale="width")+
  scale_x_discrete(labels=addline_format(levels(biodiversity.benefit.all$Scenario)),
                   expand = c(.05, .05)) +
  scale_color_manual(values = c("#294B29", "#50623A", "#76453B", "#789461", "#B19470")) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  labs(title = "", x = "", y = "Biodiversity benefit") +
  theme_minimal()+
  theme(text = element_text(size = 16, family = "sans"),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "bottom")




#
#



pgm.carbon.benefit.df.all <- as.data.frame(pgm.carbon.benefit[[1:6]], xy = TRUE)
pgm.carbon.benefit.df.all <- pgm.carbon.benefit.df.all %>% 
  pivot_longer(
    PGM_2020_avoiddeforest_carbon_benefit:PGM_2020_restor_n_avoidboth_carbon_benefit,
    names_to = "Scenario",
    values_to = "Benefit"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "PGM",
    Scenario = factor(Scenario,
                      levels = c("PGM_2020_avoiddeforest_carbon_benefit",
                                 "PGM_2020_avoiddegrad_carbon_benefit",
                                 "PGM_2020_restor_wo_avoid_carbon_benefit",
                                 "PGM_2020_avoidboth_carbon_benefit",
                                 "PGM_2020_restor_n_avoiddeforest_carbon_benefit",
                                 "PGM_2020_restor_n_avoidboth_carbon_benefit"),
                      labels = c("Avoid Deforestation", "Avoid Degradation",
                                 "Restoration without avoid", "Avoid both",
                                 "Restoration and avoid deforestation",
                                 "Restoration and avoid both")),
    #Benefit = na_if(Benefit, 0)
  ) %>% ungroup() %>% drop_na() %>% 
  left_join(pgm.area.change.df.all)
#head(pgm.carbon.benefit.df.all)


stm.carbon.benefit.df.all <- as.data.frame(stm.carbon.benefit[[1:6]], xy = TRUE)
stm.carbon.benefit.df.all <- stm.carbon.benefit.df.all %>% 
  pivot_longer(
    STM_2020_avoiddeforest_carbon_benefit:STM_2020_restor_n_avoidboth_carbon_benefit,
    names_to = "Scenario",
    values_to = "Benefit"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "STM",
    Scenario = factor(Scenario,
                      levels = c("STM_2020_avoiddeforest_carbon_benefit",
                                 "STM_2020_avoiddegrad_carbon_benefit",
                                 "STM_2020_restor_wo_avoid_carbon_benefit",
                                 "STM_2020_avoidboth_carbon_benefit",
                                 "STM_2020_restor_n_avoiddeforest_carbon_benefit",
                                 "STM_2020_restor_n_avoidboth_carbon_benefit"),
                      labels = c("Avoid Deforestation", "Avoid Degradation",
                                 "Restoration without avoid", "Avoid both",
                                 "Restoration and avoid deforestation",
                                 "Restoration and avoid both")),
    #Benefit = na_if(Benefit, 0)
  ) %>% ungroup() %>% drop_na() %>% 
  left_join(stm.area.change.df.all)
#head(stm.carbon.benefit.df.all)


#


carbon.benefit.all <- rbind(pgm.carbon.benefit.df.all, stm.carbon.benefit.df.all)
carbon.benefit.all <- carbon.benefit.all %>% mutate(Region = factor(Region, levels = c("PGM", "STM")))




##violin plot + jitter [dots colored by LULC category]
carbon.benefit.all %>% filter(Benefit>=0) %>%
  ggplot(aes(x=Scenario, y=Benefit)) +
  geom_violin(scale="width") + 
  geom_sina(data=carbon.benefit.all %>% sample_n(100000) %>% filter(Benefit>=0),
            aes(color=Cat, group=Scenario), alpha=.1, scale="width")+
  scale_x_discrete(labels=addline_format(levels(carbon.benefit.all$Scenario)),
                   expand = c(.05, .05)) +
  scale_color_manual(values = c("#294B29", "#50623A", "#76453B", "#789461", "#B19470")) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  labs(title = "", x = "", y = "Carbon benefit (MgC/ha)") +
  theme_minimal()+
  theme(text = element_text(size = 16, family = "sans"),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "bottom")




#
#



pgm.costs.df.all <- as.data.frame(pgm.costs[[1:6]], xy = TRUE)
pgm.costs.df.all <- pgm.costs.df.all %>% 
  pivot_longer(
    PGM_2020_avoiddeforest_costs:PGM_2020_restor_n_avoidboth_costs,
    names_to = "Scenario",
    values_to = "Costs"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "PGM",
    Scenario = factor(Scenario,
                      levels = c("PGM_2020_avoiddeforest_costs",
                                 "PGM_2020_avoiddegrad_costs",
                                 "PGM_2020_restor_wo_avoid_costs",
                                 "PGM_2020_avoidboth_costs",
                                 "PGM_2020_restor_n_avoiddeforest_costs",
                                 "PGM_2020_restor_n_avoidboth_costs"),
                      labels = c("Avoid Deforestation", "Avoid Degradation",
                                 "Restoration without avoid", "Avoid both",
                                 "Restoration and avoid deforestation",
                                 "Restoration and avoid both")),
    #Benefit = na_if(Benefit, 0)
  ) %>% ungroup() %>% drop_na() %>% 
  left_join(pgm.area.change.df.all)
#head(pgm.costs.df.all)


stm.costs.df.all <- as.data.frame(stm.costs[[1:6]], xy = TRUE)
stm.costs.df.all <- stm.costs.df.all %>% 
  pivot_longer(
    STM_2020_avoiddeforest_costs:STM_2020_restor_n_avoidboth_costs,
    names_to = "Scenario",
    values_to = "Costs"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "STM",
    Scenario = factor(Scenario,
                      levels = c("STM_2020_avoiddeforest_costs",
                                 "STM_2020_avoiddegrad_costs",
                                 "STM_2020_restor_wo_avoid_costs",
                                 "STM_2020_avoidboth_costs",
                                 "STM_2020_restor_n_avoiddeforest_costs",
                                 "STM_2020_restor_n_avoidboth_costs"),
                      labels = c("Avoid Deforestation", "Avoid Degradation",
                                 "Restoration without avoid", "Avoid both",
                                 "Restoration and avoid deforestation",
                                 "Restoration and avoid both")),
    #Benefit = na_if(Benefit, 0)
  ) %>% ungroup() %>% drop_na() %>% 
  left_join(stm.area.change.df.all)
#head(stm.costs.df.all)


#


costs.all <- rbind(pgm.costs.df.all, stm.costs.df.all)
costs.all <- costs.all %>% mutate(Region = factor(Region, levels = c("PGM", "STM")))




##boxplot
costs.all %>% 
  ggplot(aes(x=Scenario, y=Costs)) +
  geom_boxplot(varwidth = T, outlier.shape = NA, fill="#FF4500", color="#000000", show.legend = F) + 
  #geom_hline(yintercept=1, linetype='dashed', linewidth=1)+
  scale_x_discrete(labels=addline_format(levels(costs.all$Scenario)),
                   expand = c(.05, .05)) +
  scale_y_continuous(limits = c(0, 600)) +
  labs(title = "", x = "", y = "Costs (R$ / ha / year)") +
  theme_minimal()+
  theme(text = element_text(size = 16, family = "sans"),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "bottom")




costs.all %>% 
  left_join(biodiversity.benefit.all) %>% 
  mutate(BBenefit.year = Benefit/10,
         B.CBr = ifelse(BBenefit.year == 0, NA, Costs/BBenefit.year)) %>% 
  ggplot(aes(x=Scenario, y=B.CBr)) +
  geom_boxplot(varwidth = T, outlier.shape = NA, fill="#603a62", color="#A4A4A4", show.legend = F) + 
  #geom_hline(yintercept=1, linetype='dashed', linewidth=1)+
  scale_x_discrete(labels=addline_format(levels(costs.all$Scenario)),
                   expand = c(.05, .05)) +
  scale_y_continuous(limits = c(0, 400)) +
  labs(title = "", x = "", y = "Cost-Benefit ratio (R$ / Biodiversity)") +
  theme_minimal()+
  theme(text = element_text(size = 16, family = "sans"),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "bottom")




costs.all %>% 
  left_join(carbon.benefit.all) %>% 
  mutate(CBenefit.year = Benefit/10,
         C.CBr = ifelse(CBenefit.year == 0, NA, Costs/CBenefit.year)) %>% 
  ggplot(aes(x=Scenario, y=C.CBr)) +
  geom_boxplot(varwidth = T, outlier.shape = NA, fill="#603a62", color="#A4A4A4", show.legend = F) + 
  #geom_hline(yintercept=1, linetype='dashed', linewidth=1)+
  scale_x_discrete(labels=addline_format(levels(costs.all$Scenario)),
                   expand = c(.05, .05)) +
  scale_y_continuous(limits = c(0, 400)) +
  labs(title = "", x = "", y = "Cost-Benefit ratio (R$ / MgC)") +
  theme_minimal()+
  theme(text = element_text(size = 16, family = "sans"),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "bottom")





#
#



pgm.area.change.df.primary <- as.data.frame(pgm.area.change[[9:13]], xy = TRUE)
pgm.area.change.df.primary <- pgm.area.change.df.primary %>% 
  pivot_longer(
    area_change_2020_avoiddeforest2:area_change_2020_restor_n_avoidboth2,
    names_to = "Scenario",
    values_to = "Cat"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "PGM",
    Scenario = factor(Scenario,
                      levels = c("area_change_2020_avoiddeforest2",
                                 "area_change_2020_avoiddegrad2",
                                 "area_change_2020_avoidboth2",
                                 "area_change_2020_restor_n_avoiddeforest2",
                                 "area_change_2020_restor_n_avoidboth2"),
                      labels = c("Avoid Deforestation", "Avoid Degradation",
                                 "Avoid both",
                                 "Restoration and avoid deforestation",
                                 "Restoration and avoid both")),
    #Cat = na_if(Cat, 0),
    Cat = factor(Cat,
                 levels = c(1,10,25,100,125),
                 labels = c("UPF", "DPF", "RDPF", "SF", "DSF"))
  ) %>% ungroup() %>% drop_na()



stm.area.change.df.primary <- as.data.frame(stm.area.change[[9:13]], xy = TRUE)
stm.area.change.df.primary <- stm.area.change.df.primary %>% 
  pivot_longer(
    area_change_2020_avoiddeforest2:area_change_2020_restor_n_avoidboth2,
    names_to = "Scenario",
    values_to = "Cat"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "STM",
    Scenario = factor(Scenario,
                      levels = c("area_change_2020_avoiddeforest2",
                                 "area_change_2020_avoiddegrad2",
                                 "area_change_2020_avoidboth2",
                                 "area_change_2020_restor_n_avoiddeforest2",
                                 "area_change_2020_restor_n_avoidboth2"),
                      labels = c("Avoid Deforestation", "Avoid Degradation",
                                 "Avoid both",
                                 "Restoration and avoid deforestation",
                                 "Restoration and avoid both")),
    #Cat = na_if(Cat, 0),
    Cat = factor(Cat,
                 levels = c(1,10,25,100,125),
                 labels = c("UPF", "DPF", "RDPF", "SF", "DSF"))
  ) %>% ungroup() %>% drop_na()



pgm.biodiversity.benefit.df.primary <- as.data.frame(pgm.biodiversity.benefit[[7:11]], xy = TRUE)
pgm.biodiversity.benefit.df.primary <- pgm.biodiversity.benefit.df.primary %>% 
  pivot_longer(
    PGM_2020_avoiddeforest2_biodiversity_benefit:PGM_2020_restor_n_avoidboth2_biodiversity_benefit,
    names_to = "Scenario",
    values_to = "Benefit"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "PGM",
    Scenario = factor(Scenario,
                      levels = c("PGM_2020_avoiddeforest2_biodiversity_benefit",
                                 "PGM_2020_avoiddegrad2_biodiversity_benefit",
                                 "PGM_2020_avoidboth2_biodiversity_benefit",
                                 "PGM_2020_restor_n_avoiddeforest2_biodiversity_benefit",
                                 "PGM_2020_restor_n_avoidboth2_biodiversity_benefit"),
                      labels = c("Avoid Deforestation", "Avoid Degradation",
                                 "Avoid both",
                                 "Restoration and avoid deforestation",
                                 "Restoration and avoid both")),
    #Benefit = na_if(Benefit, 0)
  ) %>% ungroup() %>% drop_na() %>% 
  left_join(pgm.area.change.df.primary)
#head(pgm.biodiversity.benefit.df.primary)


stm.biodiversity.benefit.df.primary <- as.data.frame(stm.biodiversity.benefit[[7:11]], xy = TRUE)
stm.biodiversity.benefit.df.primary <- stm.biodiversity.benefit.df.primary %>% 
  pivot_longer(
    STM_2020_avoiddeforest2_biodiversity_benefit:STM_2020_restor_n_avoidboth2_biodiversity_benefit,
    names_to = "Scenario",
    values_to = "Benefit"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "STM",
    Scenario = factor(Scenario,
                      levels = c("STM_2020_avoiddeforest2_biodiversity_benefit",
                                 "STM_2020_avoiddegrad2_biodiversity_benefit",
                                 "STM_2020_avoidboth_biodiversity2_benefit",
                                 "STM_2020_restor_n_avoiddeforest2_biodiversity_benefit",
                                 "STM_2020_restor_n_avoidboth2_biodiversity_benefit"),
                      labels = c("Avoid Deforestation", "Avoid Degradation",
                                 "Avoid both",
                                 "Restoration and avoid deforestation",
                                 "Restoration and avoid both")),
    #Benefit = na_if(Benefit, 0)
  ) %>% ungroup() %>% drop_na() %>% 
  left_join(stm.area.change.df.primary)
#head(stm.biodiversity.benefit.df.primary)


#


biodiversity.benefit.primary <- rbind(pgm.biodiversity.benefit.df.primary, stm.biodiversity.benefit.df.primary)
biodiversity.benefit.primary <- biodiversity.benefit.primary %>% mutate(Region = factor(Region, levels = c("PGM", "STM")))




##violin plot + jitter [dots colored by LULC category]
biodiversity.benefit.primary %>% filter(Benefit>=0) %>%
  ggplot(aes(x=Scenario, y=Benefit)) +
  geom_violin(scale="width") + 
  geom_sina(data=biodiversity.benefit.primary %>% sample_n(100000) %>% filter(Benefit>=0),
            aes(color=Cat, group=Scenario), alpha=.1, scale="width")+
  scale_x_discrete(labels=addline_format(levels(biodiversity.benefit.primary$Scenario)),
                   expand = c(.05, .05)) +
  scale_color_manual(values = c("#294B29", "#50623A", "#76453B", "#789461", "#B19470")) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  labs(title = "", x = "", y = "Biodiversity benefit") +
  theme_minimal()+
  theme(text = element_text(size = 16, family = "sans"),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "bottom")




#
#



pgm.carbon.benefit.df.primary <- as.data.frame(pgm.carbon.benefit[[7:11]], xy = TRUE)
pgm.carbon.benefit.df.primary <- pgm.carbon.benefit.df.primary %>% 
  pivot_longer(
    PGM_2020_avoiddeforest2_carbon_benefit:PGM_2020_restor_n_avoidboth2_carbon_benefit,
    names_to = "Scenario",
    values_to = "Benefit"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "PGM",
    Scenario = factor(Scenario,
                      levels = c("PGM_2020_avoiddeforest2_carbon_benefit",
                                 "PGM_2020_avoiddegrad2_carbon_benefit",
                                 "PGM_2020_avoidboth2_carbon_benefit",
                                 "PGM_2020_restor_n_avoiddeforest2_carbon_benefit",
                                 "PGM_2020_restor_n_avoidboth2_carbon_benefit"),
                      labels = c("Avoid Deforestation", "Avoid Degradation",
                                 "Avoid both",
                                 "Restoration and avoid deforestation",
                                 "Restoration and avoid both")),
    #Benefit = na_if(Benefit, 0)
  ) %>% ungroup() %>% drop_na() %>% 
  left_join(pgm.area.change.df.primary)
#head(pgm.carbon.benefit.df.primary)


stm.carbon.benefit.df.primary <- as.data.frame(stm.carbon.benefit[[7:11]], xy = TRUE)
stm.carbon.benefit.df.primary <- stm.carbon.benefit.df.primary %>% 
  pivot_longer(
    STM_2020_avoiddeforest2_carbon_benefit:STM_2020_restor_n_avoidboth2_carbon_benefit,
    names_to = "Scenario",
    values_to = "Benefit"
  ) %>% 
  group_by(Scenario) %>% 
  mutate(
    Cell = row_number(),
    Region = "STM",
    Scenario = factor(Scenario,
                      levels = c("STM_2020_avoiddeforest2_carbon_benefit",
                                 "STM_2020_avoiddegrad2_carbon_benefit",
                                 "STM_2020_avoidboth_carbon2_benefit",
                                 "STM_2020_restor_n_avoiddeforest2_carbon_benefit",
                                 "STM_2020_restor_n_avoidboth2_carbon_benefit"),
                      labels = c("Avoid Deforestation", "Avoid Degradation",
                                 "Avoid both",
                                 "Restoration and avoid deforestation",
                                 "Restoration and avoid both")),
    #Benefit = na_if(Benefit, 0)
  ) %>% ungroup() %>% drop_na() %>% 
  left_join(stm.area.change.df.primary)
#head(stm.carbon.benefit.df.primary)


#


carbon.benefit.primary <- rbind(pgm.carbon.benefit.df.primary, stm.carbon.benefit.df.primary)
carbon.benefit.primary <- carbon.benefit.primary %>% mutate(Region = factor(Region, levels = c("PGM", "STM")))




##violin plot + jitter [dots colored by LULC category]
carbon.benefit.primary %>% filter(Benefit>=0) %>%
  ggplot(aes(x=Scenario, y=Benefit)) +
  geom_violin(scale="width") + 
  geom_sina(data=carbon.benefit.primary %>% sample_n(100000) %>% filter(Benefit>=0),
            aes(color=Cat, group=Scenario), alpha=.1, scale="width")+
  scale_x_discrete(labels=addline_format(levels(carbon.benefit.primary$Scenario)),
                   expand = c(.05, .05)) +
  scale_color_manual(values = c("#294B29", "#50623A", "#76453B", "#789461", "#B19470")) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  labs(title = "", x = "", y = "Carbon benefit (MgC/ha)") +
  theme_minimal()+
  theme(text = element_text(size = 16, family = "sans"),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "bottom")




#
#
























































#==============================| previous modeling approach





#filter the cells with conservation action for each scenario
# [the difference in total forests between scenario and 2020 Real]

pgm.allforest.list <- list.files("rasters/PGM/all_forest_mask/", pattern = ".tif", full.names = T, recursive = T)

pgm.allforest <- stack(pgm.allforest.list)
#names(pgm.allforest) <- unlist(strsplit(pgm.allforest.list, "/|.tif"))[seq(4,32,4)]


pgm.conservationaction.avoiddegrad <- pgm.allforest[["PGM_2020_avoiddegrad"]] - pgm.allforest[["PGM_2020_real"]]
pgm.conservationaction.avoiddegrad[pgm.conservationaction.avoiddegrad!=1]<-NA

pgm.conservationaction.avoiddeforest <- pgm.allforest[["PGM_2020_avoiddeforest"]] - pgm.allforest[["PGM_2020_real"]]
pgm.conservationaction.avoiddeforest[pgm.conservationaction.avoiddeforest!=1]<-NA

pgm.conservationaction.avoiddeforest2 <- pgm.allforest[["PGM_2020_avoiddeforest2"]] - pgm.allforest[["PGM_2020_real"]]
pgm.conservationaction.avoiddeforest2[pgm.conservationaction.avoiddeforest2!=1]<-NA

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
#names(stm.allforest) <- unlist(strsplit(stm.allforest.list, "/|.tif"))[seq(4,32,4)]


stm.conservationaction.avoiddegrad <- stm.allforest[["STM_2020_avoiddegrad"]] - stm.allforest[["STM_2020_real"]]
stm.conservationaction.avoiddegrad[stm.conservationaction.avoiddegrad!=1]<-NA

stm.conservationaction.avoiddeforest <- stm.allforest[["STM_2020_avoiddeforest"]] - stm.allforest[["STM_2020_real"]]
stm.conservationaction.avoiddeforest[stm.conservationaction.avoiddeforest!=1]<-NA

stm.conservationaction.avoiddeforest2 <- stm.allforest[["STM_2020_avoiddeforest2"]] - stm.allforest[["STM_2020_real"]]
stm.conservationaction.avoiddeforest2[stm.conservationaction.avoiddeforest2!=1]<-NA

stm.conservationaction.avoidboth <- stm.allforest[["STM_2020_avoidboth"]] - stm.allforest[["STM_2020_real"]]
stm.conservationaction.avoidboth[stm.conservationaction.avoidboth!=1]<-NA

stm.conservationaction.restor_wo_avoid <- stm.allforest[["STM_2020_restor_wo_avoid"]] - stm.allforest[["STM_2020_real"]]
stm.conservationaction.restor_wo_avoid[stm.conservationaction.restor_wo_avoid!=1]<-NA

stm.conservationaction.restor_n_avoid_deforest <- stm.allforest[["STM_2020_restor_n_avoid_deforest"]] - stm.allforest[["STM_2020_real"]]
stm.conservationaction.restor_n_avoid_deforest[stm.conservationaction.restor_n_avoid_deforest!=1]<-NA

stm.conservationaction.restor_n_avoid_both <- stm.allforest[["STM_2020_restor_n_avoid_both"]] - stm.allforest[["STM_2020_real"]]
stm.conservationaction.restor_n_avoid_both[stm.conservationaction.restor_n_avoid_both!=1]<-NA



## STM mask for big rivers (Tapajos & Amazonas)
stm.mask <- raster("rasters/STM/raw/mapbiomas-brazil-collection-70-stm-2010-100mpx.tif")
stm.mask[stm.mask==0]<-NA
stm.mask[stm.mask==11]<-NA # 11 == Wetland
stm.mask[stm.mask==33]<-NA # 33 == River

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
stm.biodiversity.benefit <- mask(stm.biodiversity.benefit, stm.mask)
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
  ggplot(aes(x=Scenario, y=Biodiversity, fill = Scenario)) +
  geom_violin(show.legend = F) + 
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange", width = 0.2, show.legend = F) +
  scale_fill_manual(values = c("#41644A", "#41644A", "#41644A", 
                               "#D3756B", "#D3756B", "#D3756B")) +
  coord_flip() +
  facet_wrap(~Region) +
  labs(title = "", x = "", y = "Biodiversity benefit values") +
  theme_minimal()+
  theme(text = element_text(size = 16, family = "sans"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 14))


#########################################
####    comparing carbon benefits    ####
####      between scenarios and      ####
####         real projections        ####
####           (landscape)           ####
#########################################
carbon.benefit.list <- list.files("models.output/carbon.benefits/", pattern = ".tif", full.names = T, recursive = T)

pgm.carbon.benefit <- stack(grep("PGM", carbon.benefit.list, value = T))
#names(pgm.carbon.benefit) <- unlist(strsplit(carbon.benefit.list, "/|.tif"))[seq(3,21,3)]
#plot(pgm.carbon.benefit, nr=2, col = terrain.colors(length(seq(0, 225, by = 25)), rev = T), breaks= seq(0, 225, by = 25)) ## res = 1673 x 881
pgm.carbon.benefit.df <- as.data.frame(pgm.carbon.benefit, xy = TRUE)
pgm.carbon.benefit.df <- pgm.carbon.benefit.df %>% 
                            pivot_longer(
                              PGM_2010_real_carbon_benefit:PGM_2020_restor_wo_avoid_carbon_benefit,
                              names_to = "ID",
                              values_to = "Carbon"
                              ) %>% 
                            mutate(
                              #across(c('x', 'y'), round, 6),
                              Region = "PGM",
                              Scenario = factor(case_when(str_detect(ID, "avoiddegrad")~ "Avoid degradation",
                                                          str_detect(ID, "avoiddeforest2")~ "Avoid deforestation primary forest",
                                                          str_detect(ID, "avoiddeforest")~ "Avoid deforestation",
                                                          str_detect(ID, "avoidboth")~ "Avoid both",
                                                          str_detect(ID, "restor_wo_avoid")~ "Restoration without avoid",
                                                          str_detect(ID, "restor_n_avoid_deforest")~ "Restoration and avoid deforestation",
                                                          str_detect(ID, "restor_n_avoid_both")~ "Restoration and avoid both",
                                                          str_detect(ID, "2010_real")~ "2010 Real",
                                                          str_detect(ID, "2020_real")~ "2020 Real"),
                                                levels = c("2010 Real",
                                                           "2020 Real",
                                                           "Avoid degradation",
                                                           "Avoid deforestation",
                                                           "Avoid deforestation primary forest",
                                                           "Avoid both",
                                                           "Restoration without avoid",
                                                           "Restoration and avoid deforestation",
                                                           "Restoration and avoid both")),
                              Carbon = na_if(Carbon, 0)
                              )


stm.carbon.benefit <- stack(grep("STM", carbon.benefit.list, value = T))
#names(stm.carbon.benefit) <- unlist(strsplit(carbon.benefit.list, "/|.tif"))[seq(24,42,3)]
stm.carbon.benefit <- mask(stm.carbon.benefit, stm.mask)
#plot(stm.carbon.benefit, nr=2, col = terrain.colors(length(seq(0, 225, by = 25)), rev = T), breaks= seq(0, 225, by = 25)) ##res = 1673 x 881
stm.carbon.benefit.df <- as.data.frame(stm.carbon.benefit, xy = TRUE)
stm.carbon.benefit.df <- stm.carbon.benefit.df %>% 
                               pivot_longer(
                                 STM_2010_real_carbon_benefit:STM_2020_restor_wo_avoid_carbon_benefit,
                                 names_to = "ID",
                                 values_to = "Carbon"
                                 ) %>% 
                               mutate(
                                 #across(c('x', 'y'), round, 6),
                                 Region = "STM",
                                 Scenario = factor(case_when(str_detect(ID, "avoiddegrad")~ "Avoid degradation",
                                                             str_detect(ID, "avoiddeforest2")~ "Avoid deforestation primary forest",
                                                             str_detect(ID, "avoiddeforest")~ "Avoid deforestation",
                                                             str_detect(ID, "avoidboth")~ "Avoid both",
                                                             str_detect(ID, "restor_wo_avoid")~ "Restoration without avoid",
                                                             str_detect(ID, "restor_n_avoid_deforest")~ "Restoration and avoid deforestation",
                                                             str_detect(ID, "restor_n_avoid_both")~ "Restoration and avoid both",
                                                             str_detect(ID, "2010_real")~ "2010 Real",
                                                             str_detect(ID, "2020_real")~ "2020 Real"),
                                                   levels = c("2010 Real",
                                                              "2020 Real",
                                                              "Avoid degradation",
                                                              "Avoid deforestation",
                                                              "Avoid deforestation primary forest",
                                                              "Avoid both",
                                                              "Restoration without avoid",
                                                              "Restoration and avoid deforestation",
                                                              "Restoration and avoid both")),
                                 Carbon = na_if(Carbon, 0)
                               )


carbon.benefit <- rbind(pgm.carbon.benefit.df, stm.carbon.benefit.df)

cowplot::plot_grid(
  
  carbon.benefit %>% filter(ID == "PGM_2020_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Real", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = c(.1,.8)),
  
  carbon.benefit %>% filter(ID == "PGM_2020_avoiddegrad_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid degradation", x = "", y = "Latitude") +
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
  
  carbon.benefit %>% filter(ID == "PGM_2020_avoiddeforest2_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid deforestation primary forest", x = "", y = "Latitude") +
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
  
  carbon.benefit %>% filter(ID == "PGM_2020_restor_wo_avoid_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration without avoid", x = "", y = "") +
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
  
  carbon.benefit %>% filter(ID == "PGM_2020_restor_n_avoid_both_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>% 
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration and avoid both", x = "Longitude", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  ncol = 4, align = "hv")



cowplot::plot_grid(
  
  carbon.benefit %>% filter(ID == "STM_2020_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Real", x = "", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = c(.1,.8)),
  
  carbon.benefit %>% filter(ID == "STM_2020_avoiddegrad_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid degradation", x = "", y = "Latitude") +
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
  
  carbon.benefit %>% filter(ID == "STM_2020_avoiddeforest2_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Avoid deforestation primary forest", x = "", y = "Latitude") +
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
  
  carbon.benefit %>% filter(ID == "STM_2020_restor_wo_avoid_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>% 
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration without avoid", x = "", y = "") +
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
  
  carbon.benefit %>% filter(ID == "STM_2020_restor_n_avoid_both_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>% 
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = Carbon)) +
    scale_fill_gradient2(low = "#F3F6F4", mid = "#8fce00", high = "#374f00", na.value = NA) +
    labs(title = "Restoration and avoid both", x = "Longitude", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  ncol = 4, align = "hv")





carbon.benefit %>% group_by(Region, Scenario) %>% summarise(mean.benefit = mean(Carbon, na.rm=T))

cowplot::plot_grid(
  
  carbon.benefit %>% filter(ID == "PGM_2020_real_carbon_benefit" | ID == "PGM_2010_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("gray", "#374f00")) +
    scale_fill_manual(values = c("gray", "#374f00")) +
    geom_hline(aes(yintercept=83.4), linewidth=1, color = "gray") +
    geom_hline(aes(yintercept=76.9), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "Real", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(ID == "PGM_2020_avoiddegrad_carbon_benefit" | ID == "PGM_2010_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
  ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("gray", "#374f00")) +
    scale_fill_manual(values = c("gray", "#374f00")) +
    geom_hline(aes(yintercept=83.4), linewidth=1, color = "gray") +
    geom_hline(aes(yintercept=87.5), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "Avoid degradation", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_avoiddeforest_carbon_benefit" | ID == "PGM_2010_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
  ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("grey", "#374f00")) +
    scale_fill_manual(values = c("grey", "#374f00")) +
    geom_hline(aes(yintercept=83.4), linewidth=1, color = "grey") +
    geom_hline(aes(yintercept=79.1), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "Avoid deforestation", x = "", y = "Carbon benefit values") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "PGM" & ID == "PGM_2020_avoiddeforest2_carbon_benefit" | ID == "PGM_2010_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("grey", "#374f00")) +
    scale_fill_manual(values = c("grey", "#374f00")) +
    geom_hline(aes(yintercept=83.4), linewidth=1, color = "grey") +
    geom_hline(aes(yintercept=82.4), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "Avoid deforestation primary forest", x = "", y = "Carbon benefit values") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(ID == "PGM_2020_avoidboth_carbon_benefit" | ID == "PGM_2010_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
  ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("grey", "#374f00")) +
    scale_fill_manual(values = c("grey", "#374f00")) +
    geom_hline(aes(yintercept=83.4), linewidth=1, color = "grey") +
    geom_hline(aes(yintercept=90.3), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "Avoid both", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(ID == "PGM_2020_restor_wo_avoid_carbon_benefit" | ID == "PGM_2010_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("grey", "#374f00")) +
    scale_fill_manual(values = c("grey", "#374f00")) +
    geom_hline(aes(yintercept=83.4), linewidth=1, color = "grey") +
    geom_hline(aes(yintercept=79), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "Restoration without avoid", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(ID == "PGM_2020_restor_n_avoid_deforest_carbon_benefit" | ID == "PGM_2010_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("grey", "#374f00")) +
    scale_fill_manual(values = c("grey", "#374f00")) +
    geom_hline(aes(yintercept=83.4), linewidth=1, color = "grey") +
    geom_hline(aes(yintercept=80.9), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "Restoration and avoid deforestation", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(ID == "PGM_2020_restor_n_avoid_both_carbon_benefit" | ID == "PGM_2010_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
  ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("grey", "#374f00"), labels = c("Real", "Scenario")) +
    scale_fill_manual(values = c("grey", "#374f00"), labels = c("Real", "Scenario")) +
    geom_hline(aes(yintercept=83.4), linewidth=1, color = "grey") +
    geom_hline(aes(yintercept=90.4), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "Restoration and avoid both", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = c(.8,.8)),
  
  ncol = 4, align = "hv")




cowplot::plot_grid(
  
  carbon.benefit %>% filter(ID == "STM_2020_real_carbon_benefit" | ID == "STM_2010_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("gray", "#374f00")) +
    scale_fill_manual(values = c("gray", "#374f00")) +
    geom_hline(aes(yintercept=94.6), linewidth=1, color = "gray") +
    geom_hline(aes(yintercept=89.4), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "Real", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(ID == "STM_2020_avoiddegrad_carbon_benefit" | ID == "STM_2010_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("gray", "#374f00")) +
    scale_fill_manual(values = c("gray", "#374f00")) +
    geom_hline(aes(yintercept=94.6), linewidth=1, color = "gray") +
    geom_hline(aes(yintercept=97.9), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "Avoid degradation", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "STM" & ID == "STM_2020_avoiddeforest_carbon_benefit" | ID == "STM_2010_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("grey", "#374f00")) +
    scale_fill_manual(values = c("grey", "#374f00")) +
    geom_hline(aes(yintercept=94.6), linewidth=1, color = "grey") +
    geom_hline(aes(yintercept=92.1), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "Avoid deforestation", x = "", y = "Carbon benefit values") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(Region == "STM" & ID == "STM_2020_avoiddeforest2_carbon_benefit" | ID == "STM_2010_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("grey", "#374f00")) +
    scale_fill_manual(values = c("grey", "#374f00")) +
    geom_hline(aes(yintercept=94.6), linewidth=1, color = "grey") +
    geom_hline(aes(yintercept=94.3), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "Avoid deforestation primary forest", x = "", y = "Carbon benefit values") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(ID == "STM_2020_avoidboth_carbon_benefit" | ID == "STM_2010_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("grey", "#374f00")) +
    scale_fill_manual(values = c("grey", "#374f00")) +
    geom_hline(aes(yintercept=94.6), linewidth=1, color = "grey") +
    geom_hline(aes(yintercept=101), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "Avoid both", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(ID == "STM_2020_restor_wo_avoid_carbon_benefit" | ID == "STM_2010_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("grey", "#374f00")) +
    scale_fill_manual(values = c("grey", "#374f00")) +
    geom_hline(aes(yintercept=94.6), linewidth=1, color = "grey") +
    geom_hline(aes(yintercept=90), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "Restoration without avoid", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(ID == "STM_2020_restor_n_avoid_deforest_carbon_benefit" | ID == "STM_2010_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("grey", "#374f00")) +
    scale_fill_manual(values = c("grey", "#374f00")) +
    geom_hline(aes(yintercept=94.6), linewidth=1, color = "grey") +
    geom_hline(aes(yintercept=94.1), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "Restoration and avoid deforestation", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  
  
  carbon.benefit %>% filter(ID == "STM_2020_restor_n_avoid_both_carbon_benefit" | ID == "STM_2010_real_carbon_benefit") %>% 
    #mutate(Carbon = na_if(Carbon, 0)) %>%
    ggplot(aes(y=Carbon, fill = Scenario, colour = Scenario)) +
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("grey", "#374f00"), labels = c("Real", "Scenario")) +
    scale_fill_manual(values = c("grey", "#374f00"), labels = c("Real", "Scenario")) +
    geom_hline(aes(yintercept=94.6), linewidth=1, color = "grey") +
    geom_hline(aes(yintercept=101), linewidth=1, linetype='dashed', color = "#374f00") +
    labs(title = "Restoration and avoid both", x = "", y = "") +
    theme_minimal() +
    theme(legend.position = c(.8,.8)),
  
  ncol = 4, align = "hv")





#carbon.benefit %>% 
#  group_by(ID) %>% drop_na() %>% 
#  summarise(d = list(density(Carbon, from = min(.$Carbon), to = max(.$Carbon)))) %>% 
#  do(rbind(data.frame(Region = "PGM",
#                      Scenario = "Real",
#                      x = .$d[[1]]$x,    
#                      y = .$d[[1]]$y - .$d[[6]]$y),
#           data.frame(Region = "PGM",
#                      Scenario = "Avoid degradation",
#                      x = .$d[[1]]$x,    
#                      y = .$d[[1]]$y - .$d[[5]]$y),
#           data.frame(Region = "PGM",
#                      Scenario = "Avoid deforestation",
#                      x = .$d[[1]]$x,    
#                      y = .$d[[1]]$y - .$d[[4]]$y),
#           data.frame(Region = "PGM",
#                      Scenario = "Avoid deforestation primary forest",
#                      x = .$d[[1]]$x,    
#                      y = .$d[[1]]$y - .$d[[3]]$y),
#           data.frame(Region = "PGM",
#                      Scenario = "Avoid both",
#                      x = .$d[[1]]$x,    
#                      y = .$d[[1]]$y - .$d[[2]]$y),
#           data.frame(Region = "PGM",
#                      Scenario = "Restoration without avoid",
#                      x = .$d[[1]]$x,    
#                      y = .$d[[1]]$y - .$d[[9]]$y),
#           data.frame(Region = "PGM",
#                      Scenario = "Restoration and avoid deforestation",
#                      x = .$d[[1]]$x,    
#                      y = .$d[[1]]$y - .$d[[8]]$y),
#           data.frame(Region = "PGM",
#                      Scenario = "Restoration and avoid both",
#                      x = .$d[[1]]$x,    
#                      y = .$d[[1]]$y - .$d[[7]]$y),
#           data.frame(Region = "STM",
#                      Scenario = "Real",
#                      x = .$d[[10]]$x,    
#                      y = .$d[[10]]$y - .$d[[15]]$y),
#           data.frame(Region = "STM",
#                      Scenario = "Avoid degradation",
#                      x = .$d[[10]]$x,    
#                      y = .$d[[10]]$y - .$d[[14]]$y),
#           data.frame(Region = "STM",
#                      Scenario = "Avoid deforestation",
#                      x = .$d[[10]]$x,    
#                      y = .$d[[10]]$y - .$d[[13]]$y),
#           data.frame(Region = "STM",
#                      Scenario = "Avoid deforestation primary forest",
#                      x = .$d[[10]]$x,    
#                      y = .$d[[10]]$y - .$d[[12]]$y),
#           data.frame(Region = "STM",
#                      Scenario = "Avoid both",
#                      x = .$d[[10]]$x,    
#                      y = .$d[[10]]$y - .$d[[11]]$y),
#           data.frame(Region = "STM",
#                      Scenario = "Restoration without avoid",
#                      x = .$d[[10]]$x,    
#                      y = .$d[[10]]$y - .$d[[18]]$y),
#           data.frame(Region = "STM",
#                      Scenario = "Restoration and avoid deforestation",
#                      x = .$d[[10]]$x,    
#                      y = .$d[[10]]$y - .$d[[17]]$y),
#           data.frame(Region = "STM",
#                      Scenario = "Restoration and avoid both",
#                      x = .$d[[10]]$x,    
#                      y = .$d[[10]]$y - .$d[[16]]$y))) %>%
#  mutate(Scenario = factor(Scenario, 
#                           levels = c("Real",
#                                      "Avoid degradation",
#                                      "Avoid deforestation",
#                                      "Avoid deforestation primary forest",
#                                      "Avoid both",
#                                      "Restoration without avoid",
#                                      "Restoration and avoid deforestation",
#                                      "Restoration and avoid both"))) %>% 
#  ggplot(aes(x=x)) +    # now plot
#  geom_area(aes(y=y), fill="#73A9AD") + 
#  geom_ribbon(aes(ymin = 0, ymax = ifelse(y >= 0, 0, y)), fill = "#F5F0BB") +
#  labs(title = '', y="Density", x = "Carbon benefit values (Mg C / ha)") +
#  theme_minimal()+
#  facet_wrap(Region~Scenario, nrow = 2) #+
#  #theme(text = element_text(size = 44, family = "sans"),
#  #      axis.title = element_text(face="bold"),
#  #      axis.text.x=element_text(size = 32))





real.tp <- data.frame(
  Period1 = pgm.carbon.benefit[["PGM_2010_real_carbon_benefit"]][],
  Period2 = pgm.carbon.benefit[["PGM_2020_real_carbon_benefit"]][]
)

real.tp <- real.tp %>% drop_na() %>% 
  mutate(Carbon_cat1 = factor(case_when(Period1 > 150 ~ ">150",
                                        Period1 > 126 & Period1 <= 150 ~ "126-150",
                                        Period1 > 101 & Period1 <= 125 ~ "101-125",
                                        Period1 > 76 & Period1 <= 100 ~ "76-100",
                                        Period1 > 51 & Period1 <= 75~ "51-75",
                                        Period1 > 25 & Period1 <= 50 ~ "26-50",
                                        Period1 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")),
         Carbon_cat2 = factor(case_when(Period2 > 150 ~ ">150",
                                        Period2 > 126 & Period2 <= 150 ~ "126-150",
                                        Period2 > 101 & Period2 <= 125 ~ "101-125",
                                        Period2 > 76 & Period2 <= 100 ~ "76-100",
                                        Period2 > 51 & Period2 <= 75~ "51-75",
                                        Period2 > 25 & Period2 <= 50 ~ "26-50",
                                        Period2 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")))


avoiddegrad.tp <- data.frame(
  Period1 = pgm.carbon.benefit[["PGM_2010_real_carbon_benefit"]][],
  Period2 = pgm.carbon.benefit[["PGM_2020_avoiddegrad_carbon_benefit"]][]
)

avoiddegrad.tp <- avoiddegrad.tp %>% drop_na() %>% 
  mutate(Carbon_cat1 = factor(case_when(Period1 > 150 ~ ">150",
                                        Period1 > 126 & Period1 <= 150 ~ "126-150",
                                        Period1 > 101 & Period1 <= 125 ~ "101-125",
                                        Period1 > 76 & Period1 <= 100 ~ "76-100",
                                        Period1 > 51 & Period1 <= 75~ "51-75",
                                        Period1 > 25 & Period1 <= 50 ~ "26-50",
                                        Period1 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")),
         Carbon_cat2 = factor(case_when(Period2 > 150 ~ ">150",
                                        Period2 > 126 & Period2 <= 150 ~ "126-150",
                                        Period2 > 101 & Period2 <= 125 ~ "101-125",
                                        Period2 > 76 & Period2 <= 100 ~ "76-100",
                                        Period2 > 51 & Period2 <= 75~ "51-75",
                                        Period2 > 25 & Period2 <= 50 ~ "26-50",
                                        Period2 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")))


avoiddeforest.tp <- data.frame(
  Period1 = pgm.carbon.benefit[["PGM_2010_real_carbon_benefit"]][],
  Period2 = pgm.carbon.benefit[["PGM_2020_avoiddeforest_carbon_benefit"]][]
)

avoiddeforest.tp <- avoiddeforest.tp %>% drop_na() %>% 
  mutate(Carbon_cat1 = factor(case_when(Period1 > 150 ~ ">150",
                                        Period1 > 126 & Period1 <= 150 ~ "126-150",
                                        Period1 > 101 & Period1 <= 125 ~ "101-125",
                                        Period1 > 76 & Period1 <= 100 ~ "76-100",
                                        Period1 > 51 & Period1 <= 75~ "51-75",
                                        Period1 > 25 & Period1 <= 50 ~ "26-50",
                                        Period1 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")),
         Carbon_cat2 = factor(case_when(Period2 > 150 ~ ">150",
                                        Period2 > 126 & Period2 <= 150 ~ "126-150",
                                        Period2 > 101 & Period2 <= 125 ~ "101-125",
                                        Period2 > 76 & Period2 <= 100 ~ "76-100",
                                        Period2 > 51 & Period2 <= 75~ "51-75",
                                        Period2 > 25 & Period2 <= 50 ~ "26-50",
                                        Period2 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")))


avoiddeforest2.tp <- data.frame(
  Period1 = pgm.carbon.benefit[["PGM_2010_real_carbon_benefit"]][],
  Period2 = pgm.carbon.benefit[["PGM_2020_avoiddeforest2_carbon_benefit"]][]
)

avoiddeforest2.tp <- avoiddeforest2.tp %>% drop_na() %>% 
  mutate(Carbon_cat1 = factor(case_when(Period1 > 150 ~ ">150",
                                        Period1 > 126 & Period1 <= 150 ~ "126-150",
                                        Period1 > 101 & Period1 <= 125 ~ "101-125",
                                        Period1 > 76 & Period1 <= 100 ~ "76-100",
                                        Period1 > 51 & Period1 <= 75~ "51-75",
                                        Period1 > 25 & Period1 <= 50 ~ "26-50",
                                        Period1 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")),
         Carbon_cat2 = factor(case_when(Period2 > 150 ~ ">150",
                                        Period2 > 126 & Period2 <= 150 ~ "126-150",
                                        Period2 > 101 & Period2 <= 125 ~ "101-125",
                                        Period2 > 76 & Period2 <= 100 ~ "76-100",
                                        Period2 > 51 & Period2 <= 75~ "51-75",
                                        Period2 > 25 & Period2 <= 50 ~ "26-50",
                                        Period2 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")))


avoidboth.tp <- data.frame(
  Period1 = pgm.carbon.benefit[["PGM_2010_real_carbon_benefit"]][],
  Period2 = pgm.carbon.benefit[["PGM_2020_avoidboth_carbon_benefit"]][]
)

avoidboth.tp <- avoidboth.tp %>% drop_na() %>% 
  mutate(Carbon_cat1 = factor(case_when(Period1 > 150 ~ ">150",
                                        Period1 > 126 & Period1 <= 150 ~ "126-150",
                                        Period1 > 101 & Period1 <= 125 ~ "101-125",
                                        Period1 > 76 & Period1 <= 100 ~ "76-100",
                                        Period1 > 51 & Period1 <= 75~ "51-75",
                                        Period1 > 25 & Period1 <= 50 ~ "26-50",
                                        Period1 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")),
         Carbon_cat2 = factor(case_when(Period2 > 150 ~ ">150",
                                        Period2 > 126 & Period2 <= 150 ~ "126-150",
                                        Period2 > 101 & Period2 <= 125 ~ "101-125",
                                        Period2 > 76 & Period2 <= 100 ~ "76-100",
                                        Period2 > 51 & Period2 <= 75~ "51-75",
                                        Period2 > 25 & Period2 <= 50 ~ "26-50",
                                        Period2 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")))


restor_wo_avoid.tp <- data.frame(
  Period1 = pgm.carbon.benefit[["PGM_2010_real_carbon_benefit"]][],
  Period2 = pgm.carbon.benefit[["PGM_2020_restor_wo_avoid_carbon_benefit"]][]
)

restor_wo_avoid.tp <- restor_wo_avoid.tp %>% drop_na() %>% 
  mutate(Carbon_cat1 = factor(case_when(Period1 > 150 ~ ">150",
                                        Period1 > 126 & Period1 <= 150 ~ "126-150",
                                        Period1 > 101 & Period1 <= 125 ~ "101-125",
                                        Period1 > 76 & Period1 <= 100 ~ "76-100",
                                        Period1 > 51 & Period1 <= 75~ "51-75",
                                        Period1 > 25 & Period1 <= 50 ~ "26-50",
                                        Period1 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")),
         Carbon_cat2 = factor(case_when(Period2 > 150 ~ ">150",
                                        Period2 > 126 & Period2 <= 150 ~ "126-150",
                                        Period2 > 101 & Period2 <= 125 ~ "101-125",
                                        Period2 > 76 & Period2 <= 100 ~ "76-100",
                                        Period2 > 51 & Period2 <= 75~ "51-75",
                                        Period2 > 25 & Period2 <= 50 ~ "26-50",
                                        Period2 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")))


restor_n_avoid_deforest.tp <- data.frame(
  Period1 = pgm.carbon.benefit[["PGM_2010_real_carbon_benefit"]][],
  Period2 = pgm.carbon.benefit[["PGM_2020_restor_n_avoid_deforest_carbon_benefit"]][]
)

restor_n_avoid_deforest.tp <- restor_n_avoid_deforest.tp %>% drop_na() %>% 
  mutate(Carbon_cat1 = factor(case_when(Period1 > 150 ~ ">150",
                                        Period1 > 126 & Period1 <= 150 ~ "126-150",
                                        Period1 > 101 & Period1 <= 125 ~ "101-125",
                                        Period1 > 76 & Period1 <= 100 ~ "76-100",
                                        Period1 > 51 & Period1 <= 75~ "51-75",
                                        Period1 > 25 & Period1 <= 50 ~ "26-50",
                                        Period1 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")),
         Carbon_cat2 = factor(case_when(Period2 > 150 ~ ">150",
                                        Period2 > 126 & Period2 <= 150 ~ "126-150",
                                        Period2 > 101 & Period2 <= 125 ~ "101-125",
                                        Period2 > 76 & Period2 <= 100 ~ "76-100",
                                        Period2 > 51 & Period2 <= 75~ "51-75",
                                        Period2 > 25 & Period2 <= 50 ~ "26-50",
                                        Period2 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")))


restor_n_avoid_both.tp <- data.frame(
  Period1 = pgm.carbon.benefit[["PGM_2010_real_carbon_benefit"]][],
  Period2 = pgm.carbon.benefit[["PGM_2020_restor_n_avoid_both_carbon_benefit"]][]
)

restor_n_avoid_both.tp <- restor_n_avoid_both.tp %>% drop_na() %>% 
  mutate(Carbon_cat1 = factor(case_when(Period1 > 150 ~ ">150",
                                        Period1 > 126 & Period1 <= 150 ~ "126-150",
                                        Period1 > 101 & Period1 <= 125 ~ "101-125",
                                        Period1 > 76 & Period1 <= 100 ~ "76-100",
                                        Period1 > 51 & Period1 <= 75~ "51-75",
                                        Period1 > 25 & Period1 <= 50 ~ "26-50",
                                        Period1 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")),
         Carbon_cat2 = factor(case_when(Period2 > 150 ~ ">150",
                                        Period2 > 126 & Period2 <= 150 ~ "126-150",
                                        Period2 > 101 & Period2 <= 125 ~ "101-125",
                                        Period2 > 76 & Period2 <= 100 ~ "76-100",
                                        Period2 > 51 & Period2 <= 75~ "51-75",
                                        Period2 > 25 & Period2 <= 50 ~ "26-50",
                                        Period2 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")))







addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

cowplot::plot_grid(
  
  real.tp %>% sample_n(size = 100000, replace = T) %>%
    ggplot(aes(axis1 = Carbon_cat1, axis2 = Carbon_cat2)) +
    geom_flow(aes(fill = Carbon_cat1), width = .15, curve_type = "quintic") +
    geom_stratum(width = .15) +
    scale_x_discrete(limits = c("Carbon_cat1", "Carbon_cat2"), 
                     breaks=c("Carbon_cat1", "Carbon_cat2"), 
                     labels=addline_format(c("2010 Real", "2020 Real")),
                     expand = c(.05, .05)) +
    scale_fill_manual(values = c("#65451F", "#83764F", "#263A29", "#7C9070", "#A4D0A4", "#FEE8B0", "#F97B22")) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal()+
    theme(axis.text.y= element_blank(), legend.position = "none"),
  
  
  avoiddegrad.tp %>% sample_n(size = 100000, replace = T) %>%
    ggplot(aes(axis1 = Carbon_cat1, axis2 = Carbon_cat2)) +
    geom_flow(aes(fill = Carbon_cat1), width = .15, curve_type = "quintic") +
    geom_stratum(width = .15) +
    scale_x_discrete(limits = c("Carbon_cat1", "Carbon_cat2"), 
                     breaks=c("Carbon_cat1", "Carbon_cat2"), 
                     labels=addline_format(c("2010 Real", "Avoid degradation")),
                     expand = c(.05, .05)) +
    scale_fill_manual(values = c("#65451F", "#83764F", "#263A29", "#7C9070", "#A4D0A4", "#FEE8B0", "#F97B22")) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal()+
    theme(axis.text.y= element_blank(), legend.position = "none"),
  
  
  avoiddeforest.tp %>% sample_n(size = 100000, replace = T) %>%
    ggplot(aes(axis1 = Carbon_cat1, axis2 = Carbon_cat2)) +
    geom_flow(aes(fill = Carbon_cat1), width = .15, curve_type = "quintic") +
    geom_stratum(width = .15) +
    scale_x_discrete(limits = c("Carbon_cat1", "Carbon_cat2"), 
                     breaks=c("Carbon_cat1", "Carbon_cat2"), 
                     labels=addline_format(c("2010 Real", "Avoid deforestation")),
                     expand = c(.05, .05)) +
    scale_fill_manual(values = c("#65451F", "#83764F", "#263A29", "#7C9070", "#A4D0A4", "#FEE8B0", "#F97B22")) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal()+
    theme(axis.text.y= element_blank(), legend.position = "none"),
  
  
  avoiddeforest2.tp %>% sample_n(size = 100000, replace = T) %>%
    ggplot(aes(axis1 = Carbon_cat1, axis2 = Carbon_cat2)) +
    geom_flow(aes(fill = Carbon_cat1), width = .15, curve_type = "quintic") +
    geom_stratum(width = .15) +
    scale_x_discrete(limits = c("Carbon_cat1", "Carbon_cat2"), 
                     breaks=c("Carbon_cat1", "Carbon_cat2"), 
                     labels=addline_format(c("2010 Real", "Avoid deforestation primary forest")),
                     expand = c(.05, .05)) +
    scale_fill_manual(values = c("#65451F", "#83764F", "#263A29", "#7C9070", "#A4D0A4", "#FEE8B0", "#F97B22")) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal()+
    theme(axis.text.y= element_blank(), legend.position = "none"),
  
  
  avoidboth.tp %>% sample_n(size = 100000, replace = T) %>%
    ggplot(aes(axis1 = Carbon_cat1, axis2 = Carbon_cat2)) +
    geom_flow(aes(fill = Carbon_cat1), width = .15, curve_type = "quintic") +
    geom_stratum(width = .15) +
    scale_x_discrete(limits = c("Carbon_cat1", "Carbon_cat2"), 
                     breaks=c("Carbon_cat1", "Carbon_cat2"), 
                     labels=addline_format(c("2010 Real", "Avoid both")),
                     expand = c(.05, .05)) +
    scale_fill_manual(values = c("#65451F", "#83764F", "#263A29", "#7C9070", "#A4D0A4", "#FEE8B0", "#F97B22")) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal()+
    theme(axis.text.y= element_blank(), legend.position = "none"),
  
  
  restor_wo_avoid.tp %>% sample_n(size = 100000, replace = T) %>%
    ggplot(aes(axis1 = Carbon_cat1, axis2 = Carbon_cat2)) +
    geom_flow(aes(fill = Carbon_cat1), width = .15, curve_type = "quintic") +
    geom_stratum(width = .15) +
    scale_x_discrete(limits = c("Carbon_cat1", "Carbon_cat2"), 
                     breaks=c("Carbon_cat1", "Carbon_cat2"), 
                     labels=addline_format(c("2010 Real", "Restoration without avoid")),
                     expand = c(.05, .05)) +
    scale_fill_manual(values = c("#65451F", "#83764F", "#263A29", "#7C9070", "#A4D0A4", "#FEE8B0", "#F97B22")) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal()+
    theme(axis.text.y= element_blank(), legend.position = "none"),
  
  
  restor_n_avoid_deforest.tp %>% sample_n(size = 100000, replace = T) %>%
    ggplot(aes(axis1 = Carbon_cat1, axis2 = Carbon_cat2)) +
    geom_flow(aes(fill = Carbon_cat1), width = .15, curve_type = "quintic") +
    geom_stratum(width = .15) +
    scale_x_discrete(limits = c("Carbon_cat1", "Carbon_cat2"), 
                     breaks=c("Carbon_cat1", "Carbon_cat2"), 
                     labels=addline_format(c("2010 Real", "Restoration and avoid deforestation")),
                     expand = c(.05, .05)) +
    scale_fill_manual(values = c("#65451F", "#83764F", "#263A29", "#7C9070", "#A4D0A4", "#FEE8B0", "#F97B22")) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal()+
    theme(axis.text.y= element_blank(), legend.position = "none"),
  
  
  restor_n_avoid_both.tp %>% sample_n(size = 100000, replace = T) %>%
    ggplot(aes(axis1 = Carbon_cat1, axis2 = Carbon_cat2)) +
    geom_flow(aes(fill = Carbon_cat1), width = .15, curve_type = "quintic") +
    geom_stratum(width = .15) +
    scale_x_discrete(limits = c("Carbon_cat1", "Carbon_cat2"), 
                     breaks=c("Carbon_cat1", "Carbon_cat2"), 
                     labels=addline_format(c("2010 Real", "Restoration and avoid both")),
                     expand = c(.05, .05)) +
    scale_fill_manual(values = c("#65451F", "#83764F", "#263A29", "#7C9070", "#A4D0A4", "#FEE8B0", "#F97B22")) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal()+
    theme(axis.text.y= element_blank(), legend.position = "none"),

ncol = 4, align = "hv")





real.tp <- data.frame(
  Period1 = stm.carbon.benefit[["STM_2010_real_carbon_benefit"]][],
  Period2 = stm.carbon.benefit[["STM_2020_real_carbon_benefit"]][]
)

real.tp <- real.tp %>% drop_na() %>% 
  mutate(Carbon_cat1 = factor(case_when(Period1 > 150 ~ ">150",
                                        Period1 > 126 & Period1 <= 150 ~ "126-150",
                                        Period1 > 101 & Period1 <= 125 ~ "101-125",
                                        Period1 > 76 & Period1 <= 100 ~ "76-100",
                                        Period1 > 51 & Period1 <= 75~ "51-75",
                                        Period1 > 25 & Period1 <= 50 ~ "26-50",
                                        Period1 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")),
         Carbon_cat2 = factor(case_when(Period2 > 150 ~ ">150",
                                        Period2 > 126 & Period2 <= 150 ~ "126-150",
                                        Period2 > 101 & Period2 <= 125 ~ "101-125",
                                        Period2 > 76 & Period2 <= 100 ~ "76-100",
                                        Period2 > 51 & Period2 <= 75~ "51-75",
                                        Period2 > 25 & Period2 <= 50 ~ "26-50",
                                        Period2 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")))


avoiddegrad.tp <- data.frame(
  Period1 = stm.carbon.benefit[["STM_2010_real_carbon_benefit"]][],
  Period2 = stm.carbon.benefit[["STM_2020_avoiddegrad_carbon_benefit"]][]
)

avoiddegrad.tp <- avoiddegrad.tp %>% drop_na() %>% 
  mutate(Carbon_cat1 = factor(case_when(Period1 > 150 ~ ">150",
                                        Period1 > 126 & Period1 <= 150 ~ "126-150",
                                        Period1 > 101 & Period1 <= 125 ~ "101-125",
                                        Period1 > 76 & Period1 <= 100 ~ "76-100",
                                        Period1 > 51 & Period1 <= 75~ "51-75",
                                        Period1 > 25 & Period1 <= 50 ~ "26-50",
                                        Period1 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")),
         Carbon_cat2 = factor(case_when(Period2 > 150 ~ ">150",
                                        Period2 > 126 & Period2 <= 150 ~ "126-150",
                                        Period2 > 101 & Period2 <= 125 ~ "101-125",
                                        Period2 > 76 & Period2 <= 100 ~ "76-100",
                                        Period2 > 51 & Period2 <= 75~ "51-75",
                                        Period2 > 25 & Period2 <= 50 ~ "26-50",
                                        Period2 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")))


avoiddeforest.tp <- data.frame(
  Period1 = stm.carbon.benefit[["STM_2010_real_carbon_benefit"]][],
  Period2 = stm.carbon.benefit[["STM_2020_avoiddeforest_carbon_benefit"]][]
)

avoiddeforest.tp <- avoiddeforest.tp %>% drop_na() %>% 
  mutate(Carbon_cat1 = factor(case_when(Period1 > 150 ~ ">150",
                                        Period1 > 126 & Period1 <= 150 ~ "126-150",
                                        Period1 > 101 & Period1 <= 125 ~ "101-125",
                                        Period1 > 76 & Period1 <= 100 ~ "76-100",
                                        Period1 > 51 & Period1 <= 75~ "51-75",
                                        Period1 > 25 & Period1 <= 50 ~ "26-50",
                                        Period1 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")),
         Carbon_cat2 = factor(case_when(Period2 > 150 ~ ">150",
                                        Period2 > 126 & Period2 <= 150 ~ "126-150",
                                        Period2 > 101 & Period2 <= 125 ~ "101-125",
                                        Period2 > 76 & Period2 <= 100 ~ "76-100",
                                        Period2 > 51 & Period2 <= 75~ "51-75",
                                        Period2 > 25 & Period2 <= 50 ~ "26-50",
                                        Period2 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")))


avoiddeforest2.tp <- data.frame(
  Period1 = stm.carbon.benefit[["STM_2010_real_carbon_benefit"]][],
  Period2 = stm.carbon.benefit[["STM_2020_avoiddeforest2_carbon_benefit"]][]
)

avoiddeforest2.tp <- avoiddeforest2.tp %>% drop_na() %>% 
  mutate(Carbon_cat1 = factor(case_when(Period1 > 150 ~ ">150",
                                        Period1 > 126 & Period1 <= 150 ~ "126-150",
                                        Period1 > 101 & Period1 <= 125 ~ "101-125",
                                        Period1 > 76 & Period1 <= 100 ~ "76-100",
                                        Period1 > 51 & Period1 <= 75~ "51-75",
                                        Period1 > 25 & Period1 <= 50 ~ "26-50",
                                        Period1 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")),
         Carbon_cat2 = factor(case_when(Period2 > 150 ~ ">150",
                                        Period2 > 126 & Period2 <= 150 ~ "126-150",
                                        Period2 > 101 & Period2 <= 125 ~ "101-125",
                                        Period2 > 76 & Period2 <= 100 ~ "76-100",
                                        Period2 > 51 & Period2 <= 75~ "51-75",
                                        Period2 > 25 & Period2 <= 50 ~ "26-50",
                                        Period2 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")))


avoidboth.tp <- data.frame(
  Period1 = stm.carbon.benefit[["STM_2010_real_carbon_benefit"]][],
  Period2 = stm.carbon.benefit[["STM_2020_avoidboth_carbon_benefit"]][]
)

avoidboth.tp <- avoidboth.tp %>% drop_na() %>% 
  mutate(Carbon_cat1 = factor(case_when(Period1 > 150 ~ ">150",
                                        Period1 > 126 & Period1 <= 150 ~ "126-150",
                                        Period1 > 101 & Period1 <= 125 ~ "101-125",
                                        Period1 > 76 & Period1 <= 100 ~ "76-100",
                                        Period1 > 51 & Period1 <= 75~ "51-75",
                                        Period1 > 25 & Period1 <= 50 ~ "26-50",
                                        Period1 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")),
         Carbon_cat2 = factor(case_when(Period2 > 150 ~ ">150",
                                        Period2 > 126 & Period2 <= 150 ~ "126-150",
                                        Period2 > 101 & Period2 <= 125 ~ "101-125",
                                        Period2 > 76 & Period2 <= 100 ~ "76-100",
                                        Period2 > 51 & Period2 <= 75~ "51-75",
                                        Period2 > 25 & Period2 <= 50 ~ "26-50",
                                        Period2 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")))


restor_wo_avoid.tp <- data.frame(
  Period1 = stm.carbon.benefit[["STM_2010_real_carbon_benefit"]][],
  Period2 = stm.carbon.benefit[["STM_2020_restor_wo_avoid_carbon_benefit"]][]
)

restor_wo_avoid.tp <- restor_wo_avoid.tp %>% drop_na() %>% 
  mutate(Carbon_cat1 = factor(case_when(Period1 > 150 ~ ">150",
                                        Period1 > 126 & Period1 <= 150 ~ "126-150",
                                        Period1 > 101 & Period1 <= 125 ~ "101-125",
                                        Period1 > 76 & Period1 <= 100 ~ "76-100",
                                        Period1 > 51 & Period1 <= 75~ "51-75",
                                        Period1 > 25 & Period1 <= 50 ~ "26-50",
                                        Period1 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")),
         Carbon_cat2 = factor(case_when(Period2 > 150 ~ ">150",
                                        Period2 > 126 & Period2 <= 150 ~ "126-150",
                                        Period2 > 101 & Period2 <= 125 ~ "101-125",
                                        Period2 > 76 & Period2 <= 100 ~ "76-100",
                                        Period2 > 51 & Period2 <= 75~ "51-75",
                                        Period2 > 25 & Period2 <= 50 ~ "26-50",
                                        Period2 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")))


restor_n_avoid_deforest.tp <- data.frame(
  Period1 = stm.carbon.benefit[["STM_2010_real_carbon_benefit"]][],
  Period2 = stm.carbon.benefit[["STM_2020_restor_n_avoid_deforest_carbon_benefit"]][]
)

restor_n_avoid_deforest.tp <- restor_n_avoid_deforest.tp %>% drop_na() %>% 
  mutate(Carbon_cat1 = factor(case_when(Period1 > 150 ~ ">150",
                                        Period1 > 126 & Period1 <= 150 ~ "126-150",
                                        Period1 > 101 & Period1 <= 125 ~ "101-125",
                                        Period1 > 76 & Period1 <= 100 ~ "76-100",
                                        Period1 > 51 & Period1 <= 75~ "51-75",
                                        Period1 > 25 & Period1 <= 50 ~ "26-50",
                                        Period1 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")),
         Carbon_cat2 = factor(case_when(Period2 > 150 ~ ">150",
                                        Period2 > 126 & Period2 <= 150 ~ "126-150",
                                        Period2 > 101 & Period2 <= 125 ~ "101-125",
                                        Period2 > 76 & Period2 <= 100 ~ "76-100",
                                        Period2 > 51 & Period2 <= 75~ "51-75",
                                        Period2 > 25 & Period2 <= 50 ~ "26-50",
                                        Period2 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")))


restor_n_avoid_both.tp <- data.frame(
  Period1 = stm.carbon.benefit[["STM_2010_real_carbon_benefit"]][],
  Period2 = stm.carbon.benefit[["STM_2020_restor_n_avoid_both_carbon_benefit"]][]
)

restor_n_avoid_both.tp <- restor_n_avoid_both.tp %>% drop_na() %>% 
  mutate(Carbon_cat1 = factor(case_when(Period1 > 150 ~ ">150",
                                        Period1 > 126 & Period1 <= 150 ~ "126-150",
                                        Period1 > 101 & Period1 <= 125 ~ "101-125",
                                        Period1 > 76 & Period1 <= 100 ~ "76-100",
                                        Period1 > 51 & Period1 <= 75~ "51-75",
                                        Period1 > 25 & Period1 <= 50 ~ "26-50",
                                        Period1 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")),
         Carbon_cat2 = factor(case_when(Period2 > 150 ~ ">150",
                                        Period2 > 126 & Period2 <= 150 ~ "126-150",
                                        Period2 > 101 & Period2 <= 125 ~ "101-125",
                                        Period2 > 76 & Period2 <= 100 ~ "76-100",
                                        Period2 > 51 & Period2 <= 75~ "51-75",
                                        Period2 > 25 & Period2 <= 50 ~ "26-50",
                                        Period2 <= 25 ~ "0-25"), 
                              levels = c(">150",
                                         "126-150",
                                         "101-125",
                                         "76-100",
                                         "51-75",
                                         "26-50",
                                         "0-25")))





cowplot::plot_grid(
  
  real.tp %>% sample_n(size = 100000, replace = T) %>%
    ggplot(aes(axis1 = Carbon_cat1, axis2 = Carbon_cat2)) +
    geom_flow(aes(fill = Carbon_cat1), width = .15, curve_type = "quintic") +
    geom_stratum(width = .15) +
    scale_x_discrete(limits = c("Carbon_cat1", "Carbon_cat2"), 
                     breaks=c("Carbon_cat1", "Carbon_cat2"), 
                     labels=addline_format(c("2010 Real", "2020 Real")),
                     expand = c(.05, .05)) +
    scale_fill_manual(values = c("#65451F", "#83764F", "#263A29", "#7C9070", "#A4D0A4", "#FEE8B0", "#F97B22")) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal()+
    theme(axis.text.y= element_blank(), legend.position = "none"),
  
  
  avoiddegrad.tp %>% sample_n(size = 100000, replace = T) %>%
    ggplot(aes(axis1 = Carbon_cat1, axis2 = Carbon_cat2)) +
    geom_flow(aes(fill = Carbon_cat1), width = .15, curve_type = "quintic") +
    geom_stratum(width = .15) +
    scale_x_discrete(limits = c("Carbon_cat1", "Carbon_cat2"), 
                     breaks=c("Carbon_cat1", "Carbon_cat2"), 
                     labels=addline_format(c("2010 Real", "Avoid degradation")),
                     expand = c(.05, .05)) +
    scale_fill_manual(values = c("#65451F", "#83764F", "#263A29", "#7C9070", "#A4D0A4", "#FEE8B0", "#F97B22")) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal()+
    theme(axis.text.y= element_blank(), legend.position = "none"),
  
  
  avoiddeforest.tp %>% sample_n(size = 100000, replace = T) %>%
    ggplot(aes(axis1 = Carbon_cat1, axis2 = Carbon_cat2)) +
    geom_flow(aes(fill = Carbon_cat1), width = .15, curve_type = "quintic") +
    geom_stratum(width = .15) +
    scale_x_discrete(limits = c("Carbon_cat1", "Carbon_cat2"), 
                     breaks=c("Carbon_cat1", "Carbon_cat2"), 
                     labels=addline_format(c("2010 Real", "Avoid deforestation")),
                     expand = c(.05, .05)) +
    scale_fill_manual(values = c("#65451F", "#83764F", "#263A29", "#7C9070", "#A4D0A4", "#FEE8B0", "#F97B22")) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal()+
    theme(axis.text.y= element_blank(), legend.position = "none"),
  
  
  avoiddeforest2.tp %>% sample_n(size = 100000, replace = T) %>%
    ggplot(aes(axis1 = Carbon_cat1, axis2 = Carbon_cat2)) +
    geom_flow(aes(fill = Carbon_cat1), width = .15, curve_type = "quintic") +
    geom_stratum(width = .15) +
    scale_x_discrete(limits = c("Carbon_cat1", "Carbon_cat2"), 
                     breaks=c("Carbon_cat1", "Carbon_cat2"), 
                     labels=addline_format(c("2010 Real", "Avoid deforestation primary forest")),
                     expand = c(.05, .05)) +
    scale_fill_manual(values = c("#65451F", "#83764F", "#263A29", "#7C9070", "#A4D0A4", "#FEE8B0", "#F97B22")) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal()+
    theme(axis.text.y= element_blank(), legend.position = "none"),
  
  
  avoidboth.tp %>% sample_n(size = 100000, replace = T) %>%
    ggplot(aes(axis1 = Carbon_cat1, axis2 = Carbon_cat2)) +
    geom_flow(aes(fill = Carbon_cat1), width = .15, curve_type = "quintic") +
    geom_stratum(width = .15) +
    scale_x_discrete(limits = c("Carbon_cat1", "Carbon_cat2"), 
                     breaks=c("Carbon_cat1", "Carbon_cat2"), 
                     labels=addline_format(c("2010 Real", "Avoid both")),
                     expand = c(.05, .05)) +
    scale_fill_manual(values = c("#65451F", "#83764F", "#263A29", "#7C9070", "#A4D0A4", "#FEE8B0", "#F97B22")) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal()+
    theme(axis.text.y= element_blank(), legend.position = "none"),
  
  
  restor_wo_avoid.tp %>% sample_n(size = 100000, replace = T) %>%
    ggplot(aes(axis1 = Carbon_cat1, axis2 = Carbon_cat2)) +
    geom_flow(aes(fill = Carbon_cat1), width = .15, curve_type = "quintic") +
    geom_stratum(width = .15) +
    scale_x_discrete(limits = c("Carbon_cat1", "Carbon_cat2"), 
                     breaks=c("Carbon_cat1", "Carbon_cat2"), 
                     labels=addline_format(c("2010 Real", "Restoration without avoid")),
                     expand = c(.05, .05)) +
    scale_fill_manual(values = c("#65451F", "#83764F", "#263A29", "#7C9070", "#A4D0A4", "#FEE8B0", "#F97B22")) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal()+
    theme(axis.text.y= element_blank(), legend.position = "none"),
  
  
  restor_n_avoid_deforest.tp %>% sample_n(size = 100000, replace = T) %>%
    ggplot(aes(axis1 = Carbon_cat1, axis2 = Carbon_cat2)) +
    geom_flow(aes(fill = Carbon_cat1), width = .15, curve_type = "quintic") +
    geom_stratum(width = .15) +
    scale_x_discrete(limits = c("Carbon_cat1", "Carbon_cat2"), 
                     breaks=c("Carbon_cat1", "Carbon_cat2"), 
                     labels=addline_format(c("2010 Real", "Restoration and avoid deforestation")),
                     expand = c(.05, .05)) +
    scale_fill_manual(values = c("#65451F", "#83764F", "#263A29", "#7C9070", "#A4D0A4", "#FEE8B0", "#F97B22")) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal()+
    theme(axis.text.y= element_blank(), legend.position = "none"),
  
  
  restor_n_avoid_both.tp %>% sample_n(size = 100000, replace = T) %>%
    ggplot(aes(axis1 = Carbon_cat1, axis2 = Carbon_cat2)) +
    geom_flow(aes(fill = Carbon_cat1), width = .15, curve_type = "quintic") +
    geom_stratum(width = .15) +
    scale_x_discrete(limits = c("Carbon_cat1", "Carbon_cat2"), 
                     breaks=c("Carbon_cat1", "Carbon_cat2"), 
                     labels=addline_format(c("2010 Real", "Restoration and avoid both")),
                     expand = c(.05, .05)) +
    scale_fill_manual(values = c("#65451F", "#83764F", "#263A29", "#7C9070", "#A4D0A4", "#FEE8B0", "#F97B22")) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal()+
    theme(axis.text.y= element_blank(), legend.position = "none"),
  
  ncol = 4, align = "hv")





##checking areas that turned low carbon stock after the conservation strategy
#above150 <- stm.carbon.benefit[["STM_2010_real_carbon_benefit"]]
#above150[] <- ifelse(stm.carbon.benefit[["STM_2010_real_carbon_benefit"]][] > 150, 1, NA)
#teste <- mask(stm.carbon.benefit[["STM_2020_restor_n_avoid_both_carbon_benefit"]], above150)
#teste[] <- ifelse(teste[] <= 150, 1, NA)
#plot(teste, col="black")
#plot(stm.shp, add=T)
#length(which(!is.na(teste[])))
#
#lulc_who_the_decrasing <- mask(LULC2020_restor_n_avoid_both, teste)
#table(lulc_who_the_decrasing[])
#
#upf_who_the_decrasing <- mask(raster("rasters/STM/2020_restor_n_avoid_both/UPFls.tif"), teste)
#tsd_who_the_decrasing <- mask(raster("rasters/STM/2020_restor_n_avoid_both/TSDls.tif"), teste)
#ed_who_the_decrasing <- mask(raster("rasters/STM/2020_restor_n_avoid_both/edgedist.tif"), teste)
#
#par(mfrow=c(2,2))
#
#hist(upf_who_the_decrasing[])
#  
#hist(tsd_who_the_decrasing[])
#  
#hist(ed_who_the_decrasing[])
#
#par(mfrow=c(1,1))





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
pgm.conservact.avoiddeforest2.carbbenefitmask <- mask(pgm.carbon.benefit[["PGM_2020_avoiddeforest2_carbon_benefit"]], pgm.conservationaction.avoiddeforest2)
pgm.conservact.avoidboth.carbbenefitmask <- mask(pgm.carbon.benefit[["PGM_2020_avoidboth_carbon_benefit"]], pgm.conservationaction.avoidboth)
pgm.conservact.restor_wo_avoid.carbbenefitmask <- mask(pgm.carbon.benefit[["PGM_2020_restor_wo_avoid_carbon_benefit"]], pgm.conservationaction.restor_wo_avoid)
pgm.conservact.restor_n_avoid_deforest.carbbenefitmask <- mask(pgm.carbon.benefit[["PGM_2020_restor_n_avoid_deforest_carbon_benefit"]], pgm.conservationaction.restor_n_avoid_deforest)
pgm.conservact.restor_n_avoid_both.carbbenefitmask <- mask(pgm.carbon.benefit[["PGM_2020_restor_n_avoid_both_carbon_benefit"]], pgm.conservationaction.restor_n_avoid_both)

pgm.conservact.carbbenefitmask <- stack(pgm.conservact.avoiddegrad.carbbenefitmask,
                                        pgm.conservact.avoiddeforest.carbbenefitmask,
                                        pgm.conservact.avoiddeforest2.carbbenefitmask,
                                        pgm.conservact.avoidboth.carbbenefitmask,
                                        pgm.conservact.restor_wo_avoid.carbbenefitmask,
                                        pgm.conservact.restor_n_avoid_deforest.carbbenefitmask,
                                        pgm.conservact.restor_n_avoid_both.carbbenefitmask)

names(pgm.conservact.carbbenefitmask) <- c("pgm_avoiddegrad_carbbenefitmask",
                                           "pgm_avoiddeforest_carbbenefitmask",
                                           "pgm_avoiddeforest2_carbbenefitmask",
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
                                str_detect(ID, "avoiddeforest2")~ "Avoid deforestation primary forest",
                                str_detect(ID, "avoiddeforest")~ "Avoid deforestation",
                                str_detect(ID, "avoidboth")~ "Avoid both",
                                str_detect(ID, "restor_wo_avoid")~ "Restoration without avoid",
                                str_detect(ID, "restor_n_avoid_deforest")~ "Restoration and avoid deforestation",
                                str_detect(ID, "restor_n_avoid_both")~ "Restoration and avoid both",
                                str_detect(ID, "real")~ "Real"),
                      levels = c("Real",
                                 "Avoid degradation",
                                 "Avoid deforestation",
                                 "Avoid deforestation primary forest",
                                 "Avoid both",
                                 "Restoration without avoid",
                                 "Restoration and avoid deforestation",
                                 "Restoration and avoid both")),
    Carbon = na_if(Carbon, 0)
  )


stm.conservact.avoiddegrad.carbbenefitmask <- mask(stm.carbon.benefit[["STM_2020_avoiddegrad_carbon_benefit"]], stm.conservationaction.avoiddegrad)
stm.conservact.avoiddeforest.carbbenefitmask <- mask(stm.carbon.benefit[["STM_2020_avoiddeforest_carbon_benefit"]], stm.conservationaction.avoiddeforest)
stm.conservact.avoiddeforest2.carbbenefitmask <- mask(stm.carbon.benefit[["STM_2020_avoiddeforest2_carbon_benefit"]], stm.conservationaction.avoiddeforest2)
stm.conservact.avoidboth.carbbenefitmask <- mask(stm.carbon.benefit[["STM_2020_avoidboth_carbon_benefit"]], stm.conservationaction.avoidboth)
stm.conservact.restor_wo_avoid.carbbenefitmask <- mask(stm.carbon.benefit[["STM_2020_restor_wo_avoid_carbon_benefit"]], stm.conservationaction.restor_wo_avoid)
stm.conservact.restor_n_avoid_deforest.carbbenefitmask <- mask(stm.carbon.benefit[["STM_2020_restor_n_avoid_deforest_carbon_benefit"]], stm.conservationaction.restor_n_avoid_deforest)
stm.conservact.restor_n_avoid_both.carbbenefitmask <- mask(stm.carbon.benefit[["STM_2020_restor_n_avoid_both_carbon_benefit"]], stm.conservationaction.restor_n_avoid_both)

stm.conservact.carbbenefitmask <- stack(stm.conservact.avoiddegrad.carbbenefitmask,
                                        stm.conservact.avoiddeforest.carbbenefitmask,
                                        stm.conservact.avoiddeforest2.carbbenefitmask,
                                        stm.conservact.avoidboth.carbbenefitmask,
                                        stm.conservact.restor_wo_avoid.carbbenefitmask,
                                        stm.conservact.restor_n_avoid_deforest.carbbenefitmask,
                                        stm.conservact.restor_n_avoid_both.carbbenefitmask)

names(stm.conservact.carbbenefitmask) <- c("stm_avoiddegrad_carbbenefitmask",
                                           "stm_avoiddeforest_carbbenefitmask",
                                           "stm_avoiddeforest2_carbbenefitmask",
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
                                str_detect(ID, "avoiddeforest2")~ "Avoid deforestation primary forest",
                                str_detect(ID, "avoiddeforest")~ "Avoid deforestation",
                                str_detect(ID, "avoidboth")~ "Avoid both",
                                str_detect(ID, "restor_wo_avoid")~ "Restoration without avoid",
                                str_detect(ID, "restor_n_avoid_deforest")~ "Restoration and avoid deforestation",
                                str_detect(ID, "restor_n_avoid_both")~ "Restoration and avoid both",
                                str_detect(ID, "real")~ "Real"),
                      levels = c("Real",
                                 "Avoid degradation",
                                 "Avoid deforestation",
                                 "Avoid deforestation primary forest",
                                 "Avoid both",
                                 "Restoration without avoid",
                                 "Restoration and avoid deforestation",
                                 "Restoration and avoid both")),
    Carbon = na_if(Carbon, 0)
  )


conservact.carbbenefitmask.df <- rbind(pgm.conservact.carbbenefitmask.df, stm.conservact.carbbenefitmask.df)

conservact.carbbenefitmask.df %>% group_by(Region, Scenario) %>% drop_na() %>% 
  summarise(nvalues=n(), mean.benefit = mean(Carbon, na.rm=T))

conservact.carbbenefit.overview <- conservact.carbbenefitmask.df %>% group_by(Region, Scenario) %>% drop_na() %>%
  sample_n(size = 20000, replace = T) %>% 
  summarise(
    Carbon_mean = mean(Carbon, na.rm=T),
    Carbon_sd = sd(Carbon, na.rm=T)
  ) %>% 
  ungroup()






conservact.carbbenefit.overview %>%
  pivot_longer(Carbon_mean:Carbon_sd,
               names_to = c("set", ".value"),
               names_pattern = "(.+)_(.+)"
  ) %>% 
  ggplot(aes(x = Scenario, y = mean, fill = Region)) +
  geom_errorbar(aes(ymin = mean, ymax = mean + sd),
                width = 0.2, colour = "black",
                position = position_dodge(0.9)
  ) + 
  geom_bar(
    stat = "identity", position = position_dodge(width = .9), 
    color="black", alpha=.7#, show_guide=FALSE
  ) +
  scale_fill_manual(values = c("#440154", "#374f00")) +
  scale_x_discrete(breaks=unique(conservact.carbbenefit.overview$Scenario), 
                   labels=addline_format(c("Avoid degradation", 
                                           "Avoid deforestation",
                                           "Avoid deforestation primary forest",
                                           "Avoid both",
                                           "Restoration without avoid", 
                                           "Restoration and avoid deforestation", 
                                           "Restoration and avoid both"))) + 
  #geom_text(
  #  aes(y= mean + sd, x= Scenario, group = Region, 
  #      label = paste0("(", min, " - ", max, ")")), 
  #  position = position_dodge(width = .9), vjust = -1
  #) +
  labs(x="", y="Carbon (Mg C / ha)") +
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        text = element_text(size = 22, family = "sans"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 18))


#conservact.carbbenefitmask.df %>% 
#  ggplot(aes(x=Scenario, y=Carbon, fill = Scenario)) +
#  geom_violin(show.legend = F) + 
#  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange", width = 0.2, show.legend = F) +
#  scale_fill_manual(values = c("#41644A", "#41644A", "#41644A", 
#                               "#D3756B", "#D3756B", "#D3756B")) +
#  coord_flip() +
#  facet_wrap(~Region) +
#  labs(title = "", x = "", y = "Carbon benefit values") +
#  theme_minimal()+
#  theme(text = element_text(size = 16, family = "sans"),
#        axis.title = element_text(face="bold"),
#        axis.text.x=element_text(size = 14))


#########################################
####          comparing costs        ####
####      between scenarios and      ####
####         real projections        ####
####           (landscape)           ####
#########################################
#dir.create("models.output/biodiversity.cost.benefits", recursive = T)
#dir.create("models.output/carbon.cost.benefits", recursive = T)

##PGM
#import costs
pgm.opportunity.cost <- raster("models.output/opportunity.costs/PGM_2010_real_base_opportunity_cost.tif")
#pgm.opportunity.cost[is.na(pgm.opportunity.cost)]<-0
pgm.harvest.cost <- raster("models.output/opportunity.costs/PGM_2010_real_base_harvest_cost.tif")
pgm.stack.opportunity.cost <- stack(pgm.opportunity.cost,pgm.harvest.cost)
pgm.max.opportunity.cost <- max(pgm.stack.opportunity.cost, na.rm = T)
pgm.firecontrol.cost <- raster("models.output/opportunity.costs/PGM_2010_real_base_firecontrol.tif")
pgm.passiverestor.cost <- raster("models.output/opportunity.costs/PGM_2010_real_base_passiverestoration.tif")


#setting scenario costs
pgm.avoiddegrad.cost <- pgm.harvest.cost + pgm.firecontrol.cost
pgm.avoiddeforest.cost <- pgm.opportunity.cost + pgm.harvest.cost
pgm.avoidboth.cost <- pgm.max.opportunity.cost + pgm.firecontrol.cost
pgm.restor_wo_avoid.cost <- pgm.passiverestor.cost + pgm.opportunity.cost
pgm.restor_n_avoid_deforest.cost <- pgm.passiverestor.cost + pgm.opportunity.cost + pgm.harvest.cost
pgm.restor_n_avoid_both.cost <- pgm.passiverestor.cost + pgm.max.opportunity.cost + pgm.firecontrol.cost

#converting in dataframe
pgm.costs <- stack(c(pgm.avoiddegrad.cost, pgm.avoiddeforest.cost, pgm.avoidboth.cost, pgm.restor_wo_avoid.cost,
                     pgm.restor_n_avoid_deforest.cost, pgm.restor_n_avoid_both.cost))

names(pgm.costs) <- c("pgm_avoiddegrad_cost", "pgm_avoiddeforest_cost", "pgm_avoidboth_cost", "pgm_restor_wo_avoid_cost",
                    "pgm_restor_n_avoid_deforest_cost", "pgm_restor_n_avoid_both_cost")

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
stm.opportunity.cost <- mask(stm.opportunity.cost, stm.mask)
#stm.opportunity.cost[is.na(stm.opportunity.cost)]<-0
stm.harvest.cost <- raster("models.output/opportunity.costs/STM_2010_real_base_harvest_cost.tif")
stm.harvest.cost <- mask(stm.harvest.cost, stm.mask)
stm.stack.opportunity.cost <- stack(stm.opportunity.cost,stm.harvest.cost)
stm.max.opportunity.cost <- max(stm.stack.opportunity.cost, na.rm = T)
stm.firecontrol.cost <- raster("models.output/opportunity.costs/STM_2010_real_base_firecontrol.tif")
stm.firecontrol.cost <- mask(stm.firecontrol.cost, stm.mask)
stm.passiverestor.cost <- raster("models.output/opportunity.costs/STM_2010_real_base_passiverestoration.tif")
stm.passiverestor.cost <- mask(stm.passiverestor.cost, stm.mask)


#setting scenario costs
stm.avoiddegrad.cost <- stm.harvest.cost + stm.firecontrol.cost
stm.avoiddeforest.cost <- stm.opportunity.cost + stm.harvest.cost
stm.avoidboth.cost <- stm.max.opportunity.cost + stm.firecontrol.cost
stm.restor_wo_avoid.cost <- stm.passiverestor.cost + stm.opportunity.cost
stm.restor_n_avoid_deforest.cost <- stm.passiverestor.cost + stm.opportunity.cost + stm.harvest.cost
stm.restor_n_avoid_both.cost <- stm.passiverestor.cost + stm.max.opportunity.cost + stm.firecontrol.cost

#converting in dataframe
stm.costs <- stack(c(stm.avoiddegrad.cost, stm.avoiddeforest.cost, stm.avoidboth.cost, stm.restor_wo_avoid.cost,
                     stm.restor_n_avoid_deforest.cost, stm.restor_n_avoid_both.cost))

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



costs.overview <- costs %>% group_by(Region, Scenario) %>% 
  summarise(
    Costs_mean = mean(Costs, na.rm=T),
    Costs_sd = sd(Costs, na.rm=T),
    Costs_min = round(min(Costs, na.rm=T),1),
    Costs_max = round(max(Costs, na.rm=T),1)
  ) %>% 
  ungroup()




costs.overview %>%
  pivot_longer(Costs_mean:Costs_max,
               names_to = c("set", ".value"),
               names_pattern = "(.+)_(.+)"
  ) %>% 
  ggplot(aes(x = Scenario, y = mean, fill = Region)) +
  geom_errorbar(aes(ymin = mean, ymax = mean + sd),
                width = 0.2, colour = "black",
                position = position_dodge(0.9)
  ) + 
  geom_bar(
    stat = "identity", position = position_dodge(width = .9), 
    color="black", alpha=.7#, show_guide=FALSE
  ) +
  scale_fill_manual(values = c("#440154", "#374f00")) +
  scale_x_discrete(breaks=unique(costs.overview$Scenario), 
                   labels=addline_format(c("Avoid degradation", 
                                           "Avoid deforestation", 
                                           "Avoid both",
                                           "Restoration without avoid", 
                                           "Restoration and avoid deforestation", 
                                           "Restoration and avoid both"))) + 
  geom_text(
    aes(y= mean + sd, x= Scenario, group = Region, 
        label = paste0("(", min, " - ", max, ")")), 
    position = position_dodge(width = .9), vjust = -1
  ) +
  labs(x="", y="Costs (R$ / ha / year)") +
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        text = element_text(size = 22, family = "sans"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 18))


#par(mfrow=c(1,2)) ##res = 1326 x 876
#plot(pgm.opportunity.cost, col = terrain.colors(length(seq(0, 1500, by = 125)), rev = T), breaks= seq(0, 1500, by = 125))
#plot(stm.opportunity.cost, col = terrain.colors(length(seq(0, 1500, by = 125)), rev = T), breaks= seq(0, 1500, by = 125))

#plot(pgm.harvest.cost, col = terrain.colors(length(seq(0, 200, by = 20)), rev = T), breaks= seq(0, 200, by = 20))
#plot(stm.harvest.cost, col = terrain.colors(length(seq(0, 200, by = 20)), rev = T), breaks= seq(0, 200, by = 20))

#plot(pgm.firecontrol.cost, col = terrain.colors(length(seq(0, 600, by = 50)), rev = T), breaks= seq(0, 600, by = 50))
#plot(stm.firecontrol.cost, col = terrain.colors(length(seq(0, 600, by = 50)), rev = T), breaks= seq(0, 600, by = 50))

#plot(pgm.passiverestor.cost, col = terrain.colors(length(seq(0, 8000, by = 500)), rev = T), breaks= seq(0, 8000, by = 500))
#plot(stm.passiverestor.cost, col = terrain.colors(length(seq(0, 8000, by = 500)), rev = T), breaks= seq(0, 8000, by = 500))

##res = 1673 x 881
#plot(pgm.costs)
#plot(stm.costs)


#######################################################
##### analyzing the benefit-cost ratio (landscape) ####
#######################################################
#benefit.cost.ratio.df <- costs %>% dplyr::select(-ID) %>% 
#  left_join(
#    (biodiversity.benefit %>% dplyr::select(-ID) %>% filter(Scenario!="Real")),
#    by = c("x", "y", "Region", "Scenario")
#  )
#
#benefit.cost.ratio.df <- benefit.cost.ratio.df %>% 
#  left_join(
#    (carbon.benefit %>% dplyr::select(-ID)%>% filter(Scenario!="Real")),
#    by = c("x", "y", "Region", "Scenario")
#  )
#
##rm(list= ls()[!(ls() %in% c("pgm.biodiversity.benefit", "stm.biodiversity.benefit", "biodiversity.benefit", 
##                            "pgm.carbon.benefit", "stm.carbon.benefit", "carbon.benefit", 
##                            "pgm.costs", "stm.costs", "costs", "benefit.cost.ratio.df"))])
##gc()
##
#
#
#benefit.cost.ratio.df <- benefit.cost.ratio.df %>% 
#  mutate(Biodiversity.BCR = ifelse(is.na(Costs) | is.na(Biodiversity), NA, Biodiversity/Costs),
#         Carbon.BCR = ifelse(is.na(Costs) | is.na(Carbon), NA, Carbon/Costs))
#
##benefit.cost.ratio.df <- read.csv("data/benefit_cost_ratio.csv")
#
#
#benefit.cost.ratio.overview <- benefit.cost.ratio.df %>% group_by(Region, Scenario) %>% 
#  summarise(
#    Biodiversity_mean = mean(Biodiversity, na.rm=T),
#    Carbon_mean = mean(Carbon, na.rm=T),
#    Costs_mean = mean(Costs, na.rm=T),
#    Biodiversity_sd = sd(Biodiversity, na.rm=T),
#    Carbon_sd = sd(Carbon, na.rm=T),
#    Costs_sd = sd(Costs, na.rm=T),
#  )
#
#
#
#library(ggrepel)
#
#empty_theme <- theme(                              
#  plot.background = element_blank(), 
#  panel.grid.major = element_blank(), 
#  panel.grid.minor = element_blank(), 
#  panel.border = element_blank(), 
#  panel.background = element_blank(),
#  axis.line = element_blank(),
#  axis.ticks = element_blank(),
#  axis.text.y = element_text(angle = 90)
#)
#
#
#cowplot::plot_grid(
#
#ggplot(data = benefit.cost.ratio.overview, 
#       aes(x = Costs_mean, y = Biodiversity_mean, 
#           fill = Scenario, 
#           shape = Scenario,
#           label = Scenario)) + 
#  geom_point(aes(size = 8)) +
#  #geom_errorbar(aes(ymin = Biodiversity_mean-Biodiversity_sd, ymax = Biodiversity_mean+Biodiversity_sd)) + 
#  #geom_errorbarh(aes(xmin = Costs_mean-Costs_sd, xmax = Costs_mean+Costs_sd)) +
#  geom_hline(data=filter(benefit.cost.ratio.df, Region=="PGM"), 
#             aes(yintercept=163.0597), linetype="dashed") + 
#  geom_hline(data=filter(benefit.cost.ratio.df, Region=="STM"), 
#             aes(yintercept=201.396), linetype="dashed") +
#  geom_vline(data=filter(benefit.cost.ratio.df, Region=="PGM"), 
#             aes(xintercept=225.7966), linetype="dashed") + 
#  geom_vline(data=filter(benefit.cost.ratio.df, Region=="STM"), 
#             aes(xintercept=329.7966), linetype="dashed") +
#  scale_x_continuous(breaks = c(160,400), 
#                     labels=c("160" = "Low", "400" = "High")) +
#  scale_y_continuous(breaks = c(160,210), 
#                     labels=c("160" = "Low", "210" = "High")) +
#  scale_shape_manual(values = c(21, 21, 21, 24, 24, 24)) +
#  scale_fill_manual(values = c("#41644A", "#41644A", "#41644A", 
#                               "#D3756B", "#D3756B", "#D3756B")) +
#  geom_label_repel(size = 5,
#                   fill = NA,
#                   label.size = NA,
#                   min.segment.length = 0,
#                   segment.size = .5,
#                   segment.color=NA,
#                   box.padding = unit(.25, "lines"),
#                   nudge_y = 1.0E-6) +
#  facet_grid(~Region)+
#  labs(title = "Biodiversity Benefit x Cost",
#       x = "Costs",
#       y = "Benefits") +
#  empty_theme +
#  theme(legend.position = "none"),
#
#
#
#ggplot(data = benefit.cost.ratio.overview, 
#       aes(x = Costs_mean, y = Carbon_mean, 
#           fill = Scenario, 
#           shape = Scenario,
#           label = Scenario)) + 
#  geom_point(aes(size = 8)) +
#  #geom_errorbar(aes(ymin = Carbon_mean-Carbon_sd, ymax = Carbon_mean+Carbon_sd)) + 
#  #geom_errorbarh(aes(xmin = Costs_mean-Costs_sd, xmax = Costs_mean+Costs_sd)) +
#  geom_hline(data=filter(benefit.cost.ratio.df, Region=="PGM"), 
#             aes(yintercept=85.04721), linetype="dashed") + 
#  geom_hline(data=filter(benefit.cost.ratio.df, Region=="STM"), 
#             aes(yintercept=80.79686), linetype="dashed") +
#  geom_vline(data=filter(benefit.cost.ratio.df, Region=="PGM"), 
#             aes(xintercept=225.7966), linetype="dashed") + 
#  geom_vline(data=filter(benefit.cost.ratio.df, Region=="STM"), 
#             aes(xintercept=329.7966), linetype="dashed") +
#  scale_x_continuous(breaks = c(160,400), 
#                     labels=c("160" = "Low", "400" = "High")) +
#  scale_y_continuous(breaks = c(81,90), 
#                     labels=c("81" = "Low", "90" = "High")) +
#  scale_shape_manual(values = c(21, 21, 21, 24, 24, 24)) +
#  scale_fill_manual(values = c("#41644A", "#41644A", "#41644A", 
#                               "#D3756B", "#D3756B", "#D3756B")) +
#  geom_label_repel(size = 5,
#                   fill = NA,
#                   label.size = NA,
#                   min.segment.length = 0,
#                   segment.size = .5,
#                   segment.color=NA,
#                   box.padding = unit(.25, "lines"),
#                   nudge_y = 1.0E-6) +
#  facet_grid(~Region)+
#  labs(title = "Carbon Benefit x Cost",
#       x = "Costs",
#       y = "Benefits") +
#  empty_theme +
#  theme(legend.position = "none"),
#
#ncol = 1, align = "hv")
#  
#
#
#########################################
####         comparing costs         ####
####      between scenarios and      ####
####         real projections        ####
####  (conservation action areas)    ####
#########################################
pgm.conservact.avoiddegrad.costmask <- mask(pgm.costs[["pgm_avoiddegrad_cost"]], pgm.conservationaction.avoiddegrad)
pgm.conservact.avoiddeforest.costmask <- mask(pgm.costs[["pgm_avoiddeforest_cost"]], pgm.conservationaction.avoiddeforest)
pgm.conservact.avoidboth.costmask <- mask(pgm.costs[["pgm_avoidboth_cost"]], pgm.conservationaction.avoidboth)
pgm.conservact.restor_wo_avoid.costmask <- mask(pgm.costs[["pgm_restor_wo_avoid_cost"]], pgm.conservationaction.restor_wo_avoid)
pgm.conservact.restor_n_avoid_deforest.costmask <- mask(pgm.costs[["pgm_restor_n_avoid_deforest_cost"]], pgm.conservationaction.restor_n_avoid_deforest)
pgm.conservact.restor_n_avoid_both.costmask <- mask(pgm.costs[["pgm_restor_n_avoid_both_cost"]], pgm.conservationaction.restor_n_avoid_both)

pgm.conservact.costmask <- stack(pgm.conservact.avoiddegrad.costmask,
                                 pgm.conservact.avoiddeforest.costmask,
                                 pgm.conservact.avoidboth.costmask,
                                 pgm.conservact.restor_wo_avoid.costmask,
                                 pgm.conservact.restor_n_avoid_deforest.costmask,
                                 pgm.conservact.restor_n_avoid_both.costmask)

names(pgm.conservact.costmask) <- c("pgm_avoiddegrad_costmask",
                                    "pgm_avoiddeforest_costmask",
                                    "pgm_avoidboth_costmask",
                                    "pgm_restor_wo_avoid_costmask",
                                    "pgm_restor_n_avoid_deforest_costmask",
                                    "pgm_restor_n_avoid_both_costmask")

pgm.conservact.costmask.df <- as.data.frame(pgm.conservact.costmask, xy = TRUE)
pgm.conservact.costmask.df <- pgm.conservact.costmask.df %>% 
  pivot_longer(
    pgm_avoiddegrad_costmask:pgm_restor_n_avoid_both_costmask,
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


stm.conservact.avoiddegrad.costmask <- mask(stm.costs[["stm_avoiddegrad_cost"]], stm.conservationaction.avoiddegrad)
stm.conservact.avoiddeforest.costmask <- mask(stm.costs[["stm_avoiddeforest_cost"]], stm.conservationaction.avoiddeforest)
stm.conservact.avoidboth.costmask <- mask(stm.costs[["stm_avoidboth_cost"]], stm.conservationaction.avoidboth)
stm.conservact.restor_wo_avoid.costmask <- mask(stm.costs[["stm_restor_wo_avoid_cost"]], stm.conservationaction.restor_wo_avoid)
stm.conservact.restor_n_avoid_deforest.costmask <- mask(stm.costs[["stm_restor_n_avoid_deforest_cost"]], stm.conservationaction.restor_n_avoid_deforest)
stm.conservact.restor_n_avoid_both.costmask <- mask(stm.costs[["stm_restor_n_avoid_both_cost"]], stm.conservationaction.restor_n_avoid_both)

stm.conservact.costmask <- stack(stm.conservact.avoiddegrad.costmask,
                                 stm.conservact.avoiddeforest.costmask,
                                 stm.conservact.avoidboth.costmask,
                                 stm.conservact.restor_wo_avoid.costmask,
                                 stm.conservact.restor_n_avoid_deforest.costmask,
                                 stm.conservact.restor_n_avoid_both.costmask)

names(stm.conservact.costmask) <- c("stm_avoiddegrad_costmask",
                                    "stm_avoiddeforest_costmask",
                                    "stm_avoidboth_costmask",
                                    "stm_restor_wo_avoid_costmask",
                                    "stm_restor_n_avoid_deforest_costmask",
                                    "stm_restor_n_avoid_both_costmask")

stm.conservact.costmask.df <- as.data.frame(stm.conservact.costmask, xy = TRUE)
stm.conservact.costmask.df <- stm.conservact.costmask.df %>% 
  pivot_longer(
    stm_avoiddegrad_costmask:stm_restor_n_avoid_both_costmask,
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
                                str_detect(ID, "restor_n_avoid_both")~ "Restoration and avoid both",
                                str_detect(ID, "real")~ "Real"),
                      levels = c("Avoid degradation",
                                 "Avoid deforestation",
                                 "Avoid both",
                                 "Restoration without avoid",
                                 "Restoration and avoid deforestation",
                                 "Restoration and avoid both",
                                 "Real")),
    Costs = na_if(Costs, 0)
  )


conservact.costmask.df <- rbind(pgm.conservact.costmask.df, stm.conservact.costmask.df)

conservact.costmask.df %>% group_by(Region, Scenario) %>% drop_na() %>% 
  summarise(nvalues=n(), mean.benefit = mean(Costs, na.rm=T))

conservact.cost.overview <- conservact.costmask.df %>% group_by(Region, Scenario) %>% drop_na() %>%
  sample_n(size = 20000, replace = T) %>% 
  summarise(
    Costs_mean = mean(Costs, na.rm=T),
    Costs_sd = sd(Costs, na.rm=T),
    Costs_min = round(min(Costs, na.rm=T),1),
    Costs_max = round(max(Costs, na.rm=T),1)
  ) %>% 
  ungroup()




conservact.cost.overview %>%
  pivot_longer(Costs_mean:Costs_max,
               names_to = c("set", ".value"),
               names_pattern = "(.+)_(.+)"
  ) %>% 
  ggplot(aes(x = Scenario, y = mean, fill = Region)) +
  geom_errorbar(aes(ymin = mean, ymax = mean + sd),
                width = 0.2, colour = "black",
                position = position_dodge(0.9)
  ) + 
  geom_bar(
    stat = "identity", position = position_dodge(width = .9), 
    color="black", alpha=.7#, show_guide=FALSE
  ) +
  scale_fill_manual(values = c("#440154", "#374f00")) +
  scale_x_discrete(breaks=unique(conservact.cost.overview$Scenario), 
                   labels=addline_format(c("Avoid degradation", 
                                           "Avoid deforestation", 
                                           "Avoid both",
                                           "Restoration without avoid", 
                                           "Restoration and avoid deforestation", 
                                           "Restoration and avoid both"))) + 
  geom_text(
    aes(y= mean + sd, x= Scenario, group = Region, 
        label = paste0("(", min, " - ", max, ")")), 
    position = position_dodge(width = .9), vjust = -1
  ) +
  labs(x="", y="Costs (R$ / ha / year)") +
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        text = element_text(size = 22, family = "sans"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 18))


#####################################################################
#### analyzing the benefit-cost ratio (conservation action area) ####
#####################################################################
#conservact.benefit.cost.ratio.df <- conservact.costmask.df %>% dplyr::select(-ID) %>% 
#  left_join(
#    (conservact.biodbenefitmask.df %>% dplyr::select(-ID) %>% filter(Scenario!="Real")),
#    by = c("x", "y", "Region", "Scenario")
#  )
#
#conservact.benefit.cost.ratio.df <- conservact.benefit.cost.ratio.df %>% 
#  left_join(
#    (conservact.carbbenefitmask.df %>% dplyr::select(-ID)%>% filter(Scenario!="Real")),
#    by = c("x", "y", "Region", "Scenario")
#  )
#
#
#
#conservact.benefit.cost.ratio.df <- conservact.benefit.cost.ratio.df %>% 
#  mutate(Biodiversity.BCR = ifelse(is.na(Costs) | is.na(Biodiversity), NA, Biodiversity/Costs),
#         Carbon.BCR = ifelse(is.na(Costs) | is.na(Carbon), NA, Carbon/Costs))
#
##conservact.benefit.cost.ratio.df <- read.csv("data/conservact_benefit_cost_ratio.csv")
#
#
#conservact.benefit.cost.ratio.overview <- conservact.benefit.cost.ratio.df %>% group_by(Region, Scenario) %>% 
#  summarise(
#    Biodiversity_mean = mean(Biodiversity, na.rm=T),
#    Carbon_mean = mean(Carbon, na.rm=T),
#    Costs_mean = mean(Costs, na.rm=T),
#    Biodiversity_sd = sd(Biodiversity, na.rm=T),
#    Carbon_sd = sd(Carbon, na.rm=T),
#    Costs_sd = sd(Costs, na.rm=T),
#  )
#
#
#
#library(ggrepel)
#
#empty_theme <- theme(                              
#  plot.background = element_blank(), 
#  panel.grid.major = element_blank(), 
#  panel.grid.minor = element_blank(), 
#  panel.border = element_blank(), 
#  panel.background = element_blank(),
#  axis.line = element_blank(),
#  axis.ticks = element_blank(),
#  axis.text.y = element_text(angle = 90)
#)
#
#
#cowplot::plot_grid(
#  
#  ggplot(data = conservact.benefit.cost.ratio.overview, 
#         aes(x = Costs_mean, y = Biodiversity_mean, 
#             fill = Scenario, 
#             shape = Scenario,
#             label = Scenario)) + 
#    geom_point(aes(size = 8)) +
#    #geom_errorbar(aes(ymin = Biodiversity_mean-Biodiversity_sd, ymax = Biodiversity_mean+Biodiversity_sd)) + 
#    #geom_errorbarh(aes(xmin = Costs_mean-Costs_sd, xmax = Costs_mean+Costs_sd)) +
#    #geom_hline(data=filter(conservact.benefit.cost.ratio.overview, Region=="PGM"), 
#    #           aes(yintercept=150), linetype="dashed") + 
#    #geom_hline(data=filter(conservact.benefit.cost.ratio.overview, Region=="STM"), 
#    #           aes(yintercept=150), linetype="dashed") +
#    #geom_vline(data=filter(conservact.benefit.cost.ratio.overview, Region=="PGM"), 
#    #           aes(xintercept=625), linetype="dashed") + 
#    #geom_vline(data=filter(conservact.benefit.cost.ratio.overview, Region=="STM"), 
#    #           aes(xintercept=750), linetype="dashed") +
#    #scale_x_continuous(breaks = c(160,400), 
#    #                   labels=c("160" = "Low", "400" = "High")) +
#    #scale_y_continuous(breaks = c(160,210), 
#    #                   labels=c("160" = "Low", "210" = "High")) +
#    scale_shape_manual(values = c(21, 21, 21, 24, 24, 24)) +
#    scale_fill_manual(values = c("#41644A", "#41644A", "#41644A", 
#                                 "#D3756B", "#D3756B", "#D3756B")) +
#    geom_label_repel(size = 5,
#                     fill = NA,
#                     label.size = NA,
#                     min.segment.length = 0,
#                     segment.size = .5,
#                     segment.color=NA,
#                     box.padding = unit(.25, "lines"),
#                     nudge_y = 1.0E-6) +
#    facet_grid(~Region)+
#    labs(title = "Biodiversity Benefit x Cost",
#         x = "Costs",
#         y = "Benefits") +
#    empty_theme +
#    theme(legend.position = "none"),
#  
#  
#  
#  ggplot(data = conservact.benefit.cost.ratio.overview, 
#         aes(x = Costs_mean, y = Carbon_mean, 
#             fill = Scenario, 
#             shape = Scenario,
#             label = Scenario)) + 
#    geom_point(aes(size = 8)) +
#    #geom_errorbar(aes(ymin = Carbon_mean-Carbon_sd, ymax = Carbon_mean+Carbon_sd)) + 
#    #geom_errorbarh(aes(xmin = Costs_mean-Costs_sd, xmax = Costs_mean+Costs_sd)) +
#    #geom_hline(data=filter(conservact.benefit.cost.ratio.overview, Region=="PGM"), 
#    #           aes(yintercept=85.04721), linetype="dashed") + 
#    #geom_hline(data=filter(conservact.benefit.cost.ratio.overview, Region=="STM"), 
#    #           aes(yintercept=80.79686), linetype="dashed") +
#    #geom_vline(data=filter(conservact.benefit.cost.ratio.overview, Region=="PGM"), 
#    #           aes(xintercept=225.7966), linetype="dashed") + 
#    #geom_vline(data=filter(conservact.benefit.cost.ratio.overview, Region=="STM"), 
#    #           aes(xintercept=329.7966), linetype="dashed") +
#    #scale_x_continuous(breaks = c(160,400), 
#    #                   labels=c("160" = "Low", "400" = "High")) +
#    #scale_y_continuous(breaks = c(81,90), 
#    #                   labels=c("81" = "Low", "90" = "High")) +
#    scale_shape_manual(values = c(21, 21, 21, 24, 24, 24)) +
#    scale_fill_manual(values = c("#41644A", "#41644A", "#41644A", 
#                                 "#D3756B", "#D3756B", "#D3756B")) +
#    geom_label_repel(size = 5,
#                     fill = NA,
#                     label.size = NA,
#                     min.segment.length = 0,
#                     segment.size = .5,
#                     segment.color=NA,
#                     box.padding = unit(.25, "lines"),
#                     nudge_y = 1.0E-6) +
#    facet_grid(~Region)+
#    labs(title = "Carbon Benefit x Cost",
#         x = "Costs",
#         y = "Benefits") +
#    empty_theme +
#    theme(legend.position = "none"),
#  
#  ncol = 1, align = "hv")
#
#
#
#



#################################################
#### the optimal mix of conservation actions ####
#################################################
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

















  
  


























##next steps
#
#
#
#
##step 3 - mask with costs
#
###biodiversity
#pgm.conservact.avoiddegrad.biodcostmask <- mask(pgm.costs[["pgm_avoiddegrad_cost"]], pgm.conservact.avoiddegrad.biodbenefitmask)
#pgm.conservact.avoiddeforest.biodcostmask <- mask(pgm.costs[["pgm_avoiddeforest_cost"]], pgm.conservact.avoiddeforest.biodbenefitmask)
#pgm.conservact.avoidboth.biodcostmask <- mask(pgm.costs[["pgm_avoidboth_cost"]], pgm.conservact.avoidboth.biodbenefitmask)
#pgm.conservact.restor_wo_avoid.biodcostmask <- mask(pgm.costs[["pgm_restor_wo_avoid_cost"]], pgm.conservact.restor_wo_avoid.biodbenefitmask)
#pgm.conservact.restor_n_avoid_deforest.biodcostmask <- mask(pgm.costs[["pgm_restor_n_avoid_deforest.cost"]], pgm.conservact.restor_n_avoid_deforest.biodbenefitmask)
#pgm.conservact.restor_n_avoid_both.biodcostmask <- mask(pgm.costs[["pgm_restor_n_avoid_both_cost"]], pgm.conservact.restor_n_avoid_both.biodbenefitmask)
#
#pgm.conservact.biodcostmask <- stack(pgm.conservact.avoiddegrad.biodcostmask,
#                                     pgm.conservact.avoiddeforest.biodcostmask,
#                                     pgm.conservact.avoidboth.biodcostmask,
#                                     pgm.conservact.restor_wo_avoid.biodcostmask,
#                                     pgm.conservact.restor_n_avoid_deforest.biodcostmask,
#                                     pgm.conservact.restor_n_avoid_both.biodcostmask)
#
#names(pgm.conservact.carbbenefitmask) <- c("pgm_avoiddegrad_biodcostmask",
#                                           "pgm_avoiddeforest_biodcostmask",
#                                           "pgm_avoidboth_biodcostmask",
#                                           "pgm_restor_wo_avoid_biodcostmask",
#                                           "pgm_restor_n_avoid_deforest_biodcostmask",
#                                           "pgm_restor_n_avoid_both_biodcostmask")
#
#
#
#
#stm.conservact.avoiddegrad.biodcostmask <- mask(stm.costs[["stm_avoiddegrad_cost"]], stm.conservact.avoiddegrad.biodbenefitmask)
#stm.conservact.avoiddeforest.biodcostmask <- mask(stm.costs[["stm_avoiddeforest_cost"]], stm.conservact.avoiddeforest.biodbenefitmask)
#stm.conservact.avoidboth.biodcostmask <- mask(stm.costs[["stm_avoidboth_cost"]], stm.conservact.avoidboth.biodbenefitmask)
#stm.conservact.restor_wo_avoid.biodcostmask <- mask(stm.costs[["stm_restor_wo_avoid_cost"]], stm.conservact.restor_wo_avoid.biodbenefitmask)
#stm.conservact.restor_n_avoid_deforest.biodcostmask <- mask(stm.costs[["stm_restor_n_avoid_deforest_cost"]], stm.conservact.restor_n_avoid_deforest.biodbenefitmask)
#stm.conservact.restor_n_avoid_both.biodcostmask <- mask(stm.costs[["stm_restor_n_avoid_both_cost"]], stm.conservact.restor_n_avoid_both.biodbenefitmask)
#
#stm.conservact.biodcostmask <- stack(stm.conservact.avoiddegrad.biodcostmask,
#                                     stm.conservact.avoiddeforest.biodcostmask,
#                                     stm.conservact.avoidboth.biodcostmask,
#                                     stm.conservact.restor_wo_avoid.biodcostmask,
#                                     stm.conservact.restor_n_avoid_deforest.biodcostmask,
#                                     stm.conservact.restor_n_avoid_both.biodcostmask)
#
#names(stm.conservact.carbbenefitmask) <- c("stm_avoiddegrad_biodcostmask",
#                                           "stm_avoiddeforest_biodcostmask",
#                                           "stm_avoidboth_biodcostmask",
#                                           "stm_restor_wo_avoid_biodcostmask",
#                                           "stm_restor_n_avoid_deforest_biodcostmask",
#                                           "stm_restor_n_avoid_both_biodcostmask")
#
#
#
###carbon
#pgm.conservact.avoiddegrad.carbcostmask <- mask(pgm.costs[["pgm_avoiddegrad_cost"]], pgm.conservact.avoiddegrad.carbbenefitmask)
#pgm.conservact.avoiddeforest.carbcostmask <- mask(pgm.costs[["pgm_avoiddeforest_cost"]], pgm.conservact.avoiddeforest.carbbenefitmask)
#pgm.conservact.avoidboth.carbcostmask <- mask(pgm.costs[["pgm_avoidboth_cost"]], pgm.conservact.avoidboth.carbbenefitmask)
#pgm.conservact.restor_wo_avoid.carbcostmask <- mask(pgm.costs[["pgm_restor_wo_avoid_cost"]], pgm.conservact.restor_wo_avoid.carbbenefitmask)
#pgm.conservact.restor_n_avoid_deforest.carbcostmask <- mask(pgm.costs[["pgm_restor_n_avoid_deforest.cost"]], pgm.conservact.restor_n_avoid_deforest.carbbenefitmask)
#pgm.conservact.restor_n_avoid_both.carbcostmask <- mask(pgm.costs[["pgm_restor_n_avoid_both_cost"]], pgm.conservact.restor_n_avoid_both.carbbenefitmask)
#
#pgm.conservact.carbcostmask <- stack(pgm.conservact.avoiddegrad.carbcostmask,
#                                     pgm.conservact.avoiddeforest.carbcostmask,
#                                     pgm.conservact.avoidboth.carbcostmask,
#                                     pgm.conservact.restor_wo_avoid.carbcostmask,
#                                     pgm.conservact.restor_n_avoid_deforest.carbcostmask,
#                                     pgm.conservact.restor_n_avoid_both.carbcostmask)
#
#names(pgm.conservact.carbbenefitmask) <- c("pgm_avoiddegrad_carbcostmask",
#                                           "pgm_avoiddeforest_carbcostmask",
#                                           "pgm_avoidboth_carbcostmask",
#                                           "pgm_restor_wo_avoid_carbcostmask",
#                                           "pgm_restor_n_avoid_deforest_carbcostmask",
#                                           "pgm_restor_n_avoid_both_carbcostmask")
#
#
#
#stm.conservact.avoiddegrad.carbcostmask <- mask(stm.costs[["stm_avoiddegrad_cost"]], stm.conservact.avoiddegrad.carbbenefitmask)
#stm.conservact.avoiddeforest.carbcostmask <- mask(stm.costs[["stm_avoiddeforest_cost"]], stm.conservact.avoiddeforest.carbbenefitmask)
#stm.conservact.avoidboth.carbcostmask <- mask(stm.costs[["stm_avoidboth_cost"]], stm.conservact.avoidboth.carbbenefitmask)
#stm.conservact.restor_wo_avoid.carbcostmask <- mask(stm.costs[["stm_restor_wo_avoid_cost"]], stm.conservact.restor_wo_avoid.carbbenefitmask)
#stm.conservact.restor_n_avoid_deforest.carbcostmask <- mask(stm.costs[["stm_restor_n_avoid_deforest_cost"]], stm.conservact.restor_n_avoid_deforest.carbbenefitmask)
#stm.conservact.restor_n_avoid_both.carbcostmask <- mask(stm.costs[["stm_restor_n_avoid_both_cost"]], stm.conservact.restor_n_avoid_both.carbbenefitmask)
#
#stm.conservact.carbcostmask <- stack(stm.conservact.avoiddegrad.carbcostmask,
#                                     stm.conservact.avoiddeforest.carbcostmask,
#                                     stm.conservact.avoidboth.carbcostmask,
#                                     stm.conservact.restor_wo_avoid.carbcostmask,
#                                     stm.conservact.restor_n_avoid_deforest.carbcostmask,
#                                     stm.conservact.restor_n_avoid_both.carbcostmask)
#
#names(stm.conservact.carbbenefitmask) <- c("stm_avoiddegrad_carbcostmask",
#                                           "stm_avoiddeforest_carbcostmask",
#                                           "stm_avoidboth_carbcostmask",
#                                           "stm_restor_wo_avoid_carbcostmask",
#                                           "stm_restor_n_avoid_deforest_carbcostmask",
#                                           "stm_restor_n_avoid_both_carbcostmask")
#
#
##step 4 - choose cells ranked by costs until a budget constraint
##and calculate the proportion of land of each scenario in each simulation
##[simulations!]
#
#constraint.sim <- data.frame(Regiao=NA, Scenario=NA, Benefit=NA, N_cells=NA, Area=NA, Constraint=NA)
#
###biodiversity
#r <- pgm.conservact.biodcostmask[[c(1,2,4)]]
##plot(r)
##values(r)[1,]
#
#r.ord <- calc(r, fun=function(x, na.rm) x[order(x, decreasing = F)])
##plot(r.ord)
##values(r.ord)[1,]
#
#r.scen <- calc(r, function(x, na.rm) order(x, decreasing = F))
#r.scen <- mask(r.scen, r.ord[[1]])
##plot(r.scen)
##values(r.scen)[1,]
#table(values(r.scen[[1]]))
#
#ss <- sort(as.vector(r.ord[[1]]), decreasing = F)
#
#for (i in seq(0,100000000,2000000)) {
#  
#  s_lim <- ss[cumsum(ss) <= i]
#  
#  constraint <- r.ord[[1]] 
#  constraint[!constraint %in% s_lim ] <- NA
#  #cellStats(constraint, sum)
#  #plot(constraint)
#  
#  constraint_scen <- mask(r.scen[[1]], constraint)
#  #plot(constraint_scen)
#  constraint_df <- as.data.frame(table(values(constraint_scen)))
#  colnames(constraint_df) <- c("Scenario", "N_cells")
#  constraint_df$Regiao <- "PGM"
#  constraint_df$Benefit <- "Biodiversity"
#  constraint_df$Area <- constraint_df$N_cells/length(r.scen[[1]][Which(!is.na(r.scen[[1]]))])
#  constraint_df$Constraint <- i
#  
#  
#  constraint.sim <- rbind(constraint.sim, constraint_df)
#  
#}
#
#constraint.sim <- constraint.sim[-1,]
#
#
#r2 <- stm.conservact.biodcostmask[[c(1,2,4)]]
##plot(r2)
##values(r2)[1,]
#
#r2.ord <- calc(r2, fun=function(x, na.rm) x[order(x, decreasing = F)])
##plot(r2.ord)
##values(r2.ord)[1,]
#
#r2.scen <- calc(r2, function(x, na.rm) order(x, decreasing = F))
#r2.scen <- mask(r2.scen, r2.ord[[1]])
##plot(r2.scen)
##values(r2.scen)[1,]
#table(values(r2.scen[[1]]))
#
#ss2 <- sort(as.vector(r2.ord[[1]]), decreasing = F)
#
#for (i in seq(0,100000000,2000000)) {
#  
#  s2_lim <- ss2[cumsum(ss2) <= i]
#  
#  constraint <- r2.ord[[1]] 
#  constraint[!constraint %in% s2_lim ] <- NA
#  #cellStats(constraint, sum)
#  #plot(constraint)
#  
#  constraint_scen <- mask(r2.scen[[1]], constraint)
#  #plot(constraint_scen)
#  constraint_df <- as.data.frame(table(values(constraint_scen)))
#  colnames(constraint_df) <- c("Scenario", "N_cells")
#  constraint_df$Regiao <- "STM"
#  constraint_df$Benefit <- "Biodiversity"
#  constraint_df$Area <- constraint_df$N_cells/length(r2.scen[[1]][Which(!is.na(r2.scen[[1]]))])
#  constraint_df$Constraint <- i
#  
#  
#  constraint.sim <- rbind(constraint.sim, constraint_df)
#  
#}
#
#
#r3 <- pgm.conservact.carbcostmask[[c(1,2,4)]]
##plot(r3)
##values(r3)[1,]
#
#r3.ord <- calc(r3, fun=function(x, na.rm) x[order(x, decreasing = F)])
##plot(r3.ord)
##values(r3.ord)[1,]
#
#r3.scen <- calc(r3, function(x, na.rm) order(x, decreasing = F))
#r3.scen <- mask(r3.scen, r3.ord[[1]])
##plot(r3.scen)
##values(r3.scen)[1,]
#table(values(r3.scen[[1]]))
#
#ss3 <- sort(as.vector(r3.ord[[1]]), decreasing = F)
#
#for (i in seq(0,100000000,2000000)) {
#  
#  s3_lim <- ss3[cumsum(ss3) <= i]
#  
#  constraint <- r3.ord[[1]] 
#  constraint[!constraint %in% s3_lim ] <- NA
#  #cellStats(constraint, sum)
#  #plot(constraint)
#  
#  constraint_scen <- mask(r3.scen[[1]], constraint)
#  #plot(constraint_scen)
#  constraint_df <- as.data.frame(table(values(constraint_scen)))
#  colnames(constraint_df) <- c("Scenario", "N_cells")
#  constraint_df$Regiao <- "PGM"
#  constraint_df$Benefit <- "Carbon"
#  constraint_df$Area <- constraint_df$N_cells/length(r3.scen[[1]][Which(!is.na(r3.scen[[1]]))])
#  constraint_df$Constraint <- i
#  
#  
#  constraint.sim <- rbind(constraint.sim, constraint_df)
#  
#}
#
#
#r4 <- stm.conservact.carbcostmask[[c(1,2,4)]]
##plot(r4)
##values(r4)[1,]
#
#r4.ord <- calc(r4, fun=function(x, na.rm) x[order(x, decreasing = F)])
##plot(r4.ord)
##values(r4.ord)[1,]
#
#r4.scen <- calc(r4, function(x, na.rm) order(x, decreasing = F))
#r4.scen <- mask(r4.scen, r4.ord[[1]])
##plot(r4.scen)
##values(r4.scen)[1,]
#table(values(r4.scen[[1]]))
#
#ss4 <- sort(as.vector(r4.ord[[1]]), decreasing = F)
#
#for (i in seq(0,100000000,2000000)) {
#  
#  s4_lim <- ss4[cumsum(ss4) <= i]
#  
#  constraint <- r4.ord[[1]] 
#  constraint[!constraint %in% s4_lim ] <- NA
#  #cellStats(constraint, sum)
#  #plot(constraint)
#  
#  constraint_scen <- mask(r4.scen[[1]], constraint)
#  #plot(constraint_scen)
#  constraint_df <- as.data.frame(table(values(constraint_scen)))
#  colnames(constraint_df) <- c("Scenario", "N_cells")
#  constraint_df$Regiao <- "STM"
#  constraint_df$Benefit <- "Carbon"
#  constraint_df$Area <- constraint_df$N_cells/length(r4.scen[[1]][Which(!is.na(r4.scen[[1]]))])
#  constraint_df$Constraint <- i
#  
#  
#  constraint.sim <- rbind(constraint.sim, constraint_df)
#  
#}
#
#
#
#constraint.sim$Scenario <- ifelse(constraint.sim$Scenario == 1, "Degradation", 
#                                  ifelse(constraint.sim$Scenario == 2, "Deforestation", "Restoration"))
#
#
#
#constraint.sim %>% mutate(Scenario = factor(Scenario, levels = c("Degradation", 
#                                                                 "Deforestation", 
#                                                                 "Restoration"))) %>% 
#  ggplot(aes(x=Constraint, y=Area, color=Scenario))+
#  geom_smooth(se=F) +
#  scale_x_continuous(breaks = seq(0,100000000,10000000),
#                     labels = unit_format(unit = "M", scale = 1e-6)) +
#  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
#  facet_grid(Benefit~Regiao) +
#  labs(x="Budget constraint", y="Proportional landscape area")+
#  theme_bw()+
#  theme(plot.title = element_text(size = 22, family = "sans", face = "bold"),
#        text = element_text(size = 22, family = "sans"),
#        axis.title = element_text(face="bold"),
#        axis.text.x=element_text(angle = 45, hjust=1, size = 20),
#        panel.spacing = unit(2, "lines"),
#        legend.position = "top")
#
#
#
##
























##### Global (simple) benefits & costs bar plot ####
### benefits
#conservact.benefit.cost.ratio.overview2 <- conservact.benefit.cost.ratio.df %>% group_by(Scenario) %>% 
#  summarise(
#    Biodiversity_mean = mean(Biodiversity, na.rm=T),
#    Carbon_mean = mean(Carbon, na.rm=T),
#    Costs_mean = mean(Costs, na.rm=T),
#    Biodiversity_sd = sd(Biodiversity, na.rm=T),
#    Carbon_sd = sd(Carbon, na.rm=T),
#    Costs_sd = sd(Costs, na.rm=T),
#  ) %>% 
#  ungroup()
#
#
#
#
#addline_format <- function(x,...){
#  gsub('\\s','\n',x)
#}
#
#
#conservact.benefit.cost.ratio.overview2 %>%
#  pivot_longer(Biodiversity_mean:Costs_sd,
#               names_to = c("set", ".value"),
#               names_pattern = "(.+)_(.+)"
#  ) %>% 
#  filter(set != "Costs") %>%
#  ggplot(aes(x = Scenario, y = mean, fill = set)) +
#  geom_bar(
#    stat = "identity", position = position_dodge(), alpha=.7
#  ) +
#  scale_fill_manual(values = c("#440154", "#374f00")) +
#  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
#                width = 0.2, colour = "black",
#                position = position_dodge(0.9)
#  ) + 
#  scale_x_discrete(breaks=unique(conservact.benefit.cost.ratio.overview2$Scenario), 
#                   labels=addline_format(c("Avoid degradation", 
#                                           "Avoid deforestation", 
#                                           "Avoid both",
#                                           "Restoration without avoid", 
#                                           "Restoration and avoid deforestation", 
#                                           "Restoration and avoid both"))) +
#  labs(x="", y="Benefits") +
#  theme_bw()+
#  theme(legend.position = "top",
#        legend.title = element_blank(),
#        text = element_text(size = 22, family = "sans"),
#        axis.title = element_text(face="bold"),
#        axis.text.x=element_text(size = 18))
#
#
#
#
#
#
#conservact.benefit.cost.ratio.overview2 %>%
#  pivot_longer(Biodiversity_mean:Costs_sd,
#               names_to = c("set", ".value"),
#               names_pattern = "(.+)_(.+)"
#  ) %>% 
#  filter(set == "Costs") %>%
#  ggplot(aes(x = Scenario, y = mean)) +
#  geom_bar(
#    stat = "identity", alpha=.7
#  ) +
#  scale_fill_manual(values = "grey25") +
#  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
#                width = 0.2, colour = "grey25"
#  ) + 
#  scale_x_discrete(breaks=unique(conservact.benefit.cost.ratio.overview2$Scenario), 
#                   labels=addline_format(c("Avoid degradation", 
#                                           "Avoid deforestation", 
#                                           "Avoid both",
#                                           "Restoration without avoid", 
#                                           "Restoration and avoid deforestation", 
#                                           "Restoration and avoid both"))) +
#  labs(x="", y="Costs") +
#  theme_bw()+
#  theme(legend.position = "top",
#        legend.title = element_blank(),
#        text = element_text(size = 22, family = "sans"),
#        axis.title = element_text(face="bold"),
#        axis.text.x=element_text(size = 18))
#
#
