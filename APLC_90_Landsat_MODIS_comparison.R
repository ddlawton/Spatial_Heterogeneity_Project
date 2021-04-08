#####
# APLC
#   90
#    Landsat/MODIS dat
#     comparison
#####
rm(list=ls())
library(tidyverse)
library(raster)
library(rgdal)
library(spatialEco)
library(sf)
library(scales)
library(ggpubr)
library(gridExtra)
library(viridis)
library(lubridate)
library(zoo)
library(data.table)
library(patchwork)

modis <- fread("Data/processed/APLC_90day_MODIS_April62021.csv")
landsat <- fread("Data/processed/APLC_90day_Landsat_April72021.csv")
names(modis)
CT_modis <- modis %>%
  filter(
    Species == '10' | Species == '11'
  )

CT_landsat <- landsat %>%
  filter(
    Species == '10' | Species == '11'
  )

CT_modis$prominence <- as.numeric(CT_modis$prominence)
CT_landsat$prominence <- as.numeric(CT_landsat$prominence)


modis_ndvi <- ggplot(CT_modis, aes(x=diff_days,y=NDVIs/1000, color=as.factor(binary_outbreak))) + geom_smooth() + xlab("Days (0 = observation date)") +
  scale_color_manual(name = "", labels = c("Non-outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62")) + theme_pubr() + xlim(-90,0)  + 
  coord_cartesian(ylim=c(0, 0.5)) +
  ggtitle("MODIS")


landsat_ndvi <- ggplot(CT_landsat, aes(x=diff_days,y=NDVIs/1000, color=as.factor(binary_outbreak))) + geom_smooth() + xlab("Days (0 = observation date)") +
  scale_color_manual(name = "", labels = c("Non-outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62")) + theme_pubr() + xlim(-90,0) + 
  coord_cartesian(ylim=c(0, 0.5)) +
  ggtitle("Landsat")

modis_ndvi + landsat_ndvi


str(CT_landsat)
modis_prominence <- ggplot(CT_modis, aes(x=diff_days,y=prominence, color=as.factor(binary_outbreak))) + geom_smooth(method="gam") + xlab("Days (0 = observation date)") +
  scale_color_manual(name = "", labels = c("Non-outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62")) + theme_pubr() + xlim(-90,0)  + 
  coord_cartesian(ylim=c(-2e+14, 2e+14)) +
  ggtitle("MODIS")


landsat_prominence <- ggplot(CT_landsat, aes(x=diff_days,y=prominence, color=as.factor(binary_outbreak))) + geom_smooth(method="gam") + xlab("Days (0 = observation date)") +
  scale_color_manual(name = "", labels = c("Non-outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62")) + theme_pubr() + xlim(-90,0) + 
  coord_cartesian(ylim=c(-2e+14, 2e+14)) +
  ggtitle("Landsat")

modis_entropy <- ggplot(CT_modis, aes(x=diff_days,y=spatial_hetero, color=as.factor(binary_outbreak))) + geom_smooth(method="gam") + xlab("Days (0 = observation date)") +
  scale_color_manual(name = "", labels = c("Non-outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62")) + theme_pubr() + xlim(-90,0)  + 
  coord_cartesian(ylim=c(0,1)) +
  ggtitle("MODIS")


landsat_entropy <- ggplot(CT_landsat, aes(x=diff_days,y=entropy, color=as.factor(binary_outbreak))) + geom_smooth(method="gam") + xlab("Days (0 = observation date)") +
  scale_color_manual(name = "", labels = c("Non-outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62")) + theme_pubr() + xlim(-90,0) + 
  coord_cartesian(ylim=c(0, 5)) +
  ggtitle("Landsat")

modis_prominence + landsat_prominence + modis_entropy + landsat_entropy
names(CT_modis)
CT_modis <- CT_modis %>%
  mutate(
    contrasts_scaled = scale(contrasts,center = TRUE, scale=TRUE),
    dissimilarity_scaled = scale(dissimilarity,center = TRUE, scale=TRUE),
    entropy_scaled = scale(entropy,center = TRUE, scale=TRUE),
    prominence_scaled = scale(CT_modis$prominence,center = TRUE, scale=TRUE)) 

CT_modis <- CT_modis %>% mutate(
  contrasts_normalized = (CT_modis$contrasts_scaled - min(CT_modis$contrasts_scaled))/(max(CT_modis$contrasts_scaled)-min(CT_modis$contrasts_scaled)),
  dissimilarity_normalized = (CT_modis$dissimilarity_scaled - min(CT_modis$dissimilarity_scaled))/(max(CT_modis$dissimilarity_scaled)-min(CT_modis$dissimilarity_scaled)),
  entropy_normalized = (CT_modis$entropy_scaled - min(CT_modis$entropy_scaled))/(max(CT_modis$entropy_scaled)-min(CT_modis$entropy_scaled)),
  prominence_normalized = (CT_modis$prominence_scaled - min(CT_modis$prominence_scaled))/(max(CT_modis$prominence_scaled)-min(CT_modis$prominence_scaled))) 

names(CT_modis)

CT_modis <- CT_modis %>% mutate(
  spatial_hetero = rowMeans(.[,35:38], na.rm = TRUE)
)


ggplot(CT_modis,aes(x=diff))