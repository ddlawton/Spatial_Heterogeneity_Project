######
# Desert locust --
#  Exporatory graphing
#   
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
library(splitstackshape)
library(anytime)
library(mgcv)
library(viridis)
library(gratia)
dat <- fread("Data/processed/DL_database_hierarchical_data_Feb_15.csv")
dat$NDVIs <- dat$NDVIs/1000

dat <- dat %>% drop_na(c(contrasts,correlation,dissimilarity,entropy,homogeneity,NDVIs)) %>% filter(between(NDVIs,0,1)) %>% filter(between(diff_days,-60,0))

dat$contrasts_normalized <- (dat$contrasts - min(dat$contrasts))/(max(dat$contrasts)-min(dat$contrasts))
#dat$correlation_normalized <- (dat$correlation - min(dat$correlation))/(max(dat$correlation)-min(dat$correlation))
dat$dissimilarity_normalized <- (dat$dissimilarity - min(dat$dissimilarity))/(max(dat$dissimilarity)-min(dat$dissimilarity))
dat$entropy_normalized <- (dat$entropy - min(dat$entropy))/(max(dat$entropy)-min(dat$entropy))
#dat$homogeneity_normalized <- (dat$homogeneity - min(dat$homogeneity))/(max(dat$homogeneity)-min(dat$homogeneity))


dat2 <- dat %>% mutate(
  spatial_hetero = rowMeans(dplyr::select(.,ends_with("_normalized")), na.rm = TRUE)) 

ggplot(dat2,aes(x=as.factor(outbreak),y=spatial_hetero)) + geom_boxplot()

ggplot(dat2,aes(x=NDVIs,y=spatial_hetero)) + geom_smooth()

ggplot(dat2,aes(x=diff_days,y=NDVIs,color=as.factor(outbreak))) +geom_smooth() + xlim(-60,0)

ggplot(dat2,aes(x=diff_days,y=dissimilarity,color=as.factor(outbreak))) +geom_smooth()  + 
  ylab("Spatial Heoterogenity")

ggplot(dat2,aes(x=diff_days,y=spatial_hetero,color=as.factor(outbreak))) +geom_smooth()  + 
  ylab("Spatial Heoterogenity") + facet_wrap(~Zone)

dat2$outbreak <- as.factor(dat2$outbreak)
dat2$Zone <- as.factor(dat2$Zone)
dat2$ECO_NAME <- as.factor(dat2$ECO_NAME)

mod <- bam((outbreak) ~ te(diff_days,NDVIs) + te(diff_days,spatial_hetero) + te(Longitude,Latitude) +
           s((Zone),bs="re") + s((ECO_NAME),bs="re"), dat=dat2,discrete = TRUE,nthreads = 7,
           family=binomial(),select=TRUE)

summary(mod)
draw(mod)
appraise(mod)
dat2$preds <- predict(mod,type="response")

ggplot(dat2,aes(x=diff_days,y=NDVIs,z=(preds))) + stat_summary_2d() + scale_fill_viridis() #+ ylim(0,.75)
ggplot(dat2,aes(x=diff_days,y=spatial_hetero,z=(preds))) + stat_summary_hex() + scale_fill_viridis() + facet_wrap(~Zone)

