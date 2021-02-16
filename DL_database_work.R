######
# Desert locust --
#  Hierarchical structure extraction and
#   general formatting for modeling
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
dat <- fread("data/raw/DL_MODIS_dat.csv")

names(dat)

dat$Date <- as.Date(dat$DATE,format="%m/%d/%y")
dat$SHPAPPGREG <- as.factor(dat$SHPAPPGREG)
dat$SHPAPPSOL <- as.factor(dat$SHPAPPSOL)
dat$SHPAPPTRAN <- as.factor(dat$SHPAPPTRAN)
dat$SHPAPPUNK<- as.factor(dat$SHPAPPUNK)
dat$SHPDENGRP<- as.factor(dat$SHPDENGRP)
dat$SHPDENISOL<- as.factor(dat$SHPDENISOL)
dat$SHPDENSCAT<- as.factor(dat$SHPDENSCAT)
dat$adult<- as.factor(dat$adult)
dat$band<- as.factor(dat$band)

dat$Month <- month(dat$Date)

dat2 <- dat %>%
  mutate(Season = case_when(
    Month == 1 ~ "Winter",
    Month == 2 ~ "Winter",
    Month == 3 ~ "Winter",
    Month == 4 ~ "Spring",
    Month == 5 ~ "Spring",
    Month == 6 ~ "Spring",
    Month == 7 ~ "Summer",
    Month == 8 ~ "Summer",
    Month == 9 ~ "Summer",
    Month == 10 ~ "Fall",
    Month ==11 ~ "Fall",
    Month == 12 ~ "Fall"
  ))

flevels  <- c("Winter","Spring","Summer","Fall")
dat2$Season <- factor(dat2$Season,levels=flevels)


# Now adding Recession/Invasion zones from Cyril

Recession_area <- readOGR(dsn = "~/Dropbox (ASU)/Research/Graduate/Chapter 4 - Hierarchical scaling of locust swarms/Hierarchical_locust_modeling/Desert_locust_model/Data/Raw/SGRmap", layer = "recessionAreaClean")
Invasion_area <- readOGR(dsn = "~/Dropbox (ASU)/Research/Graduate/Chapter 4 - Hierarchical scaling of locust swarms/Hierarchical_locust_modeling/Desert_locust_model/Data/Raw/SGRmap", layer = "invasionareaclean_dl")

areas <- bind(Recession_area,Invasion_area)


coords <- as.data.frame(cbind(dat$Longitude,dat$Latitude))

NDVI_projected <- SpatialPointsDataFrame(coords,data=dat,
                                         proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))



projection(areas) == projection(NDVI_projected) #check to see if projections are the same

test <- point.in.poly(NDVI_projected,areas)

test2 <- test@data

dat2 <- test2 %>% dplyr::select(-c(area))

dat2_projected <- SpatialPointsDataFrame(coords,data=dat2,
                                         proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

######
# Adding WWF ecoregions
###

ecoregions <- readOGR(dsn = "~/Dropbox (ASU)/Research/Graduate/Chapter 4 - Hierarchical scaling of locust swarms/Hierarchical_locust_modeling/Desert_locust_model/Data/Raw/SGRmap/WFF_ecoregions", layer = "wwf_terr_ecos")

projection(ecoregions) == projection(dat2_projected)


test <- point.in.poly(dat2_projected,ecoregions)
test2 <- test@data

names(test2)
dat <- test2 %>%
  dplyr::select(1:30,ECO_NAME)





dat4 <- cSplit(setDT(dat)[, lapply(.SD, gsub, pattern = "[][}]", 
                                    replacement = "")], names(dat), sep=",", direction='long', fixed = FALSE, "long")

dat4[dat4 == ""] <- NA # define NA pattern
dat5 <- dat4[rowSums(is.na(dat4)) != ncol(dat4), ]

names(dat5)
dat6 <- dat5 %>% tidyr::fill(!c("NDVIs","contrasts","correlation","dates","dissimilarity","entropy",
                                "homogeneity"), .direction = 'down') 

dat6$diff_days <- as.numeric(difftime(anytime(dat6$dates/1000), (dat6$Date), units = "days"))



write.csv(dat6,file="Data/processed/DL_database_hierarchical_data_Feb_15.csv")

