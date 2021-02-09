######
# Australian plague locust --
#  Hierarhical structure extraction and
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

dat <- fread("data/processed/CT_outbreaks_Feb_7_MODIS.csv")

str(dat)

dat$Date <- as.Date(dat$Date,format="%m/%d/%y")
dat$Source <- as.factor(dat$Source)
dat$Species <- as.factor(dat$Species)
dat$AdultDensity <- as.factor(dat$AdultDensity)
dat$AdultStage<- as.factor(dat$AdultStage)
dat$DataQuality<- as.factor(dat$DataQuality)
dat$NymphDensity<- as.factor(dat$NymphDensity)
dat$NymphStage<- as.factor(dat$NymphStage)

dens_target <- c("0","1","2","3")

dat2 <- dat %>%
  mutate(
    spp_code = case_when(
      Species == 10 ~ "None",
      Species == 11 ~ "CT",
      Species == 12 ~ "AG",
      Species == 13 ~ "GM",
      Species == 14 ~ "LM",
      Species == 15 ~ "AC",
      Species == 16 ~ "Other",
      Species == 17 ~ "AT",
      Species == 18 ~ "BS",
      Species == 19 ~ "OA",
      Species == 20 ~ "PS",
      Species == 21 ~ "US",
      Species == 22 ~ "VS"),
    binary_outbreak = case_when(
      NymphDensity %in% dens_target ~ "0",
      NymphDensity == "4" ~ "1"
    )
  )


## Extraction seasonal rainfall zones
#Loading and projecting data

names(dat)
spdf <- SpatialPointsDataFrame(coords=dat[,8:7],data=dat2,
                               proj4string = CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs "))


#Loading in raster
rainzones <- raster("~/Dropbox (ASU)/Research/Gradaute/Chapter 4 - Hierarchical scaling of locust swarms/Hierarchical_locust_modeling/Australian_plague_locust_model/Data/Raw/seasonal_rainfall_zones.tif")

crs(rainzones) <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs "

#Check if same projection
projection(rainzones) == projection(spdf)

#Extract data from raster
spdf <- raster::extract(rainzones,spdf,method="simple",df=TRUE,sp=TRUE)


## Extraction Bioregions

bioregions <- readOGR(dsn="~/Dropbox (ASU)/Research/Gradaute/Chapter 4 - Hierarchical scaling of locust swarms/Hierarchical_locust_modeling/Australian_plague_locust_model/Data/Raw/IBRA_bioregions", layer = "IBRA_bioregions")

crs(bioregions)

bioregions2 <- spTransform(bioregions,CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs "))


projection(bioregions2) == projection(spdf)

spdf2 <- point.in.poly(spdf,bioregions2)

summary(dat3$REG_NAME_7)
#Renaming zones
dat3 <- spdf2@data %>%
  mutate(Minor_rain_zones = case_when(
    seasonal_rainfall_zones == 3 ~ "Summer Dominant (350mm-650mm)",
    seasonal_rainfall_zones == 6 ~ "Summer (650mm-1200mm)",
    seasonal_rainfall_zones == 7 ~ "Summer (350mm-650mm)",
    seasonal_rainfall_zones == 10 ~ "Uniform (500mm-800mm)",
    seasonal_rainfall_zones == 11 ~ "Uniform (250mm-500mm)",
    seasonal_rainfall_zones == 17 ~ "Winter (>800mm)",
    seasonal_rainfall_zones == 18 ~ "Winter (500mm-800mm)",
    seasonal_rainfall_zones == 19 ~ "Winter (250mm-500mm)",
    seasonal_rainfall_zones == 4 ~ "Arid (<350mm)",
    seasonal_rainfall_zones == 8 ~ "Arid (<350mm)",
    seasonal_rainfall_zones == 12 ~ "Arid (<350mm)",
    seasonal_rainfall_zones == 20 ~ "Arid (<350mm)"
  )) %>%
  mutate(Major_rain_zones = case_when(
    seasonal_rainfall_zones == 3 ~ "Summer Dominant",
    seasonal_rainfall_zones == 6 ~ "Summer",
    seasonal_rainfall_zones == 7 ~ "Summer",
    seasonal_rainfall_zones == 10 ~ "Uniform",
    seasonal_rainfall_zones == 11 ~ "Uniform",
    seasonal_rainfall_zones == 17 ~ "Winter",
    seasonal_rainfall_zones == 18 ~ "Winter",
    seasonal_rainfall_zones == 19 ~ "Winter",
    seasonal_rainfall_zones == 4 ~ "Arid",
    seasonal_rainfall_zones == 8 ~ "Arid",
    seasonal_rainfall_zones == 12 ~ "Arid",
    seasonal_rainfall_zones == 20 ~ "Arid")) %>%  
  mutate(Source_code = case_when(
    Source == 1 ~ "APLC Survey",
    Source == 2 ~ "APLC Control",
    Source == 3 ~ "NSW Data",
    Source == 4 ~ "APLC aerial Survey",
    Source == 5 ~ "SA Data",
    Source == 7 ~ "NSW DPI aerial Survey",
    TRUE ~ "APLC Survey")) %>%
  dplyr::select(!seasonal_rainfall_zones) %>% 
  dplyr::select(!STATE_CODE)


dat4 <- cSplit(setDT(dat3)[, lapply(.SD, gsub, pattern = "[][}]", 
                                    replacement = "")], names(dat3), sep=",", direction='long', fixed = FALSE, "long")

dat4[dat4 == ""] <- NA # define NA pattern
dat5 <- dat4[rowSums(is.na(dat4)) != ncol(dat4), ]

dat6 <- dat5 %>% tidyr::fill(c("AdultDensity", "Date", "AdultStage","DataQuality", "Latitude", 
                               "Longitude","ID","NymphDensity","NymphStage","Source","spp_code","Source_code","Species",
                               "binary_outbreak","REG_NAME_7","Minor_rain_zones","Major_rain_zones","unixdates","system.index"), .direction = 'down') 

dat6$diff_days <- as.numeric(difftime(anytime(dat6$dates/1000), (dat6$Date), units = "days"))


yq <- as.yearqtr(as.yearmon(dat6$Date) + 1/12)
dat6$Season <- factor(format(yq, "%q"), levels = 1:4, 
                      labels = c("Summer", "Fall", "Winter", "Spring"))

write.csv(dat6,file="Data/processed/APLC_database_hierarchical_data_Feb_8.csv")

