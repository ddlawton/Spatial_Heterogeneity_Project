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
library(rnaturalearth)
library(scales)
library(spdplyr)

rawdat <- fread("Data/raw/DL_MODIS_dat.csv")
names(rawdat)

rawdat <- rawdat %>% dplyr::select("DATE","Latitude","Longitude","outbreak")


dat <- fread("Data/processed/DL_database_hierarchical_data_Feb_24.csv")
dat$NDVIs <- dat$NDVIs/1000

dat <- dat %>% drop_na(c(contrasts,correlation,dissimilarity,entropy,homogeneity,NDVIs)) %>% 
  filter(between(NDVIs,0,1)) %>% filter(between(diff_days,-78,0)) %>% 
  dplyr::select(!c("contrasts_normalized", "dissimilarity_normalized","entropy_normalized",
            "spatial_hetero"))


dat <- dat %>%
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
dat$Season <- factor(dat$Season,levels=flevels)
dat$Date <- as.Date(dat$Date)



cols <- c("SHPAPPGREG", "SHPAPPSOL", "SHPAPPTRAN", "SHPAPPUNK","SHPDENGRP","SHPDENISOL","SHPDENSCAT",
          "adult","band","outbreak","swarm","Zone","ECO_NAME","binary_outbreak")

for(col in cols){
  set(dat, j = col, value = as.factor(dat[[col]]))
}


dat$DOY <- yday(dat$Date)
dat$Year <- year(dat$Date)

names(dat)

dat <- dat %>%
  filter(between(diff_days,-60,0)) %>% filter(NDVIs >= 0) %>%
  mutate(
    contrasts_scaled = scale(contrasts,center = TRUE, scale=TRUE),
    correlation_scaled = scale(correlation,center = TRUE, scale=TRUE),
    dissimilarity_scaled = scale(dissimilarity,center = TRUE, scale=TRUE),
    entropy_scaled = scale(entropy,center = TRUE, scale=TRUE),
    homogeneity_scaled = scale(homogeneity,center = TRUE, scale=TRUE)) 

dat <- dat %>% mutate(
  contrasts_normalized = (dat$contrasts_scaled - min(dat$contrasts_scaled))/(max(dat$contrasts_scaled)-min(dat$contrasts_scaled)),
  correlation_normalized = (dat$correlation_scaled - min(dat$correlation_scaled))/(max(dat$correlation_scaled)-min(dat$correlation_scaled)),
  dissimilarity_normalized = (dat$dissimilarity_scaled - min(dat$dissimilarity_scaled))/(max(dat$dissimilarity_scaled)-min(dat$dissimilarity_scaled)),
  entropy_normalized = (dat$entropy_scaled - min(dat$entropy_scaled))/(max(dat$entropy_scaled)-min(dat$entropy_scaled)),
  homogeneity_normalized = (dat$homogeneity_scaled - min(dat$homogeneity_scaled))/(max(dat$homogeneity_scaled)-min(dat$homogeneity_scaled))) 
names(dat)

names(dat)

dat <- dat %>% mutate(
  spatial_hetero = rowMeans(.[,43:47], na.rm = TRUE)
)

ggplot(dat,aes(x=as.factor(outbreak),y=spatial_hetero)) + geom_boxplot()


ggplot(dat,aes(x=NDVIs,y=spatial_hetero)) + geom_smooth()

ggplot(dat,aes(x=diff_days,y=NDVIs,color=as.factor(outbreak))) +geom_smooth() + xlim(-60,0)

ggplot(dat,aes(x=diff_days,y=dissimilarity,color=as.factor(outbreak))) +geom_smooth()  + 
  ylab("Dissimilarity")

ggplot(dat,aes(x=diff_days,y=contrasts,color=as.factor(outbreak))) +geom_smooth()  + 
  ylab("contrasts")

ggplot(dat,aes(x=diff_days,y=entropy,color=as.factor(outbreak))) +geom_smooth()  + 
  ylab("Entropy")

ggplot(dat,aes(x=diff_days,y=spatial_hetero,color=as.factor(outbreak))) +geom_smooth()  + 
  ylab("Spatial Heoterogenity") + ylim(0,.5) + theme_pubclean()

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



#Bioregion breakdown

date_dat <- rawdat
date_dat$DATE <- as.Date(date_dat$DATE,format="%m/%d/%y")

ecoregions <- readOGR(dsn = "/Users/ddlawton/Dropbox (ASU)/Research/Graduate/Chapter 4 - Hierarchical scaling of locust swarms/Hierarchical_locust_modeling/Desert_locust_model/Data/Raw/SGRmap/WFF_ecoregions/",
                      layer = "wwf_terr_ecos")


names(date_dat)

coords <- cbind(date_dat$Longitude,date_dat$Latitude)
date_dat_projected <- SpatialPointsDataFrame(coords,data=date_dat,
                                         proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))


plot(date_dat_projected)
projection(ecoregions) == projection(date_dat_projected) #check to see if projections are the same


test <- point.in.poly(date_dat_projected,ecoregions)
test2 <- test@data
names(test2)

date_dat_projected2 <- test2 %>%
  dplyr::select(1:4,ECO_NAME)




tab <- dat %>% group_by(ECO_NAME,outbreak) %>% tally() %>% pivot_wider(names_from = outbreak,values_from = n) %>%
  mutate(total = `0`+`1`) %>%  filter(total > 200) %>% filter(`1` >= 100)


tab_names <- unique(tab$ECO_NAME)

names <- c("Mediterranean woodlands and forests","Mediterranean dry woodlands and steppe","Atlantic coastal desert",
           "Mediterranean acacia-argania dry woodlands and succulent thickets")

dat <- dat %>%filter((ECO_NAME %in% tab_names)) %>% filter(!(ECO_NAME %in% names)) %>% droplevels()

dat$ECO_NAME <- as.factor(dat$ECO_NAME)
# Filtering out unneeded ecoregions
used_eco_regions <- unique(dat$ECO_NAME) %>% droplevels()

filter_ecoregions <- ecoregions %>%
  filter(ECO_NAME %in% used_eco_regions)

filtered_ecoregions2 <- st_as_sf(filter_ecoregions)
world <- ne_countries(scale = "medium", returnclass = "sf")
names(dat)


eco_names <- unique(dat$ECO_NAME) %>% droplevels()
date_dat_projected2$outbreak <- as.factor(date_dat_projected2$outbreak)


eco <- list()
for(i in eco_names){
  print(paste("Now computing",i))
  ecorgn <- dat %>% filter(ECO_NAME == i)
  
  ob_summary <- date_dat_projected2 %>% 
    dplyr::filter(ECO_NAME == i) %>% 
    group_by(outbreak) %>% 
    tally()
  
  no_outbreaks <- paste0("Non-outbreak (N: ",(ob_summary %>% filter(outbreak=="0"))$n,")")
  outbreaks <- paste0("Outbreak (N: ",(ob_summary %>% filter(outbreak=="1"))$n,")")
  
  g1 <- ggplot(ecorgn,aes(x=diff_days,y=spatial_hetero,color=as.factor(outbreak))) +
    geom_smooth() + ylab("Heterogeneity") + xlab("") +
    xlab("Days (0 = date of outbreak)") +
    scale_color_discrete(name = "", labels = c("Non-outbreak", "Outbreak"))+
    geom_vline(xintercept=0,color="red",linetype=2)  + 
    geom_vline(xintercept=-35,color="black",linetype=2) + 
    geom_vline(xintercept=-21,color="black",linetype=2) + 
    scale_x_continuous(limits=c(-60,0))+
    scale_y_continuous(breaks= pretty_breaks(n=4))+
    expand_limits(y=0)+
    theme_pubr(legend = "top")
  g2 <- ggplot(ecorgn,aes(x=diff_days,y=spatial_hetero,color=as.factor(outbreak))) +
    geom_smooth() + ylab("Heterogeneity") + 
    xlab("Days (0 = date of outbreak)") +
    scale_color_discrete(name = "", labels = c("Non-outbreak", "Outbreak"))+
    geom_vline(xintercept=0,color="red",linetype=2)  + 
    geom_vline(xintercept=-35,color="black",linetype=2) + 
    geom_vline(xintercept=-21,color="black",linetype=2) + 
    scale_x_continuous(limits=c(-60,0))+
    scale_y_continuous(breaks= pretty_breaks(n=4))+
    facet_wrap(vars(Season),ncol=4,drop=FALSE)+
    expand_limits(y=0)+
    theme_pubr()  + theme(legend.position = "none")
  g3 <- ggplot(data = world) + geom_sf() + 
    geom_sf(data = (filtered_ecoregions2 %>% filter(ECO_NAME == toString(i))), aes(fill = ECO_NAME)) +
    coord_sf(xlim = c(-20,80), ylim = c(0,40), expand = FALSE) +
    theme_pubr() + theme(legend.position = "none")
  g4 <- ggplot((date_dat_projected2 %>% filter(ECO_NAME == toString(i),outbreak==0)),
               aes(x=DATE,y=as.integer(outbreak),color=(outbreak))) + geom_col(color="black")+
    theme_pubr() + ylab("Number of observations") + xlab("Year") +
    scale_color_discrete(name = "", labels = c("Outbreak")) +
    ggtitle(paste(no_outbreaks))+
    scale_y_continuous(breaks= pretty_breaks(n=4))+
    scale_x_date(date_breaks = "5 years", date_minor_breaks = "1 year",
                 date_labels = "%Y",limits= as.Date(c('2000-01-01','2020-12-31')))
  
  g5 <- ggplot((date_dat_projected2 %>% filter(ECO_NAME == toString(i),outbreak==1)),
               aes(x=DATE,y=as.integer(outbreak),color=(outbreak))) + geom_col(color="black")+
    theme_pubr() + ylab("Number of observations") + xlab("Year") +
    scale_color_discrete(name = "", labels = c("Outbreak")) +
    ggtitle(paste(outbreaks))+
    scale_y_continuous(breaks= pretty_breaks(n=4))+
    scale_x_date(date_breaks = "5 years", date_minor_breaks = "1 year",
                 date_labels = "%Y",limits= as.Date(c('2000-01-01','2020-12-31')))
  
  t <- toString(i)
  p1 <- ggplotGrob(g1)
  p2 <- ggplotGrob(g2)
  p3 <- ggplotGrob(g3)
  p4 <- ggplotGrob(g4)
  p5 <- ggplotGrob(g5)
  
  lay <- rbind(c(3,3,3,3,3),
               c(3,3,3,3,3),
               c(1,1,1,4,4),
               c(1,1,1,5,5),
               c(2,2,2,2,2))
  
  eco[[i]] <- grid.arrange(
    grobs = list(p1,p2,p3,p4,p5),
    layout_matrix = lay,
    top=t
  )
}


ml <- marrangeGrob(eco, nrow=1, ncol=1)
ggsave("Figures/DL_breakdown_ecoregion_March22.pdf", ml,height=11,width=8.5,units="in")

# mod results by bioregion

mod <- readRDS("models/DL_mods/mod.rds")
dat$preds <- predict(mod,type="response")


eco <- list()
for(i in eco_names){
  print(paste("Now computing",i))
  ecorgn <- dat %>% filter(ECO_NAME == i)
  
  ob_summary <- date_dat_projected2 %>% 
    dplyr::filter(ECO_NAME == i) %>% 
    group_by(outbreak) %>% 
    tally()
  
  no_outbreaks <- paste0("Non-outbreak (N: ",(ob_summary %>% filter(outbreak=="0"))$n,")")
  outbreaks <- paste0("Outbreak (N: ",(ob_summary %>% filter(outbreak=="1"))$n,")")

  g1 <- ggplot(ecorgn,aes(x=diff_days,y=spatial_hetero,z=preds)) +
    stat_summary_hex(bins=21) + ylab("Heterogeneity") + xlab("") +
    xlab("Days (0 = date of outbreak)") +
    scale_fill_viridis() +
    scale_x_continuous(limits=c(-60,0))+
    #scale_y_continuous(breaks= pretty_breaks(n=4))+
    #expand_limits(y=0)+
    theme_pubr(legend = "top")
  g2 <- ggplot(ecorgn,aes(x=diff_days,y=spatial_hetero,z=preds)) +
    stat_summary_hex(bins=21) + ylab("Heterogeneity") + 
    xlab("Days (0 = date of outbreak)") +
    scale_fill_viridis() +
    scale_x_continuous(limits=c(-60,0))+
    #scale_y_continuous(breaks= pretty_breaks(n=4))+
    facet_wrap(vars(Season),ncol=4,drop=FALSE)+
   # expand_limits(y=0)+
    theme_pubr()  + theme(legend.position = "none")
  g3 <- ggplot(data = world) + geom_sf() + 
    geom_sf(data = (filtered_ecoregions2 %>% filter(ECO_NAME == toString(i))), aes(fill = ECO_NAME)) +
    coord_sf(xlim = c(-20,80), ylim = c(0,40), expand = FALSE) +
    theme_pubr() + theme(legend.position = "none")
  g4 <- ggplot((date_dat_projected2 %>% filter(ECO_NAME == toString(i),outbreak==0)),
               aes(x=DATE,y=as.integer(outbreak),color=(outbreak))) + geom_col(color="black")+
    theme_pubr() + ylab("Number of observations") + xlab("Year") +
    scale_color_discrete(name = "", labels = c("Outbreak")) +
    ggtitle(paste(no_outbreaks))+
    scale_y_continuous(breaks= pretty_breaks(n=4))+
    scale_x_date(date_breaks = "5 years", date_minor_breaks = "1 year",
                 date_labels = "%Y",limits= as.Date(c('2000-01-01','2020-12-31')))
  
  g5 <- ggplot((date_dat_projected2 %>% filter(ECO_NAME == toString(i),outbreak==1)),
               aes(x=DATE,y=as.integer(outbreak),color=(outbreak))) + geom_col(color="black")+
    theme_pubr() + ylab("Number of observations") + xlab("Year") +
    scale_color_discrete(name = "", labels = c("Outbreak")) +
    ggtitle(paste(outbreaks))+
    scale_y_continuous(breaks= pretty_breaks(n=4))+
    scale_x_date(date_breaks = "5 years", date_minor_breaks = "1 year",
                 date_labels = "%Y",limits= as.Date(c('2000-01-01','2020-12-31')))
  
  t <- toString(i)
  p1 <- ggplotGrob(g1)
  p2 <- ggplotGrob(g2)
  p3 <- ggplotGrob(g3)
  p4 <- ggplotGrob(g4)
  p5 <- ggplotGrob(g5)
  
  lay <- rbind(c(3,3,3,3,3),
               c(3,3,3,3,3),
               c(1,1,1,4,4),
               c(1,1,1,5,5),
               c(2,2,2,2,2))
  
  eco[[i]] <- grid.arrange(
    grobs = list(p1,p2,p3,p4,p5),
    layout_matrix = lay,
    top=t
  )
}


ml <- marrangeGrob(eco, nrow=1, ncol=1)
ggsave("Figures/DL_breakdown_ecoregion_March22_modeled.pdf", ml,height=11,width=8.5,units="in")

k.check(mod)
