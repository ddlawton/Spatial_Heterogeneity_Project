######
# Australian Plague Locust --
#  Exporatory graphing
#   Last modified: March 8th 2020
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


rawdat <- fread("Data/raw/CT_dat.csv")
names(rawdat)

dens_target <- c("0","1","2","3")

rawdat$Date <- as.Date(rawdat$Date,format="%d-%b-%y")
rawdat$Source <- as.factor(rawdat$Source)
rawdat$Species <- as.factor(rawdat$Species)
rawdat$spp_code <- as.factor(rawdat$spp_code)
rawdat$`Adult Density` <- as.factor(rawdat$`Adult Density`)
rawdat$`Adult Stage`<- as.factor(rawdat$`Adult Stage`)
rawdat$`Data Quality`<- as.factor(rawdat$`Data Quality`)
rawdat$`Nymph Density`<- as.factor(rawdat$`Nymph Density`)
rawdat$`Nymph Stage`<- as.factor(rawdat$`Nymph Stage`)


rawdat <- rawdat %>%  mutate(outbreak = case_when(
  `Nymph Density` %in% dens_target ~ "0",
  `Nymph Density` == "4" ~ "1"
)
)
names(rawdat)

rawdat <- rawdat %>% dplyr::select("Date","Latitude","Longitude","outbreak")


dat <- fread("Data/processed/APLC_database_hierarchical_data_Feb_8.csv")

dat$Date <- as.Date(dat$Date,format="%m/%d/%y")
dat$Source <- as.factor(dat$Source)
dat$Species <- as.factor(dat$Species)
dat$AdultDensity <- as.factor(dat$AdultDensity)
dat$AdultStage<- as.factor(dat$AdultStage)
dat$DataQuality<- as.factor(dat$DataQuality)
dat$NymphDensity<- as.factor(dat$NymphDensity)
dat$NymphStage<- as.factor(dat$NymphStage)
dat$spp_code<- as.factor(dat$spp_code)
dat$Season<- as.factor(dat$Season)
dat$Major_rain_zones<- as.factor(dat$Major_rain_zones)
dat$Minor_rain_zones<- as.factor(dat$Minor_rain_zones)
dat$REG_NAME_7<- as.factor(dat$REG_NAME_7)

dat <- dat %>% filter(spp_code == "None" | spp_code == "CT") %>%
  dplyr::select(6,8:13,15:29)

summary(dat)



dat$NDVIs <- dat$NDVIs/1000

dat <- dat %>% drop_na(c(contrasts,correlation,dissimilarity,entropy,homogeneity,NDVIs)) %>% filter(between(NDVIs,0,1)) %>% filter(between(diff_days,-60,0))

dat$contrasts_normalized <- (dat$contrasts - min(dat$contrasts))/(max(dat$contrasts)-min(dat$contrasts))
#dat$correlation_normalized <- (dat$correlation - min(dat$correlation))/(max(dat$correlation)-min(dat$correlation))
dat$dissimilarity_normalized <- (dat$dissimilarity - min(dat$dissimilarity))/(max(dat$dissimilarity)-min(dat$dissimilarity))
dat$entropy_normalized <- (dat$entropy - min(dat$entropy))/(max(dat$entropy)-min(dat$entropy))
#dat$homogeneity_normalized <- (dat$homogeneity - min(dat$homogeneity))/(max(dat$homogeneity)-min(dat$homogeneity))


dat <- dat %>% mutate(
  spatial_hetero = rowMeans(dplyr::select(.,ends_with("_normalized")), na.rm = TRUE)) 

summary(dat$Season)

flevels  <- c("Winter","Spring","Summer","Fall")
dat$Season <- factor(dat$Season,levels=flevels)


ggplot(dat,aes(x=as.factor(binary_outbreak),y=spatial_hetero)) + geom_boxplot()

ggplot(dat,aes(x=NDVIs,y=spatial_hetero)) + geom_smooth()

ggplot(dat,aes(x=diff_days,y=NDVIs,color=as.factor(binary_outbreak))) +geom_smooth() + xlim(-60,0)

ggplot(dat,aes(x=diff_days,y=dissimilarity,color=as.factor(binary_outbreak))) +geom_smooth()  + 
  ylab("dissimilarity")

ggplot(dat,aes(x=diff_days,y=contrasts,color=as.factor(binary_outbreak))) +geom_smooth()  + 
  ylab("contrasts")

ggplot(dat,aes(x=dissimilarity,y=contrasts)) + geom_smooth()

ggplot(dat,aes(x=diff_days,y=entropy,color=as.factor(binary_outbreak))) +geom_smooth()  + # this one is weird... 
  ylab("entropy")

ggplot(dat,aes(x=diff_days,y=spatial_hetero,color=as.factor(binary_outbreak))) +geom_smooth()  + 
  ylab("Spatial Heoterogenity")


# Bioregion breakdown

date_dat <- rawdat

ecoregions <- readOGR(dsn="/Users/ddlawton/Dropbox (ASU)/Research/Graduate/Chapter 4 - Hierarchical scaling of locust swarms/Hierarchical_locust_modeling/Australian_plague_locust_model/Data/Raw/IBRA_bioregions", layer = "IBRA_bioregions")


crs(ecoregions)

ecoregions <- spTransform(ecoregions,CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")) # something is going wrong with the CRS command...doesnt know what GDA94 is??


coords <- cbind(date_dat$Longitude,date_dat$Latitude)
date_dat_projected <- SpatialPointsDataFrame(coords,data=date_dat,
                        proj4string = crs("+proj=lcc +lat_0=-33.25 +lon_0=147 +lat_1=-30.75 +lat_2=-35.75 +x_0=9300000 +y_0=4500000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

?CRS
plot(date_dat_projected)
projection(ecoregions) == projection(date_dat_projected) #check to see if projections are the same


test <- point.in.poly(date_dat_projected,ecoregions) #make sure in the same projection!
test2 <- test@data
names(test2)

date_dat_projected2 <- test2 %>%
  dplyr::select(1:4,REG_NAME_7)




tab <- dat %>% group_by(REG_NAME_7,binary_outbreak) %>% tally() %>% pivot_wider(names_from = binary_outbreak,values_from = n) %>%
  mutate(total = `0`+`1`) %>%  filter(total > 200) %>% filter(`1` >= 100)
view(tab)

tab_names <- unique(tab$REG_NAME_7) %>% droplevels()

dat <- dat %>% filter((REG_NAME_7 %in% tab_names)) %>% droplevels()

dat$REG_NAME_7 <- as.factor(dat$REG_NAME_7)
# Filtering out unneeded ecoregions
used_eco_regions <- unique(dat$REG_NAME_7) %>% droplevels()

filter_ecoregions <- ecoregions %>%
  filter(REG_NAME_7 %in% used_eco_regions)

filtered_ecoregions2 <- st_as_sf(filter_ecoregions)
world <- ne_countries(scale = "medium", returnclass = "sf")
names(dat)

head(date_dat_projected2)
eco_names <- unique(dat$REG_NAME_7) %>% droplevels()
date_dat_projected2$outbreak <- as.factor(date_dat_projected2$binary_outbreak)

eco <- list()
for(i in eco_names){
  print(paste("Now computing",i))
  ecorgn <- dat %>% filter( REG_NAME_7 == i)
  
  #ob_summary <- date_dat_projected2 %>% 
  #  dplyr::filter( REG_NAME_7 == i) %>% 
  #  group_by(outbreak) %>% 
  #  tally()
  
  #no_outbreaks <- paste0("Non-outbreak (N: ",(ob_summary %>% filter(outbreak=="0"))$n,")")
  #outbreaks <- paste0("Outbreak (N: ",(ob_summary %>% filter(outbreak=="1"))$n,")")
  
  g1 <- ggplot(ecorgn,aes(x=diff_days,y=contrasts,color=as.factor(binary_outbreak))) +
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
  g2 <- ggplot(ecorgn,aes(x=diff_days,y=contrasts,color=as.factor(binary_outbreak))) +
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
  #g3 <- ggplot(data = world) + geom_sf() + 
  #  geom_sf(data = (filtered_ecoregions2 %>% filter( REG_NAME_7 == toString(i))), aes(fill =  REG_NAME_7)) +
  #  coord_sf(xlim = c(-20,80), ylim = c(0,40), expand = FALSE) +
  #  theme_pubr() + theme(legend.position = "none")
  #g4 <- ggplot((date_dat_projected2 %>% filter( REG_NAME_7 == toString(i),outbreak==0)),
  #             aes(x=DATE,y=as.integer(outbreak),color=(outbreak))) + geom_col(color="black")+
  #  theme_pubr() + ylab("Number of observations") + xlab("Year") +
  #  scale_color_discrete(name = "", labels = c("Outbreak")) +
  #  ggtitle(paste(no_outbreaks))+
  #  scale_y_continuous(breaks= pretty_breaks(n=4))+
  #  scale_x_date(date_breaks = "5 years", date_minor_breaks = "1 year",
  #               date_labels = "%Y",limits= as.Date(c('2000-01-01','2020-12-31')))
  
  #g5 <- ggplot((date_dat_projected2 %>% filter( REG_NAME_7 == toString(i),outbreak==1)),
  #             aes(x=DATE,y=as.integer(outbreak),color=(outbreak))) + geom_col(color="black")+
  #  theme_pubr() + ylab("Number of observations") + xlab("Year") +
  #  scale_color_discrete(name = "", labels = c("Outbreak")) +
  #  ggtitle(paste(outbreaks))+
  #  scale_y_continuous(breaks= pretty_breaks(n=4))+
  #  scale_x_date(date_breaks = "5 years", date_minor_breaks = "1 year",
  #               date_labels = "%Y",limits= as.Date(c('2000-01-01','2020-12-31')))
  
  t <- toString(i)
  p1 <- ggplotGrob(g1)
  p2 <- ggplotGrob(g2)
  #p3 <- ggplotGrob(g3)
  #p4 <- ggplotGrob(g4)
  #p5 <- ggplotGrob(g5)
  
  #lay <- rbind(c(3,3,3,3,3),
  #             c(3,3,3,3,3),
  #             c(1,1,1,4,4),
  #             c(1,1,1,5,5),
  #             c(2,2,2,2,2))
  
  lay <- rbind(c(1,1,1,1,1),
               c(1,1,1,1,1),
               c(2,2,2,2,2))
  
  #eco[[i]] <- grid.arrange(
  #  grobs = list(p1,p2,p3,p4,p5),
  #  layout_matrix = lay,
  #  top=t
  #)
  
  eco[[i]] <- grid.arrange(
    grobs = list(p1,p2),
    layout_matrix = lay,
    top=t
  )
}


ml <- marrangeGrob(eco, nrow=1, ncol=1)
ggsave("Figures/APL_breakdown_ecoregion_Feb26.pdf", ml,height=11,width=8.5,units="in")

