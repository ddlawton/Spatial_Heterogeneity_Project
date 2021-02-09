#####
# APLC Database Work
#####

rm(list=ls())
library(tidyverse)
library(anytime)
library(splitstackshape)
library(data.table)
library(mgcv)
library(gratia)
library("rnaturalearth")
library("rnaturalearthdata")
library(ggpubr)
library(viridis)
library(gridExtra)
dat <- read.csv("data/processed/APLC_database_hierarchical_data_Feb_8.csv")

tab <- dat %>% group_by(as.factor(Species), as.factor(NymphDensity)) %>% tally() %>% pivot_wider(names_from = `as.factor(NymphDensity)`, values_from=n)
view(tab)

str(dat)

cols <- c("AdultDensity", "AdultStage", "DataQuality", "NymphDensity","NymphStage","Source","Species",
          "spp_code","binary_outbreak","REG_NAME_7","Minor_rain_zones","Major_rain_zones","Source_code",
          "Season")
dat[cols] <- lapply(dat[cols], factor)



CT <- dat %>%
  filter(
    Species == '10' | Species == '11'
  )

AG <- dat %>%
  filter(
    Species == '10' | Species == '12'
  )

AC <- dat %>%
  filter(
    Species == '10' | Species == '15'
  )

OA <- dat %>%
  filter(
    Species == '10' | Species == '19'
  )

CT <- CT %>% filter(between(diff_days,-60,0)) %>% filter(NDVIs >= 0)
names(CT)

ggplot(CT,aes(x=diff_days,y=NDVIs/1000,color=binary_outbreak)) + geom_smooth() + theme_pubr()
ggplot(CT,aes(x=diff_days,y=dissimilarity,color=binary_outbreak)) + geom_smooth() + theme_pubr()



dat <- CT

date_dat <- APL_points

date_dat$Date <- as.Date(date_dat$Date, format="%Y-%m-%d")
names(date_dat)


bioregions <- readOGR(dsn = "Australian_plague_locust_model/Data/Processed",
                      layer = "Simplified_bioregion")
filtered_bioregions3 <- st_as_sf(bioregions)

world <- ne_countries(scale = "medium", returnclass = "sf")

bio_names <- unique(dat$REG_NAME_7)

date_dat$Outbreak <- as.factor(date_dat$binary_outbreak)

dat$Outbreak <- as.factor(dat$binary_outbreak)
dat$Bioregn <- (dat$REG_NAME_7)
date_dat$Bioregn <- (date_dat$REG_NAME_7)

dat <- dat %>% filter(Season != "NA")

bio <- list()
for(i in bio_names){
  print(paste("Now computing",i))
  
  biorgn <- dat %>% filter(Bioregn == i)
  
  ob_summary <- date_dat %>% 
    dplyr::filter(Bioregn == i) %>% 
    group_by(Outbreak) %>% 
    tally()
  
  no_outbreaks <- paste0("Non-outbreak (N: ",(ob_summary %>% filter(Outbreak=="0"))$n,")")
  outbreaks <- paste0("Outbreak (N: ",(ob_summary %>% filter(Outbreak=="1"))$n,")")
  
  
  g1 <- ggplot(biorgn,aes(x=days_diff,y=mean,color=(Outbreak))) +
    geom_smooth() + ylab("NDVI") + xlab("") +
    xlab("Days (0 = date of outbreak)")  +
    scale_color_manual(name = "", labels = c("Non-Outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62"))+
    geom_vline(xintercept=0,color="red",linetype=2)  + 
    geom_vline(xintercept=-35,color="black",linetype=2) + 
    geom_vline(xintercept=-21,color="black",linetype=2) + 
    scale_x_continuous(limits=c(-78,32),breaks=c(-78,-32,0,32))+
    scale_y_continuous(breaks= pretty_breaks(n=4))+
    expand_limits(y=0)+
    theme_pubr(legend = "top")
  g2 <- ggplot(biorgn,aes(x=days_diff,y=mean,color=(Outbreak))) +
    geom_smooth() + ylab("NDVI") + 
    xlab("Days (0 = date of outbreak)") +
    scale_color_manual(name = "", labels = c("Non-Outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62"))+
    geom_vline(xintercept=0,color="red",linetype=2)  + 
    geom_vline(xintercept=-35,color="black",linetype=2) + 
    geom_vline(xintercept=-21,color="black",linetype=2) + 
    scale_x_continuous(limits=c(-78,32),breaks=c(-78,-32,0,32))+
    scale_y_continuous(breaks= pretty_breaks(n=4))+
    facet_wrap(vars(Season),ncol=3,drop=FALSE)+
    expand_limits(y=0)+
    theme_pubr()  + theme(legend.position = "none")
  g3 <- ggplot(data = world) + geom_sf() + 
    geom_sf(data = (filtered_bioregions3 %>% filter(REG_NAME_7 == i)), aes(fill = REG_NAME_7)) +
    coord_sf(xlim = c(130,155), ylim = c(-40,-15), expand = FALSE)+
    theme_pubr() + theme(legend.position = "none")
  g4 <- ggplot((date_dat %>% filter(Bioregn == toString(i),Outbreak==0)),
               aes(x=Date,y=(as.integer(Outbreak)+1),color=(Outbreak))) + geom_col(color="black")+
    theme_pubr() + ylab("Number of observations") + xlab("Year") +
    scale_color_discrete(name = "", labels = c("Outbreak")) +
    ggtitle(paste(no_outbreaks))+
    scale_y_continuous(breaks= pretty_breaks(n=4))+
    scale_x_date(date_breaks = "5 years", date_minor_breaks = "1 year",
                 date_labels = "%Y",limits= as.Date(c('2000-01-01','2020-12-31')))
  
  g5 <- ggplot((date_dat %>% filter(Bioregn == toString(i),Outbreak==1)),
               aes(x=Date,y=as.integer(Outbreak),color=(Outbreak))) + geom_col(color="black")+
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
  
  bio[[i]] <- grid.arrange(
    grobs = list(p1,p2,p3,p4,p5),
    layout_matrix = lay,
    top=t
  )
}



plot <- grid.arrange(
  grobs = list(p1,p2,p3,p4,p5),
  layout_matrix = lay,
  top=t)

ml <- marrangeGrob(bio, nrow=1, ncol=1)
ggsave("MS_figures/APL_bioregion_breakdown_Dec25.pdf", ml,height=11,width=8.5,units="in")







##########
# Garaged Code
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

P <- ggplot(data = world) + geom_sf() + geom_point(data=dat,
                                                   aes(x=Longitude,y=Latitude),size=.1) +
  coord_sf(xlim = c(-20,80), ylim = c(0,40), expand = FALSE) +
  scale_color_manual(name = "", labels = c("Non-outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62"))+
  theme_pubr() + theme(legend.position = "none")  

P



CT <- dat8 %>%
  filter(
    Species == '10' | Species == '11'
  )

AG <- dat8 %>%
  filter(
    Species == '10' | Species == '12'
  )

AC <- dat8 %>%
  filter(
    Species == '10' | Species == '15'
  )

OA <- dat8 %>%
  filter(
    Species == '10' | Species == '19'
  )

 # CT
g1 <- ggplot(CT,aes(x=diff_days,y=NDVIs/1000,color=as.factor(outbreak))) + geom_smooth() + ylab("NDVIs") + xlab("Days (0 = observation date)") +
  scale_color_manual(name = "", labels = c("Non-outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62")) + theme_pubr() + xlim(-70,0) + scale_y_continuous(expand = c(0, 0))

g2 <- ggplot(CT,aes(x=diff_days,y=dissimilarity,color=as.factor(outbreak))) + geom_smooth() + ylab("Dissimilarity") + xlab("Days (0 = observation date)") +
  scale_color_manual(name = "", labels = c("Non-outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62")) + theme_pubr() + xlim(-70,0) + scale_y_continuous(expand = c(0, 0))

g3 <- ggplot(CT,aes(x=diff_days,y=contrasts,color=as.factor(outbreak))) + geom_smooth() + ylab("contrasts") + xlab("Days (0 = observation date)") +
  scale_color_manual(name = "", labels = c("Non-outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62")) + theme_pubr() + xlim(-70,0) + scale_y_continuous(expand = c(0, 0))
 
g4 <-ggplot(CT,aes(x=diff_days,y=entropy,color=as.factor(outbreak))) + geom_smooth() + ylab("entropy") + xlab("Days (0 = observation date)") +
  scale_color_manual(name = "", labels = c("Non-outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62")) + theme_pubr() + xlim(-70,0) + scale_y_continuous(expand = c(0, 0))

g5 <-ggplot(CT,aes(x=diff_days,y=correlation,color=as.factor(outbreak))) + geom_smooth() + ylab("correlation") + xlab("Days (0 = observation date)") +
  scale_color_manual(name = "", labels = c("Non-outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62")) + theme_pubr() + xlim(-70,0) + scale_y_continuous(expand = c(0, 0))
 
p1 <- ggplotGrob(g1)
p2 <- ggplotGrob(g2)
p3 <- ggplotGrob(g3)
p4 <- ggplotGrob(g4)
p5 <- ggplotGrob(g5)

lay <- rbind(c(1,2),
             c(3,4),
             c(5,NA))

grid.arrange(
  grobs = list(p1,p2,p3,p4,p5),
  layout_matrix = lay)

# AG
g6 <- ggplot(AG,aes(x=diff_days,y=NDVIs/1000,color=as.factor(NymphDensity))) + geom_smooth() + ylab("NDVIs") + xlab("Days (0 = observation date)") +
  theme_pubr() + xlim(-70,0) + scale_y_continuous(expand = c(0, 0))

g7 <- ggplot(AG,aes(x=diff_days,y=dissimilarity,color=as.factor(NymphDensity))) + geom_smooth() + ylab("Dissimilarity") + xlab("Days (0 = observation date)") +
  theme_pubr() + xlim(-70,0) + scale_y_continuous(expand = c(0, 0))

g8 <- ggplot(AG,aes(x=diff_days,y=contrasts,color=as.factor(NymphDensity))) + geom_smooth() + ylab("contrasts") + xlab("Days (0 = observation date)") +
  theme_pubr() + xlim(-70,0) + scale_y_continuous(expand = c(0, 0))

g9 <-ggplot(AG,aes(x=diff_days,y=entropy,color=as.factor(NymphDensity))) + geom_smooth() + ylab("entropy") + xlab("Days (0 = observation date)") +
  theme_pubr() + xlim(-70,0) + scale_y_continuous(expand = c(0, 0))

g10 <-ggplot(AG,aes(x=diff_days,y=correlation,color=as.factor(NymphDensity))) + geom_smooth() + ylab("correlation") + xlab("Days (0 = observation date)") +
  theme_pubr() + xlim(-70,0) + scale_y_continuous(expand = c(0, 0))

p6 <- ggplotGrob(g6)
p7 <- ggplotGrob(g7)
p8 <- ggplotGrob(g8)
p9 <- ggplotGrob(g9)
p10 <- ggplotGrob(g10)

lay <- rbind(c(1,2),
             c(3,4),
             c(5,NA))

grid.arrange(
  grobs = list(p6,p7,p8,p9,p10),
  layout_matrix = lay)


# AC
g11 <- ggplot(AC,aes(x=diff_days,y=NDVIs/1000,color=as.factor(NymphDensity))) + geom_smooth() + ylab("NDVIs") + xlab("Days (0 = observation date)") +
  theme_pubr() + xlim(-70,0) + scale_y_continuous(expand = c(0, 0))

g12 <- ggplot(AC,aes(x=diff_days,y=dissimilarity,color=as.factor(NymphDensity))) + geom_smooth() + ylab("Dissimilarity") + xlab("Days (0 = observation date)") +
  theme_pubr() + xlim(-70,0) + scale_y_continuous(expand = c(0, 0))

g13 <- ggplot(AC,aes(x=diff_days,y=contrasts,color=as.factor(NymphDensity))) + geom_smooth() + ylab("contrasts") + xlab("Days (0 = observation date)") +
  theme_pubr() + xlim(-70,0) + scale_y_continuous(expand = c(0, 0))

g14 <-ggplot(AC,aes(x=diff_days,y=entropy,color=as.factor(NymphDensity))) + geom_smooth() + ylab("entropy") + xlab("Days (0 = observation date)") +
  theme_pubr() + xlim(-70,0) + scale_y_continuous(expand = c(0, 0))

g15 <-ggplot(AC,aes(x=diff_days,y=correlation,color=as.factor(NymphDensity))) + geom_smooth() + ylab("correlation") + xlab("Days (0 = observation date)") +
  theme_pubr() + xlim(-70,0) + scale_y_continuous(expand = c(0, 0))

p11 <- ggplotGrob(g11)
p12 <- ggplotGrob(g12)
p13 <- ggplotGrob(g13)
p14 <- ggplotGrob(g14)
p15 <- ggplotGrob(g15)

lay <- rbind(c(1,2),
             c(3,4),
             c(5,NA))

grid.arrange(
  grobs = list(p11,p12,p13,p14,p15),
  layout_matrix = lay)

# OA
g16 <- ggplot(OA,aes(x=diff_days,y=NDVIs/1000,color=as.factor(NymphDensity))) + geom_smooth() + ylab("NDVIs") + xlab("Days (0 = observation date)") +
  theme_pubr() + xlim(-70,0) + scale_y_continuous(expand = c(0, 0))

g17 <- ggplot(OA,aes(x=diff_days,y=dissimilarity,color=as.factor(NymphDensity))) + geom_smooth() + ylab("Dissimilarity") + xlab("Days (0 = observation date)") +
  theme_pubr() + xlim(-70,0) + scale_y_continuous(expand = c(0, 0))

g18 <- ggplot(OA,aes(x=diff_days,y=contrasts,color=as.factor(NymphDensity))) + geom_smooth() + ylab("contrasts") + xlab("Days (0 = observation date)") +
  theme_pubr() + xlim(-70,0) + scale_y_continuous(expand = c(0, 0))

g19 <-ggplot(OA,aes(x=diff_days,y=entropy,color=as.factor(NymphDensity))) + geom_smooth() + ylab("entropy") + xlab("Days (0 = observation date)") +
  theme_pubr() + xlim(-70,0) + scale_y_continuous(expand = c(0, 0))

g20 <-ggplot(OA,aes(x=diff_days,y=correlation,color=as.factor(NymphDensity))) + geom_smooth() + ylab("correlation") + xlab("Days (0 = observation date)") +
  theme_pubr() + xlim(-70,0) + scale_y_continuous(expand = c(0, 0))

p16 <- ggplotGrob(g16)
p17 <- ggplotGrob(g17)
p18 <- ggplotGrob(g18)
p19 <- ggplotGrob(g19)
p20 <- ggplotGrob(g20)

lay <- rbind(c(1,2),
             c(3,4),
             c(5,NA))

grid.arrange(
  grobs = list(p16,p17,p18,p19,p20),
  layout_matrix = lay)





CT2 <- CT %>% filter(diff_days >= -60)

mod <- bam(outbreak~ te(diff_days,((NDVIs))) + te(diff_days,((entropy))) +
             te(diff_days,((dissimilarity))) + te(diff_days,((contrasts))) + 
             te(Longitude,Latitude), family=binomial(),select = TRUE,discrete=TRUE,nthreads=4,data=CT2)



summary(mod)
draw(mod)
appraise(mod)
k.check(mod)

CT$preds <- predict(mod,newdata = CT,type="response")



ggplot(CT,aes(x=diff_days,y=NDVIs/1000,z=(preds))) +
  stat_summary_2d() + ylab("NDVI")  +
  scale_fill_viridis() + ylim(0,.25) 

ggplot(CT,aes(x=diff_days,y=dissimilarity,z=(preds))) +
  stat_summary_2d() + ylab("")  +
  scale_fill_viridis()+ ylim(0,75)

ggplot(CT,aes(x=diff_days,y=contrasts,z=(preds))) +
  stat_summary_2d() + ylab("")  +
  scale_fill_viridis() + ylim(0,15000)

ggplot(CT,aes(x=diff_days,y=entropy,z=(preds))) +
  stat_summary_2d() + ylab("")  +
  scale_fill_viridis() 

