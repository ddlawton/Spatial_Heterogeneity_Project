#####
# APLC Database Work
#####

rm(list=ls())
library(parallel)
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
library(lubridate)
dat <- fread("data/processed/APLC_database_hierarchical_data_Feb_8.csv")

tab <- dat %>% group_by(as.factor(Species), as.factor(NymphDensity)) %>% tally() %>% pivot_wider(names_from = `as.factor(NymphDensity)`, values_from=n)
view(tab)

tab <- dat %>% group_by(REG_NAME_7,binary_outbreak) %>% tally() %>% pivot_wider(names_from = binary_outbreak,values_from = n) %>%
  mutate(total = `0`+`1`) %>%  filter(total > 200) %>% filter(`1` >= 100)
view(tab)

tab_names <- unique(as.factor(tab$REG_NAME_7)) %>% droplevels()

dat <- dat %>% filter((REG_NAME_7 %in% tab_names)) %>% droplevels()

cols <- c("AdultDensity", "AdultStage", "DataQuality", "NymphDensity","NymphStage","Source","Species",
          "spp_code","binary_outbreak","REG_NAME_7","Minor_rain_zones","Major_rain_zones","Source_code",
          "Season")

for(col in cols){
  set(dat, j = col, value = as.factor(dat[[col]]))
}

dat$Date <- as.Date(dat$Date)
dat$Year <- year(dat$Date)
dat$DOY <- (yday(dat$Date) + 183) %% 366



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


CT <- CT %>%
  filter(between(diff_days,-60,0)) %>% filter(NDVIs >= 0) %>%
  mutate(
    contrasts_scaled = scale(contrasts,center = TRUE, scale=TRUE),
    correlation_scaled = scale(correlation,center = TRUE, scale=TRUE),
    dissimilarity_scaled = scale(dissimilarity,center = TRUE, scale=TRUE),
    entropy_scaled = scale(entropy,center = TRUE, scale=TRUE),
    homogeneity_scaled = scale(homogeneity,center = TRUE, scale=TRUE)) 

CT <- CT %>% mutate(
    contrasts_normalized = (CT$contrasts_scaled - min(CT$contrasts_scaled))/(max(CT$contrasts_scaled)-min(CT$contrasts_scaled)),
    correlation_normalized = (CT$correlation_scaled - min(CT$correlation_scaled))/(max(CT$correlation_scaled)-min(CT$correlation_scaled)),
    dissimilarity_normalized = (CT$dissimilarity_scaled - min(CT$dissimilarity_scaled))/(max(CT$dissimilarity_scaled)-min(CT$dissimilarity_scaled)),
    entropy_normalized = (CT$entropy_scaled - min(CT$entropy_scaled))/(max(CT$entropy_scaled)-min(CT$entropy_scaled)),
    homogeneity_normalized = (CT$homogeneity_scaled - min(CT$homogeneity_scaled))/(max(CT$homogeneity_scaled)-min(CT$homogeneity_scaled))) 

CT <- CT %>% mutate(
    spatial_hetero = rowMeans(.[,35:39], na.rm = TRUE)
  )




summary(CT$contrasts)
summary(CT$correlation)
summary(CT$dissimilarity)
summary(CT$entropy)
summary(CT$homogeneity)

#ggplot(CT,aes(x=spatial_hetero,y=homogeneity_scaled)) + geom_smooth() + theme_pubr()



tab <- CT %>% group_by(REG_NAME_7, binary_outbreak)  %>% tally() %>% pivot_wider(names_from = binary_outbreak,values_from = n)

tab2 <- (tab %>% filter(`0` > 200, `1` > 10) %>% droplevels()) 
names_list <- levels(tab2$REG_NAME_7)


summary(CT)


#ggplot(CT,aes(x=spatial_hetero,y=NDVIs/1000)) + geom_smooth() + theme_pubr()
#ggplot(CT,aes(x=diff_days,y=spatial_hetero,color=as.factor(binary_outbreak))) + geom_smooth() + theme_pubr()

names(CT)

CT %>%
  group_by(as.factor(Year),as.factor(binary_outbreak)) %>%
  tally() %>% pivot_wider(names_from = `as.factor(binary_outbreak)`,values_from = n)

mod <- bam(as.factor(binary_outbreak) ~ 
             te(diff_days, spatial_hetero) + 
             te(Longitude, Latitude, k = 15) + 
             s((REG_NAME_7), bs = "re") + 
             s((Major_rain_zones), bs = "re") + 
             s((Season), bs = "re") +
             s(DOY,bs="cc", k=12) +
             s(Year,bs="tp", k=15) +
             ti(DOY,Year, bs = c("cc", "tp"), k = c(12, 15)),
           family=binomial(),data=CT,select = TRUE,discrete=TRUE,nthreads=8)

mod2 <- bam(as.factor(binary_outbreak) ~ 
             te(diff_days, spatial_hetero) + 
             te(Longitude, Latitude, k = 15) + 
             s((REG_NAME_7), bs = "re") + 
             s((Major_rain_zones), bs = "re") + 
             s((Season), bs = "re") +
             te(DOY,Year, bs = c("cc", "tp"), k = c(12, 15)),
           family=binomial(),data=CT,select = TRUE,discrete=TRUE,nthreads=8)


mod3 <- bam(as.factor(binary_outbreak) ~ 
              te(diff_days, spatial_hetero) + 
              s((REG_NAME_7), bs = "re") + 
              s((Major_rain_zones), bs = "re") + 
              s((Season), bs = "re") +
              te(DOY,Year, bs = c("cc", "tp"), k = c(12, 15))+
              te(Year,Longitude, Latitude, k = 15), 
            family=binomial(),data=CT,select = TRUE,chunk.size=5000,cluster=cl)

###see if you have multiple cores

detectCores()

###indicate number of cores used for parallel processing
if (detectCores()>1) {
  cl <- makeCluster(detectCores()-1)
} else cl <- NULL

cl

system.time(mod1 <- bam(as.factor(binary_outbreak) ~ 
                          te(diff_days, spatial_hetero) + 
                          te(Longitude, Latitude, k = 15) + 
                          s((REG_NAME_7), bs = "re") + 
                          s((Major_rain_zones), bs = "re") + 
                          s((Season), bs = "re") +
                          s(DOY,bs="cc", k=12) +
                          s(Year,bs="tp", k=15) +
                          ti(DOY,Year, bs = c("cc", "tp"), k = c(12, 15)),
                        family=binomial(),data=CT,select = TRUE,cluster=cl))

if (!is.null(cl)) stopCluster(cl)




summary(mod)
summary(mod2)
AIC(mod,mod2)
draw(mod2)
draw(mod)
concurvity(mod)
k.check(mod)

names(CT)
pdat <- with(CT, expand.grid(Major_rain_zones = levels(Major_rain_zones), 
                             REG_NAME_7 = levels(REG_NAME_7), Season = levels(Season), 
                             Longitude = seq(min(Longitude), max(Longitude), length = 10),
                             Latitude = seq(min(Latitude), max(Latitude), length = 10),
                             diff_days = seq(min(diff_days), max(diff_days), length = 10),
                             spatial_hetero = seq(min(spatial_hetero), max(spatial_hetero), length = 10)
)
)

pdat <- transform(pdat, pred = predict(mod, newdata = pdat, type = "response"))
names(pdat)
ggplot(pdat,aes(x=diff_days,y=spatial_hetero,z=(pred))) + stat_summary_2d(bins=15) + scale_fill_viridis()
draw(mod)




#draw(mod)
#appraise(mod)
#names(CT2)
#CT2$preds <- predict(mod,newdata = CT2,type="response")
#CT2$Season
#ggplot(CT2,aes(x=diff_days,y=NDVIs,z=preds)) +stat_summary_hex() + scale_fill_viridis() #+ ylim(0.2,0.4)

#write.csv(CT2,file="Data/Processed/CT_modeling_dat.csv")

CT3 <- CT2 %>% filter(REG_NAME_7 %in% names_list)


modI <- bam(
  binary_outbreak ~ 
    te(diff_days, spatial_hetero, by=Major_rain_zones, bs=c("tp", "tp"),k=c(5, 5), m=2)+
    te(diff_days, spatial_hetero, by=REG_NAME_7, bs=c("tp", "tp"),k=c(5, 5), m=2)+
    te(diff_days, spatial_hetero, by=Season, bs=c("tp", "tp"),k=c(5, 5), m=2)+
    te(diff_days, NDVIs, by=Major_rain_zones, bs=c("tp", "tp"),k=c(5, 5), m=2)+
    te(diff_days, NDVIs, by=REG_NAME_7, bs=c("tp", "tp"),k=c(5, 5), m=2)+
    te(diff_days, NDVIs, by=Season, bs=c("tp", "tp"),k=c(5, 5), m=2)+
    s(Latitude,Longitude,k=50),
  family=binomial(),select = TRUE,discrete=TRUE,nthreads=8,data=CT3,method="fREML",
  drop.unused.levels=FALSE)



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

