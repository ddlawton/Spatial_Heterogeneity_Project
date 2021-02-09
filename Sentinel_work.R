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
dat <- fread("data/processed/Locust_outbreaks_jan_28_S2.csv")

dat$NDVIs <- as.factor(dat$NDVIs)
names(dat2)
dat2 <- dat %>%
  select(2:18) %>%
  filter(NDVIs != "[]")


world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

P <- ggplot(data = world) + geom_sf() + geom_point(data=dat,
                                                   aes(x=Longitude,y=Latitude),size=.1) +
  coord_sf(xlim = c(-20,80), ylim = c(0,40), expand = FALSE) +
  scale_color_manual(name = "", labels = c("Non-outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62"))+
  theme_pubr() + theme(legend.position = "none")  

P

dat4 <- cSplit(setDT(dat2)[, lapply(.SD, gsub, pattern = "[][}]", 
                          replacement = "")], names(dat2), sep=",", direction='long', fixed = FALSE, "long")


dat4[dat4 == ""] <- NA # define NA pattern
dat5 <- dat4[rowSums(is.na(dat4)) != ncol(dat4), ]


dat6 <- dat5 %>% tidyr::fill(c("DATE", "ID","LOCUSTID", "Latitude", 
                             "Longitude","SHPAPPGREG","SHPAPPSOL","SHPAPPTRAN",
                             "SHPAPPUNK","SHPDENGRP", "SHPDENISOL","unixdates"), .direction = 'down')


dat6$SHPAPPGREG <- as.factor(dat6$SHPAPPGREG)
dat6$SHPAPPSOL <- as.factor(dat6$SHPAPPSOL)
dat6$SHPAPPTRAN <- as.factor(dat6$SHPAPPTRAN)
dat6$SHPAPPUNK <- as.factor(dat6$SHPAPPUNK)
dat6$SHPDENGRP <- as.factor(dat6$SHPDENGRP)
dat6$SHPDENISOL <- as.factor(dat6$SHPDENISOL)


dat6 %>%
  group_by(SHPAPPGREG,SHPAPPSOL,SHPAPPTRAN,SHPAPPUNK,SHPDENGRP,SHPDENISOL) %>%
  tally()

dat7 <- dat6 %>%
  filter(NDVIs >= 0) 


dat7$diff_days <- as.numeric(difftime(anytime(dat7$dates/1000), anytime(dat7$unixdates/1000), units = "days"))

dat8 <- dat7 %>% filter(diff_days <= 0)

dat9 <- dat8 %>% mutate(
  outbreak = case_when(
    SHPAPPGREG == 1 | SHPDENGRP == 1 | SHPAPPTRAN == 1 ~ 1,
    TRUE ~ 0
  )
)



dat9$outbreak <- as.factor(dat9$outbreak)
names(dat9)

summary(dat9$NDVIs)
g1 <- ggplot(dat9,aes(x=diff_days,y=NDVIs/1000,color=outbreak)) + geom_smooth() + ylab("NDVIs") + xlab("Days (0 = observation date)") +
  scale_color_manual(name = "", labels = c("Non-outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62")) + theme_pubr()

g2 <- ggplot(dat9,aes(x=diff_days,y=dissimilarity,color=outbreak)) + geom_smooth() + ylab("Dissimilarity") + xlab("Days (0 = observation date)") +
  scale_color_manual(name = "", labels = c("Non-outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62")) + theme_pubr()

g3 <- ggplot(dat9,aes(x=diff_days,y=contrasts,color=outbreak)) + geom_smooth() + ylab("contrasts") + xlab("Days (0 = observation date)") +
  scale_color_manual(name = "", labels = c("Non-outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62")) + theme_pubr()

g4 <-ggplot(dat9,aes(x=diff_days,y=entropy,color=outbreak)) + geom_smooth() + ylab("entropy") + xlab("Days (0 = observation date)") +
  scale_color_manual(name = "", labels = c("Non-outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62")) + theme_pubr()


names(dat9)


mod <- bam(outbreak~ te(diff_days,((NDVIs))) + te(diff_days,((entropy))) +
             te(diff_days,((dissimilarity))) + te(diff_days,((contrasts))) + 
             te(Longitude,Latitude), family=binomial(),select = TRUE,discrete=TRUE,nthreads=4,data=dat9)



summary(mod)
draw(mod)
appraise(mod)
k.check(mod)

dat9$preds <- predict(mod,newdata = dat9,type="response")



ggplot(dat9,aes(x=diff_days,y=NDVIs/1000,z=(preds))) +
  stat_summary_2d() + ylab("")  +
  scale_fill_viridis() + ylim(0,.25)

ggplot(dat9,aes(x=diff_days,y=dissimilarity,z=(preds))) +
  stat_summary_2d() + ylab("")  +
  scale_fill_viridis()+ ylim(0,75)

ggplot(dat9,aes(x=diff_days,y=contrasts,z=(preds))) +
  stat_summary_2d() + ylab("")  +
  scale_fill_viridis() + ylim(0,15000)

ggplot(dat9,aes(x=diff_days,y=entropy,z=(preds))) +
  stat_summary_2d() + ylab("")  +
  scale_fill_viridis() + ylim(2,5)



ggplot(dat9,aes(x=diff_days,y=NDVIs/1000,color=SHPAPPGREG)) + geom_smooth() + ylab("NDVIs") + xlab("Days (0 = observation date)") +
  scale_color_manual(name = "", labels = c("Non-outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62"))

ggplot(dat9,aes(x=diff_days,y=dissimilarity,color=SHPAPPGREG)) + geom_smooth() + ylab("Dissimilarity") + xlab("Days (0 = observation date)") +
  scale_color_manual(name = "", labels = c("Non-outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62"))

ggplot(dat9,aes(x=diff_days,y=contrasts,color=SHPAPPGREG)) + geom_smooth() + ylab("contrasts") + xlab("Days (0 = observation date)") +
  scale_color_manual(name = "", labels = c("Non-outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62"))

ggplot(dat9,aes(x=diff_days,y=entropy,color=SHPAPPGREG)) + geom_smooth() + ylab("entropy") + xlab("Days (0 = observation date)") +
  scale_color_manual(name = "", labels = c("Non-outbreak", "Outbreak"),values=c("#67a9cf","#ef8a62"))

p1 <- ggplotGrob(g1)
p2 <- ggplotGrob(g2)
p3 <- ggplotGrob(g3)
p4 <- ggplotGrob(g4)

lay <- rbind(c(1,2),
             c(3,4))

grid.arrange(
  grobs = list(p1,p2,p3,p4),
  layout_matrix = lay)
