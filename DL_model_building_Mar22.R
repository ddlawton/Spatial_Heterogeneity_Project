rm(list=ls())
library(parallel)
library(tidyverse)
library(anytime)
library(splitstackshape)
library(data.table)
library(mgcv)
library(lubridate)
library(gratia)
dat <- fread("Data/processed/DL_database_hierarchical_data_Feb_24.csv")
dat$NDVIs <- dat$NDVIs/1000

dat <- dat %>% drop_na(c(contrasts,correlation,dissimilarity,entropy,homogeneity,NDVIs)) %>% 
  filter(between(NDVIs,0,1)) %>% filter(between(diff_days,-78,0)) %>% 
  select(!c("contrasts_normalized", "dissimilarity_normalized","entropy_normalized",
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

dat <- dat %>% mutate(
  spatial_hetero = rowMeans(.[,41:45], na.rm = TRUE)
)

summary(dat$spatial_hetero)
str(dat)

mod1 <- readRDS("models/DL_mods/mod.rds")
summary(mod1)
draw(mod1)
appraise(mod1)
aov<-anova(mod,mod1)
summary(aov)
draw(mod)


mod2 <- bam(as.factor(binary_outbreak) ~ 
             te(diff_days, spatial_hetero, k = 5) + 
             te(ECO_NAME, Zone, bs=c("re","re")),
           family=binomial(),data=dat,select = TRUE,discrete=TRUE,nthreads=nc)
summary(mod2)
plot(mod2)
