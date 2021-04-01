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

mod3 <- readRDS("models/mod3.rds")
mod2 <- readRDS("models/mod2.rds")
mod <- readRDS("models/mod.rds")

k.check(mod3)
AIC(mod,mod3)
concurvity(mod)
summary(mod3)
draw(mod3)
appraise(mod3)
AIC(mod,mod2)
