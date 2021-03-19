######
# Australian Plague Locust --
#  PCA analyses - Heterogenity Index Construction
#   
#####
rm(list=ls())

library(tidyverse)
library(data.table)
library(factoextra)
library(FactoMineR)

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

dat$NDVIs <- dat$NDVIs/1000

dat <- dat %>% drop_na(c(contrasts,correlation,dissimilarity,entropy,homogeneity,NDVIs)) %>% filter(between(NDVIs,0,1)) %>% filter(between(diff_days,-60,0))

PCA_dat <- dat %>% dplyr::select(8,9,11,12,13)

res.pca <- prcomp(PCA_dat, scale = TRUE)
fviz_eig(res.pca)

fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# Eigenvalues
eig.val <- get_eigenvalue(res.pca)
eig.val

# Results for Variables
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

 
PCA_Cords <- as.data.frame(res.var$coord)

PCA_Cords$Variance_explained <- eig.val$variance.percent
PCA_Cords$Cum_Variance_explained <- eig.val$cumulative.variance.percent

PCA_Cords <- PCA_Cords %>%
  mutate(
    Dim.1_weighted = (abs(Dim.1*58.7349491)),
    Dim.2_weighted = (abs(Dim.2*34.2182652))
    ) %>%
  mutate(
    weights = rowSums(.[8:9])  
    ) %>%
  dplyr::select(!c(3:9))






