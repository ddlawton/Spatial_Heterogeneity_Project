########
# Conducting PCA on 
#  APLC Database
####
library(tidyverse)
library(FactoMineR)
library(factoextra)

dat <- read.csv("Data/Processed/APLC_database_hierarchical_data_Feb_8.csv")
names(dat)

dat.RS <- dat %>% select(correlation,dissimilarity,entropy,homogeneity) %>% drop_na()

res.pca <- prcomp(dat.RS, scale = TRUE)

fviz_eig(res.pca)

fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
