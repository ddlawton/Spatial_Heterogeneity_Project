########
# Conducting PCA on 
#  APLC Database
####
library(tidyverse)
library(FactoMineR)
library(factoextra)

dat <- read.csv("Data/Processed/APLC_database_hierarchical_data_Feb_8.csv")
dat <- dat  %>% drop_na()

dat.RS <- dat %>% select(correlation,dissimilarity,entropy,homogeneity)

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

# Helper function 
#::::::::::::::::::::::::::::::::::::::::
var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}
# Compute Coordinates
#::::::::::::::::::::::::::::::::::::::::
loadings <- res.pca$rotation
sdev <- res.pca$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 
head(var.coord[, 1:4])


var.cos2 <- var.coord^2
head(var.cos2[, 1:4])


# Compute contributions
#::::::::::::::::::::::::::::::::::::::::
comp.cos2 <- apply(var.cos2, 2, sum)
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
head(var.contrib[, 1:4])

dat$PC1 <- res.pca$x[,2]
dat$spp_code <- as.factor(dat$spp_code)
ggplot((dat %>% filter(spp_code == 'None' | spp_code =='CT')),aes(x=diff_days,y=PC1,color=as.factor(binary_outbreak))) +geom_smooth() + xlim(-60,0)

## Unsure if this is the correct way forward. I would like to make an index like Bateman and Merritt 2020...
