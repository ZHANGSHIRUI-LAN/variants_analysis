install.packages(c("FactoMineR", "factoextra"))

library(dplyr)
library("factoextra")
library("FactoMineR")
install.packages("RColorBrewer")
library(RColorBrewer)

devtools::install_github("kassambara/ggpubr")
library(ggpubr) 
library(scatterplot3d)

#clear memory
rm()
rm(list=ls()) 

setwd(" ")
data0 <- read.csv("Compound list POS_more 5 year.csv",header = T,encoding='UTF-8')

data_1<-data0[,-1]

data_X<-t(data_1)

colnames(data_X)<-data0$Compound

pca.data <- PCA(data_X, scale.unit = TRUE, graph = FALSE,ncp = 3)

fviz_eig(pca.data, addlabels = TRUE, ylim = c(0, 70))

fviz_pca_var(pca.data, col.var = "cos2",
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
             repel = TRUE) 

pca.data <- PCA(data_X, scale.unit = TRUE, graph = FALSE,ncp = 2)

a <- fviz_pca_ind(pca.data, col.ind = "cos2", 
                  gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"), 
                  repel = TRUE)

ggpar(a,
      title = "Principal Component Analysis",
      xlab = "PC1", ylab = "PCA2",
      legend.title = "Cos2", legend.position = "top",
      ggtheme = theme_minimal())


