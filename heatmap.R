library(RColorBrewer)
library(pheatmap)

setwd("")
data0 <- read.csv("Compound list POS_more 5 year.csv",header = T,encoding='UTF-8')


data3<-data.matrix(data0)

rownames(data3)<-data0$Compound

data3<-data3[,-(1:1)]

p<-pheatmap(data3,fontsize_row=8,fontsize_col=8,height=10000,cellheight=10,cellwidth =10,cluster_cols=T,cluster_rows=F)

p<-pheatmap(data3,color=colorRampPalette(brewer.pal(9, "Blues"))(20),fontsize_row=8,fontsize_col=8,height=10000,cellheight=10,cellwidth =10,cluster_cols=T,cluster_rows=F)



