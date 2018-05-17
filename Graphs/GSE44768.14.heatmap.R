library(gplots)
library("RColorBrewer")

GSE<-read.csv("GSE44768.14.csv",h=T)
row.names(GSE)<-GSE$X
library("RColorBrewer")
GSE<-GSE[,3:232]
GSE_Matrix<-data.matrix(GSE)
heatmap.2(GSE_Matrix, margins=c(5,10),trace="none",labCol = FALSE,dendrogram="both",col=bluered(4096))
