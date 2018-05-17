library(gplots)
library("RColorBrewer")
GSE4835020<-read.csv("GSE48350.20.csv",h=T)
row.names(GSE4835020)<-GSE4835020$X
GSE4835020_Matrix<-data.matrix(GSE4835020)
GSE4835020_Matrix<-GSE4835020_Matrix[,-1]
GSE4835020_Matrix<-GSE4835020_Matrix[,-1]
A<-GSE4835020_Matrix[1,]
A[A == 1] = "red"
A[A == 0] = "blue"
I = order(GSE4835020_Matrix[2,])
GSE4835020_Matrix = GSE4835020_Matrix[, I]
A=A[I]
I = order(A)
GSE4835020_Matrix = GSE4835020_Matrix[, I]
A=A[I]
GSE4835020_Matrix<-GSE4835020_Matrix[-1,]
GSE4835020_Matrix<-GSE4835020_Matrix[-2,]
heatmap.2(GSE4835020_Matrix,margins=c(5,10),trace="none",scale="row",main ="GSE48350.20\nHuman[Brain Region] from 20 cases & 43 controls", Rowv=NA,Colv=NA, dendrogram="both",labCol = FALSE,ColSideColors=A,col=bluered(4096))

