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


library(gplots)

# create some data
d <- matrix(rnorm(120),12,10)

# cluster it
hr <- hclust(as.dist(1-cor(t(d), method="pearson")), method="complete")

# define some clusters
mycl <- cutree(hr, h=max(hr$height/1.5))

# get a color palette equal to the number of clusters
clusterCols <- rainbow(length(unique(mycl)))

# create vector of colors for side bar
myClusterSideBar <- clusterCols[mycl]

# choose a color palette for the heat map
myheatcol <- rev(redgreen(75))

# draw the heat map
heatmap.2(d, main="Hierarchical Cluster", Rowv=as.dendrogram(hr), Colv=NA, dendrogram="row", scale="row", col=myheatcol, density.info="none", trace="none", RowSideColors= myClusterSideBar)

# cutree returns a vector of cluster membership
# in the order of the original data rows
# examine it
mycl

# examine the cluster membership by it's order
# in the heatmap
mycl[hr$order]

# grab a cluster
cluster1 <- d[mycl == 1,]

# or simply add the cluster ID to your data
foo <- cbind(d, clusterID=mycl)

# examine the data with cluster ids attached, and ordered like the heat map
foo[hr$order,]
