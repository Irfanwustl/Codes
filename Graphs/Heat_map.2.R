library(gplots)
library("RColorBrewer")

#read file
d<-read.table("normalized_by_TMM.tsv",h=T,row.names=1)

#read gene list
gene<-read.table("top200.txt",h=F,row.names=1)

#keep interested rows
list_gene<- which(row.names(d) %in% row.names(gene))
d<-as.matrix(d[list_gene,])

# cluster it
hr <- hclust(as.dist(1-cor(t(d), method="pearson")), method="complete")

# define some clusters
mycl <- cutree(hr, h=max(hr$height/1.5))

# get a color palette equal to the number of clusters
clusterCols <- rainbow(length(unique(mycl)))

# create vector of colors for side bar
myClusterSideBar <- clusterCols[mycl]

# choose a color palette for the heat map
mypalette <- rev(redgreen(75))

# draw the heat map
heatmap.2(dv, main="DE Analysis\nCondition1 vs Condition2", Rowv=as.dendrogram(hr), Colv=NA ,dendrogram="row", scale="row", col=mypalette, density.info="none", trace="none", lhei=c(1,2.5), lwid=c(0.5,2.5), keysize=0.75, key.par = list(cex=0.5),margins=c(8,8),RowSideColors= myClusterSideBar)
