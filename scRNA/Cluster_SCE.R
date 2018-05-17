library(scater)
library(cellity)
counts<-read.table("Subsets/HPAP013.tsv",h=T,row.names = 1)
pca<-read.csv("HPAP013.PCA_allF.csv",h=T,row.names = 1)
row.names(pca)<-pca[,1]
I<-pca[which(pca$quality == "1"),]
l_i <- which(colnames(counts) %in% row.names(I))
counts2<-counts[,l_i]
counts2 = counts2[which(rowSums(counts2) > 0),which(colSums(counts2) > 0)]
counts2<-as.matrix(counts2)
sce <- SingleCellExperiment(assays = list(counts = counts2,logcounts=log2(counts2+1)), colData = I)
pdf("HPAP013_AllF.pdf")
plotPCA(sce, colour_by = "quality")
dev.off()

Marker<-read.table("Markers.txt")
l_gene <- which(row.names(counts2) %in% Marker[,1])
counts3<-counts2[l_gene,]
sce <- SingleCellExperiment(assays = list(counts = counts3,logcounts=log2(counts3+1)), colData = I)
pdf("HPAP013_Marker_AllF.pdf")
plotPCA(sce, colour_by = "quality")
dev.off()

pca<-read.csv("HPAP013.PCA_comF.csv",h=T,row.names = 1)
row.names(pca)<-pca[,1]
I<-pca[which(pca$quality == "1"),]
l_i <- which(colnames(counts) %in% row.names(I))
counts2<-counts[,l_i]
counts2 = counts2[which(rowSums(counts2) > 0),which(colSums(counts2) > 0)]
counts2<-as.matrix(counts2)
sce <- SingleCellExperiment(assays = list(counts = counts2,logcounts=log2(counts2+1)), colData = I)
pdf("HPAP013_ComF.pdf")
plotPCA(sce, colour_by = "quality")
dev.off()

Marker<-read.table("Markers.txt")
l_gene <- which(row.names(counts2) %in% Marker[,1])
counts3<-counts2[l_gene,]
sce <- SingleCellExperiment(assays = list(counts = counts3,logcounts=log2(counts3+1)), colData = I)
pdf("HPAP013_Marker_comF.pdf")
plotPCA(sce, colour_by = "quality")
dev.off()
