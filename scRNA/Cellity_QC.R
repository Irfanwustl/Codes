library(cellity)
library(org.Hs.eg.db)
data("extra_human_genes")
data("feature_info")
GO_terms <- feature_info[[1]]
common_features <- feature_info[[2]]
counts<-read.table("HTSEQ_Count/Count_HPAP001.txt",h=T,row.names = 1)
stats<-read.table("Stats_Celloine/SC_HPAP001",h=T)
stats<-stats[,-12]
counts = counts[which(rowSums(counts) > 0),which(colSums(counts) > 0)]
counts2 = counts[which(rowSums(counts) > 3000),which(colSums(counts) > 1000)]
counts_nm <- normalise_by_factor(counts2, colSums(counts2))
row.names(stats)<-stats[,1]
I<-colnames(counts2)
stats2<-stats[I,]
sample_features <- extract_features(counts_nm, stats2,common_features = common_features,GO_terms = GO_terms, extra_genes = extra_human_genes,organism = "human")
sample_features_all <- sample_features[[1]]
training_quality_PCA_allF <- assess_cell_quality_PCA(sample_features_all, file = "./HPAP001_QC_PCA_allF.pdf")
write.csv(training_quality_PCA_allF,"HPAP001.PCA_allF.csv")
sample_features_all <- sample_features[[2]]
training_quality_PCA_allF <- assess_cell_quality_PCA(sample_features_all, file = "./HPAP001_QC_PCA_comF.pdf")
write.csv(training_quality_PCA_allF,"HPAP001.PCA_comF.csv")
write.csv(sample_features[[1]],"HPAP001.All_Features.csv")
