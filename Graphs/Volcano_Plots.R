library(dplyr)
library(ggplot2)
library(ggrepel)

results<-read.table("filename.tsv",h=T)

results["group"] <- "NotSignificant"
results[which(results['PValue'] < 0.05 & abs(results['logFC']) < 1.5 ),"group"] <- "Significant"
results[which(results['PValue'] > 0.05 & abs(results['logFC']) > 1.5 ),"group"] <- "FoldChange"
results[which(results['PValue'] < 0.05 & abs(results['logFC']) > 1.5 ),"group"] <- "Significant&FoldChange"
top_peaks <- results[with(results, order(results$logFC, results$PValue)),][1:5,]
top_peaks <- rbind(top_peaks, results[with(results, order(-results$logFC, results$PValue)),][1:5,])
a <- list()



p = ggplot(results, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(col=group)) +
  scale_color_manual(values=c("col1", "col2","col3","col4"))
p

p+geom_text(data=filter(results, padj<0.05), aes(label=Gene))
p+geom_text_repel(data=filter(results, padj<0.05), aes(label=Gene))


p = ggplot(results, aes(logFC, -log10(PValue))) +
+     geom_point(aes(col=sig)) +
+     scale_color_manual(values=c("red", "green"))
> p+geom_text_repel(data=filter(top_peaks, top_peaks$PValue<0.05), aes(label=X))
