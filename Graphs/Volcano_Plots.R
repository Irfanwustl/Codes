library(dplyr)
library(ggplot2)
library(ggrepel)

results<-read.csv("filename.tsv",h=T)
results<-results[which(results['padj'] != 'NA'),]
  results["group"] <- "NotSignificant"
  results[which(results['padj'] < 0.05 & abs(results['log2FoldChange']) < 1.5 ),"group"] <- "Significant"
  results[which(results['padj'] > 0.05 & abs(results['log2FoldChange']) > 1.5 ),"group"] <- "FoldChange"
  results[which(results['padj'] < 0.05 & abs(results['log2FoldChange']) > 1.5 ),"group"] <- "Significant&FoldChange"
  top_peaks <- results[with(results, order(results$log2FoldChange, results$padj)),][1:10,]
  top_peaks <- rbind(top_peaks, results[with(results, order(-results$log2FoldChange, results$padj)),][1:10,])
  a <- list()
results = mutate(results, sig=ifelse(results$padj<0.05, "padj<0.05", "Not Sig"))
p= ggplot(results, aes(log2FoldChange, -log10(padj))) +
     geom_point(aes(col=group)) +
     scale_color_manual(values=c("pink", "black","green","red"))
png("volcano_1.png")
p+geom_text_repel(data=filter(top_peaks, top_peaks$padj<0.05), aes(label=X)) + ggtitle("Volcano Plot of Comparison \n v1 (Left) vs v9 (Right)") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()





p= ggplot(results, aes(log2FoldChange, -log10(padj))) +
+     geom_point(aes(col=sig)) +
+     scale_color_manual(values=c("red", "black"))

p = ggplot(results, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(col=group)) +
  scale_color_manual(values=c("col1", "col2","col3","col4"))
p

p+geom_text(data=filter(results, padj<0.05), aes(label=Gene))
p+geom_text_repel(data=filter(results, padj<0.05), aes(label=Gene))


p = ggplot(results, aes(log2FoldChange, -log10(padj))) +
     geom_point(aes(col=results$group)) +
     scale_color_manual(values=c("black","red", "green"))
> p+geom_text_repel(data=filter(top_peaks, top_peaks$padj<0.05), aes(label=X))
