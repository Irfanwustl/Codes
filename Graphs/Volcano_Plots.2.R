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
     scale_color_manual(values=c("grey", "black","green","red"))
png("volcano_1.png",height=1200,width=1600)
p+geom_text_repel(data=filter(top_peaks, top_peaks$padj<0.05), aes(label=X)) + ggtitle("Volcano Plot of Comparison \n ICOS-enriched (Left) vs Mem1-enriched (Right)") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


results["group"] <- "NotSignificant"
results[which(results['padj'] < 0.05 & (results['log2FoldChange']) > 1.5 ),"group"] <- "196 Genes enriched in Mem1"
results[which(results['padj'] < 0.05 & (results['log2FoldChange']) < -1.5 ),"group"] <- "148 Genes enriched in ICOS"
gene_8<-read.table("scripts/REVIGO_8.txt",h=F,row.names=1)
l_gene<-which(results$X %in% row.names(gene_8))
top_peaks <-results[l_gene,]
p= ggplot(results, aes(log2FoldChange, -log10(padj))) +
     geom_point(aes(col=group)) +
     scale_color_manual(values=c("#0011d1", "#d10000","#4f4f4f"))

p+
geom_text_repel(data=filter(top_peaks, top_peaks$log2FoldChange>0), aes(label=X),nudge_x = 15,color="black" )+geom_text_repel(data=filter(top_peaks, top_peaks$log2FoldChange<0), aes(label=X),nudge_x = -15,color="black" ) + ggtitle("Volcano Plot of Comparison \n ICOS-enriched (Left) vs Mem1-enriched (Right)") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(), panel.border = element_rect(colour = "black"))



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
