library(dplyr)
library(ggplot2)
library(ggrepel)

results<-read.csv("filename.tsv",h=T)
results<-results[which(results['padj'] != 'NA'),]
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
