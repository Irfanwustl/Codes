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



library("plotly")
diff_df["group"] <- "NotSignificant"
diff_df[which(diff_df['PValue'] < 0.05 & abs(diff_df['logFC']) < 1.5 ),"group"] <- "Significant"
diff_df[which(diff_df['PValue'] > 0.05 & abs(diff_df['logFC']) > 1.5 ),"group"] <- "FoldChange"
diff_df[which(diff_df['PValue'] < 0.05 & abs(diff_df['logFC']) > 1.5 ),"group"] <- "Significant&FoldChange"
top_peaks <- diff_df[with(diff_df, order(diff_df$logFC, diff_df$PValue)),][1:5,]
top_peaks <- rbind(top_peaks, diff_df[with(diff_df, order(-diff_df$logFC, diff_df$PValue)),][1:5,])
a <- list()
for (i in seq_len(nrow(top_peaks))) {
    m <- top_peaks[i, ]
    a[[i]] <- list(
    x = m[["logFC"]],
    y = -log10(m[["PValue"]]),
    text = m[["X"]],
    xref = "x",
    yref = "y",
    showarrow = TRUE,
    arrowhead = 0.5,
    ax = 20,
    ay = -40
    )
}
p <- plot_ly(data = diff_df, x = diff_df$logFC, y = -log10(diff_df$PValue), text = diff_df$X, mode = "markers", color = diff_df$group) %>%
layout(title ="Clinical Outcomes: Stable vs Responsive",titlefont = list(size=40)) %>%
layout(annotations = a,xaxis = list(title="logFC\n<0:Enriched in Responsive, >0:Enriched in Stable)",titlefont=list(size=40)), yaxis = list(title="-log10(PValue)",titlefont = list(size=40)),margin =list(l = 60, r = 60, b = 120, t = 120, pad = 5))
p




library(Rtsne)
train<-
Labels<-train$Gene.ID
train$label<-as.factor(train$Gene.ID)
colors = rainbow(length(unique(train$Gene.ID)))
names(colors) = unique(train$Gene.ID)
tsne <- Rtsne(train[,-1], dims = 2, perplexity=17, verbose=TRUE, max_iter = 500)
exeTimeTsne<- system.time(Rtsne(train[,-1], dims = 2, perplexity=17, verbose=TRUE, max_iter = 500))
## Plotting
plot(tsne$Y, t='n', main="tsne")
text(tsne$Y, labels=train$Gene.ID, col=colors[train$Gene.ID])




library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(gplots)
library(calibrate)
library(data.table)
##set data table as matrix and fix names
t.names <- paste(dat.d.3.t$T,dat.d.3.t$count,sep=',')
dat.d.3.t <- as.matrix(subset(dat.d.3.t,select=-c(T,count)))
rownames(dat.d.3.t) <- t.names
##create correlation matrix and pvalue matrix
num.t <- nrow(dat.d.3.t)
cor.mat <- mat.or.vec(num.t,num.t)
cor.p <- mat.or.vec(num.t,num.t)
for (i in 1:num.t){
    for (j in 1:num.t){
        r1 <- dat.d.3.t[i,]
        r2 <- dat.d.3.t[j,]
        c <- cor.test(r1,r2,alternative='two.sided',method='pearson')
        cor.mat[i,j] <- c$estimate
        cor.p[i,j] <- c$p.value
    }
}
rownames(cor.mat) <- t.names
colnames(cor.mat) <- t.names
colors <- rev(colorRampPalette(brewer.pal(8,'RdBu'))(20))
heatmap.2(cor.mat,trace='none', Rowv=FALSE, Colv=FALSE,
col=colors,
cellnote=ifelse(cor.p<=2e-16,'***',ifelse(cor.p<5e-12,'**',ifelse(cor.p<5e-6,'*',''))),notecol='white',notecex=1,
cexRow=1.5, cexCol=1.5,
#margins=c(20,20),
keysize=0.5,lwid=c(0.3,3),lhei=c(0.3,3)
)





cor(Q_HTG, use ="everything")






| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/

cat Angiogenesis.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/Angiogenesis.txt
cat Apoptosis.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/Apoptosis.txt
cat Cardio_Toxicity.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/Cardio_Toxicity.txt
cat Cell_Cycle.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/Cell_Cycle.txt
cat Cluster_of_Differentiation.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/Cluster_of_Differentiation.txt
cat Controls.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/Controls.txt
cat DMPK.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/DMPK.txt
cat DNA_Repair.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/DNA_Repair.txt
cat Drug_Targets.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/Drug_Targets.txt
cat EGFR_HER_Pathway.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/EGFR_HER_Pathway.txt
cat EGF_PDGF_Pathway.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/EGF_PDGF_Pathway.txt
cat FGFR_Pathway.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/FGFR_Pathway.txt
cat HK.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/HK.txt
cat Hedgehog_Pathway.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/Hedgehog_Pathway.txt
cat Housekeepers.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/Housekeepers.txt
cat Hypoxia.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/Hypoxia.txt
cat IO.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/IO.txt
cat Immuno_Oncology.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/Immuno_Oncology.txt
cat JAK_STAT_Pathway.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/JAK_STAT_Pathway.txt
cat List.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/List.txt
cat Lymphoma_Genes.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/Lymphoma_Genes.txt
cat MAP_Kinase_Pathway.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/MAP_Kinase_Pathway.txt
cat Markers.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/Markers.txt
cat NFkB_Pathway.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/NFkB_Pathway.txt
cat Other_Genes_of_Interest.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/Other_Genes_of_Interest.txt
cat PI3K_AKT_Pathway.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/PI3K_AKT_Pathway.txt
cat Stem_cells.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/Stem_cells.txt
cat Stress_Toxity.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/Stress_Toxity.txt
cat Tissue_Specific_Genes.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/Tissue_Specific_Genes.txt
cat WNT_Pathway.txt| sed 's/VLP139/NSCLC_VLP139/g' | sed 's/VLP151/NSCLC_VLP151/g' | sed 's/VLP168_18752/NSCLC_VLP168_18752/g' | sed 's/VLP168_25118/NSCLC_VLP168_25118/g' | sed 's/\tVLP168/\tPDAC_VLP168/g' > New/WNT_Pathway.txt

cat



cut -f 1-795 ../normalized_by_TMM.tsv    >    HPAP001.tsv &
cut -f 1,796-1595 ../normalized_by_TMM.tsv    >    HPAP002.tsv &
cut -f 1,1596-2240 ../normalized_by_TMM.tsv    >    HPAP003.tsv &
cut -f 1,2241-3000 ../normalized_by_TMM.tsv    >    HPAP004.tsv &
cut -f 1,3001-3654 ../normalized_by_TMM.tsv    >    HPAP005.tsv &
cut -f 1,3655-4454 ../normalized_by_TMM.tsv    >    AEBK009.tsv &
cut -f 1,4455-5254 ../normalized_by_TMM.tsv    >    HPAP006.tsv &
cut -f 1,5255-6054 ../normalized_by_TMM.tsv    >    HPAP007.tsv &
cut -f 1,6055-6854 ../normalized_by_TMM.tsv    >    HPAP008.tsv &
cut -f 1,6855-7654 ../normalized_by_TMM.tsv    >    HPAP009.tsv &
cut -f 1,7655-8334 ../normalized_by_TMM.tsv    >    HPAP010.tsv &
cut -f 1,8335-9125 ../normalized_by_TMM.tsv    >    HPAP011.tsv &
cut -f 1,9126-9925 ../normalized_by_TMM.tsv    >    HPAP012.tsv &
cut -f 1,9926-10725 ../normalized_by_TMM.tsv    >    HPAP013.tsv &
cut -f 1,10726-11427 ../normalized_by_TMM.tsv    >    HPAP014.tsv &
cut -f 1,11428-12227 ../normalized_by_TMM.tsv    >    HPAP015.tsv &
cut -f 1,12228-12587 ../normalized_by_TMM.tsv    >    HPAP016.tsv &
cut -f 1,12588-12906 ../normalized_by_TMM.tsv    >    HPAP017.tsv &
cut -f 1,12588-12907 ../normalized_by_TMM.tsv    >    HPAP017.tsv &




python2.7 /project/ibilab/tools/bin/htseq-count -f bam -r pos -a 0 -s no sorted_bam_files /project/ibilab/library/genes/GTCh38/Homo_sapiens.GRCh38.90.gtf > output_counts
python2.7 /project/ibilab/tools/bin/htseq-count -f bam -r pos -a 0 -s no ../mapped_reads/sc43690.bam /project/ibilab/library/genes/GTCh38/Homo_sapiens.GRCh38.90_revised_ERCC.gtf > sc43690_counts



commands.append("python")
commands.append(os.path.dirname(sys.argv[0])+ "/lib/stats.py") # change  PATH
commands.append("-i")
commands.append(sorted_sam_files)
commands.append("-o")
commands.append(stats_output)
commands.append("-n")
commands.append(sample_name + "_{#}.counts")
commands.append("-g")
commands.append(GTF_file)

python /lib/stats.py -i ../mapped_reads/sc43690.sam -o sc43690_stats -n sc43690.counts -g /project/ibilab/library/genes/GTCh38/Homo_sapiens.GRCh38.90_revised_ERCC.gtf



TT10 TT11 曾經有error。要檢查sam檔跟bam有沒有出入


cat Header_HPAP004 HTSEQ_HPAP004.txt > Count_HPAP004.txt
cat Header_HPAP005 HTSEQ_HPAP005.txt > Count_HPAP005.txt
cat Header_HPAP006 HTSEQ_HPAP006.txt > Count_HPAP006.txt
cat Header_HPAP007 HTSEQ_HPAP007.txt > Count_HPAP007.txt
cat Header_HPAP008 HTSEQ_HPAP008.txt > Count_HPAP008.txt
cat Header_HPAP009 HTSEQ_HPAP009.txt > Count_HPAP009.txt
cat Header_HPAP010 HTSEQ_HPAP010.txt > Count_HPAP010.txt
cat Header_HPAP011 HTSEQ_HPAP011.txt > Count_HPAP011.txt
cat Header_HPAP012 HTSEQ_HPAP012.txt > Count_HPAP012.txt
cat Header_HPAP013 HTSEQ_HPAP013.txt > Count_HPAP013.txt
cat Header_HPAP014 HTSEQ_HPAP014.txt > Count_HPAP014.txt


out <- queryMany(Genes, scopes="symbol", fields="ensembl.gene", species="human")


cut -f 1 HPAP015.lst  |sed ':a;N;$!ba;s/\n/ /g'  >Header_HPAP015






HPAP001<-read.csv("HPAP001.All_Features.csv",row.names=1)
HPAP002<-read.csv("HPAP002.All_Features.csv",row.names=1)
HPAP003<-read.csv("HPAP003.All_Features.csv",row.names=1)
HPAP004<-read.csv("HPAP004.All_Features.csv",row.names=1)
HPAP005<-read.csv("HPAP005.All_Features.csv",row.names=1)
HPAP006<-read.csv("HPAP006.All_Features.csv",row.names=1)
HPAP007<-read.csv("HPAP007.All_Features.csv",row.names=1)
HPAP008<-read.csv("HPAP008.All_Features.csv",row.names=1)
HPAP009<-read.csv("HPAP009.All_Features.csv",row.names=1)
HPAP010<-read.csv("HPAP010.All_Features.csv",row.names=1)
HPAP011<-read.csv("HPAP011.All_Features.csv",row.names=1)
HPAP012<-read.csv("HPAP012.All_Features.csv",row.names=1)
HPAP013<-read.csv("HPAP013.All_Features.csv",row.names=1)
HPAP014<-read.csv("HPAP014.All_Features.csv",row.names=1)
HPAP015<-read.csv("HPAP015.All_Features.csv",row.names=1)
HPAP016<-read.csv("HPAP016.All_Features.csv",row.names=1)
HPAP017<-read.csv("HPAP017.All_Features.csv",row.names=1)
AEBK009<-read.csv("AEBK009.All_Features.csv",row.names=1)
HPAP001<-read.csv("HPAP001.All_Features.csv",row.names=1)



setwd("~/BIC/Klaus_Kaestner_HPAP_Jonathan_2017_9/data/F800HT/QC")
HPAP001<-read.csv("HPAP001.All_Features.csv",row.names=1)
HPAP002<-read.csv("HPAP002.All_Features.csv",row.names=1)
HPAP003<-read.csv("HPAP003.All_Features.csv",row.names=1)
HPAP004<-read.csv("HPAP004.All_Features.csv",row.names=1)
HPAP005<-read.csv("HPAP005.All_Features.csv",row.names=1)
HPAP006<-read.csv("HPAP006.All_Features.csv",row.names=1)
HPAP007<-read.csv("HPAP007.All_Features.csv",row.names=1)
HPAP008<-read.csv("HPAP008.All_Features.csv",row.names=1)
HPAP009<-read.csv("HPAP009.All_Features.csv",row.names=1)
HPAP010<-read.csv("HPAP010.All_Features.csv",row.names=1)
HPAP011<-read.csv("HPAP011.All_Features.csv",row.names=1)
HPAP012<-read.csv("HPAP012.All_Features.csv",row.names=1)
HPAP013<-read.csv("HPAP013.All_Features.csv",row.names=1)
HPAP014<-read.csv("HPAP014.All_Features.csv",row.names=1)
HPAP015<-read.csv("HPAP015.All_Features.csv",row.names=1)
HPAP016<-read.csv("HPAP016.All_Features.csv",row.names=1)
HPAP017<-read.csv("HPAP017.All_Features.csv",row.names=1)
AEBK009<-read.csv("AEBK009.All_Features.csv",row.names=1)
HPAP001<-read.csv("HPAP001.All_Features.csv",row.names=1)


AllF_PCA<-read.csv("ALL.PCA_allF.csv")
ComF_PCA<-read.csv("ALL.PCA_comF.csv")
AllF_PCA<-AllF_PCA[-c(1:17),]
ComF_PCA<-ComF_PCA[-c(1:17),]

row.names(AllF_PCA)<-AllF_PCA$cell
row.names(ComF_PCA)<-ComF_PCA$cell

IA<-AllF_PCA[which(AllF_PCA$quality == "1"),]
IC<-ComF_PCA[which(ComF_PCA$quality == "1"),]
IC0<-ComF_PCA[which(ComF_PCA$quality == "0"),]
IA0<-AllF_PCA[which(AllF_PCA$quality == "0"),]



l_A <- which(row.names(HPAP001) %in% row.names(IA))
l_C <- which(row.names(HPAP001) %in% row.names(IC))
l_A0 <- which(row.names(HPAP001) %in% row.names(IA0))
l_C0 <- which(row.names(HPAP001) %in% row.names(IC0))
HPAP001_A<-HPAP001[l_A,]
HPAP001_C<-HPAP001[l_C,]
HPAP001_A0<-HPAP001[l_A0,]
HPAP001_C0<-HPAP001[l_C0,]

l_A <- which(row.names(HPAP002) %in% row.names(IA))
l_C <- which(row.names(HPAP002) %in% row.names(IC))
l_A0 <- which(row.names(HPAP002) %in% row.names(IA0))
l_C0 <- which(row.names(HPAP002) %in% row.names(IC0))
HPAP002_A<-HPAP002[l_A,]
HPAP002_C<-HPAP002[l_C,]
HPAP002_A0<-HPAP002[l_A0,]
HPAP002_C0<-HPAP002[l_C0,]

l_A <- which(row.names(HPAP003) %in% row.names(IA))
l_C <- which(row.names(HPAP003) %in% row.names(IC))
l_A0 <- which(row.names(HPAP003) %in% row.names(IA0))
l_C0 <- which(row.names(HPAP003) %in% row.names(IC0))
HPAP003_A<-HPAP003[l_A,]
HPAP003_C<-HPAP003[l_C,]
HPAP003_A0<-HPAP003[l_A0,]
HPAP003_C0<-HPAP003[l_C0,]

l_A <- which(row.names(HPAP004) %in% row.names(IA))
l_C <- which(row.names(HPAP004) %in% row.names(IC))
l_A0 <- which(row.names(HPAP004) %in% row.names(IA0))
l_C0 <- which(row.names(HPAP004) %in% row.names(IC0))
HPAP004_A<-HPAP004[l_A,]
HPAP004_C<-HPAP004[l_C,]
HPAP004_A0<-HPAP004[l_A0,]
HPAP004_C0<-HPAP004[l_C0,]

l_A <- which(row.names(HPAP005) %in% row.names(IA))
l_C <- which(row.names(HPAP005) %in% row.names(IC))
l_A0 <- which(row.names(HPAP005) %in% row.names(IA0))
l_C0 <- which(row.names(HPAP005) %in% row.names(IC0))
HPAP005_A<-HPAP005[l_A,]
HPAP005_C<-HPAP005[l_C,]
HPAP005_A0<-HPAP005[l_A0,]
HPAP005_C0<-HPAP005[l_C0,]

l_A <- which(row.names(HPAP006) %in% row.names(IA))
l_C <- which(row.names(HPAP006) %in% row.names(IC))
l_A0 <- which(row.names(HPAP006) %in% row.names(IA0))
l_C0 <- which(row.names(HPAP006) %in% row.names(IC0))
HPAP006_A<-HPAP006[l_A,]
HPAP006_C<-HPAP006[l_C,]
HPAP006_A0<-HPAP006[l_A0,]
HPAP006_C0<-HPAP006[l_C0,]

l_A <- which(row.names(HPAP007) %in% row.names(IA))
l_C <- which(row.names(HPAP007) %in% row.names(IC))
l_A0 <- which(row.names(HPAP007) %in% row.names(IA0))
l_C0 <- which(row.names(HPAP007) %in% row.names(IC0))
HPAP007_A<-HPAP007[l_A,]
HPAP007_C<-HPAP007[l_C,]
HPAP007_A0<-HPAP007[l_A0,]
HPAP007_C0<-HPAP007[l_C0,]

l_A <- which(row.names(HPAP008) %in% row.names(IA))
l_C <- which(row.names(HPAP008) %in% row.names(IC))
l_A0 <- which(row.names(HPAP008) %in% row.names(IA0))
l_C0 <- which(row.names(HPAP008) %in% row.names(IC0))
HPAP008_A<-HPAP008[l_A,]
HPAP008_C<-HPAP008[l_C,]
HPAP008_A0<-HPAP008[l_A0,]
HPAP008_C0<-HPAP008[l_C0,]

l_A <- which(row.names(HPAP009) %in% row.names(IA))
l_C <- which(row.names(HPAP009) %in% row.names(IC))
l_A0 <- which(row.names(HPAP009) %in% row.names(IA0))
l_C0 <- which(row.names(HPAP009) %in% row.names(IC0))
HPAP009_A<-HPAP009[l_A,]
HPAP009_C<-HPAP009[l_C,]
HPAP009_A0<-HPAP009[l_A0,]
HPAP009_C0<-HPAP009[l_C0,]

l_A <- which(row.names(HPAP010) %in% row.names(IA))
l_C <- which(row.names(HPAP010) %in% row.names(IC))
l_A0 <- which(row.names(HPAP010) %in% row.names(IA0))
l_C0 <- which(row.names(HPAP010) %in% row.names(IC0))
HPAP010_A<-HPAP010[l_A,]
HPAP010_C<-HPAP010[l_C,]
HPAP010_A0<-HPAP010[l_A0,]
HPAP010_C0<-HPAP010[l_C0,]

l_A <- which(row.names(HPAP011) %in% row.names(IA))
l_C <- which(row.names(HPAP011) %in% row.names(IC))
l_A0 <- which(row.names(HPAP011) %in% row.names(IA0))
l_C0 <- which(row.names(HPAP011) %in% row.names(IC0))
HPAP011_A<-HPAP011[l_A,]
HPAP011_C<-HPAP011[l_C,]
HPAP011_A0<-HPAP011[l_A0,]
HPAP011_C0<-HPAP011[l_C0,]

l_A <- which(row.names(HPAP012) %in% row.names(IA))
l_C <- which(row.names(HPAP012) %in% row.names(IC))
l_A0 <- which(row.names(HPAP012) %in% row.names(IA0))
l_C0 <- which(row.names(HPAP012) %in% row.names(IC0))
HPAP012_A<-HPAP012[l_A,]
HPAP012_C<-HPAP012[l_C,]
HPAP012_A0<-HPAP012[l_A0,]
HPAP012_C0<-HPAP012[l_C0,]

l_A <- which(row.names(HPAP013) %in% row.names(IA))
l_C <- which(row.names(HPAP013) %in% row.names(IC))
l_A0 <- which(row.names(HPAP013) %in% row.names(IA0))
l_C0 <- which(row.names(HPAP013) %in% row.names(IC0))
HPAP013_A<-HPAP013[l_A,]
HPAP013_C<-HPAP013[l_C,]
HPAP013_A0<-HPAP013[l_A0,]
HPAP013_C0<-HPAP013[l_C0,]

l_A <- which(row.names(HPAP014) %in% row.names(IA))
l_C <- which(row.names(HPAP014) %in% row.names(IC))
l_A0 <- which(row.names(HPAP014) %in% row.names(IA0))
l_C0 <- which(row.names(HPAP014) %in% row.names(IC0))
HPAP014_A<-HPAP015[l_A,]
HPAP014_C<-HPAP015[l_C,]
HPAP014_A0<-HPAP015[l_A0,]
HPAP014_C0<-HPAP015[l_C0,]

l_A <- which(row.names(HPAP015) %in% row.names(IA))
l_C <- which(row.names(HPAP015) %in% row.names(IC))
l_A0 <- which(row.names(HPAP015) %in% row.names(IA0))
l_C0 <- which(row.names(HPAP015) %in% row.names(IC0))
HPAP015_A<-HPAP015[l_A,]
HPAP015_C<-HPAP015[l_C,]
HPAP015_A0<-HPAP015[l_A0,]
HPAP015_C0<-HPAP015[l_C0,]

l_A <- which(row.names(HPAP016) %in% row.names(IA))
l_C <- which(row.names(HPAP016) %in% row.names(IC))
l_A0 <- which(row.names(HPAP016) %in% row.names(IA0))
l_C0 <- which(row.names(HPAP016) %in% row.names(IC0))
HPAP016_A<-HPAP016[l_A,]
HPAP016_C<-HPAP016[l_C,]
HPAP016_A0<-HPAP016[l_A0,]
HPAP016_C0<-HPAP016[l_C0,]

l_A <- which(row.names(HPAP017) %in% row.names(IA))
l_C <- which(row.names(HPAP017) %in% row.names(IC))
l_A0 <- which(row.names(HPAP017) %in% row.names(IA0))
l_C0 <- which(row.names(HPAP017) %in% row.names(IC0))
HPAP017_A<-HPAP017[l_A,]
HPAP017_C<-HPAP017[l_C,]
HPAP017_A0<-HPAP017[l_A0,]
HPAP017_C0<-HPAP017[l_C0,]

l_A <- which(row.names(AEBK009) %in% row.names(IA))
l_C <- which(row.names(AEBK009) %in% row.names(IC))
l_A0 <- which(row.names(AEBK009) %in% row.names(IA0))
l_C0 <- which(row.names(AEBK009) %in% row.names(IC0))
AEBK009_A<-AEBK009[l_A,]
AEBK009_C<-AEBK009[l_C,]
AEBK009_A0<-AEBK009[l_A0,]
AEBK009_C0<-AEBK009[l_C0,]

l_A <- which(row.names(HPAP001) %in% row.names(IA))
l_C <- which(row.names(HPAP001) %in% row.names(IC))
l_A0 <- which(row.names(HPAP001) %in% row.names(IA0))
l_C0 <- which(row.names(HPAP001) %in% row.names(IC0))
HPAP001_A<-HPAP001[l_A,]
HPAP001_C<-HPAP001[l_C,]
HPAP001_A0<-HPAP001[l_A0,]
HPAP001_C0<-HPAP001[l_C0,]





write.csv(summary(HPAP001),"Summary_HPAP001.csv")
write.csv(summary(HPAP001_A),"Summary_HPAP001_A.csv")
write.csv(summary(HPAP001_C),"Summary_HPAP001_C.csv")
write.csv(summary(HPAP001_A0),"Summary_HPAP001_A0.csv")
write.csv(summary(HPAP001_C0),"Summary_HPAP001_C0.csv")

write.csv(summary(HPAP002),"Summary_HPAP002.csv")
write.csv(summary(HPAP002_A),"Summary_HPAP002_A.csv")
write.csv(summary(HPAP002_C),"Summary_HPAP002_C.csv")
write.csv(summary(HPAP002_A0),"Summary_HPAP002_A0.csv")
write.csv(summary(HPAP002_C0),"Summary_HPAP002_C0.csv")

write.csv(summary(HPAP003),"Summary_HPAP003.csv")
write.csv(summary(HPAP003_A),"Summary_HPAP003_A.csv")
write.csv(summary(HPAP003_C),"Summary_HPAP003_C.csv")
write.csv(summary(HPAP003_A0),"Summary_HPAP003_A0.csv")
write.csv(summary(HPAP003_C0),"Summary_HPAP003_C0.csv")

write.csv(summary(HPAP004),"Summary_HPAP004.csv")
write.csv(summary(HPAP004_A),"Summary_HPAP004_A.csv")
write.csv(summary(HPAP004_C),"Summary_HPAP004_C.csv")
write.csv(summary(HPAP004_A0),"Summary_HPAP004_A0.csv")
write.csv(summary(HPAP004_C0),"Summary_HPAP004_C0.csv")

write.csv(summary(HPAP005),"Summary_HPAP005.csv")
write.csv(summary(HPAP005_A),"Summary_HPAP005_A.csv")
write.csv(summary(HPAP005_C),"Summary_HPAP005_C.csv")
write.csv(summary(HPAP005_A0),"Summary_HPAP005_A0.csv")
write.csv(summary(HPAP005_C0),"Summary_HPAP005_C0.csv")

write.csv(summary(HPAP006),"Summary_HPAP006.csv")
write.csv(summary(HPAP006_A),"Summary_HPAP006_A.csv")
write.csv(summary(HPAP006_C),"Summary_HPAP006_C.csv")
write.csv(summary(HPAP006_A0),"Summary_HPAP006_A0.csv")
write.csv(summary(HPAP006_C0),"Summary_HPAP006_C0.csv")

write.csv(summary(HPAP007),"Summary_HPAP007.csv")
write.csv(summary(HPAP007_A),"Summary_HPAP007_A.csv")
write.csv(summary(HPAP007_C),"Summary_HPAP007_C.csv")
write.csv(summary(HPAP007_A0),"Summary_HPAP007_A0.csv")
write.csv(summary(HPAP007_C0),"Summary_HPAP007_C0.csv")

write.csv(summary(HPAP008),"Summary_HPAP008.csv")
write.csv(summary(HPAP008_A),"Summary_HPAP008_A.csv")
write.csv(summary(HPAP008_C),"Summary_HPAP008_C.csv")
write.csv(summary(HPAP008_A0),"Summary_HPAP008_A0.csv")
write.csv(summary(HPAP008_C0),"Summary_HPAP008_C0.csv")

write.csv(summary(HPAP009),"Summary_HPAP009.csv")
write.csv(summary(HPAP009_A),"Summary_HPAP009_A.csv")
write.csv(summary(HPAP009_C),"Summary_HPAP009_C.csv")
write.csv(summary(HPAP009_A0),"Summary_HPAP009_A0.csv")
write.csv(summary(HPAP009_C0),"Summary_HPAP009_C0.csv")

write.csv(summary(HPAP010),"Summary_HPAP010.csv")
write.csv(summary(HPAP010_A),"Summary_HPAP010_A.csv")
write.csv(summary(HPAP010_C),"Summary_HPAP010_C.csv")
write.csv(summary(HPAP010_A0),"Summary_HPAP010_A0.csv")
write.csv(summary(HPAP010_C0),"Summary_HPAP010_C0.csv")

write.csv(summary(HPAP011),"Summary_HPAP011.csv")
write.csv(summary(HPAP011_A),"Summary_HPAP011_A.csv")
write.csv(summary(HPAP011_C),"Summary_HPAP011_C.csv")
write.csv(summary(HPAP011_A0),"Summary_HPAP011_A0.csv")
write.csv(summary(HPAP011_C0),"Summary_HPAP011_C0.csv")

write.csv(summary(HPAP012),"Summary_HPAP012.csv")
write.csv(summary(HPAP012_A),"Summary_HPAP012_A.csv")
write.csv(summary(HPAP012_C),"Summary_HPAP012_C.csv")
write.csv(summary(HPAP012_A0),"Summary_HPAP012_A0.csv")
write.csv(summary(HPAP012_C0),"Summary_HPAP012_C0.csv")

write.csv(summary(HPAP013),"Summary_HPAP013.csv")
write.csv(summary(HPAP013_A),"Summary_HPAP013_A.csv")
write.csv(summary(HPAP013_C),"Summary_HPAP013_C.csv")
write.csv(summary(HPAP013_A0),"Summary_HPAP013_A0.csv")
write.csv(summary(HPAP013_C0),"Summary_HPAP013.csv")

write.csv(summary(HPAP014),"Summary_HPAP014.csv")
write.csv(summary(HPAP014_A),"Summary_HPAP014_A.csv")
write.csv(summary(HPAP014_C),"Summary_HPAP014_C.csv")
write.csv(summary(HPAP014_A0),"Summary_HPAP014_A0.csv")
write.csv(summary(HPAP014_C0),"Summary_HPAP014_C0.csv")

write.csv(summary(HPAP015),"Summary_HPAP015.csv")
write.csv(summary(HPAP015_A),"Summary_HPAP015_A.csv")
write.csv(summary(HPAP015_C),"Summary_HPAP015_C.csv")
write.csv(summary(HPAP015_A0),"Summary_HPAP015_A0.csv")
write.csv(summary(HPAP015_C0),"Summary_HPAP015_C0.csv")


write.csv(summary(HPAP016),"Summary_HPAP016.csv")
write.csv(summary(HPAP016_A),"Summary_HPAP016_A.csv")
write.csv(summary(HPAP016_C),"Summary_HPAP016_C.csv")
write.csv(summary(HPAP016_A0),"Summary_HPAP016_A0.csv")
write.csv(summary(HPAP016_C0),"Summary_HPAP016_C0.csv")


write.csv(summary(HPAP017),"Summary_HPAP017.csv")
write.csv(summary(HPAP017_A),"Summary_HPAP017_A.csv")
write.csv(summary(HPAP017_C),"Summary_HPAP017_C.csv")
write.csv(summary(HPAP017_A0),"Summary_HPAP017_A0.csv")
write.csv(summary(HPAP017_C0),"Summary_HPAP017_C0.csv")

write.csv(summary(AEBK009),"Summary_AEBK009.csv")
write.csv(summary(AEBK009_A),"Summary_AEBK009_A.csv")
write.csv(summary(AEBK009_C),"Summary_AEBK009_C.csv")
write.csv(summary(AEBK009_A0),"Summary_AEBK009_A0.csv")
write.csv(summary(AEBK009_C0),"Summary_AEBK009_C0.csv")



cor.test(Folate$WB_5.MTH, Folate$SRCASE..Case...Control..1..Case..2...Control,  method = "pearson", use = "complete.obs")
cor.test(Folate$WB_THF, Folate$SRCASE..Case...Control..1..Case..2...Control,  method = "pearson", use = "complete.obs")
cor.test(Folate$WB_5.10.MT, Folate$SRCASE..Case...Control..1..Case..2...Control,  method = "pearson", use = "complete.obs")
cor.test(Folate$Total.THFs.nM., Folate$SRCASE..Case...Control..1..Case..2...Control,  method = "pearson", use = "complete.obs")
cor.test(Folate$Ratio_5.MTHF.THF, Folate$SRCASE..Case...Control..1..Case..2...Control,  method = "pearson", use = "complete.obs")
cor.test(Folate$Ratio_THF.5.MTHF, Folate$SRCASE..Case...Control..1..Case..2...Control,  method = "pearson", use = "complete.obs")


cor.test(Folate_M$WB_5.MTH, Folate_M$SRCASE..Case...Control..1..Case..2...Control,  method = "pearson", use = "complete.obs")
cor.test(Folate_M$WB_THF, Folate_M$SRCASE..Case...Control..1..Case..2...Control,  method = "pearson", use = "complete.obs")
cor.test(Folate_M$WB_5.10.MT, Folate_M$SRCASE..Case...Control..1..Case..2...Control,  method = "pearson", use = "complete.obs")
cor.test(Folate_M$Total.THFs.nM., Folate_M$SRCASE..Case...Control..1..Case..2...Control,  method = "pearson", use = "complete.obs")
cor.test(Folate_M$Ratio_5.MTHF.THF, Folate_M$SRCASE..Case...Control..1..Case..2...Control,  method = "pearson", use = "complete.obs")
cor.test(Folate_M$Ratio_THF.5.MTHF, Folate_M$SRCASE..Case...Control..1..Case..2...Control,  method = "pearson", use = "complete.obs")

cor.test(Folate_F$WB_5.MTH, Folate_F$SRCASE..Case...Control..1..Case..2...Control,  method = "pearson", use = "complete.obs")
cor.test(Folate_F$WB_THF, Folate_F$SRCASE..Case...Control..1..Case..2...Control,  method = "pearson", use = "complete.obs")
cor.test(Folate_F$WB_5.10.MT, Folate_F$SRCASE..Case...Control..1..Case..2...Control,  method = "pearson", use = "complete.obs")
cor.test(Folate_F$Total.THFs.nM., Folate_F$SRCASE..Case...Control..1..Case..2...Control,  method = "pearson", use = "complete.obs")
cor.test(Folate_F$Ratio_5.MTHF.THF, Folate_F$SRCASE..Case...Control..1..Case..2...Control,  method = "pearson", use = "complete.obs")
cor.test(Folate_F$Ratio_THF.5.MTHF, Folate_F$SRCASE..Case...Control..1..Case..2...Control,  method = "pearson", use = "complete.obs")


cor.test(Folate$WB_5.MTH, Folate$SRCASE..Case...Control..1..Case..2...Control,  method = "spearman", use = "complete.obs")
cor.test(Folate$WB_THF, Folate$SRCASE..Case...Control..1..Case..2...Control,  method = "spearman", use = "complete.obs")
cor.test(Folate$WB_5.10.MT, Folate$SRCASE..Case...Control..1..Case..2...Control,  method = "spearman", use = "complete.obs")
cor.test(Folate$Total.THFs.nM., Folate$SRCASE..Case...Control..1..Case..2...Control,  method = "spearman", use = "complete.obs")
cor.test(Folate$Ratio_5.MTHF.THF, Folate$SRCASE..Case...Control..1..Case..2...Control,  method = "spearman", use = "complete.obs")
cor.test(Folate$Ratio_THF.5.MTHF, Folate$SRCASE..Case...Control..1..Case..2...Control,  method = "spearman", use = "complete.obs")


cor.test(Folate_M$WB_5.MTH, Folate_M$SRCASE..Case...Control..1..Case..2...Control,  method = "spearman", use = "complete.obs")
cor.test(Folate_M$WB_THF, Folate_M$SRCASE..Case...Control..1..Case..2...Control,  method = "spearman", use = "complete.obs")
cor.test(Folate_M$WB_5.10.MT, Folate_M$SRCASE..Case...Control..1..Case..2...Control,  method = "spearman", use = "complete.obs")
cor.test(Folate_M$Total.THFs.nM., Folate_M$SRCASE..Case...Control..1..Case..2...Control,  method = "spearman", use = "complete.obs")
cor.test(Folate_M$Ratio_5.MTHF.THF, Folate_M$SRCASE..Case...Control..1..Case..2...Control,  method = "spearman", use = "complete.obs")
cor.test(Folate_M$Ratio_THF.5.MTHF, Folate_M$SRCASE..Case...Control..1..Case..2...Control,  method = "spearman", use = "complete.obs")

cor.test(Folate_F$WB_5.MTH, Folate_F$SRCASE..Case...Control..1..Case..2...Control,  method = "spearman", use = "complete.obs")
cor.test(Folate_F$WB_THF, Folate_F$SRCASE..Case...Control..1..Case..2...Control,  method = "spearman", use = "complete.obs")
cor.test(Folate_F$WB_5.10.MT, Folate_F$SRCASE..Case...Control..1..Case..2...Control,  method = "spearman", use = "complete.obs")
cor.test(Folate_F$Total.THFs.nM., Folate_F$SRCASE..Case...Control..1..Case..2...Control,  method = "spearman", use = "complete.obs")
cor.test(Folate_F$Ratio_5.MTHF.THF, Folate_F$SRCASE..Case...Control..1..Case..2...Control,  method = "spearman", use = "complete.obs")
cor.test(Folate_F$Ratio_THF.5.MTHF, Folate_F$SRCASE..Case...Control..1..Case..2...Control,  method = "spearman", use = "complete.obs")


t.test(FDB[,"WB_5_MTH"][FDB$SR_Case == 1],FDB[,"WB_5_MTH"][FDB$SR_Case == 2])
t.test(FDB[,"WB_THF"][FDB$SR_Case == 1],FDB[,"WB_THF"][FDB$SR_Case == 2])
t.test(FDB[,"WB_5.10.MT"][FDB$SR_Case == 1],FDB[,"WB_5.10.MT"][FDB$SR_Case == 2])

t.test(FDB[,"WBnM_5.MTH"][FDB$SR_Case == 1],FDB[,"WBnM_5.MTH"][FDB$SR_Case == 2])
t.test(FDB[,"WBnM_THF"][FDB$SR_Case == 1],FDB[,"WBnM_THF"][FDB$SR_Case == 2])
t.test(FDB[,"WBnM_5.10.MT"][FDB$SR_Case == 1],FDB[,"WBnM_5.10.MT"][FDB$SR_Case == 2])

t.test(FDB[,"Plasma_5.MTH"][FDB$SR_Case == 1],FDB[,"Plasma_5.MTH"][FDB$SR_Case == 2])
t.test(FDB[,"Plasma_THF"][FDB$SR_Case == 1],FDB[,"Plasma_THF"][FDB$SR_Case == 2])
t.test(FDB[,"Plasma_5.10.MT"][FDB$SR_Case == 1],FDB[,"Plasma_5.10.MT"][FDB$SR_Case == 2])

t.test(FDB[,"FRBC_5.MTH"][FDB$SR_Case == 1],FDB[,"FRBC_5.MTH"][FDB$SR_Case == 2])
t.test(FDB[,"FRBC_THF"][FDB$SR_Case == 1],FDB[,"FRBC_THF"][FDB$SR_Case == 2])
t.test(FDB[,"FRBC_5.10.MT"][FDB$SR_Case == 1],FDB[,"FRBC_5.10.MT"][FDB$SR_Case == 2])

t.test(FDB[,"Ratio_5.MTHF.THF"][FDB$SR_Case == 1],FDB[,"Ratio_5.MTHF.THF"][FDB$SR_Case == 2])
t.test(FDB[,"Ratio_THF.5.MTHF"][FDB$SR_Case == 1],FDB[,"Ratio_THF.5.MTHF"][FDB$SR_Case == 2])
t.test(FDB[,"Total.THFs.nM."][FDB$SR_Case == 1],FDB[,"Total.THFs.nM."][FDB$SR_Case == 2])


> t.test(FDB[,"WB_5_MTH"][FDB$SR_Case == 1],FDB[,"WB_5_MTH"][FDB$SR_Case == 2])

Welch Two Sample t-test

data:  FDB[, "WB_5_MTH"][FDB$SR_Case == 1] and FDB[, "WB_5_MTH"][FDB$SR_Case == 2]
t = -0.22054, df = 394.89, p-value = 0.8256
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
-24.39310  19.47239
sample estimates:
mean of x mean of y
224.9429  227.4033

> t.test(FDB[,"WB_THF"][FDB$SR_Case == 1],FDB[,"WB_THF"][FDB$SR_Case == 2])

Welch Two Sample t-test

data:  FDB[, "WB_THF"][FDB$SR_Case == 1] and FDB[, "WB_THF"][FDB$SR_Case == 2]
t = 1.033, df = 362.42, p-value = 0.3023
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
-3.889655 12.497051
sample estimates:
mean of x mean of y
21.8932   17.5895

> t.test(FDB[,"WB_5.10.MT"][FDB$SR_Case == 1],FDB[,"WB_5.10.MT"][FDB$SR_Case == 2])

Welch Two Sample t-test

data:  FDB[, "WB_5.10.MT"][FDB$SR_Case == 1] and FDB[, "WB_5.10.MT"][FDB$SR_Case == 2]
t = 1.6763, df = 336.13, p-value = 0.09461
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
-0.5362141  6.7196004
sample estimates:
mean of x mean of y
10.56954   7.47785

>
> t.test(FDB[,"WBnM_5.MTH"][FDB$SR_Case == 1],FDB[,"WBnM_5.MTH"][FDB$SR_Case == 2])

Welch Two Sample t-test

data:  FDB[, "WBnM_5.MTH"][FDB$SR_Case == 1] and FDB[, "WBnM_5.MTH"][FDB$SR_Case == 2]
t = -0.22004, df = 394.89, p-value = 0.826
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
-53.08071  42.39489
sample estimates:
mean of x mean of y
489.5726  494.9155

> t.test(FDB[,"WBnM_THF"][FDB$SR_Case == 1],FDB[,"WBnM_THF"][FDB$SR_Case == 2])

Welch Two Sample t-test

data:  FDB[, "WBnM_THF"][FDB$SR_Case == 1] and FDB[, "WBnM_THF"][FDB$SR_Case == 2]
t = 1.0331, df = 362.41, p-value = 0.3023
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
-8.731373 28.057607
sample estimates:
mean of x mean of y
49.14112  39.47800

> t.test(FDB[,"WBnM_5.10.MT"][FDB$SR_Case == 1],FDB[,"WBnM_5.10.MT"][FDB$SR_Case == 2])

Welch Two Sample t-test

data:  FDB[, "WBnM_5.10.MT"][FDB$SR_Case == 1] and FDB[, "WBnM_5.10.MT"][FDB$SR_Case == 2]
t = 1.6765, df = 336.1, p-value = 0.09457
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
-1.17386 14.72046
sample estimates:
mean of x mean of y
23.1533   16.3800

>
> t.test(FDB[,"Plasma_5.MTH"][FDB$SR_Case == 1],FDB[,"Plasma_5.MTH"][FDB$SR_Case == 2])

Welch Two Sample t-test

data:  FDB[, "Plasma_5.MTH"][FDB$SR_Case == 1] and FDB[, "Plasma_5.MTH"][FDB$SR_Case == 2]
t = -0.8747, df = 378.01, p-value = 0.3823
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
-3.463107  1.330610
sample estimates:
mean of x mean of y
21.97646  23.04271

> t.test(FDB[,"Plasma_THF"][FDB$SR_Case == 1],FDB[,"Plasma_THF"][FDB$SR_Case == 2])

Welch Two Sample t-test

data:  FDB[, "Plasma_THF"][FDB$SR_Case == 1] and FDB[, "Plasma_THF"][FDB$SR_Case == 2]
t = 2.4212, df = 360.45, p-value = 0.01596
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
0.05738569 0.55388785
sample estimates:
mean of x mean of y
1.2217172 0.9160804

> t.test(FDB[,"Plasma_5.10.MT"][FDB$SR_Case == 1],FDB[,"Plasma_5.10.MT"][FDB$SR_Case == 2])

Welch Two Sample t-test

data:  FDB[, "Plasma_5.10.MT"][FDB$SR_Case == 1] and FDB[, "Plasma_5.10.MT"][FDB$SR_Case == 2]
t = 7.9872, df = 281.48, p-value = 3.554e-14
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
0.2515058 0.4160136
sample estimates:
mean of x mean of y
0.5181818 0.1844221

>
> t.test(FDB[,"FRBC_5.MTH"][FDB$SR_Case == 1],FDB[,"FRBC_5.MTH"][FDB$SR_Case == 2])

Welch Two Sample t-test

data:  FDB[, "FRBC_5.MTH"][FDB$SR_Case == 1] and FDB[, "FRBC_5.MTH"][FDB$SR_Case == 2]
t = 2.1973, df = 374.33, p-value = 0.02861
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
11.54791 208.13131
sample estimates:
mean of x mean of y
983.7306  873.8910

> t.test(FDB[,"FRBC_THF"][FDB$SR_Case == 1],FDB[,"FRBC_THF"][FDB$SR_Case == 2])

Welch Two Sample t-test

data:  FDB[, "FRBC_THF"][FDB$SR_Case == 1] and FDB[, "FRBC_THF"][FDB$SR_Case == 2]
t = 1.5049, df = 340.82, p-value = 0.1333
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
-8.542592 64.192320
sample estimates:
mean of x mean of y
100.38265  72.55779

> t.test(FDB[,"FRBC_5.10.MT"][FDB$SR_Case == 1],FDB[,"FRBC_5.10.MT"][FDB$SR_Case == 2])

Welch Two Sample t-test

data:  FDB[, "FRBC_5.10.MT"][FDB$SR_Case == 1] and FDB[, "FRBC_5.10.MT"][FDB$SR_Case == 2]
t = 2.1271, df = 315.39, p-value = 0.03419
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
1.255931 32.223168
sample estimates:
mean of x mean of y
47.02347  30.28392

>
> t.test(FDB[,"Ratio_5.MTHF.THF"][FDB$SR_Case == 1],FDB[,"Ratio_5.MTHF.THF"][FDB$SR_Case == 2])

Welch Two Sample t-test

data:  FDB[, "Ratio_5.MTHF.THF"][FDB$SR_Case == 1] and FDB[, "Ratio_5.MTHF.THF"][FDB$SR_Case == 2]
t = 1.2746, df = 197.97, p-value = 0.204
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
-19.96337  92.93070
sample estimates:
mean of x mean of y
75.05000  38.56633

> t.test(FDB[,"Ratio_THF.5.MTHF"][FDB$SR_Case == 1],FDB[,"Ratio_THF.5.MTHF"][FDB$SR_Case == 2])

Welch Two Sample t-test

data:  FDB[, "Ratio_THF.5.MTHF"][FDB$SR_Case == 1] and FDB[, "Ratio_THF.5.MTHF"][FDB$SR_Case == 2]
t = 0.45558, df = 381.56, p-value = 0.649
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
-0.09330237  0.14957904
sample estimates:
mean of x mean of y
0.1876429 0.1595045

> t.test(FDB[,"Total.THFs.nM."][FDB$SR_Case == 1],FDB[,"Total.THFs.nM."][FDB$SR_Case == 2])

Welch Two Sample t-test

data:  FDB[, "Total.THFs.nM."][FDB$SR_Case == 1] and FDB[, "Total.THFs.nM."][FDB$SR_Case == 2]
t = 2.8814, df = 374.13, p-value = 0.004187
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
48.97967 259.46546
sample estimates:
mean of x mean of y
1131.1301  976.9075






python2.7 /project/ibilab/tools/miniconda3/envs/py27/bin/geneBody_coverage.py -r /project/ibilab/library/annotation/hg19/hg19_refseq.bed12 -i Archives/QS_1008.bam -o /project/ibilab/projects/Jos_Melenhorst_RNASeq_Jos_2017_10/large_data/RseQC/Genebody_Coverage/QS_1008
