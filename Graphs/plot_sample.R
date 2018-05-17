#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

disease = "PD"
gse_id = "GSE19587" 
type = c("Dorsal Motor Nucleus of the Vagus","Inferior Olivary Nucleus")
case = "Parkinson's disease"
ctl = "Control"

#gse_id="GSE8397"

workdir = "."
gene_list <- toupper(c("TNF", "IL6", "HAMP", "Slc39a14", "TF", "Tfrc", "Slc11a2"))

# load series and platform data from GEO
gset <- getGEO(gse_id, GSEMatrix =TRUE, getGPL=FALSE)

pdf(paste(gse_id,".pdf",sep=''))
par(mfrow=c(4,2))

for(idx in 1:length(gset)){  #one GSE can contains multiple GPLs
  gset1 = gset[[idx]]
  gpl_id = annotation(gset1) 
  samples = sampleNames(gset1)
  titles = pData(phenoData(gset1))$title

  #download GPL annotation file if not exists already
  gpl_file = paste(workdir,"/",gpl_id,".soft", sep='')
  if(!file.exists(gpl_file)) gpl = getGEO(gpl_id, destdir=".")
  
  #match gene symbols with platform IDs
  if(gpl_id == "GPL571"){
    for(i in 1:length(gene_list)){
      #Gene symbols are called by $"Gene Symbol"
      #gpl_idx = which(Table(gpl)$"Gene Symbol" == gene_list[i])
      rowid = Table(gpl)$ID[Table(gpl)$"Gene Symbol" == gene_list[i]]
      exp = exprs(gset1)[as.character(rowid),]
      if(is.matrix(exp)){
        exp = colMeans(exp) #take means for all transcripts
      }
      
      #match GSE samples
      for(t in 1:length(type)){
        case_idx <- grep(paste(type[t],".*",case,sep=''),titles, ignore.case=TRUE)
        ctl_idx <- grep(paste(type[t],".*",ctl,sep=''),titles, ignore.case=TRUE)
        title <- paste(disease, type[t], "(Gene:",gene_list[i],")")
        boxplot(exp[case_idx],exp[ctl_idx], main=title, xlab="Disease vs. Control", ylab="Expression Level")
      }    
    }
  }
}

dev.off()

