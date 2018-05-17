#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

get_gset <- function(gse_id) {
  gset <- getGEO(gse_id, GSEMatrix =TRUE, getGPL=FALSE)
  return(gset)
}

disease = "AD"
ADlist = read.table("AD_list.txt", header = T, sep='\t', stringsAsFactors=F)

workdir = "."
gene_list <- toupper(c("TNF", "IL6", "HAMP", "Slc39a14", "TF", "Tfrc", "Slc11a2"))
GPLlist = read.table("GPL.txt", header = T, sep='\t', stringsAsFactors=F, row.names = 1)

#for(i in 1:length(row.names(ADlist))){
for(i in c(25)){
  gse_id = ADlist[i,1]
  gpl_id = ADlist[i,2]
  field =  ADlist[i,3]
  case =  ADlist[i,4]
  ctl =  ADlist[i,5]
  cat("Running... i=",i,"GSE=",gse_id, "GPL=",gpl_id)
  pdf(paste(disease,"_",i,"_",gse_id,".pdf",sep=''))
  par(mfrow=c(4,2))
  
  gset <- NULL
  attempt <- 1
  while( is.null(gset) && attempt <= 3 ) {
    attempt <- attempt + 1
    try(
      gset <- get_gset(gse_id)
    )
    if(is.null(gset)) Sys.sleep(10)
  } 
  
  for(idx in 1:length(gset)){  #one GSE can contains multiple GPLs
    gset1 = gset[[idx]]
    gpl_id1 = annotation(gset1) 
    if(gpl_id1 == gpl_id){
      samples = sampleNames(gset1)
      types = pData(phenoData(gset1))[,eval(field)]
      
      #download GPL annotation file if not exists already
      gpl = getGEO(gpl_id, destdir=".") 

      #match gene symbols with platform IDs
      for(j in 1:length(gene_list)){
        if(is.na(GPLlist[gpl_id, 2])){
          rowid = Table(gpl)$ID[Table(gpl)[,eval(GPLlist[gpl_id, 1])] == gene_list[j]]
        } else {
          geneTag1 = paste("^", gene_list[j], " ",GPLlist[gpl_id, 2], sep='')
          geneTag2 = paste(GPLlist[gpl_id, 2], gene_list[j], GPLlist[gpl_id, 2])
          geneTag3 = paste(GPLlist[gpl_id, 2], " ",gene_list[j],"$",sep='')
          geneTag4 = paste("^", gene_list[j], "$", sep='')
          geneTag = paste(geneTag1, geneTag2, geneTag3, geneTag4, sep="|")
          rowid = Table(gpl)$ID[grepl(geneTag, Table(gpl)[,eval(GPLlist[gpl_id, 1])])]
        }
        rowid = intersect(rowid, rownames(exprs(gset1)))
        if(length(rowid)==0) next
        exp = exprs(gset1)[as.character(rowid),]
        if(is.matrix(exp)){
          exp = colMeans(exp) #take means for all transcripts
        }
        
        #match GSE samples
        case_idx <- grep(case,types, ignore.case=TRUE)
        ctl_idx <- grep(ctl,types, ignore.case=TRUE)
        pval <- NULL
        tryCatch( {
          pval <- signif(wilcox.test(exp[case_idx],exp[ctl_idx])$p.value,4)
        }, error=function(e){})
        if(is.null(pval)) {
          pval=signif(t.test(exp[case_idx],exp[ctl_idx])$p.value,4)
          pval <- paste(pval,"(t.test)")
        }  
        title <- paste(gse_id, " Gene:",gene_list[j]," p=",pval,sep='')
        boxplot(exp[case_idx],exp[ctl_idx], main=title, xlab="Disease vs. Control", ylab="Expression Level")
      }
    }
  }
  dev.off()
}


