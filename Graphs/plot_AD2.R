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
gene_list <- toupper(c("IL6", "TNF", "HAMP", "Slc39a14", "TF", "Tfrc", "Slc11a2", "SLC40A1", "SLC39A8", "FTL"))
gse_list <- c("GSE44770","GSE44771","GSE44768","GSE48350")
#gse_list <- c("GSE48350")
#gse_list <- c("GSE44770","GSE44771","GSE44768")
GPLlist = read.table("GPL.txt", header = T, sep='\t', stringsAsFactors=F, row.names = 1)

for(i in 1:length(row.names(ADlist))){
  #if(i<17) next
  gse_id = ADlist[i,1]
  if(length(intersect(gse_id, gse_list)) < 1) next
  
  gpl_id = ADlist[i,2]
  field =  ADlist[i,3]
  case =  ADlist[i,4]
  ctl =  ADlist[i,5]
  
  cat("Running... i=",i,"GSE=",gse_id, "GPL=",gpl_id)
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

      case_idx <- grep(case,types, ignore.case=TRUE)
      ctl_idx <- grep(ctl,types, ignore.case=TRUE)
      
      rel_idx = c(case_idx, ctl_idx)
      
      allexp <- matrix(0,length(gene_list)+1,length(rel_idx)+1, dimnames = list(c("type",gene_list), c("Cor_w_IL6",samples[rel_idx])))
      allexp[1,] = c("-", rep(1,length(case_idx)), rep(0,length(ctl_idx)))
      
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
        exp = exprs(gset1)[as.character(rowid),rel_idx]
        if(is.matrix(exp)){
          exp = colMeans(exp) #take means for all transcripts
        }
        
        if(j == 1){
          cor_gene_exp = exp
        }
        
        cor = cor(exp, cor_gene_exp)
        
        allexp[gene_list[j],] = c(cor,exp)
        
        #match GSE samples
        #case_idx <- grep(case,types, ignore.case=TRUE)
        #ctl_idx <- grep(ctl,types, ignore.case=TRUE)
        #exp[case_idx],exp[ctl_idx]
        
      }
      
    fname = paste(gse_id,i,"csv",sep=".")  
    write.csv(allexp, fname)
      
    }
  }
}


