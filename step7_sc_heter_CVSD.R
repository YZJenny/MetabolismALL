### CV_SD+GSEA:tumor cells by batch + normal cells in ALL
rm(list=ls())
library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(dplyr)
options(stringsAsFactors=FALSE)
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

pathway2gene <- read.gmt('/remote-home/yanzijun/CRU/TALL_M/data/GeneSet7.gmt')
genelst <- unique(pathway2gene$gene)
print(length(genelst))

####1. Loading the data ####
pbmc <- readRDS('pbmc_metab.RDS')
meta <- readRDS('meta.RDS')
print(all(colnames(pbmc)==rownames(meta)))
pbmc$celltype_abbr <- meta$celltype_abbr


type='ALL' 
if(type=='all'){
  pbmc_metab <- pbmc[genelst,]
}else if(type=='ALL'){
  pbmc_metab <- pbmc[genelst,Cells(pbmc)[-grep('^PBMMC',pbmc$orig.ident)]]
}

print(table(pbmc_metab$orig.ident))


df <- pbmc_metab@meta.data[,c('orig.ident','celltype_abbr')]
df <- df %>% mutate(celltype_abbr=if_else(celltype_abbr=='ALL',orig.ident,celltype_abbr))
print(all(rownames(df)==colnames(pbmc_metab)))
pbmc_metab$celltype_3 <- df$celltype_abbr

samples <- c(sort(unique(pbmc_metab$orig.ident)),
             sort(unique(as.character(pbmc_metab$celltype_abbr))[-1]))
print(samples)


#### 2.CV and SD ####
get_gsea <- function(s,type,pbmc_metab){
  print(s)
  if(s=='allsample'){ 
    each_metabolic_sce <- pbmc_metab[,pbmc_metab$celltype_abbr=='ALL']
  }else{
    each_metabolic_sce <- pbmc_metab[,pbmc_metab$celltype_3==s]
  }
  
  each_metabolic_tpm <- each_metabolic_sce@assays$RNA@data
  each_metabolic_tpm <- each_metabolic_tpm[rowSums(each_metabolic_tpm)>0,]
  x <- each_metabolic_tpm
  
  # calculate SD
  gene_sd <- apply(x, 1, sd)
  gene_sd <- sort(gene_sd, decreasing = T)
  print(head(gene_sd))
  
  # calculate CV
  gene_mean <- apply(x, 1, mean)
  gene_cv <- gene_sd / gene_mean
  gene_cv <- sort(gene_cv, decreasing = T)
  print(head(gene_cv))
  
  #### GSEA of SD
  gsea_sd <- GSEA(geneList = gene_sd, TERM2GENE=pathway2gene, verbose=FALSE,maxGSSize = 2000,minGSSize = 10,
                  pvalueCutoff = 1,pAdjustMethod='fdr',seed = FALSE,nPermSimple=1000) #有些通路为NA
  gsea_sd <- as.data.frame(gsea_sd)
  
  #### GSEA of CV
  gsea_cv <- GSEA(geneList = gene_cv, TERM2GENE=pathway2gene, verbose=FALSE,maxGSSize = 2000,minGSSize = 10,
                  pvalueCutoff = 1,pAdjustMethod='fdr',seed = FALSE,nPermSimple=1000) #有些通路为NA
  gsea_cv <- as.data.frame(gsea_cv)
  
  gsea <- list(gsea_sd=gsea_sd,gsea_cv=gsea_cv)
}

##### 2.1 inter-heterogeneity #####
enrich_sd <- data.frame(x=NULL,y=NULL,NES=NULL,PVAL=NULL)
enrich_cv <- data.frame(x=NULL,y=NULL,NES=NULL,PVAL=NULL)
for(s in samples){
  gsea  <- get_gsea(s,type,pbmc_metab)
  
  gsea_sd <- gsea$gsea_sd
  saveRDS(gsea_sd,paste('../CV_SD/',s,'_',type,'_SD.RDS',sep=''))
  
  # get the result
  gsea_pathways <- rownames(gsea_sd)
  enrich_sd <- rbind(enrich_sd,
                     data.frame(x=s,y=gsea_pathways,NES=gsea_sd$NES,PVAL=gsea_sd$p.adjust))
  
  gsea_cv <- gsea$gsea_sd
  saveRDS(gsea_cv,paste('../CV_SD/',s,'_',type,'_CV.RDS',sep=''))
  
  # get the result
  gsea_pathways <- rownames(gsea_cv)
  enrich_cv <- rbind(enrich_cv,
                     data.frame(x=s,y=gsea_pathways,NES=gsea_cv$NES,PVAL=gsea_cv$p.adjust))
  
}
saveRDS(enrich_sd,paste('../CV_SD/enriched_pathway','_',type,'_SD.RDS',sep=''))
saveRDS(enrich_cv,paste('../CV_SD/enriched_pathway','_',type,'_CV.RDS',sep=''))

##### 2.2 intra-heterogeneity #####
gsea  <- get_gsea('allsample',type,pbmc_metab)
gsea_sd <- gsea$gsea_sd
saveRDS(gsea_sd,paste('../CV_SD/allsample_',type,'_SD.RDS',sep=''))

gsea_cv <- gsea$gsea_cv
saveRDS(gsea_sd,paste('../CV_SD/allsample_',type,'_CV.RDS',sep=''))