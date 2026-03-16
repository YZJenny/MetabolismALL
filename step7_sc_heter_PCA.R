### PCA+GSEA: tumor cells by batch + normal cells in ALL
rm(list=ls())
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
library(scater)
library(stringr)
library(scran)
library(gtools)
library(dplyr)
options(stringsAsFactors=FALSE)
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

pathway2gene <- read.gmt('/remote-home/yanzijun/CRU/TALL_M/data/GeneSet7.gmt')
genelst <- unique(pathway2gene$gene)
print(length(genelst))

#### 1. Loading the data ####
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

#### 2.PCA ####
get_gsea <- function(s,type,pbmc_metab){
  print(s)
  if(s=='allsample'){ ## heterogeneity of tumor cell
    each_metabolic_sce <- pbmc_metab[,pbmc_metab$celltype_abbr=='ALL']
  }else{
    each_metabolic_sce <- pbmc_metab[,pbmc_metab$celltype_3==s]
  }
  
  each_metabolic_tpm <- each_metabolic_sce@assays$RNA@data
  each_metabolic_tpm <- each_metabolic_tpm[rowSums(each_metabolic_tpm)>0,]
  x <- each_metabolic_tpm
  ntop <- nrow(x)
  rv <- rowVars(x)
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(x[select,]))
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  ### select PCs that explain at least 80% of the variance
  cum_var <- cumsum(percentVar)
  select_pcs <- which(cum_var>0.8)[1]
  
  ### plot the PCA and explained variances
  tmp_plotdata <- data.frame(x=1:length(percentVar),y=percentVar,
                             sel=c(rep("y",select_pcs),rep("n",length(percentVar)-select_pcs)),
                             celltype=rep(s,length(percentVar)))
  
  ### GSEA
  pre_rank_matrix <- as.matrix(rowSums(abs(pca$rotation[,1:select_pcs])))
  print(dim(pre_rank_matrix))
  
  genes <- as.numeric(pre_rank_matrix[,1])
  names(genes) <- rownames(pre_rank_matrix)
  genes <- sort(genes, decreasing = T)
  
  gsea <- GSEA(geneList = genes, TERM2GENE=pathway2gene, verbose=FALSE,maxGSSize = 2000,minGSSize = 10,
               pvalueCutoff = 1,pAdjustMethod='fdr',seed = FALSE,nPermSimple=1000) #有些通路为NA
  res <- list(gsea=gsea,tmp_plotdata=tmp_plotdata)
  return(res)
}

##### 2.1 inter-heterogeneity ##### 
pc_plotdata <- data.frame(x=numeric(),y=numeric(),sel=character(),celltype=character())
enrich_data_df <- data.frame(x=NULL,y=NULL,NES=NULL,PVAL=NULL)

for(s in samples){
  res  <- get_gsea(s,type,pbmc_metab)
  gsea_result <- as.data.frame(res$gsea)
  saveRDS(gsea_result,paste('../PC/',s,'_',type,'.RDS',sep=''))
  
  gsea_pathways <- rownames(gsea_result)
  enrich_data_df <- rbind(enrich_data_df,
                          data.frame(x=s,y=gsea_pathways,NES=gsea_result$NES,PVAL=gsea_result$p.adjust))
  
  tmp_plotdata <- res$tmp_plotdata
  pc_plotdata <- rbind(pc_plotdata,tmp_plotdata)
}
saveRDS(enrich_data_df,paste('../PC/enriched_pathway_',type,'.RDS',sep=''))
saveRDS(pc_plotdata,paste('../PC/pc_plotdata_',type,'.RDS',sep=''))

##### 2.2 intra-heterogeneity ##### 
res <- get_gsea('allsample',type,pbmc_metab)
gsea_result <- res$gsea
saveRDS(gsea_result,paste('../PC/allsample_',type,'.RDS',sep=''))

pc_plotdata <- res$tmp_plotdata
saveRDS(pc_plotdata,paste('../PC/allsample_',type,'_pc_plot.RDS',sep=''))