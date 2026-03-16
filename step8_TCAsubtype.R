## classified tumor cells accroding the quantail of activity score
#### 1. classified tumor cells into three states####
rm(list=ls())
library(Seurat)
library(future)
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

plan("multiprocess", workers = 20)

scale_rows<- function(x){ #row：gene, col: sample
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m)/s)
}

score_irGSEA <- readRDS('../irGSEA/pbmc_all_superpathway.RDS')
ssgsea <- score_irGSEA$assays$ssgsea
meta <- score_irGSEA$meta

ssgsea_tumor <- scale_rows(as.matrix(ssgsea[,meta$celltype_2=='leukemia cells']))

## 1/4 as cutoff
new_meta <- as.data.frame(t(ssgsea_tumor[rownames(ssgsea_tumor)=='TCA cycle',,drop=FALSE]))
new_meta$subtype <- 'Neutral'
new_meta$subtype[new_meta$`TCA cycle`< quantile(new_meta$`TCA cycle`,0.25)] <- 'L-A'
new_meta$subtype[new_meta$`TCA cycle`> quantile(new_meta$`TCA cycle`,0.75)] <- 'H-A'
saveRDS(new_meta,'../TCAsubtype_new/subtype_info.RDS')


#### 2. findmk of three states ####
pbmc <- readRDS('pbmc.RDS')
sub_meta <- new_meta[new_meta$subtype != 'M-A',]
cells <- rownames(sub_meta)

pbmc_tumor <- subset(pbmc,cells = cells)
print(all(colnames(pbmc_tumor)==rownames(sub_meta)))

pbmc_tumor@meta.data$subtype <- sub_meta$subtype 

Idents(pbmc_tumor) <- 'subtype'
mk <- FindMarkers(pbmc_tumor,ident.1 = 'H-A',ident.2 = 'L-A',only.pos = FALSE,logfc.threshold = 0,
                  min.pct = 0,min.cells.feature=0,min.cells.group =0)
saveRDS(mk,'../TCAsubtype_new/mk_HighvsLow_2.RDS')
saveRDS(pbmc_tumor,'../TCAsubtype_new/pbmc_tumor.RDS')


## KEGG/Hallmarks enrichment
rm(list=ls())
library(clusterProfiler)
library(org.Hs.eg.db)
#setwd('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/')
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

mk <- readRDS('../TCAsubtype_new/mk_HighvsLow.RDS')
mk <- na.omit(mk)
#mk <- mk[mk$p_val_adj<0.05,]

gs <-bitr(rownames(mk), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db",drop = TRUE)
markers1<-cbind(mk[gs[,1],],gs)
markers1 <- markers1[!duplicated(markers1$SYMBOL),]
geneList = markers1$avg_log2FC
names(geneList) = markers1$SYMBOL
geneList = sort(geneList,decreasing = T)


gmt_KEGG <- read.gmt('/remote-home/yanzijun/publicData/GMT/KEGG_clusterprofile.gmt')
gmt_HM <- read.gmt('/remote-home/yanzijun/publicData/GMT/h.all.v6.2.symbols.gmt')

term='KEGG'
if(term=='KEGG'){
  gmt <- gmt_KEGG
}else if(term=='HM'){
  gmt <- gmt_HM
}

egmt <- GSEA(geneList, TERM2GENE=gmt,verbose=F,maxGSSize=1000,pvalueCutoff = 0.05)
y=data.frame(egmt)
print(y[,1:7])

saveRDS(egmt,paste('../TCAsubtype_new/GSEA_',term,'.RDS',sep=''))