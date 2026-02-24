#### 用上下1/4分位点代替v1版本的median划分
#### 1. 根据TCA cycle的activity划分tumor cells into high and low subtype ###
rm(list=ls())
library(Seurat)
library(future)
#setwd('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/')
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

plan("multiprocess", workers = 20)

scale_rows<- function(x){ #行为基因，列为样本
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
new_meta$subtype[new_meta$`TCA cycle`< quantile(new_meta$`TCA cycle`,0.25)] <- 'Low'
new_meta$subtype[new_meta$`TCA cycle`> quantile(new_meta$`TCA cycle`,0.75)] <- 'High'
saveRDS(new_meta,'../TCAsubtype_new/subtype_info.RDS')


## findmk
pbmc <- readRDS('pbmc.RDS')
sub_meta <- new_meta[new_meta$subtype != 'Neutral',]
cells <- rownames(sub_meta)

pbmc_tumor <- subset(pbmc,cells = cells)
print(all(colnames(pbmc_tumor)==rownames(sub_meta)))

pbmc_tumor@meta.data$subtype <- sub_meta$subtype 

Idents(pbmc_tumor) <- 'subtype'
mk <- FindMarkers(pbmc_tumor,ident.1 = 'High',ident.2 = 'Low',only.pos = FALSE,logfc.threshold = 0,
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


## Fig2F: 可视化
rm(list=ls())
library(ggplot2)
#setwd('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/')
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

term='KEGG'
egmt <- readRDS(paste('../TCAsubtype_new/GSEA_',term,'.RDS',sep=''))
y=data.frame(egmt)

plot_df <- y[order(y$NES,decreasing = T),]

if(term=='KEGG'){
  # plot_df <- plot_df[c(1,2,9,12,14:20,22,24,25,27,34,37,40,48,50,51,53,55,61,65:67,76,82),] #GSE_nc
  plot_df <- plot_df[c(2,12,14:19,24,25,27,32,35,40:41,46,51,54,57:58,69:71),] 
}

plot_df$Description <- factor(plot_df$Description,levels = rev(plot_df$Description ))

library(RColorBrewer)
cols <- rev(colorRampPalette(c("grey","#E88482"))(10))
p <-ggplot(data=plot_df, aes(x=Description , y=NES ,fill=p.adjust))+
  geom_bar(stat="identity",width = 0.7)+coord_flip()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(color='black',size = 8),
        axis.title  = element_text(color='black',size = 8))+
  scale_fill_gradientn(colours = cols)

ggsave(paste("../TCAsubtype/Fig2F_GSEAbarplot_",term,"_red.pdf",sep=''),
       width = 5.5,height = 4.5)


cols <- rev(colorRampPalette(c("grey","#63A8D2"))(10))
p <-ggplot(data=plot_df, aes(x=Description , y=NES ,fill=p.adjust))+
  geom_bar(stat="identity",width = 0.7)+coord_flip()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(color='black',size = 8),
        axis.title  = element_text(color='black',size = 8))+
  scale_fill_gradientn(colours = cols)

ggsave(paste("../TCAsubtype/Fig2F_GSEAbarplot_",term,"_blue.pdf",sep=''),width = 5.5,height = 4.5)

