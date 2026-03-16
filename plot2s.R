### sFig2
rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(RColorBrewer)
library(clusterProfiler)
col1 <-colorRampPalette(brewer.pal(8,'Paired'))(17)

setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

pbmc_normal <- readRDS('pbmc_normal_metab.RDS')
meta <- readRDS('meta.RDS') #要用最新的meta文件
map <- meta$celltype_abbr
names(map) <- rownames(meta)

#### sFig2A: 代谢基因的TSNE of normal cells in ALL (Related Fig2A) ####
pbmc_normal$celltype_abbr <- map[colnames(pbmc_normal)]
CT <- c('B','T','Mono.','NK','Ery.')
print(CT)
Color <-c("#33A02C", "#B2DF8A","#FB9A99","#FDBF6F","#1F78B4")
names(Color) <- CT

cells=Cells(pbmc_normal)[-grep('PBMMC',pbmc_normal$orig.ident)] # normal cells in ALL
tmp <- subset(pbmc_normal,cells=cells)
tmp <- RunTSNE(tmp, dims = 1:30,check_duplicates = FALSE,seed.use = 1)

pdf('sFig_TSNE_normal_metab.pdf',width = 6.2,height = 5)
TSNEPlot(pbmc_normal, group.by = "celltype_abbr",cols = Color)
TSNEPlot(tmp, group.by = "celltype_abbr",cols = Color)
dev.off()


#### sFig2B. variance of PC: step7_intra_ALL_heter.R (Related to Fig2B)#####
rm(list=ls())
library(ggplot2)
library(clusterProfiler)
library(scater)
library(stringr)
library(scran)
library(gtools)
options(stringsAsFactors=FALSE)
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

type='ALL'
method='CV'
plot_df <- readRDS(paste('../CV_SD/enriched_pathway','_',type,'_',method,'.RDS',sep=''))
pvals <- plot_df$PVAL
plot_df$PVAL <- -log10(pvals)

#sort
pathway_pv_sum <- by(plot_df$PVAL,plot_df$y,FUN=sum)
pathway_order <- names(pathway_pv_sum)[order(pathway_pv_sum,decreasing = T)]

plot_df <- plot_df[plot_df$y %in% pathway_order,]
samples <-  unique(plot_df$x)
plot_df$x <- factor(plot_df$x, levels = samples)
plot_df$y <- factor(plot_df$y,levels = pathway_order)

##buble plot
p <- ggplot(plot_df, aes(x = x, y = y, size = PVAL, color = NES)) +
  geom_point(shape=19) +
  #ggtitle("pathway heterogeneity") +
  labs(x = NULL, y = NULL,
       size = "-log10 adjp", color = "NES") +
  scale_size(range = c(0, 2.5)) +
  scale_color_gradient( low = "white", high = "red") +
  #scale_color_gradient2(low="red",mid="white",high="blue",midpoint = 1) +
  theme(legend.position = "bottom", legend.direction = "vertical",
        legend.box = "horizontal",
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(colour="black",size=8),
        legend.title =  element_text(colour="black",size=8),
        axis.line = element_line(size=0.3, colour = "black"),
        #panel.grid.major = element_line(colour = "#d3d3d3"),
        #panel.grid.minor = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.3),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 8,angle=90,hjust=1,vjust=0.5),
        axis.text.y=element_text(colour="black", size = 8)) +
  theme(plot.margin = unit(rep(1,4),"lines"))
ggsave(paste("../CV_SD/sFig_enriched_pathway_",type,'_',method,".pdf",sep=''),p,
       width = 5.5,height=3.3,units="in",device="pdf",useDingbats=FALSE)



#### sFig2C. dotplot of SD: step7_intra_ALL_heter_CVSD.R(Related to Fig2B)####


#### sFig2D. PCA of sc TCA subtype (Related to Fig2D)####
rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(sva)
library(FactoMineR)
library(factoextra)
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

pbmc_tumor <- readRDS('pbmc_tumor_metab.RDS')
meta <- readRDS('meta.RDS') #用最新meta
map <- meta$celltype_abbr
names(map) <- rownames(meta)
pbmc_tumor$celltype_abbr <- map[colnames(pbmc_tumor)]

new_meta <- readRDS('../TCAsubtype_new/subtype_info.RDS') ## TCA subtype info from step10_TCAsubtype_v2.R
print(table(new_meta$subtype))

print(all(colnames(pbmc_tumor)==rownames(new_meta)))
pbmc_tumor$subtype <- new_meta$subtype


pca.plot = function(dat,col,pal){
  df.pca <- PCA(t(dat), graph = FALSE)
  fviz_pca_ind(df.pca,
               geom.ind = "point",
               col.ind = col ,
               addEllipses = FALSE,
               legend.title = "Groups",
               palette = pal
  )
}

gmt <- read.gmt('/remote-home/yanzijun/CRU/TALL_M/data/GeneSet7.gmt')
geneset <- gmt$gene[gmt$term=='TCA cycle']

## 整理sc数据
sc_EXP <- pbmc_tumor@assays$SCT@data
sc_EXP <- sc_EXP[rownames(sc_EXP) %in% geneset,]

sc_Subtype <- as.data.frame(pbmc_tumor$subtype)
colnames(sc_Subtype) <- 'subtype'
print(all(colnames(sc_EXP)==sc_Subtype$sampleID))

##保存数据
saveRDS(sc_EXP, '../TCAsubtype_new/exp_tumor_TCAcycle.RDS')
saveRDS(sc_Subtype, '../TCAsubtype_new/meta_tumor_TCAcycle.RDS')

p1 <- pca.plot(sc_EXP,factor(sc_Subtype$subtype),pal=c('#E88482','#63A8D2','#DADADA'))  
pdf('sFig2_PCA_sc_subtype.pdf',width = 5,height = 4)
print(p1)
dev.off()


#### sFig2E. barplot of sc TCA subtype in each sample (Related Fig2D)####
phylog_df <- pbmc_tumor@meta.data[,c('orig.ident',"subtype")]
phylog_df <- table(phylog_df$orig.ident,phylog_df[,"subtype"])
phylog_df <- data.frame(phylog_df)
colnames(phylog_df) <- c('SampleID','subtype','Freq')
phylog_df$subtype <- factor(phylog_df$subtype,levels = c('Low','Neutral','High'))
phylog_df$SampleID <- factor(phylog_df$SampleID)

p1 <- ggplot(phylog_df,aes(x=SampleID,y=Freq,fill=subtype))+
  geom_col(position = "fill", width = 0.6)+
  #coord_flip()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(face = 'bold',size = 8),
        axis.text.x = element_text(face = 'bold',size = 8,hjust = 0.5),
        axis.line = element_line(size=1, colour = "black"))+
  labs(x='',y='Proportion of cells')+theme(legend.position="right")+
  scale_fill_manual(values = c('#63A8D2','#DADADA','#E88482'))
ggsave('/remote-home/yanzijun/CRU/ped_M/res_ALL/TCAsubtype_new/sFig_subtyperatio_bysample.pdf',p1,width = 9,height = 5)




#### sFig 计算metaGene ssgsea score与superpathway ssgsea score的correlation ####
scale_rows<- function(x){ #行为基因，列为样本
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m)/s)
}

score_irGSEA <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/irGSEA/pbmc_all_superpathway.RDS')
ssgsea <- score_irGSEA$assays$ssgsea
meta <- score_irGSEA$meta

ssgsea_tumor <- scale_rows(as.matrix(ssgsea[,meta$celltype_2=='leukemia cells']))
cor.res <- round(cor(t(ssgsea_tumor),method = 'spearman'), 3)

plot.cor <- pheatmap(cor.res,scale = 'none',na_col="white",cluster_rows = F,cluster_cols = F,
                     colorRampPalette(c("navy", "white", "firebrick3"))(length(seq(0,1,by = 0.01))),
                     breaks = seq(0,1,by = 0.01),legend_breaks = seq(0,1,0.2))
ggsave('/remote-home/yanzijun/CRU/ped_M/res_ALL/irGSEA/corplot_superpathway.pdf',plot.cor,width = 5,height = 4)



