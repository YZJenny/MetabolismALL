### sFig1
rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(RColorBrewer)
library(clusterProfiler)
col <- colorRampPalette(brewer.pal(8,'Paired'))(8)
col <- col[-6] #删除深红色
#setwd('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/')
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

#### sFig1D. markergene的featureplot (Related Fig1A)####
pbmc <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/pbmc.RDS')
## ALL:'CD7'
## Tcell:'CD3E','CD3D'
## B:'CD79A','CD19', 'CD79B','MS4A1'
## Mono:'CST3','LYZ'
## NK:'GNLY','NKG7', 'FGFBP2', 'PRF1'
## Ery: HBA1
marker.genes <- rev(c('CD7','CD3E','CD3D','CD79A','CD19','CST3','LYZ','GNLY','NKG7', 'FGFBP2', 'PRF1', 'HBA1' ))
df.gene <- data.frame(stringsAsFactors = F)
for (gene in marker.genes) {
  df.sub <- data.frame(expvalue = pbmc@assays$RNA@data[gene,],
                       gene = rep(gene, ncol(pbmc@assays$RNA@data)),
                       celltype = pbmc$celltype_abbr)
  df.gene <- rbind(df.gene, df.sub)
}
df.plot <- df.gene
df.plot$gene <- factor(df.gene$gene, levels = marker.genes)
df.plot$celltype <- factor(df.gene$celltype, 
                           levels = c('ALL', 'T', 'B', 'Mono.','NK','Ery.'))

color.cell <- c("#A6CEE3" ,"#B2DF8A","#33A02C", "#FB9A99","#FDBF6F","#1F78B4")
plot.vln <- 
  ggplot(data = df.plot, aes(x = gene, y = expvalue, color = celltype, fill = celltype)) + 
  geom_violin(trim = T, scale = 'width') + 
  scale_color_manual(labels = c('ALL', 'T', 'B', 'Mono.','NK','Ery.'),values = color.cell) + 
  scale_fill_manual(labels = c('ALL', 'T', 'B', 'Mono.','NK','Ery.'),values = color.cell) + 
  facet_grid( ~ celltype) + 
  theme_classic() + coord_flip() +
  stat_summary(fun= mean, geom = "point",shape = 23, size = 2, color = "black") + 
  labs(x = 'Gene', y = 'Expression Level') + 
  theme(axis.text.y = element_text(size = 8, color = "black", face = 'italic'), 
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8,color = "black"), 
        strip.text.x = element_text(size = 8, color = "black"), legend.position = 'none')
ggsave("/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/sFig1_Vln_markergene.pdf",plot.vln,
       height = 8, width = 16, units = 'cm')


#### sFig1E. Ratio of celltype by samples (Related Fig1A)####
phylog_df <- pbmc@meta.data[,c('orig.ident',"celltype_abbr")]
phylog_df <- table(phylog_df$orig.ident,phylog_df[,"celltype_abbr"])
phylog_df <- data.frame(phylog_df)
colnames(phylog_df) <- c('SampleID','CellType','Freq')
phylog_df$CellType <- factor(phylog_df$CellType)
phylog_df$SampleID <- factor(phylog_df$SampleID,
                             levels = c(unique(phylog_df$SampleID)[c(10:12,1:6,13:20,7:9)]))

p1 <- ggplot(phylog_df,aes(x=SampleID,y=Freq,fill=CellType))+
  geom_col(position = "fill", width = 0.6)+
  #coord_flip()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(face = 'bold',size = 8),
        axis.text.x = element_text(face = 'bold',size = 8,hjust = 0.5),
        axis.line = element_line(size=1, colour = "black"))+
  labs(x='',y='Proportion of cells')+theme(legend.position="right")+
  scale_fill_manual(values = color.cell)
ggsave('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/sFig_celltyperatio_bysample.pdf',p1,width = 7,height = 5)


#### sFig1FC. Barplot of number of genes of 7 superpathway ####
pathway2gene <- read.gmt('/remote-home/yanzijun/CRU/TALL_M/data/GeneSet7.gmt')
plot_data <- as.data.frame(table(pathway2gene$term))
colnames(plot_data) <- c('Superpathway','Num')
plot_data <- plot_data[order(plot_data$Num,decreasing = T),]
plot_data$Superpathway <- factor(plot_data$Superpathway,levels = rev(plot_data$Superpathway))
p <-ggplot(data=plot_data, aes(x=Superpathway, y=Num,fill=Superpathway))+
  geom_bar(stat="identity")+coord_flip()+
  theme_classic()+scale_fill_manual(values = rep('grey',7))
ggsave('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/sFig_barplot_superpathway.pdf',p,width = 7,height = 5)



#### sFig1F-G. step6.3 enrichment of each celltype in BMMC and ALL (related to Fig1B)####


#### sFig1. Heatmap of metab genes in ALL vs BMMC (Related Fig1D-F) ####
library(pheatmap)

pbmc_metab <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/pbmc_metab.RDS')
pathway2gene <- read.gmt('/remote-home/yanzijun/CRU/TALL_M/data/GeneSet7.gmt')
genelst <- unique(pathway2gene$gene)
print(length(genelst))

Idents(pbmc_metab) <- 'type'
tmp <- subset(pbmc_metab,downsample=500)
markers <- as.data.frame(genelst)
markerdata <- ScaleData(tmp, features = as.character(unique(markers$genelst)), assay = "RNA")

anno <- data.frame(type=markerdata$type)
rownames(anno) <- colnames(markerdata)
anno <- anno[order(anno$type),,drop=FALSE]

exp <- markerdata@assays$RNA@scale.data
exp <- as.data.frame(exp)
exp <- dplyr::select(exp,rownames(anno))

## 注释颜色
anno_col <- list(
  type=c(ALL='#E88482',BMMC='#63A8D2'))

bk <- c(seq(-3,-0.1,by=0.01),seq(0,3,by=0.01))
p <- pheatmap::pheatmap(
  exp,
  color=c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
  scale = "none",
  gaps_row = c(1,2),gaps_col = c(length(which(anno$type=='BMMC'))),
  annotation_col = anno,annotation_colors = anno_col,
  show_rownames = F,show_colnames = F,
  cluster_rows = F,cluster_cols = F,
  legend_breaks=seq(-2,2,1),breaks=bk
)
ggsave('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/sFig_Heatmap.pdf',p,width = 4,height = 5)

