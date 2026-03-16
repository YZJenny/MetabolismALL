#### Fig2D: Boxplot of M5-13 from scEA between sc TCA subtype ####
rm(list=ls())
library(ggplot2)
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/scFEA/')

scale_rows<- function(x){ #row:gene, col:sample
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m)/s)
}

scFEA <- read.csv('ALL_sc_sample5k_module168.csv',row.names = 1)
scFEA <- t(scale_rows(t(scFEA)))
M2C <- read.csv("/remote-home/yanzijun/software/scFEA-master/data/Human_M168_information.symbols.csv",row.names = 1)
M2C$name <- paste(M2C$Compound_IN_name,M2C$Compound_OUT_name,sep=' > ')

map <- M2C$name
names(map) <- rownames(map)


subtype <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/TCAsubtype_new/subtype_info.RDS')

ol <- intersect(rownames(scFEA),rownames(subtype))
sub_scFEA <- scFEA[rownames(scFEA) %in% ol,]
sub_subtype <- subtype[match(rownames(sub_scFEA),rownames(subtype)),]

print(all(rownames(sub_scFEA)==rownames(sub_subtype)))

plot.df <- cbind(sub_subtype,sub_scFEA[,paste('M_',5:13,sep='')])
plot.df$`TCA cycle`=NULL
plot.df <- reshape2::melt(plot.df)
plot.df$subtype <- factor(plot.df$subtype,levels = c('Low','Neutral','High'))
plot.df$name <- map[plot.df$variable]
plot.df$variable <- factor(plot.df$variable,levels = paste('M_',5:13,sep=''))

p <- ggplot(plot.df, aes(x=variable, y=value,fill=subtype)) +
  geom_boxplot(width = 0.5,outlier.shape = NA)+
  theme_classic()+labs(y='Metabolic flux',x='')+
  theme(plot.title = element_text(hjust = 0.5,size = 12),
        axis.text = element_text(size=12,face="plain",color="black"),
        axis.title  = element_text(size=12,face="plain",color="black"),
        axis.text.x = element_text(hjust = 1,angle =90),
        #legend.position="none"
  )+
  scale_fill_manual(values = c('#63A8D2','#DADADA','#E88482'))+
  ylim(c(0,6))+
  stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", method.args = list(alternative = "two.sided"))

ggsave('Fig2D_boxplot_TCAflux_subtype.pdf',p, width = 7,height = 3)


#### Fig2E. pheatmap+ KEGG enrichment ####
rm(list=ls())
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(future)
library(pheatmap)
plan("multiprocess", workers = 20)
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/TCAsubtype_new/')

pbmc <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/pbmc.RDS')
subtype <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/TCAsubtype_new/subtype_info.RDS')

pbmc_tumor <- subset(pbmc,cells=rownames(subtype))
print(all(colnames(pbmc_tumor)==rownames(subtype)))
pbmc_tumor$subtype <- factor(subtype$subtype,levels = c('High','Medium','Low'))

##### 1. Heatmap #####
mk <- readRDS('mk.RDS')
markerdata <- ScaleData(pbmc_tumor, features = as.character(unique(mk$markers)), assay = "RNA")
Idents(markerdata) <- markerdata$subtype

top.markers <- mk %>% group_by(cluster) %>% top_n(n=50,wt=avg_log2FC)

p1 <- DoHeatmap(markerdata,features = mk$gene,assay = 'RNA',slot = 'scale.data')+NoLegend()
ggsave('Heatmap_cells.pdf',p1,width = 10,height = 10)

tmp <- AverageExpression(markerdata, return.seurat = TRUE)
plot_data <- tmp@assays$RNA@scale.data[mk$gene,]
p2 <- pheatmap(plot_data,scale = 'none',cluster_cols = F,cluster_rows = F, 
               show_colnames = T,show_rownames = F, treeheight_row=0,treeheight_col=0,
               colorRampPalette(c("navy", "white", "firebrick3"))(length(seq(-1.2,1.2,by = 0.1))),
               breaks = seq(-1.2,1.2,by = 0.1),legend_breaks = seq(-1,1,1))
ggsave('Heatmap_subtypes.pdf',p2,width = 2,height = 7)


##### 2. KEGG/Hallmarks enrichment #####
###### 2.1 Enrichment ######
rm(list=ls())
library(clusterProfiler)
library(org.Hs.eg.db)
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/TCAsubtype_new/')

gmt_KEGG <- read.gmt('/remote-home/yanzijun/publicData/GMT/KEGG_clusterprofile.gmt')
gmt_HM <- read.gmt('/remote-home/yanzijun/publicData/GMT/h.all.v6.2.symbols.gmt')
gmt_GO <- read.gmt('/remote-home/yanzijun/publicData/GMT/c5.all.v6.2.symbols.gmt')

mk <- readRDS('mk.RDS')
mk <- na.omit(mk)
#mk <- mk[mk$p_val_adj<0.05,]

type='High'
gene_lst <- mk$gene[mk$cluster==type]
print(length(gene_lst))

gs <-bitr(gene_lst, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db",drop = TRUE)
markers1<-cbind(mk[gs[,1],],gs)
markers1 <- markers1[!duplicated(markers1$SYMBOL),]
geneList = markers1$avg_log2FC
names(geneList) = markers1$SYMBOL
geneList = sort(geneList,decreasing = T)


term='GO'
if(term=='KEGG'){
  gmt <- gmt_KEGG
}else if(term=='HM'){
  gmt <- gmt_HM
}else if(term=='GO'){
  gmt <- gmt_GO
}

egmt <- GSEA(geneList, TERM2GENE=gmt,verbose=F,maxGSSize=1000,pvalueCutoff = 0.05)
y=data.frame(egmt)
print(y[,1:7])
saveRDS(egmt,paste('GSEA/',type,'_',term,'.RDS',sep=''))

## enrichment: enrichKEGG需要联网，本地跑
library(enrichplot)
# GO
ego <- enrichGO(gs$ENTREZID,keyType = 'ENTREZID',OrgDb="org.Hs.eg.db",ont='BP',
                pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = TRUE)

y=data.frame(ego)

# KEGG需要联网，本地跑
ekegg <- enrichKEGG(gs$ENTREZID,organism = 'hsa',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
y=data.frame(ekegg)
saveRDS(ekegg,paste('Enrich/',type,'_KEGG.RDS',sep=''))

show_H <- c(1,5,9,12,13,15,18:20,23,26,27,31,32)
show_y <- y[show_H,]
saveRDS(show_y,paste('Enrich/',type,'_KEGG_pick.RDS',sep=''))

show_N <- c(1,6,12,15,18,28)
show_y <- y[show_N,]
saveRDS(show_y,paste('Enrich/',type,'_KEGG_pick.RDS',sep=''))

show_y <- y
saveRDS(show_y,paste('Enrich/',type,'_KEGG_pick.RDS',sep=''))

# 查看ID对应的Symbol
cc <- unlist(strsplit(c,'/'))
print(gs$SYMBOL[gs$ENTREZID %in% cc]) 

####### 2.2 plot #######
### 可视化
rm(list=ls())
get_barplot <- function(df,col){
  df$logfdr <- -log(df$p.adjust,10)
  df$Description=factor(df$Description,levels = rev(df$Description))
  p <- ggplot(df,aes(Description,logfdr))+
    #coord_cartesian(expand = F)+
    geom_bar(stat='identity',width = 0.8,colour=col,fill=col)+
    coord_flip() + #horizontal
    labs(y='-log10(padj)',x="")+
    theme(panel.grid=element_blank(),
          panel.background = element_blank(),
          axis.ticks.y= element_blank())+ #delete backgroud
    guides(fill=FALSE)+ #delete legend
    theme(axis.text.x = element_text(size = 8,colour = 'black'),
          axis.text.y = element_text(size = 8,colour = 'black'), 
          axis.title.x = element_text(size = 8,colour = 'black'), 
          axis.title.y = element_text(size = 6,colour = 'black'))+
    ylim(c(0,19))
  return(p)
}

mycol <- c('#E88482','grey','#63A8D2')
names(mycol) <- c('H-A','M-A','L-A')

df_H <- readRDS('Enrich/High_KEGG_pick.RDS')
df_N <- readRDS('Enrich/Middle_KEGG_pick.RDS')
df_L <- readRDS('Enrich/Low_KEGG_pick.RDS')

plot_H <- get_barplot(df=df_H,col = mycol[1])
ggsave('Enrich/High_KEGG_pick.pdf',plot_H,width = 5,height = 3)

plot_N <- get_barplot(df=df_N,col = mycol[2])
ggsave('Enrich/Middle_KEGG_pick.pdf',plot_N,width = 5,height = 1.8)

plot_L <- get_barplot(df=df_L,col = mycol[3])
ggsave('Enrich/Low_KEGG_pick.pdf',plot_L,width = 5,height = 1.1)
