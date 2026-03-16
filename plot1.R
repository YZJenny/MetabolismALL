### Fig1 
rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(RColorBrewer)
library(clusterProfiler)
col <- colorRampPalette(brewer.pal(8,'Paired'))(8)
col <- col[-6] #删除深红色
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')


### Fig1A
pbmc <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/pbmc.RDS')

CT <- c('ALL','Ery.','T','B','Mono.','NK')
print(CT)
Color <- col[1:length(CT)]
names(Color) <- CT


#### Fig1A by celltype ####
p <- DimPlot(pbmc, group.by='celltype_abbr',label=F,pt.size = 0.01,
             cols = Color)
ggsave("/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/Fig1A_UMAP_CellType.pdf",p,width = 7,height = 6)

p <- DimPlot(pbmc, group.by='celltype_abbr',label=F,reduction = 'tsne',pt.size = 0.01,
             cols = Color)
ggsave("/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/Fig1A_TSNE_CellType.pdf",p,width = 7,height = 6)


#### Fig1A by study/type/sample/tisuue ####
## study
p <- DimPlot(pbmc, group.by='study',label=F,pt.size = 0.01,
             cols = c('#A6CEE3','#B2DF8A','#FB9A99'))
ggsave("/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/Fig1A_UMAP_study.pdf",p,width = 7,height = 6)

p <- DimPlot(pbmc, group.by='study',label=F,reduction = 'tsne',pt.size = 0.01,
             cols = c('#A6CEE3','#B2DF8A','#FB9A99'))
ggsave("/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/Fig1A_TSNE_study.pdf",p,width = 7,height = 6)

## type
p <- DimPlot(pbmc, group.by='type',label=F,pt.size = 0.01,
             cols = c('#A6CEE3',"#FB9A99"))
ggsave("/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/Fig1A_UMAP_type.pdf",p,width = 7,height = 6)

p <- DimPlot(pbmc, group.by='type',label=F,reduction = 'tsne',pt.size = 0.01,
             cols = c('#A6CEE3',"#FB9A99"))
ggsave("/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/Fig1A_TSNE_type.pdf",p,width = 7,height = 6)


## tissue
p <- DimPlot(pbmc, group.by='tissue',label=F,pt.size = 0.01,
             cols = c('#A6CEE3','#B2DF8A','#FB9A99'))
ggsave("/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/Fig1A_UMAP_tissue.pdf",p,width = 7,height = 6)

p <- DimPlot(pbmc, group.by='tissue',label=F,reduction = 'tsne',pt.size = 0.01,
             cols = c('#A6CEE3','#B2DF8A','#FB9A99'))
ggsave("/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/Fig1A_TSNE_tissue.pdf",p,width = 7,height = 6)

## orig.ident
col1 <- c(colorRampPalette(brewer.pal(8,'Paired'))(17),colorRampPalette(brewer.pal(8,'Set2'))(3))
p <- DimPlot(pbmc, group.by='orig.ident',label=F,pt.size = 0.01,
             cols = col1)
ggsave("/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/Fig1A_UMAP_sample.pdf",p,width = 7.5,height = 6)

p <- DimPlot(pbmc, group.by='orig.ident',label=F,reduction = 'tsne',pt.size = 0.01,
             cols = col1)
ggsave("/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/Fig1A_TSNE_sample.pdf",p,width = 7.5,height = 6)


#### Fig1A celltype number cycle plot ####
df<-as.data.frame(table(pbmc$celltype_abbr))
colnames(df) <- c('category','count')
df$fraction<-df$count/sum(df$count)
df$ymax<-cumsum(df$fraction)
df$ymin<-c(0,head(df$ymax,n=-1))

# df$labelPosition<-(df$ymax + df$ymin)/2
# df$label<-paste0(df$category,"\n cells: ",df$count)
p <- ggplot(df,aes(ymax=ymax,ymin=ymin,
              xmax=4,xmin=3))+
  geom_rect(aes(fill=category))+
  # geom_label(x=3.5,aes(y=labelPosition,label=label),size=4)+
  #scale_fill_brewer(palette = 4)+
  scale_fill_manual(values = Color)+
  coord_polar(theta = "y")+
  xlim(-6,4)+
  theme_void()+
  theme(legend.position = "none")
ggsave("/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/Fig1A_Pieplot_cellNum.pdf",p,width = 7,height = 7)
  


#### Fig1B violin plot of overall metab activity in ALL/BMMC的celltype ####
## Fig1B
rm(list=ls())
library(irGSEA)
library(Seurat)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(RColorBrewer)
library(ggpubr)
library(ggsignif)

assay <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/irGSEA/pbmc_all_superpathway.RDS')

score <- assay$meta
score$metabGene1 <- assay$assays$ssgsea@data['metabGene',]

score_tumor <- score[score$celltype=='leukemia cells',]
score_normal <- score[score$celltype !='leukemia cells',]

score_ALL <- score[-grep('^PBMMC',score$orig.ident),]
score_BMMC <- score[grep('^PBMMC',score$orig.ident),]


#### Fig1B (left pannel). ALL:boxplot of overall activity by celltype #####
data=score_ALL
CT=c('ALL cells','B cells',"Erythrocytes", "Monocytes" , "NK cells" , "pDCs" )
mycol = colorRampPalette(brewer.pal(8, "Set2"))(8)

plot.df <- data[,c('celltype_2','metabGene1')]
colnames(plot.df) <- c('celltype','Activity')
plot.df$celltype[plot.df$celltype=='leukemia cells'] <- 'ALL cells'

med <- aggregate(round(plot.df$logActivity,5),list(plot.df$celltype),median)
rank <- med$Group.1[order(med$x,decreasing = T)]

plot.df$celltype=factor(plot.df$celltype,levels = rank)
plot.df$logActivity <- log(plot.df$Activity,10)

my_comparisons <- list( c("ALL cells", "pDCs"), c("ALL cells", "Monocytes"), c("ALL cells", "B cells"),
                        c("ALL cells", "NK cells"), c("ALL cells", "Erythrocytes"))

p <- ggplot(plot.df, aes(x=celltype, y=logActivity,fill=celltype)) + 
  geom_violin(trim = FALSE,width = 1)+
  #geom_boxplot(width = 0.5,outlier.shape = NA)+
  #stat_boxplot(geom = "errorbar",width=0.2)+
  stat_summary(fun= median, geom = "point",shape = 23, size = 2, color = "black")+ #添加中位值
  theme_classic()+labs(y='Metabolic activity score',x='')+
  theme(plot.title = element_text(hjust = 0.5,size = 12),
        axis.text = element_text(size=12,face="plain",color="black"),
        axis.title  = element_text(size=12,face="plain",color="black"),
        axis.text.x = element_text(hjust = 1,angle =20),
        legend.position="none")+
  # geom_signif(comparisons = my_comparisons,
  #             test = "wilcox.test", test.args = c("great"),
  #             y_position=seq(from=max(plot.df$Activity),
  #                            to=max(plot.df$Activity)+4*0.04,
  #                            by=0.04),
  #             map_signif_level=TRUE)+
  scale_fill_manual(values =mycol[1:length(CT)])

pdf('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/irGSEA/Fig1_vlnplot_ssgsea_ALL.pdf',width = 5,height = 4)
print(p)
dev.off()


#### Fig1B (right pannel). BMMC:boxplot of overall activity by celltype #####
data=score_BMMC
CT=c('B cells',"Erythrocytes", "Monocytes" , "NK cells" )
mycol = c("#E78AC3","#FFD92F","#8DA0CB", "#A6D854")
names(mycol) <- CT

med <- aggregate(round(data$metabGene1,2),list(data$celltype_2),median)
rank <- med$Group.1[order(med$x,decreasing = T)]

plot.df <- data[,c('celltype_2','metabGene1')]
colnames(plot.df) <- c('celltype','Activity')
plot.df$celltype=factor(plot.df$celltype,levels = rank)

p <- ggplot(plot.df, aes(x=celltype, y=Activity,fill=celltype)) + 
  geom_violin(trim = FALSE,width = 1)+
  #geom_boxplot(width = 0.5,outlier.shape = NA)+
  #stat_boxplot(geom = "errorbar",width=0.2)+
  stat_summary(fun= median, geom = "point",shape = 23, size = 2, color = "black")+ #添加中位值
  theme_classic()+labs(y='Metabolic activity score',x='')+
  theme(plot.title = element_text(hjust = 0.5,size = 12),
        axis.text = element_text(size=12,face="plain",color="black"),
        axis.title  = element_text(size=12,face="plain",color="black"),
        axis.text.x = element_text(hjust = 1,angle = 20),
        legend.position="none")+
  stat_compare_means()+
  scale_fill_manual(values =mycol[rank])

pdf('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/irGSEA/Fig1_vlnplot_addModulescore_BMMC.pdf',
    width = 4,height = 4)
print(p)
dev.off()

#### Fig1C pheatmap of celltype between ALL vs BMMC ####
rm(list=ls())
library(ggplot2)
library(dplyr)
library(pheatmap)
#setwd('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/')
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

## 导入数据
irGSEA <- readRDS('../irGSEA/pbmc_all_superpathway.RDS')
meta <- readRDS('meta.RDS') #要用最新的meta文件
assays <- irGSEA$assays

get_plotdf <- function(assays,meta,method,metabGene='FALSE',type){
  df <- as.data.frame(t(as.matrix(assays[[method]]@counts)))
  print(all(rownames(df)==rownames(meta)))
  df$celltype <- meta$celltype_abbr
  df$celltype <- as.factor(df$celltype)
  
  cells <- rownames(meta)[meta$type==type]
  df <- df[rownames(df) %in% cells,]
  # 根据类型分组，计算每个基因在每个类别中的平均中位值
  plot.df <- aggregate(.~ celltype, data = df, median)
  rownames(plot.df) <- plot.df$celltype;plot.df$celltype=NULL
  plot.df <- t(plot.df) 
  if(metabGene=='FALSE'){
    plot.df <- plot.df[rownames(plot.df) !='metabGene',]
  }
  return(plot.df)
}



#### in ALL
bk <- c(seq(-1.5,-0.1,by=0.01),seq(0,1.5,by=0.01))
df1 <- get_plotdf(assays,meta,method='ssgsea',metabGene='FALSE',type='ALL')
p1 <- pheatmap(df1,scale = 'row',cluster_cols = T,cluster_rows = F,
               color=c(colorRampPalette(colors = c("#63A8D2","white"))(length(bk)/2),colorRampPalette(colors = c("white","#E88482"))(length(bk)/2)),
               legend_breaks=seq(-1.5,1.5,1),breaks=bk)

color <- colorRampPalette(c('grey',"white","red"))(100)
df2 <- get_plotdf(assays,meta,method='addMS',metabGene='FALSE',type='ALL')
p2 <- pheatmap(df2,scale = 'row',cluster_cols = T,cluster_rows = F,color=color)

p <- cowplot::plot_grid(p1$gtable, p2$gtable, ncol= 2)
pdf('../irGSEA/Fig1_Heatmap_ALL_superpathway.pdf', width = 8,height = 4)
print(p)
dev.off()


#### in BMMC
bk <- c(seq(-1.5,-0.1,by=0.01),seq(0,1.5,by=0.01))
df1 <- get_plotdf(assays,meta,method='ssgsea',metabGene='FALSE',type='BMMC')
p1 <- pheatmap(df1,scale = 'row',cluster_cols = T,cluster_rows = F,
               color=c(colorRampPalette(colors = c("#63A8D2","white"))(length(bk)/2),colorRampPalette(colors = c("white","#E88482"))(length(bk)/2)),
               legend_breaks=seq(-1.5,1.5,1),breaks=bk)

color <- colorRampPalette(c('grey',"white","red"))(100)
df2 <- get_plotdf(assays,meta,method='addMS',metabGene='FALSE',type='BMMC')
p2 <- pheatmap(df2,scale = 'row',cluster_cols = T,cluster_rows = F,color=color)

p <- cowplot::plot_grid(p1$gtable, p2$gtable, ncol= 2)
pdf('../irGSEA/Fig1_Heatmap_BMMC_superpathway.pdf',width = 7.5,height = 4)
print(p)
dev.off()

#### Fig1D TSNE of metablism ####
rm(list=ls())
pbmc_metab <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/pbmc_metab.RDS')
tmp <- RunUMAP(pbmc_metab,dims = 1:50,reduction = 'pca',n.neighbors = 50,min.dist = 0.5,seed.use = 123) 
tmp <- RunTSNE(pbmc_metab,dims = 1:50,seed.use = 123,check_duplicates = FALSE) 

#saveRDS(pbmc_metab,'/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/pbmc_metab.RDS')

pathway2gene <- read.gmt('/remote-home/yanzijun/CRU/TALL_M/data/GeneSet7.gmt')
genelst <- unique(pathway2gene$gene)
print(length(genelst))

## all cells
mycol <- c('#E88482','#63A8D2')
names(mycol) <- c('ALL','BMMC')

pdf('Fig1D_Dimplot_metab.pdf',width = 6,height =5)
TSNEPlot(pbmc_metab, group.by = "type",cols=mycol)
UMAPPlot(pbmc_metab, group.by = "type",cols=mycol)
DimPlot(pbmc_metab, group.by = "type",cols=mycol,reduction = "pca")
dev.off()


#### Fig1E violin plot of overall metab activity in ALL vs BMMC ####
rm(list=ls())
library(ggplot2)
library(dplyr)
library(RColorBrewer)

#setwd('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/')
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

## 导入数据
irGSEA <- readRDS('../irGSEA/pbmc_all_superpathway.RDS')
meta <- meta <- readRDS('meta.RDS') #要用最新的meta文件
assays <- irGSEA$assays

method='ssgsea'
df <- as.data.frame(t(as.matrix(assays[[method]]@counts)))
df_norm <- as.data.frame(scale(df,scale = TRUE,center = TRUE)) #因为是行为样本，列为pathway，所以进行列标准化
print(all(rownames(df_norm)==rownames(meta)))
df_norm$type <- meta$type
df_norm$type <- as.factor(df_norm$type)

data <- df_norm[,c('metabGene','type')]
med <- aggregate(round(data$metabGene,2),list(data$type),median)
rank <- med$Group.1[order(med$x,decreasing = T)]
print(rank)

plot.df <- data[,c('type','metabGene')]
colnames(plot.df) <- c('type','Activity')
plot.df$type=factor(plot.df$type,levels = rank)

p <- ggplot(plot.df, aes(x=type, y=Activity,fill=type)) + 
  geom_violin(trim = FALSE,width = 1)+
  #geom_boxplot(width = 0.5,outlier.shape = NA)+
  #stat_boxplot(geom = "errorbar",width=0.2)+
  stat_summary(fun= median, geom = "point",shape = 23, size = 2, color = "black")+ #添加中位值
  theme_classic()+labs(y='Metabolic activity score',x='')+
  theme(plot.title = element_text(hjust = 0.5,size = 8),
        axis.text = element_text(size=8,face="plain",color="black"),
        axis.title  = element_text(size=8,face="plain",color="black"),
        #axis.text.x = element_text(hjust = 1,angle =20),
        legend.position="none")+
  stat_compare_means()+
  scale_fill_manual(values = c('#E88482','#63A8D2'))

ggsave(paste('../irGSEA/Fig1_vlnplot_',method,'_all.pdf',sep=''),p,width = 4,height = 6,unit='cm')

#### Fig1F. bulk DEG ratio step5_2 ####

#### Fig1G violin plot of overall metab activity in celltype between ALL vs BMMC ####
rm(list=ls())
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggsignif)
library(ggpubr)
#setwd('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/')
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

## 导入数据
irGSEA <- readRDS('../irGSEA/pbmc_all_superpathway.RDS')
meta <- meta <- readRDS('meta.RDS') #要用最新的meta文件
assays <- irGSEA$assays

method='ssgsea'
df <- as.data.frame(t(as.matrix(assays[[method]]@counts)))
df_norm <- as.data.frame(scale(df,scale = TRUE,center = TRUE)) #因为是行为样本，列为pathway，所以进行列标准化
print(all(rownames(df_norm)==rownames(meta)))
df_norm$celltype <- meta$celltype_abbr
df_norm$celltype <- as.factor(df_norm$celltype)
df_norm$type <- meta$type
df_norm$type <- as.factor(df_norm$type)

## 将所有BMMC中的B/T，合并成ALL的ctrol数据
ctrl <- df_norm$metabGene[df_norm$celltype %in% c('T','B') & df_norm$type=='BMMC']
ctrl_data <- data.frame(metabGene=ctrl,type='BMMC',celltype='ALL')

data <- rbind(df_norm[,c('metabGene','type','celltype')],ctrl_data)
med <- aggregate(round(data$metabGene,2),list(data$celltype),median)
rank <- med$Group.1[order(med$x,decreasing = T)]
print(rank)

plot.df <- data[,c('type','celltype','metabGene')]
colnames(plot.df) <- c('type','celltype','Activity')
plot.df$celltype=factor(plot.df$celltype,levels = rank)

p <- ggplot(plot.df, aes(x=celltype, y=Activity,fill=type)) + 
  geom_violin(trim = FALSE,width = 1)+
  stat_summary(fun= median, aes(group=type),geom = "point",position=position_dodge(.95),
               shape = 23, size = 2, color = "black")+ #添加中位值
  theme_classic()+labs(y='Metabolic activity score',x='')+
  theme(plot.title = element_text(hjust = 0.5,size = 8),
        axis.text = element_text(size=8,face="plain",color="black"),
        axis.title  = element_text(size=8,face="plain",color="black"),
        #axis.text.x = element_text(hjust = 1,angle =20),
        legend.position="none")+
  stat_compare_means( aes(label = ..p.signif..))+
  scale_fill_manual(values = c('#E88482','#63A8D2'))

ggsave(paste('../irGSEA/Fig1_vlnplot_',method,'_ALLvsBMMC.pdf',sep=''),p,width = 12,height = 6,unit='cm')


#### Fig1H heatmap plot of superpathway metab activity in celltype between ALL vs BMMC ####
rm(list=ls())
library(ggplot2)
library(dplyr)
library(pheatmap)

#setwd('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/')
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

## 导入数据
irGSEA <- readRDS('../irGSEA/pbmc_all_superpathway.RDS')
meta <- meta <- readRDS('meta.RDS') #要用最新的meta文件
assays <- irGSEA$assays

method='ssgsea'
df <- as.data.frame(t(as.matrix(assays[[method]]@counts)))
df_norm <- as.data.frame(scale(df,scale = TRUE,center = FALSE)) #因为是行为样本，列为pathway，所以进行列标准化

print(all(rownames(df_norm)==rownames(meta)))
df_norm$celltype <- meta$celltype_abbr
df_norm$celltype <- as.factor(df_norm$celltype)
df_norm$type <- meta$type
df_norm$type <- as.factor(df_norm$type)
df_norm$metabGene=NULL

## 将所有BMMC中的B/T，合并成ALL的ctrol数据
ctrl <- df_norm[df_norm$celltype %in% c('T','B') & df_norm$type=='BMMC',]
ctrl$celltype <- 'ALL'

data <- rbind(df_norm,ctrl)

CT <- c('ALL','Ery.','T','B','Mono.','NK')
pathway <- colnames(df)[-8]

pvalue <- matrix(NA,nrow = length(CT),ncol=length(pathway))
fc <- matrix(NA,nrow = length(CT),ncol=length(pathway))

for(i in 1:length(CT)){
  print(i)
  ct=CT[i]
  for(j in 1:length(pathway) ){
    print(j)
    term=pathway[j]
    all <- data[data$type=='ALL' & data$celltype==ct,term]
    bmmc <- data[data$type=='BMMC' & data$celltype==ct,term]
    p <- wilcox.test(all,bmmc)$p.value
    pvalue[i,j] <- p
    
    fc[i,j] <- log(mean(all)/mean(bmmc),2)
  }
}

colnames(pvalue)=colnames(fc)=pathway
rownames(pvalue)=rownames(fc)=CT
print(fc)

if (!is.null(pvalue)){
  ssmt <- pvalue< 0.01
  pvalue[ssmt] <-'**'
  smt <- pvalue >0.01& pvalue <0.05
  pvalue[smt] <- '*'
  pvalue[!ssmt&!smt]<- ''
} else {
  pvalue <- F
}
#可视化
color <- colorRampPalette(c('#63A8D2',"white","#E88482"))(50)

p <- pheatmap(t(fc),scale = "none",cluster_row = T, cluster_col = T, border=NA,
         display_numbers = t(pvalue),fontsize_number = 12, number_color = "black",color=color)
ggsave(paste('../irGSEA/Fig1_heatmap_',method,'_ALLvsBMMC.pdf',sep=''),p,width = 5,height = 4)



