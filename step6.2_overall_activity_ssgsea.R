####1. 所有代谢通路基因的总的代谢活性(by ssgsea for Fig1)

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

#### 1.1 可视化：所有细胞，分type的celltype boxplot #####
plot.df <- score[,c('type','celltype_2','metabGene1')]
colnames(plot.df) <- c('type','celltype','Activity')
plot.df$celltype[plot.df$celltype=='leukemia cells'] <- 'ALL cells'
plot.df$logActivity <- log(plot.df$Activity,10)

mycol <- c('#E88482','#63A8D2')
names(mycol) <- c('ALL','BMMC')

p <- ggplot(plot.df, aes(x = celltype, y = logActivity)) + 
  geom_violin(aes(fill = type), trim = FALSE)+
  stat_summary(fun= median, aes(group=type),geom = "point",position=position_dodge(.9),
               shape = 23, size = 1, color = "black")+ #添加中位值
  theme_classic()+labs(y='Metabolic activity score',x='')+
  theme(plot.title = element_text(hjust = 0.5,size = 8),
        axis.text = element_text(size=8,face="plain",color="black"),
        axis.title  = element_text(size=8,face="plain",color="black"),
        axis.text.x = element_text(hjust = 1,angle =20))+
  scale_fill_manual(values = mycol)

pdf('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/irGSEA/Fig1_vlnplot_ssgsea_allcell.pdf',
    width = 4.2,height = 2.1)
#VlnPlot(pbmc,features = 'metabGene1',split.by = 'type',group.by = 'celltype_2',pt.size = 0)
print(p)
dev.off()

## 做显著性检验
CT <- unique(plot.df$celltype)
for(ct in CT){
  df <- plot.df[plot.df$celltype==ct,]
  if(length(unique(df$type))==2){
    print(ct)
    ALL <- as.numeric(df$Activity[df$type=='ALL'])
    BMMC <- as.numeric(df$Activity[df$type=='BMMC'])
    res <- wilcox.test(ALL,BMMC,alternative = 'less')
    print(res$p.value)
    #print(summary(ALL))
    #print(summary(BMMC))
  }
}
## 只有NK/Mono是ALL显著大于BMMC


#### 1.2 ALL样本中 celltype boxplot #####
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


#### 1.3 BMMC 样本中 celltype boxplot #####
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



####### 初始画图
pdf('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/irGSEA/vlnplot_addModulescore.pdf',width = 10,height = 4)
VlnPlot(pbmc,features = 'metabGene1',split.by = 'type',group.by = 'celltype_2')
VlnPlot(pbmc_tumor,features = 'metabGene1',group.by = 'orig.ident')
VlnPlot(pbmc_normal,features = 'metabGene1',split.by = 'type',group.by = 'celltype_2')
VlnPlot(pbmc_ALL,features = 'metabGene1',group.by = 'celltype_2')
VlnPlot(pbmc_BMMC,features = 'metabGene1',group.by = 'celltype_2')
dev.off()

