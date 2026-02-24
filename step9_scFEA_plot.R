####1. 所有代谢通路基因的总的代谢通量
rm(list=ls())
library(Seurat)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(RColorBrewer)
library(ggpubr)
library(ggsignif)

col <-colorRampPalette(brewer.pal(8,'Paired'))(8)
#setwd('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/scFEA/')
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/scFEA/')

CT <- c('ALL cells','B cells',"Erythrocytes ", "Monocytes" , "NK cells" , "pDCs" )
print(CT)
mycol <- c('#A6CEE3','#33A02C','#B2DF8A','#FB9A99','#1F78B4','#E31A1C')
names(mycol) <- CT

scFEA <- read.csv('ALL_sc_sample5k_module168.csv',
                  row.names = 1)
scFEA[1:3,1:4]
all.meta <- read.csv('ALL_meta_sample5k.csv',row.names =1)
M2C <- read.csv("/remote-home/yanzijun/software/scFEA-master/data/Human_M168_information.symbols.csv",row.names = 1)


########################可视化
score_tumor <- scFEA[all.meta$celltype_2 =='leukemia cells',]
score_normal <- scFEA[all.meta$celltype_2 !='leukemia cells',]

score_ALL <- scFEA[-grep('^PBMMC',all.meta$orig.ident),]
score_BMMC <- scFEA[grep('^PBMMC',all.meta$orig.ident),]

#### 1.1 可视化：所有细胞，分type的celltype boxplot #####
data=data.frame(flux=apply(scFEA,1,mean),type=all.meta$type,celltype=all.meta$celltype_2)
data$celltype[data$celltype=='leukemia cells'] <- 'ALL cells'

plot.df <- data
p <- ggplot(plot.df, aes(x = celltype, y = flux)) + 
  geom_violin(aes(fill = type), trim = FALSE)+
  stat_summary(fun= median, aes(group=type),geom = "point",position=position_dodge(.9),
               shape = 23, size = 2, color = "black")+ #添加中位值
  theme_classic()+labs(y='Metabolic flux',x='')+
  theme(plot.title = element_text(hjust = 0.5,size = 12),
        axis.text = element_text(size=12,face="plain",color="black"),
        axis.title  = element_text(size=12,face="plain",color="black"),
        axis.text.x = element_text(hjust = 1,angle =20))

pdf('vlnplot_flux_allcell.pdf', width = 7,height = 4.1)
print(p)
dev.off()


#### 1.2 可视化：所有细胞，分type的celltype barplot #####
df <- plot.df
print(kruskal.test(flux~celltype,data=df)$p.value)

# 统计 数据数据均值、标准差、标准误
mean <- aggregate(df$flux, by=list(df$type,df$celltype), FUN=mean)
sd <- aggregate(df$flux, by=list(df$type,df$celltype), FUN=sd)

df_res <- data.frame(mean, sd=sd$x)
colnames(df_res) = c("type",'celltype',  "Mean", "Sd")


p <-  ggplot(df_res, aes(x=celltype, y=Mean, fill=type)) +
  geom_bar(stat="identity", position=position_dodge(),color="black", width=.8) +
  #scale_fill_manual(values =mycol) +
  geom_errorbar(aes(ymin=Mean-Sd, ymax=Mean +Sd),position=position_dodge(.8), width=.4) +
  theme_classic()+
  labs(y='Metabolic flux')+
  theme(plot.title = element_text(hjust = 0.5,size = 8),
        axis.text = element_text(size=8,face="plain",color="black"),
        axis.title  = element_text(size=8,face="plain",color="black"),
        axis.text.x = element_text(hjust = 1,angle =90),
        legend.position="none")
pdf('barplot_flux_allcell.pdf',width = 5,height = 3.1)
print(p)
dev.off()


## 做显著性检验
data <- plot.df
CT <- unique(plot.df$celltype)
for(ct in CT){
  df <- data[data$celltype==ct,]
  if(length(unique(df$type))==2){
    print(ct)
    ALL <- as.numeric(df$flux[df$type=='ALL'])
    BMMC <- as.numeric(df$flux[df$type=='BMMC'])
    res <- wilcox.test(ALL,BMMC,alternative = 'greater')
    print(res$p.value)
    # print(summary(ALL))
    # print(summary(BMMC))
  }
}
## 只有Bcell是ALL显著小于BMMC


#### 1.2 ALL样本中 celltype boxplot #####
plot.df <- data[data$type=='ALL',]

med <- aggregate(round(plot.df$flux,6),list(plot.df$celltype),median)
rank <- med$Group.1[order(med$x,decreasing = T)]

plot.df$celltype=factor(plot.df$celltype,levels = rank)

my_comparisons <- list( c("ALL cells", "pDCs"), c("ALL cells", "Monocytes"), c("ALL cells", "B cells"),
                        c("ALL cells", "NK cells"), c("ALL cells", "Erythrocytes"))

p <- ggplot(plot.df, aes(x=celltype, y=flux,fill=celltype)) + 
  #geom_violin(trim = FALSE,width = 5)+
  geom_boxplot(width = 0.5,outlier.shape = NA)+
  #stat_boxplot(geom = "errorbar",width=0.2)+
  #stat_summary(fun= median, geom = "point",shape = 23, size = 2, color = "black")+ #添加中位值
  theme_classic()+labs(y='Metabolic flux',x='')+
  theme(plot.title = element_text(hjust = 0.5,size = 12),
        axis.text = element_text(size=12,face="plain",color="black"),
        axis.title  = element_text(size=12,face="plain",color="black"),
        axis.text.x = element_text(hjust = 1,angle =20),
        legend.position="none")+
  geom_signif(comparisons = my_comparisons,
              test = "wilcox.test", test.args = c("great"),
              y_position=seq(from=max(plot.df$flux),
                             to=max(plot.df$flux)+4*0.004,
                             by=0.04),
              map_signif_level=TRUE)+
  scale_fill_manual(values =mycol)

pdf('sFig_boxplot_flux_ALL.pdf',width = 5,height = 4)
print(p)
dev.off()


#### 1.3 BMMC 样本中 celltype boxplot #####
plot.df <- data[data$type=='BMMC',]

med <- aggregate(round(plot.df$flux,6),list(plot.df$celltype),median)
rank <- med$Group.1[order(med$x,decreasing = T)]

plot.df$celltype=factor(plot.df$celltype,levels = rank)

CT=c('B cells',"Erythrocytes", "Monocytes" , "NK cells" )
mycol = c('#33A02C','#B2DF8A','#FB9A99','#1F78B4','#E31A1C')
names(mycol) <- CT

p <- ggplot(plot.df, aes(x=celltype, y=flux,fill=celltype)) + 
  #geom_violin(trim = FALSE,width = 1)+
  geom_boxplot(width = 0.5,outlier.shape = NA)+
  #stat_boxplot(geom = "errorbar",width=0.2)+
  #stat_summary(fun= median, geom = "point",shape = 23, size = 2, color = "black")+ #添加中位值
  theme_classic()+labs(y='Metabolic activity score',x='')+
  theme(plot.title = element_text(hjust = 0.5,size = 12),
        axis.text = element_text(size=12,face="plain",color="black"),
        axis.title  = element_text(size=12,face="plain",color="black"),
        axis.text.x = element_text(hjust = 1,angle = 20),
        legend.position="none")+
  stat_compare_means()+
  scale_fill_manual(values =mycol[rank])

pdf('vlnplot_flux_BMMC.pdf',width = 4,height = 4)
print(p)
dev.off()
