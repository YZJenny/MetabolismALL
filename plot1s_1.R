
rm(list=ls())
library(Seurat)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(gridExtra)
require(RColorBrewer)

pbmc <- readRDS('/remote-home/yanzijun/CRU/scTALL/res/Seurat/RDS/pbmc_QC.RDS') ##/remote-home/yanzijun/CRU/scTALL/scr/step1_QC.R
#### sFig1A: vlnbox of nFeature/nCount/percent.mt of each batch after QC ####

print(mean(pbmc$nCount_RNA))
print(mean(pbmc$nFeature_RNA))

print(median(pbmc$nCount_RNA))
print(median(pbmc$nFeature_RNA))

plot_QC_bybatch<- function(feature){
  p <- VlnPlot(pbmc,features =feature ,group.by = 'batch',pt.size = 0)+
    theme(axis.text=element_text(size=8,colour="black"),axis.title=element_text(size = 8,colour="black"),
          legend.position = 'none',
          plot.margin =unit(c(-0.4,0,-0.2,0.1), "cm"),
          axis.line = element_line(colour = "black", size = 0.3))+
    labs(x="",y="Numbers",title = "")
  return(p)
}

plot.nFeature_RNA <- plot_QC_bybatch('nFeature_RNA')
plot.nCount_RNA <- plot_QC_bybatch('nCount_RNA')
plot.percent.mt <-plot_QC_bybatch('percent.mt')
plot.QC <- grid.arrange(plot.nFeature_RNA,plot.nCount_RNA,plot.percent.mt,nrow = 1, ncol = 3)
ggsave('/remote-home/yanzijun/CRU/scTALL/res/Seurat/FIG/sFig1_vlnplot.pdf',plot.QC,width = 13.5,height = 4,units = 'cm')


#### sFig1B: vlnbox of nFeature/nCount/percent.mt of all cells after QC ####
mycol=brewer.pal(3, "Set3")
df <- data.frame(cells=colnames(pbmc),Gene='Gene',UMI='UMI',MT='MT',
                 nFeature_RNA=pbmc$nFeature_RNA,nCount_RNA=pbmc$nCount_RNA,perMT=pbmc$percent.mt)

plot.gene <- ggplot(df, aes(x=Gene, y=nFeature_RNA,fill=Gene)) + 
  geom_violin(trim=FALSE,color='black',size=0.3) + 
  geom_boxplot(width=0.3,position=position_dodge(0.9),size=0.3)+
  scale_fill_manual(values = mycol[1])+ 
  theme_classic()+ 
  theme(axis.text.x=element_blank(),axis.ticks.x =element_blank(),
        plot.margin =unit(c(0.2,0,0.2,0.1), "cm"),
        panel.grid = element_blank(),
        axis.text=element_text(size=8,colour="black"),axis.title=element_text(size = 8,colour="black"),
        axis.ticks.y =element_line(colour="black",size=0.25),
        axis.line = element_line(colour = "black", size = 0.3),legend.position = 'none')+ 
  labs(x="All",y="Numbers")

plot.UMI <- ggplot(df, aes(x=Gene, y=nCount_RNA,fill=UMI)) + 
  geom_violin(trim=FALSE,color='black',size=0.3) + 
  geom_boxplot(width=0.3,position=position_dodge(0.9),size=0.3)+
  scale_fill_manual(values = mycol[2])+ 
  theme_classic()+ 
  theme(axis.text.x=element_blank(),axis.ticks.x =element_blank(),
        plot.margin = unit(c(0.2,0,0.2,0.1), "cm"),
        panel.grid = element_blank(),
        axis.text=element_text(size=8,colour="black"),axis.title=element_text(size = 8,colour="black"),
        axis.ticks.y =element_line(colour="black",size=0.25),
        axis.line = element_line(colour = "black", size = 0.3),legend.position = 'none')+ 
  labs(x="All",y="Numbers")

plot.MT <- ggplot(df, aes(x=MT, y=perMT,fill=MT)) + 
  geom_violin(trim=FALSE,color='black',size=0.3) + 
  geom_boxplot(width=0.3,position=position_dodge(0.9),size=0.3)+
  scale_fill_manual(values = mycol[3])+ 
  theme_classic()+ 
  theme(axis.text.x=element_blank(),axis.ticks.x =element_blank(),
        plot.margin = unit(c(0.2,0,0.2,0.1), "cm"),
        panel.grid = element_blank(),
        axis.text=element_text(size=8,colour="black"),axis.title=element_text(size = 8,colour="black"),
        axis.ticks.y =element_line(colour="black",size=0.25),
        axis.line = element_line(colour = "black", size = 0.3),
        legend.position = 'none')+ 
  labs(x="All",y="Numbers")

library(gridExtra)
plot.all <- grid.arrange(plot.gene,plot.UMI, plot.MT, nrow = 1, ncol = 3)
ggsave('/remote-home/yanzijun/CRU/scTALL/res/Seurat/FIG/sFig1D.pdf',plot.gene,width = 6,height = 4,units = 'cm')
