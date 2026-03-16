rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(tidyverse)
library(future)
plan('multiprocess',workers=80)
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')


pbmc <- readRDS('pbmc.RDS')
meta <- readRDS('meta.RDS')
pbmc@meta.data <- meta

Idents(pbmc) <- 'orig.ident'
case_ID <- unique(pbmc$orig.ident)[-grep('^PBMMC',unique(pbmc$orig.ident))]
control_ID <- unique(pbmc$orig.ident)[grep('^PBMMC',unique(pbmc$orig.ident))]

mk <- FindMarkers(pbmc,ident.1 = case_ID,ident.2 = control_ID,
                  logfc.threshold = 0,min.pct = 0,min.cells.feature = 0,
                  min.cells.group = 0,only.pos = FALSE)
saveRDS(mk,'mk_ALL_BMMC_bysample.RDS')


#### 1. barplot of DEG ratio ####
rm(list=ls())
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(tidyverse)
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

pathway2gene <- read.gmt('/remote-home/yanzijun/CRU/TALL_M/data/GeneSet7.gmt')
mk <- readRDS('mk_ALL_BMMC_bysample.RDS')


p_val_adj <- 0.05
logfc <- 0
data <- rownames_to_column(mk,'gene')
data <- merge(data,pathway2gene)

data$group <- 'Nosig'
data$group[data$avg_log2FC > logfc & data$p_val_adj < p_val_adj] <- 'Up'
data$group[data$avg_log2FC < -logfc & data$p_val_adj < p_val_adj] <- 'Down'

plot_df <- as.data.frame(table(data$group,data$term))
colnames(plot_df) <- c('type','term','Freq')
plot_df$type <- factor(plot_df$type,levels = c('Nosig','Down','Up'))
plot_df$term <- factor(plot_df$term,levels = rev(c('TCA cycle','Nucleotide','Amino acid','Carbohydrate',
                                                   'Vitamin cofactor','Lipid','Energy')))


mycol <- c('#9cb0c3','#DADADA','#eab080')
names(mycol) <- c('Down','Nosig','Up')

p <- ggplot(plot_df,aes(x=term,y=Freq,fill=type))+
  geom_col(position = "fill", width = 0.6)+
  coord_flip()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(color='black',size = 12))+
  scale_fill_manual(values = mycol)+
  labs(x='',y='Proportion of cells')+theme(legend.position="right")+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

ggsave('sFig1F_DEGratio_ALL_BMMC_bysample.pdf',p,width = 4,height = 2.5,units = 'in')

