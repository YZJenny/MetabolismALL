#### 方式2. ALL vs BMMC (区分sample:用bulk方法分析DEG)
rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(tidyverse)
library(future)
plan('multiprocess',workers=80)
#setwd('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/')
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


#### Fig1E. ALL vs BMMC 条形ratio柱状图 ####
rm(list=ls())
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(tidyverse)
#setwd('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/')
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


mycol <- c('#63A8D2','#DADADA','#E88482')
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

ggsave('Fig1E_DEGratio_ALL_BMMC_bysample.pdf',p,width = 4,height = 2.5,units = 'in')


#### FIG1F. ALL vs BMMC GSEA图 ####
rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
#setwd('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/')
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

pathway2gene <- read.gmt('/remote-home/yanzijun/CRU/TALL_M/data/GeneSet7.gmt')
mk <- readRDS('mk_ALL_BMMC_bysample.RDS')
mk <- na.omit(mk)
mk <- mk[mk$p_val_adj<0.05,]

gs <-bitr(rownames(mk), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db",drop = TRUE)
markers1<-cbind(mk[gs[,1],],gs)

# markers1 <- mk
# markers1$SYMBOL <- rownames(markers1)

markers1 <- markers1[!duplicated(markers1$SYMBOL),]
geneList = markers1$avg_log2FC
names(geneList) = markers1$SYMBOL
geneList = sort(geneList,decreasing = T)

egmt <- GSEA(geneList, TERM2GENE=pathway2gene,verbose=F,maxGSSize=1000,pvalueCutoff = 1)
#egmt <- fgsea::fgsea(pathways = split(pathway2gene$gene, pathway2gene$term),geneList)
y=data.frame(egmt)
print(y[,1:7])

saveRDS(egmt,'GSEA_ALL_BMMC_bysample.RDS') #GSE_nc用GSEA_ALL_BMMC_bysample_2.RDS


egmt <- readRDS('GSEA_ALL_BMMC_bysample.RDS')
y=data.frame(egmt)

library(enrichplot)
pdf('GSEAplot_ALL_BMMC_bysample.pdf',width = 4, height = 3.6)
for(i in seq_along(egmt@result$ID)){
  p <- gseaplot2(egmt, geneSetID = i, title = egmt@result$ID[i],pvalue_table = TRUE,color = 'black')
  print(p)
}
dev.off()


library(RColorBrewer)
cols <- rev(colorRampPalette(c("grey","#E88482"))(10))
y$Description <- factor(y$Description,levels = rev(y$Description))
p <-ggplot(data=y, aes(x=Description , y=NES ,fill=pvalue))+
  geom_bar(stat="identity",width = 0.7)+coord_flip()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(color='black',size = 8),
        axis.title  = element_text(color='black',size = 8))+
  scale_fill_gradientn(colours = cols)
ggsave("Fig1F_GSEAbarplot_ALL_BMMC_bysample.pdf",width = 4,height = 2.5)


pdf("GSEAdotplot_ALL_BMMC_bysample.pdf",
    width = 4,height = 3)#气泡图，展示geneset被激活还是抑制
dotplot(egmt,split=".sign")+facet_grid(~.sign)
dev.off()

