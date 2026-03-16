## Fig2
rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(RColorBrewer)
library(clusterProfiler)
col1 <-colorRampPalette(brewer.pal(8,'Paired'))(17)

setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

pbmc_tumor <- readRDS('pbmc_tumor_metab.RDS')
meta <- readRDS('meta.RDS') #要用最新的meta文件
map <- meta$celltype_abbr
names(map) <- rownames(meta)
pbmc_tumor$celltype_abbr <- map[colnames(pbmc_tumor)]


#### Fig2A (upper pannel). tumor cells in ALL ####
pdf('Fig2A_TSNE_tumor_metab.pdf', width = 6.5,height =5)
DimPlot(pbmc_tumor,group.by = "orig.ident",cols = col1,reduction = 'tsne')
dev.off()


#### Fig2A (lower pannel). TSNE of metaGene ssgsea score####
scale_rows<- function(x){ #row: gene, col: sample
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m)/s)
}

score_irGSEA <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/irGSEA/pbmc_all_superpathway.RDS')
ssgsea <- score_irGSEA$assays$ssgsea
meta <- score_irGSEA$meta
ssgsea_tumor <- scale_rows(as.matrix(ssgsea[,meta$celltype_2=='leukemia cells']))

all(colnames(ssgsea_tumor)==colnames(pbmc_tumor))
df <- as.data.frame(pbmc_tumor@reductions$tsne@cell.embeddings)
df$metabGene <- as.numeric(ssgsea_tumor[which(rownames(ssgsea_tumor)=='metabGene'),])
df$TCA <- as.numeric(ssgsea_tumor[which(rownames(ssgsea_tumor)=='TCA cycle'),])

p1 <- ggplot(df,aes(x=tSNE_1,y=tSNE_2,color=metabGene))+
  geom_point(size=0.5)+
  scale_colour_gradientn(colours = 
                           (rev(c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", 
                                  "#66C2A5", "#5E4FA2"))), 
                         na.value = "transparent", 
                         limits = c(-4, 4), 
                         breaks = seq(-4, 4, by = 2), 
                         guide = guide_colorbar(frame.colour = "black", 
                                                ticks.colour = "black"))+
  theme_classic()


pdf('/remote-home/yanzijun/CRU/ped_M/res_ALL/irGSEA/TSNE_ssgsea.pdf',width = 8,height = 7)
print(p1)
dev.off()


#### Fig2B. dotplot of heterogenity: related to step7_ALL_heter_PCA.R ####
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
plot_df <- readRDS(paste('../PC/enriched_pathway','_',type,'.RDS',sep=''))
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
ggsave(paste("../PC/Fig2B_enriched_pathway_",type,".pdf"),p,
       width = 5.5,height=3.3,units="in",device="pdf",useDingbats=FALSE)


#### Fig2C. TSNE ALL samples by TCA subtype ####
new_meta <- readRDS('../TCAsubtype_new/subtype_info.RDS') 
print(table(new_meta$subtype))

print(all(colnames(pbmc_tumor)==rownames(new_meta)))
pbmc_tumor$subtype <- new_meta$subtype

pdf('Fig2_TSNE_tumor_metab_bySubtype.pdf', width = 6,height =5)
TSNEPlot(pbmc_tumor,group.by = "subtype",cols = c('#E88482','#63A8D2','#DADADA'))
dev.off()

