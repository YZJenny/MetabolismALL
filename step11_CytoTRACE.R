rm(list=ls())
library(Seurat)
library(CytoTRACE2)
library(ggplot2)
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/CytoTRACE2/')

## load data
pbmc_tumor <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/pbmc_tumor_metab.RDS')
meta <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/meta.RDS') #用最新meta
map <- meta$celltype_abbr
names(map) <- rownames(meta)
pbmc_tumor$celltype_abbr <- map[colnames(pbmc_tumor)]

subtype <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/TCAsubtype_new/subtype_info.RDS')
print(all(colnames(pbmc_tumor)==rownames(subtype)))
pbmc_tumor$subtype <- subtype$subtype

#### 1.run CytoTRACE ####
Cyto <- cytotrace2(pbmc_tumor, is_seurat = TRUE, slot_type = "counts", species = 'human',seed = 1234)
Cyto
saveRDS(Cyto,'CytoTRACE.RDS')

#### 2.plot ####
Cyto <- readRDS('CytoTRACE.RDS')
##### Fig3D (right).boxplot #####
library(ggpubr)
print(all(colnames(Cyto)==rownames(subtype)))
Cyto@meta.data$subtype <- subtype$subtype
Cyto@meta.data$subtype <- factor(Cyto@meta.data$subtype,levels = c('High','Medium','Low'))
p1 <- ggboxplot(Cyto@meta.data, x="subtype", y="CytoTRACE2_Score", width = 0.6, 
                color = "black",
                fill="subtype",
                xlab = F, 
                bxp.errorbar=F,
                outlier.shape=NA, 
                legend = "none")+
  #scale_fill_manual(values = c('#E88482','grey','#63A8D2'))+
  scale_fill_manual(values = c('grey','grey','grey'))+
  ylab('Potency score')
###指定组比较
my_comparisons <- list(c("High", "Medium"), c("High", "Low"))
p2 <- p1+stat_compare_means(comparisons = my_comparisons,
                            method = "wilcox.test")
ggsave('CytoTRACE_boxplot.pdf',p2,width = 3,height = 3)



##### Fig3D (left).CytoTRACE2_Score #####
library(ggplot2)
p1 <- ggplot(df,aes(x=tSNE_1,y=tSNE_2,color=CytoTRACE2_Score))+
  geom_point(size=0.5)+
  #scale_colour_gradient2(low = 'grey100',mid='grey', high = 'darkblue')+
  scale_colour_gradientn(colours = 
                           (rev(c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", 
                              "#66C2A5", "#5E4FA2"))), 
                         na.value = "transparent", 
                         limits = c(0, 0.8), 
                         breaks = seq(0, 0.8, by = 0.2), 
                         labels = c("0.0 (More diff.)", 
                                    "0.2", "0.4", "0.6",  "0.8 (Less diff.)"), 
                         name = "Relative\norder \n", 
                         guide = guide_colorbar(frame.colour = "black", 
                                                ticks.colour = "black"))+
  theme_classic()
pdf('CytoTRACE_TSNE.pdf',width = 8,height = 7)
print(p1)
dev.off()

