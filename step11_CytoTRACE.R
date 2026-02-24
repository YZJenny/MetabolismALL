rm(list=ls())
library(Seurat)
library(CytoTRACE2)
library(ggplot2)
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/CytoTRACE2/')

## 导入sc数据
pbmc_tumor <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/pbmc_tumor_metab.RDS')
meta <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/meta.RDS') #用最新meta
map <- meta$celltype_abbr
names(map) <- rownames(meta)
pbmc_tumor$celltype_abbr <- map[colnames(pbmc_tumor)]

subtype <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/TCAsubtype_new/subtype_info_kmeans.RDS')
print(all(colnames(pbmc_tumor)==rownames(subtype)))
pbmc_tumor$subtype <- subtype$subtype

#### 1.run CytoTRACE ####
Cyto <- cytotrace2(pbmc_tumor, is_seurat = TRUE, slot_type = "counts", species = 'human',seed = 1234)
Cyto
saveRDS(Cyto,'CytoTRACE.RDS')

#### 2.可视化 ####
Cyto <- readRDS('CytoTRACE.RDS')
##### Fig2F.boxplot #####
library(ggpubr)
print(all(colnames(Cyto)==rownames(subtype)))
Cyto@meta.data$subtype <- subtype$subtype
Cyto@meta.data$subtype <- factor(Cyto@meta.data$subtype,levels = c('High','Medium','Low'))
p1 <- ggboxplot(Cyto@meta.data, x="subtype", y="CytoTRACE2_Score", width = 0.6, 
                color = "black",#轮廓颜色
                fill="subtype",#填充
                xlab = F, #不显示x轴的标签
                bxp.errorbar=F,#显示误差条
                outlier.shape=NA, #不显示outlier
                legend = "none")+
  #scale_fill_manual(values = c('#E88482','grey','#63A8D2'))+
  scale_fill_manual(values = c('grey','grey','grey'))+
  ylab('Potency score')
###指定组比较
my_comparisons <- list(c("High", "Medium"), c("High", "Low"))
p2 <- p1+stat_compare_means(comparisons = my_comparisons,
                            method = "wilcox.test")
ggsave('CytoTRACE_boxplot.pdf',p2,width = 3,height = 3)


##### sFig2G. CytoTRACE与TCA score的相关性 #####
scale_rows<- function(x){ #行为基因，列为样本
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
df$CytoTRACE2_Score <- Cyto@meta.data$CytoTRACE2_Score
df$metabGene <- as.numeric(ssgsea_tumor[which(rownames(ssgsea_tumor)=='metabGene'),])
df$TCA <- as.numeric(ssgsea_tumor[which(rownames(ssgsea_tumor)=='TCA cycle'),])
df$subtype <- subtype$subtype

print(cor.test(df$CytoTRACE2_Score,df$metabGene,method = 'spearman')) #0.306,p-value < 2.2e-16
print(cor.test(df$CytoTRACE2_Score,df$TCA,method = 'spearman')) #0.313,p-value < 2.2e-16
aggregate(df$CytoTRACE2_Score, by=list(type=df$subtype),mean)
aggregate(df$CytoTRACE2_Score, by=list(type=df$subtype),median)

p1 <- ggplot(df, aes(x = CytoTRACE2_Score, y = metabGene)) +
  geom_point()+
  geom_density_2d(alpha = 0.5)+
  geom_density_2d_filled()
p2 <- ggplot(df, aes(x = CytoTRACE2_Score, y = TCA)) +
  geom_point()+
  geom_density_2d(alpha = 0.5)+
  geom_density_2d_filled()
pdf('CytoTRACE_densitycor.pdf',width = 8,height = 7)
print(p1)
print(p2)
dev.off()





p1 <- p.cor <- ggplot(data=df, aes(x=CytoTRACE2_Score, y=metabGene)) + 
  geom_point()+ 
  geom_smooth(method = 'lm', formula = y ~ x, se = T) + ##se= T意思为画出置信区间, se=F意思则为不画出置信区间
  stat_cor(data=df, method = "spearman")+
  theme_classic(base_size = 10)+
  theme(plot.title = element_text(hjust = 0.5,color='black'),
        axis.text = element_text(size = 10,color='black'))
pdf('CytoTRACE_densitycor.pdf',width = 8,height = 7)
print(p1)
dev.off()




##### sFig2F.CytoTRACE2_Score #####
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

