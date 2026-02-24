### CV_SDA+GSEA肿瘤间和肿瘤内异质性: tumor cells by batch + normal cells in ALL
rm(list=ls())
library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(dplyr)
options(stringsAsFactors=FALSE)
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

pathway2gene <- read.gmt('/remote-home/yanzijun/CRU/TALL_M/data/GeneSet7.gmt')
genelst <- unique(pathway2gene$gene)
print(length(genelst))

####1. Loading the data ####
pbmc <- readRDS('pbmc_metab.RDS')
meta <- readRDS('meta.RDS')
print(all(colnames(pbmc)==rownames(meta)))
pbmc$celltype_abbr <- meta$celltype_abbr


type='ALL' ##肿瘤样本
if(type=='all'){
  pbmc_metab <- pbmc[genelst,]
}else if(type=='ALL'){
  pbmc_metab <- pbmc[genelst,Cells(pbmc)[-grep('^PBMMC',pbmc$orig.ident)]]
}

print(table(pbmc_metab$orig.ident))


df <- pbmc_metab@meta.data[,c('orig.ident','celltype_abbr')]
df <- df %>% mutate(celltype_abbr=if_else(celltype_abbr=='ALL',orig.ident,celltype_abbr))
print(all(rownames(df)==colnames(pbmc_metab)))
pbmc_metab$celltype_3 <- df$celltype_abbr

samples <- c(sort(unique(pbmc_metab$orig.ident)),
             sort(unique(as.character(pbmc_metab$celltype_abbr))[-1]))
print(samples)


#### 2.CV and SD ####
get_gsea <- function(s,type,pbmc_metab){
  print(s)
  if(s=='allsample'){ ## 所有样本的tumor cell合在一起，看肿瘤间异质性
    each_metabolic_sce <- pbmc_metab[,pbmc_metab$celltype_abbr=='ALL']
  }else{
    each_metabolic_sce <- pbmc_metab[,pbmc_metab$celltype_3==s]
  }
  
  each_metabolic_tpm <- each_metabolic_sce@assays$RNA@data
  each_metabolic_tpm <- each_metabolic_tpm[rowSums(each_metabolic_tpm)>0,]
  x <- each_metabolic_tpm
  
  # 计算SD
  gene_sd <- apply(x, 1, sd)
  gene_sd <- sort(gene_sd, decreasing = T)
  print(head(gene_sd))
  
  # 计算CV
  gene_mean <- apply(x, 1, mean)
  gene_cv <- gene_sd / gene_mean
  gene_cv <- sort(gene_cv, decreasing = T)
  print(head(gene_cv))
  
  #### GSEA of SD
  gsea_sd <- GSEA(geneList = gene_sd, TERM2GENE=pathway2gene, verbose=FALSE,maxGSSize = 2000,minGSSize = 10,
                  pvalueCutoff = 1,pAdjustMethod='fdr',seed = FALSE,nPermSimple=1000) #有些通路为NA
  gsea_sd <- as.data.frame(gsea_sd)
  
  #### GSEA of CV
  gsea_cv <- GSEA(geneList = gene_cv, TERM2GENE=pathway2gene, verbose=FALSE,maxGSSize = 2000,minGSSize = 10,
                  pvalueCutoff = 1,pAdjustMethod='fdr',seed = FALSE,nPermSimple=1000) #有些通路为NA
  gsea_cv <- as.data.frame(gsea_cv)
  
  gsea <- list(gsea_sd=gsea_sd,gsea_cv=gsea_cv)
}

##### 2.1 每个肿瘤内的异质性 #####
enrich_sd <- data.frame(x=NULL,y=NULL,NES=NULL,PVAL=NULL)
enrich_cv <- data.frame(x=NULL,y=NULL,NES=NULL,PVAL=NULL)
for(s in samples){
  gsea  <- get_gsea(s,type,pbmc_metab)
  
  gsea_sd <- gsea$gsea_sd
  saveRDS(gsea_sd,paste('../CV_SD/',s,'_',type,'_SD.RDS',sep=''))
  
  # get the result
  gsea_pathways <- rownames(gsea_sd)
  enrich_sd <- rbind(enrich_sd,
                     data.frame(x=s,y=gsea_pathways,NES=gsea_sd$NES,PVAL=gsea_sd$p.adjust))
  
  gsea_cv <- gsea$gsea_sd
  saveRDS(gsea_cv,paste('../CV_SD/',s,'_',type,'_CV.RDS',sep=''))
  
  # get the result
  gsea_pathways <- rownames(gsea_cv)
  enrich_cv <- rbind(enrich_cv,
                     data.frame(x=s,y=gsea_pathways,NES=gsea_cv$NES,PVAL=gsea_cv$p.adjust))
  
}
saveRDS(enrich_sd,paste('../CV_SD/enriched_pathway','_',type,'_SD.RDS',sep=''))
saveRDS(enrich_cv,paste('../CV_SD/enriched_pathway','_',type,'_CV.RDS',sep=''))

##### 2.2 肿瘤间的异质性 #####
gsea  <- get_gsea('allsample',type,pbmc_metab)
gsea_sd <- gsea$gsea_sd
saveRDS(gsea_sd,paste('../CV_SD/allsample_',type,'_SD.RDS',sep=''))

gsea_cv <- gsea$gsea_cv
saveRDS(gsea_sd,paste('../CV_SD/allsample_',type,'_CV.RDS',sep=''))



#### 3. 可视化 (sFig2C-D, Related to Fig2C) ####
rm(list=ls())
library(ggplot2)
library(clusterProfiler)
library(scater)
library(stringr)
library(scran)
library(gtools)
options(stringsAsFactors=FALSE)
#setwd('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/')
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

## 只挑选p<0.05的pathway可视化
# min_pval <- by(enrich_data_df$PVAL, enrich_data_df$y, FUN=min)
# select_pathways <- names(min_pval)[(min_pval<=0.05)] #
# print(select_pathways)
# plot_df <- enrich_data_df[enrich_data_df$y%in% select_pathways,]

##### 3.1 肿瘤内异质性 ##### 
type='ALL'
method='CV'
plot_df <- readRDS(paste('../CV_SD/enriched_pathway','_',type,'_',method,'.RDS',sep=''))
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
ggsave(paste("../CV_SD/sFig_enriched_pathway_",type,'_',method,".pdf",sep=''),p,
       width = 5.5,height=3.3,units="in",device="pdf",useDingbats=FALSE)

##### 3.2 肿瘤间异质性 #####
type='ALL'
method='SD'
plot_df <- readRDS(paste('../CV_SD/allsample','_',type,'_',method,'.RDS',sep=''))
pvals <- plot_df$p.adjust
plot_df$PVAL <- -log10(pvals)
plot_df$ID <- factor(plot_df$ID,levels = rev(plot_df$ID))

##barplot
p <- ggplot(plot_df, aes(x = ID, y = PVAL, fill=PVAL)) +
  geom_bar(stat = "identity")+
  coord_flip()+
  labs(x = NULL, y = NULL,size = "-log10 adjp") +
  scale_fill_gradient(low = "white", high = "red") +
  theme(axis.text  =element_text(colour="black", size = 8))+
  theme_classic()+
  geom_hline(aes(yintercept=-log(0.05,10)), colour="#990000", linetype="dashed")

ggsave(paste("../CV_SD/sFig_allsample_",type,'_',method,".pdf",sep=''),p,
       width = 5.5,height=3.3,units="in",device="pdf",useDingbats=FALSE)