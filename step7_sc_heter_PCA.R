### PCA+GSEA肿瘤间和肿瘤内异质性: tumor cells by batch + normal cells in ALL
rm(list=ls())
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
library(scater)
library(stringr)
library(scran)
library(gtools)
library(dplyr)
options(stringsAsFactors=FALSE)
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

pathway2gene <- read.gmt('/remote-home/yanzijun/CRU/TALL_M/data/GeneSet7.gmt')
genelst <- unique(pathway2gene$gene)
print(length(genelst))

#### 1. Loading the data ####
pbmc <- readRDS('pbmc_metab.RDS')
meta <- readRDS('meta.RDS')
print(all(colnames(pbmc)==rownames(meta)))
pbmc$celltype_abbr <- meta$celltype_abbr

type='ALL'
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

#### 2.PCA ####
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
  ntop <- nrow(x)
  rv <- rowVars(x)
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(x[select,]))
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  ### select PCs that explain at least 80% of the variance
  cum_var <- cumsum(percentVar)
  select_pcs <- which(cum_var>0.8)[1]
  
  ### plot the PCA and explained variances
  tmp_plotdata <- data.frame(x=1:length(percentVar),y=percentVar,
                             sel=c(rep("y",select_pcs),rep("n",length(percentVar)-select_pcs)),
                             celltype=rep(s,length(percentVar)))
  
  ### GSEA
  pre_rank_matrix <- as.matrix(rowSums(abs(pca$rotation[,1:select_pcs])))
  print(dim(pre_rank_matrix))
  
  genes <- as.numeric(pre_rank_matrix[,1])
  names(genes) <- rownames(pre_rank_matrix)
  genes <- sort(genes, decreasing = T)
  
  gsea <- GSEA(geneList = genes, TERM2GENE=pathway2gene, verbose=FALSE,maxGSSize = 2000,minGSSize = 10,
               pvalueCutoff = 1,pAdjustMethod='fdr',seed = FALSE,nPermSimple=1000) #有些通路为NA
  res <- list(gsea=gsea,tmp_plotdata=tmp_plotdata)
  return(res)
}

##### 2.1 每个肿瘤内的异质性 ##### 
pc_plotdata <- data.frame(x=numeric(),y=numeric(),sel=character(),celltype=character())
enrich_data_df <- data.frame(x=NULL,y=NULL,NES=NULL,PVAL=NULL)

for(s in samples){
  res  <- get_gsea(s,type,pbmc_metab)
  gsea_result <- as.data.frame(res$gsea)
  saveRDS(gsea_result,paste('../PC/',s,'_',type,'.RDS',sep=''))
  
  gsea_pathways <- rownames(gsea_result)
  enrich_data_df <- rbind(enrich_data_df,
                          data.frame(x=s,y=gsea_pathways,NES=gsea_result$NES,PVAL=gsea_result$p.adjust))
  
  tmp_plotdata <- res$tmp_plotdata
  pc_plotdata <- rbind(pc_plotdata,tmp_plotdata)
}
saveRDS(enrich_data_df,paste('../PC/enriched_pathway_',type,'.RDS',sep=''))
saveRDS(pc_plotdata,paste('../PC/pc_plotdata_',type,'.RDS',sep=''))

##### 2.2 肿瘤间异质性 ##### 
res <- get_gsea('allsample',type,pbmc_metab)
gsea_result <- res$gsea
saveRDS(gsea_result,paste('../PC/allsample_',type,'.RDS',sep=''))

pc_plotdata <- res$tmp_plotdata
saveRDS(pc_plotdata,paste('../PC/allsample_',type,'_pc_plot.RDS',sep=''))


#### 3. 可视化 ####
rm(list=ls())
library(ggplot2)
library(clusterProfiler)
library(scater)
library(stringr)
library(scran)
library(gtools)
options(stringsAsFactors=FALSE)
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

##### 3.1 肿瘤内异质性 #####
## 只挑选p<0.05的pathway可视化
# min_pval <- by(enrich_data_df$PVAL, enrich_data_df$y, FUN=min)
# select_pathways <- names(min_pval)[(min_pval<=0.05)] #
# print(select_pathways)
# plot_df <- enrich_data_df[enrich_data_df$y%in% select_pathways,]

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
ggsave(paste("../PC/Fig2C_enriched_pathway_",type,".pdf"),p,
       width = 5.5,height=3.3,units="in",device="pdf",useDingbats=FALSE)


## plot variance
pc_plotdata <- readRDS(paste('../PC/pc_plotdata_',type,'.RDS',sep=''))
pc_plotdata$celltype <- factor(pc_plotdata$celltype,levels=mixedsort(samples))
p <- ggplot(pc_plotdata) + geom_point(aes(x,y,colour=factor(sel)),size=0.5) +
  scale_color_manual(values=c("gray","#ff4000")) +
  facet_wrap(~celltype,scales="free",ncol = 4) + theme_bw() + 
  labs(x="Principal components", y="Explained variance (%)") +
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor= element_blank(),
        axis.line=element_line(size=0.2,colour="black"),
        axis.ticks = element_line(colour = "black",size=0.2),
        axis.text.x=element_text(colour="black", size = 6),
        axis.text.y=element_text(colour="black", size = 6),
        strip.background = element_rect(fill="white",size=0.2,colour = NULL),
        strip.text=element_text(size=6))

ggsave(paste("../PC/variance_plot_",type,".pdf",sep=""),p,
       device="pdf",useDingbats=FALSE)


##### 3.2 肿瘤间异质性 #####
type='ALL'
egmt <- readRDS(paste('../PC/allsample_',type,'.RDS',sep=''))
plot_df <- as.data.frame(egmt)

pdf(paste("../PC/Fig3_allsample_",type,".pdf",sep=''))
for(i in 1:7){
  p <- enrichplot::gseaplot2(egmt,geneSetID = i,title=plot_df$ID[i],pvalue_table = TRUE,color = '#E31A1C')
  print(p)
}
dev.off()
