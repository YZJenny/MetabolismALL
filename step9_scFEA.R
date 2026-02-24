##### scFEA using node9
rm(list=ls())
library(Seurat)
library(dplyr)
#setwd('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/')
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')


#### 1. make input csv file ####
pbmc <- readRDS('pbmc.RDS')

## sample cells
Idents(pbmc) <- pbmc$celltype_2
sub.pbmc <- subset(pbmc,downsample=5000) 
data <- as.matrix(sub.pbmc@assays$RNA@data)
write.csv(data,'../scFEA/ALL_sc_sample5k.csv',row.names = T)
meta <- sub.pbmc@meta.data
write.csv(meta,'../scFEA/ALL_meta_sample5k.csv',row.names = T)


#### 2. move to node9 /local/yanzijun/software/scFEA-master/input_addZou ####
# cd ~/software/scFEA-master/
# 
# nohup python src/scFEA.py --data_dir data --input_dir input_addZou \
# --test_file ALL_sc_sample5k.csv \
# --moduleGene_file module_gene_m168.csv \
# --stoichiometry_matrix cmMat_c70_m168.csv \
# --res_dir output_addZou \
# --sc_imputation True &

  
#### 3. 可视化move to CRU ../scFEA/
rm(list=ls())
library(ggplot2)
library(Seurat)
library(dplyr)
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

## 3.1 input
scFEA <- read.csv('../scFEA/ALL_sc_sample5k_module168.csv',row.names = 1)
scFEA[1:3,1:4]
all.meta <- readRDS('meta.RDS')
all.meta[1:3,]

M2C <- read.csv("/remote-home/yanzijun/software/scFEA-master/data/Human_M168_information.symbols.csv",row.names = 1)

## ALL samples
ALLcells <- rownames(all.meta)[-grep('^PBMMC',all.meta$orig.ident)]
olcells <- intersect(ALLcells,rownames(scFEA))

meta <- all.meta[which(rownames(all.meta)%in%olcells),]

flux <- scFEA[rownames(scFEA) %in% rownames(meta),]

scale_rows<- function(x){ #行为基因，列为样本
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m)/s)
}

flux_norm <- t(scale_rows(t(as.matrix(flux))))

## 3.2 findmk
counts <- as.data.frame(t(flux))
info <- meta[,'celltype_abbr',drop=FALSE]
colnames(info) <- 'celltype'
info$celltype <- factor(info$celltype,levels = sort(unique(info$celltype)))
print(all(colnames(counts)==rownames(info)))

## 细胞类型特异性的pathway
obj <- CreateSeuratObject(counts = counts,meta.data = info)
Idents(obj)=obj$celltype
mk <- FindAllMarkers(obj,slot = "data",logfc.threshold = 0,only.pos = TRUE)
print(table(mk$cluster[mk$avg_log2FC>0.25]))
print(table(mk$cluster))
  
print(table(mk$cluster[mk$avg_log2FC>0.01]))
mk[mk$avg_log2FC>0.01,]



#### 可视化
### 3.3 pheatmap by celltype in ALL
df <- dplyr::select(flux,gsub('-','_',unique(mk$gene)))
df$celltype <- info$celltype
df$celltype <- factor(df$celltype,levels=sort(unique(df$celltype)))
print(levels(df$celltype))

# 根据类型分组，计算每个基因在每个类别中的平均表达值
plot.df <- aggregate(.~ celltype, data = df, mean)
rownames(plot.df) <- plot.df$celltype;plot.df$celltype=NULL
plot.df <- t(plot.df)

## 只show cycle的metabolite
tag <- rownames(plot.df)
tag[!rownames(plot.df) %in% paste('M_',5:14,sep='')] <- ''
rownames(plot.df) <- tag

color <- colorRampPalette(c("blue","white","red"))(100)
library(pheatmap)
p <- pheatmap(plot.df,scale = 'row',cluster_cols = F,cluster_rows = F,color=color)
pdf('../scFEA/pheatmap_DEM_celltype.pdf',
    width = 2,height = 6.6) #width = 2,height = 6.6
print(p) 
dev.off()


## 3.4 only M5-M14 boxplot
library(ggpubr)
library(ggsignif)
library(RColorBrewer)
col = colorRampPalette(brewer.pal(12, "Set3"))(12)

CT=levels(info$celltype)
mycol <- col[1:length(CT)]
names(mycol) <- CT

p.lst <- list()
for(M in paste('M',5:13,sep='_')){
  print(M)
  df <- data.frame(flux=as.numeric(flux_norm[,colnames(flux_norm)==M]),celltype=meta$celltype_abbr)
  print(kruskal.test(flux~celltype,data=df)$p.value)
  
  p <-  ggplot(df, aes(x=celltype, y=flux, fill=celltype)) +
    geom_boxplot(position = position_dodge(0.05),width=0.5,outlier.shape = NA) +
    scale_fill_manual(values =mycol) +
    theme_classic()+
    labs(y='Metabolic flux',title = M,x='')+
    theme(plot.title = element_text(hjust = 0.5,size = 8),
          axis.text = element_text(size=8,face="plain",color="black"),
          axis.title  = element_text(size=8,face="plain",color="black"),
          axis.text.x = element_text(hjust = 1,angle =90),
          legend.position="none")+
    coord_cartesian(ylim=c(-1,3))
    
  p.lst[[M]] <- p
}

p <- cowplot::plot_grid(plotlist = p.lst,ncol = 5)
ggsave(filename = '../scFEA/boxplot_TCA_celltype.pdf',plot = p,
       width = 10,height = 4)


## 3.5 Fig2D-F only M5-M14 pheatmap by celltype in ALL
Module <- paste('M',5:14,sep='_')
df <- flux[,colnames(flux) %in% Module]
df$celltype <- meta$celltype_abbr

plot_df <- aggregate(. ~ celltype,data=df,FUN = median)
rownames(plot_df) <- plot_df$celltype;plot_df$celltype=NULL

p1 <- pheatmap(t(plot_df),scale = 'row',cluster_cols = T,cluster_rows = F,
               color = colorRampPalette(c("grey","white","red"))(100))
ggsave(filename = '../scFEA/Fig2D_heatmap_TCA_celltype.pdf',plot = p1,
       width = 3,height = 4)


## 3.6 Fig2D-F pheatmap by batch in ALL cells
Module <- paste('M',5:14,sep='_')
df <- flux[,colnames(flux) %in% Module]
print(all(rownames(df)==rownames(meta)))
df$orig.ident <- meta$orig.ident
df <- df[ meta$celltype_abbr=='ALL',]

plot_df <- aggregate(. ~ orig.ident,data=df,FUN = median)
rownames(plot_df) <- plot_df$orig.ident;plot_df$orig.ident=NULL

p1 <- pheatmap(t(plot_df),scale = 'row',cluster_cols = T,cluster_rows = F,
               color = colorRampPalette(c("grey","white","red"))(100),clustering_method = "ward.D")
ggsave(filename = '../scFEA/Fig2E_heatmap_TCA_bybatch.pdf',plot = p1,
       width = 5,height = 4)


## 3.7 pheatmap by batch in celltype in BMMC
Module <- paste('M',5:14,sep='_')
df <- scFEA[,colnames(scFEA) %in% Module]
df$celltype <- all.meta$celltype_abbr
df <- df[grep('PBMMC',all.meta$orig.ident),]

plot_df <- aggregate(. ~ celltype,data=df,FUN = median)
rownames(plot_df) <- plot_df$celltype;plot_df$celltype=NULL

p1 <- pheatmap(t(plot_df),scale = 'row',cluster_cols = T,cluster_rows = F,
               color = colorRampPalette(c("grey","white","red"))(100))
ggsave(filename = '../scFEA/heatmap_TCA_BMMC.pdf',plot = p1,
       width = 2.5,height = 4)
