### metabolic genes的single-cell landscape
rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(harmony)
library(RColorBrewer)
library(clusterProfiler)
#setwd('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/')
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

pathway2gene <- read.gmt('/remote-home/yanzijun/CRU/TALL_M/data/GeneSet7.gmt')
genelst <- unique(pathway2gene$gene)
print(length(genelst)) #1804

pbmc_harmony <- readRDS('pbmc.RDS')

#### 1.所有细胞
DefaultAssay(pbmc_harmony) <- 'RNA'
pbmc_metab <- ScaleData(object = pbmc_harmony, features = genelst)
pbmc_metab <- RunPCA(pbmc_harmony,features =genelst )
pbmc_metab <- RunTSNE(pbmc_metab, dims = 1:30,check_duplicates = FALSE)
pbmc_metab <- RunUMAP(pbmc_metab, dims = 1:30)
saveRDS(pbmc_metab,'pbmc_metab.RDS')


#### 2. 肿瘤和正常细胞分开
DefaultAssay(pbmc_harmony) <- 'RNA'
tumorcell <- Cells(pbmc_harmony)[pbmc_harmony$celltype_2=='leukemia cells']
pbmc_tumor <- subset(pbmc_harmony,cells = tumorcell)
print(table(pbmc_tumor$celltype_2))

normalcell <- Cells(pbmc_harmony)[pbmc_harmony$celltype_2 !='leukemia cells']
pbmc_normal <- subset(pbmc_harmony,cells = normalcell)
print(table(pbmc_normal$celltype_2))

print(table(pbmc_normal$orig.ident[pbmc_normal$celltype_2=='leukemia cells']))


pbmc_tumor <- ScaleData(object = pbmc_tumor, features = genelst)
pbmc_tumor <- RunPCA(pbmc_tumor,features =genelst)
pbmc_tumor <- RunTSNE(pbmc_tumor, dims = 1:30,check_duplicates = FALSE)
pbmc_tumor <- RunUMAP(pbmc_tumor, dims = 1:30)

pbmc_normal <- ScaleData(object = pbmc_normal, features = genelst)
pbmc_normal <- RunPCA(pbmc_normal,features =genelst)
pbmc_normal <- RunTSNE(pbmc_normal, dims = 1:30,check_duplicates = FALSE)
pbmc_normal <- RunUMAP(pbmc_normal, dims = 1:30)

saveRDS(pbmc_tumor,'pbmc_tumor_metab.RDS')
saveRDS(pbmc_normal,'pbmc_normal_metab.RDS')


#### 4.可视化
#### 4.1 所有细胞
pdf('TSNE_harmony_metab.pdf',width = 18,height = 10)
  TSNEPlot(pbmc_metab, group.by = "celltype")+
  TSNEPlot(pbmc_metab, group.by = "celltype_2")+
  TSNEPlot(pbmc_metab, group.by = "tissue")+
  TSNEPlot(pbmc_metab, group.by = "study")+
  TSNEPlot(pbmc_metab, group.by = "type")
dev.off()

pdf('UMAP_harmony_metab.pdf',width = 18,height = 10)
  UMAPPlot(pbmc_metab, group.by = "celltype")+
  UMAPPlot(pbmc_metab, group.by = "celltype_2")+
  UMAPPlot(pbmc_metab, group.by = "tissue")+
  UMAPPlot(pbmc_metab, group.by = "study")+
  UMAPPlot(pbmc_metab, group.by = "type")
dev.off()


#### 4.2tumor和正常细胞分开画
## tumor cells
pdf('TSNE_harmony_tumor_metab.pdf',width = 18,height =5)
TSNEPlot(pbmc_tumor, group.by = "orig.ident")+
  TSNEPlot(pbmc_tumor, group.by = "tissue")+
  TSNEPlot(pbmc_tumor, group.by = "study")
dev.off()

pdf('UMAP_harmony_tumor_metab.pdf',width = 18,height = 5)
UMAPPlot(pbmc_tumor, group.by = "orig.ident")+
  UMAPPlot(pbmc_tumor, group.by = "tissue")+
  UMAPPlot(pbmc_tumor, group.by = "study")
dev.off()


## normal cells
pdf('TSNE_harmony_normal_metab.pdf',width = 18,height = 10)
TSNEPlot(pbmc_normal, group.by = "orig.ident")+
  TSNEPlot(pbmc_normal, group.by = "celltype")+
  TSNEPlot(pbmc_normal, group.by = "celltype_2")+
  TSNEPlot(pbmc_normal, group.by = "tissue")+
  TSNEPlot(pbmc_normal, group.by = "study")+
  TSNEPlot(pbmc_normal, group.by = "type")
dev.off()

pdf('UMAP_harmony_normal_metab.pdf',width = 18,height = 10)
UMAPPlot(pbmc_normal, group.by = "orig.ident")+
  UMAPPlot(pbmc_normal, group.by = "celltype")+
  UMAPPlot(pbmc_normal, group.by = "celltype_2")+
  UMAPPlot(pbmc_normal, group.by = "tissue")+
  UMAPPlot(pbmc_normal, group.by = "study")+
  UMAPPlot(pbmc_normal, group.by = "type")
dev.off()
