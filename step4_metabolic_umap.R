### single-cell landscape of metabolic genes
rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(harmony)
library(RColorBrewer)
library(clusterProfiler)
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

pathway2gene <- read.gmt('/remote-home/yanzijun/CRU/TALL_M/data/GeneSet7.gmt')
genelst <- unique(pathway2gene$gene)
print(length(genelst)) #1804

pbmc_harmony <- readRDS('pbmc.RDS')

#### 1.all cells
DefaultAssay(pbmc_harmony) <- 'RNA'
pbmc_metab <- ScaleData(object = pbmc_harmony, features = genelst)
pbmc_metab <- RunPCA(pbmc_harmony,features =genelst )
pbmc_metab <- RunTSNE(pbmc_metab, dims = 1:30,check_duplicates = FALSE)
pbmc_metab <- RunUMAP(pbmc_metab, dims = 1:30)
saveRDS(pbmc_metab,'pbmc_metab.RDS')


#### 2. separate leukemia and normal cells
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