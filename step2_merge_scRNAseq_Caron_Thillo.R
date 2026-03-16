#######step2. merge two published scRNA-seq data: Caron et al. (GSE132509) and Thillo et al.(nc2021)
rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(harmony)
library(RColorBrewer)

## 1. input data
pbmc_GSE <- readRDS('/remote-home/yanzijun/CRU/ped_M/data/GSE132509/pbmc_batch.RDS')
pbmc_nc <- readRDS('/remote-home/yanzijun/CRU/scTALL/public/NC_2021/single_cell_rnaseq_input/tall_filtered_merged_umap.Rds')

pbmc_GSE@meta.data <- pbmc_GSE@meta.data[,c('orig.ident','cluster')]
pbmc_GSE@meta.data$study <- 'sciRep'

pbmc_nc@meta.data <- pbmc_nc@meta.data[,c('orig.ident','class')]
pbmc_nc@meta.data$study <- 'nc'

colnames(pbmc_nc@meta.data) <- colnames(pbmc_GSE@meta.data)

pbmc_nc$cluster <- as.character(pbmc_nc$cluster)
pbmc_nc$cluster[as.character(pbmc_nc$cluster)=='Malignant_T_Cells'] <- 'leukemia cells'

## 2. merge data 
pbmc <- merge(x=pbmc_GSE,y=list(pbmc_nc))
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)
pbmc <- RunTSNE(pbmc, dims = 1:30)
pbmc <- RunUMAP(pbmc, dims = 1:30)

pbmc <- FindNeighbors(pbmc,dims = 1:30)
pbmc <- FindClusters(pbmc)

pbmc$celltype <- pbmc$cluster
pbmc$celltype[pbmc$cluster %in% c('B_Cells','Monocytes')] <- 'B/Mono cells'
pbmc$celltype[pbmc$cluster %in% c('NK_Cells','T_Cells')] <- 'T/NK cells'
table(pbmc$celltype)

pbmc$tissue <- 'TALL'
pbmc$tissue[grep('^PBMMC',pbmc$orig.ident)] <- 'PBMMC'
pbmc$tissue[grep('^ETV|HHD',pbmc$orig.ident)] <- 'BALL'
table(pbmc$tissue)

pbmc$type <- 'ALL'
pbmc$type[pbmc$tissue=='PBMMC'] <- 'BMMC'
table(pbmc$type)

## 3. harmony
pbmc_harmony <- RunHarmony(pbmc,group.by.vars = 'study')
pbmc_harmony <- RunTSNE(pbmc_harmony, dims = 1:30,reduction = "harmony")
pbmc_harmony <- RunUMAP(pbmc_harmony,dims=1:30, reduction = "harmony", seed.use = 123,n.components=2)
pbmc_harmony <- FindNeighbors(pbmc_harmony,dims = 1:30)
pbmc_harmony <- FindClusters(pbmc_harmony)

## 4. re-name cluster
pbmc_harmony$celltype_2 <- pbmc_harmony$celltype

ori_tnse  <- as.data.frame(Embeddings(pbmc_harmony,reduction = 'tsne'))
ori_tnse$celltype <- pbmc_harmony$celltype
head(ori_tnse)

Monocytes <- rownames(ori_tnse)[ori_tnse$tSNE_1 < -20 & ori_tnse$tSNE_2 > -20 &  ori_tnse$celltype=='B/Mono cells']
pbmc_harmony$celltype_2[Cells(pbmc_harmony) %in% Monocytes]='Monocytes'

Bcells <-  rownames(ori_tnse)[ori_tnse$tSNE_1> -20 & ori_tnse$tSNE_2 < -20 &  ori_tnse$celltype=='B/Mono cells']
pbmc_harmony$celltype_2[Cells(pbmc_harmony) %in% Bcells]='B cells'

pbmc_harmony$celltype_2[pbmc_harmony$celltype_2=='T/NK cells']='NK cells'
pbmc_harmony$celltype_2[pbmc_harmony$celltype_2=='B/Mono cells']='B cells'

## 5. save
saveRDS(pbmc,'/remote-home/yanzijun/CRU/Metabolism_ALL/res_ALL/Seurat/pbmc_GSE_nc_batch.RDS')
saveRDS(pbmc_harmony,'/remote-home/yanzijun/CRU/Metabolism_ALL/res_ALL/Seurat/pbmc_GSE_nc.RDS')
