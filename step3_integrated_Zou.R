#### 整合nc、GSE、zou数据
rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)

#### step1. 给zou数据匹配细胞类型，merge好的nc和GSE数据当成reference ####
pbmc_ref <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/pbmc.RDS')
DefaultAssay(pbmc_ref) <- 'RNA'

pbmc_query <- readRDS('/remote-home/yanzijun/CRU/scTALL/res/Seurat/RDS/pbmc_harmony.RDS')
DefaultAssay(pbmc_query) <- 'RNA'

## return.model = TRUE for MapQuery
pbmc_ref <-  RunUMAP(pbmc_ref,dims=1:30, reduction = "harmony",seed.use = 123,n.components=2,return.model = TRUE)
pbmc_query <- RunUMAP(pbmc_query,dims=1:30, reduction = "harmony", seed.use = 123,n.components=2,return.model = TRUE)

##
anchors <- FindTransferAnchors(reference = pbmc_ref, query = pbmc_query, dims = 1:30,reference.reduction = "pca")

predictions_1 <- TransferData(anchorset = anchors, refdata = pbmc_ref$celltype, dims = 1:30,k.weight = 30)
pbmc_query <- AddMetaData(pbmc_query, metadata = predictions_1[,c('predicted.id'),drop=FALSE])
colnames(pbmc_query@meta.data)[colnames(pbmc_query@meta.data)=='predicted.id'] <- 'celltype'

predictions_2 <- TransferData(anchorset = anchors, refdata = pbmc_ref$celltype_2, dims = 1:30,k.weight = 30)
pbmc_query <- AddMetaData(pbmc_query, metadata = predictions_2[,c('predicted.id'),drop=FALSE])
colnames(pbmc_query@meta.data)[colnames(pbmc_query@meta.data)=='predicted.id'] <- 'celltype_2'

pbmc_query$celltype_2[pbmc_query$celltype=='B/Mono cells' & pbmc_query$celltype_2 %in% c('leukemia cells','NK cells')] <- 'Monocytes'

pdf('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/sFig_pbmc_Zou_transferLabel.pdf',width = 14,height = 6)
DimPlot(pbmc_query,group.by='celltype')+DimPlot(pbmc_query,group.by='celltype_2')
dev.off()

### 与reference数据集拥有共同的UMAP图
pbmc_query <- MapQuery(anchorset = anchors, reference = pbmc_ref, query = pbmc_query,
                       refdata = list(celltype = "celltype_2"), reference.reduction = "pca", 
                       reduction.model = "umap")

p1 <- DimPlot(pbmc_ref, reduction = "umap", group.by = "celltype_2", label = TRUE, label.size = 3,
              repel = TRUE)  + ggtitle("Reference annotations")
p2 <- DimPlot(pbmc_query, reduction = "ref.umap", group.by = "celltype_2", label = TRUE,
              label.size = 3, repel = TRUE) + ggtitle("Query transferred labels")
pdf('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/sFig_pbmc_harmony_transferLable.pdf',width = 12,height = 5)
p1 + p2
dev.off()

pbmc_query$cluster <- pbmc_query$celltype
pbmc_query$study <- 'Zou'
pbmc_query$tissue <- 'TALL'
pbmc_query$type <- 'ALL'

saveRDS(pbmc_query,'/remote-home/yanzijun/CRU/scTALL/res/Seurat/RDS/pbmc_transfer.RDS')


####
#### step2.整合zou数据到GSE132509/nc2021中####
rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(harmony)
library(RColorBrewer)

## 1. 输入数据
pbmc_ref <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/pbmc.RDS')
pbmc_query <- readRDS('/remote-home/yanzijun/CRU/scTALL/res/Seurat/RDS/pbmc_transfer.RDS')

ol <- intersect(colnames(pbmc_query@meta.data),colnames(pbmc_ref@meta.data))

pbmc_ref@meta.data <- dplyr::select(pbmc_ref@meta.data,ol)
pbmc_query@meta.data <- dplyr::select(pbmc_query@meta.data,ol)
print(all(colnames(pbmc_query@meta.data)==colnames(pbmc_ref@meta.data)))


## 2. 合并数据
pbmc <- merge(x=pbmc_ref,y=list(pbmc_query))
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)
pbmc <- RunTSNE(pbmc, dims = 1:30)
pbmc <- RunUMAP(pbmc, dims = 1:30)

pbmc <- FindNeighbors(pbmc,dims = 1:30)
pbmc <- FindClusters(pbmc)

## 3. harmony
pbmc_harmony <- RunHarmony(pbmc,group.by.vars = 'study')
# pbmc_harmony <- RunTSNE(pbmc_harmony, dims = 1:30,reduction = "harmony")
# pbmc_harmony <- RunUMAP(pbmc_harmony,dims=1:30, reduction = "harmony", seed.use = 123,n.components=2)
pbmc_harmony <- RunUMAP(pbmc_harmony,dims = 1:50,reduction = "harmony",n.neighbors = 50,min.dist = 0.5,seed.use = 123) #10,0.5
pbmc_harmony <- RunTSNE(pbmc_harmony,dims = 1:50,reduction = "harmony",seed.use = 123) 
pbmc_harmony <- FindNeighbors(pbmc_harmony,dims = 1:30)
pbmc_harmony <- FindClusters(pbmc_harmony,resolution = 0.1)


#### step3. 重新匹配细胞类型: 1. Mono和pDC合并/2.高表达NK marker的为NK，其他都为T cell ####
all(colnames(pbmc)==colnames(pbmc_harmony))
pbmc_harmony$celltype_2 <- pbmc$celltype_2
table(pbmc_harmony$celltype_2[pbmc_harmony$celltype=='B/Mono cells' & pbmc_harmony$celltype_2 %in% c('leukemia cells','NK cells')])
pbmc_harmony$celltype_2[pbmc_harmony$celltype=='B/Mono cells' & pbmc_harmony$celltype_2 %in% c('leukemia cells','NK cells')] <- 'Monocytes'

pbmc$celltype_2[pbmc$celltype=='B/Mono cells' & pbmc$celltype_2 %in% c('leukemia cells','NK cells')] <- 'Monocytes'

## 
Idents(pbmc_harmony) <- 'celltype_2'
NK.pbmc_harmony <- subset(pbmc_harmony,idents='NK cells')
NK.pbmc_harmony <- RunUMAP(NK.pbmc_harmony,dims=1:30, reduction = "harmony", seed.use = 123,n.components=2)
NK.pbmc_harmony <- FindNeighbors(NK.pbmc_harmony,dims = 1:30)
NK.pbmc_harmony <- FindClusters(NK.pbmc_harmony,resolution = 0.1)
NKcells <- colnames(NK.pbmc_harmony)[NK.pbmc_harmony$seurat_clusters=='2']
Tcells <- colnames(NK.pbmc_harmony)[NK.pbmc_harmony$seurat_clusters !='2']

celltype_3 <-  pbmc_harmony$celltype_2
celltype_3[celltype_3=='pDCs'] <- 'Monocytes'
celltype_3[colnames(pbmc_harmony) %in% NKcells] <- 'NK cells'
celltype_3[colnames(pbmc_harmony) %in% Tcells] <- 'T cells'

pbmc_harmony$celltype_3 <- celltype_3
pbmc_harmony$celltype_abbr <- ''
pbmc_harmony$celltype_abbr[pbmc_harmony$celltype_3=='leukemia cells'] <- 'ALL'
pbmc_harmony$celltype_abbr[pbmc_harmony$celltype_3=='NK cells'] <- 'NK'
pbmc_harmony$celltype_abbr[pbmc_harmony$celltype_3=='Erythrocytes'] <- 'Ery.'
pbmc_harmony$celltype_abbr[pbmc_harmony$celltype_3=='B cells'] <- 'B'
pbmc_harmony$celltype_abbr[pbmc_harmony$celltype_3=='Monocytes'] <- 'Mono.'
pbmc_harmony$celltype_abbr[pbmc_harmony$celltype_3=='T cells'] <- 'T'
pbmc_harmony$celltype_abbr <- factor(pbmc_harmony$celltype_abbr,levels = c('ALL','B','T','Mono.','NK','Ery.'))

saveRDS(pbmc,'/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/pbmc_batch.RDS')
saveRDS(pbmc_harmony,'/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/pbmc.RDS')

meta <- pbmc@meta.data
saveRDS(meta,'/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/meta_batch.RDS')

meta_harmony <- pbmc_harmony@meta.data
saveRDS(meta_harmony,'/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/meta.RDS')


#### step5. Find MK ####
library(future)
plan(multisession, workers=40)

pbmc_harmony <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/pbmc.RDS')
print(table(pbmc_harmony$celltype_abbr))

Idents(pbmc_harmony) <- 'celltype_abbr'
MK <- FindAllMarkers(pbmc_harmony,only.pos = T)
saveRDS(MK,'/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/MK_celltype.RDS')

print(table(MK$cluster)) #no ALL mk
