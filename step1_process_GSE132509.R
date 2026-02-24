#######step1.整理单细胞公共数据：cALL and PBMC from GSE132509
rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
setwd('/remote-home/yanzijun/CRU/ped_M/data/GSE132509/')
meta <- read.csv(file = '/remote-home/yanzijun/CRU/ped_M/data/GSE132509/GSE132509_cell_annotations.tsv',sep='\t')
print(table(meta$orig.ident,meta$celltype)) ##PBMMC.1/3里分别有8/2tumor cells,删掉
meta <- meta[!(meta$orig.ident %in% c('PBMMC.1','PBMMC.3') & meta$celltype %in% c('ETV6.RUNX1.1','ETV6.RUNX1.4','HHD.1')),]

PBMC1 <- Read10X('PBMC1/')
ctr1=CreateSeuratObject(counts = PBMC1, min.cells = 0, min.features = 0, project = 'PBMC1')

PBMC2 <- Read10X('PBMC2/')
ctr2=CreateSeuratObject(counts = PBMC2, min.cells = 0, min.features = 0, project = 'PBMC2')

PBMC3 <- Read10X('PBMC3/')
ctr3=CreateSeuratObject(counts = PBMC3, min.cells = 0, min.features = 0, project = 'PBMC3')

BALL1 <- Read10X('BALL1/')
pbmc1=CreateSeuratObject(counts = BALL1, min.cells = 0, min.features = 0, project = 'ETV6.RUNX1.1')

BALL2 <- Read10X('BALL2/')
pbmc2=CreateSeuratObject(counts = BALL2, min.cells = 0, min.features = 0, project = 'ETV6.RUNX1.2')

BALL3 <- Read10X('BALL3/')
pbmc3=CreateSeuratObject(counts = BALL3, min.cells = 0, min.features = 0, project = 'ETV6.RUNX1.3')

BALL4 <- Read10X('BALL4/')
pbmc4=CreateSeuratObject(counts = BALL4, min.cells = 0, min.features = 0, project = 'ETV6.RUNX1.4')

BALL5 <- Read10X('BALL5/')
pbmc5=CreateSeuratObject(counts = BALL5, min.cells = 0, min.features = 0, project = 'HHD.1')

BALL6 <- Read10X('BALL6/')
pbmc6=CreateSeuratObject(counts = BALL6, min.cells = 0, min.features = 0, project = 'HHD.2')

TALL1 <- Read10X('TALL1/')
pbmc7=CreateSeuratObject(counts = TALL1, min.cells = 0, min.features = 0, project = 'PRE-T.1')

TALL2 <- Read10X('TALL2/')
pbmc8=CreateSeuratObject(counts = TALL2, min.cells = 0, min.features = 0, project = 'PRE-T.2')

tmp1 <- merge(ctr1, y = list(ctr2,ctr3,pbmc1,pbmc2,pbmc3,pbmc4,pbmc5,pbmc6,pbmc7,pbmc8),
             add.cell.ids = c('PBMMC.1','PBMMC.2','PBMMC.3',"ETV6.RUNX1.1", "ETV6.RUNX1.2", "ETV6.RUNX1.3", "ETV6.RUNX1.4",
                              "HHD.1" ,"HHD.2","PRE-T.1", "PRE-T.2"), project = "cALL")

new.name <- c(Cells(tmp1)[grep('PBMMC.1',Cells(tmp1))],
              sub("-[0-9]+$", "", Cells(tmp1)[-grep('PBMMC.1',Cells(tmp1))])) ### 细胞删除后面-[0-9]
tmp2 <- RenameCells(tmp1,new.names = new.name) 

pbmc_QC <- subset(tmp2,cells = meta$cell_id) ### only QC pass cALL data
print(all(colnames(pbmc_QC)==meta$cell_id))

## add meta
rownames(meta) <- meta$cell_id
meta$cluster <- 'leukemia cells'
meta$cluster[meta$celltype=='B cells + Mono'] <- 'B/Mono cells'
meta$cluster[meta$celltype=='Erythrocytes'] <- 'Erythrocytes'
meta$cluster[meta$celltype=='T cells + NK'] <- 'T/NK cells'

pbmc_QC@meta.data <- meta
head(pbmc_QC@meta.data)
saveRDS(pbmc_QC,'/remote-home/yanzijun/CRU/ped_M/data/GSE132509/pbmc_QC.RDS')

####### process
pbmc_batch <- NormalizeData(object = pbmc_QC, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_batch <- FindVariableFeatures(object = pbmc_batch, selection.method = "vst", nfeatures = 4000)
pbmc_batch <- ScaleData(object = pbmc_batch, features = VariableFeatures(object = pbmc_batch))
pbmc_batch <- RunPCA(object = pbmc_batch, seed.use=123, npcs=150,
                     features = VariableFeatures(object = pbmc_batch), ndims.print=1,nfeatures.print=1)
pbmc_batch <- RunTSNE(pbmc_batch, dims = 1:50, seed.use = 123,n.components=2)
pbmc_batch <- RunUMAP(pbmc_batch, dims = 1:50, seed.use = 123,n.components=2)

pdf('pbmc_batch.pdf')
DimPlot(pbmc_batch,group.by = 'celltype')
DimPlot(pbmc_batch,group.by = 'cluster')
ggplot(meta, aes(x = UMAP1, y = UMAP2, col = celltype)) + geom_point() #from article
ggplot(meta, aes(x = UMAP1, y = UMAP2, col = cluster)) + geom_point()#from article
dev.off()

## re-set TSNE/UMAP coordinates
pbmc_batch@reductions$tsne@cell.embeddings <- as.matrix(meta[,c('tSNE_1','tSNE_2')])

pbmc_batch@reductions$umap@cell.embeddings <- as.matrix(meta[,c('UMAP1','UMAP2')])
colnames(pbmc_batch@reductions$umap@cell.embeddings) <- c('UMAP_1','UMAP_2')

TSNEPlot(pbmc_batch,group.by = 'cluster')
dev.off()

## split into cALL
ALL.cells <- Cells(pbmc_batch)[-grep('^PBMMC',Cells(pbmc_batch))]
pbmc_ALL <- subset(pbmc_batch,cells=ALL.cells)

pdf('pbmc_ALL.pdf')
DimPlot(pbmc_ALL,group.by = 'celltype')
DimPlot(pbmc_ALL,group.by = 'cluster')
dev.off()

saveRDS(pbmc_batch,'/remote-home/yanzijun/CRU/ped_M/data/GSE132509/pbmc_batch.RDS')
saveRDS(pbmc_ALL,'/remote-home/yanzijun/CRU/ped_M/data/GSE132509/pbmc_ALL.RDS')

