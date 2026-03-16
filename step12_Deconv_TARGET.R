#### Deconv bulk RNA-seq from TARGET cohort
rm(list=ls())
library(MuSiC)
library(Biobase)
library(Seurat)
library(xbioc)
library(BisqueRNA)
library(ggplot2)
library(SummarizedExperiment)
library(RColorBrewer)
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/')

#### 1.load data ####
## 1.1 load sc data
pbmc <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/pbmc.RDS')
meta <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/meta.RDS')
pbmc@meta.data <- meta
subtype <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/TCAsubtype_new/subtype_info.RDS')
low_cells <- rownames(subtype)[subtype$subtype=='Low'] 
Middle_cells <- rownames(subtype)[subtype$subtype=='Middle'] 
high_cells <- rownames(subtype)[subtype$subtype=='High']

pbmc@meta.data$subtype <- as.character(pbmc$celltype_abbr)
pbmc@meta.data$subtype[which(colnames(pbmc) %in% high_cells)] <- 'High'
pbmc@meta.data$subtype[which(colnames(pbmc) %in% Middle_cells)] <- 'Middle'
pbmc@meta.data$subtype[which(colnames(pbmc) %in% low_cells)] <- 'Low'
print(table(pbmc@meta.data$subtype))

Idents(pbmc) <- 'type'
pbmc_ALL <- subset(pbmc,idents='ALL')


## 1.2 load bulk data
## ALL
df <- readRDS('/remote-home/yanzijun/CRU/ped_M/res/RNAseq/exp_ALL.count.rds')
bulk.mtx_ALL <- as.matrix(na.omit(df))


#### 2. Bisque ####
get_Bisque <- function(bulk.mtx,subpbmc){
  featureData <- data.frame(featureNames=rownames(bulk.mtx))
  rownames(featureData) <- rownames(bulk.mtx)
  phenoData <- data.frame(sampleNames =colnames(bulk.mtx))
  rownames(phenoData) <- colnames(bulk.mtx)
  
  bulk.eset <- ExpressionSet(assayData = bulk.mtx,
                             featureData = new("AnnotatedDataFrame", data = featureData),
                             phenoData = new("AnnotatedDataFrame", data = phenoData))
  
  featureData <- data.frame(featureNames=rownames(subpbmc))
  rownames(featureData) <- rownames(subpbmc)
  phenoData <- data.frame(SubjectName =subpbmc$orig.ident,cellType=subpbmc$subtype)
  rownames(phenoData) <- colnames(subpbmc)
  
  sc.eset <- ExpressionSet(assayData = as.matrix(subpbmc@assays$RNA@counts),
                           featureData = new("AnnotatedDataFrame", data = featureData),
                           phenoData = new("AnnotatedDataFrame", data = phenoData))
  
  Bisque <- ReferenceBasedDecomposition(bulk.eset = bulk.eset,sc.eset = sc.eset,
                                        markers = intersect(rownames(sc.eset),rownames(bulk.mtx)),
                                        cell.types = 'cellType',
                                        subject.names = 'SubjectName',use.overlap = FALSE)
  return(Bisque)
}

Bisque_ALL <- get_Bisque(bulk.mtx = bulk.mtx_ALL,subpbmc = pbmc_ALL)
saveRDS(Bisque_ALL,'Bisque/TARGET_ALL_Bisque.rds')