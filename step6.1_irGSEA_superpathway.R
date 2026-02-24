####2. 细胞类型特异性的代谢通路活性(7个superpathway)
rm(list=ls())
library(irGSEA)
library(Seurat)
library(ggplot2)
library(dplyr)
#setwd('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/')
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')

pathway2gene <- read.gmt('/remote-home/yanzijun/CRU/TALL_M/data/GeneSet7.gmt')
genelst <- split(pathway2gene$gene,pathway2gene$term)
print(length(genelst))
metabGene <- list(metabGene=unique(pathway2gene$gene))
genelst  <- c(genelst,metabGene)
print(length(genelst))

pbmc <- readRDS('pbmc.RDS')
meta <- readRDS('meta.RDS')
pbmc@meta.data <- meta


## irGSEA中没有addModuleScore
tmp <- AddModuleScore(pbmc,features = genelst,name = names(genelst))
addMS <- tmp@meta.data[,(ncol(tmp@meta.data) - length(genelst)+1) : ncol(tmp@meta.data)]
colnames(addMS) <- gsub('[0-9]','',colnames(addMS))

obj <- irGSEA.score(object = pbmc, assay = "RNA", 
                          slot = "data", seeds = 123, ncores = 40,
                          min.cells = 3, min.feature = 0,
                          maxGSSize=2000,
                          custom = T, geneset = genelst, msigdb = F, 
                          species = "Homo sapiens", 
                          subcategory = NULL, 
                          method = c("AUCell","ssgsea"),
                          aucell.MaxRank = NULL, ucell.MaxRank = 2000, 
                          kcdf = 'Gaussian')
## 保存结果
assays <- obj@assays
assays$RNA=NULL
assays$SCT=NULL

addMS.obj <- CreateAssayObject(counts = t(addMS))
assays$addMS <- addMS.obj

meta <- obj@meta.data
irGSEA <- list(assays=assays,meta=meta)
saveRDS(irGSEA,'../irGSEA/pbmc_all_superpathway.RDS')

