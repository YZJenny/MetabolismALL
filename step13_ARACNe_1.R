######### 跑ARACEN：制作输入文件
rm(list=ls())
library(corto)
library(Seurat)
library(clusterProfiler)
library(FactoMineR)
library(ggplot2)
source('/remote-home/yanzijun/CRU/ped_M/scr/cortoplot_myself_v3.R')
source('/remote-home/yanzijun/software/InferLoop-main/inferloop/InferLoop.R')

setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/ARACNe/')
print(getwd())

#### 1. merge cells into 100 bin based on PCA ####
generateBin <- function(indata, used_coords, n=100, seed=123, type=NULL){
  DATA=indata
  used_coords=used_coords
  if(is.null(type)==FALSE){
    type_coords=as.numeric(as.factor(type))*10
    used_coords=cbind(used_coords, type_coords)
  }
  n=n
  set.seed(seed)
  KM=kmeans(used_coords,centers=n)
  CLST=KM$cluster
  names(CLST)=colnames(indata)
  OUT=.generate_mean(DATA, CLST)
  ##########################
  RETURN=list()
  RETURN$mat=OUT
  RETURN$clst=CLST
  return(RETURN)
}

indata <- as.matrix(readRDS('../TCAsubtype_new/exp_tumor_TCAcycle.RDS'))
inmeta <- readRDS('../TCAsubtype_new/meta_tumor_TCAcycle.RDS')

df.pca <- PCA(t(indata), graph = FALSE)
used_coords=as.data.frame(df.pca$ind$coord[,1:2])

print(all(rownames(used_coords)==inmeta$sampleID))
used_coords$subtype <- inmeta$subtype_new

# ggplot(used_coords,aes(x=Dim.1,y=Dim.2,color=subtype))+geom_point()
# dev.off()

indata_High <- indata[,colnames(indata) %in% inmeta$sampleID[inmeta$subtype_new=='sc_High']]
indata_Low <- indata[,colnames(indata) %in% inmeta$sampleID[inmeta$subtype_new=='sc_Low']]

used_coords_High <- used_coords[used_coords$subtype=='sc_High',]
used_coords_Low <- used_coords[used_coords$subtype=='sc_Low',]

print(all(colnames(indata_High)==rownames(used_coords_High)))
print(all(colnames(indata_Low)==rownames(used_coords_Low)))

used_coords_High$subtype=used_coords_Low$subtype=NULL

BIN_High=generateBin(indata_High,used_coords_High, n=500)
BIN_Low=generateBin(indata_Low,used_coords_Low, n=500)

mat_High <- BIN_High$mat
colnames(mat_High) <- paste('High_',colnames(mat_High),sep='')
mat_Low <- BIN_Low$mat
colnames(mat_Low) <- paste('Low_',colnames(mat_Low),sep='')
print(all(rownames(mat_High)==rownames(mat_Low)))
mat <- as.data.frame(cbind(mat_High,mat_Low))
mat <- tibble::rownames_to_column(mat,'Gene')

write.table(mat,'mat.txt',row.names = F,col.names = T,sep='\t',quote = F)

saveRDS(BIN_High, file='BIN_High.rds')
saveRDS(BIN_Low, file='BIN_Low.rds')

