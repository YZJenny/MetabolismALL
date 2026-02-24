#######step2.合并单细胞公共数据：GSE132509/nc2021
rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(harmony)
library(RColorBrewer)

## 1. 输入数据
pbmc_GSE <- readRDS('/remote-home/yanzijun/CRU/ped_M/data/GSE132509/pbmc_batch.RDS')
pbmc_nc <- readRDS('/remote-home/yanzijun/CRU/scTALL/public/NC_2021/single_cell_rnaseq_input/tall_filtered_merged_umap.Rds')

pbmc_GSE@meta.data <- pbmc_GSE@meta.data[,c('orig.ident','cluster')]
pbmc_GSE@meta.data$study <- 'sciRep'

pbmc_nc@meta.data <- pbmc_nc@meta.data[,c('orig.ident','class')]
pbmc_nc@meta.data$study <- 'nc'

colnames(pbmc_nc@meta.data) <- colnames(pbmc_GSE@meta.data)

pbmc_nc$cluster <- as.character(pbmc_nc$cluster)
pbmc_nc$cluster[as.character(pbmc_nc$cluster)=='Malignant_T_Cells'] <- 'leukemia cells'

## 2. 合并数据
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


### 5. 可视化
Color1 <- c("#A6CEE3","#1F78B4","#CAB2D6" ,"#33A02C","#E31A1C" ,"#FF7F00","#6A3D9A",'#F00000')
names(Color1) <-unique(pbmc_harmony$cluster)

Color2 <- c("#A6CEE3","#1F78B4","#CAB2D6" ,"#33A02C","#E31A1C")
names(Color2) <- unique(pbmc_harmony$celltype)

Color3 <- c("#A6CEE3","#1F78B4","#CAB2D6" ,"#33A02C","#E31A1C" ,"#FF7F00")
names(Color3) <- unique(pbmc_harmony$celltype_2)

Color4 <- brewer.pal(3, 'Set3')
names(Color4) <- unique(pbmc_harmony$tissue)

pdf('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/UMAP_GSE_nc.pdf',width = 18,height = 10)
DimPlot(pbmc, group.by = "cluster",cols = Color1)+
  DimPlot(pbmc, group.by = "celltype",cols = Color2)+
  DimPlot(pbmc, group.by = "tissue",cols = Color4)+
  DimPlot(pbmc, group.by = "study")+
  DimPlot(pbmc, group.by = "type")
dev.off()

pdf('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/UMAP_GSE_nc_harmony.pdf',width = 18,height = 10)
DimPlot(pbmc_harmony, group.by = "cluster",cols = Color1)+
  DimPlot(pbmc_harmony, group.by = "celltype",cols = Color2)+
  DimPlot(pbmc_harmony, group.by = "celltype_2",cols = Color3)+
  DimPlot(pbmc_harmony, group.by = "tissue",cols = Color4)+
  DimPlot(pbmc_harmony, group.by = "study")+
  DimPlot(pbmc_harmony, group.by = "type")
dev.off()


pdf('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/TSNE_GSE_nc.pdf',width = 18,height = 10)
TSNEPlot(pbmc, group.by = "cluster",cols = Color1)+
  TSNEPlot(pbmc, group.by = "celltype",cols = Color2)+
  TSNEPlot(pbmc, group.by = "tissue",cols = Color4)+
  TSNEPlot(pbmc, group.by = "study")+
  TSNEPlot(pbmc, group.by = "type")
dev.off()

pdf('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/TSNE_GSE_nc_harmony.pdf',width = 18,height = 10)
TSNEPlot(pbmc_harmony, group.by = "cluster",cols = Color1)+
  TSNEPlot(pbmc_harmony, group.by = "celltype",cols = Color2)+
  TSNEPlot(pbmc_harmony, group.by = "celltype_2",cols = Color3)+
  TSNEPlot(pbmc_harmony, group.by = "tissue",cols = Color4)+
  TSNEPlot(pbmc_harmony, group.by = "study")+
  TSNEPlot(pbmc_harmony, group.by = "type")
dev.off()

## 区分B和Mono
pdf('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/B_Mono.pdf',width = 12,height = 6)
FeaturePlot(pbmc_harmony,c('CD19','CD14'),reduction = 'tsne')
dev.off()

saveRDS(pbmc,'/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/pbmc_batch.RDS')
saveRDS(pbmc_harmony,'/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/pbmc.RDS')


#### all clusters ratio in each sample
rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
pbmc <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/pbmc.RDS')

phylog_df <- pbmc@meta.data[,c('orig.ident',"celltype_2")]
phylog_df <- table(phylog_df$orig.ident,phylog_df[,"celltype_2"])
phylog_df <- data.frame(phylog_df)
colnames(phylog_df) <- c('SampleID','CellType','Freq')
phylog_df$CellType <- factor(phylog_df$CellType)

p1 <- ggplot(phylog_df,aes(x=SampleID,y=Freq,fill=CellType))+
  geom_col(position = "fill", width = 0.6)+
  coord_flip()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(face = 'bold',size = 8),
        axis.line = element_line(size=1, colour = "black"))+
  labs(x='',y='Proportion of cells')+theme(legend.position="right")+
  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#CAB2D6" ,"#33A02C","#E31A1C" ,"#FF7F00"))
ggsave('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/seurat/celltyperatio_bysample.pdf',p1)

