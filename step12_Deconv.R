#### Deconv bulk RNA-seq
rm(list=ls())
library(MuSiC)
library(Biobase)
library(Seurat)
library(xbioc)
library(BisqueRNA)
library(ggplot2)
library(SummarizedExperiment)
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Bisque/')

#### 1.load data ####
## 1.1 load sc data
pbmc <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/pbmc.RDS')
meta <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/meta.RDS')
pbmc@meta.data <- meta
subtype <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/TCAsubtype_new/subtype_info.RDS')
low_cells <- rownames(subtype)[subtype$subtype=='Low'] 
neutral_cells <- rownames(subtype)[subtype$subtype=='Neutral'] 
high_cells <- rownames(subtype)[subtype$subtype=='High']

Idents(pbmc) <- 'type'
subpbmc <- subset(pbmc,idents='ALL')
table(subpbmc$celltype_abbr)

subpbmc@meta.data$subtype <- as.character(subpbmc$celltype_abbr)
subpbmc@meta.data$subtype[which(colnames(subpbmc) %in% high_cells)] <- 'High'
subpbmc@meta.data$subtype[which(colnames(subpbmc) %in% neutral_cells)] <- 'Neutral'
subpbmc@meta.data$subtype[which(colnames(subpbmc) %in% low_cells)] <- 'Low'
print(table(subpbmc@meta.data$subtype))

## 1.2 load bulk data
df <- readRDS('/remote-home/yanzijun/CRU/ped_M/res/RNAseq/exp_ALL.count.rds')
bulk.mtx <- as.matrix(na.omit(df))

#### 2. Deconv ####
study='TARGET'
##### 2.1 Bisque #####
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
saveRDS(Bisque,paste(study,'_Bisque.rds',sep=''))


#### 3. plot estimated cell type proportions #####
study <- 'TARGET'
Bisque <- readRDS(paste(study,'_Bisque.rds',sep=''))

library(RColorBrewer)
col = colorRampPalette(brewer.pal(8, "Paired"))(8)
Color <- c('#E88482','#BEBEBE','#63A8D2',col[2:length(CT)])

CT <- c('High','Neutral','Low','Ery.','T','B','Mono.','NK')
names(Color) <- CT
print(Color)

## 样本按照subtype排列
subtype_bulk <- read.csv('/remote-home/yanzijun/CRU/ped_M/res/GSEA/subtype_ALL.csv',row.names = 1)
subtype_bulk <- as.data.frame(t(subtype_bulk['TCA cycle',]))
subtype_bulk <- subtype_bulk[order(subtype_bulk$`TCA cycle`),,drop=FALSE]
subtype_bulk <- na.omit(subtype_bulk,drop=FALSE)

plot_df <- as.data.frame(Bisque$bulk.props)

ol <- intersect(rownames(subtype_bulk),colnames(plot_df))

subtype_bulk <- subtype_bulk[ol,,drop=FALSE]
plot_df <- dplyr::select(plot_df,rownames(subtype_bulk))
print(all(rownames(subtype_bulk)==colnames(plot_df)))

plot_data <- reshape2::melt(as.matrix(plot_df))
colnames(plot_data) <- c('subtype','SampleID','Freq')
plot_data$subtype <- factor(plot_data$subtype,levels = rev(CT))


p2 <- ggplot(plot_data,aes(x=SampleID,y=Freq,fill=subtype))+
  geom_col(stat = "identity", width = 0.5,position=position_stack())+
  #scale_y_continuous(labels = scales::percent_format(scale = 1))+
  #coord_flip()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(colour = 'black',size = 8),
        axis.text.x = element_text(colour = 'black',size = 8,hjust = 0.5),
        axis.line = element_line(size=1, colour = "black"))+
  labs(x='',y='Proportion of cells')+theme(legend.position="right")+
  scale_fill_manual(values = Color)+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_blank())

ggsave(paste(study,'_bulk_ratio.pdf',sep=''),p2,width = 16,height = 4)


celltype='High'
tmp <- plot_data$Freq[plot_data$subtype==celltype]
cor.test(1:length(tmp),tmp) #0.311294 

celltype='Neutral'
tmp <- plot_data$Freq[plot_data$subtype==celltype]
cor.test(1:length(tmp),tmp) #0.2549766

celltype='Low'
tmp <- plot_data$Freq[plot_data$subtype==celltype]
cor.test(1:length(tmp),tmp) #0.6109643
