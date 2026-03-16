#### sFig4A. PCA of subtype based on bulk RNAseq #### 
rm(list=ls())
library(Seurat)
library(clusterProfiler)
library(sva)
library(FactoMineR)
library(ggplot2)
library(factoextra)
pca.plot = function(dat,col,pal){
  df.pca <- PCA(t(dat), graph = FALSE)
  fviz_pca_ind(df.pca,
               geom.ind = "point",
               col.ind = col ,
               addEllipses = FALSE,
               legend.title = "Groups",
               palette = pal
  )
}

gmt <- read.gmt('/remote-home/yanzijun/CRU/TALL_M/data/GeneSet7.gmt')
geneset <- gmt$gene[gmt$term=='TCA cycle']


## load bulk
bulk_EXP <- readRDS('/remote-home/yanzijun/CRU/ped_M/res/RNAseq/exp_ALL.count.rds')
bulk_EXP <- bulk_EXP[rownames(bulk_EXP) %in% geneset,]

subtype <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Bisque/kmeans_result_ALL_subtype.RDS')
bulk_Subtype <- tibble::rownames_to_column(bulk_Subtype,'sampleID')
bulk_Subtype <- na.omit(bulk_Subtype)
bulk_Subtype$subtype_new <- paste('bulk_',bulk_Subtype$subtype,sep='')
bulk_Subtype$batch <- 'bulk'

bulk_EXP <- dplyr::select(bulk_EXP,bulk_Subtype$sampleID)
bulk_EXP <- na.omit(bulk_EXP)
print(all(colnames(bulk_EXP)==bulk_Subtype$sampleID))


## load sc
sc_EXP <- pbmc_tumor@assays$SCT@data
sc_EXP <- sc_EXP[rownames(sc_EXP) %in% geneset,]


sc_Subtype <- as.data.frame(pbmc_tumor$subtype)
colnames(sc_Subtype) <- 'subtype'
sc_Subtype$subtype_new <- paste('sc_',sc_Subtype$subtype,sep = '')
sc_Subtype <- tibble::rownames_to_column(sc_Subtype,'sampleID')
sc_Subtype$batch <- 'sc'
sc_Subtype$subtype[sc_Subtype$subtype=='High'] <- 'Up'
sc_Subtype$subtype[sc_Subtype$subtype=='Low'] <- 'Down'

print(all(colnames(sc_EXP)==sc_Subtype$sampleID))

## save data
saveRDS(sc_EXP, '../TCAsubtype_new/exp_tumor_TCAcycle.RDS')
saveRDS(sc_Subtype, '../TCAsubtype_new/meta_tumor_TCAcycle.RDS')


p1 <- pca.plot(sc_EXP,factor(sc_Subtype$subtype),pal=c('#63A8D2','#E88482'))  
p2 <- pca.plot(log(na.omit(bulk_EXP)+1,2),factor(bulk_Subtype$subtype),pal=c('#63A8D2','grey','#E88482')) 

pdf('/remote-home/yanzijun/CRU/ped_M/res_ALL/bulkALL/sFig4A_PCA_sc_bulk.pdf',width = 5,height = 4)
print(p1)
print(p2)
dev.off()

#### sFig4B. forest tree of cox #### 
rm(list=ls())
library(ggplot2)
setwd('/remote-home/yanzijun/CRU/ped_M')

##
subtype <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Bisque/kmeans_result_ALL_subtype.RDS')
subtype$sampleID <- apply(as.matrix(rownames(subtype)),1,
                          function(x) paste(unlist(strsplit(x,'\\.'))[1:3],collapse ='-'))
subtype <- tibble::rownames_to_column(subtype,'ID')
subtype$subtype <- ifelse(subtype$subtype == "Low", 0,  ifelse(subtype$subtype == "Middle", 1, 2))

## load clinical data
raw_clin <- read.csv('/remote-home/yanzijun/CRU/ped_M/data/RNAseq_clinical/TARGET_ALL_clinical.tsv',sep='\t',header = TRUE)
clin <- raw_clin[,c('case_submitter_id','primary_diagnosis','ethnicity','gender','vital_status','days_to_last_follow_up')]
colnames(clin)[1] <- 'sampleID'

clin <- clin[clin$vital_status %in% c("Dead", "Alive"), ]
clin$status <- ifelse(clin$vital_status == "Dead", 1, 0)
clin$days_to_last_follow_up <- as.numeric(clin$days_to_last_follow_up)
clin$time <- round(clin$days_to_last_follow_up /12, 2)

clin$group <- NA
clin$group[clin$primary_diagnosis=='T lymphoblastic leukemia/lymphoma'] <- 'TALL'
clin$group[clin$primary_diagnosis %in% c('Precursor B-cell lymphoblastic leukemia',
                                         'B lymphoblastic leukemia/lymphoma, NOS')] <- 'BALL'
clin$group[grep('^Mixed phenotype acute leukemia',clin$primary_diagnosis)] <- 'MPAL'
clin <- na.omit(clin);
clin <- clin[clin$group %in% c('BALL','TALL'),]
clin$group <- ifelse(clin$group == "TALL", 1, 0)

clin$primary_diagnosis=NULL

clin$gender[clin$gender=='unknown'] <- 'Unknown'
clin$gender <- ifelse(clin$gender == "female", 1, 0)

clin$ethnicity <-  ifelse(clin$ethnicity == "not hispanic or latino", 0,  ifelse(clin$ethnicity == "hispanic or latino", 1, 2))

## merge data
anno <- na.omit(merge(subtype,clin,all.x = TRUE))
rownames(anno) <- anno$ID
anno <- anno[order(anno$subtype),]
head(anno)
anno$sampleID=anno$ID=anno$vital_status=anno$days_to_last_follow_up=NULL
anno <- anno[,c('time','status','subtype','group','gender','ethnicity')]

library(survival)
pfilter <- 0.05   
uniresult <- data.frame()  
df <- anno
for(i in colnames(df[,3:ncol(df)])){   
  unicox <- coxph(Surv(time = time, event = status) ~ df[,i], data = df)  
  unisum<- summary(unicox)   
  pvalue <- round(unisum$coefficients[,5],3) 

    uniresult <- rbind(uniresult,
                       cbind(gene=i,
                             HR=unisum$coefficients[,2],
                             L95CI=unisum$conf.int[,3],
                             H95CI=unisum$conf.int[,4],
                             pvalue=unisum$coefficients[,5]
                       ))
  
}   
## save
write.csv(uniresult,file = "res_ALL/bulkALL/uniCox.csv",row.names = F)

## plot 
rt <- uniresult
gene=rt$gene
hr=sprintf("%.3f",as.numeric(rt$HR))
hrLow=sprintf("%.3f",as.numeric(rt$L95CI))
hrHigh=sprintf("%.3f",as.numeric(rt$H95CI))
Hazard.ratio=paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal=ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", as.numeric(rt$pvalue)))

pdf(file='res_ALL/bulkALL/forest.pdf', width = 5.5, height =3.5)
n=nrow(rt)
nRow=n+1
ylim=c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))


xlim = c(0,3)
par(mar=c(4,2,1.5,1.5))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)


par(mar=c(4,1,1.5,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.03,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'blue')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()