rm(list=ls())
library(ssmarina)
library(mixtools)
library(parallel)
library(igraph)
library(ggraph)
library(tidygraph)
source('/remote-home/yanzijun/CRU/ped_M/scr/aracne2regulon_self.R')

## step1. construct regulon obj from network.txt
get_regulon <- function(cancer){
  dset <- readRDS(paste('/remote-home/yanzijun/CRU/ped_M/res/RNAseq/exp_',cancer,'.count.zscore.rds',sep=''))
  regulon <- aracne2regulon_self(paste('/remote-home/yanzijun/CRU/ped_M/res/ARACNe/',cancer,'/',ntfile,sep=''), dset, format = '3col')
  saveRDS(regulon,paste('/remote-home/yanzijun/CRU/ped_M/res/ARACNe/',cancer,'/regulon.RDS',sep=''))
}

## step2. run MARINA
get_MR <- function(pathway,cancer){
  dset <- readRDS(paste('/remote-home/yanzijun/CRU/ped_M/res/RNAseq/exp_',cancer,'.count.zscore.rds',sep=''))
  regulon <- readRDS(paste('/remote-home/yanzijun/CRU/ped_M/res/ARACNe/',cancer,'/regulon.RDS',sep=''))

  subtype <- as.data.frame(t(readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Bisque/kmeans_result_ALL_subtype.RDS')))
  
  High <- which(colnames(dset) %in%
                na.omit(rownames(subtype)[subtype[,pathway]=='High']))
  Low <- which(colnames(dset) %in%
                  na.omit(rownames(subtype)[subtype[,pathway]=='Low']))
  
  if(length(High) >2 & length(Low) >2){ ##
    signature <- ssmarina::rowTtest(dset[, High], dset[, Low])
    signature <- (qnorm(signature$p.value/2, lower.tail=F) * sign(signature$statistic))[, 1]
    
    nullmodel <- ssmarina::ttestNull(dset[, High], dset[, Low], per=1000, repos=T)
    mrs <- ssmarina::marina(signature, regulon, nullmodel)

    mrshadow <- ssmarina::shadow(mrs, pval=0.05,targets = 10)
    mrshadow <- ssmarina::ledge(mrshadow)
    saveRDS(mrshadow,
            paste('ssMARINA/All/',cancer,'_',pathway,'.RDS',sep=''))
    return(mrshadow)
  }else{
    warnings='No records'
    print(warnings)
  }
}

## step3. pick significant MR
get_sigMR <- function(pathway,cancer,FDR=0.1){
  file=paste('ssMARINA/All/',cancer,'_',pathway,'.RDS',sep='')
  if(file.exists(file)){
    mrshadow <- readRDS(file)
    mrshadow.summary <- summary(mrshadow)
    MARINA.results <- mrshadow.summary$MARINA.results
    if(length(which(MARINA.results$FDR < FDR))>0){
      sigMR <- data.frame(pathway=pathway,TF=MARINA.results$Regulon[MARINA.results$FDR<fdr_MARINA],
                          cancer=cancer)
      saveRDS(sigMR,
              paste('ssMARINA/FDR',FDR,'/',cancer,'_',pathway,'.RDS',sep=''))
    }else{
      sigMR <- c()
    }
    return(sigMR)
  }else{
    warnings='no File'
    print(warnings)
  }
}


############################
setwd('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25')
print(getwd())

#setwd('/remote-home/yanzijun/CRU/ped_M/res/ARACNe')
ntfile <- 'network.txt'

pathway.lst <- c("Amino acid","Carbohydrate","Energy","Lipid",
                "Nucleotide","TCA cycle","Vitamin cofactor")
cancer.lst <- c('ALL','AML','CCSK','NBL','OS','RT','WT')

# 生成两个向量的两两组合
combinations <- expand.grid(pathway.lst, cancer.lst)


#### 1. construct regulon obj from network.txt #####
print('step1')
tmp <- lapply(cancer.lst, get_regulon)


#### 2. run MARINA ####
print('step2')
# cl <- makeCluster(detectCores())
library("future")
cl <- makeClusterPSOCK(124, revtunnel = TRUE, outfile = "", verbose = TRUE)
clusterExport(cl, "get_MR")
chunk_size <- ceiling(nrow(combinations) / length(cl))
chunks <- split(combinations, rep(1:length(cl), each = chunk_size, length.out = nrow(combinations)))

MR <- clusterApply(cl, chunks, function(chunk) {
  apply(chunk, 1, function(row) get_MR(row[[1]], row[[2]]))
})
stopCluster(cl)

#### 3.pick significant MR ####
print('step3')
fdr_MARINA=0.1
sigMR <- apply(combinations, 1, function(row) get_sigMR(row[[1]], row[[2]], FDR=fdr_MARINA))