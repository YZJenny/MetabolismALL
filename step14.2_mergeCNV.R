
rm(list=ls())
library(data.table)
library(dplyr)
setwd('/mdshare/node9/yanzijun/CRU/ped_M/data/CNV/gdc_download/')

#### 0. load marker file ####
marker <- read.table('/mdshare/node9/yanzijun/CRU/ped_M/data/CNV/annoProbe/Marker.txt',
                     header = F,sep='\t')[,-1]
head(marker)
marker$V2=paste('chr',marker$V2,sep='')
marker.lst <- paste(marker$V2,marker$V3,sep=':')
result <- marker %>%
  group_by(V2) %>%
  summarise(Max_Value = max(V3),
            Min_Value = min(V3))
print(result)


#### 1. merge CNV files ####
CNV <- read.table('../gdc_sample_sheet.2023-07-27.tsv',header = T,sep='\t')
CNV[1:2,]

CNV$Sample.ID.new <- NULL
for(i in 1:length(CNV$Sample.Type)){
  element = CNV$Sample.Type[i]
  split_element <- strsplit(element, ', ')[[1]]
  index <- grep('Cancer|Tumor',split_element)
  
  CNV$Sample.ID.new[i] <- strsplit(CNV$Sample.ID[i], ', ')[[1]][index]
}
CNV[1:2,]
table(CNV$Project.ID)


cancer.lst <- c('TARGET-ALL','TARGET-AML','TARGET-CCSK','TARGET-OS')

for(i in 1:length(cancer.lst)){
  cancer.type=cancer.lst[i]
  cancer.name=unlist(strsplit(cancer.type,'-'))[2]
  print(cancer.name)
  
  sample.mtx <- CNV[grep(cancer.type,CNV$Project.ID),]
  
  arraylist <- data.frame(array=sample.mtx[,'Sample.ID.new'])
  write.table(arraylist,file = paste('/mdshare/node10/yzj/CRU/ped_M/res/GISTIC2/arraylistfile_',cancer.name,'.txt',sep=''),
              sep='\t',col.names = TRUE,row.names = FALSE,quote = FALSE)
  
  File2Case <- sample.mtx$Sample.ID.new
  names(File2Case) <- apply(as.matrix(sample.mtx$File.Name),1,function(x) unlist(strsplit(x,'\\.'))[2])
  
  file_list <- paste(sample.mtx$File.ID,sample.mtx$File.Name,sep='/')
  print(length(file_list))
  
  if(length(file_list)>0){ 
    mat_list <- lapply(file_list, function(x) {
      maf <- fread(x,fill = TRUE,skip = 1) 
    })
    #print(lapply(mat_list,dim))
    merged_cnv <- do.call(rbind, mat_list)
    
    col_1 <- File2Case[merged_cnv$V1]
    merged_cnv$V1=NULL

    write.table(merged_cnv,
                file = paste('/mdshare/node10/yzj/CRU/ped_M/res/GISTIC2/segmentationfile_',cancer.name,'.bed',sep=''),
                sep='\t',col.names = FALSE,row.names = FALSE,quote = FALSE)
    

    ### 2. calculate probe number  ####

    system(
      paste("~/Downloads/bedtools2/bin/bedtools intersect -a /mdshare/node10/yzj/CRU/ped_M/res/GISTIC2/segmentationfile_",cancer.name,
            ".bed -b /mdshare/node9/yanzijun/CRU/ped_M/data/CNV/annoProbe/GPL6801_hg38.bed -c > /mdshare/node10/yzj/CRU/ped_M/res/GISTIC2/segmentationfile_",cancer.name,'.tmp',sep=''))
    
    tmp <- read.table(paste('/mdshare/node10/yzj/CRU/ped_M/res/GISTIC2/segmentationfile_',cancer.name,'.tmp',sep=''),
                      header = FALSE,sep='\t')
    colnames(tmp) <- c('Chromosome','Start','End','Copy_Number','Major_CN','Minor_CN','Num_Probes')
    tmp$Segment_mean <- log(tmp$Copy_Number/2,2)
    tmp <- tmp[,c('Chromosome','Start','End','Num_Probes','Segment_mean')]
    
    tmp <- as.data.frame(cbind(col_1,tmp))
    
    start.lst=paste(tmp$Chromosome,tmp$Start,sep=':')
    end.lst=paste(tmp$Chromosome,tmp$End,sep=':')
    CNV.df  <- tmp[start.lst %in% marker.lst & end.lst %in% marker.lst,]

    CNV.df$Chromosome <- gsub('chr','',CNV.df$Chromosome) 
    
    write.table(CNV.df,
                paste('/mdshare/node10/yzj/CRU/ped_M/res/GISTIC2/segmentationfile_',cancer.name,'.txt',sep=''),
                row.names = F,col.names = F,quote = F,sep='\t')
  }
}

