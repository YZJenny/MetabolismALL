#### Fig4A. plot estimated cell type proportions #### 
get_ratioplot <- function(study,subtype_bulk,method='Bisque'){
  CT <- c('High','Middle','Low','Ery.','T','B','Mono.','NK')
  Color <- c('#EAB080','#d6d6d6','#9cb0c3',
             "#1F78B4" ,"#B2DF8A","#33A02C", "#FB9A99","#FF7F00")
  names(Color) <- CT
  print(Color)
  
  if(method=='Bisque'){
    res <- readRDS(paste('Bisque/',study,'_Bisque.rds',sep=''))
    plot_df <- as.data.frame(res$bulk.props)
  }else if(method=='CIBERSORTx'){
    res <- read.table(paste('CIBERSORTx/',study,'_CIBERSORTx.txt',sep=''),sep='\t',header = T)
    rownames(res) <- res$Mixture
    res$Mixture=res$P.value=res$Correlation=res$RMSE=NULL
    plot_df <- as.data.frame(t(res))
    rownames(plot_df)[7:8] <- c('Mono.','Ery.')
  }
  
  subtype_bulk <- data.frame(subtype = subtype_bulk, row.names = names(subtype_bulk))
  ol <- intersect(rownames(subtype_bulk),colnames(plot_df))
  
  subtype_bulk <- subtype_bulk[ol,,drop=FALSE]
  plot_df <- dplyr::select(plot_df,rownames(subtype_bulk))
  print(all(rownames(subtype_bulk)==colnames(plot_df)))
  
  ### 总图
  plot_data <- reshape2::melt(as.matrix(plot_df))
  colnames(plot_data) <- c('subtype','SampleID','Freq')
  plot_data$subtype <- factor(plot_data$subtype,levels = rev(CT)) 
  
  p1 <- ggplot(plot_data,aes(x=SampleID,y=Freq,fill=subtype))+
    geom_col(stat = "identity", width = 0.5,position=position_stack())+
    scale_y_continuous(expand = c(0,0))+
    #coord_flip()+
    theme_bw()+
    theme(panel.grid.major=element_line(colour=NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.minor = element_blank())+
    theme(axis.text = element_text(colour = 'black',size = 8))+
    labs(x='',y='% of cells')+theme(legend.position="none")+
    scale_fill_manual(values = Color)+
    theme(axis.text.x = element_blank())+
    theme(axis.ticks.x = element_blank())+
    geom_vline(xintercept = as.numeric(table(subtype_bulk))[1])+
    geom_vline(xintercept = as.numeric(table(subtype_bulk))[1]+as.numeric(table(subtype_bulk))[2])
  ggsave(paste(method,'/',study,'_bulk_ratio.pdf',sep=''),p1,width = 16,height = 4)
  
  
  ### 分图：正常细胞
  normal <- plot_df[!rownames(plot_df) %in% c('High','Middle','Low'),]
  plot_normal <- reshape2::melt(as.matrix(normal))
  colnames(plot_normal) <- c('subtype','SampleID','Freq')
  plot_normal$subtype <- factor(plot_normal$subtype,levels = rev(CT)) 
  
  p_normal <- ggplot(plot_normal,aes(x=SampleID,y=Freq,fill=subtype))+
    geom_col(stat = "identity", width = 0.5,position=position_stack())+
    scale_y_continuous(expand = c(0,0))+
    #coord_flip()+
    theme_bw()+
    theme(panel.grid.major=element_line(colour=NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.minor = element_blank())+
    theme(axis.text = element_text(colour = 'black',size = 8))+
    labs(x='',y='% of cells')+theme(legend.position="none")+
    scale_fill_manual(values = Color)+
    theme(axis.text.x = element_blank())+
    theme(axis.ticks.x = element_blank())+
    geom_vline(xintercept = as.numeric(table(subtype_bulk))[1])+
    geom_vline(xintercept = as.numeric(table(subtype_bulk))[1]+as.numeric(table(subtype_bulk))[2])
  ggsave(paste(method,'/',study,'_bulk_ratio_normal.pdf',sep=''),p_normal,width = 16,height = 4)
  
  ### 分图：肿瘤细胞
  tumor <- plot_df[c('High','Middle','Low'),]
  
  # 找出总和为0的样本并删除
  row_sums <- colSums(tumor)
  zero_sum_samples <- names(row_sums[row_sums == 0])
  print(length(zero_sum_samples))
  tumor <- tumor[,!colnames(tumor) %in% zero_sum_samples]
  
  plot_tumor <- reshape2::melt(as.matrix(tumor))
  colnames(plot_tumor) <- c('subtype','SampleID','Freq')
  plot_tumor$subtype <- factor(plot_tumor$subtype,levels = rev(CT))
  
  p_tumor <- ggplot(plot_tumor,aes(x=SampleID,y=Freq,fill=subtype))+
    geom_col(stat = "identity", width = 0.5,position='fill')+
    scale_y_continuous(expand = c(0,0))+
    #coord_flip()+
    theme_bw()+
    theme(panel.grid.major=element_line(colour=NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.minor = element_blank())+
    theme(axis.text = element_text(colour = 'black',size = 8))+
    labs(x='',y='% of tumor cells')+theme(legend.position="none")+
    scale_fill_manual(values = Color)+
    theme(axis.text.x = element_blank())+
    theme(axis.ticks.x = element_blank())+
    geom_vline(xintercept = as.numeric(table(subtype_bulk))[1])+
    geom_vline(xintercept = as.numeric(table(subtype_bulk))[1]+as.numeric(table(subtype_bulk))[2])
  ggsave(paste(method,'/',study,'_bulk_ratio_tumor.pdf',sep=''),p_tumor,width = 6,height = 4)
  return(plot_data)
}

kmeans_result <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Bisque/kmeans_result_ALL.RDS')
subtype_bulk <- kmeans_result$cluster
plot_ALL <- get_ratioplot(study = 'TARGET_ALL',subtype_bulk = subtype_bulk,method = 'Bisque')


#### Fig4B. survplot of TCA cycle subtype #### 
rm(list=ls())
{
  library(readxl)
  library(dplyr)
  library(survival)
  library(survminer)
  library(survivalROC)
  setwd('/remote-home/yanzijun/CRU/ped_M/')
}

## load clinical data
{
  raw_clin <- read.csv('data/RNAseq_clinical/TARGET_ALL_clinical.tsv',sep='\t',row.names = 1)
  tmp1 <- read_excel("data/RNAseq_clinical/TARGET_ALL_ClinicalData_Phase_II_Validation_20211118.xlsx")
  tmp2 <- read_excel("data/RNAseq_clinical/TARGET_ALL_ClinicalData_Phase_II_Discovery_20211118.xlsx")
  print(all(colnames(tmp1)==colnames(tmp2)))
  
  tmp <- rbind(tmp1,tmp2)
  rm(tmp1);rm(tmp2) 
}

clin <- merge(raw_clin[,c('case_submitter_id','primary_diagnosis')],tmp,by.x = 'case_submitter_id',by.y = 'TARGET USI')
head(clin) 


subtype <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Bisque/kmeans_result_ALL_subtype.RDS')
subtype$class=NULL
subtype$sampleID <- apply(as.matrix(rownames(subtype)),1,
                          function(x) paste(unlist(strsplit(x,'\\.'))[1:3],collapse ='-'))
subtype <- tibble::rownames_to_column(subtype,'ID')


merge_df <- merge(clin,subtype,by.x='case_submitter_id',by.y = 'sampleID')
merge_df$subtype <- factor(merge_df$subtype,levels = c('High','Middle','Low'))

print(table(merge_df$primary_diagnosis))
print(table(merge_df$`Cell of Origin`))

get_surv <- function(df,type){
  cox.res <- coxph(Surv(time, status) ~ subtype, data = df)
  HR <- signif(exp(cox.res$coefficients[2]), digits = 2)
  
  fit <- survfit(Surv(time = time, event = status) ~ subtype, data = df)
  pvalue <- surv_pvalue(fit)$pval
  print(pvalue)
  
  plt <- ggsurvplot(fit,size=0.3,
                    censor.size=1,
                    pval = TRUE, pval.size = 2,
                    ggtheme = theme_survminer(
                      base_size = 5,
                      font.main = c(6, "plain", "black"),
                      font.submain = c(6, "plain", "black"),
                      font.caption = c(6, "plain", "black"),
                      font.x = c(5, "plain", "black"),
                      font.y = c(6, "plain", "black"),
                      font.tickslab = c(5, "plain", "black")) %+replace%
                      theme(plot.title=element_text(hjust=0.5)), 
                    title = type,
                    palette = c("#63A8D2", "grey", "#E88482"),
                    xlab = "Months", legend = "none",
                    linetype =1)
  return(plt)
}

EFS_df <- merge_df[,c('subtype','First Event','Event Free Survival Time in Days')]
colnames(EFS_df)[2:3] <- c('status','time')
print(table(EFS_df$status))
EFS_df <- na.omit(EFS_df[EFS_df$status != 'Censored',])
EFS_df$status <- ifelse(EFS_df$status == "None", 0, 1)
EFS_df$time <- round(EFS_df$time/12, 2)

OS_df <- merge_df[,c('subtype','Vital Status','Overall Survival Time in Days')]
colnames(OS_df)[2:3] <- c('status','time')
print(table(OS_df$status))
OS_df <- na.omit(OS_df[OS_df$status %in% c('Alive','Dead'),])
OS_df$status <- ifelse(OS_df$status == "Alive", 0, 1)
OS_df$time <- round(OS_df$time/12, 2)

p_EFS <- get_surv(EFS_df,'EFS')
p_OS <- get_surv(OS_df,'OS')

p <- arrange_ggsurvplots(list(p_EFS,p_OS), ncol = 2, print = FALSE)  
ggsave(filename = paste("res_FDR0.25/Surv/surv_ALL_TCAcycle.pdf", sep = ""), p,
       width = 6, height = 3, units = "cm")


#### Fig4C. dotplot of PROGENy  #### 
## https://cloud.tencent.com/developer/article/2206142
rm(list=ls())
library(progeny)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(forcats)
library(ggstance)
library(pheatmap)

setwd("/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/")
exp <- readRDS('../res/RNAseq/exp_ALL.count.zscore.rds')
exp <- na.omit(exp)
group <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Bisque/kmeans_result_ALL_subtype.RDS')
group$class=NULL

group <- group[order(group$subtype),,drop=FALSE]
sub.group <- group[group$subtype != 'Middle',,drop=FALSE]

# 1. progeny
res <- progeny(exp, scale=TRUE,
               organism="Human",
               top = 100, perm = 1)
saveRDS(res,'PROGENy/ALL_TCA cycle.RDS')
#res <- readRDS('PROGENy/ALL_TCA cycle.RDS')

sub.res <- res[match(rownames(sub.group),rownames(res)),]
print(all(rownames(sub.res)==rownames(sub.group)))


# 2. claulate difference between High and Low subtype
controls =ifelse(sub.group$subtype=='High',FALSE,TRUE) 

result = apply(sub.res, 2, function(x) {
  broom::tidy(lm(x ~ !controls)) %>%
    filter(term == "!controlsTRUE") %>%
    dplyr::select(-term)
})
res_lm <- mutate(bind_rows(result), signaling=names(result))
res_lm$p.adjust <- p.adjust(res_lm$p.value,method = 'fdr')
res_lm$signif <- ' '
res_lm$signif[res_lm$p.value <0.05 & res_lm$p.value > 0.01] <- '*'
res_lm$signif[res_lm$p.value < 0.01 & res_lm$p.value > 0.001] <- '**'
res_lm$signif[res_lm$p.value < 0.001 & res_lm$p.value >0.0001] <- '***'
res_lm$signif[res_lm$p.value < 0.0001] <- '****'
res_lm$y_location <- ifelse(res_lm$statistic>0,res_lm$statistic+0.5,res_lm$statistic-0.5) 
print(res_lm)
saveRDS(res_lm,'PROGENy/lm_ALL_TCA cycle.RDS')


#res_lm <- readRDS('PROGENy/lm_ALL_TCA cycle.RDS')
# 3.plot
p1 <- ggplot(as.data.frame(res_lm),aes(x=signaling,y=statistic))+
  geom_point(aes(size=-log10(p.value),color=signaling))+
  geom_segment(aes(x=signaling,xend=signaling,
                   y=0,yend=statistic,
                   color=signaling),
               cex=1)+
  theme_classic()+
  scale_color_manual(values = paletteer::paletteer_d('ggsci::category20_d3',n=14))+
  coord_flip()+
  theme(legend.position = 'left', 
        plot.title = element_text(size=12,hjust=0.5), 
        text = element_text(colour = 'black'),
        axis.text = element_text(color = 'black'))+
  geom_hline(yintercept = 0)+
  labs(title = cancer)+
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+ 
  xlab("")+ 
  ylab('')+
  annotate('text',x=1:14,y=res_lm$y_location,label=res_lm$signif) 
# theme(legend.text = element_text(size = 6))  
ggsave('PROGENy/ALL_TCA cycle.pdf',p1)


#### Fig4D. barplot of SMG/CNV ####
rm(list=ls())
library(ggplot2)

## SMG info
pvalue.SMG <- read.table("/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/MutSigCV_2018nature/SMGpvalue_ALL.txt", header = T, sep = "\t")
sigSMG <- rownames(pvalue.SMG)[which(pvalue.SMG$TCA.cycle <0.05)] 
sigSMG.p <- pvalue.SMG[which(pvalue.SMG$TCA.cycle <0.05),'TCA.cycle'] 
names(sigSMG.p) <- sigSMG
sigSMG.p <- sort(sigSMG.p)
print(sigSMG.p)

## CNV info
pvalue.CNV <- read.table("/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/GISTIC2/CNVpvalue_ALL.txt", header = T, sep = "\t")
sigCNV <- rownames(pvalue.CNV)[which(pvalue.CNV$TCA.cycle <0.05)] #
sigCNV.p <- pvalue.CNV[which(pvalue.CNV$TCA.cycle <0.05),'TCA.cycle'] 
names(sigCNV.p) <- sigCNV
sigCNV.p <- sort(sigCNV.p)
print(sigCNV.p)


plot_df <- data.frame(gene=c(names(sigCNV.p),names(sigSMG.p)),
                      pvalue=c(sigCNV.p,sigSMG.p),
                      type=c(rep('CNV',length(sigCNV.p)),rep('SMG',length(sigSMG.p))))
plot_df$gene <- factor(plot_df$gene,levels = plot_df$gene)

p <- ggplot(plot_df, aes(x = gene, y = pvalue, fill = type)) +
  geom_bar(stat='identity')+
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = c('#E88482','#63A8D2'))+
  theme( axis.line = element_line(size=0.3, colour = "black"),
         axis.ticks = element_line(colour = "black", size = 0.3),
         panel.border = element_blank(), panel.background = element_blank(),
         axis.text.x=element_text(colour="black", size = 8,angle=90,hjust=1,vjust=0.5),
         axis.text.y=element_text(colour="black", size = 8)) +
  theme(plot.margin = unit(rep(1,4),"lines"))+
  ylim(c(0,0.05))
ggsave('/remote-home/yanzijun/CRU/ped_M/res_ALL/bulkALL/Fig4D_SMG_CNV.pdf',p,
       width = 5.5,height=2.3,units="in",device="pdf",useDingbats=FALSE)
