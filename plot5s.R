rm(list=ls())
library(ssmarina)
library(mixtools)
library(parallel)
library(tidygraph)
library(ggrepel)

#### sFig5A. volcano of sigMR ####
mrs.es <- read.csv('/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/allMR.csv',header = T,row.names = 1)
sig <- read.csv('/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/sigMR.csv',header = T,row.names = 1)

tag_lst=rownames(sig)[1:10]
tag_gene=as.character(rownames(mrs.es))
index_tag=which(!tag_gene %in% tag_lst)
tag_gene[index_tag]=""

figure <- ggplot(mrs.es,aes(x=nes,y=log10FDR))+geom_point(aes(color=group))+
  scale_color_manual(values = c('#E88482','grey','#63A8D2'))+
  labs(x="Normalized enrichment score (NES)",y = "-Log10(FDR)")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=10,color='black'),panel.grid.major=element_line(colour=NA))+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

p <- figure+geom_text_repel(label=tag_gene,max.overlaps = 50)+
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)

ggsave(plot=p,
       filename = '/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/volcano.pdf',width = 5.5,height = 4)


#### sFig5B-C. boxplot of sigMR and target genes ####
rm(list=ls())
library(ggplot2)

ALL <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Bisque/kmeans_result_ALL_subtype.RDS')
ALL_High <- rownames(ALL)[ALL$subtype=='High']
ALL_Low <- rownames(ALL)[ALL$subtype=='Low']
head(ALL_High)

## process exp file
EXP <- readRDS('/remote-home/yanzijun/CRU/ped_M/res/RNAseq/exp_TARGAT.fpkm.rds')
EXP <- na.omit( as.data.frame(t(EXP)))
EXP[1:3,1:3]
EXP <- log(EXP+1,2)

## load sigMR
sig <- read.csv('/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/sigMR.csv',header = T,row.names = 1)
TF_high <- rownames(sig)[sig$es>0]
TF_low <- rownames(sig)[sig$es<0]
TF <- rownames(sig)

## load target
TF_target <- read.csv('/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/network_ALL_TCA cycle.csv')
target_high <- TF_target$target[TF_target$ledge=='Yes'&TF_target$TF %in% TF_high]
target_low <- TF_target$target[TF_target$ledge=='Yes'&TF_target$TF %in% TF_low]


get_boxplot <- function(EXP1,EXP2){
  df <- rbind(data.frame(group='High',TF=apply(EXP1, 1, mean)),
              data.frame(group='Low',TF=apply(EXP2, 1, mean)))
  mycol <- c('#E88482','#63A8D2')
  names(mycol) <- c('High','Low')
  
  res <- t.test(apply(EXP1, 1, mean),apply(EXP2, 1, mean))
  print(res$p.value)
  
  
  p <- ggplot(data = df) + 
    geom_boxplot(aes(x = group, y = TF,
                     fill = factor(group))) + 
    scale_fill_manual(values = mycol)+
    scale_color_manual(values = mycol)+
    theme_classic()+
    labs(title='',y='Expression',x='')+
    theme(legend.position ='none', axis.text = element_text(color = 'black'))
  return(p)
}

## HighTF
EXP_high_h <- EXP[rownames(EXP) %in% ALL_High,colnames(EXP) %in% TF_high]
EXP_low_h <- EXP[rownames(EXP) %in% ALL_Low,colnames(EXP) %in% TF_high]
p_H <- get_boxplot(EXP_high_h,EXP_low_h)
ggsave(paste('/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/boxplot_highTF.pdf',sep=''),p_H,
       width = 2.5,height = 3.5)

## LowTF
EXP_high_l <- EXP[rownames(EXP) %in% ALL_High,colnames(EXP) %in% TF_low]
EXP_low_l <- EXP[rownames(EXP) %in% ALL_Low,colnames(EXP) %in% TF_low]
p_L <- get_boxplot(EXP_high_l,EXP_low_l)
ggsave(paste('/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/boxplot_lowTF.pdf',sep=''),p_L,
       width = 2.5,height = 3.5)


## HighTarget
EXP_high_h <- EXP[rownames(EXP) %in% ALL_High,colnames(EXP) %in% target_high]
EXP_low_h <- EXP[rownames(EXP) %in% ALL_Low,colnames(EXP) %in% target_high]
p_H <- get_boxplot(EXP_high_h,EXP_low_h)
ggsave(paste('/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/boxplot_highTarget.pdf',sep=''),p_H,
       width = 2.5,height = 3.5)

## LowTarget
EXP_high_l <- EXP[rownames(EXP) %in% ALL_High,colnames(EXP) %in% target_low]
EXP_low_l <- EXP[rownames(EXP) %in% ALL_Low,colnames(EXP) %in% target_low]
p_L <- get_boxplot(EXP_high_l,EXP_low_l)
ggsave(paste('/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/boxplot_lowTarget.pdf',sep=''),p_L,
       width = 2.5,height = 3.5)



#### sFig5D-E. survplot sigMR/target ####
rm(list=ls())
library(ggplot2)
library(survival)
library(survminer)
library(survivalROC)


clin <- read.csv('/remote-home/yanzijun/CRU/ped_M/data/RNAseq_clinical/TARGET_clinical.csv')
table(clin$project_id)
sub.clin <- clin[clin$vital_status %in% c('Dead','Alive'),
                 c('case_submitter_id','project_id','vital_status','days_to_last_follow_up')]
sub.clin$status <- ifelse(sub.clin$vital_status=='Dead',1,0)
sub.clin$days_to_last_follow_up <- as.numeric(sub.clin$days_to_last_follow_up)
sub.clin$time <- round(sub.clin$days_to_last_follow_up/12,2)
sub.clin$vital_status=sub.clin$days_to_last_follow_up=NULL
sub.clin <- na.omit(sub.clin)

## load subtype
ALL <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Bisque/kmeans_result_ALL_subtype.RDS')
ALL_High <- rownames(ALL)[ALL$subtype=='High']
ALL_Low <- rownames(ALL)[ALL$subtype=='Low']
head(ALL_High)

## load exp
EXP <- readRDS('/remote-home/yanzijun/CRU/ped_M/res/RNAseq/exp_TARGAT.fpkm.rds')
EXP <- as.data.frame(t(EXP))
EXP[1:3,1:3]

sig <- read.csv('/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/sigMR.csv',header = T,row.names = 1)
TF_high <- rownames(sig)[sig$es>0]
TF_low <- rownames(sig)[sig$es<0]
TF <-  rownames(sig)

## load targets
TF_target <- read.csv('/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/network_ALL_TCA cycle.csv')
target_high <- TF_target$target[TF_target$ledge=='Yes'&TF_target$TF %in% TF_high]
target_low <- TF_target$target[TF_target$ledge=='Yes'&TF_target$TF %in% TF_low]

## exp of sigMR
EXP_high <- EXP[,colnames(EXP) %in% TF_high]
EXP_low <- EXP[,colnames(EXP) %in% TF_low]
EXP_all <- EXP[,colnames(EXP) %in% TF]

## exp of targets
EXP_high_target <- EXP[,colnames(EXP) %in% target_high]
EXP_low_target <- EXP[,colnames(EXP) %in% target_low]

## exp of GRN
EXP_high_module <- EXP[,colnames(EXP) %in% c(TF_high,target_high)]
EXP_low_module <- EXP[,colnames(EXP) %in% c(TF_low,target_low)]

## survival anaysis
subtype <- 'Hightarget'
if(subtype=='LowTF'){
  EXP_TF <- EXP_low  
}else if(subtype=='HighTF'){
  EXP_TF <- EXP_high 
}else if(subtype=='AllTF'){
  EXP_TF <- EXP_all 
}else if(subtype=='Lowtarget'){
  EXP_TF <- EXP_low_target
}else if(subtype=='Hightarget'){
  EXP_TF <- EXP_high_target
}else if(subtype=='Lowmodule'){
  EXP_TF <- EXP_low_module
}else if(subtype=='Highmodule'){
  EXP_TF <- EXP_high_module
}

EXP_tmp <- data.frame(case_submitter_id=rownames(EXP_TF),TF=apply(EXP_TF, 1, median))
head(EXP_tmp)

new_ID <-   apply(as.matrix(rownames(EXP_tmp)),1,
                  function(x) paste(unlist(strsplit(x,'\\.'))[1:3],collapse ='-')) 
EXP_tmp$case_submitter_id <- new_ID
EXP_tmp <- EXP_tmp[!duplicated(EXP_tmp$case_submitter_id),] 
head(EXP_tmp)

## merge clin info and exp
data <- merge(sub.clin,EXP_tmp,by='case_submitter_id')
print(table(data$project_id))
dim(data)

candiGenes <-'TF'

## calculate cutoff
res.cut <- surv_cutpoint(data,time = "time", event = "status",variables = candiGenes)
res.cat <- surv_categorize(res.cut)

clin_info=res.cat[,1:2]
surtype <- res.cat[,3:ncol(res.cat)]

clin_info$surtype <- surtype
fit <- survfit(Surv(time =time, event =status) ~ surtype,data=clin_info)
#pvalue <- surv_pvalue(fit)$pval

p <- ggsurvplot(fit,data = clin_info, 
                #risk.table = TRUE,
                pval = T,palette=c("#E88482","#63A8D2"),
                title = subtype,
                ggtheme = theme_survminer(size = 1))

pdf(paste('/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/survplot_TARGET_',subtype,'.pdf',sep=''),
    width = 3,height = 3.5)
print(p)
dev.off()


#### sFig5F-G. boxplot of TF in bulkRNA-seq subtypes####
rm(list=ls())
library(ggplot2)
library(ggsignif)
library(ggpubr)
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL')

sig <- read.csv('/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/sigMR.csv',header = T,row.names = 1)
Act <- rownames(sig)[sig$group=='Activation']
Rep <-  rownames(sig)[sig$group=='Repression']

type='Act'
if(type=='Act'){
  candiTF <- Act
}else if(type=='Rep'){
  candiTF <- Rep
}

## bulk expression
TALL <- read.csv('/remote-home/yanzijun/CRU/TALL/data/2303_RNAseq_CQ/res/onlyTALL_Batch_after_gene_logTPM.csv',row.names = 1)
load('/remote-home/yanzijun/CRU/TALL/data/2303_RNAseq_CQ/res/noTALL_Batch_after_gene_logTPM.RData')

exp=cbind(TALL[candiTF,],Tcell[candiTF,])
df <- as.data.frame(t(exp))
df$tissue <- c(rep('TALL',ncol(TALL)),rep('Tcell',ncol(Tcell)))
df$tissue <- factor(df$tissue,levels = c('TALL','Tcell'))

plot_data <- reshape2::melt(df)
colnames(plot_data) <- c('group','TF','expression')
plot_data <- na.omit(plot_data)
plot_data$group <- factor(plot_data$group,levels = c('TALL','Tcell'))

print(aggregate(plot_data$expression, list(plot_data$group,plot_data$TF),summary))

p <- ggplot(plot_data, aes(x = TF, y = expression))+ 
  geom_boxplot(aes(fill = group),position=position_dodge(0.8),width=0.6)+
  theme_classic()+
  theme(axis.text = element_text(size=10,face="plain",color="black"),
        axis.title  = element_text(size=10,face="plain",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position=c(0.75, 0.9),
        legend.direction = "horizontal")+
  scale_fill_manual(values = c('#E88482','#63A8D2'))+
  stat_compare_means(aes(group = group), method ='t.test',label = "p.signif")
ggsave(paste('/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/boxplot_',type,'_bulk.pdf',sep=''),p,width = 10,height = 4)

## print p-value
for(tf in candiTF){
  high <- plot_data$expression[plot_data$TF==tf & plot_data$group=='TALL']
  low <- plot_data$expression[plot_data$TF==tf & plot_data$group=='Tcell']
  res <- t.test(high,low,alternative='greater')
  print(res$p.value)
}


#### sFig5H. ALL celline in DepMap and boxplot of DepMap of TF####
rm(list=ls())
library(RColorBrewer)
library(tibble)
mycol=c('#E88482','#63A8D2')

sig <- read.csv('/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/sigMR.csv',header = T,row.names = 1)
sig <- rownames_to_column(sig,'TF')
candiTF <- sig$TF
print(candiTF)

## Depmap score
## from scr/step11_candiTF_v2.R
eff.ALL <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/DepMap/CRISPR_gene_effect_ALL.RDS')
META <- read.csv('/remote-home/yanzijun/CRU/TALL_FM/data/DepMap/22Q1/sample_info.csv')

subMETA <- META[META$DepMap_ID %in% rownames(eff.ALL),]
write.table(subMETA,'/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/CRISPR_ALL_celline.csv',row.names = FALSE,sep=',')

df <- eff.ALL[,colnames(eff.ALL) %in% candiTF]
med <- apply(df,2,median)
rank <- names(sort(med,decreasing = F))

plot.df <- reshape2::melt(df)
colnames(plot.df) <- c('TF','EffectScore')
plot.df$TF=factor(plot.df$TF,levels = rank)
plot.df <- merge(plot.df,sig[,c('TF','group')],by='TF')

p.flip <- ggplot(data = plot.df) + 
  geom_boxplot(aes(x = TF, y = EffectScore,
                   fill = factor(group),
                   #alpha = 0.8
  )) + 
  scale_fill_manual(values = mycol)+
  scale_color_manual(values = mycol)+
  theme_classic()+
  theme(legend.position ='none',
        axis.text = element_text(color = 'black'),
        axis.text.x = element_text(hjust = 1,angle = 90)
  )+ geom_hline(yintercept = 0,lty=3,col="black",lwd=0.5)
#+coord_flip()

ggsave('/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/DepMap_candiTF_all.pdf',
       p.flip,width = 21,height = 7,units = 'cm')

#### sTable2. sigMR ####
mrs <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/ssMARINA/All/ALL_TCA cycle.RDS')
mrs.es <- do.call(data.frame,mrs$es[1:4])
mrs.es$fdr <- p.adjust(mrs.es$p.value,method = 'fdr')

mrs.es$group='notSig'
mrs.es$group[mrs.es$nes<0 & mrs.es$fdr <0.05] <- 'Repression' #18
mrs.es$group[mrs.es$nes>0 & mrs.es$fdr <0.05] <- 'Activation' #20

mrs.es$log10FDR=-log10(mrs.es$fdr)
write.csv(mrs.es,'/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/allMR.csv')


sig <- mrs.es[mrs.es$group !='notSig',] 
sig <- sig[order(abs(sig$nes),decreasing = T),]
write.csv(sig,'/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/sigMR.csv')

