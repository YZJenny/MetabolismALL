#### sFig6A. from HPA ####

#### sFig6B.  survplot of PA2G4/YBX1 in TARGET ####
rm(list=ls())
library(dplyr)
library(survival)
library(survminer)
library(survivalROC)

## top TF
sig <- read.csv('/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/sigMR.csv',header = T,row.names = 1)
candiTF <- rownames(sig)[1:10]
print(candiTF)
candiTF <- candiTF[c(6:9,1:5,10)]
print(candiTF)

## clin. info.
clin <- read.csv('/remote-home/yanzijun/CRU/ped_M/data/RNAseq_clinical/TARGET_clinical.csv')
table(clin$project_id)

sub.clin <- clin[clin$vital_status %in% c('Dead','Alive'),
                 c('case_submitter_id','vital_status','days_to_last_follow_up')]
sub.clin$status <- ifelse(sub.clin$vital_status=='Dead',1,0)
sub.clin$days_to_last_follow_up <- as.numeric(sub.clin$days_to_last_follow_up)
sub.clin$time <- round(sub.clin$days_to_last_follow_up/12,2)
sub.clin$vital_status=sub.clin$days_to_last_follow_up=NULL
sub.clin <- na.omit(sub.clin)

##
cutoff.df <- read.csv('/remote-home/yanzijun/CRU/ped_M/res/progGene/cutpoint_ALL.csv',header = T)
exp.df <- na.omit(as.data.frame(t(readRDS('/remote-home/yanzijun/CRU/ped_M/res/RNAseq/exp_ALL.count.rds'))))
exp.df <- log(exp.df+1,2)
exp.df$case_submitter_id <-  apply(as.matrix(rownames(exp.df)),1,function(x) paste(unlist(strsplit(x,'\\.'))[1:3],collapse ='-'))

splots <- list()
for(i in 1:length(candiTF)){
  TF=candiTF[i]
  cutoff <- cutoff.df$cutpoint[cutoff.df$gene==TF]
  
  exp <- exp.df[,c('case_submitter_id',TF)]
  print(dim(exp))
  
  plot_data <- merge(sub.clin,exp,by='case_submitter_id')
  colnames(plot_data)[ncol(plot_data)] <- 'exp'
  plot_data$surtype <- 'High'
  plot_data$surtype[plot_data$exp <= cutoff] <- 'Low'
  print(table(plot_data$surtype))
  
  fit <- survminer::surv_fit(Surv(time = time, event = status) ~ surtype,data=plot_data)
  print(surv_pvalue(fit)$pval)
  p <- ggsurvplot(fit, data=plot_data,pval = TRUE,
                  palette = c('#E88482','#63A8D2'),
                  ggtheme = theme_survminer(base_size = 10),title=TF,
  ) 
  splots[[i]] <- p
}

p <- arrange_ggsurvplots(splots, ncol = 5, nrow = 1) 
ggsave('/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/survplot_candiTF.pdf',p,width = 12,height = 3)