rm(list=ls())
library(ggplot2)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(readr)
library(ggbreak)

setwd("/remote-home/yanzijun/CRU/ped_M/res_FDR0.25")
print(getwd())

celline <- 'PEER'
df <- read_csv('LCMS/SU056_TCAcycle.csv')
colnames(df) <- gsub(' ','',colnames(df))
colnames(df) <- gsub('\\(ug/ml\\)','',colnames(df))
df$group <- gsub('-[1-6]','',df$SampleName)
df <- df[df$group != 'QC',]
table(df$group)
df <- df[grep(celline,df$group),]
table(df$group)

map <- c('CIR','cisAA','ICIR','OXG','SUC','FUM','MAL')
names(map) <- c("Citric_acid","Cis-aconitic_acid", "Isocitric_acid","a-Ketoglutaric_acid","Succinic_acid","Fumaric_acid","Malic_acid")
df$Metabolism <- map[df$CompoundName]
df$Metabolism <- factor(df$Metabolism,levels = map)

df$group <- factor(df$group,levels = paste(celline,c('CTR','VEN','SU056','VEN-SU056'),sep='-'))
mycol= c('#DCDDDD','#808080')

mean <- aggregate(df$Concentration, by=list(df$group, df$Metabolism), FUN=mean)
sd <- aggregate(df$Concentration, by=list(df$group, df$Metabolism), FUN=sd)
len <- aggregate(df$Concentration, by=list(df$group, df$Metabolism), FUN=length)

df_res <- data.frame(mean, sd=sd$x, len=len$x)
colnames(df_res) = c("group", "Metabolism", "Mean", "Sd", "Count")
str(df_res)
df_res$Se <- df_res$Sd/sqrt(df_res$Count) 

### CTRL vs SU056
plot.data <- df_res[df_res$group %in% paste(celline,c('-CTR','-SU056'),sep=''),]
p1 <- ggplot(plot.data, aes(x=Metabolism, y=Mean, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(),color="black", width=.8) +
  scale_fill_manual(values = mycol) +
  geom_errorbar(aes(ymin=Mean-Sd, ymax=Mean +Sd),position=position_dodge(.8), width=.4) +
  theme_classic()+
  labs(y='Concentration (ug/ml)')+
  theme(axis.text = element_text(color='black',size=10))

ggsave(paste('LCMS/',celline,'_TCAcycle_SU056.pdf',sep=''),p1,width = 5,height = 3)


### VEN vs SU056+VEN
plot.data <- df_res[df_res$group %in% paste(celline,c('-VEN','-VEN-SU056'),sep=''),]
p1 <- ggplot(plot.data, aes(x=Metabolism, y=Mean, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(),color="black", width=.8) +
  scale_fill_manual(values = mycol) +
  geom_errorbar(aes(ymin=Mean-Sd, ymax=Mean +Sd),position=position_dodge(.8), width=.4) +
  theme_classic()+
  labs(y='Concentration (ug/ml)')+
  theme(axis.text = element_text(color='black',size=10))
ggsave(paste('LCMS/',celline,'_TCAcycle_VEN_SU056.pdf',sep=''),p1,width = 5.5,height = 3)


## cal. p-value
print_stars <- function(p) {
  if (is.na(p)) {
    print("p值未提供或无法计算")
  } else {
    if (p < 0.0001) {
      cat("****\n")
    } else if (p < 0.001) {
      cat("***\n")
    } else if (p < 0.01) {
      cat("**\n")
    } else if (p < 0.05) {
      cat("*\n")
    } else {
      print("ns")
    }
  }
}

type='1'
Comp <- levels(df$Metabolism)
group <- unique(df$group)
for(i in 1:length(Comp)){
  c <- Comp[i]
  print(c)
  ctr <- as.numeric(df$Concentration[df$Metabolism==c & df$group == paste(celline,'-CTR',sep='')])
  ven <- as.numeric(df$Concentration[df$Metabolism==c & df$group ==  paste(celline,'-VEN',sep='')])
  SU056 <- as.numeric(df$Concentration[df$Metabolism==c & df$group ==  paste(celline,'-SU056',sep='')])
  ven_SU056 <- as.numeric(df$Concentration[df$Metabolism==c & df$group ==  paste(celline,'-VEN-SU056',sep='')])
  
  if(type=='1'){
    p <- t.test(ctr,SU056,alternative = 'greater')$p.value
  }else if(type=='2'){
    p <- t.test(ven,ven_SU056,alternative = 'greater')$p.value
  }
  
  signif <- print_stars(p)
  print(signif)
}