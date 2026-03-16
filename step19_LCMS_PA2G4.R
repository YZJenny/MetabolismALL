rm(list=ls())
library(ggplot2)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(ggbreak)
library(ggsignif)
library(readr)
setwd("/remote-home/yanzijun/CRU/ped_M/res_FDR0.25")
print(getwd())

celline <- 'REH'
df <- read_csv(paste("LCMS/",celline,"_TCAcycle.csv",sep=''))
colnames(df) <- gsub(' ','',colnames(df))
colnames(df) <- gsub('\\(ug/ml\\)','',colnames(df))
df$group <- gsub('-[1-6]','',df$SampleName)
df <- df[df$group != 'QC',]
table(df$group)

map <- c('CIR','cisAA','ICIR','OXG','SUC','FUM','MAL')
names(map) <- c("Citric_acid","Cis-aconitic_acid", "Isocitric_acid","a-Ketoglutaric_acid","Succinic_acid","Fumaric_acid","Malic_acid")
if(celline=='PEER'){
  df$Concentration <- log(df$Concentration+1,2)
  df$Metabolism <- apply(as.matrix(df$CompoundName),1,function(x) unlist(strsplit(x,'  \\['))[1]) 
  df$Metabolism[df$Metabolism=='Cis-Aconitic acid'] <- 'cisAA'
  df$Metabolism <- factor(df$Metabolism,levels = c('CIR','cisAA','ICIR','OXG','SUC','FUM','MAL'))
}else if(celline=='REH'){
  df$Metabolism <- map[df$CompoundName]
  df$Metabolism <- factor(df$Metabolism,levels = map)
}


df$group <- factor(df$group,levels = paste(celline,c('CTR','VEN','WS6','VEN-WS6'),sep='-'))
mycol= c('#63A8D2','#E88482')

mean <- aggregate(df$Concentration, by=list(df$group, df$Metabolism), FUN=mean)
sd <- aggregate(df$Concentration, by=list(df$group, df$Metabolism), FUN=sd)
len <- aggregate(df$Concentration, by=list(df$group, df$Metabolism), FUN=length)

df_res <- data.frame(mean, sd=sd$x, len=len$x)
colnames(df_res) = c("group", "Metabolism", "Mean", "Sd", "Count")
str(df_res)
df_res$Se <- df_res$Sd/sqrt(df_res$Count) 


### CTRL vs WS6
plot.data <- df_res[df_res$group %in% paste(celline,c('CTR','WS6'),sep='-'),]
plot.data$group <- factor(plot.data$group,levels = paste(celline,c('CTR','WS6'),sep='-'))

p1 <- ggplot(plot.data, aes(x=Metabolism, y=Mean, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(),color="black", width=.8) +
  scale_fill_manual(values = mycol) +
  geom_errorbar(aes(ymin=Mean-Sd, ymax=Mean +Sd),position=position_dodge(.8), width=.4) +
  theme_classic()+
  labs(y='Concentration (ug/ml)')+
  theme(axis.text = element_text(color='black',size=10))
ggsave(paste('LCMS/',celline,'_TCAcycle_WS6_1.pdf',sep=''),p1,width = 5,height = 3)


### VEN vs WS6+VEN
plot.data <- df_res[df_res$group %in% paste(celline,c('CTR','WS6'),sep='-'),]
p1 <- ggplot(plot.data, aes(x=Metabolism, y=Mean, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(),color="black", width=.8) +
  scale_fill_manual(values = mycol) +
  geom_errorbar(aes(ymin=Mean-Sd, ymax=Mean +Sd),position=position_dodge(.8), width=.4) +
  theme_classic()+
  labs(y='Concentration (ug/ml)')+
  theme(axis.text = element_text(color='black',size=10))
ggsave(paste('LCMS/',celline,'_TCAcycle_VEN_WS6.pdf',sep=''),p1,width = 5,height = 3)


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
  ws6 <- as.numeric(df$Concentration[df$Metabolism==c & df$group ==  paste(celline,'-WS6',sep='')])
  ven_ws6 <- as.numeric(df$Concentration[df$Metabolism==c & df$group ==  paste(celline,'-VEN-WS6',sep='')])
  
  if(type=='1'){
    p <- t.test(ctr,ws6,alternative = 'greater')$p.value
  }else if(type=='2'){
    p <- t.test(ven,ven_ws6,alternative = 'greater')$p.value
  }
  
  signif <- print_stars(p)
  print(signif)
}