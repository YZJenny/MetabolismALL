rm(list = ls())
library(pheatmap)
library(ggplot2)
library(readxl)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(dplyr)
setwd("/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/")

## load drug response data from nc2022
DS_raw <- as.data.frame(read_excel("celline/nc2022/sDSS.xlsx"))
rownames(DS_raw) <- DS_raw$...1
DS_raw$...1 <- NULL
print(colnames(DS_raw))
DS <- as.data.frame(na.omit(t(scale(t(DS_raw)))))

## load cell line subtype
RS='nc2022'
cohort='CCLE'
cancer <- "ALL"
pathway <- "TCA cycle"

subtype.df <- read.csv(paste("celline/",cohort,"/subtype_", cancer, ".csv", sep = ""), row.names = 1)
colnames(subtype.df) <- gsub("X", "", colnames(subtype.df))
colnames(subtype.df) <- gsub("\\.", "-", colnames(subtype.df))
colnames(subtype.df) <- toupper(colnames(subtype.df))
print(colnames(subtype.df))

ol <- intersect(colnames(subtype.df), colnames(DS))
subtype.df <- subtype.df[, ol]
print(colnames(subtype.df))

High <- colnames(subtype.df)[which(subtype.df[pathway, ] == "High")]
Middle <- colnames(subtype.df)[which(subtype.df[pathway, ] == "Middle")]
Low <- colnames(subtype.df)[which(subtype.df[pathway, ] == "Low")]
group <- c(rep("Low", length(Low)), rep("Middle", length(Middle)), rep("High", length(High)))
data <- select(DS,c(Low,Middle,High))


##### sFig7A ##### 
coldata <- data.frame(group=group);rownames(coldata)<- colnames(data)
coldata$group[coldata$group=='Middle'] <- 'Low'

Group <- coldata$group
design = model.matrix(~Group)
fit = lmFit(data,design)
fit = eBayes(fit)
pathway_res = topTable(fit,coef = 2,number = Inf)
pathway_res[1:5,1:5]

pathway_sig <- pathway_res[pathway_res$P.Value<0.05 & abs(pathway_res$logFC) > 0.2,]

library(pheatmap)
library(stringr)
sub.data <- data[match(rownames(pathway_sig),rownames(data)),]

sort_data <- order(Group)
n <- sub.data[,sort_data]
pd <- coldata[sort_data,]

# 调整颜色梯度
range(sub.data)
breaksList = seq(-3, 3, by = 0.1)
colors <- colorRampPalette(c("#336699", "white", "tomato"))(length(breaksList))

#创建列和行注释
annCol <- data.frame(group = Group,
                     row.names = colnames(sub.data),
                     stringsAsFactors = FALSE)
p <- pheatmap(sub.data,
              annotation_col = annCol,
              #annotation_row = annRow,
              color = colors,
              breaks = breaksList,
              cluster_rows = T,
              cluster_cols = FALSE,
              show_rownames = TRUE,
              show_colnames = TRUE,
              fontsize_row = 8,
              fontsize_col = 8,
              gaps_col = c(length(Low), length(Low) + length(Middle)),
              annotation_names_row = FALSE
)
pdf(paste("celline/",cohort,"/pheatmap_",RS,'_', cancer, "_", pathway, "_all.pdf",sep=''), height = 3.5, width = 5)
print(p)
dev.off()


#### Fig6A. dotplot of cell line response to Venetoclax ####
df <- data.frame(celline=colnames(data),sDSS = as.numeric(data['Venetoclax', ]), subtype = group)
df <- df[order(df$sDSS,df$subtype),]
df <- df[df$subtype != 'Middle',]
df$order <- 1:nrow(df)


p <- ggplot(df, aes(x =order , y = sDSS, col = subtype 
                    #, shape = subtype 
)) + #col 颜色分组， shape 形状分组， size 大小分组
  geom_point()+
  scale_color_manual(values = c("#63A8D2","#E88482"))+
  theme_classic()+
  theme(axis.text = element_text(colour = 'black'))
ggsave('celline/CCLE/dotplot_nc2022_Venetoclax_ALL.pdf',p,width = 4,height = 2)
