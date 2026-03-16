rm(list=ls())
library(GSVA)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(clusterProfiler)

# ===================== 1. load data =====================
expr_matrix <- read.csv('/remote-home/yanzijun/CRU/ped_M/res/celline/CCLE/RNAseq/exp_ALL.csv',row.names = 1)

gene.set <- read.gmt('/remote-home/yanzijun/CRU/TALL_M/data/GeneSet7.gmt')
gene.set <- gene.set$gene[gene.set$term=='TCA cycle']
tca_gene_set <- list("TCA_cycle" = gene.set)

# ===================== 2. run ssGSEA =====================
ssgsea_result <- gsva(
  expr = as.matrix(expr_matrix),
  gset.idx.list = tca_gene_set,
  method = "ssgsea", 
  kcdf = "Gaussian",
  abs.ranking = TRUE,
  ssgsea.norm = TRUE
)


ssgsea_score <- as.data.frame(t(ssgsea_result)) %>% 
  tibble::rownames_to_column("CellLine")
ssgsea_score <- ssgsea_score[order(ssgsea_score$TCA_cycle,decreasing = T),]
ssgsea_score <- ssgsea_score %>%
  mutate(
    q1 = quantile(TCA_cycle, 0.25, na.rm = TRUE),  
    q3 = quantile(TCA_cycle, 0.75, na.rm = TRUE),  
    state = case_when(
      TCA_cycle >= q3 ~ "H-A",          
      TCA_cycle <= q1 ~ "L-A",          
      TRUE ~ "M-A"                 
    )
  ) %>%
  dplyr::select(-q1, -q3) 

## sort
df_sorted <- ssgsea_score %>%
  arrange(TCA_cycle) %>%  
  mutate(CellLine = factor(CellLine, levels = .$CellLine))

# ===================== 3. plot =====================
p1 <- ggplot(df_sorted, aes(x = CellLine, y = TCA_cycle,fill=state)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("L-A" = "#9cb0c3", "M-A" = "#d6d6d6", "H-A" = "#eab080")) +
  labs(x = "Cell Line",y = "Enrichment Score", title = "") +
  theme_classic()+ 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,colour = 'black'), 
    legend.position = "none"
  )
ggsave("/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/Fig4_TCA_ssGSEA_barplot.pdf", p1, width = 7, height = 3, dpi = 300)

# ===================== 4. sTable1 =====================
write.csv(ssgsea_score, "/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/sTable_cellinestate.csv", row.names = FALSE)

