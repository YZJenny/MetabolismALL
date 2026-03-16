#### Fig5A-B.network ####
rm(list=ls())
library(igraph)
library(ggraph)
library(ssmarina)
library(mixtools)
library(parallel)
library(tidygraph)
library(ggrepel)

## 1. load sigMR-target
mrs <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_FDR0.25/ssMARINA/All/ALL_TCA cycle.RDS')
mrs.es <- read.csv('/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/allMR.csv',header = T,row.names = 1)
sig <- read.csv('/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/sigMR.csv',header = T,row.names = 1)

TFs <- rownames(sig)
TF_target <- c()
for(tf in TFs){
  df <- data.frame(TF=tf,target=names(mrs$regulon[[tf]]$tfmode),tfmode=mrs$regulon[[tf]]$tfmode)
  df$ledge <-'No'
  df$ledge[df$target %in% mrs$ledge[[tf]]] <- 'Yes'
  
  TF_target <- rbind(TF_target,df)
}

write.csv(TF_target,'/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/network_ALL_TCA cycle.csv',
          row.names = FALSE,quote = FALSE)

## 2. plot
nes2TF <- abs(sig$nes)
names(nes2TF) <- rownames(sig)

get_network <- function(show.TF,col){
  plot_data <- read.csv('/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/network_ALL_TCA cycle.csv')
  plot_data <- plot_data[plot_data$ledge=='Yes',]
  plot_data$ledge=NULL
  plot_data <- plot_data[plot_data$TF %in% show.TF,]
  
  mat_weight <- reshape2::dcast(plot_data,TF ~ target, value.var = 'tfmode')
  dim(mat_weight)
  
  rownames(mat_weight) <- mat_weight$TF
  mat_weight$TF <- NULL
  
  row_add <- setdiff(colnames(mat_weight), rownames(mat_weight))
  df.row_add <- data.frame(matrix(rep(NA, length(row_add)*ncol(mat_weight)), 
                                  nrow = length(row_add), ncol = ncol(mat_weight)),
                           row.names = row_add)
  colnames(df.row_add) <- colnames(mat_weight)
  mat_weight <- rbind(mat_weight, df.row_add)
  col_add <- setdiff(rownames(mat_weight), colnames(mat_weight))
  df.col_add <- data.frame(matrix(rep(NA, length(col_add)*nrow(mat_weight)), 
                                  nrow = nrow(mat_weight), ncol = length(col_add)),
                           row.names = rownames(mat_weight))
  colnames(df.col_add) <- col_add
  mat_weight <- cbind(mat_weight, df.col_add)
  mat_weight[is.na(mat_weight)] <- 0
  
  TFgene <- sort(colnames(mat_weight))
  mat_weight <- mat_weight[TFgene, TFgene]
  
  igraph_mtx <- 
    graph_from_adjacency_matrix(t(as.matrix(mat_weight)), 
                                mode = 'directed', 
                                weighted = T, diag = T)
  V(igraph_mtx)$degree <- degree(igraph_mtx, normalized = T)
  V(igraph_mtx)$weight_degree <- strength(igraph_mtx)
  V(igraph_mtx)$page_rank <- page_rank(igraph_mtx)$vector

  nes <- rep(0.1,nrow(mat_weight))
  names(nes) <- names(strength(igraph_mtx))
  nes[show.TF] <- nes2TF[show.TF]
  V(igraph_mtx)$nes <-nes
  
  node_list <- data.frame(
    node_id = V(igraph_mtx)$name,
    degree = V(igraph_mtx)$degree,
    weight_degree = V(igraph_mtx)$weight_degree,
    page_rank = V(igraph_mtx)$page_rank,
    nes <- V(igraph_mtx)$nes)
  
  ## set nodes color
  node_type <- V(igraph_mtx)$name
  node_group <- rep('2', length(node_type))
  node_group[node_type %in% show.TF] <- '1' 
  V(igraph_mtx)$group <- node_group
  
  #saveRDS(igraph_up,'ssMARINA/FDR0.1/network.RDS')
  
  ggraph <- as_tbl_graph(igraph_mtx)
  
  p <- 
    ggraph(ggraph, layout = 'stress') + 
    geom_edge_link(aes(edge_width=weight, alpha = weight),color="gray97",
                   arrow = arrow(length = unit(0, 'mm')), 
                   end_cap = circle(1, 'mm'), linejoin = 'bevel', 
                   start_cap = circle(0.3, 'mm')) +
    scale_edge_width(range=c(0.6,1)) + 
    scale_size_continuous(range = c(3,10)) + 
    geom_node_point(aes(size = page_rank,fill = group, color = group),#
                    shape=21)+
    scale_color_manual(values = c(col, 'gray')) +
    scale_fill_manual(values = c(col, 'gray')) +
    # scale_alpha_manual(values = c( rep(1,length(vec_desc)+2), 0.1)) +
    geom_node_text(aes(filter = group == 1,label=name),size=6, repel = T) +
    geom_node_text(aes(filter = group == 2,label=name),size=0) +
    theme_void() + theme(legend.position = 'none')
  return(p)
}


show.TF <- rownames(sig)[sig$nes>0]
print(show.TF)
p_High <- get_network(show.TF = show.TF,col='#E88482')
ggsave(plot = p_High, path = '/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/',
       filename = 'pathway_TFnet_High.pdf',height = 17, width = 17)

show.TF <-  rownames(sig)[sig$nes<0]
print(show.TF)
p_Low <- get_network(show.TF = show.TF,col='#63A8D2')
ggsave(plot = p_Low, path = '/remote-home/yanzijun/CRU/ped_M/res_ALL/candiTF/',
       filename = 'pathway_TFnet_Low.pdf',height = 14, width = 15)

#### Fig5C. ssmarina plot ####
rm(list=ls())
library(corto)
library(parallel)

cancer.lst <- c("ALL")
pathway.lst <- c("TCA cycle")

library(plotrix)
source('/remote-home/yanzijun/CRU/ped_M/scr/cortoplot_myself_v3.R')

for(i in 1:length(cancer.lst)){
  cancer=cancer.lst[i]
  print(cancer)
  for(j in 1:length(pathway.lst)){
    pathway=pathway.lst[j]
    print(pathway)
    
    input <- paste('ssMARINA/All/',cancer,'_',pathway,'.RDS',sep='')
    if(file.exists(input)){
      mrs <- readRDS(input)
      
      ## marina to corto object
      new_mrs <- list()
      new_mrs$nes <- mrs$es$nes
      new_mrs$pvalue <- mrs$es$p.value
      new_mrs$sig <- as.numeric(mrs$signature)
      names(new_mrs$sig) <- rownames(mrs$signature)
      new_mrs$regulon <- mrs$regulon
      
      pdf(paste('ssMARINA/cortoplot/',cancer,'_',pathway,'_top10.pdf',sep=''),width = 8,height = 7)
      cortoplot_myself(new_mrs,10)
      dev.off()
    }
  }
}


for(i in 1:length(cancer.lst)){
  cancer=cancer.lst[i]
  print(cancer)
  for(j in 1:length(pathway.lst)){
    pathway=pathway.lst[j]
    print(pathway)
    
    input <- paste('ssMARINA/All/',cancer,'_',pathway,'.RDS',sep='')
    if(file.exists(input)){
      mrs <- readRDS(input)
      pdf(paste('ssMARINA/ssMARINAplot/',cancer,'_',pathway,'_top10.pdf',sep=''),
          width = 7,height = 7)
      plot(x=mrs,10,cex=0.7)
      dev.off()
    }
  }
}
