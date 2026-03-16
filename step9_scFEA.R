##### scFEA using node9
rm(list=ls())
library(Seurat)
library(dplyr)
setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/')


#### 1. make input csv file ####
pbmc <- readRDS('pbmc.RDS')

## sample cells
Idents(pbmc) <- pbmc$celltype_3
sub.pbmc <- subset(pbmc,downsample=5000) 
data <- as.matrix(sub.pbmc@assays$RNA@data)
write.csv(data,'../scFEA/ALL_sc_sample5k.csv',row.names = T)
meta <- sub.pbmc@meta.data
write.csv(meta,'../scFEA/ALL_meta_sample5k.csv',row.names = T)


#### 2. move to node9 /local/yanzijun/software/scFEA-master/input_addZou ####
# cd ~/software/scFEA-master/
# 
# nohup python src/scFEA.py --data_dir data --input_dir input_addZou \
# --test_file ALL_sc_sample5k.csv \
# --moduleGene_file module_gene_m168.csv \
# --stoichiometry_matrix cmMat_c70_m168.csv \
# --res_dir output_addZou \
# --sc_imputation True &