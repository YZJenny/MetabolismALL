# This script maps the leukemic cells to the T/B-lineage map. https://zenodo.org/records/14346457
{
  rm(list = ls())
  gc(full = T)
  library(Seurat)
  library(stringr)
  library(patchwork)
  library(anndata)
  library(Matrix)
  library(ggplot2)
  library(tidyverse)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(future)
  library(qs)
  library(reshape2)
  set.seed(42)
  source("/remote-home/yanzijun/CRU/TALL/data/STM_2025/scr/symphony_utils_seurat.R")
  source("/remote-home/yanzijun/CRU/TALL/data/STM_2025/scr/utils.R")
  setwd('/remote-home/yanzijun/CRU/ped_M/res_ALL')
}
set_output_dir("TB_atlas_mapping/")

{ # color coding
  ref_map_colors <- c("B Cells" = "#B3CDE3", "B Prog" = "#99CC99",
                      "DC and DC Prog" = "#7A5AA3", "Erythroid" = "#FDCEA1",
                      "HSPCs" = "#D43F27",
                      "Mature Myeloid" = "#0E8243", "Myeloid Prog" = "#D8E6B7",
                      "NK Cells" = "#4444AA", "pDC and pDC Prog" = "#A3ACD7",
                      "Thymocyte" = "#4389C8", "T Cells" = "#DECBE4",
                      "Plasma Cells" = "#6475B6") #"Mast Cells" = "#AC6FAD",
  tbcell_map_colors <- c("abT agonist" = "#D8E6B7", "abT entry" = "#6475B6", "CD4 T" = "#6DCCDA",
                        "CD4 Tmem" = "#F1615F", "CD8 T" = "#DECBE4",
                        "CD8 Tmem" = "#F7B6D2", "CD8aa T" = "#4444AA", "DN" = "#AC6FAD", "DP" = "#A3ACD7", "ETP" = "#FDCEA1",
                        "gdT" = "#7A5AA3", "HSPCs" = "#D43F27", "NKT" = "#98DF8A", "pDC" = "#B3CDE3",
                        "Treg" = "#4389C8","B Cells" = "#33CCCC", "B Prog" = "#99CC99")
  
  
}

#### 1. subsetted T-lineage cells from BMTH atlas ####
merged.ref.map <- qs::qread(file='/remote-home/yanzijun/CRU/TALL/data/STM_2025/merged.ref.map_final_266003_cells_cd4_cd8_removed_Harmony_integrated_annotated_noUMAP.qs',
                            nthreads = 70)
tbcell.ref.map_pDC <- subset(merged.ref.map, subset = Mapped_Annotation_high_lvl %in%
                              c("HSPCs",'B Prog', "Thymocyte",'B Cells', "T Cells", "pDC and pDC Prog"))
tbcell.ref.map_pDC_No_ILC_Cycling <- subset(tbcell.ref.map_pDC, cells = colnames(tbcell.ref.map_pDC)[!tbcell.ref.map_pDC$Mapped_Annotation_fine %in% c("ILC", "Cycling")])

tbcell.ref.map_pDC_No_ILC_Cycling$Mapped_Annotation_high_lvl <- factor(tbcell.ref.map_pDC_No_ILC_Cycling$Mapped_Annotation_high_lvl, 
                                                                      levels =c("HSPCs",'B Prog', "Thymocyte",'B Cells', "T Cells", "pDC and pDC Prog")) ##
qs::qsave(x = tbcell.ref.map_pDC_No_ILC_Cycling, file = file.path(RESULT_PATH, "tbcell.ref.map_pDC_No_ILC_Cycling_165101_of_266003.qs"),
          nthreads = 70)

## plot
# tbcell.ref.map <- qs::qread(file = file.path(RESULT_PATH,"tbcell.ref.map_pDC_No_ILC_Cycling_165101_of_266003.qs"),
#                            nthreads = 70)
tbcell.ref.map <- tbcell.ref.map_pDC_No_ILC_Cycling
rm(tbcell.ref.map_pDC_No_ILC_Cycling)
tbcell.ref.map$Mapped_Annotation_fine <- factor(tbcell.ref.map$Mapped_Annotation_fine,
                                               levels = c("HSPCs", 
                                                          'B Prog',
                                                          "ETP","DN", "DP",
                                                          'B Cells',
                                                          "abT entry", "abT agonist",
                                                          "CD4 T", "CD8 T",
                                                          "Treg",
                                                          "CD8aa T", "gdT", "CD4 Tmem", "CD8 Tmem",
                                                          "NKT",
                                                          "pDC")) ##
{
  pdf(paste0(RESULT_PATH, "subsetted_fine_annotations.pdf"), width = 18, height = 12)
  print(DimPlot(tbcell.ref.map, group.by = "Mapped_Annotation_high_lvl", label = T,
                raster = FALSE, label.size = 5, shuffle = T, pt.size = .1) +
          ggplot2::scale_color_manual(values = scales::alpha(ref_map_colors, .3)))
  print(DimPlot(tbcell.ref.map, group.by = "Mapped_Annotation_fine", label = T,
                raster = FALSE, label.size = 5, shuffle = T, pt.size = .1) +
          ggplot2::scale_color_manual(values =
                                        rgb(runif(length(unique(tbcell.ref.map$Mapped_Annotation_fine)), 0.2, 0.8),
                                            runif(length(unique(tbcell.ref.map$Mapped_Annotation_fine)), 0.2, 0.8),
                                            runif(length(unique(tbcell.ref.map$Mapped_Annotation_fine)), 0.2, 0.8))
          ))
  graphics.off()
}

##### 1.1 findmarkers #####
{
  plan(multisession, workers=55)
  options(future.globals.maxSize=1000*1024^3)
  Idents(tbcell.ref.map) <- tbcell.ref.map$Mapped_Annotation_fine
  lymphoid_DEGs <- FindAllMarkers(object = tbcell.ref.map, test.use = "wilcox",
                                  random.seed = 42,
                                  # latent.vars = "DonorID",
                                  assay = "RNA", verbose = T, min.pct = 0.1)
  
  lymphoid_DEGs$cluster <- factor(lymphoid_DEGs$cluster,
                                  levels = c("HSPCs", 
                                             'B Prog',
                                             "ETP","DN", "DP",
                                             'B Cells',
                                             "abT entry", "abT agonist",
                                             "CD4 T", "CD8 T",
                                             "Treg",
                                             "CD8aa T", "gdT", "CD4 Tmem", "CD8 Tmem",
                                             "NKT",
                                             "pDC")) ##
  
  lymphoid_DEGs <- lymphoid_DEGs[lymphoid_DEGs$p_val_adj < 0.05, ]
  lymphoid_DEGs_up <- lymphoid_DEGs[lymphoid_DEGs$avg_log2FC > 0, ]
  
  cell_types <- unique(lymphoid_DEGs_up$cluster)
  select.features <- c()
  cell_types
  top5.features <- c()
  cell_type_marker_list <- list()
  for (idx in 1:length(cell_types)){
    tmp.table <- lymphoid_DEGs_up[lymphoid_DEGs_up$cluster == cell_types[idx], ]
    tmp.table <- tmp.table[order(tmp.table$avg_log2FC, decreasing = T), ]
    cell_type_marker_list[[as.character(cell_types[idx])]] <- tmp.table$gene
    
    tmp.table <- head(tmp.table, 5)
    top5.features <- c(top5.features, tmp.table$gene)
  }
  top5.features <- unique(top5.features)
  print(length(top5.features))
}
saveRDS(lymphoid_DEGs,paste0(RESULT_PATH,'lymphoid_DEGs.RDS'))

{
  pdf(paste0(RESULT_PATH, "top5_features_dotplot.pdf"), width = 8, height = 18)
  print(DotPlot(tbcell.ref.map, features = top5.features, group.by = "Mapped_Annotation_fine") +
          coord_flip() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))
  graphics.off()
}

##### 1.2 re-UMAP accroding DEG from findmk #####
{
  variable_features <- unique(lymphoid_DEGs_up$gene)
  length(variable_features)
  VariableFeatures(tbcell.ref.map) <- variable_features
  tbcell.ref.map <- ScaleData(tbcell.ref.map)
  tbcell.ref.map <- RunPCA(tbcell.ref.map, verbose = T, npcs = 50)
  tbcell.ref.map <- RunHarmony.Seurat(tbcell.ref.map, group.by.vars = "DonorID")
  tbcell.ref.map[["umap"]] <- RunUMAP2(Embeddings(tbcell.ref.map, 'harmony')[, 1:50], assay = "RNA", verbose = T,
                                      umap.method='uwot', return.model=TRUE)
}
qs::qsave(x = tbcell.ref.map, file = file.path(RESULT_PATH, "tbcell.ref.map_pDC_No_ILC_Cycling_165101_of_266003_reUMAP.qs"), nthreads = 80)

##plot
{
  pdf(paste0(RESULT_PATH, "tbcell_atlas.pdf"), width = 18, height = 12)
  print(DimPlot(tbcell.ref.map, group.by = "Mapped_Annotation_high_lvl", label = T,
                raster = FALSE, label.size = 5, shuffle = T, pt.size = .1) +
          ggplot2::scale_color_manual(values = scales::alpha(ref_map_colors, .3)))
  
  print(DimPlot(tbcell.ref.map, group.by = "Mapped_Annotation_fine", label = T,
                raster = FALSE, label.size = 5, shuffle = T, pt.size = 0.1) +
          ggplot2::scale_color_manual(values = scales::alpha(tbcell_map_colors, .3)))
  graphics.off()
  
}

#### 2. run Symphony ####
tbcell.ref.map <- qread(file.path(RESULT_PATH, "tbcell.ref.map_pDC_No_ILC_Cycling_165101_of_266003_reUMAP.qs"), nthreads = 80)
print(tbcell.ref.map@reductions$umap@misc)
# tbcell.ref.map_subsetted.features_harmonyRef <- buildReferenceFromSeurat(tbcell.ref.map, verbose = TRUE, 
#                                                                          save_umap = TRUE,
#                                                                          save_uwot_path = paste0('/remote-home/yanzijun/CRU/ped_M/res_ALL/',RESULT_PATH,'/cache_symphony2.uwot'))
# qs::qsave(tbcell.ref.map_subsetted.features_harmonyRef,
#           file = file.path(RESULT_PATH, "tbcell.ref.map_subsetted.features_harmonyRef.qs"), nthreads = 80)

tbcell.ref.map_subsetted.features_harmonyRef <- qread(file.path(RESULT_PATH, "tbcell.ref.map_subsetted.features_harmonyRef.qs"), nthreads = 80)

in_house_samples_leukemic_subset <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/Seurat/pbmc_tumor_metab.RDS')
new_meta <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/TCAsubtype_new/subtype_info_kmeans.RDS')
print(all(colnames(in_house_samples_leukemic_subset)==rownames(new_meta)))
in_house_samples_leukemic_subset$subtype <- new_meta$subtype

in_house_samples_leukemic_subset_mapped <- mapQuery(in_house_samples_leukemic_subset@assays$RNA@counts,
                                                    in_house_samples_leukemic_subset@meta.data,
                                                    tbcell.ref.map_subsetted.features_harmonyRef,
                                                    return_type = 'Seurat')

{
  in_house_samples_leukemic_subset_mapped <- knnPredict.Seurat(in_house_samples_leukemic_subset_mapped,
                                                               tbcell.ref.map_subsetted.features_harmonyRef,
                                                               'Mapped_Annotation_high_lvl', confidence = TRUE)
  in_house_samples_leukemic_subset_mapped <- knnPredict.Seurat(in_house_samples_leukemic_subset_mapped,
                                                               tbcell.ref.map_subsetted.features_harmonyRef,
                                                               'Mapped_Annotation_fine', confidence = TRUE)
  
  in_house_samples_leukemic_subset_mapped$Mapped_Annotation_high_lvl <- factor(in_house_samples_leukemic_subset_mapped$Mapped_Annotation_high_lvl,
                                                                               levels = c("HSPCs",
                                                                                          'B Prog',"Thymocyte",
                                                                                          'B Cells',"T Cells",
                                                                                          "pDC and pDC Prog"))
  in_house_samples_leukemic_subset_mapped$Mapped_Annotation_fine <- factor(in_house_samples_leukemic_subset_mapped$Mapped_Annotation_fine,
                                                                           levels = c("HSPCs", 
                                                                                      'B Prog',
                                                                                      "ETP","DN", "DP",
                                                                                      'B Cells',
                                                                                      "abT entry", "abT agonist",
                                                                                      "CD4 T", "CD8 T",
                                                                                      "Treg",
                                                                                      "CD8aa T", "gdT", "CD4 Tmem", "CD8 Tmem",
                                                                                      "NKT",
                                                                                      "pDC"))#"NK", "ILC", "Cycling ILC", "Cycling TNK", "Cycling NK"
}
qs::qsave(in_house_samples_leukemic_subset_mapped,
          file = file.path(RESULT_PATH, "in_house_samples_leukemic_subset_mapped.qs"), nthreads = 70)

## merge ref and query
{
  seurat_obj_overlay_tcell_map <- merge(x = tbcell.ref.map,
                                        y = in_house_samples_leukemic_subset_mapped, merge.dr = "umap", add.cell.ids = c("ref", "query"))
  seurat_obj_overlay_tcell_map$Mapped_Annotation_high_lvl <- factor(seurat_obj_overlay_tcell_map$Mapped_Annotation_high_lvl,
                                                                               levels = c("HSPCs",
                                                                                          'B Prog',"Thymocyte",
                                                                                          'B Cells',"T Cells",
                                                                                          "pDC and pDC Prog"))
  seurat_obj_overlay_tcell_map$Mapped_Annotation_fine <- factor(seurat_obj_overlay_tcell_map$Mapped_Annotation_fine,
                                                                           levels = c("HSPCs", 
                                                                                      'B Prog',
                                                                                      "ETP","DN", "DP",
                                                                                      'B Cells',
                                                                                      "abT entry", "abT agonist",
                                                                                      "CD4 T", "CD8 T",
                                                                                      "Treg",
                                                                                      "CD8aa T", "gdT", "CD4 Tmem", "CD8 Tmem",
                                                                                      "NKT",
                                                                                      "pDC"))#"NK", "ILC", "Cycling ILC", "Cycling TNK", "Cycling NK"
  
}

qs::qsave(seurat_obj_overlay_tcell_map,file = file.path(RESULT_PATH, "seurat_obj_overlay_tbcell_map.qs"), nthreads = 70)

## plot
seurat_obj_overlay_tcell_map <- qs::qread(paste0(RESULT_PATH,'seurat_obj_overlay_tbcell_map.qs'),nthreads = 70)
new_meta <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/TCAsubtype_new/subtype_info_kmeans.RDS')

High_cells <- list("High" = paste0('query_',rownames(new_meta)[new_meta$subtype %in% "High"] ))
Medium_cells <- list("Medium" = paste0('query_',rownames(new_meta)[new_meta$subtype %in% "Medium"] ))
Low_cells <- list("Low" = paste0('query_',rownames(new_meta)[new_meta$subtype %in% "Low"] ))

{
  pdf(file = file.path(RESULT_PATH, "overlay_tcell_map.pdf"), width = 8, height = 6)
  print(DimPlot(seurat_obj_overlay_tcell_map,group.by = "Mapped_Annotation_high_lvl", label = T,
                raster = FALSE, label.size = 0, shuffle = T, pt.size = .1) +
          ggplot2::scale_color_manual(values = scales::alpha(ref_map_colors, .3))) 
  print(DimPlot(seurat_obj_overlay_tcell_map,group.by = "Mapped_Annotation_fine", label = T,
                raster = FALSE, label.size = 0, shuffle = T, pt.size = 0.1) +
          ggplot2::scale_color_manual(values = scales::alpha(tbcell_map_colors, .3)))
  graphics.off()
}

{
  pdf(file = file.path(RESULT_PATH, "cells_mapped_to_tbcell_atlas.pdf"), width = 6, height = 4)
  print(DimPlot(seurat_obj_overlay_tcell_map,
                raster = FALSE, label.size = 2, shuffle = T, pt.size = .1, cols = scales::alpha("gray90", .3),
                cells.highlight = High_cells, sizes.highlight = .1, cols.highlight = scales::alpha("#F41117", .2)) 
  )
  print(DimPlot(seurat_obj_overlay_tcell_map,
                raster = FALSE, label.size = 2, shuffle = T, pt.size = .1, cols = scales::alpha("gray90", .3),
                cells.highlight = Medium_cells, sizes.highlight = .1, cols.highlight = scales::alpha("#CCCC66", .2)) 
  )
  print(DimPlot(seurat_obj_overlay_tcell_map,
                raster = FALSE, label.size = 2, shuffle = T, pt.size = .1, cols = scales::alpha("gray90", .3),
                cells.highlight = Low_cells, sizes.highlight = .1, cols.highlight = scales::alpha("#66CCCC", .2)) 
  )
  graphics.off()
}
rm(seurat_obj_overlay_tcell_map)
gc(full = T)

#### 3. plot for statistic ####
in_house_samples_leukemic_subset_mapped <- qs::qread(paste0(RESULT_PATH,'in_house_samples_leukemic_subset_mapped.qs'),
                                                     nthreads = 70)
new_meta <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/TCAsubtype_new/subtype_info_kmeans.RDS')

print(all(colnames(in_house_samples_leukemic_subset_mapped)==rownames(new_meta)))
in_house_samples_leukemic_subset_mapped$subtype <- factor(new_meta$subtype,
                                                          levels = c('High','Medium','Low'))

##### 3.1 barplot of ratio of celltype ####
{
  pdf(file.path(RESULT_PATH, "composition_plots_high_lvl_annotation.pdf"), width = 4, height = 4)
  df.plot <- in_house_samples_leukemic_subset_mapped@meta.data[, c("subtype", "Mapped_Annotation_high_lvl", "Mapped_Annotation_high_lvl_prob")]
  df.plot <- df.plot %>%
    group_by(Mapped_Annotation_high_lvl, subtype, .drop = F) %>%
    summarize(Freq = length(Mapped_Annotation_high_lvl_prob)) %>%
    group_by(subtype) %>%
    mutate(Freq_in_Sample = Freq/sum(Freq))
  print(
    ggplot(df.plot, aes(x = subtype, y = Freq_in_Sample, fill = Mapped_Annotation_high_lvl)) +
      geom_bar(stat="identity", width = 0.8) +
      scale_y_continuous(labels = scales::percent) + ylab("Precentage of composition") + xlab("") +
      theme_minimal() + labs(fill = "") +
      scale_fill_manual(values = alpha(ref_map_colors, alpha = 1))+
      theme(axis.text = element_text(colour = 'black'))
    #theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  graphics.off()
}

{
  pdf(file.path(RESULT_PATH, "composition_plots_fine_lvl_annotation.pdf"), width = 4, height = 4)
  df.plot <- in_house_samples_leukemic_subset_mapped@meta.data[, c("subtype", "Mapped_Annotation_fine", "Mapped_Annotation_fine_prob")]
  df.plot <- df.plot %>%
    group_by(Mapped_Annotation_fine, subtype, .drop = F) %>%
    summarize(Freq = length(Mapped_Annotation_fine_prob)) %>%
    group_by(subtype) %>%
    mutate(Freq_in_Sample = Freq/sum(Freq))
  
  print(
    ggplot(df.plot, aes(x = subtype, y = Freq_in_Sample, fill = Mapped_Annotation_fine)) +
      geom_bar(stat="identity", width = 0.8) +
      scale_y_continuous(labels = scales::percent) + ylab("Precentage of composition") + xlab("") +
      theme_minimal() + labs(fill = "") +
      scale_fill_manual(values = alpha(tbcell_map_colors, alpha = 1)) +
      theme(axis.text = element_text(colour = 'black'))
    #theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  graphics.off()
}

##### 3.2 boxplot of TCA cycle score #####
irGSEA <- readRDS('irGSEA/pbmc_all_superpathway.RDS')
ssgsea <- data.frame(cell=colnames(irGSEA$assays$ssgsea@counts),
                     ssgsea=log(as.numeric(irGSEA$assays$ssgsea@counts['TCA cycle',]),10))


{
  pdf(file.path(RESULT_PATH, "ssGSEA_boxplots_high_lvl_annotation.pdf"), width = 3, height = 4)
  celltype <- data.frame(cell=colnames(in_house_samples_leukemic_subset_mapped),
                         celltype=in_house_samples_leukemic_subset_mapped$Mapped_Annotation_high_lvl)
  df.plot <- merge(ssgsea,celltype)
  print(
    ggplot(df.plot,aes(x=celltype,y=ssgsea,fill=celltype)) +
      geom_boxplot() + geom_jitter(width = 0,alpha = 0)+
      ylab("TCA cycle score") + xlab("") +
      theme_minimal() + labs(fill = "") +
      scale_fill_manual(values = alpha(ref_map_colors, alpha = 1)) +
      theme(axis.text = element_text(colour = 'black'),
            legend.position = 'none',
            axis.text.x = element_text(angle = 45, hjust = 1))
  )
  
  graphics.off()
}

{
  pdf(file.path(RESULT_PATH, "ssGSEA_boxplots_fine_lvl_annotation.pdf"), width = 5, height = 3)
  celltype <- data.frame(cell=colnames(in_house_samples_leukemic_subset_mapped),
                         celltype=in_house_samples_leukemic_subset_mapped$Mapped_Annotation_fine)
  df.plot <- merge(ssgsea,celltype)
  
  print( df.plot %>% ggplot(aes(x=celltype,y=ssgsea,fill=celltype)) +
           geom_boxplot() + geom_jitter(width = 0,alpha = 0)+
           ylab("TCA cycle score") + xlab("") +
           theme_minimal() + labs(fill = "") +
           scale_fill_manual(values = alpha(tbcell_map_colors, alpha = 1)) +
           theme(axis.text = element_text(colour = 'black'),
                 legend.position = 'none')+
           theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  graphics.off()
}

##### 3.3 heatmap: gene exp of mk of subtypes (不太符合) #####
in_house_samples_leukemic_subset_mapped <- qs::qread(paste0(RESULT_PATH,'in_house_samples_leukemic_subset_mapped.qs'),nthread=70)
mk <- readRDS('/remote-home/yanzijun/CRU/ped_M/res_ALL/TCAsubtype_new/mk_kmeans.RDS') 
top_mk <- mk %>% group_by(cluster) %>% top_n(n=20,wt=avg_log2FC)

mklst <-  split(top_mk$gene,top_mk$cluster)
pbmc <- AddModuleScore(in_house_samples_leukemic_subset_mapped,
                       features = mklst,name = "mk")

{
  pdf(file.path(RESULT_PATH, "Heatmap_high_lvl_annotation.pdf"), width = 5, height = 3)
  plot.df <- pbmc@meta.data[,c('Mapped_Annotation_high_lvl','mk1','mk2','mk3')]
  plot.df <- aggregate(. ~ Mapped_Annotation_high_lvl, data = plot.df, median)
  print(plot.df)
  
  melted_cormat <- melt(plot.df)
  
  print(ggplot(data = melted_cormat, aes(x=Mapped_Annotation_high_lvl, y=variable, fill=value)) + 
          geom_tile(color = "white")+
          scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                               midpoint = 0, limit = c(-2,2), space = "Lab", 
                               name=" ") +
          theme_minimal()+ 
          theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                           size = 12, hjust = 1))+
          coord_fixed()
  )
  graphics.off()
}


{
  pdf(file.path(RESULT_PATH, "Heatmap_fine_lvl_annotation.pdf"), width = 6, height = 3)
  plot.df <- pbmc@meta.data[,c('Mapped_Annotation_fine','mk1','mk2','mk3')]
  plot.df <- aggregate(. ~ Mapped_Annotation_fine, data = plot.df, median)
  print(plot.df)
  
  melted_cormat <- melt(plot.df)
  
  print(ggplot(data = melted_cormat, aes(x=Mapped_Annotation_fine, y=variable, fill=value)) + 
          geom_tile(color = "white")+
          scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                               midpoint = 0, limit = c(-2,2), space = "Lab", 
                               name=" ") +
          theme_minimal()+
          theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                           size = 12, hjust = 1))+
          coord_fixed()
  )
  graphics.off()
}


