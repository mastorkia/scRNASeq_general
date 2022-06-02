#!/usr/bin/env Rscript
message("Starting job...")

#.libPaths("/share/home/maider/R_libraries/R/x86_64-redhat-linux-gnu-library/4.1")
#Uncomment this if using in the server


message("Loading required libraries")

#Libraries needed for the analysis
library(dplyr)
library(Seurat)
library(patchwork)
library(plyr)
library(data.table)
library(ggplot2)
library(Matrix)
library(config)
library(matrixStats)
#require(knitr)
#require(markdown)
print(Sys.time())

rm(list = ls())
#Read the config file
config <- config::get(file= "config_seurat.yml")

message("**Loading Integrated data**")

#Load the craeted Seurat object.
load(paste0(config$seurat_object,config$seurat_object_dim))

message("**Running marker identification for each cluster..**")

#Set default assay to normalized data: RNA
DefaultAssay(Integrated_data) <- "RNA"

Markers <- FindAllMarkers(Integrated_data, min.pct = config$min.pct, logfc.threshold = config$log.fc_threshold_markers)

message(paste0("**Saving markers results in ",config$integration_folder,config$markers_outfile),"**")

#Save results
write.table(Markers, paste0(config$integration_folder,config$markers_outfile), quote = F, row.names = F, col.names = T, sep = ",")

#Plot top markers in a heatmap
Markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

png(file =paste0(config$integration_folder,"Top_markers_heatmap.png"),  
width = 1600,
height = 1600)
message("**Saving plot of Heatmap with top marker expression for each cluster**")

Heatmap <- DoHeatmap(Integrated_data, features = top10$gene, size = 8, assay="integrated") + NoLegend()
print(Heatmap)
dev.off()

#Plot cell type markers in vln plot and heatmap. We add the markers information from config file:
png(file =paste0(config$integration_folder,"Cell_type_markers_heatmap.png"),  
width = 1600,
height = 1600)
message("**Saving plot of Heatmap with cell type markers expression**")

Heatmap2 <- DoHeatmap(Integrated_data, features =strsplit(config$Cell_type_markers, ",")[[1]], size = 8,  assay="integrated") + NoLegend()
print(Heatmap2)
dev.off()

png(file =paste0(config$integration_folder,"Cell_type_markers_Violin.png"),  
width = 1600,
height = 1600)
message("**Saving Violin plot with marker expression**")

cell_markers_vln <- VlnPlot(Integrated_data, features =strsplit(config$Cell_type_markers, ",")[[1]], pt.size = 0, assay="integrated")
print(cell_markers_vln)
dev.off()

message("**Starting DE calculation...**")
metadata <- read.csv(config$metadata_file, header = T)

if(length(unique(Integrated_data@meta.data$orig.ident)) >1) {
  #Get DE for each cluster between conditions. Here we can select which conditions to compare.
  DE0 <- function(cluster0, Reference, KO){
    cell1 = rownames(Integrated_data@meta.data)[Integrated_data@meta.data$seurat_clusters == cluster0 & Integrated_data@meta.data$Condition %in% Reference]
    cell2 = rownames(Integrated_data@meta.data)[Integrated_data@meta.data$seurat_clusters == cluster0 & Integrated_data@meta.data$Condition %in% KO]
    DE0 = FindMarkers(Integrated_data, ident.1 = cell1, ident.2 = cell2,logfc.threshold=config$log.fc_threshold_DE, test.use = "MAST",latent.vars="orig.ident")
    DE0$cluster = cluster0
    DE0$Symbol <- rownames(DE0)
    return(DE0)
  }
  
  
  cluster_list <- as.numeric(levels(Integrated_data@meta.data$seurat_clusters))
  #Select which comparison you want to do.
  conditions <- levels(Integrated_data@meta.data$Condition)
  
  #Here a create a df with all possible comparison depending of how many conditions we have.
  possible_conditions <- t(combn(conditions,2))
  
  #I run the DE analysis plus I create the files for enrichment analysis.
  for (row in 1:nrow(possible_conditions)) {
    print(paste0("Calculating DEG for ", possible_conditions[row, ][1], "_vs_", possible_conditions[row, ][2]))
    DE_list1 <- lapply(FUN=DE0,cluster_list, possible_conditions[row, ][1], possible_conditions[row, ][2])
    DE_results1 = do.call(rbind,DE_list1)
    write.table(DE_results1, 
                file = paste0(config$DE_folder,possible_conditions[row, ][1], "_vs_", possible_conditions[row, ][2], "_DE_results.tsv"),
                sep = ",", 
                col.names = T, 
                row.names = F, 
                quote = F)
    message(paste0("**Saving DE results in ", config$DE_folder),"**")
    
    message(paste0("**Saving files for enrichments analysis in  ", config$enrichment_folder),"**")
    
    
    #Here I start to create the files for the enrichment analysis.
    if(!file.exists(paste0(config$enrichment_folder,possible_conditions[row, ][1],"_",possible_conditions[row, ][2]))){
      dir.create(paste0(config$enrichment_folder,possible_conditions[row, ][1],"_",possible_conditions[row, ][2])) 
    } 
    
    print(paste0("Creating directory for ",possible_conditions[row, ][1],"_",possible_conditions[row, ][2]))
      cluster_names <- unique(DE_results1$cluster)
    DE_results1$Symbol <- rownames(DE_results1)
    for (i in cluster_names){
      DE_cluster <- DE_results1 %>% filter(cluster==i)
      DE_cluster$value <- as.numeric(DE_cluster$avg_log2FC > 0)
      DE_cluster2 <- DE_cluster %>% 
        dplyr::arrange(p_val) %>%
        dplyr::mutate(log_pvalue = ifelse(value == 0, (-log10(DE_cluster$p_val)*-1), (-log10(DE_cluster$p_val)*1))) %>% dplyr::select(Symbol, log_pvalue)
      write.table(DE_cluster2, 
                  paste0(config$enrichment_folder,possible_conditions[row, ][1],"_",possible_conditions[row, ][2],"/Cluster_",i, "_rank.rnk"), 
                  sep = "\t", 
                  row.names = F, 
                  col.names = F, 
                  quote = F)
                  }
    }
}

print(Sys.time())
#message("Creating html report")                                    
#knitr::spin("scRNAseq_downstream_analysis.R") 
