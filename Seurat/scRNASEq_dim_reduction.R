#!/usr/bin/env Rscript
message("Starting job...")

#.libPaths("/share/home/maider/R_libraries/R/x86_64-redhat-linux-gnu-library/4.1")
#Uncomment this if using in the server

message("Loading required libraries")

#Libraries needed for the analysis
require(dplyr)
require(Seurat)
require(patchwork)
require(stringr)
require(plyr)
require(data.table)
require(ggplot2)
require(Matrix)
require(config)
#require(knitr)
#require(markdown)

print(Sys.time())

rm(list = ls())
#Read the config file
config <- config::get(file ="config_seurat.yml")

message("**Loading integrated dataset**")

#Load the craeted Seurat object.
load(paste0(config$seurat_object,config$seurat_object_name))

message("**Starting with dimension reduction...**")

#Dimensionality reduction
Integrated_data <- RunUMAP(Integrated_data, reduction = "pca", dims = 1:config$npcs)
Integrated_data <- FindNeighbors(Integrated_data, reduction = "pca", dims = 1:config$npcs)
Integrated_data <- FindClusters(Integrated_data, resolution = config$resolution)

message(paste0("**Saving data in ", config$seurat_object),"**")

#Plot integration UMAP and clustering UMAP.
message("**Saving plot of UMAP based on different conditions**")

png(file=paste0(config$integration_folder,"UMAP_integration_samples.png"),  
width = 1200,
height = 1200)
plot1 <- DimPlot(Integrated_data, reduction = "umap", group.by = "Condition")
print(plot1)
dev.off()

message("**Saving plot of UMAP clustering**")

png(file =paste0(config$integration_folder,"UMAP_clustering.png"),  
width = 1200,
height = 1200)
p2 <- DimPlot(Integrated_data, reduction = "umap", label = TRUE, repel = TRUE)

print(p2)
dev.off()

#Save the integrated object so we can play with it.
save(Integrated_data, file=paste0(config$seurat_object,config$seurat_object_dim))

print(Sys.time())
#message("Creating html report")                                    
#knitr::spin("scRNASEq_dim_reduction.R")   
