#!/usr/bin/env Rscript
message("Starting job...")

#.libPaths("/share/home/maider/R_libraries/R/x86_64-redhat-linux-gnu-library/4.1")

message("Loading required libraries")
#Libraries needed for the analysis
require(dplyr)
require(base)
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

message("**Reading samples...**")
#Create a function to run through each sample.
# Load dataset with count matrix and gene and barcode matrices.
#These directory have to have inside folder for each samples and inside the three necessary files.
cell_ranger_output <- list.dirs(config$data_location,recursive = F)
#print(cell_ranger_output)
#Now I get the names to call the objects and that have to be the same as the folder.
sample_names <- sub('.*\\/', '', cell_ranger_output)
print(sample_names)

#FUNCTION_1: Create demultiplexed objects from each multiplexed files. This is coming from Cellranger COUNT mode
if (config$cellranger_option=="COUNT"){
create_object <- function(sample_list,data_dir){

    Sample_dataset <- Read10X(data_dir)
    Sample_expression <- CreateSeuratObject(counts = Sample_dataset[["Gene Expression"]], project = "Sample_test")
    Sample_antibody <- CreateSeuratObject(counts = Sample_dataset[["Antibody Capture"]], project = "Sample_test")

    # filtered the cells for you, but perform this step for clarity.
    joint.bcs <- intersect(colnames(Sample_expression), colnames(Sample_antibody))

    # Subset RNA and HTO counts by joint cell barcodes
    Sample_expression <- Sample_expression[, joint.bcs]
    Sample_antibody <- as.matrix(Sample_antibody@assays$RNA@counts[, joint.bcs])

    # Normalize RNA data with log normalization
    Sample.hashtag <- NormalizeData(Sample_expression)
    # Find and scale variable features
    Sample.hashtag <- FindVariableFeatures(Sample.hashtag, selection.method = "mean.var.plot")
    Sample.hashtag <- ScaleData(Sample.hashtag, features = VariableFeatures(Sample.hashtag))

    # Add HTO data as a new assay independent from RNA
    Sample.hashtag[["HTO"]] <- CreateAssayObject(counts = Sample_antibody, project = "Sample_test")


    #Demultiplex cells based on HTO enrichment
    # Normalize HTO data, here we use centered log-ratio (CLR) transformation
    Sample.hashtag <- NormalizeData(Sample.hashtag, assay = "HTO", normalization.method = "CLR")

    # If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid
    # clustering function for large applications You can also play with additional parameters (see
    # documentation for HTODemux()) to adjust the threshold for classification Here we are using the
    # default settings
    Sample.hashtag <- HTODemux(Sample.hashtag, assay = "HTO", positive.quantile = 0.99)
    #The sample assigment is saved in the metadata dataframe inside the object Sample.hahstag.



    #Visualize demultiplexing results
    # Global classification results
    #table(Sample.hashtag$HTO_classification.global)

    # Group cells based on the max HTO signal
    Idents(Sample.hashtag) <- "HTO_maxID"
    RidgePlot(Sample.hashtag, assay = "HTO", features = rownames(Sample.hashtag[["HTO"]])[1:2], ncol = 2)


    Idents(Sample.hashtag) <- "HTO_classification.global"

    png(file =paste0(config$quality_location,sample_list,"nCount_violin.png"),  
      width = 1000,
      height = 1200)

    nCount_violin <- VlnPlot(Sample.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
    print(nCount_violin)
    dev.off()


    # First, we will remove negative cells from the object
    Sample.hashtag.subset <- subset(Sample.hashtag, idents = "Negative", invert = TRUE)

    # Calculate a distance matrix using HTO
    hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = Sample.hashtag.subset, assay = "HTO"))))

    # Calculate tSNE embeddings with a distance matrix
    Sample.hashtag.subset <- RunTSNE(Sample.hashtag.subset, distance.matrix = hto.dist.mtx, perplexity = 100)
    #DimPlot(Sample.hashtag.subset)

    # Extract the singlets  
    Sample.singlet <- subset(Sample.hashtag, idents = "Singlet")
    #Normalize data and find variable genes.
    Samples <- SplitObject(Sample.singlet, split.by = "HTO_maxID")
}

Object_list <- mapply(FUN=create_object,sample_names, cell_ranger_output)
Object_list <- unlist(Object_list)
names(Object_list) <- gsub(".*\\.", "",names(Object_list))
samples_names <- names(Object_list)

  } else if (config$cellranger_option=="MULTI") {
    
  create_object <- function(sample_list,data_dir){
    Sample_dataset <- Read10X(data_dir)
    # Initialize the Sample_object object with the raw (non-normalized data).
    Sample_object <- CreateSeuratObject(counts = Sample_dataset[["Gene Expression"]], project = paste0(sample_list), min.cells = config$min_cells, min.features = config$min_features) 
  } 
  
Object_list <- mapply(FUN=create_object,sample_names, cell_ranger_output)

} else {
    print("cellranger option is not correct")
  }

#After creating obejcts for each sample start the filtering steps:
filtering_object <- function(Sample_object,sample_list){  
  Sample_object[["percent.mt"]] <- PercentageFeatureSet(Sample_object, pattern = "^MT-")
  #Filtering
  Sample_object <- subset(Sample_object, subset = nFeature_RNA > config$nFeature_min & nFeature_RNA < config$nFeature_max & nCount_RNA < config$nCount_max & nCount_RNA > config$nCount_min & percent.mt < config$percent_mit)
  
  png(file =paste0(config$quality_location,sample_list,"_Violin_data_quality.png"),  
      width = 1000,
      height = 1200)
  
  quality_plot <- VlnPlot(Sample_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
  print(quality_plot)
  dev.off()

  #Check features correlations.
  png(file =paste0(config$quality_location,sample_list,"_feature_correlations.png"),  
      width = 1000,
      height = 1200)

  plot1 <- FeatureScatter(Sample_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(Sample_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- plot1 + plot2

  print(plot3)
  dev.off()

  
  #Normalize data and find variable genes.
  Sample_object <- NormalizeData(Sample_object)
  Sample_object <- FindVariableFeatures(Sample_object, selection.method = "vst", nfeatures = config$HVF_selection)
}

Object_list <- lapply(FUN=filtering_object, Object_list, sample_names)
message("**Starting with QC and data normalization...**")
#Run the function through all samples and get one object per sample. 



print(Sys.time())

#message("**Subseting data**")
#In order to make a test lets subset all objects. This one is optional.
#subset_object <- function(Seurat_object){
# Seurat_object <- subset(x = Seurat_object, downsample = 500)
#}

#Object_list <- mapply(FUN=subset_object,Object_list)
features <- SelectIntegrationFeatures(object.list = Object_list, nfeatures = config$HVF_selection)

dim_red <- function(Object_list){
  Object_list <- ScaleData(Object_list, features = features, verbose = FALSE)
  Object_list <- RunPCA(Object_list, features=features, verbose = FALSE)
}

Object_list <- lapply(FUN=dim_red, Object_list)
print(Sys.time())

options(future.global.maxSize= 8000 * 1024^2)

message("**Starting integration...**")
if (length(Object_list) > 1) {
  #In order to integrate data select integration features through the objects. 
  
  anchors <- FindIntegrationAnchors(object.list = Object_list, anchor.features = features, reduction = "rpca")
  
  #Create 'integrated' data assay through all the samples. 
  Integrated_data <- IntegrateData(anchorset =anchors)
  
  DefaultAssay(Integrated_data) <- "integrated"
} else {
  Integrated_data <- Object_list[[1]]
}
print(Sys.time())

message("**Load metadata and add information to object**")
#Add to metadata condition information about each sample.
metadata <- read.csv(config$metadata_file, header = T)
cell_names <- rownames(Integrated_data@meta.data)
Integrated_data@meta.data <- base::merge(Integrated_data@meta.data, metadata, by="orig.ident", all=T)
rownames(Integrated_data@meta.data) <- cell_names
Integrated_data@meta.data$Condition <- relevel(factor(Integrated_data@meta.data$Condition), ref=config$condition_reference)



#Dimensionality determination
#Scale data for visualization.
Integrated_data <- ScaleData(Integrated_data, verbose = FALSE)

Integrated_data <- RunPCA(Integrated_data, npcs = config$npcs, verbose = FALSE)
png(file =paste0(config$quality_location,"elbow_plot.png"),  
      width = 1000,
      height = 1200)
elbowplot <- ElbowPlot(Integrated_data)
print(elbowplot)
dev.off()


message(paste0("**Saving integrated data in ", config$seurat_object),"**")
save(Integrated_data, file=paste0(config$seurat_object,config$seurat_object_name))



#Save integrated raw count matrix
message("**Saving integrated raw count matrix**")
Raw_counts_mat <- as.matrix(Integrated_data[["RNA"]]@counts)

write.table(Raw_counts_mat, file=paste0(config$seurat_object,"Raw_count_mat.csv"), quote=F, row.names=T, col.names=T, sep=",")


#Save integrated normalized count matrix
message("**Saving integrated normalized count matrix**")
Normalized_counts_mat <- as.matrix(Integrated_data[["RNA"]]@data)

write.table(Normalized_counts_mat, file=paste0(config$seurat_object,"Normalized_count_mat.csv"), quote=F, row.names=T, col.names=T, sep=",")




print(Sys.time())
#message("Creating html report") 
#knitr::spin("scRNASeq_Integration.R")
