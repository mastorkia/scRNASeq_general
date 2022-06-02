# scRNASeq
Main pipeline to analyze 10xGenomics Chromium Single Cell Gene Expression and TotalSeq single cell RNASeq data

### Cellranger

#### Running Cellranger for TotalSeq data
In TotalSeq sequencing antibody hashing was performed to multiplex the data. We identified two ways of demultiplexing this data:
1. Use Cellranger multi option
2. Use Cellranger count option

In both ways the main script for Cellranger will be TotalSeq_cellranger.py but the output will be different. In MULTI mode the samples will be demultiplexed and in COUNT mode further R script steps are necessary. 

Usage of ```TotalSeq_cellranger.py``` is as follows :

```python
 -c/--config-file                 - The configuration file (workflow order and software parameters)
 -f/--fastq-file-list             - File list containing complete path names for fastq.gz files.
 -a/--aim                         - Cellranger mode to run the analysis. Could be MULTI/COUNT
```
Since we are getting the whole path in -f option our analysis will be done in that directory. So we can run our command in our own home but the analysis will be performed in the path directory.

In the path direcotry as long with fastq files we will need a CMO_reference.csv file (similar to Illumina Samplesheet) and config files for each of the samples. These files will change a little bit depending on MULTI or COUNT purposes.

Information about HTO - sample link and HTO sequence have to be provided by the company.

#### Cellranger MULTI option

In the TotalSeq/MULTI directory you will have the main script ```TotalSeq_cellranger.py```:

Example:
```
python TotalSeq_cellranger.py -f sample.lst -c config_cellranger.yml -a MULTI
```

where, ```sample.lst``` file contains the absolute path to fastq file(s). Each sample will be comprised of 8 fastq files. 
```bash
$ cat ./cellranger/MULTI/sample.list 
/share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022/test_cellranger/MMC-79-89-93-A_Ab_S1_L001_I1_001.fastq.gz
/share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022/test_cellranger/MMC-79-89-93-A_Ab_S1_L001_I2_001.fastq.gz
/share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022/test_cellranger/MMC-79-89-93-A_Ab_S1_L001_R1_001.fastq.gz
/share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022/test_cellranger/MMC-79-89-93-A_Ab_S1_L001_R2_001.fastq.gz
/share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022/test_cellranger/MMC-79-89-93-A_GEX_S1_L002_I1_001.fastq.gz
/share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022/test_cellranger/MMC-79-89-93-A_GEX_S1_L002_I2_001.fastq.gz
/share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022/test_cellranger/MMC-79-89-93-A_GEX_S1_L002_R1_001.fastq.gz
/share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022/test_cellranger/MMC-79-89-93-A_GEX_S1_L002_R2_001.fastq.gz
```

And where, ```CMO_reference.csv``` and ```MMC-79-89-93-A.config``` will look in the ```/share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022/test_cellranger/``` location:


```bash
$ cat ./cellranger/MULTI/CMO_reference.csv 
id,name,read,pattern,sequence,feature_type
HTO1,HTO1,R2,5PNNNNNNNNNN(BC),GTCAACTCTTTAGCG,Multiplexing Capture
HTO2,HTO2,R2,5PNNNNNNNNNN(BC),TGATGGCCTATTGGG,Multiplexing Capture
HTO3,HTO3,R2,5PNNNNNNNNNN(BC),TTCCGCCTCTCTTTG,Multiplexing Capture
HTO4,HTO4,R2,5PNNNNNNNNNN(BC),AGTAAGTTCAGCGTA,Multiplexing Capture

$ cat MULTI/MMC-79-89-93-A.config
[gene-expression]
reference,/share/data/RNA_Seq/reference_genomes/10X_genomics/human/refdata-cellranger-GRCh38-3.0.0/
cmo-set,/share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022/CMO_reference.csv

[libraries]
fastq_id,fastqs,feature_types
MMC-79-89-93-A_GEX,/share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022,Gene Expression
MMC-79-89-93-A_Ab,/share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022,Multiplexing Capture #created in step1 

[samples]
sample_id,cmo_ids
MMC-0079,HTO1
MMC-0089,HTO2
MMC-0093,HTO3
```

#### Cellranger COUNT option

In the TotalSeq/COUNT directory you will have the main script ```TotalSeq_cellranger.py```:

Example:
```
python TotalSeq_cellranger.py -f sample.lst -c config_cellranger.yml -a COUNT
```

where, ```sample.lst``` file contains the absolute path to fastq file(s). Each sample will be comprised of 8 fastq files. 
```bash
$ cat ./cellranger/COUNT/sample.list 
/share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022/test_cellranger/MMC-79-89-93-A_Ab_S1_L001_I1_001.fastq.gz
/share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022/test_cellranger/MMC-79-89-93-A_Ab_S1_L001_I2_001.fastq.gz
/share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022/test_cellranger/MMC-79-89-93-A_Ab_S1_L001_R1_001.fastq.gz
/share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022/test_cellranger/MMC-79-89-93-A_Ab_S1_L001_R2_001.fastq.gz
/share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022/test_cellranger/MMC-79-89-93-A_GEX_S1_L002_I1_001.fastq.gz
/share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022/test_cellranger/MMC-79-89-93-A_GEX_S1_L002_I2_001.fastq.gz
/share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022/test_cellranger/MMC-79-89-93-A_GEX_S1_L002_R1_001.fastq.gz
/share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022/test_cellranger/MMC-79-89-93-A_GEX_S1_L002_R2_001.fastq.gz
```
And where, ```CMO_reference.csv``` and ```MMC-79-89-93-A_count_Abs.config``` will look in the ```/share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022/test_cellranger/``` location:


```bash
$ cat ./cellranger/COUNT/CMO_reference.csv 
id,name,read,pattern,sequence,feature_type
HTO1,HTO1,R2,5PNNNNNNNNNN(BC),GTCAACTCTTTAGCG,Antibody Capture
HTO2,HTO2,R2,5PNNNNNNNNNN(BC),TGATGGCCTATTGGG,Antibody Capture
HTO3,HTO3,R2,5PNNNNNNNNNN(BC),TTCCGCCTCTCTTTG,Antibody Capture
HTO4,HTO4,R2,5PNNNNNNNNNN(BC),AGTAAGTTCAGCGTA,Antibody Capture

$ cat ./cellranger/COUNT/MMC-79-89-93-A_count_Abs.config
fastqs,sample,library_type
/share/data/RNA_Seq/MECFS_scRNAseq/Working/test_April25_2022,MMC-79-89-93-A_GEX,Gene Expression
/share/data/RNA_Seq/MECFS_scRNAseq/Working/test_April25_2022,MMC-79-89-93-A_Ab,Antibody Capture 
```
#### Running Cellranger for Chromium Single Cell Gene Expression data

In the 10xCOUNT directory you will have the main script ```cellranger_count.py```

Usage:

```python
 -c/--config-file                 - The configuration file (workflow order and software parameters)
 -f/--fastq-file-list             - File list containing complete path names for fastq.gz files.
```

Example:
```
python cellranger_count.py -f sample.lst -c config_cellranger.yml
```

where, ```sample.lst``` file contains the absolute path to fastq file(s).

```bash
$ cat ./cellranger/10xCOUNT/sample.list 
pbmc_1k_v3_S1_L001_I1_001.fastq.gz
pbmc_1k_v3_S1_L001_R1_001.fastq.gz
pbmc_1k_v3_S1_L001_R2_001.fastq.gz
pbmc_1k_v3_S1_L002_I1_001.fastq.gz
pbmc_1k_v3_S1_L002_R1_001.fastq.gz
pbmc_1k_v3_S1_L002_R2_001.fastq.gz
```


#### Configuration file
The configuration file is a simple text file that plays a central role in the pipeline. This text file is
composed of different parts allowing to build the workflow,
to provide software parameters, to compress or remove
temporary/intermediate data, and to set up the scheduler if needed (see figure below for more details).

This file works for any of the options above and will be located in each directory.

```bash
$ cat ./cellranger/config_cellranger.yml
transcriptome_loc:/share/data/RNA_Seq/reference_genomes/10X_genomics/human/refdata-cellranger-GRCh38-3.0.0/
software_loc:/share/data/software/cellranger/cellranger-6.1.2/bin/cellranger
expect_cells:1000
local_cores:16
local_mem:10G
```

### Seurat

Once count matrices are obtained from Cellranger we will start working with Seurat package. 

#### Create the desired environment previous to Seurat

First we will create the desired folder structure using ```create_proj_str.py```.

Usage:

```python
 -p/--path  <Project path>
 
```
Example:
```
python ./Seurat/create_proj_str.py -p /share/home/maider/pipelines/scRNAseq
```
Then we will create a symbolic link to all the matrix ouputs from Cellranger using ```symbolic_cellranger_out.py```. 

Usage:

```python
 -cp/--cellranger_path <path containing cellranger output for all samples> 
 -m/--metadata csv file containing sample name and condition information:q!

```
Example:
```
python ./Seurat/symbolic_cellranger_out.py -cp /share/data/RNA_Seq/MECFS_scRNAseq/Working/April19_2022/Cellranger -m seurat_sample_metadata.csv
```

where, ```seurat_sample_metadata.csv``` file contains phenotipic information about the samples. 

```bash
$ cat ./Seurat/seurat_sample_metadata.csv 
porig.ident,MMA ID,Record ID,SAM ID,Current Age,Sex,Collection date,Condition
MMC_0079,MMA_801,463,SAM000801,38,Female,4_22_2021,PASC
MMC_0089,MMA_807,254,SAM000807,41,Female,4_22_2021,HC
MMC_0093,MMA_808,43,SAM000808,30,Female,4_22_2021,ME_CFS
MMC_0097,MMA_809,74,SAM000809,25,Female,4_22_2021,ME_CFS
MMC_0107,MMA_810,503,SAM000810,26,Female,4_22_2021,PASC
MMC_0117,MMA_800,436,SAM000800,30,Female,4_22_2021,PASC
MMC_0125,MMA_802,230,SAM000802,40,Female,4_22_2021,ME_CFS
MMC_0135,MMA_811,45,SAM000811,57,Female,4_22_2021,ME_CFS
MMC_0144,MMA_803,454,SAM000803,47,Female,4_22_2021,PASC
```
### Run Seurat through different scripts. 
Seurat will run through three different scripts. The configuration file will determine the input name, ouput name or the type of analysis we used in cellranger. 

The ```config_seurat.yml``` provides the  modification of the following options:
```bash
$ cat ./Seurat/config_seurat.yml
default:
  data_location: "data"
  cellranger_out: "/share/data/RNA_Seq/10X/Raw_fastq/" #output location for all samples
  cellranger_option: "COUNT" #can be MULTI or COUNT
  metadata_file: "seurat_sample_metadata.csv"
  condition_reference: "HC" #condition to use as reference
  quality_location: "Quality_figures/" #output of quality figures
  seurat_object: "Seurat_objects/" #output of Seurat object
  seurat_object_name: "Seurat_object.RData" #name of the output Seurat object
  seurat_object_dim: "dim_Seurat_object.RData" #name of the output Seurat obejct with dimensio reduction
  integration_folder: "Integration_results/" #outout of downstream integration figures and tables
  DE_folder: "DE_Analysis/" #output of DE analysis
  enrichment_folder: "Enrichment_files/" #output of files to load in GSEA
  
#These all are for filtering steps
  min_cells: 3 
  min_features: 200
  nFeature_max: 8000
  nFeature_min: 50
  nCount_max: 60000
  nCount_min: 20
  percent_mit: 10
  HVF_selection: 2000

#Number of dimensions to use for reduction
  npcs: 30
  resolution: 0.3
  Cell_type_markers: D3D,KLRF1,FCGR3A,CD14,CD1C,IL3RA,IGKC,MYL9,HBB,CD79A

  markers_outfile: "Markers.csv" #output name of markers

#Downstream analysis
  min.pct: 0.1
  log.fc_threshold_markers: 0.25
  log.fc_threshold_DE: 0.25

```

Once the configuration file is set we will run ```scRNASeq_Integration_TotalSeq.R``` script in which data will be filtered, normalized and samples will be integrated in a unique matrix.

```
Rscript scRNASeq_Integration_TotalSeq.R
```

Next, we will perform dimensionality reduction and cell type identification by ```scRNASEq_dim_reduction.R```.

```
Rscript scRNASEq_dim_reduction.R
```

To finish, we will perform DEG and create files to use as input in GSEA through ```scRNAseq_downstream_analysis.R```.

```
Rscript scRNAseq_downstream_analysis.R
```




#### Workflows Management

single cell RNASeq pipeline:
* mapp the reads, demultiplex and get the count matrices, 
* data normalization and dimensionality reduction 
* cell type identification
* downstream analysis: 
  * identification of dysregulated genes between conditions
  * creation of files to use in GSEA enrichment analysis
  


#### Platforms, Installation and Customization
The pipeline currently runs on any recent GNU/Linux system
(Debian Lenny and more, Ubuntu 12.04 and more, and
CentOS 6 and more were tested). 

The pipeline was developed to be straightforward to install in several ways : manually (git clone)
A unique file (localConfig.pm) needs to be filled at installation to ensure the integration of the whole software list
(path and version)



# scRNASeq
