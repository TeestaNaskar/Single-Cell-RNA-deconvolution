#seuratObj = readRDS("sc.NormByLocationRep.Harmony.final.rds")
#this script is for doing single cell RNA deconvolution of human placenta with the template data from WashU paper
setwd("/sc/arion/projects/MetaDope/Teesta/WashU.singlecelldata")

library(scCustomize)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(readxl)
library(openxlsx)
library(tibble)
library(magrittr)
library(dplyr)
seuratObj = readRDS("sc.NormByLocationRep.Harmony.final.rds")
seuratObj = FindClusters(seuratObj, resolution = 1)

#next step is to export the count data from the seuratObject
#and metadata

metadata = seuratObj@meta.data
count = seuratObj@assays$RNA@counts
print(dim(count))
#since its taking too much time to finish so subsetting the object as per random sampling
Idents(seuratObj) = "assig2"
object_subset <- subset(seuratObj, idents = c("HPL20874-M", "HPL20877-F", "HPL20878-M")) #specifying selected subjects, we can add more or all depending upon the question
count = object_subset@assays$RNA$counts
print(dim(count))
data = count %>% as.matrix %>% t %>% as.data.frame
print('Loaded cells')
print(Sys.time())
print(data[1:10, 1:10])
print(head(colnames(data)))
print(head(rownames(data)))

metadata = object_subset@meta.data
#remove low-expressed genes
gene_names <- colnames(data)[2:ncol(data)]
sample_names <- data[,1]
print(head(gene_names))
print(head(sample_names))
data = data[, 2:ncol(data)] #remove sample column
colsums = colSums(data)
print('Removing low-expressed genes')
print(Sys.time())
dim(data)
data = data[, colsums > 10] #remove low-expressed genes
print('Finished removing non-expressed genes')
print(dim(data))
print(Sys.time())
print('Transposing cells')
print(Sys.time())
data <- data.table::transpose(data)
print('Transposed cells')
print(Sys.time())

#set column and row names to sample names and gene names appropriately
rownames(data) <- gene_names[which(colsums > 10)]
#Extract the celltype labels
labels <- metadata$cluster_name
#load required package for doing dwls
library(DWLS)
library(stringr)
library(data.table)
library(openxlsx)
require(Matrix)

#Build signature matrix: takes several hours and requires large memory
#allocated 320GB for this dataset (70K cells). DWLS doesn't support multi-threading
#so this was run on one core
print('Building signature')
Signature<-buildSignatureMatrixMAST(scdata=data,
                                   id=labels,
                                  path="../RNAdeconvol/script/DWLS.results",
                                 diff.cutoff=0.5,
                                   pval.cutoff=0.01)
save(Signature,file="washU_DWLS.RData")
load("washU_DWLS.RData")
