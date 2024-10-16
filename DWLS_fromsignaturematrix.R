#this script is for estimating the cell fraction after having the signature matrix
setwd("/sc/arion/projects/MetaDope/Teesta/WashU.singlecelldata")
library(data.table, lib='/sc/arion/projects/MetaDope/Randy/R_4.0.3/')
library(DWLS, lib='/sc/arion/projects/MetaDope/Randy/R_4.0.3/')
library(stringr)
library(Seurat)
library(dplyr)
library(openxlsx)
require(Matrix)
library(biomaRt)

load("washU_DWLS.RData")

#import bulk data
#counts = read.csv("../../Proj_Placenta/RNAseq/ALL_COUNTS.csv", check.names = F)
counts <- read.csv("/sc/arion/projects/MetaDope/Teesta/Proj_Placenta/RNAseq/humanplacenta.withEnsembleID.VST_counts.csv", check.names = F)
rownames(counts) <- make.names(counts[,1], unique = TRUE) #gene symbols
counts <- counts[,2:ncol(counts)]
print(head(colnames(counts)))
##loading meta data and consider counts only for 93 subjects that are without SANDY
meta = read.xlsx("/sc/arion/projects/MetaDope/Teesta/Proj_Placenta/RNAseq/HUMAN.PLACENTA.METADATA/WorkingMetadata.updatedbyTeesta.consolidatedinfo.Anissa.Greg.Placenta.inventory.xlsx", sheet = 3)

meta = meta[1:93,]
print(dim(meta))
counts = counts[,colnames(counts) %in% meta$Placenta_Seq_ID]
print(dim(counts))

samples<- colnames(counts) # grabbing a vector of your subject names
print(head(rownames(counts)))
print(head(colnames(counts)))

#converting to CPM 
#CPM <-sweep(counts,2,as.numeric(colSums(counts)),FUN="/")*1000000 
#print(head(rownames(CPM)))
#print(head(colnames(CPM)))
#print(head(rownames(sig)))
#rownames of sig contain .1, .2 like that so need to remove after . using gsub
rownames(Signature) = sub("\\..*$", "", rownames(Signature))

print(head(rownames(Signature)))
#converting rownames of sig from ensemble id to gene names by using biomart
 
#value = rownames(sig)
#ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

#genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "chromosome_name"), values = value,mart = ensembl)

#idx <- match(value, genemap$ensembl_gene_id)
#sig$gene_name <- genemap$hgnc_symbol[idx]
#rownames(sig) = make.names(sig$gene_name, unique = TRUE)
#print(sig[1:5, 1:5])



#CPM.trim <- CPM[rownames(sig),]
#samples<- colnames(CPM.trim)
#print(head(samples))
#print(dim(CPM.trim))

# Initializing Results Data Frames for Each Algorithm that DWLS employs (DWLS being the superstar)
DWLS <- data.frame()
SVR <- data.frame()
OLS <- data.frame()

# Computing Cell Fractions for Each Sample and Adding to Results Data Frames
# Computing Cell Fractions for Each Sample and Adding to Results Data Frames
for(sample in samples){
  bulk <- counts[,sample]
  names(bulk) <- rownames(counts) #subsetting column deletes rownames; add gene symbols back!
  tr<-trimData(Signature,bulk)
  solDWLS<-solveDampenedWLS(tr$sig,tr$bulk)
  solSVR <- solveSVR(tr$sig,tr$bulk)
  solOLS <- solveOLS(tr$sig,tr$bulk)
  DWLS <-rbind(DWLS,solDWLS)
  SVR <- rbind(SVR,solSVR)
  OLS <- rbind(OLS,solOLS)
}

model<- list(DWLS,SVR,OLS) # List of dataframes

# Adding the rownames and columns back to each dataframe in list
L <- lapply(model, function(df)
            {
              colnames(df) <- colnames(Signature)
              rownames(df) <- samples
              df
            }
           )

# Writing each dataframe in list to file
write.csv(L[1],file="placentahuman_CellFractions_tees_ref_DWLS.csv")
write.csv(L[2],file="placentahuman_CellFractions_tees_pDMS_ref_SVR.csv")
write.csv(L[3],file="placentahuman_CellFractions_tees_pDMS_ref_OLS.csv")
