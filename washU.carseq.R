#####important note: this script is subject to test and under modification so please do not consider this script to run for your own experiment
#this script is for running Carseq
setwd('/sc/arion/projects/MetaDope/Teesta/WashU.singlecelldata')

list.of.packages = c("dplyr", "readxl", "tidyr", "stringr", "tibble")

new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

library(data.table)
require(Matrix)
library(openxlsx)


counts = read.csv("/sc/arion/projects/MetaDope/Teesta/Proj_Placenta/RNAseq/ALL_COUNTS.csv", check.names=F)
rownames(counts) <- counts[,1]
counts = counts[,2:ncol(counts)]

meta = read.xlsx('/sc/arion/projects/MetaDope/Teesta/Proj_Placenta/RNAseq/HUMAN.PLACENTA.METADATA/WorkingMetadata.updatedbyTeesta.consolidatedinfo.Anissa.Greg.Placenta.inventory.xlsx', sheet= 2)
#meta$condition = paste0(meta$Rat.ID,'_',meta$Group)

#cellfracDWLS = read.csv('CellFractions_Allen_PL_IL_Alexprl_ref_DWLS.csv')
cellfrac = read.csv("CellFractions_Allen_PL_IL_Alexprl_ref_SVR.csv")
cellfrac=cellfrac[order(match(cellfrac$X, colnames(counts))),]
#CellfracSVR = read.csv("CellFractions_Allen_PL_IL_Alexprl_ref_SVR.csv")
#Remove cells with low fractions as recommended by the developers: https://github.com/Sun-lab/CARseq/issues/3
#According to them, keeping in cells with low fractions can lead to convergence issues and many genes
#having NAs in their results

print(paste('Starting number of cell types:', dim(cellfrac)[2]-1))
#print(paste('Starting number of cell types:', dim(cellfracOLS)[2]-1))
#print(paste('Starting number of cell types:', dim(CellfracSVR)[2]-1))
#make sure you remove the ID column as this will make apply treat all values as characters instead of numerics
max_val_all_celltypes = apply(cellfrac[,2:ncol(cellfrac)], 2, max)
#max_val_all_celltypesOLS = apply(cellfracOLS [,2:ncol(cellfracOLS)], 2, max)
#max_val_all_celltypesSVR = apply(CellfracSVR [,2:ncol(CellfracSVR)], 2, max)
#remove cell types with max values below a threshold
threshold = 0.015
print(threshold)
cells_equalTo_orAbove_threshold = cellfrac[,2:ncol(cellfrac)][, which(max_val_all_celltypes >= threshold)]
cellfrac = cbind(cellfrac$X, cells_equalTo_orAbove_threshold)
colnames(cellfrac)[1] = 'X'
print(paste('After cell type removal:', dim(cellfrac)[2]-1))

print(all(colnames(counts)==cellfrac$X))

#cells_equalTo_orAbove_thresholdOLS = cellfracOLS [,2:ncol(cellfracOLS)][, which(max_val_all_celltypesOLS >= threshold)]
#cellfracOLS = cbind(cellfracOLS$X, cells_equalTo_orAbove_thresholdOLS)
#colnames(cellfracOLS)[1] = 'X'

#cells_equalTo_orAbove_thresholdSVR = CellfracSVR [,2:ncol(CellfracSVR)][, which(max_val_all_celltypesSVR >= threshold)]
#CellfracSVR = cbind(CellfracSVR$X, cells_equalTo_orAbove_thresholdSVR)
#colnames(CellfracSVR)[1] = 'X'



# cellfrac = cellfrac[ , -which(names(cellfrac) %in% c("L2_IT_RHP", "L5_IT_CTX", "L2_3_IT_PPP",
#                                     "L5_NP_CTX","Lamp5","L4_5_IT_CTX", "Sncg","L6b_CTX","Sst","Vip"
#                                     ,"Sst_Chodl"
#                                     ,"L2_3_IT_ENTl"
#                                     ,"L5_IT_TPE_ENT"
#                                     #,"L6_CT_CTX"
#                                     #,"L5_PT_CTX"
#                                     ))]
#print(paste('After cell type removal:', dim(cellfrac)[2]-1))

#print(all(colnames(counts)==cellfrac$X))
meta_sub = meta[meta$Group %in% c("Control", "Cannabis"),]

counts_sub = counts[,colnames(counts) %in% as.character(meta_sub$Rat.ID)]
# counts_sub = counts_sub[c("Vehicle", "10 mg/kg CBD "),]
cellfrac_sub = cellfrac[cellfrac$X %in% meta_sub$Placenta_Seq_ID,]
print(all(colnames(counts_sub)==cellfrac_sub$X))

#cellfracOLS_sub = cellfracOLS[cellfracOLS$X %in% meta_sub$Rat.ID,]
#print(all(colnames(counts_sub)==cellfracOLS_sub$X))

#cellfracSVR_sub = CellfracSVR[CellfracSVR$X %in% meta_sub$Rat.ID,]
#print(all(colnames(counts_sub)==cellfracSVR_sub$X))

print('Confirmed samples are in same order in counts and cell fractions.')


groups=c()
for (sample in colnames(counts_sub)){
  groups = c(groups, meta_sub[meta_sub$Placenta_Seq_ID==sample, 'Group'])
}
print('Groups made.')


library(CARseq)

#groups = c("Vehicle", "`10 mg/kg CBD`")
#resDWLS = run_CARseq(counts, cellfracDWLS_sub[,2:ncol(cellfracDWLS_sub)], groups)
#save(resDWLS, file = "DWLScarseq_results_rawcounts.RData")

#resOLS = run_CARseq(counts, cellfracOLS_sub[,2:ncol(cellfracOLS_sub)], groups)
#save(resOLS, file = "OLScarseq_results_rawcounts.RData")

#resSVR = run_CARseq(counts, cellfracSVR_sub[,2:ncol(cellfracSVR_sub)], groups)
#save(resSVR, file = "SVRcarseq_results_rawcounts.RData")

#res = run_CARseq(counts_sub, cellfrac_sub[,2:ncol(cellfrac_sub)], groups)
  #print(res$p)
  #save(res, file = paste0("new.deconvolution.Alex/carseq_Alex.heroin.prl_",Vehicle,'_',10 mg/kg CBD ,"_OLS.RData"))

res = run_CARseq(counts_sub, cellfrac_sub[,2:ncol(cellfrac_sub)], groups)
  #res.p = data.frame(print(res$p))
  save(res, file = paste0("carseq_rawcounts_SVR.RData"))
