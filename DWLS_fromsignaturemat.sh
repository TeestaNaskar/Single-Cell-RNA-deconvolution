#BSUB -L /bin/bash
#BSUB -P acc_DADisorders
#BSUB -q premium
#BSUB -n 1
#BSUB -R span[ptile=4]
#BSUB -R himem
#BSUB -R rusage[mem=550000]
#BSUB -W 144:00
#BSUB -o dwlsstdout.%J.%I
#BSUB -e dwlsstderr.%J.%I

module load  R/4.0.3
module load openssl
module load udunits

R CMD BATCH DWLS_fromsignaturematrix.R

##note: this shell script is for runnning the R script DWLS_fromsignaturematrix.R
