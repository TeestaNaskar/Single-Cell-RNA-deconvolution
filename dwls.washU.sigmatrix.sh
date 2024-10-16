#BSUB -L /bin/bash
#BSUB -P acc_DADisorders
#BSUB -q premium
#BSUB -n 1
#BSUB -R span[ptile=4]
#BSUB -R himem
#BSUB -R rusage[mem=14000000]
#BSUB -W 144:00
#BSUB -o dwlsstdout.%J.%I
#BSUB -e dwlsstderr.%J.%I

module load R
module load openssl
module load udunits

R CMD BATCH washU.dwls.R
