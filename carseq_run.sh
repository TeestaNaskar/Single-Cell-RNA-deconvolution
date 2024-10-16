#BSUB -L /bin/bash
#BSUB -P acc_DADisorders
#BSUB -q premium
#BSUB -n 1
#BSUB -R span[ptile=14]
#BSUB -R rusage[mem=320000]
#BSUB -W 48:00
#BSUB -o stdout.%J.%I
#BSUB -e stderr.%J.%I

ml R/4.1.0
module load openssl
module load udunits


R CMD BATCH washU.carseq.R
