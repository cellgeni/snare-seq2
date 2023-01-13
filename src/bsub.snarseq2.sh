#! /bin/bash -e
#BSUB -G cellgeni
#BSUB -J snar.seq[2-8]
#BSUB -o %J.%I.snar.seq.out
#BSUB -e %J.%I.snar.seq.err
#BSUB -q normal
#BSUB -n 16
#BSUB -M64000
#BSUB -R "span[hosts=1] select[mem>64000] rusage[mem=64000]"
#BSUB -q normal

cd /lustre/scratch126/cellgen/cellgeni/pasham/raw/2110.SNARE-seq/v04.fq/pipeline.out

PARS=`head -n $LSB_JOBINDEX samples.txt | tail -n 1`

/nfs/cellgeni/pasham/projects/2110.SNARE-seq/src/pipeline/snareseq2.sh $PARS
