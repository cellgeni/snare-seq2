#!/bin/bash -e

# example parameters
# TAG=test
# path to fq folders
# ATACDIR=/lustre/scratch126/cellgen/cellgeni/pasham/raw/2110.SNARE-seq/v04.fq/SNAREsq313165052
# CDNADIR=/lustre/scratch126/cellgen/cellgeni/pasham/raw/2110.SNARE-seq/v04.fq/SNAREsq313165053

# string to filter (grep on full path) fqs for specific sample
# ATACTAG='SNAREsq313165052'
# CDNATAG='SNAREsq313165053'
# expects following files for each run
# ATAC: R1 and R3 - bio; R2 - barcodes;  
# cDNA: R1 - bio; R2 - barcodes; R3 UMI (pos 0-10)
# barcodes shouls be at 15-23,53-61,91-99
# SPECIES='human'


# parameters
TAG=$1
ATACDIR=$2
CDNADIR=$3
ATACTAG=$4
CDNATAG=$5
SPECIES=$6

# paths to bins and scripts
STAR=/nfs/users/nfs_p/pm19/nfs/ghub/STAR/source/STAR
BWA="/nfs/cellgeni/pasham/bin/bwa-mem2/bwa-mem2 mem"
PIPELINEHOME=`echo "$(cd "$(dirname "$0")"; pwd)"`
#/nfs/cellgeni/pasham/projects/2110.SNARE-seq/src/pipeline/src

CPUS=16
cd $TAG

# reference
if [ $SPECIES == 'human' ]
then
  STARREF=/nfs/cellgeni/STAR/human/2020A/index
  BWAREF=/nfs/cellgeni/pasham/bin/2020A.bwa.index/GRCh38_v32_modified
elif [ $SPECIES == 'mouse' ]
then
  STARREF=/nfs/cellgeni/STAR/mouse/2020A/index
  BWAREF=/nfs/cellgeni/pasham/bin/2020A.bwa.index/mm10_vM23_modified
else
  echo "only human or mouse are allowed"
  exit
fi

STARBC=/nfs/cellgeni/pasham/projects/2110.SNARE-seq/src/v04/barcodes/bc.white.list.txt
BARCODES=$PIPELINEHOME/barcodes/barcode



# init
source activate ngs
extractbc="python $PIPELINEHOME/extractBCnUMI.py"
bindfa="python $PIPELINEHOME/bind.fastqs.py"
addcom="python $PIPELINEHOME/fastqs2comment.py"

mkdir $TAG
cd $TAG

# preprocess barcodes and UMIs
# _extract barcodes from R2
# _ATAC
for i in `ls -1 $ATACDIR/*R2_001.fastq.gz | grep $ATACTAG`
do
  f=`basename $i | sed s/_R2_001.fastq.gz//`
  echo $f
  $extractbc $ATACDIR/${f}_R2_001.fastq.gz 15-23,53-61,91-99 ${BARCODES}1.txt,${BARCODES}2.txt,${BARCODES}3.txt 2> atac_${f}_bc_001.json | gzip > atac_${f}_bc_001.fastq.gz
done

# cDNA
for i in `ls -1 $CDNADIR/*R2_001.fastq.gz | grep $CDNATAG`
do
  f=`basename $i | sed s/_R2_001.fastq.gz//`
  echo $f
  $extractbc $CDNADIR/${f}_R2_001.fastq.gz 15-23,53-61,91-99 ${BARCODES}1.txt,${BARCODES}2.txt,${BARCODES}3.txt 2> cdna_${f}_bc_001.json | gzip > cdna_${f}_bc_001.fastq.gz
done

# cDNA: extract umis and join with bcs
for i in `ls -1 $CDNADIR/*R3_001.fastq.gz | grep $CDNATAG`
do
  f=`basename $i | sed s/_R3_001.fastq.gz//`
  echo $f
  $extractbc  $CDNADIR/${f}_R3_001.fastq.gz 0-10 '' 2> cdna_${f}_umi_001.json | gzip > cdna_${f}_umi_001.fastq.gz
  $bindfa cdna_${f}_bc_001.fastq.gz cdna_${f}_umi_001.fastq.gz | gzip > cdna_${f}_bcumi_001.fastq.gz
done

# map RNA
TAG=cDNA
UMILEN=10   ## 10x v2
STR=Forward ## 3' 10x
BAM="--outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 2 --limitBAMsortRAM 60000000000 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM CB UB GX GN CR CY UR UY"
CBLEN=24 # default 16
UMIstart=`expr $CBLEN + 1`

mkdir $TAG
## for multiple fastq files; change grep options according to your file format
R1=`find ./* | grep cdna | grep "_bcumi_" | sort | tr '\n' ','`
R2=`find $CDNADIR/* | grep $CDNATAG | grep "_R1_" | sort | tr '\n' ','`

echo "cDNA Reads:"
echo $R1
echo $R2

$STAR --runThreadN $CPUS --genomeDir $STARREF --readFilesIn $R2 $R1 --runDirPerm All_RWX --readFilesCommand zcat \
     --soloType CB_UMI_Simple --soloCBwhitelist $STARBC --soloBarcodeReadLength 1 --soloUMIlen $UMILEN --soloStrand $STR \
     --soloUMIdedup 1MM_CR --soloCBmatchWLtype Exact --soloUMIfiltering MultiGeneUMI_CR \
     --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
     --soloCBlen $CBLEN --soloUMIstart $UMIstart\
     $BAM \
     --outFileNamePrefix $TAG/ \
     --soloFeatures Gene GeneFull Velocyto --soloOutFileNames output/ genes.tsv barcodes.tsv matrix.mtx


# add commentary to ATAC fq to preserve BC
for i in `ls -1 $ATACDIR/*R1_001.fastq.gz | grep $ATACTAG | sed s/_R1_001.fastq.gz//`
do
  f=`basename $i`
  echo $f
  $addcom  $ATACDIR/${f}_R1_001.fastq.gz atac_${f}_bc_001.fastq.gz 'CB:Z:' | gzip > atac_${f}_R1c_001.fastq.gz
  $addcom  $ATACDIR/${f}_R3_001.fastq.gz atac_${f}_bc_001.fastq.gz 'CB:Z:' | gzip > atac_${f}_R3c_001.fastq.gz
done

# map ATAC
mkdir ATAC

R1=`find ./* | grep atac | grep "_R1c_" | sort | tr '\n' ' '`
R2=`find ./* | grep atac | grep "_R3c_" | sort | tr '\n' ' '`

R1="'""< zcat $R1""'"
R2="'""< zcat $R2""'"
echo "ATAC Reads:"
echo $R1
echo $R2


cmd="$BWA -t $CPUS -C $BWAREF $R1 $R2 2> ATAC/bwa.log | samtools view -b > ATAC/out.bam"
eval $cmd

samtools sort --threads $CPUS -m 3G -o ATAC/out.sorted.bam ATAC/out.bam
samtools index ATAC/out.sorted.bam

# run archR
$PIPELINEHOME/run.archr.R . $SPECIES $CPUS $PIPELINEHOME/barcodes/bc.white.list.txt

cd ..