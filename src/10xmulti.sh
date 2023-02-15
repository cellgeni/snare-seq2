#!/bin/bash -e

# some random facts
# BWA uses random hit for multimapers

# example parameters
# TAG=E7_75_rep1
# path to fq folders
# ATACDIR=/lustre/scratch126/cellgen/cellgeni/alexp/GSE205117/fastqs
# CDNADIR=/lustre/scratch126/cellgen/cellgeni/alexp/GSE205117/fastqs

# string to filter (grep on full path) fqs for specific sample
# ATACTAG='GSM6205429'
# CDNATAG='GSM6205418'
# expects following files for each run
# ATAC: R1 and R3 - bio; R2 - barcodes;  
# cDNA: R1 - bio; R2 - barcodes; R3 UMI (pos 0-10)
# barcodes shouls be at 15-23,53-61,91-99
# SPECIES='mouse'

# parameters
TAG=$1
ATACDIR=$2
CDNADIR=$3
ATACTAG=$4
CDNATAG=$5
SPECIES=$6

#rm -r $TAG
mkdir $TAG
cd $TAG

exec > Log

echo "SNARE-seq2 pipeline by Pasha Mazin"
echo "email: pm19@sanger.ac.uk"
echo Time: `date`
echo "Parameters: $1 $2 $3 $4 $5 $6"



# paths to bins and scripts
STAR=/nfs/users/nfs_p/pm19/nfs/ghub/STAR/source/STAR
BWA="/nfs/cellgeni/pasham/bin/bwa-mem2/bwa-mem2 mem"
PIPELINEHOME=`echo "$(cd "$(dirname "$0")"; pwd)"`
#/nfs/cellgeni/pasham/projects/2110.SNARE-seq/src/pipeline/src

CPUS=16

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

# these two barcodes are parallele, nth from cnda matches nth from atac
CDNABC=$PIPELINEHOME/../barcodes/barcode.10x.gex.737K-arc-v1.txt
ATACBC=$PIPELINEHOME/../barcodes/barcode.10x.737K-arc-v1.rev.txt

# init binaries
source activate ngs # should include samtools
extractbc="python $PIPELINEHOME/extractBCnUMI.py"
bindfa="python $PIPELINEHOME/bind.fastqs.py"
addcom="python $PIPELINEHOME/fastqs2comment.py"
picard="java -jar /nfs/cellgeni/pasham/bin/picard.jar"

# preprocess barcodes and UMIs
# _extract barcodes from R2
# _ATAC
echo "process ATAC barcodes:"
echo Time: `date`
for i in `ls -1 $ATACDIR/*R2_001.fastq.gz | grep $ATACTAG`
do
  f=`basename $i | sed s/_R2_001.fastq.gz//`
  echo $f
  $extractbc $ATACDIR/${f}_R2_001.fastq.gz 0-16 ${ATACBC} 2> atac_${f}_bc_001.json | gzip > atac_${f}_bc_001.fastq.gz
done

# map RNA
echo "MAP cDNA:"
echo Time: `date`
CDNAOUT=cDNA
UMILEN=12   ## 10x v2
STR=Forward #Reverse ## 3' 10x
BAM="--outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 2 --limitBAMsortRAM 60000000000 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM CB UB GX GN CR CY UR UY"
CBLEN=16 
UMIstart=`expr $CBLEN + 1`

mkdir $CDNAOUT
## for multiple fastq files; change grep options according to your file format
R1=`find $CDNADIR/* | grep $CDNATAG | grep "_R1_" | sort | tr '\n' ','`
R2=`find $CDNADIR/* | grep $CDNATAG | grep "_R2_" | sort | tr '\n' ','`

echo "cDNA Reads:"
echo $R1
echo $R2

$STAR --runThreadN $CPUS --genomeDir $STARREF --readFilesIn $R2 $R1 --runDirPerm All_RWX --readFilesCommand zcat \
     --soloType CB_UMI_Simple --soloCBwhitelist $CDNABC --soloBarcodeReadLength 1 --soloUMIlen $UMILEN --soloStrand $STR \
     --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
     --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
     --soloCBlen $CBLEN --soloUMIstart $UMIstart\
     $BAM \
     --outFileNamePrefix $CDNAOUT/ \
     --soloFeatures Gene GeneFull Velocyto --soloOutFileNames output/ genes.tsv barcodes.tsv matrix.mtx


# add commentary to ATAC fq to preserve BC
echo "add ATAC barcodes to fq:"
echo Time: `date`
for i in `ls -1 $ATACDIR/*R1_001.fastq.gz | grep $ATACTAG | sed s/_R1_001.fastq.gz//`
do
  f=`basename $i`
  echo $f
  $addcom  $ATACDIR/${f}_R1_001.fastq.gz atac_${f}_bc_001.fastq.gz 'CB:Z:' | gzip > atac_${f}_R1c_001.fastq.gz &
  $addcom  $ATACDIR/${f}_R3_001.fastq.gz atac_${f}_bc_001.fastq.gz 'CB:Z:' | gzip > atac_${f}_R3c_001.fastq.gz &
  wait
done

# map ATAC
echo "map ATAC:"
echo Time: `date`
mkdir ATAC

R1=`find ./* | grep atac | grep "_R1c_" | sort | tr '\n' ' '`
R2=`find ./* | grep atac | grep "_R3c_" | sort | tr '\n' ' '`

R1="'""< zcat $R1""'"
R2="'""< zcat $R2""'"
echo "ATAC Reads:"
echo $R1
echo $R2

# consider to use -M (see https://bio-bwa.sourceforge.net/bwa.shtml)
cmd="$BWA -t $CPUS -C $BWAREF $R1 $R2 2> ATAC/bwa.log | samtools view -b > ATAC/out.bam"
eval $cmd

samtools sort --threads $CPUS -m 3G -o ATAC/out.sorted.bam ATAC/out.bam
samtools index ATAC/out.sorted.bam

$picard MarkDuplicates -I ATAC/out.sorted.bam -O ATAC/out.sorted_markdupl.bam -M ATAC/markdupl_metrics.txt --VERBOSITY WARNING --BARCODE_TAG CB
samtools index ATAC/out.sorted_markdupl.bam

samtools idxstats ATAC/out.sorted_markdupl.bam > ATAC/out.sorted_markdupl.idxstats.txt
samtools stats -@ $CPUS ATAC/out.sorted_markdupl.bam > ATAC/out.sorted_markdupl.stats.txt
rm -f ATAC/out.bam ATAC/out.sorted.bam
rm -f ./*.fastq.gz

# run archR
echo "run archR:"
echo Time: `date`
$PIPELINEHOME/run.archr.R . $SPECIES $CPUS $ATACBC $CDNABC 

cd ..

echo "All done!"
echo Time: `date`
