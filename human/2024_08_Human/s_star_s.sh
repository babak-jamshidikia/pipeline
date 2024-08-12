#!/bin/bash

# @Author: babak jamshidikia

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH -J "STAR genome alignment single"


# check if we have 6 arguments
if [ ! $# == 4 ]; then
  echo "Usage: $0 [STAR index] [Read 1 file] [Read 2 file] [target dir e.g. /tmp/] [Read 1 marker, e.g. R1] [GTF file]"
  exit
fi

echo "----- START ----- "

echo " Genome index ------> " $1
echo " Read 1 ----------->  " $2
echo " Target directory ->  " $3
echo " GTFfile ---------->  " $4


# $1 -> Genome index
# $2 -> Read 1
# $3 -> Target directory
# $4 -> GTFfile

# remove the file extension and potential "R1" markings
# (works for double extension, e.g. .fastq.gz)
#target=`expr ${2/$5/} : \(.*\)\..*\.'`

# create the target directory, STAR will not do that for us
#mkdir $3/$target
STAR --genomeDir $1 --runThreadN 5 --readFilesIn $2 $3  --sjdbGTFfile $5 --readFilesCommand\
 zcat --outFileNamePrefix $4 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --genomeLoad NoSharedMemory 
echo "----- END --------------- "
