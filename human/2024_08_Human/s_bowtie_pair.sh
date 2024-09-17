#!/bin/bash

# @Author: Babak jamshidikia <bjamshidikia>
# @Email:  jamshidikia@gmaiil.com
# @Last modified by:   bjamshidikia
# @Last modified time: Friday, May 6, 2016 4:17 PM
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -J "bowtie2 rRNA filtering single"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=bjamshidkia@arizona.edu



# check if we have 5 arguments
if [ ! $# == 5 ]; then
  echo "Usage: $0 [rRNA index argument] [Read 1 file] [target dir e.g. /awesome/project/]"
  exit
fi

# $1 -> rRNA index
# $2 -> Read 1
# $3 -> Target directory

echo "start"

echo "index  = " $1
echo "read1  = " $2
echo "read2  = " $3
echo "target directory = " $4
echo "target file = " $5


# remove the file extension and potential "R1" markings
# (works for double extension, e.g. .fastq.gz)
#target = $5 #`expr ${2//} : \(.*\)\..*\.'`
#echo $target
#exit

# run on 20 CPUs
# set fixed seed
# memory mapped IO for multiple instances
# display timing information
# write gz unmapping reads [== no rRNA] to target dir

#bowtie2 -x $1 -U $2  -S /dev/null --no-unal --omit-sec-seq --threads 20 --mm --seed 1337 --time --un-gz $3/bwt.fastq.gz 2> $3/bwt.log
#bowtie2 -x $1 -U $2  -S /dev/null --no-unal --omit-sec-seq --threads 20 --mm --seed 1337 --time --un-gz $3 2> ${3}.log
#bowtie2 -x $1 -1 $2 -2 $3 -S /dev/null --no-unal --omit-sec-seq --threads 20 --mm --seed 1337 --time --un-conc-gz $4/$target.fastq.gz  2> $4/bwt.log
bowtie2 -x $1 -1 $2 -2 $3 -S /dev/null --no-unal --omit-sec-seq --threads 16 --mm --seed 1337 --time --un-conc-gz ${4}/$5_%.fastq.gz 2> $4/$5.log


echo "end"