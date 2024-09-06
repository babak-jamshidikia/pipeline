#!/bin/bash

# @Author: Babak jamshidikia
# @Email:  b.jamshidikia@gmail.com


#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=10G
#SBATCH -J "flexbar single"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=bjamshidkia@arizona.edu


# check if we have 3 arguments
if [ ! $# == 4 ]; then
  echo "Usage: $0 [Read 1 file] [target dir e.g. /tmp/]"
  exit
fi

echo "start"
# $1 -> Read 1
# $2 -> Read 2
# $3 -> Target directory

echo  "$1  read1 =  " $1
echo  "$2  read2 = " $2
echo  "target directory = " $3
echo  "target name =  " $4  

target=$4
# run on 10 CPUs
# compress with bz2
# only 30nt or longer
# no uncalled bases
# quality min phred 28
# use sanger quality values (i.e. Illumina 1.9+ encoding)

flexbar -r $1 -p $2 -t $3/$target -n 16 -z GZ -m 30 -u 0  -q TAIL -qt 28 -as "AGATCGGAAGAG" -qf sanger -j

echo "end"
