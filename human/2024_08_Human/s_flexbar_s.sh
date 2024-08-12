#!/bin/bash

# @Author: Babak jamshidikia
# @Email:  b.jamshidikia@gmail.com


#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH -J "flexbar single"


# check if we have 3 arguments
if [ ! $# == 3 ]; then
  echo "Usage: $0 [Read 1 file] [target dir e.g. /tmp/]"
  exit
fi

echo "start"
# $1 -> Read 1
# $2 -> Target directory

echo  "$1  read =  " $1
echo  "$2 target = " $2

# run on 10 CPUs
# compress with bz2
# only 30nt or longer
# no uncalled bases
# quality min phred 28
# use sanger quality values (i.e. Illumina 1.9+ encoding)

flexbar -r $1 -p $2-t $3 -n 15 -z GZ -m 30 -u 0  -q TAIL -qt 28 -as "AGATCGGAAGAG" -qf sanger -j

echo "end"
