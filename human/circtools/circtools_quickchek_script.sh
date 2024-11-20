#! /bin/bash

#@Babak jamshidikia
#@email : bjamshidkia@arizona.edu
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=60G
#SBATCH -J "circtools"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=bjamshidkia@arizona.edu 


#circtools quickcheck -d 01_detect/ -s ../star  -l minus,plus -g 1,2,1,2,1,2,1,2  -o 02_quickcheck/  -C .Chimeric.out.junction


circtools quickcheck -d /prj/2024_08_Human/circtools_detect/01_detect/ -s  /prj/2024_08_Human/circresult/  -l minus,plus -g 1,1,1,1,2,2,2,2\
  -o /prj/2024_08_Human/circtools_detect/02_quickcheck  -C .Chimeric.out.junction
