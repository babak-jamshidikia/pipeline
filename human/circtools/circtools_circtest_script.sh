#! /bin/bash

#Babak jamshidikia
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=60G
#SBATCH -J "circtools"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=bjamshidkia@arizona.edu 


#circtools quickcheck -d 01_detect/ -s ../star  -l minus,plus -g 1,2,1,2,1,2,1,2  -o 02_quickcheck/  -C .Chimeric.out.junction

circtools circtest -d /prj/2024_08_Human/circtools_detect/01_detect/ -p 0.01 -s 2 -r 4 -C 2 -g 1,1,1,1,2,2,2,2\
 -l SiC,SiC_TM -c 4,5,6,7,8,9,10,11 -o /prj/2024_08_Human/circtools_detect/04_circtest/



