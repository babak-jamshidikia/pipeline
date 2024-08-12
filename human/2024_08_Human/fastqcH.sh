#!/bin/bash
#SBATCH  --job-name=fastqcAT1
#SBATCH  --ntasks=1
#SBATCH  --cpus-per-task=4
#SBATCH  --mem-per-cpu=64g


fastqc -o /home/bjamshidkia/bprj/2024_08_Human/siC/fastqc/siC_TM_3  /prj/2024_08_Human/raw_data/SiC_TM_3_2.fq.gz

