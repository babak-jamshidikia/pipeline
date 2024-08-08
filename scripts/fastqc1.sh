#!/bin/bash
#SBATCH  --job-name=fastqcAT1
#SBATCH  --ntasks=1
#SBATCH  --cpus-per-task=4
#SBATCH  --mem-per-cpu=64g


fastqc -o /home/bjamshidkia/bprj/2024_05_heart/vent/vent_3/fastqc  /prj/2024_05_Heart/raw_data/WT_Vent_3.fastq.gz

