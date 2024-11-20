#! /bin/bash


#/home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_star_pair.sh  /biodb/genomes/homo_sapiens >
#        /prj/2024_08_Human/bowtie/SiC_1_1.fastq.gz  /prj/2024_08_Human/bowtie/SiC_1_2.fastq.gz  /prj/2024_08_Human/star/SiC_1  /biodb/genomes/homo_sapiens/GRCh38_107/GRCh38.107.gtf

sbatch --parsable /home/bjamshidkia/repositories/pipeline/human/circtools/slurm_circtools_mapping.sh /biodb/genomes/homo_sapiens/GRCh38_107/star /SiC_1_1.fastq.gz\
  SiC_1_2.fastq.gz /prj/2024_08_Human/circresult SiC_1 /biodb/genomes/homo_sapiens/GRCh38_107/GRCh38.107.gtf

