#! /bin/bash

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=60G
#SBATCH -J "circtools detect"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=bjamshidkia@arizona.edu 


#circtools detect @samplesheet \
#      -mt1 @mate1 \
#      -mt2 @mate2 \
#      -B @bam_files.txt \
#      -D \
#      -R mouse_repeats.gtf \
#      -an ../../Mus_musculus.GRCm38.90.gtf \
#      -Pi \
#      -F \
#      -M \
#      -Nr 5 6 \
#      -fg \
#      -G \
#      -A Mus_musculus.GRCm38.dna.primary_assembly.fa


circtools detect @/prj/2024_08_Human/circtools_detect/samplesheet \
      -mt1 @/prj/2024_08_Human/circtools_detect/mate1 \
      -mt2 @/prj/2024_08_Human/circtools_detect/mate2 \
      -B @/prj/2024_08_Human/circtools_detect/bam_files.txt \
      -D \
      -R /prj/2024_08_Human/GTFrepeats.gtf \
      -an /biodb/genomes/homo_sapiens/GRCh38_107/GRCh38.107.gtf \
      -Pi \
      -F \
      -M \
      -Nr 5 6 \
      -fg \
      -G \
      -A /biodb/genomes/homo_sapiens/GRCh38_107/GRCh38_107.fa\
      -T 20\
      -O /prj/2024_08_Human/circtools_detect/01_detect
