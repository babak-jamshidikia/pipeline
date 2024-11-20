#! /bin/bash

# @author : Babak Jamshidikia
# @email  : bjamshidkia@arizona.edu
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=60G
#SBATCH -J "circtools"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=bjamshidkia@arizona.edu 



#circtools primex -d DCC/CircCoordinates 
#		  -f Mus_musculus.GRCm38.dna.primary_assembly.fa 
#		  -g Mus_musculus.GRCm38.90.gtf 
#		  -O mm 
#		  -G Ryr2 
#		  -T "Ryr2 primer"

circtools primex -d /prj/2024_08_Human/circtools_detect/01_detect/CircCoordinates\
		 -f /biodb/genomes/homo_sapiens/GRCh38_107/GRCh38_107.fa\
		 -g /biodb/genomes/homo_sapiens/GRCh38_107/GRCh38.107.gtf\
		 -O hs\
		 -G FARSA\
		 -T "FARSA primer"\
		 -o /prj/2024_08_Human/circtools_detect/05_primex/
