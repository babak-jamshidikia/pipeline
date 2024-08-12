

all: flexbar bowtie  STAR
#all: bowtie

#Rules for fastqc
#fastqc: makingFastqcdirectory /mnt/c/row_data/fastqc/sub35_1_fastqc.html /mnt/c/row_data/fastqc/sub35_2_fastqc.html
#makingFastqcdirectory:
#	echo "making fastqc directory"
#	mkdir /mnt/c/row_data/fastqc
#/mnt/c/row_data/fastqc: /mnt/c/row_data/fastqc/sub35_1_fastqc.html /mnt/c/row_data/fastqc/sub35_2_fastqc.html
#/mnt/c/row_data/fastqc/sub35_1_fastqc.html:
#	echo "making fastqc for the row file output is in /mnt/c/row_data/fastqc/sub35_1_fastqc.html  "
#	fastqc -t 2 /mnt/c/row_data/sub35_1.fastq -o /mnt/c/row_data/fastqc 
#/mnt/c/row_data/fastqc/sub35_2_fastqc.html:
#	echo "making fastqc for the row file output is in /mnt/c/row_data/fastqc/sub35_2_fastqc.html  "
#	fastqc -t 2 /mnt/c/row_data/sub35_2.fastq -o /mnt/c/row_data/fastqc 

#Rules for flexbar


#----

#----



flexbar: making_flexbar_directory /prj/2024_05_Heart/flexbar/flb_WT_Atria_1.fastq.gz /prj/2024_05_Heart/flexbar/flb_WT_Atria_2.fastq.gz /prj/2024_05_Heart/flexbar/flb_WT_Atria_3.fastq.gz \
 /prj/2024_05_Heart/flexbar/flb_WT_Atria_1.fastq.gz /prj/2024_05_Heart/flexbar/flb_WT_Atria_2.fastq.gz /prj/2024_05_Heart/flexbar/flb_WT_Atria_3.fastq.gz \
	 /prj/2024_05_Heart/flexbar/flb_WT_Vent_1.fastq.gz /prj/2024_05_Heart/flexbar/flb_WT_Vent_2.fastq.gz /prj/2024_05_Heart/flexbar/flb_WT_Vent_3.fastq.gz

making_flexbar_directory: 
	
	echo "making directory for flexbar"
	echo "     "
	
	mkdir /prj/2024_05_Heart/flexbar 

#Atria
/prj/2024_05_Heart/flexbar/flb_WT_Atria_1.fastq.gz:
	$(eval jidAt1fl:= $(jidAt1fl)$(shell sbatch --parsable ~/bioscripts/2024_05_Heart/s_flexbar_s.sh /prj/2024_05_Heart/raw_data/WT_Atria_1.fastq.gz  /prj/2024_05_Heart/flexbar/flb_WT_Atria_1))
	echo " jid_Atriom_1fl ->  $(jidAt1fl)"
/prj/2024_05_Heart/flexbar/flb_WT_Atria_2.fastq.gz:
	$(eval jidAt2fl:=$(jidAt2fl)$(shell sbatch  --parsable ~/bioscripts/2024_05_Heart/s_flexbar_s.sh /prj/2024_05_Heart/raw_data/WT_Atria_2.fastq.gz  /prj/2024_05_Heart/flexbar/flb_WT_Atria_2))
	echo " jid_Atriom_2fl -> $(jidAt2fl)"
/prj/2024_05_Heart/flexbar/flb_WT_Atria_3.fastq.gz:
	$(eval jidAt3fl:=$(jidAt3fl)$(shell sbatch  --parsable ~/bioscripts/2024_05_Heart/s_flexbar_s.sh /prj/2024_05_Heart/raw_data/WT_Atria_3.fastq.gz  /prj/2024_05_Heart/flexbar/flb_WT_Atria_3))
	echo " jid_Atriom_3fl -> $(jidAt3fl)"
# Vents
/prj/2024_05_Heart/flexbar/flb_WT_Vent_1.fastq.gz:
	$(eval jidVe1fl:=$(jidVe1fl)$(shell sbatch --parsable ~/bioscripts/2024_05_Heart/s_flexbar_s.sh /prj/2024_05_Heart/raw_data/WT_Vent_1.fastq.gz /prj/2024_05_Heart/flexbar/flb_WT_Vent_1))
	echo " jid_Vent_1fl -> $(jidVe1fl)"
/prj/2024_05_Heart/flexbar/flb_WT_Vent_2.fastq.gz:
	$(eval jidVe2fl:=$(jidVe2fl)$(shell sbatch  --parsable ~/bioscripts/2024_05_Heart/s_flexbar_s.sh /prj/2024_05_Heart/raw_data/WT_Vent_2.fastq.gz /prj/2024_05_Heart/flexbar/flb_WT_Vent_2))
	echo " jid_Vent_2fl -> $(jidVe2fl)"
/prj/2024_05_Heart/flexbar/flb_WT_Vent_3.fastq.gz:
	$(eval jidVe3fl:=$(jidVe3fl)$(shell sbatch --parsable ~/bioscripts/2024_05_Heart/s_flexbar_s.sh /prj/2024_05_Heart/raw_data/WT_Vent_3.fastq.gz /prj/2024_05_Heart/flexbar/flb_WT_Vent_3))
	echo " jid_Vent_3fl -> $(jidVe3fl)"


#Rules for bowtie
#-------------
#sbatch s_bowtie_s.sh /biodb/genomes/mus_musculus/rRNA/mmusculus.rRNA ~/bprj/2024_05_heart/Atria/Atria_1/flexbar/flb_WT_Atria_1.fastq.gz ~/bprj/2024_05_heart/Atria/Atria_1/bwt/bwt_WT_Atria_1

#-------------





bowtie: makingbowtiefolder /prj/2024_05_Heart/bowtie/bwt_WT_Atria_1.fastq.gz  /prj/2024_05_Heart/bowtie/bwt_WT_Atria_2.fastq.gz /prj/2024_05_Heart/bowtie/bwt_WT_Atria_3.fastq.gz \
 	/prj/2024_05_Heart/bowtie/bwt_WT_Vent_1.fastq.gz /prj/2024_05_Heart/bowtie/bwt_WT_Vent_2.fastq.gz /prj/2024_05_Heart/bowtie/bwt_WT_Vent_3.fastq.gz

makingbowtiefolder:
	mkdir /prj/2024_05_Heart/bowtie

#Atria
/prj/2024_05_Heart/bowtie/bwt_WT_Atria_1.fastq.gz:
	$(eval jidAt1bw:=$(jidAt1bw)$(shell sbatch  --parsable --dependency=afterok:$(jidAt1fl) ~/bioscripts/2024_05_Heart/s_bowtie_s.sh /biodb/genomes/mus_musculus/rRNA/mmusculus.rRNA \
 	/prj/2024_05_Heart/flexbar/flb_WT_Atria_1.fastq.gz /prj/2024_05_Heart/bowtie/bwt_WT_Atria_1.fastq.gz))
	echo " jid_Atria_1bw -> $(jidAt1bw)"

/prj/2024_05_Heart/bowtie/bwt_WT_Atria_2.fastq.gz:
	$(eval jidAt2bw:=$(jidAt2bw)$(shell sbatch  --parsable --dependency=afterok:$(jidAt2fl) ~/bioscripts/2024_05_Heart/s_bowtie_s.sh /biodb/genomes/mus_musculus/rRNA/mmusculus.rRNA \
	/prj/2024_05_Heart/flexbar/flb_WT_Atria_2.fastq.gz /prj/2024_05_Heart/bowtie/bwt_WT_Atria_2.fastq.gz))
	echo " jid_Atria_2bw -> $(jidAt2bw)"
/prj/2024_05_Heart/bowtie/bwt_WT_Atria_3.fastq.gz:
	$(eval jidAt3bw :=$(jidAt3bw)$(shell sbatch  --parsable --dependency=afterok:$(jidAt3fl) ~/bioscripts/2024_05_Heart/s_bowtie_s.sh /biodb/genomes/mus_musculus/rRNA/mmusculus.rRNA \
	 /prj/2024_05_Heart/flexbar/flb_WT_Atria_3.fastq.gz /prj/2024_05_Heart/bowtie/bwt_WT_Atria_3.fastq.gz))
	echo " jid_Atria_3bw -> $(jidAt3bw)"
# Vents -----------------------------------

/prj/2024_05_Heart/bowtie/bwt_WT_Vent_1.fastq.gz:
	$(eval jidVe1bw:=$(jidVe1bw)$(shell sbatch  --parsable --dependency=afterok:$(jidVe1fl) ~/bioscripts/2024_05_Heart/s_bowtie_s.sh /biodb/genomes/mus_musculus/rRNA/mmusculus.rRNA \
	 /prj/2024_05_Heart/flexbar/flb_WT_Vent_1.fastq.gz /prj/2024_05_Heart/bowtie/bwt_WT_Vent_1.fastq.gz))
	echo " jid_Vent_1bw -> $(jidVe1bw)"
/prj/2024_05_Heart/bowtie/bwt_WT_Vent_2.fastq.gz:
	$(eval jidVe2bw:=$(jidVe2bw)$(shell sbatch  --parsable --dependency=afterok:$(jidVe2fl)  ~/bioscripts/2024_05_Heart/s_bowtie_s.sh /biodb/genomes/mus_musculus/rRNA/mmusculus.rRNA \
	 /prj/2024_05_Heart/flexbar/flb_WT_Vent_2.fastq.gz /prj/2024_05_Heart/bowtie/bwt_WT_Vent_2.fastq.gz))
	echo " jid_Vent_2bw -> $(jidVe2bw)"
/prj/2024_05_Heart/bowtie/bwt_WT_Vent_3.fastq.gz:
	$(eval jidVe3bw:=$(jidVe3bw)$(shell sbatch  --parsable --dependency=afterok:$(jidVe3fl)  ~/bioscripts/2024_05_Heart/s_bowtie_s.sh /biodb/genomes/mus_musculus/rRNA/mmusculus.rRNA \
	 /prj/2024_05_Heart/flexbar/flb_WT_Vent_3.fastq.gz /prj/2024_05_Heart/bowtie/bwt_WT_Vent_3.fastq.gz))
	echo " jid_Vent_3bw -> $(jidVe3bw)"



#Rules for STAR
#----------------------
#sbatch s_star_s.sh /biodb/genomes/mus_musculus/GRCm39_107/star ~/bprj/2024_05_heart/Atria/Atria_1/bwt/bwt.fastq.gz ~/bprj/2024_05_heart/Atria/Atria_1/star /biodb/genomes/mus_musculus/GRCm39_107/GRCm39.107.gtf

#---------------------



STAR: makingSTARfolder /prj/2024_05_Heart/star/str_WT_Atria_1.sortedBycoord.out.bam /prj/2024_05_Heart/star/str_WT_Atria_2.sortedBycoord.out.bam /prj/2024_05_Heart/star/str_WT_Atria_3.sortedBycoord.out.bam \
	/prj/2024_05_Heart/star/str_WT_Vent_1.sortedBycoord.out.bam /prj/2024_05_Heart/star/str_WT_Vent_2.sortedBycoord.out.bam /prj/2024_05_Heart/star/str_WT_Vent_3.sortedBycoord.out.bam
makingSTARfolder:
	mkdir  /prj/2024_05_Heart/STAR

/prj/2024_05_Heart/star/str_WT_Atria_1.sortedBycoord.out.bam:
	$(eval jidAt1st:=$(jidAt1st)$(shell sbatch --parsable --dependency=afterok:$(jidAt1bw)  ~/bioscripts/2024_05_Heart/s_star_s.sh /biodb/genomes/mus_musculus/GRCm39_107/star \
	 /prj/2024_05_Heart/bowtie/bwt_WT_Atria_1.fastq.gz  /prj/2024_05_Heart/star/str_WT_Atria_1 /biodb/genomes/mus_musculus/GRCm39_107/GRCm39.107.gtf))
	echo " jid_Atria_1st -> $(jidAt1st)"

/prj/2024_05_Heart/star/str_WT_Atria_2.sortedBycoord.out.bam:
	$(eval jidAt2st:=$(jidAt2st)$(shell sbatch --parsable --dependency=afterok:$(jidAt2bw) ~/bioscripts/2024_05_Heart/s_star_s.sh /biodb/genomes/mus_musculus/GRCm39_107/star \
	 /prj/2024_05_Heart/bowtie/bwt_WT_Atria_2.fastq.gz  /prj/2024_05_Heart/star/str_WT_Atria_2 /biodb/genomes/mus_musculus/GRCm39_107/GRCm39.107.gtf))
	echo " jid_Atria_2st -> $(jidAt2st)"

/prj/2024_05_Heart/star/str_WT_Atria_3.sortedBycoord.out.bam:
	$(eval jidAt3st:=$(jidAt3s)$(shell sbatch  --parsable --dependency=afterok:$(jidAt3bw) ~/bioscripts/2024_05_Heart/s_star_s.sh /biodb/genomes/mus_musculus/GRCm39_107/star \
	 /prj/2024_05_Heart/bowtie/bwt_WT_Atria_3.fastq.gz  /prj/2024_05_Heart/star/str_WT_Atria_3 /biodb/genomes/mus_musculus/GRCm39_107/GRCm39.107.gtf))
	echo " jid_Atria_3st -> $(jidAt3st)"

/prj/2024_05_Heart/star/str_WT_Vent_1.sortedBycoord.out.bam:
	$(eval jidVe1st:=$(jidVe1st)$(shell sbatch  --parsable --dependency=afterok:$(jidVe1bw) ~/bioscripts/2024_05_Heart/s_star_s.sh /biodb/genomes/mus_musculus/GRCm39_107/star \
	 /prj/2024_05_Heart/bowtie/bwt_WT_Vent_1.fastq.gz  /prj/2024_05_Heart/star/str_WT_Vent_1 /biodb/genomes/mus_musculus/GRCm39_107/GRCm39.107.gtf))
	echo " jid_Vent_1st -> $(jidVe1st)"

/prj/2024_05_Heart/star/str_WT_Vent_2.sortedBycoord.out.bam:
	$(evak  jidVe2st:=$(jidVe2st)$(shell sbatch  --parsable --dependency=afterok:$(jidVe2bw) ~/bioscripts/2024_05_Heart/s_star_s.sh /biodb/genomes/mus_musculus/GRCm39_107/star \
	 /prj/2024_05_Heart/bowtie/bwt_WT_Vent_2.fastq.gz  /prj/2024_05_Heart/star/str_WT_Vent_2 /biodb/genomes/mus_musculus/GRCm39_107/GRCm39.107.gtf))
	echo " jid_Vent_2st -> $(jidVe2st)"

/prj/2024_05_Heart/star/str_WT_Vent_3.sortedBycoord.out.bam:
	$(evak  jidVe3st:=$(jidVe3st)$(shell sbatch  --parsable --dependency=afterok:$(jidVe3bw) ~/bioscripts/2024_05_Heart/s_star_s.sh /biodb/genomes/mus_musculus/GRCm39_107/star \
	 /prj/2024_05_Heart/bowtie/bwt_WT_Vent_3.fastq.gz  /prj/2024_05_Heart/star/str_WT_Vent_3 /biodb/genomes/mus_musculus/GRCm39_107/GRCm39.107.gtf))
	echo " jid_Vent_3st -> $(jidVe3st)"

#TRACING  


# clean:
# 	rm -r /mnt/c/row_data/fastqc
# 	rm -r /mnt/c/row_data/flexbar 
# 	rm -r /mnt/c/row_data/bowtie 
# 	sudo rm -r /home/STOUT/temp1
# 	sudo rm -r /home/STOUT/temp2
