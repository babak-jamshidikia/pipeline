
all:   flexbar  bowtie  STAR
#all: STAR

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


flexbar: making_flexbar_directory  /prj/2024_08_human/flexbar/SiC_1.fastq.gz 
 /prj/2024_08_human/flexbar/SiC_2.fastq.gz /prj/2024_08_human/flexbar/SiC_3.fastq.gz   /prj/2024_08_human/flexbar/SiC_4.fastq.gz  \
        /prj/2024_08_human/flexbar/SiC_TM_1.fastq.gz /prj/2024_08_human/flexbar/SiC_TM_2.fastq.gz \
 /prj/2024_08_human/flexbar/SiC_TM_3.fastq.gz  /prj/2024_08_human/flexbar/SiC_TM_4.fastq.gz

making_flexbar_directory: 
	
	echo "making directory for flexbar"
	echo "     "
	
	mkdir /prj/2024_08_Human/flexbar 

#SiC
/prj/2024_08_human/flexbar/SiC_1.fastq.gz:
	$(eval jidSiC1fl:= $(jidSiC1fl)$(shell sbatch --parsable /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_flexbar_pair.sh /prj/2024_08_Human/raw_data/SiC_1_1.fq.gz \
                                                     /prj/2024_08_Human/raw_data/SiC_1_2.fq.gz  /prj/2024_08_Human/flexbar/SiC_1))
	echo " jid_SiC_1fl ->  $(jidSiC1fl)"
/prj/2024_08_human/flexbar/SiC_2.fastq.gz:
	$(eval jidSiC2fl:=$(jidSiC2fl)$(shell sbatch  --parsable  /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_flexbar_pair.sh /prj/2024_08_Human/raw_data/SiC_2_1.fq.gz \
                                                     /prj/2024_08_Human/raw_data/SiC_2_2.fq.gz  /prj/2024_08_Human/flexbar/SiC_2))
	echo " jid_SiC_2fl -> $(jidSiC2fl)"
/prj/2024_08_human/flexbar/SiC_3.fastq.gz:
	$(eval jidSiC3fl:=$(jidSiC3fl)$(shell sbatch  --parsable /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_flexbar_pair.sh /prj/2024_08_Human/raw_data/SiC_3_1.fq.gz \
                                                    /prj/2024_08_Human/raw_data/SiC_3_2.fq.gz  /prj/2024_08_Human/flexbar/SiC_3))
	echo " jid_SiC_3fl -> $(jidSiC3fl)"

/prj/2024_08_human/flexbar/SiC_4.fastq.gz:
	$(eval jidSiC4fl:=$(jidSiC4fl)$(shell sbatch  --parsable /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_flexbar_pair.sh /prj/2024_08_Human/raw_data/SiC_4_1.fq.gz \
                                                    /prj/2024_08_Human/raw_data/SiC_4_2.fq.gz  /prj/2024_08_Human/flexbar/SiC_4))
	echo " jid_SiC_4fl -> $(jidSiC4fl)"

# SiC_TM

/prj/2024_08_human/flexbar/SiC_TM_1.fastq.gz:
	$(eval jidSiCTM1fl:=$(jidSiCTM1fl)$(shell sbatch --parsable /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_flexbar_pair.sh /prj/2024_08_Human/raw_data/SiC_TM_1_1.fq.gz \
                                                    /prj/2024_08_Human/raw_data/Sic/SiC_TM_1_2.fq.gz  /prj/2024_08_Human/flexbar/SiC/SiC_TM_1))
	echo " jid_SIC_TM_1fl -> $(jidSiCTM1fl)"
/prj/2024_08_human/flexbar/SiC_TM_2.fastq.gz:
	$(eval jidSiCTM2fl:=$(jidSiCTM2fl)$(shell sbatch  --parsable  /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_flexbar_pair.sh /prj/2024_08_Human/raw_data/SiC_TM_2_1.fq.gz \
                                                    /prj/2024_08_Human/raw_data/Sic/SiC_TM_2_2.fq.gz  /prj/2024_08_human/flexbar/SiC/SiC_TM_2))
	echo " jid_SiC_TM_2fl -> $(jidSiCTM2fl)"
/prj/2024_08_human/flexbar/SiC_TM_3.fastq.gz:
	$(eval jidSiCTM3fl:=$(jidSiCTM3fl)$(shell sbatch --parsable  /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_flexbar_pair.sh /prj/2024_08_Human/raw_data/SiC_TM_3_1.fq.gz \
                                                     /prj/2024_08_Human/raw_data/Sic/SiC_TM_3_2.fq.gz  /prj/2024_08_human/flexbar/SiC/SiC_TM_3))
	echo " jid_SiC_TM_3fl -> $(jidSiCTM3fl)"

/prj/2024_08_human/flexbar/SiC_TM_4.fastq.gz:

	$(eval jidSiCTM4fl:=$(jidSiCTM4fl)$(shell sbatch --parsable  /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_flexbar_pair.sh /prj/2024_08_Human/raw_data/SiC_TM_4_1.fq.gz \                                   
                  /prj/2024_08_Human/raw_data/Sic/SiC_TM_4_2.fq.gz  /prj/2024_08_human/flexbar/SiC/SiC_TM_4))
	echo " jid_SiC_TM_4fl -> $(jidSiCTM4fl)"



#Rules for bowtie
#-------------

bowtie: makingbowtiefolder /prj/2024_08_Human/bowtie/bwt_SiC_1_1.fastq.gz /prj/2024_08_Human/bowtie/bwt_SiC_2_1.fastq.gz /prj/2024_08_Human/bowtie/bwt_SiC_3_1.fastq.gz /prj/2024_08_Human/bowtie/bwt_SiC_4_1.fastq.gz \
 /prj/2024_08_Human/bowtie/bwt_SiC_TM_1_1.fastq.gz  /prj/2024_08_Human/bowtie/bwt_SiC_TM_2_1.fastq.gz /prj/2024_08_Human/bowtie/bwt_SiC_TM_3_1.fastq.gz /prj/2024_08_Human/bowtie/bwt_SiC_TM_3_1.fastq.gz

makingbowtiefolder:
	mkdir /prj/2024_08_Human/bowtie

#SiC
/prj/2024_08_Human/bowtie/bwt_SiC_1_1.fastq.gz:  #/prj/2024_08_Human/bowtie/bwt_SiC_1_2.fastq.gz: 

	$(eval jidSiC1bw:=$(jidSiC1bw)$(shell sbatch  --parsable --dependency=afterok:$(jidSiC1fl)   /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_bowtie_pair.sh \
         /biodb/genomes/homo_sapiens/rRNA_cluster/rRNA_precursor \
 	/prj/2024_08_Human/flexbar/SiC_1_1.fastq.gz  /prj/2024_08_Human/flexbar/SiC_1_2.fastq.gz  /prj/2024_08_Human/bowtie/  bwt_SiC_1))
	echo " jid_Atria_1bw -> $(jidSiC1bw)"

/prj/2024_08_Human/bowtie/bwt_SiC_2_1.fastq.gz: #/prj/2024_08_Human/bowtie/bwt_SiC_2_2.fastq.gz

	$(eval jidSiC2bw:=$(jidSiC2bw)$(shell sbatch  --parsable --dependency=afterok:$(jidSiC2fl)   /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_bowtie_pair.sh \
         /biodb/genomes/homo_sapiens/rRNA_cluster/rRNA_precursor \
 	/prj/2024_08_Human/flexbar/SiC_2_1.fastq.gz  /prj/2024_08_Human/flexbar/SiC_2_2.fastq.gz  /prj/2024_08_Human/bowtie/  bwt_SiC_2))
	echo " jid_Atria_1bw -> $(jidSiC2bw)"

/prj/2024_08_Human/bowtie/bwt_SiC_3_1.fastq.gz: # /prj/2024_08_Human/bowtie/bwt_SiC_3_2.fastq.gz

	$(eval jidSiC3bw:=$(jidSiC3bw)$(shell sbatch  --parsable --dependency=afterok:$(jidSiC3fl)   /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_bowtie_pair.sh \
         /biodb/genomes/homo_sapiens/rRNA_cluster/rRNA_precursor \
 	/prj/2024_08_Human/flexbar/SiC_3_1.fastq.gz  /prj/2024_08_Human/flexbar/SiC_3_2.fastq.gz  /prj/2024_08_Human/bowtie/  bwt_SiC_3))
	echo " jid_Atria_1bw -> $(jidSiC3bw)"
/prj/2024_08_Human/bowtie/bwt_SiC_4_1.fastq.gz: #/prj/2024_08_Human/bowtie/bwt_SiC_4_2.fastq.gz

	$(eval jidSiC4bw:=$(jidSiC4bw)$(shell sbatch  --parsable --dependency=afterok:$(jidSiC4fl)   /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_bowtie_pair.sh \
         /biodb/genomes/homo_sapiens/rRNA_cluster/rRNA_precursor \
 	/prj/2024_08_Human/flexbar/SiC_4_1.fastq.gz  /prj/2024_08_Human/flexbar/SiC_4_2.fastq.gz  /prj/2024_08_Human/bowtie/  bwt_SiC_4))
	echo " jid_Atria_1bw -> $(jidSiC4bw)"



# just for test ------------


#	$(eval jidAt1bw:=$(jidAt1bw)$(shell sbatch  --parsable    /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_bowtie_pair.sh \
#         /biodb/genomes/homo_sapiens/rRNA_cluster/rRNA_precursor \
# 	/prj/2024_08_Human/flexbar/SiC_1_1.fastq.gz  /prj/2024_08_Human/flexbar/SiC_1_2.fastq.gz  /prj/2024_08_Human/bowtie/  bwt_SiC_1))

#	echo " jid_Atria_1bw -> $(jidAt1bw)"


# ------------------


# SiC_TM -----------------------------------

/prj/2024_08_Human/bowtie/bwt_SiC_TM_1_1.fastq.gz: #/prj/2024_08_Human/bowtie/bwt_SiC_TM_1_2.fastq.gz: 

	$(eval jidSiCTM1bw:=$(jidSiCTM1bw)$(shell sbatch  --parsable --dependency=afterok:$(jidSiCTM1fl)   /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_bowtie_pair.sh \
         /biodb/genomes/homo_sapiens/rRNA_cluster/rRNA_precursor \
 	/prj/2024_08_Human/flexbar/SiC_TM_1_1.fastq.gz  /prj/2024_08_Human/flexbar/SiC_TM_1_2.fastq.gz  /prj/2024_08_Human/bowtie/  bwt_SiC_TM_1))
	echo " jid_SiC_TM_1bw -> $(jiSiCTM1bw)"

/prj/2024_08_Human/bowtie/bwt_SiC_TM_2_1.fastq.gz: #/prj/2024_08_Human/bowtie/bwt_SiC_TM_2_2.fastq.gz: 

	$(eval jidSiCTM2bw:=$(jidSiCTM2bw)$(shell sbatch  --parsable --dependency=afterok:$(jidSiCTM2fl)   /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_bowtie_pair.sh \
         /biodb/genomes/homo_sapiens/rRNA_cluster/rRNA_precursor \
 	/prj/2024_08_Human/flexbar/SiC_TM_2_1.fastq.gz  /prj/2024_08_Human/flexbar/SiC_TM_2_2.fastq.gz  /prj/2024_08_Human/bowtie/  bwt_SiC_TM_2))
	echo " jid_SiC_TM_1bw -> $(jiSiCTM2bw)"

/prj/2024_08_Human/bowtie/bwt_SiC_TM_3_1.fastq.gz: #/prj/2024_08_Human/bowtie/bwt_SiC_TM_3_2.fastq.gz: 

	$(eval jidSiCTM3bw:=$(jidSiCTM3bw)$(shell sbatch  --parsable --dependency=afterok:$(jidSiCTM3fl)   /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_bowtie_pair.sh \
         /biodb/genomes/homo_sapiens/rRNA_cluster/rRNA_precursor \
 	/prj/2024_08_Human/flexbar/SiC_TM_3_1.fastq.gz  /prj/2024_08_Human/flexbar/SiC_TM_3_2.fastq.gz  /prj/2024_08_Human/bowtie/  bwt_SiC_TM_3))
	echo " jid_SiC_TM_1bw -> $(jiSiCTM3bw)"

/prj/2024_08_Human/bowtie/bwt_SiC_TM_4_1.fastq.gz: #/prj/2024_08_Human/bowtie/bwt_SiC_TM_3_2.fastq.gz: 

	$(eval jidSiCTM4bw:=$(jidSiCTM4bw)$(shell sbatch  --parsable --dependency=afterok:$(jidSiCTM4fl)   /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_bowtie_pair.sh \
         /biodb/genomes/homo_sapiens/rRNA_cluster/rRNA_precursor \
 	/prj/2024_08_Human/flexbar/SiC_TM_4_1.fastq.gz  /prj/2024_08_Human/flexbar/SiC_TM_4_2.fastq.gz  /prj/2024_08_Human/bowtie/  bwt_SiC_TM_4))
	echo " jid_SiC_TM_1bw -> $(jiSiCTM4bw)"


#Rules for STAR
#----------------------
#sbatch s_star_s.sh /biodb/genomes/mus_musculus/GRCm39_107/star ~/bprj/2024_05_heart/Atria/Atria_1/bwt/bwt.fastq.gz ~/bprj/2024_05_heart/Atria/Atria_1/star /biodb/genomes/mus_musculus/GRCm39_107/GRCm39.107.gtf

#---------------------



STAR: makingSTARfolder  /prj/2024_08_Human/star/str_SiC_1Aligned.sortedByCoord.out.bam /prj/2024_08_Human/star/str_SiC_2Aligned.sortedByCoord.out.bam \
/prj/2024_08_Human/star/str_SiC_3Aligned.sortedByCoord.out.bam /prj/2024_08_Human/star/str_SiC_4Aligned.sortedByCoord.out.bam \
/prj/2024_08_Human/star/str_SiC_TM_1Aligned.sortedByCoord.out.bam /prj/2024_08_Human/star/str_SiC_TM_2Aligned.sortedByCoord.out.bam \
/prj/2024_08_Human/star/str_SiC_TM_3Aligned.sortedByCoord.out.bam /prj/2024_08_Human/star/str_SiC_TM_4Aligned.sortedByCoord.out.bam

makingSTARfolder:
	mkdir  /prj/2024_08_Human/star

#---------------- just for test

#/prj/2024_08_Human/star/str_SiC_1.sortedBycoord.out.bam:
#	$(eval jidSiC1_st:=$(jidSiC1st)$(shell sbatch --parsable   /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_star_pair.sh  /biodb/genomes/homo_sapiens/GRCh38_107/star \
#	/prj/2024_08_Human/bowtie/bwt_SiC_1_1.fastq.gz  /prj/2024_08_Human/bowtie/bwt_SiC_1_2.fastq.gz  /prj/2024_08_Human/star/str_SiC_1  /biodb/genomes/homo_sapiens/GRCh38_107/GRCh38.107.chr_patch_hapl_scaff.gtf))
#	echo " jid_Atria_1st -> $(jidSiC1st)"

#-----------------

/prj/2024_08_Human/star/str_SiC_1Aligned.sortedByCoord.out.bam:
	$(eval jidSiC1st:=$(jidSiC1st)$(shell sbatch --parsable --dependency=afterok:$(jidSiC1bw)  /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_star_pair.sh  /biodb/genomes/homo_sapiens/GRCh38_107/star \
	/prj/2024_08_Human/bowtie/bwt_SiC_1_1.fastq.gz  /prj/2024_08_Human/bowtie/bwt_SiC_1_2.fastq.gz  /prj/2024_08_Human/star/str_SiC_1  /biodb/genomes/homo_sapiens/GRCh38_107/GRCh38.107.chr_patch_hapl_scaff.gtf))
	echo " jid_Atria_1st -> $(jidSiC1st)"


/prj/2024_08_Human/star/str_SiC_2Aligned.sortedByCoord.out.bam:
	$(eval jidSiC2st:=$(jidSiC2st)$(shell sbatch --parsable --dependency=afterok:$(jidSiC2bw)  /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_star_pair.sh  /biodb/genomes/homo_sapiens/GRCh38_107/star \
	/prj/2024_08_Human/bowtie/bwt_SiC_2_1.fastq.gz  /prj/2024_08_Human/bowtie/bwt_SiC_2_2.fastq.gz  /prj/2024_08_Human/star/str_SiC_2  /biodb/genomes/homo_sapiens/GRCh38_107/GRCh38.107.chr_patch_hapl_scaff.gtf))
	echo " jid_Atria_1st -> $(jidSiC2st)"

/prj/2024_08_Human/star/str_SiC_3Aligned.sortedByCoord.out.bam:
	$(eval jidSiC3st:=$(jidSiC3st)$(shell sbatch --parsable --dependency=afterok:$(jidSiC3bw)  /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_star_pair.sh  /biodb/genomes/homo_sapiens/GRCh38_107/star \
	/prj/2024_08_Human/bowtie/bwt_SiC_3_1.fastq.gz  /prj/2024_08_Human/bowtie/bwt_SiC_3_2.fastq.gz  /prj/2024_08_Human/star/str_SiC_3  /biodb/genomes/homo_sapiens/GRCh38_107/GRCh38.107.chr_patch_hapl_scaff.gtf))
	echo " jid_Atria_1st -> $(jidSiC3st)"

/prj/2024_08_Human/star/str_SiC_4Aligned.sortedByCoord.out.bam:
	$(eval jidSiC4st:=$(jidSiC4st)$(shell sbatch --parsable --dependency=afterok:$(jidSiC4bw)  /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_star_pair.sh  /biodb/genomes/homo_sapiens/GRCh38_107/star \
	/prj/2024_08_Human/bowtie/bwt_SiC_4_1.fastq.gz  /prj/2024_08_Human/bowtie/bwt_SiC_4_2.fastq.gz  /prj/2024_08_Human/star/str_SiC_3  /biodb/genomes/homo_sapiens/GRCh38_107/GRCh38.107.chr_patch_hapl_scaff.gtf))
	echo " jid_Atria_1st -> $(jidSiC4st)"


#SiC_TM

/prj/2024_08_Human/star/str_SiC_TM_1Aligned.sortedByCoord.out.bam:
	$(eval jidSiCTM1st:=$(jidSiCTM1st)$(shell sbatch --parsable --dependency=afterok:$(jidSiCTM1bw)  /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_star_pair.sh  /biodb/genomes/homo_sapiens/GRCh38_107/star \
	/prj/2024_08_Human/bowtie/bwt_SiC_TM_1_1.fastq.gz  /prj/2024_08_Human/bowtie/bwt_SiC_TM_1_2.fastq.gz  /prj/2024_08_Human/star/str_SiC_TM_1  /biodb/genomes/homo_sapiens/GRCh38_107/GRCh38.107.chr_patch_hapl_scaff.gtf))
	echo " jid_Atria_1st -> $(jidSiCTM1st)"

/prj/2024_08_Human/star/str_SiC_TM_2Aligned.sortedByCoord.out.bam:
	$(eval jidSiCTM2st:=$(jidSiCTM2st)$(shell sbatch --parsable --dependency=afterok:$(jidSiCTM2bw)  /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_star_pair.sh  /biodb/genomes/homo_sapiens/GRCh38_107/star \
	/prj/2024_08_Human/bowtie/bwt_SiC_TM_2_1.fastq.gz  /prj/2024_08_Human/bowtie/bwt_SiC_TM_2_2.fastq.gz  /prj/2024_08_Human/star/str_SiC_TM_2  /biodb/genomes/homo_sapiens/GRCh38_107/GRCh38.107.chr_patch_hapl_scaff.gtf))
	echo " jid_Atria_1st -> $(jidSiCTM2st)"

/prj/2024_08_Human/star/str_SiC_TM_3Aligned.sortedByCoord.out.bam:
	$(eval jidSiCTM3st:=$(jidSiCTM3st)$(shell sbatch --parsable --dependency=afterok:$(jidSiCTM3bw)  /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_star_pair.sh  /biodb/genomes/homo_sapiens/GRCh38_107/star \
	/prj/2024_08_Human/bowtie/bwt_SiC_TM_3_1.fastq.gz  /prj/2024_08_Human/bowtie/bwt_SiC_TM_3_2.fastq.gz  /prj/2024_08_Human/star/str_SiC_TM_3  /biodb/genomes/homo_sapiens/GRCh38_107/GRCh38.107.chr_patch_hapl_scaff.gtf))
	echo " jid_Atria_1st -> $(jidSiCTM3st)"

/prj/2024_08_Human/star/str_SiC_TM_4Aligned.sortedByCoord.out.bam:
	$(eval jidSiCTM4st:=$(jidSiCTM4st)$(shell sbatch --parsable --dependency=afterok:$(jidSiCTM4bw)  /home/bjamshidkia/repositories/pipeline/human/2024_08_Human/s_star_pair.sh  /biodb/genomes/homo_sapiens/GRCh38_107/star \
	/prj/2024_08_Human/bowtie/bwt_SiC_TM_4_1.fastq.gz  /prj/2024_08_Human/bowtie/bwt_SiC_TM_4_2.fastq.gz  /prj/2024_08_Human/star/str_SiC_TM_4  /biodb/genomes/homo_sapiens/GRCh38_107/GRCh38.107.chr_patch_hapl_scaff.gtf))
	echo " jid_Atria_1st -> $(jidSiCTM4st)"







# clean:
# 	rm -r /mnt/c/row_data/fastqc
# 	rm -r /mnt/c/row_data/flexbar 
# 	rm -r /mnt/c/row_data/bowtie 
# 	sudo rm -r /home/STOUT/temp1
# 	sudo rm -r /home/STOUT/temp2
