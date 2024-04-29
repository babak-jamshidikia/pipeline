

all: fastqc flexbar bowtie STAR
#all: STAR

#Rules for fastqc
fastqc: makingFastqcdirectory /mnt/c/row_data/fastqc/sub35_1_fastqc.html /mnt/c/row_data/fastqc/sub35_2_fastqc.html
makingFastqcdirectory:
	echo "making fastqc directory"
	mkdir /mnt/c/row_data/fastqc
#/mnt/c/row_data/fastqc: /mnt/c/row_data/fastqc/sub35_1_fastqc.html /mnt/c/row_data/fastqc/sub35_2_fastqc.html
/mnt/c/row_data/fastqc/sub35_1_fastqc.html:
	echo "making fastqc for the row file output is in /mnt/c/row_data/fastqc/sub35_1_fastqc.html  "
	fastqc -t 2 /mnt/c/row_data/sub35_1.fastq -o /mnt/c/row_data/fastqc 
/mnt/c/row_data/fastqc/sub35_2_fastqc.html:
	echo "making fastqc for the row file output is in /mnt/c/row_data/fastqc/sub35_2_fastqc.html  "
	fastqc -t 2 /mnt/c/row_data/sub35_2.fastq -o /mnt/c/row_data/fastqc 

#Rules for flexbar

flexbar: making_flexbar_directory /mnt/c/row_data/flexbar/flex_35_1.fastq.gz /mnt/c/row_data/flexbar/flex_35_2.fastq.gz
making_flexbar_directory: 
	
	tput setaf 206;echo "making directory for flexbar"
	tput setaf 255;echo "     "
	
	mkdir /mnt/c/row_data/flexbar 
/mnt/c/row_data/flexbar/flex_35_1.fastq.gz:
	
	echo "making flexbar for the row file output is  /mnt/c/row_data/flexbar/flex_35_1.fastq   "
	
	flexbar -r /mnt/c/row_data/sub35_1.fastq  -t  /mnt/c/row_data/flexbar/flex_35_1   -n 15 -z GZ -m 30 -u 0  -q TAIL -qt 28 -as "AGATCGGAAGAG" -qf sanger -j
	
/mnt/c/row_data/flexbar/flex_35_2.fastq.gz:
	
	echo "making flexbar for the row file output is  /mnt/c/row_data/flexbar/flex_35_2.fastq   "
	
	flexbar -r /mnt/c/row_data/sub35_2.fastq   -t  /mnt/c/row_data/flexbar/flex_35_2   -n 15 -z GZ -m 30 -u 0  -q TAIL -qt 28 -as "AGATCGGAAGAG" -qf sanger -j




#Rules for bowtie 
bowtie: makingbowtiefolder /mnt/c/row_data/bowtie/bwt_35_1.fastq.gz /mnt/c/row_data/bowtie/bwt_35_2.fastq.gz 
makingbowtiefolder:
	mkdir /mnt/c/row_data/bowtie
/mnt/c/row_data/bowtie/bwt_35_1.fastq.gz:
	echo "making bowtie for the row file output is in /mnt/c/row_data/bowtie /bwt_35_1.fastq   "
	bowtie2 -x /mnt/c/row_data/bwtindex/35 -U /mnt/c/row_data/flexbar/flex_35_1.fastq.gz   -S /dev/null --no-unal --omit-sec-seq --threads 20 --mm --seed 1337 --time --un-gz /mnt/c/row_data/bowtie/bwt_35_1.fastq.gz 2> /mnt/c/row_data/bowtie/bwt_35_1.log
/mnt/c/row_data/bowtie/bwt_35_2.fastq.gz:
	echo "making bowtie for the row file output is in /mnt/c/row_data/bowtie/bwt_35_2.fastq   "
	bowtie2 -x /mnt/c/row_data/bwtindex/35 -U /mnt/c/row_data/flexbar/flex_35_2.fastq.gz   -S /dev/null --no-unal --omit-sec-seq --threads 20 --mm --seed 1337 --time --un-gz /mnt/c/row_data/bowtie/bwt_35_2.fastq.gz 2> /mnt/c/row_data/bowtie/bwt_35_2.log




#Rules for STAR
STAR: makingSTARfolder /mnt/c/row_data/STAR/ST_35_1.bam /mnt/c/row_data/STAR/ST_35_2.bam
makingSTARfolder:
	mkdir /mnt/c/row_data/STAR
/mnt/c/row_data/STAR/ST_35_1.bam:
	echo "making BAM File for the first file"
	sudo STAR --genomeDir /mnt/c/row_data/stindex/ --runThreadN 4 --readFilesIn /mnt/c/row_data/bowtie/bwt_35_1.fastq.gz  --sjdbGTFfile /mnt/c/row_data/stindex/genomic.gtf --readFilesCommand zcat --outFileNamePrefix /mnt/c/row_data/STAR/ST_35_1 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --quantMode GeneCounts --outReadsUnmapped Fastx --outWigType bedGraph --genomeLoad NoSharedMemory --outTmpDir /home/STOUT/temp1
/mnt/c/row_data/STAR/ST_35_2.bam:
	echo "making BAM File for the second file"
	sudo STAR --genomeDir /mnt/c/row_data/stindex/ --runThreadN 4 --readFilesIn /mnt/c/row_data/bowtie/bwt_35_2.fastq.gz --sjdbGTFfile /mnt/c/row_data/stindex/genomic.gtf --readFilesCommand zcat --outFileNamePrefix /mnt/c/row_data/STAR/ST_35_2 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --quantMode GeneCounts --outReadsUnmapped Fastx --outWigType bedGraph --genomeLoad NoSharedMemory --outTmpDir /home/STOUT/temp2

#TRACING  


# clean:
# 	rm -r /mnt/c/row_data/fastqc
# 	rm -r /mnt/c/row_data/flexbar 
# 	rm -r /mnt/c/row_data/bowtie 
# 	sudo rm -r /home/STOUT/temp1
# 	sudo rm -r /home/STOUT/temp2