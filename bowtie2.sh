#! /bin/bash
# bowtie2 -x /mnt/c/Applications/data/index/ -1 /mnt/c/Applications/data/flexbar/35_1.fastq.gz -2/mnt/c/Applications/data/flexbar/35_2.fastq.gz  --no-unal --omit-sec-seq --threads 4  --time --un-conc-gz /mnt/c/Applications/data/bwt/35.fastq.gz 2> /mnt/c/Applications/data/bwt/35.log
  bowtie2 -x /mnt/c/Applications/bwt/index/35 -1 /mnt/c/Applications/data/flexbar/35_1.fastq.gz -2 /mnt/c/applications/data/flexbar/35_2.fastq.gz --no-unal --omit-sec-seq --threads 4 --mm --seed 1337 --time --un-conc-gz /mnt/c/applications/bwt/35.fastq.gz 2> /mnt/c/applications/bwt/35.log

