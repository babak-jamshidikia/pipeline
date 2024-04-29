#! /bin/bash

#flexbar -r $1 -p $2 -t $3/$target  -n 15 -z GZ -m 30 -u 0  -q TAIL -qt 28 -as "AGATCGGAAGAG" -qf sanger -j

#flexbar -r /mnt/c/Applications/data/SRR27959136_1.fastq.gz -p  /mnt/c/Applications/data/SRR27959136_2.fastq.gz  -t  /mnt/c/Applications/data/flexbar/36/36  -n 10 -z GZ -m 30 -u 0  -q TAIL -qt 28 -aa TruSeq -ap ON  -qf sanger -j
#flexbar -r /mnt/c/Applications/data/SRR27959136_1.fastq.gz -p  /mnt/c/Applications/data/SRR27959136_2.fastq.gz  -t  /mnt/c/Applications/data/flexbar/36/36\  -n 4 -z GZ -m 151 -aa TruSeq -ap ON  -qf sanger 
flexbar -r /mnt/c/Applications/data/SRR27959135_1.fastq.gz -p  /mnt/c/Applications/data/SRR27959135_2.fastq.gz  -t  /mnt/c/Applications/data/flexbar/35   -n 15 -z GZ -m 30 -u 0  -q TAIL -qt 28 -as "AGATCGGAAGAG" -qf sanger -j
