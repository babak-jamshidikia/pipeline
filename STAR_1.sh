#! /bin/bash

STAR --genomeDir /mnt/c/Applications/data/genome/ref --runThreadN 4 --readFilesIn /mnt/c/Applications/bwt/35.fastq.1.gz /mnt/c/Applications/bwt/35.fastq.2.gz --sjdbGTFfile /mnt/c/Applications/data/genome/genomic.gtf --readFilesCommand zcat --outFileNamePrefix /home/STAROUT/output --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --quantMode GeneCounts --outReadsUnmapped Fastx --outWigType bedGraph --genomeLoad NoSharedMemory --outTmpDir /home/STAROUT/tmpout
