#! /bin/bash

STAR --genomeDir $1 --runThreadN 8  --readFilesIn $2 $3 --sjdbGTFfile $6 --readFilesCommand zcat --outFileNamePrefix $4/$target/ --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --quantMode GeneCounts --outReadsUnmapped Fastx --outWigType bedGraph --genomeLoad NoSharedMemory  