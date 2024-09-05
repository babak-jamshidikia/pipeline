library(limma)
library(edgeR)
library(Rsubread)
library(org.Mm.eg.db)

bam.files <- list.files(path = "/prj/2024_08_Human/star",pattern=".bam$",full.name = TRUE)
bam.files
annotfile <-"/biodb/genomes/homo_sapiens/GRCh38_107/GRCh38.107.gtf"
print(annotfile)

fc_H <-featureCounts(files = bam.files,annot.ext = annotfile,isGTFAnnotationFile=TRUE,GTF.featureType='gene',GTF.attrType = "gene_id",nthreads=5 )
names(fc_H)

saveRDS(fc_H,"fcdataHuman.rds")

rr <-"End scriptR"
print(rr)
