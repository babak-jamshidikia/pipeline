library(limma)
library(edgeR)
library(Rsubread)
library(org.Mm.eg.db)

bam.files <- list.files(path = "/prj/2024_05_Heart/star",pattern=".bam$",full.name = TRUE)
bam.files
annotfile <-"/biodb/genomes/mus_musculus/GRCm39_107/star/GRCm39.107.gtf"
print(annotfile)

fc <-featureCounts(files = bam.files,annot.ext = annotfile,isGTFAnnotationFile=TRUE,GTF.featureType='gene',GTF.attrType = "gene_id",nthreads=5 )
names(fc)

saveRDS(fc,"fcdata3.rds")

rr <-"End scriptR"
print(rr)
