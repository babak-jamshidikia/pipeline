library(tidyverse)
library(limma)
library(edgeR)
library(org.Mm.eg.db)
library(xlsx)
fc <- readRDS("/home/bjamshidikia/srvrdata/fcdata3.rds")
head(fc$counts)
names(fc)
head(fc$counts)
y <-DGEList(counts =  fc$counts,genes = fc$annotation )
dim(y)
#y<-y[,c(1,2,3)]
dim(y)


o <-order(rowSums(y$counts),decreasing = TRUE)
head(o)
y <-y[o,]
dim(y)
d<-duplicated(y$genes$GeneID)
dim(d)
y <-y[!d,]

nrow(y)
#y$counts <- y$counts[,c(1,2,3)]
y$samples$lib.size <- colSums(y$counts)
y<-normLibSizes(y)
y$samples
plotMDS(y)
sname <- c("Atria","Atria","Atria","Vent","Vent","Vent")
snum <-c(1,2,3,1,2,3)
data.frame(Sample = colnames(y),sname,snum)
design <- model.matrix(~sname+snum)
rownames(design) <-colnames(y)
design

y<-estimateDisp(y,design,robust = TRUE)
y$common.dispersion
plotBCV(y)
fit <-glmFit(y,design)
lrt <-glmLRT(fit)

tt <- topTags(lrt,n=nrow(lrt), p.value=0.05)
colnames(design)
o <-order(lrt$table$PValue)
cpm(y)[o[1:10],]
summary(decideTests(lrt))
plotMD(lrt)
abline(h=c(-2,2),col = "blue")
genename <- mapIds(org.Mm.eg.db, keys=rownames(lrt), keytype="ENSEMBL",column="SYMBOL", multiVals = "first")
lrt$genes$GENENAME <-genename
#write.csv(lrt,file="~/Desktop/BIOINFORMATICS_jacobi/LRT-1.csv")
#writxlsx(tt,"~/Desktop/BIOINFORMATICS_jacobi/tt-6.xlsx",asTable = FALSE, overwrite = TRUE)

wb <- createWorkbook()
s <- createSheet(wb,'result')
addDataFrame(tt, s, row.names=FALSE, startRow=3)

saveWorkbook(wb,'~/Desktop/BIOINFORMATICS_jacobi/tt-7.xlsx')

