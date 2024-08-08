library(tidyverse)
library(limma)
library(edgeR)
library(org.Mm.eg.db)
library(xlsx)
library(biomaRt)
library(dplyr)
fc <- readRDS("/home/bjamshidikia/srvrdata/fcdata.rds")
head(fc$counts)
names(fc)
head(fc$counts)
y <-DGEList(counts =  fc$counts,genes = fc$annotation )
#idfound <-y$genes$GeneID %in% mappedRkeys(org.Mm.eg)
#y<-y[idfound,]
#dim(y$)
#egREFSEQ <- toTable(org.Mm.egREFSEQ)
#head(egREFSEQ)
#m<-  match(y$genes$GeneID,egREFSEQ$accession)
#head(m)
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

data.frame(Sample = colnames(y),sname,snum)
design <- model.matrix(~sname)
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
ensemblM <- useEnsembl( biomart = "genes",dataset= "mmusculus_gene_ensembl")

#attribM <- c("external_gene_name","description","uniprot_gn_id","go_id","ensembl_gene_id")
attribM <- c("external_gene_name","description","ensembl_gene_id")


filtersM <- "ensembl_gene_id"
valueM <-list(tt$table$GeneID)

dataM <-getBM(attributes = attribM ,filters = filtersM , values = valueM ,mart =  ensemblM)

DdataM <- duplicated(dataM$ensembl_gene_id)
dataM <-dataM[!DdataM,]

jointD <- inner_join(tt$table,dataM,by = c("GeneID" = "ensembl_gene_id"))



#dataM2 <-dataM2[!dataM2$ensembl_gene_id == "null",]




#creating EXCEL FILE
wb <- createWorkbook()
s <- createSheet(wb,'result')
addDataFrame(jointD, s, row.names=FALSE, startRow=1)

saveWorkbook(wb,'~/Desktop/BIOINFORMATICS_jacobi/jointD2.xlsx')

