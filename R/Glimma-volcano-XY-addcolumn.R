library(tidyverse)
library(limma)
library(edgeR)
library(org.Mm.eg.db)
#library(xlsx)
library(dplyr)
#library(package = "mouse4302.db",character.only = TRUE)
library(biomaRt)
library(topGO)
library(openxlsx)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(pheatmap)
library(dendextend)
library(viridis)
library(Glimma)

fc <- readRDS("/home/bjamshidikia/srvrdata/fcdata3.rds")
head(fc$counts)

ensembly2 <- useEnsembl( biomart = "genes",dataset= "mmusculus_gene_ensembl")
attriby2 <- c("external_gene_name","uniprot_gn_symbol","ensembl_gene_id")


filtersy2 <- "ensembl_gene_id"
valuey2 <-list(fc$annotation$GeneID)
datay2 <-getBM(attributes = attriby2 ,filters = filtersy2 , values = valuey2 ,mart =  ensembly2)
datay2 <- datay2[!duplicated(datay2$ensembl_gene_id),]

fc$annotation <- left_join(fc$annotation,datay2,by = c("GeneID" = "ensembl_gene_id"))
dim(fc)
dim(fc$annotation)
dim(fc$counts)


head(fc)
y <- DGEList(counts = fc$counts,genes = fc$annotation,group = c("atria","atria","atria","vent","vent","vent"))

RPKM <-rpkm(y)
RPKM1 <- as.data.frame(RPKM)
RPKM4 <- filter_all(RPKM1,any_vars(. > 0))



RPKM4

RPKM5 <- RPKM4
RPKM6 <- log10(RPKM5)





#heatmap(matRPKM)


keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y$samples
head(y)


y$samples$lib.size <- colSums(y$counts)
y<-normLibSizes(y)
y$samples
plotMDS(y)

sname <- c("Atria","Atria","Atria","Vent","Vent","Vent")

data.frame(Sample = colnames(y),sname)
design <- model.matrix(~sname)
rownames(design) <-colnames(y)
design

y<-estimateDisp(y,design,robust = TRUE)
y$common.dispersion

plotBCV(y)

fit <-glmFit(y,design)
lrt <-glmLRT(fit)



tt <- topTags(lrt,n=nrow(lrt), p.value=0.05)
#tt2 <- topTags(lrt,n=nrow(lrt),p.value = 1)



#lrt$genes <-left_join(lrt$genes,datay2,by = c("GeneID" = "ensembl_gene_id"))
#y$genes <- left_join(y$genes,datay3,by = c("GeneID" = "ensembl_gene_id"))







glimmaVolcano(lrt , dge = y,main = paste("Atria", "vs", "Vent"),status.cols = c("dodgerblue", "gray", "firebrick") )

glimmaVolcano(
  lrt,
  dge = NULL,
  counts = y$counts,
  groups = y$samples$group,
  status = edgeR::decideTestsDGE(lrt),
  anno =lrt$genes ,
  display.columns = NULL,
  status.cols = c("dodgerblue", "gray", "firebrick"),
  sample.cols = NULL,
  p.adj.method = "BH",
  transform.counts = c("logcpm", "cpm", "rpkm", "none"),
  main = paste("Atria", "vs", "Vent"),
  xlab = "logFC",
  ylab = "negLog10PValue",
  html = "~/Desktop/BIOINFORMATICS_jacobi/graphs/glimma/b-1.html",
  width = 1200,
  height = 1200,

)

plothtml <- glimmaMA(lrt, dge=y, status.colours=c("#3977db","#3d3f42","#db0d4e"))
htmlwidgets::saveWidget(
  plothtml,
  "~/Desktop/BIOINFORMATICS_jacobi/graphs/glimma/b-2.html",
  selfcontained = TRUE,
  libdir = NULL,
  background = "white",
  title = "MA plot"  #class(widget)[[1]],
  #knitrOptions = list()
)

#glimmaXY(x=lrt$table$logFC, y=lrt$table$logCPM , dge=y,status.cols = c("dodgerblue", "silver", "firebrick"))




