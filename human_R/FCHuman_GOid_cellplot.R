library(tidyverse)
library(limma)
library(edgeR)
library(dplyr)
library(biomaRt)
library(topGO)
library(openxlsx)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(dendextend)
library(viridis)
library(CellPlot)
library(data.table)
library(parallel)
library(stringr)
# Bioconductor
library(BiocParallel)
library(DESeq2)
library(annotate)



fc <- readRDS("/home/bjamshidikia/srvrdata/fcdataHuman.rds")
head(fc$counts)
head(fc)
y <- DGEList(counts = fc$counts,genes = fc$annotation,group = c("SiC","SiC","SiC","SiC","SiC_TM","SiC_TM","SiC_TM","SiC_TM"))
RPKM <-rpkm(y)
RPKM1 <- as.data.frame(RPKM)
RPKM4 <- filter_all(RPKM1,any_vars(. > 0))



RPKM4

# removing duplicate data 


# ------------------




keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y$samples
head(y)


y$samples$lib.size <- colSums(y$counts)
y<-normLibSizes(y)
y$samples
plotMDS(y)

sname <- c("SiC","SiC","SiC","SiC","SiC_TM","SiC_TM","SiC_TM","SiC_TM")
#snum <-c(1,2,3,1,2,3)
data.frame(Sample = colnames(y),sname)
design <- model.matrix(~sname)
rownames(design) <-colnames(y)
design

y<-estimateDisp(y,design,robust = TRUE)
y$common.dispersion

plotBCV(y)

fit <-glmFit(y,design)
fit1 <- glmQLFit(y,design,robust = TRUE)
plotQLDisp(fit1)
lrt <-glmLRT(fit)

tt <- topTags(lrt,n=nrow(lrt), p.value=0.05)
#ensembly1 <- useEnsembl( biomart = "genes",dataset= "mmusculus_gene_ensembl")
#attriby1 <- c("external_gene_name","description","ensembl_gene_id")


#filtersy1 <- "ensembl_gene_id"
#valuey1 <-list(tt$table$GeneID)

#datay1 <-getBM(attributes = attribM ,filters = filtersM , values = valueM ,mart =  ensembly1)

#Ddatay1 <- duplicated(datay1$ensembl_gene_id)
#datay1 <-datay1[!Ddatay1,]

#jointy1 <- left_join(tt$table,datay1,by = c("GeneID" = "ensembl_gene_id"))




colnames(design)
head(tt)

geneList1<- rep(0, times=nrow(RPKM4))
geneList1

str(geneList1)
names(geneList1) <- row.names(RPKM4) 
geneList1
geneList2 <- replace(geneList1,names(geneList1) %in% tt$table$GeneID, tt$table$logFC)
geneList2
#github code



if (!file.exists("~/Desktop/BIOINFORMATICS_jacobi/go_mapping_tmp")){
  ensembl <- useEnsembl( biomart = "genes",dataset= "hsapiens_gene_ensembl")
  go_annotations <- getBM(attributes = c("ensembl_gene_id", "go_id") , values =row.names(RPKM4), mart = ensembl)
  go_annotations
  dim(go_annotations)
  go_annotations <- go_annotations %>%
    group_by(ensembl_gene_id) %>%
    mutate(go_id = paste0(go_id, collapse = ", ")) %>%
    distinct()
  
  write.table(go_annotations, file = "~/Desktop/BIOINFORMATICS_jacobi/go_mapping_tmp", row.names = F, na = "", col.names = F, sep = "\t", quote = FALSE)
}

# read GO gene<>GO mappings from file
geneID2GO <- readMappings(file = "~/Desktop/BIOINFORMATICS_jacobi/go_mapping_tmp")

topDiffGenesUp <- function(allScore) {
  return(allScore > 0)
}

topDiffGenesdown <- function(allScore) {
  return(allScore < 0)
}

topDiffGenes <- function(allScore) {
  return(allScore)
}

#funcGO <-function(ont,TopN,P_valcutoff,updown){
ont = "BP"
TopN = 20
P_valcutoff = 0.05
updown = "both"

 
  
  if(updown=="up" ){gTitle <-"Result UP"
  graphfill<-"green"
  Gsel <- topDiffGenesUp
  }else if(updown=="down" ){
    gTitle <-"Result Down"
   graphfill<-"red"
    Gsel <- topDiffGenesdown
  }else if(updown=="both" ){
    gTitle <-"Result both"       
    graphfill<-"blue"
    Gsel <- topDiffGenes
} 
 
  
  go_up_bp <- new("topGOdata",
                  nodeSize = 10,
                  ontology = ont,
                  geneSel = Gsel,
                  allGenes = geneList2,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO
  )
  go_up_bp
  test.stat <- new("elimCount", testStatistic = GOFisherTest, cutOff = P_valcutoff)
  resultElimUp <- getSigGroups(go_up_bp,test.stat)
  allResUp <-GenTable(go_up_bp,Elim = resultElimUp,
                      orderBy= "Elim",topNodes = TopN,numChar = 100)
  allResUp             
  allResup2 <- allResUp
  #allResup2$Elim <- -log2(as.numeric(allResup2$Elim))
  allResup2$LogEnrich <-  log2(allResup2$Significant/allResup2$Expected)
  #tt$table$logFC <- -log2(tt$table$logFC)
  colnames(allResup2) <- c("GO.ID","Term","Annotated","Significant","Expected","pvalCutOff","LogEnrich")
  
  
  
  
  
  ga <- genesInTerm(go_up_bp) # GenesAnnotated | list of genes per go-terms
  ga <- ga[allResup2$GO.ID] # eliminate missing terms
  ga2 <- ga
  names(ga) <- NULL
  allResup2$GenesAnnotated <- ga

  xs <- TTgoid [,c("FDR", "logFC")] # significant stats subset
  xs <- subset(xs, FDR < 0.05)
  allResup2$GenesSignificant <- lapply(allResup2$GenesAnnotated, intersect, rownames(xs)) # extract genes
  
  gadf <- ldply (ga2, data.frame)
  colnames(gadf) <- c("GO.ID","GeneID")
  ttdf <- as.data.frame(tt)
  TTgoid <- inner_join(gadf,ttdf,c("GeneID" = "GeneID" ))
  
  goidlist <- allResup2$GO.ID
  
  # TTgoidnd <- TTgoid[!d,]
  
  
  lsfc <- list()
  lspadj <- list()
  i <- 1
  for(Gid in goidlist) {
    
    print(Gid)
    Fgoid <-  filter(TTgoid,TTgoid$GO.ID == Gid)
    print(Fgoid)
    lsfc[i] <- list(Fgoid$logFC)
    lspadj[i] <- list(Fgoid$FDR)
  #  lsfc <- list(lsfc)
    #print(lsfc)
    #print(lspadj)
  #  ind <- match(Gid,allResup2$GO.ID)
  #  allResup2[ind,"log2Foldchange"] <- lsfc
    
    #rows_update(allResup2, tibble( "log2Foldchange"= lsfc, "GO.ID" = Gid), by = "GO.ID")
    i <- i +1
    }
  #print(i)
  
  allResup2$log2Foldchange <- lsfc
  allResup2$padj <- lspadj
  # 
  
  #allResup2$log2Foldchange <- -log2(allResup2$log2Foldchange)
  #return(TTtoAllres)
  
  #TTtoAllres %>% relocate(go_id,Term,Annotated,Significant,Expected, .before = logFC)
  #TTtoAllres <- TTtoAllres[-c(3)]
  
  #ei.rows <- mclapply(allResup2$GenesSignificant, function (y) {
  #  if (length(y)) as.list(xs[y,,drop=FALSE])
  #  else as.list(rep(NA_real_, length(xs)))
  #}, mc.cores = 10)
  #ei <- mclapply(names(xs), function(z) {
  #  lapply(ei.rows, "[[", z)
  #}, mc.cores = 10)
  #ei <- structure(ei, names = names(xs), row.names = seq(nrow(r)), class = "data.frame")
  #row.names(ei) <- NULL
  #allResup2 <- data.frame(allResup2, ei, stringsAsFactors = FALSE, check.names = FALSE)
  
  
  
  
  
  x <- subset(allResup2,pvalCutOff <= 0.05 & Significant > 20)
  x <- x[order(-x$LogEnrich),]
  pdf(file = "/home/bjamshidikia/Desktop/BIOINFORMATICS_jacobi/graphs/cellplot/p2.pdf")
  cell.plot(x = setNames(x$LogEnrich,x$Term),
            cells = x$log2Foldchange ,
            main = " GO enrichment (Not vs CLL)  ",
            x.mar = c(0.4,0),
            key.n = 7,
            y.mar = c(0.1,0),
            cex = 1.6,
            cell.outer = 3,
            bar.scale = 0.7,
            space = 0.2,
  )
  
  dev.off()

  
  pdf(file = "/home/bjamshidikia/Desktop/BIOINFORMATICS_jacobi/graphs/cellplot/S3.pdf")
  sym.plot( x= setNames(x$LogEnrich,x$Term),
            cells = x$log2Foldchange,
            x.annotated = x$Annotated,
            main = "GO enrichment (NoT vs CLL)",
            x.mar = c(0.47,0),
            key.n = 7,
            cex = 1.6,
            axis.cex = 0.8,
            group.cex = 0.7
  )
  dev.off()  
  
  
  




