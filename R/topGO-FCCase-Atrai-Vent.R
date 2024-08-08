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
fc <- readRDS("/home/bjamshidikia/srvrdata/fcdata3.rds")
head(fc$counts)
head(fc)
y <- DGEList(counts = fc$counts,genes = fc$annotation,group = c("atria","atria","atria","vent","vent","vent"))
RPKM <-rpkm(y)
RPKM1 <- as.data.frame(RPKM)
#RPKM3 <- RPKM1 %>% 
 # filter(if_any(everything(), ~ .x>0))
RPKM4 <- filter_all(RPKM1,any_vars(. > 0))


#RPKM2<-filter(RPKM1,RPKM1$str_WT_Atria_1Aligned.sortedByCoord.out.bam>0 |
 #             RPKM1$str_WT_Atria_2Aligned.sortedByCoord.out.bam>0|
  #            RPKM1$str_WT_Atria_3Aligned.sortedByCoord.out.bam>0|
   #           RPKM1$str_WT_Vent_1Aligned.sortedByCoord.out.bam>0|
    #          RPKM1$str_WT_Vent_2Aligned.sortedByCoord.out.bam>0|
     #         RPKM1$str_WT_Vent_3Aligned.sortedByCoord.out.bam>0
#)

RPKM4

keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y$samples
#y1 <- y[,c(1,2,3)]
#y2 <-y[,c(4,5,6)]
#head(y1)
#head(y2)
head(y)
# for Atria ->y1 
#oy <-order(rowSums(y$counts),decreasing = TRUE)
#head(oy)
#y <-y[oy,]
#dim(y)
#d<-duplicated(y$genes$GeneID)
#dim(d)
#y <-y[!d,]
#nrow(y)


y$samples$lib.size <- colSums(y$counts)
y<-normLibSizes(y)
y$samples
plotMDS(y)

sname <- c("Atria","Atria","Atria","Vent","Vent","Vent")
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

plogfc <- filter(tt$table, logFC>0)
dim(plogfc)
nlogfc <- filter(tt$table,logFC<0)
dim(nlogfc)
# for positive logFC
#allGenesy1 <- jointy1$PValue
#names(allGenesy1) <- jointy1$external_gene_name
#VallGenesy1 <- as.vector(allGenesy1)




#gnames <- names(geneID2GO)
#geneList <- factor(as.integer(gnames %in% y$genes$GeneID))
#names(geneList) <-  gnames #jointy1$external_gene_name
`#genesely1 <-plogfc$PValue
#names(genesely1) <- plogfc$external_gene_name
#geneList1 <- c(1:10)
geneList1<- rep(0, times=nrow(RPKM4))
geneList1

str(geneList1)
names(geneList1) <- row.names(RPKM4) 
geneList1
geneList2 <- replace(geneList1,names(geneList1) %in% tt$table$GeneID, tt$table$logFC)
geneList2
#github code



if (!file.exists("~/Desktop/BIOINFORMATICS_jacobi/go_mapping_tmp")){
  ensembl <- useEnsembl( biomart = "genes",dataset= "mmusculus_gene_ensembl")
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



go_up_bp <- new("topGOdata",
                nodeSize = 10,
                ontology = "BP",
                geneSel = topDiffGenesUp,
                allGenes = geneList2,
                annot = annFUN.gene2GO,
                gene2GO = geneID2GO
)
go_up_bp
test.stat <- new("elimCount", testStatistic = GOFisherTest, cutOff = 0.05)
resultElimUp <- getSigGroups(go_up_bp,test.stat)
allResUp <-GenTable(go_up_bp,Elim = resultElimUp,
                  orderBy= "Term",topNodes = 20,numChar = 100)
allResUp             

go_down_bp <- new("topGOdata",
                nodeSize = 10,
                ontology = "BP",
                geneSel = topDiffGenesdown,
                allGenes = geneList2,
                annot = annFUN.gene2GO,
                gene2GO = geneID2GO
)
go_down_bp
test.stat <- new("elimCount", testStatistic = GOFisherTest, cutOff = 0.05)
resultElimdown <- getSigGroups(go_down_bp,test.stat)
allResdown <-GenTable(go_down_bp,Elim = resultElimdown,
                    orderBy= "Term",ranksOf = "Elim",topNodes = 20,numChar=100)
allResdown

go_bp <- new("topGOdata",
                  nodeSize = 10,
                  ontology = "BP",
                  geneSel = topDiffGenes,
                  allGenes = geneList2,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO
)
go_bp
test.stat <- new("elimCount", testStatistic = GOFisherTest, cutOff = 0.05)
resultElim <- getSigGroups(go_bp,test.stat)

allResboth <-GenTable(go_bp,Elim = resultElim,
                      orderBy= "Term",ranksOf = "Elim",topNodes = 20,numChar = 100)
allResboth
# creating the plotwith ggplot

basic_plot_Up <- ggplot(allResUp, aes(x = reorder(Term,Elim), y = -log2(as.numeric(Elim)))) +
  geom_bar(stat = "identity", fill = "skyblue") + coord_flip() +
  labs(title = "Result UP")
basic_plot_Up

basic_plot_down <- ggplot(allResdown, aes(x = reorder(Term,Elim), y = -log2(as.numeric(Elim)))) +
  geom_bar(stat = "identity", fill = "red") + coord_flip() +
  labs(title = "Result Down")
basic_plot_down

basic_plot_Up | basic_plot_down

#creating EXCEL FILE

wb <- createWorkbook()

addWorksheet(wb, sheetName = "allresUP");
writeDataTable(wb, sheet = 1, x = allResUp, rowNames = TRUE);

addWorksheet(wb, sheetName = "allresdwn");
writeDataTable(wb, sheet = 2, x = allResdown, rowNames = TRUE);

addWorksheet(wb, sheetName = "allresBoth");
writeDataTable(wb, sheet = 3, x = allResboth, rowNames = TRUE);

#addWorksheet(wb, sheetName = "RPKM4");
#writeDataTable(wb, sheet = 4, x = RPKM4, rowNames = TRUE);

saveWorkbook(wb, "~/Desktop/BIOINFORMATICS_jacobi/ALLRES.xlsx", overwrite = TRUE)

