library(tidyverse)
library(limma)
library(edgeR)
#library(xlsx)
library(dplyr)
library(biomaRt)
library(topGO)
library(openxlsx)
#library(ggplot2)
#library(ggrepel)
library(dendextend)

fc <- readRDS("/home/bjamshidikia/srvrdata/fcdata3.rds")
head(fc$counts)
head(fc)
y <- DGEList(counts = fc$counts,genes = fc$annotation,group = c("atria","atria","atria","vent","vent","vent"))
RPKM <-rpkm(y)
RPKM1 <- as.data.frame(RPKM)
#RPKM3 <- RPKM1 %>% 
 # filter(if_any(everything(), ~ .x>0))
RPKM4 <- filter_all(RPKM1,any_vars(. > 0))



RPKM4
#RPKM5 <- head(RPKM4,1000)
#RPKM5 <- RPKM4
#RPKM6 <- log10(RPKM5)

#colnames(RPKM6) <-c("Atria_1","Atria_2","Atria_3","Vent_1","Vent_2","Vent_3") #substr(colnames(RPKM6),8,14)

#matRPKM <- as.matrix(RPKM6)
#my_heatmap <- pheatmap(matRPKM,
#                       cluster_rows = FALSE,
#                       color = viridis(n = 2000,option = "magma"),
#                       legend_breaks = c(-2,0,2),
#                       show_rownames = FALSE,
#                       show_colnames = TRUE,
#                       border_color = NA,
#                       main = "Atria vs Vent  \n Heatmap based on Log10(RPKM) values ",
#                       scale = "row",
                       
#                      )
#my_heatmap


#class(my_heatmap)
#names(my_heatmap)
#my_heatmap$tree_row %>%
#  as.dendrogram() %>%
#  plot(horiz = TRUE)
  


 #heatmap(matRPKM)


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
#tt2 <- topTags(lrt,n=nrow(lrt),p.value = 1)
# for tt
#ensembly1 <- useEnsembl( biomart = "genes",dataset= "mmusculus_gene_ensembl")
#attriby1 <- c("external_gene_name","description","ensembl_gene_id")


#filtersy1 <- "ensembl_gene_id"
#valuey1 <-list(tt$table$GeneID)

#datay1 <-getBM(attributes = attriby1 ,filters = filtersy1 , values = valuey1 ,mart =  ensembly1)

#Ddatay1 <- duplicated(datay1$ensembl_gene_id)
#datay1 <-datay1[!Ddatay1,]

#tt <- left_join(tt$table,datay1,by = c("GeneID" = "ensembl_gene_id"))
#tt2 <- left_join(tt2$table,datay1,by = c("GeneID" = "ensembl_gene_id"))



colnames(design)
head(tt)


#dfTT <- as.data.frame(tt)
#dfTT$FDR <- -log10(dfTT$FDR)
#dfTT
#dfTT$diffexpressed <- "NO"
#dfTT
#dfTT$diffexpressed[dfTT$logFC > 0 & dfTT$FDR > -log10(0.05)] <- "UP"
#dfTT$diffexpressed[dfTT$logFC < 0 & dfTT$FDR > -log10(0.05)] <- "DOWN"


#dfTT$delabel <- ifelse(dfTT$external_gene_name %in% head(dfTT[order(dfTT$logFC), "external_gene_name"], 10), dfTT$external_gene_name , NA)

#p1 <- ggplot(dfTT,aes(x= logFC,y= FDR,col = diffexpressed,label = delabel)) +
#  geom_vline(xintercept = 0,col = "green",linetype = "dashed")+
#  geom_hline(yintercept = -log10(0.05),col= "green",linetype = "dashed")+
#  geom_point(size = 2)+
#  scale_color_manual(values = c("#00AFBB", "grey", "red"), 
#                     labels = c("Downregulated", "Not significant", "Upregulated"))+
#  coord_cartesian(ylim = c(0, 350), xlim = c(-15, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
#  labs(color = 'Severe', #legend_title, 
#       x = expression("logFC"), y = expression("-log"[10]*"(FDR)")) + 
#  scale_x_continuous(breaks = seq(-15, 10, 5))+
#  scale_y_continuous(breaks = seq(0,350,50))+
#  ggtitle("Volcano Plot") +
#  geom_text_repel(max.overlaps = Inf)
  


#p1
#dfTT
# for tt2

#dfTT2 <- as.data.frame(tt2)
#dfTT2$LOGFDR <- -log10(dfTT2$FDR)
#dfTT2
#dfTT2$diffexpressed <- "NO"
#dfTT2
#dfTT2$diffexpressed[dfTT2$logFC > 0 & dfTT2$FDR < (0.05)] <- "UP"
#dfTT2$diffexpressed[dfTT2$logFC < 0 & dfTT2$FDR < (0.05)] <- "DOWN"

#nUpdfTT2 <- length(dfTT2$diffexpressed[dfTT2$logFC > 0 & dfTT2$FDR < (0.05)])
#nDwdfTT2 <- length(dfTT2$diffexpressed[dfTT2$logFC < 0 & dfTT2$FDR < (0.05)])




#if(length(head(dfTT2[order(dfTT2$logFC >0), "external_gene_name"], 20))==20){
#dfTT2$delabel <- ifelse(dfTT2$external_gene_name %in% head(dfTT2[order(dfTT2$logFC,decreasing =  TRUE), "external_gene_name"], 10), dfTT2$external_gene_name , NA)
#}

#dfTT2$delabel <- ifelse(dfTT2$external_gene_name %in% head(dfTT2[order(dfTT2$logFC,decreasing =  TRUE), "external_gene_name"], 10), dfTT2$external_gene_name , NA)


#if(length(head(dfTT2[order(dfTT2$logFC < 0), "external_gene_name"], 10))==10){
#dfTT2$delabel <- ifelse(dfTT2$external_gene_name %in% head(dfTT2[order(dfTT2$logFC ,decreasing = FALSE), "external_gene_name"], 13), dfTT2$external_gene_name ,dfTT2$delabel)
#}




#p2 <- ggplot(dfTT2,aes(x= logFC,y= LOGFDR,col = diffexpressed,label = delabel)) +
#  geom_vline(xintercept = 0,col = "black",linetype = "dashed")+
#  geom_hline(yintercept = -log10(0.05),col= "black",linetype = "dashed")+
#  geom_point(size = 2,alpha = 0.2)+
#  scale_color_manual(values = c("#00AFBB", "black", "red"), 
#                     labels = c("Downregulated", "Not significant", "Upregulated"))+
#  coord_cartesian(ylim = c(-20, 350), xlim = c(-15, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
#  labs(color = 'Severe', #legend_title, 
#       x = expression("log"[2]* "(Fold Change)"),
#       y = expression("-log"[10]*"(FDR)")) + 
  #scale_x_continuous(breaks = seq(-15, 10, 5))+
  #scale_y_continuous(breaks = seq(0,5,1))+
#  ggtitle("Volcano Plot") +
  #geom_text_repel(color = "black", max.overlaps = Inf,xlim = c(-Inf, Inf), ylim = c(-Inf, Inf))+
#  geom_label_repel(fill = "white",max.overlaps = Inf, xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),color = "black",box.padding = 0.25,alpha = 0.7)
  

#p2
#dfTT2

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

funcGO <-function(ont,TopN,P_valcutoff,updown){

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
allResup2$Elim <- -log2(as.numeric(allResup2$Elim))


basic_plot_Up <- ggplot(allResup2, aes(x = reorder(Term,Elim),y= Elim)) +
  geom_bar(stat = "identity", fill = graphfill) + coord_flip() +
  labs(title = paste("Result ",updown,ont,".pdf"),
       x = "Term",
       y = "-log2(Elim)")

basic_plot_Up
#ggsave(
 # filename = paste("Result ",updown,ont,".pdf"),
#  plot = basic_plot_Up,
#  device = "pdf",
#  path = "~/Desktop/BIOINFORMATICS_jacobi/graphs"
#)

} # end of function

pdf()
funcGO("BP",20,0.05,"up")
funcGO("BP",20,0.05,"down")
funcGO("BP",20,0.05,"both")

funcGO("MF",20,0.05,"up")
funcGO("MF",20,0.05,"down")
funcGO("MF",20,0.05,"both")

funcGO("CC",20,0.05,"up")
funcGO("CC",20,0.05,"down")
funcGO("CC",20,0.05,"both")
dev.off()



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

