library(tidyverse)
library(limma)
library(edgeR)
library(org.Mm.eg.db)
library(dplyr)
library(biomaRt)
library(topGO)
library(openxlsx)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(dendextend)
library(viridis)

fc <- readRDS("/home/bjamshidikia/srvrdata/fcdataHuman.rds")
head(fc$counts)
head(fc)
y <- DGEList(counts = fc$counts,genes = fc$annotation,group = c("SiC_1","SiC_2","SiC_3","SiC_4","SiC_TM_1","SiC_TM_2","SiC_TM_3","SiC_TM_4"))
RPKM <-rpkm(y)
RPKM1 <- as.data.frame(RPKM)
RPKM4 <- filter_all(RPKM1,any_vars(. > 0))



RPKM4
# removing duplicate data 
#row1 <- filter(RPKM4,row.names(RPKM4) == "ENSMUSG00000116933" )
#row2 <- filter(RPKM4,row.names(RPKM4) == "ENSMUSG00000022956" )
#row <- row1+row2


#RPKM4 <- filter(RPKM4,row.names(RPKM4) != "ENSMUSG00000116933")
#RPKM4 <- filter(RPKM4,row.names(RPKM4) != "ENSMUSG00000022956")
# ------------------
RPKM4[nrow(RPKM4) + 1,] <- row
RPKM4



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

funcGO <-function(ont,TopN,P_valcutoff,updown){
#ont = "BP"
#TopN = 20
#P_valcutoff = 0.05
#updown = "up"

  if(!dir.exists(paste0("~/Desktop/BIOINFORMATICS_jacobi/graphs/HeatMaps/",ont))){
  system(paste0("mkdir ~/Desktop/BIOINFORMATICS_jacobi/graphs/HeatMaps/",ont))
  }
  
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
  
  
  #basic_plot_Up <- ggplot(allResup2, aes(x = reorder(Term,Elim),y= Elim)) +
  #  geom_bar(stat = "identity", fill = graphfill) + coord_flip() +
  #  labs(title = paste("Result ",updown,ont,".pdf"),
  #       x = "Term",
  #       y = "-log2(Elim)")
  
#  basic_plot_Up
#------------------------------------  
  
  
#---------------------
  
  #print(updown)
  if(!dir.exists(paste0("~/Desktop/BIOINFORMATICS_jacobi/graphs/HeatMaps/",ont,"/",updown))){
  system(paste0("mkdir ~/Desktop/BIOINFORMATICS_jacobi/graphs/HeatMaps/",ont,"/",updown))
  }

  gendf <- as.data.frame(geneList2)
  
  ensembl <- useEnsembl( biomart = "genes",dataset= "hsapiens_gene_ensembl")
  godata2 <- getBM(attributes = c("ensembl_gene_id","external_gene_name"), filters = "ensembl_gene_id", values =row.names(gendf), mart = ensembl)
  
  wb <- createWorkbook()
  
for (i in 1:nrow(allResUp) ){
#i = 1

#print(allResUp$GO.ID[i])
  gt <- as.data.frame(genesInTerm(go_up_bp,allResUp$GO.ID[i]))
  colnames(gt) <- c("E_ID")
  Termdata <- filter(godata2,ensembl_gene_id %in% gt$E_ID)
  #print(nrow(Termdata))
  TermDataEXL <- inner_join(Termdata,tt$table,c("ensembl_gene_id" = "GeneID" ))
  
  #print(nrow(TermDataEXL))
  #print("--------------")
  
  TermDataEXL <- TermDataEXL[-c(3,4,5,6,7,9,10,11)]
  #writing to worksheet
  
  
  addWorksheet(wb, sheetName = paste0(i,"_", substr(allResUp$Term[i],1,28) ))
  writeDataTable(wb, sheet = i, x = TermDataEXL, rowNames = FALSE)
  #print(i)
  #print(allResUp$GO.ID[i])
  #print(allResUp$Annotated[i])
  #print(allResUp$Significant[i])
  #print("================")
  #-----------------------
  #preparing data for making heatmap
  #RPKMHeat <-filter(RPKM4,row.names(RPKM4) %in% gt$E_ID)
  RPKMHeat <-filter(RPKM4,row.names(RPKM4) %in% TermDataEXL$ensembl_gene_id)
  RPKMHeat <- log10(RPKMHeat)
  colnames(RPKMHeat) <-c("SiC_1","SiC_2","SiC_3","SiC_4","SiC_TM_1","SiC_TM_2","SiC_TM_3","SiC_TM_4")
  RPKMHeatGenename <- RPKMHeat
  RPKMHeatGenename$GID <- row.names(RPKMHeatGenename)
  RPKMHeatGenename <- inner_join(RPKMHeatGenename,Termdata,c("GID" = "ensembl_gene_id"))
  RPKMHeatGenename <- filter(RPKMHeatGenename,RPKMHeatGenename$external_gene_name != "")
  rownames(RPKMHeatGenename) <- RPKMHeatGenename$external_gene_name
  RPKMHeatGenename <- RPKMHeatGenename[-c(9,10,11,12)] 
  matRPKMHeat <- as.matrix(RPKMHeatGenename)
  # preparing PDF file
  pdf(file = paste0("~/Desktop/BIOINFORMATICS_jacobi/graphs/HeatMaps/",ont,"/",updown,"/",sprintf("%004d", i),"_" , str_replace_all(allResup2$GO.ID[i],":","_"),"_",str_replace_all(allResup2$Term[i]," ","_"),".pdf"),height = 2 + (nrow(matRPKMHeat)*0.08))
  #making the heatmap visual
  my_heatmap <- pheatmap(matRPKMHeat,
                         cluster_rows = FALSE,
                         color = viridis(n = 2000,option = "magma"),
                                            legend_breaks = c(-2,0,2),
                         show_rownames = TRUE,
                         show_colnames = TRUE,
                         border_color = NA,
                         main = paste(gTitle," ",ont,"\n","GO ID : ",allResup2$GO.ID[i],"\n", allResUp$Term[i]),
                         scale = "row",
                         fontsize_row = 5,
                         cellheight = 5,
                         width=20000, 
                         height=20000
                         
  )
  
  
  
 my_heatmap
 dev.off() # printing PDF file
  
  # ------
  i <- i+1
} #end of for loop 
  #saving excel file 
  saveWorkbook(wb, paste0("~/Desktop/BIOINFORMATICS_jacobi/graphs/HeatMaps/",ont , "/",updown ,"/genes_per_GO_ID.xlsx"), overwrite = TRUE)

}  # end of function

funcGO("BP",20,0.05,"up")
funcGO("BP",20,0.05,"down")
funcGO("BP",20,0.05,"both")


funcGO("MF",20,0.05,"up")
funcGO("MF",20,0.05,"down")
funcGO("MF",20,0.05,"both")


funcGO("CC",20,0.05,"up")
funcGO("CC",20,0.05,"down")
funcGO("CC",20,0.05,"both")





