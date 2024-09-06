#----------------------------------------------------------------------------------------------
#                           function list 
#   funchaetmap()           makes heatmap graph 
#   funcGlimavolc(varfilename)                          makes interactive volcano plot
#   funcVolcano()                                       makes volcano graph
#   funcGOBasicplot(ont,TopN,P_valcutoff,updown)        makes bbasic plot based on GO data
#   funcGOGroupGene(ont,TopN,P_valcutoff,updown)        makeg heatmap and excel file based on  GOid in allresults table
#
#-----------------------------------------------------------------------------------------------

library(tidyverse)
library(edgeR)
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

fc <- readRDS("/home/bjamshidikia/srvrdata/fcdataHuman.rds")
head(fc$counts)
head(fc)
y <- DGEList(counts = fc$counts,genes = fc$annotation,group = c("SiC","SiC","SiC","SiC","SiC_TM","SiC_TM","SiC_TM","SiC_TM"))
RPKM <-rpkm(y)
RPKM1 <- as.data.frame(RPKM)
#RPKM3 <- RPKM1 %>% 
 # filter(if_any(everything(), ~ .x>0))
RPKM4 <- filter_all(RPKM1,any_vars(. > 0))



RPKM4

RPKM5 <- RPKM4

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
tt2 <- topTags(lrt,n=nrow(lrt),p.value = 1)
# for tt
ensembly1 <- useEnsembl( biomart = "genes",dataset= "hsapiens_gene_ensembl")
attriby1 <- c("external_gene_name","description","ensembl_gene_id")



#tt <- left_join(tt$table,datay1,by = c("GeneID" = "ensembl_gene_id"))
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


funchaetmap <- function(){
  RPKM6 <- log10(RPKM5)
  colnames(RPKM6) <-c("SiC_1","SiC_2","SiC_3","SiC_4","SiC_TM_1","SiC_TM_2","SiC_TM_3","SiC_TM_4") #substr(colnames(RPKM6),8,14)
  matRPKM <- as.matrix(RPKM6)
  my_heatmap <- pheatmap(matRPKM,
                         cluster_rows = FALSE,
                         color = viridis(n = 2000,option = "magma"),
                         #                       legend_breaks = c(-2,0,2),
                         show_rownames = FALSE,
                         show_colnames = TRUE,
                         border_color = NA,
                         main = "Atria vs Vent  \n Heatmap based on Log10(RPKM) values ",
                         scale = "row",
                         
  )
  my_heatmap
}

funcGlimavolc <- function(varfilename){
  
  #glimmaVolcano(lrt , dge = y,main = paste("Atria", "vs", "Vent"),status.cols = c("dodgerblue", "gray", "firebrick") )
  
  glimmaVolcano(
    lrt,
    dge = y,
    counts = y$counts,
    groups = y$samples$group,
    #status = edgeR::decideTestsDGE(lrt),
    status = edgeR::decideTests.DGEExact(lrt),
    anno =lrt$genes ,
    display.columns = NULL,
    status.cols = c("dodgerblue", "gray", "firebrick"),
    sample.cols = NULL,
    p.adj.method = "BH",
    transform.counts = c("logcpm", "cpm", "rpkm", "none"),
    main = paste("SiC", "vs", "SiC_TM"),
    xlab = "logFC",
    ylab = "negLog10PValue",
    html = paste0("~/Desktop/BIOINFORMATICS_jacobi/graphs/glimma/",varfilename,".html"),
    width = 1200,
    height = 1200,
    
  )
  
  #plothtml <- glimmaMA(lrt, dge=y, status = edgeR::decideTests.DGEExact(lrt),status.colours=c("#3977db","#3d3f42","#db0d4e"))
  #htmlwidgets::saveWidget(
  #  plothtml,
  #  paste0("~/Desktop/BIOINFORMATICS_jacobi/graphs/glimma/",varfilename,"MA.html"),
  #  selfcontained = TRUE,
  #  libdir = NULL,
  #  background = "white",
  #  title = "MA plot"  #class(widget)[[1]],
  #knitrOptions = list()
  #)
  #glimmaXY(x=lrt$table$logFC, y=lrt$table$logCPM , dge=y,status.cols = c("dodgerblue", "silver", "firebrick"))
} #end of function Glima 




funcVolcano <- function(){
  filtersy1 <- "ensembl_gene_id"
  valuey1 <-list(tt$table$GeneID)
  
  datay1 <-getBM(attributes = attriby1 ,filters = filtersy1 , values = valuey1 ,mart =  ensembly1)
  
  Ddatay1 <- duplicated(datay1$ensembl_gene_id)
  datay1 <-datay1[!Ddatay1,]
  
tt2 <- left_join(tt2$table,datay1,by = c("GeneID" = "ensembl_gene_id"))
dfTT2 <- as.data.frame(tt2)
dfTT2$LOGFDR <- -log10(dfTT2$FDR)
dfTT2
dfTT2$diffexpressed <- "NO"
dfTT2
dfTT2$diffexpressed[dfTT2$logFC > 0 & dfTT2$FDR < (0.05)] <- "UP"
dfTT2$diffexpressed[dfTT2$logFC < 0 & dfTT2$FDR < (0.05)] <- "DOWN"

nUpdfTT2 <- length(dfTT2$diffexpressed[dfTT2$logFC > 0 & dfTT2$FDR < (0.05)])
nDwdfTT2 <- length(dfTT2$diffexpressed[dfTT2$logFC < 0 & dfTT2$FDR < (0.05)])
dfTT2$delabel <- ifelse(dfTT2$external_gene_name %in% head(dfTT2[order(dfTT2$logFC,decreasing =  TRUE), "external_gene_name"], 10), dfTT2$external_gene_name , NA)
dfTT2$delabel <- ifelse(dfTT2$external_gene_name %in% head(dfTT2[order(dfTT2$logFC ,decreasing = FALSE), "external_gene_name"], 13), dfTT2$external_gene_name ,dfTT2$delabel)
p2 <- ggplot(dfTT2,aes(x= logFC,y= LOGFDR,col = diffexpressed,label = delabel)) +
  geom_vline(xintercept = 0,col = "black",linetype = "dashed")+
  geom_hline(yintercept = -log10(0.05),col= "black",linetype = "dashed")+
  geom_point(size = 2,alpha = 0.2)+
  scale_color_manual(values = c("#00AFBB", "black", "red"), 
                     labels = c("Downregulated", "Not significant", "Upregulated"))+
  coord_cartesian(ylim = c(-20, 350), xlim = c(-15, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Severe', #legend_title, 
       x = expression("log"[2]* "(Fold Change)"),
       y = expression("-log"[10]*"(FDR)")) + 
  #scale_x_continuous(breaks = seq(-15, 10, 5))+
  #scale_y_continuous(breaks = seq(0,5,1))+
  ggtitle("Volcano Plot") +
  #geom_text_repel(color = "black", max.overlaps = Inf,xlim = c(-Inf, Inf), ylim = c(-Inf, Inf))+
  geom_label_repel(fill = "white",max.overlaps = Inf, xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),color = "black",box.padding = 0.25,alpha = 0.7)
  

p2
#dfTT2
}


funcGOBasicplot <-function(ont,TopN,P_valcutoff,updown){

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

#wb <- createWorkbook()

#addWorksheet(wb, sheetName = "allresUP");
#writeDataTable(wb, sheet = 1, x = allResUp, rowNames = TRUE);

#addWorksheet(wb, sheetName = "allresdwn");
#writeDataTable(wb, sheet = 2, x = allResdown, rowNames = TRUE);

#addWorksheet(wb, sheetName = "allresBoth");
#writeDataTable(wb, sheet = 3, x = allResboth, rowNames = TRUE);
#saveWorkbook(wb, "~/Desktop/BIOINFORMATICS_jacobi/ALLRES.xlsx", overwrite = TRUE)

}

funcGOGroupGene <-function(ont,TopN,P_valcutoff,updown){
#  ont = "BP"
#  TopN = 20
#  P_valcutoff = 0.05
#  updown = "up"
  
  if(!dir.exists(paste0("~/Desktop/BIOINFORMATICS_jacobi/graphs/HumanHeatmaps/",ont))){
    system(paste0("mkdir ~/Desktop/BIOINFORMATICS_jacobi/graphs/HumanHeatmaps/",ont))
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
  
  wbxl <- createWorkbook()  
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
  
  
  #------------------------------------  
  
  
  #---------------------
  
  #print(updown)
  if(!dir.exists(paste0("~/Desktop/BIOINFORMATICS_jacobi/graphs/HumanHeatmaps/",ont,"/",updown))){
    system(paste0("mkdir ~/Desktop/BIOINFORMATICS_jacobi/graphs/HumanHeatmaps/",ont,"/",updown))
  }
  
  gendf <- as.data.frame(geneList2)
  
  ensembl <- useEnsembl( biomart = "genes",dataset= "hsapiens_gene_ensembl")
  godata2 <- getBM(attributes = c("ensembl_gene_id","external_gene_name"), filters = "ensembl_gene_id", values =row.names(gendf), mart = ensembl)
  
 
  
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
    
    
    addWorksheet(wbxl, sheetName = paste0(i,"_", substr(allResUp$Term[i],1,28) ))
    writeDataTable(wbxl, sheet = i, x = TermDataEXL, rowNames = FALSE)
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
    RPKMHeatGenename <- RPKMHeatGenename[-c(7,8)] 
    matRPKMHeat <- as.matrix(RPKMHeatGenename)
    # preparing PDF file
    pdf(file = paste0("~/Desktop/BIOINFORMATICS_jacobi/graphs/HumanHeatmaps/",ont,"/",updown,"/",sprintf("%004d", i),"_" , str_replace_all(allResup2$GO.ID[i],":","_"),"_",str_replace_all(allResup2$Term[i]," ","_"),".pdf"),height = 2 + (nrow(matRPKMHeat)*0.08))
    #making the heatmap visual
    my_heatmap <- pheatmap(matRPKMHeat,
                           cluster_rows = FALSE,
                           color = viridis(n = 2000,option = "magma"),
                           #                   legend_breaks = c(-2,0,2),
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
  saveWorkbook(wbxl, paste0("~/Desktop/BIOINFORMATICS_jacobi/graphs/HumanHeatmaps/",ont , "/",updown ,"/genes_per_GO_ID.xlsx"), overwrite = TRUE)
  
}  # end of function

funchaetmap()
funcGlimavolc(varfilename = "H1")
funcGOBasicplot("BP",20,0.05,"up")
funcGOGroupGene("BP",20,0.05,"up")
