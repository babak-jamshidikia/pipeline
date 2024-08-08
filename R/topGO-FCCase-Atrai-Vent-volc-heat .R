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
RPKM5 <- RPKM4
RPKM6 <- log10(RPKM5)

colnames(RPKM6) <-c("Atria_1","Atria_2","Atria_3","Vent_1","Vent_2","Vent_3") #substr(colnames(RPKM6),8,14)

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
fit1 <- glmQLFit(y,design,robust = TRUE)
plotQLDisp(fit1)
lrt <-glmLRT(fit)

tt <- topTags(lrt,n=nrow(lrt), p.value=0.05)
tt2 <- topTags(lrt,n=nrow(lrt),p.value = 1)
# for tt
ensembly1 <- useEnsembl( biomart = "genes",dataset= "mmusculus_gene_ensembl")
attriby1 <- c("external_gene_name","description","ensembl_gene_id")


filtersy1 <- "ensembl_gene_id"
valuey1 <-list(tt$table$GeneID)

datay1 <-getBM(attributes = attriby1 ,filters = filtersy1 , values = valuey1 ,mart =  ensembly1)

Ddatay1 <- duplicated(datay1$ensembl_gene_id)
datay1 <-datay1[!Ddatay1,]

tt <- left_join(tt$table,datay1,by = c("GeneID" = "ensembl_gene_id"))
tt2 <- left_join(tt2$table,datay1,by = c("GeneID" = "ensembl_gene_id"))



colnames(design)
head(tt)


dfTT <- as.data.frame(tt)
dfTT$FDR <- -log10(dfTT$FDR)
dfTT
dfTT$diffexpressed <- "NO"
dfTT
dfTT$diffexpressed[dfTT$logFC > 0 & dfTT$FDR > -log10(0.05)] <- "UP"
dfTT$diffexpressed[dfTT$logFC < 0 & dfTT$FDR > -log10(0.05)] <- "DOWN"


dfTT$delabel <- ifelse(dfTT$external_gene_name %in% head(dfTT[order(dfTT$logFC), "external_gene_name"], 10), dfTT$external_gene_name , NA)

# for tt2

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
dfTT2

