devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")
BiocManager::install("TCGAbiolinks")
setwd("~/PhD")

sort( sapply(ls(),function(x){object.size(get(x))}))


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::valid()

#BiocManager::install(c("TCGAbiolinks","RTCGAToolbox","downloader","readr","gaia","DMRcate","WGCNA","clusterProfiler","pathview","SummarizedExperiment","VennDiagram","stringr","tictoc","GOSemSim","pheatmap","edgeR","DESeq2","ELMER","parallel","MultiAssayExperiment","seqlm","IlluminaHumanMethylation450kanno.ilmn12.hg19","biomaRt","fishplot","FDb.InfiniumMethylation.hg19","circlize"))

library(TCGAbiolinks)
library(SummarizedExperiment)
library(RTCGAToolbox)
library(downloader)
library(readr)
library(gaia)
library(DMRcate)
library(WGCNA)
library(clusterProfiler)
library(pathview)
library(VennDiagram)
library(stringr)
library(tictoc)
library(GOSemSim)
library(pheatmap)
library(edgeR)
library(DESeq2)
library(ELMER)
library(parallel)
library(MultiAssayExperiment)
library(seqlm)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(biomaRt)
library(fishplot)
library(FDb.InfiniumMethylation.hg19)
library(circlize)
library(EDASeq)
library(DescTools)
library(CeRNASeek)
library(ggpubr)

# DMR

p_h <- "TCGA-ESCA"
#p_h <- "TCGA-HNSC"

#t_c <- "Primary solid Tumor"
t_c <- "Primary Tumor"
t_n <- "Solid Tissue Normal"

stage_n <- c("not reported")
stage_1 <- c("stage i","stage ia","stage ib")
stage_2 <- c("stage ii","stage iia","stage iib")
stage_3 <- c("stage iii","stage iiia","stage iiib","stage iiic")
stage_4 <- c("stage iv","stage iva","stage ivb","stage ivc")

################
load("meth.RData")
load("genelist.RData")
load("ansEA.RData")
load("dataDEG.RData")
load("rse.RData")

load(paste0(p_h,"_meth.RData"))
load(paste0(p_h,"_rse.RData"))
load(paste0(p_h,"dataDEG.RData"))
load(paste0(p_h,"_genelist.RData"))
load(paste0(p_h,"_ansEA.RData"))
load(paste0(p_h,"_datafile_PCC.RData"))

############### Data Clinical

dataClin_h <- GDCquery_clinic(p_h, "Clinical")

table(dataClin_h$ajcc_clinical_stage)
table(dataClin_h$tumor_stage)

dataClin_h_n <- subset(dataClin_h, tumor_stage %in% stage_n)
dataClin_h_1 <- subset(dataClin_h, tumor_stage %in% stage_1)
dataClin_h_2 <- subset(dataClin_h, tumor_stage %in% stage_2)
dataClin_h_3 <- subset(dataClin_h, tumor_stage %in% stage_3)
dataClin_h_4 <- subset(dataClin_h, tumor_stage %in% stage_4)
dataClin_h_low <- subset(dataClin_h, tumor_stage %in% c(stage_1,stage_2,stage_3))
dataClin_h_s <- subset(dataClin_h, tumor_stage %in% c(stage_1,stage_2,stage_3,stage_4))

dataClin_G1 <- subset(dataClin_h, tumor_stage %in% c(stage_1,stage_2))
dataClin_G2 <- subset(dataClin_h, tumor_stage %in% c(stage_3,stage_4))
######## CNV

lastRunDate <- getFirehoseRunningDates()[1]
lastAnalyseDate <- getFirehoseAnalyzeDates(1)

ESCA.data <-  getFirehoseData(dataset = "ESCA",runDate = lastRunDate, gistic2_Date = lastAnalyseDate,GISTIC = TRUE,
                              fileSizeLimit = 10000)

ESCA.data@GISTIC@ThresholdedByGene$AMP <- apply(ESCA.data@GISTIC@ThresholdedByGene,1, function(y) sum(length(which(y>0))))
ESCA.data@GISTIC@ThresholdedByGene$DEL <- apply(ESCA.data@GISTIC@ThresholdedByGene,1, function(y) sum(length(which(y<0))))

cnv_h <- ESCA.data@GISTIC@ThresholdedByGene
rownames(cnv_h) <- cnv_h$Gene.Symbol
colnames(cnv_h) <- str_replace_all(substr(colnames(cnv_h),1,12),"[.]","-")
cnv_h_n <- cnv_h[,colnames(cnv_h) %in% substr(colnames(rse_h_n),1,12)]
cnv_h_c <- cnv_h[,colnames(cnv_h) %in% substr(colnames(rse_h_c),1,12)]
cnv_h_c_1 <- cnv_h[,colnames(cnv_h) %in% substr(colnames(rse_h_c_1),1,12)]
cnv_h_c_2 <- cnv_h[,colnames(cnv_h) %in% substr(colnames(rse_h_c_2),1,12)]
cnv_h_c_3 <- cnv_h[,colnames(cnv_h) %in% substr(colnames(rse_h_c_3),1,12)]
cnv_h_c_4 <- cnv_h[,colnames(cnv_h) %in% substr(colnames(rse_h_c_4),1,12)]
cnv_h_c_n <- cnv_h[,colnames(cnv_h) %in% substr(colnames(rse_h_c_n),1,12)]
cnv_h_c_low <- cnv_h[,colnames(cnv_h) %in% substr(colnames(rse_h_c_low),1,12)]
cnv_g1 <- cnv_h[,colnames(cnv_h) %in% substr(colnames(rse_G1),1,12)]
cnv_g2 <- cnv_h[,colnames(cnv_h) %in% substr(colnames(rse_G1),1,12)]

cnv_h_n$`Gene-Symbol` <- cnv_h$`Gene-Symbol`
cnv_h_n$`Locus-ID` <- cnv_h$`Locus-ID`
cnv_h_n$Cytoband <- cnv_h$Cytoband
cnv_h_c$`Gene-Symbol` <- cnv_h$`Gene-Symbol`
cnv_h_c$`Locus-ID` <- cnv_h$`Locus-ID`
cnv_h_c$Cytoband <- cnv_h$Cytoband
cnv_h_c_1$`Gene-Symbol` <- cnv_h$`Gene-Symbol`
cnv_h_c_1$`Locus-ID` <- cnv_h$`Locus-ID`
cnv_h_c_1$Cytoband <- cnv_h$Cytoband
cnv_h_c_2$`Gene-Symbol` <- cnv_h$`Gene-Symbol`
cnv_h_c_2$`Locus-ID` <- cnv_h$`Locus-ID`
cnv_h_c_2$Cytoband <- cnv_h$Cytoband
cnv_h_c_3$`Gene-Symbol` <- cnv_h$`Gene-Symbol`
cnv_h_c_3$`Locus-ID` <- cnv_h$`Locus-ID`
cnv_h_c_3$Cytoband <- cnv_h$Cytoband
cnv_h_c_4$`Gene-Symbol` <- cnv_h$`Gene-Symbol`
cnv_h_c_4$`Locus-ID` <- cnv_h$`Locus-ID`
cnv_h_c_4$Cytoband <- cnv_h$Cytoband
cnv_h_c_n$`Gene-Symbol` <- cnv_h$`Gene-Symbol`
cnv_h_c_n$`Locus-ID` <- cnv_h$`Locus-ID`
cnv_h_c_n$Cytoband <- cnv_h$Cytoband
cnv_h_c_low$`Gene-Symbol` <- cnv_h$`Gene-Symbol`
cnv_h_c_low$`Locus-ID` <- cnv_h$`Locus-ID`
cnv_h_c_low$Cytoband <- cnv_h$Cytoband

cnv_g1$`Gene-Symbol` <- cnv_h$`Gene-Symbol`
cnv_g1$`Locus-ID` <- cnv_h$`Locus-ID`
cnv_g1$Cytoband <- cnv_h$Cytoband
cnv_g2$`Gene-Symbol` <- cnv_h$`Gene-Symbol`
cnv_g2$`Locus-ID` <- cnv_h$`Locus-ID`
cnv_g2$Cytoband <- cnv_h$Cytoband

cnv_h_n$AMP <- apply(cnv_h_n[,colnames(cnv_h_n) %in% substr(colnames(rse_h_n),1,12)],1, function(y) sum(length(which(y>0))))
cnv_h_n$DEL <- apply(cnv_h_n[,colnames(cnv_h_n) %in% substr(colnames(rse_h_n),1,12)],1, function(y) sum(length(which(y<0))))
cnv_h_n$NC <- apply(cnv_h_n[,colnames(cnv_h_n) %in% substr(colnames(rse_h_n),1,12)],1, function(y) sum(length(which(y==0))))
cnv_h_c$AMP <- apply(cnv_h_c[,colnames(cnv_h_c) %in% substr(colnames(rse_h_c),1,12)],1, function(y) sum(length(which(y>0))))
cnv_h_c$DEL <- apply(cnv_h_c[,colnames(cnv_h_c) %in% substr(colnames(rse_h_c),1,12)],1, function(y) sum(length(which(y<0))))
cnv_h_c$NC <- apply(cnv_h_c[,colnames(cnv_h_c) %in% substr(colnames(rse_h_c),1,12)],1, function(y) sum(length(which(y==0))))
cnv_h_c_1$AMP <- apply(cnv_h_c_1[,colnames(cnv_h_c_1) %in% substr(colnames(rse_h_c_1),1,12)],1, function(y) sum(length(which(y>0))))
cnv_h_c_1$DEL <- apply(cnv_h_c_1[,colnames(cnv_h_c_1) %in% substr(colnames(rse_h_c_1),1,12)],1, function(y) sum(length(which(y<0))))
cnv_h_c_1$NC <- apply(cnv_h_c_1[,colnames(cnv_h_c_1) %in% substr(colnames(rse_h_c_1),1,12)],1, function(y) sum(length(which(y==0))))
cnv_h_c_2$AMP <- apply(cnv_h_c_2[,colnames(cnv_h_c_2) %in% substr(colnames(rse_h_c_2),1,12)],1, function(y) sum(length(which(y>0))))
cnv_h_c_2$DEL <- apply(cnv_h_c_2[,colnames(cnv_h_c_2) %in% substr(colnames(rse_h_c_2),1,12)],1, function(y) sum(length(which(y<0))))
cnv_h_c_2$NC <- apply(cnv_h_c_2[,colnames(cnv_h_c_2) %in% substr(colnames(rse_h_c_2),1,12)],1, function(y) sum(length(which(y==0))))
cnv_h_c_3$AMP <- apply(cnv_h_c_3[,colnames(cnv_h_c_3) %in% substr(colnames(rse_h_c_3),1,12)],1, function(y) sum(length(which(y>0))))
cnv_h_c_3$DEL <- apply(cnv_h_c_3[,colnames(cnv_h_c_3) %in% substr(colnames(rse_h_c_3),1,12)],1, function(y) sum(length(which(y<0))))
cnv_h_c_3$NC <- apply(cnv_h_c_3[,colnames(cnv_h_c_3) %in% substr(colnames(rse_h_c_3),1,12)],1, function(y) sum(length(which(y==0))))
cnv_h_c_4$AMP <- apply(cnv_h_c_4[,colnames(cnv_h_c_4) %in% substr(colnames(rse_h_c_4),1,12)],1, function(y) sum(length(which(y>0))))
cnv_h_c_4$DEL <- apply(cnv_h_c_4[,colnames(cnv_h_c_4) %in% substr(colnames(rse_h_c_4),1,12)],1, function(y) sum(length(which(y<0))))
cnv_h_c_4$NC <- apply(cnv_h_c_4[,colnames(cnv_h_c_4) %in% substr(colnames(rse_h_c_4),1,12)],1, function(y) sum(length(which(y==0))))
cnv_h_c_n$AMP <- apply(cnv_h_c_n[,colnames(cnv_h_c_n) %in% substr(colnames(rse_h_c_n),1,12)],1, function(y) sum(length(which(y>0))))
cnv_h_c_n$DEL <- apply(cnv_h_c_n[,colnames(cnv_h_c_n) %in% substr(colnames(rse_h_c_n),1,12)],1, function(y) sum(length(which(y<0))))
cnv_h_c_n$NC <- apply(cnv_h_c_n[,colnames(cnv_h_c_n) %in% substr(colnames(rse_h_c_n),1,12)],1, function(y) sum(length(which(y==0))))
cnv_h_c_low$AMP <- apply(cnv_h_c_low[,colnames(cnv_h_c_low) %in% substr(colnames(rse_h_c_low),1,12)],1, function(y) sum(length(which(y>0))))
cnv_h_c_low$DEL <- apply(cnv_h_c_low[,colnames(cnv_h_c_low) %in% substr(colnames(rse_h_c_low),1,12)],1, function(y) sum(length(which(y<0))))
cnv_h_c_low$NC <- apply(cnv_h_c_low[,colnames(cnv_h_c_low) %in% substr(colnames(rse_h_c_low),1,12)],1, function(y) sum(length(which(y==0))))

cnv_g1$AMP <- apply(cnv_g1[,colnames(cnv_g1) %in% substr(colnames(rse_G1),1,12)],1, function(y) sum(length(which(y>0))))
cnv_g1$DEL <- apply(cnv_g1[,colnames(cnv_g1) %in% substr(colnames(rse_G1),1,12)],1, function(y) sum(length(which(y<0))))
cnv_g1$NC <- apply(cnv_g1[,colnames(cnv_g1) %in% substr(colnames(rse_G1),1,12)],1, function(y) sum(length(which(y==0))))
cnv_g2$AMP <- apply(cnv_g2[,colnames(cnv_g2) %in% substr(colnames(rse_G2),1,12)],1, function(y) sum(length(which(y>0))))
cnv_g2$DEL <- apply(cnv_g2[,colnames(cnv_g2) %in% substr(colnames(rse_G2),1,12)],1, function(y) sum(length(which(y<0))))
cnv_g2$NC <- apply(cnv_g2[,colnames(cnv_g2) %in% substr(colnames(rse_G2),1,12)],1, function(y) sum(length(which(y==0))))

###############

# Old version
#
#hyper <- met.dmr.h.n.c1[values(met.dmr.h.n.c1)[,"status.Primary.Tumor.Solid.Tissue.Normal"]=="Hypermethylated"]
#hypo <- met.dmr.h.n.c1[values(met.dmr.h.n.c1)[,"status.Primary.Tumor.Solid.Tissue.Normal"]=="Hypomethylated"]
#
#hyper_h_n_c1 <- sort(unique(unlist(strsplit(capture.output(cat(values(hyper)[,"Gene_Symbol"],sep=";")),";"))))
#hypo_h_n_c1 <- sort(unique(unlist(strsplit(capture.output(cat(values(hypo)[,"Gene_Symbol"],sep=";")),";"))))
#

met <- met.dmr.h.n.c1
rows <- rownames(met[met$status=="Hypermethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_c)[rowData(met_h_c)$probeID %in% rows,]$Gene_Symbol
hyper_h_n_c1 <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))
rows <- rownames(met[met$status=="Hypomethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_c)[rowData(met_h_c)$probeID %in% rows,]$Gene_Symbol
hypo_h_n_c1 <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))

met <- met.dmr.h.n.c2
rows <- rownames(met[met$status=="Hypermethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_c)[rowData(met_h_c)$probeID %in% rows,]$Gene_Symbol
hyper_h_n_c2 <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))
rows <- rownames(met[met$status=="Hypomethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_c)[rowData(met_h_c)$probeID %in% rows,]$Gene_Symbol
hypo_h_n_c2 <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))

met <- met.dmr.h.n.c3
rows <- rownames(met[met$status=="Hypermethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_c)[rowData(met_h_c)$probeID %in% rows,]$Gene_Symbol
hyper_h_n_c3 <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))
rows <- rownames(met[met$status=="Hypomethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_c)[rowData(met_h_c)$probeID %in% rows,]$Gene_Symbol
hypo_h_n_c3 <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))

met <- met.dmr.h.n.c4
rows <- rownames(met[met$status=="Hypermethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_c)[rowData(met_h_c)$probeID %in% rows,]$Gene_Symbol
hyper_h_n_c4 <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))
rows <- rownames(met[met$status=="Hypomethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_c)[rowData(met_h_c)$probeID %in% rows,]$Gene_Symbol
hypo_h_n_c4 <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))

met <- met.dmr.h.n.c.low
rows <- rownames(met[met$status=="Hypermethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_c)[rowData(met_h_c)$probeID %in% rows,]$Gene_Symbol
hyper_h_n_c_low <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))
rows <- rownames(met[met$status=="Hypomethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_c)[rowData(met_h_c)$probeID %in% rows,]$Gene_Symbol
hypo_h_n_c_low <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))

hyper_probes_c4 <- rownames(met.dmr.h.n.c4[met.dmr.h.n.c4$status=="Hypermethylated in Solid Tissue Normal",])
hypo_probes_c4 <- rownames(met.dmr.h.n.c4[met.dmr.h.n.c4$status=="Hypomethylated in Solid Tissue Normal",])
hyper_probes_low <- rownames(met.dmr.h.n.c.low[met.dmr.h.n.c.low$status=="Hypermethylated in Solid Tissue Normal",])
hypo_probes_low <- rownames(met.dmr.h.n.c.low[met.dmr.h.n.c.low$status=="Hypomethylated in Solid Tissue Normal",])


up_h_c1 <- rownames(dataDEGs_h_n_c1[dataDEGs_h_n_c1$logFC>0,])
down_h_c1 <- rownames(dataDEGs_h_n_c1[dataDEGs_h_n_c1$logFC<0,])

up_h_c2 <- rownames(dataDEGs_h_n_c2[dataDEGs_h_n_c2$logFC>0,])
down_h_c2 <- rownames(dataDEGs_h_n_c2[dataDEGs_h_n_c2$logFC<0,])

up_h_c3 <- rownames(dataDEGs_h_n_c3[dataDEGs_h_n_c3$logFC>0,])
down_h_c3 <- rownames(dataDEGs_h_n_c3[dataDEGs_h_n_c3$logFC<0,])

up_h_c4 <- rownames(dataDEGs_h_n_c4[dataDEGs_h_n_c4$logFC>0,])
down_h_c4 <- rownames(dataDEGs_h_n_c4[dataDEGs_h_n_c4$logFC<0,])

up_h_c_low <- rownames(dataDEGs_h_n_c_low[dataDEGs_h_n_c_low$logFC>0,])
down_h_c_low <- rownames(dataDEGs_h_n_c_low[dataDEGs_h_n_c_low$logFC<0,])

up_h_c1_c2 <- rownames(dataDEGs_h_c1_c2[dataDEGs_h_c1_c2$logFC>0,])
down_h_c1_c2 <- rownames(dataDEGs_h_c1_c2[dataDEGs_h_c1_c2$logFC<0,])

up_h_c2_c3 <- rownames(dataDEGs_h_c2_c3[dataDEGs_h_c2_c3$logFC>0,])
down_h_c2_c3 <- rownames(dataDEGs_h_c2_c3[dataDEGs_h_c2_c3$logFC<0,])

up_h_c3_c4 <- rownames(dataDEGs_h_c3_c4[dataDEGs_h_c3_c4$logFC>0,])
down_h_c3_c4 <- rownames(dataDEGs_h_c3_c4[dataDEGs_h_c3_c4$logFC<0,])

up_h_c_low_c4 <- rownames(dataDEGs_h_c_low_c4[dataDEGs_h_c_low_c4$logFC>0,])
down_h_c_low_c4 <- rownames(dataDEGs_h_c_low_c4[dataDEGs_h_c_low_c4$logFC<0,])

ups <- unique(c(up_h_c1,up_h_c2,up_h_c3,up_h_c4))
downs <- unique(c(down_h_c1,down_h_c2,down_h_c3,down_h_c4))
ups_s <- unique(c(up_h_c1_c2,up_h_c2_c3,up_h_c3_c4))
downs_s <- unique(c(down_h_c1_c2,down_h_c2_c3,down_h_c3_c4))



ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes UPs"), RegulonList = ups)
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = ups,
                        nBar = 20,filename = paste0("DEA_genes_UPs.pdf"))

ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes DOWNs"), RegulonList = downs)
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = downs,
                        nBar = 20,filename = paste0("DEA_genes_DOWNs.pdf"))

ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes UP_Ss"), RegulonList = ups_s)
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = ups_s,
                        nBar = 20,filename = paste0("DEA_genes_UP_Ss.pdf"))

ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes DOWN_Ss"), RegulonList = downs_s)
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = downs_s,
                        nBar = 20,filename = paste0("DEA_genes_DOWN_Ss.pdf"))

ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes stage"), RegulonList = union(ups_s,downs_s))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = union(ups_s,downs_s),
                        nBar = 20,filename = paste0("DEA_genes_Ss.pdf"))

hypers <- unique(c(hyper_h_n_c1,hyper_h_n_c2,hyper_h_n_c3,hyper_h_n_c4))
hypos <- unique(c(hypo_h_n_c1,hypo_h_n_c2,hypo_h_n_c3,hypo_h_n_c4))

up_h_c1_c2[up_h_c1_c2 %in% down_h_c2]
up_h_c1_c2[up_h_c1_c2 %in% up_h_c2]
up_h_c1_c2[!up_h_c1_c2 %in% c(up_h_c2,down_h_c2)]

up_h_c_low_c4[up_h_c_low_c4 %in% down_h_c4]
up_h_c_low_c4[up_h_c_low_c4 %in% up_h_c4]
up_h_c_low_c4[!up_h_c_low_c4 %in% c(up_h_c4,down_h_c4)]

up_h_c3_c4[up_h_c3_c4 %in% down_h_c4]
up_h_c3_c4[up_h_c3_c4 %in% up_h_c4]
up_h_c3_c4[!up_h_c3_c4 %in% c(up_h_c4,down_h_c4)]

up_h_c1_c2[up_h_c1_c2 %in% down_h_c1]
up_h_c1_c2[up_h_c1_c2 %in% up_h_c1]
up_h_c1_c2[!up_h_c1_c2 %in% c(up_h_c1,down_h_c1)]

up_h_c1_c2[up_h_c1_c2 %in% c(down_h_c1,down_h_c2)]
up_h_c1_c2[up_h_c1_c2 %in% c(up_h_c1,up_h_c2)]
up_h_c1_c2[!up_h_c1_c2 %in% c(up_h_c1,down_h_c1,up_h_c2,down_h_c2)]

up_h_c1_c2[up_h_c1_c2 %in% intersect(down_h_c1,down_h_c2)]
up_h_c1_c2[up_h_c1_c2 %in% intersect(up_h_c1,up_h_c2)]

down_h_c1_c2[down_h_c1_c2 %in% down_h_c2]
down_h_c1_c2[down_h_c1_c2 %in% up_h_c2]
down_h_c1_c2[!down_h_c1_c2 %in% c(up_h_c2,down_h_c2)]

down_h_c1_c2[down_h_c1_c2 %in% down_h_c1]
down_h_c1_c2[down_h_c1_c2 %in% up_h_c1]
down_h_c1_c2[!down_h_c1_c2 %in% c(up_h_c1,down_h_c1)]

down_h_c1_c2[down_h_c1_c2 %in% c(down_h_c1,down_h_c2)]
down_h_c1_c2[down_h_c1_c2 %in% c(up_h_c1,up_h_c2)]
down_h_c1_c2[!down_h_c1_c2 %in% c(up_h_c1,down_h_c1,up_h_c2,down_h_c2)]

down_h_c1_c2[down_h_c1_c2 %in% intersect(down_h_c1,down_h_c2)]
down_h_c1_c2[down_h_c1_c2 %in% intersect(up_h_c1,up_h_c2)]

hsGO2 <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont="MF", computeIC=FALSE)
genes <- up_h_c_low_c4
mgeneSim(genes, semData=hsGO2, measure="Wang", combine="BMA", verbose=FALSE)

dataDEGs_h_c_low_c4[selectedGenes[selectedGenes %in% hypo_h_n_c4 & !(selectedGenes%in% hypo_h_n_c_low)],]

cnv_h_c_low["BMP3",]




########### Circos HG19

query.maf.hg19 <- GDCquery(project = p_h, 
                           data.category = "Simple nucleotide variation", 
                           data.type = "Simple somatic mutation",
                           access = "open", 
                           file.type = "bcgsc.ca_ESCA.IlluminaHiSeq_DNASeq.1.somatic.maf",
                           legacy = TRUE)

GDCdownload(query.maf.hg19)
maf <- GDCprepare(query.maf.hg19)

mut <- maf
# Select only potentially damaging mutations
mut <- mut[mut$Variant_Classification %in% c("Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Frame_Shift_Del","Frame_Shift_Ins"),]
# Select recurrent mutations (identified in at least two samples)
mut.id <- paste0(mut$Chromosome,":",mut$Start_Position,"–",mut$End_Position,"|",mut$Reference_Allele)
mut <- cbind(mut.id, mut)
numSamples <- table(mut.id)
s.mut <- names(which(numSamples>=1))
# Prepare selected mutations data for circos plot
s.mut <- mut[mut$mut.id %in% s.mut,]
s.mut <- s.mut[,c("Chromosome","Start_Position","End_Position","Variant_Classification","Hugo_Symbol")]
s.mut <- unique(s.mut)
s.mut[,1] <- as.character(s.mut[,1])
s.mut[,4] <- as.character(s.mut[,4])
s.mut[,5] <- as.character(s.mut[,5])
typeNames <- unique(s.mut[,4])
type <- c(4:1)
names(type) <- typeNames[1:4]
Type <- type[s.mut[,4]]
s.mut <- cbind(s.mut,Type)
s.mut <- s.mut[,c(1:3,6,4,5)]
Chromosome <- sapply(s.mut[,1],function(x) paste0("chr",x))
s.mut <- cbind(Chromosome, s.mut[,-1])

# Load recurrent CNV data for selected cancer (e.g. "ESCA")
load("ESCA_CNV_results.rda")

# Prepare selected sample CNV data for circos plot
s.cnv <- as.data.frame(RecCNV[RecCNV[,"q-value"]<=10^-4,c(1:4,6)])
s.cnv <- s.cnv[,c(1,3,4,2)]
xidx <- which(s.cnv$Chromosome==23)
yidx <- which(s.cnv$Chromosome==24)
s.cnv[xidx,"Chromosome"] <- "X"
s.cnv[yidx,"Chromosome"] <- "Y"
Chromosome <- sapply(s.cnv[,1],function(x) paste0("chr",x))
s.cnv <- cbind(Chromosome, s.cnv[,-1])
s.cnv[,1] <- as.character(s.cnv[,1])
s.cnv[,4] <- as.character(s.cnv[,4])
s.cnv <- cbind(s.cnv,CNV=1)
colnames(s.cnv) <- c("Chromosome","Start_position","End_position","Aberration_Kind","CNV")

# Draw genomic circos plot
pdf("CircosPlot.pdf",width=15,height=15)
par(mar=c(1,1,1,1), cex=1)
circos.initializeWithIdeogram()
# Add CNV results
colors <- c("forestgreen","firebrick")
names(colors)  <- c(0,1)
circos.genomicTrackPlotRegion(s.cnv,  ylim = c(0,1.2),
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region, value, ytop.column = 2, ybottom = 0,
                                                   col = colors[value[[1]]],
                                                   border="white")
                                cell.xlim = get.cell.meta.data("cell.xlim")
                                circos.lines(cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                              })
# Add mutation results
#colors <- c("blue","green","red","gold")
#names(colors)  <- typeNames[1:4]
#circos.genomicTrackPlotRegion(s.mut, ylim = c(1.2,4.2),
#                              panel.fun = function(region, value, ...) {
#                                circos.genomicPoints(region, value, cex = 0.8, pch = 16, col = 
#                                                       colors[value[[2]]], ...)
#                              })

circos.clear()

legend(-0.2, 0.2, bty="n", y.intersp=1, c("Amp","Del"), pch=15, col=c("firebrick","forestgreen"), 
       title="CNVs", text.font=3, cex=1.2, title.adj=0)
#legend(-0.2, 0, bty="n", y.intersp=1, names(colors), pch=16, col=colors, title="Mutations", text.font
#       =3, cex=1.2, title.adj=0)
dev.off()


############# Circus plot G1

RecCNV <- RecCNV_G1

# Prepare selected sample CNV data for circos plot
s.cnv <- as.data.frame(RecCNV[RecCNV[,"q-value"]<=10^-4,c(1:4,6)])
s.cnv <- s.cnv[,c(1,3,4,2)]
xidx <- which(s.cnv$Chromosome==23)
yidx <- which(s.cnv$Chromosome==24)
s.cnv[xidx,"Chromosome"] <- "X"
s.cnv[yidx,"Chromosome"] <- "Y"
Chromosome <- sapply(s.cnv[,1],function(x) paste0("chr",x))
s.cnv <- cbind(Chromosome, s.cnv[,-1])
s.cnv[,1] <- as.character(s.cnv[,1])
s.cnv[,4] <- as.character(s.cnv[,4])
s.cnv <- cbind(s.cnv,CNV=1)
colnames(s.cnv) <- c("Chromosome","Start_position","End_position","Aberration_Kind","CNV")

# Draw genomic circos plot

pdf("CircosPlot_ESCA_G1.pdf",width=8,height=8)
par(mar=c(1,1,1,1), cex=1)
circos.initializeWithIdeogram()
# Add CNV results
colors <- c("forestgreen","firebrick")
names(colors)  <- c(0,1)
circos.genomicTrackPlotRegion(s.cnv,  ylim = c(0,1.2),
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region, value, ytop.column = 2, ybottom = 0,
                                                   col = colors[value[[1]]],
                                                   border="white")
                                cell.xlim = get.cell.meta.data("cell.xlim")
                                circos.lines(cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                              })
circos.clear()

legend(-0.2, 0.2, bty="n", y.intersp=1, c("Amplifications","Deletions"), pch=15, col=c("firebrick","forestgreen"), 
       title="ESCA Early Stage Samples\nCopy Number Variations", text.font=3, cex=1.2, title.adj=0)
dev.off()


RecCNV <- RecCNV_G2

# Prepare selected sample CNV data for circos plot
s.cnv <- as.data.frame(RecCNV[RecCNV[,"q-value"]<=10^-4,c(1:4,6)])
s.cnv <- s.cnv[,c(1,3,4,2)]
xidx <- which(s.cnv$Chromosome==23)
yidx <- which(s.cnv$Chromosome==24)
s.cnv[xidx,"Chromosome"] <- "X"
s.cnv[yidx,"Chromosome"] <- "Y"
Chromosome <- sapply(s.cnv[,1],function(x) paste0("chr",x))
s.cnv <- cbind(Chromosome, s.cnv[,-1])
s.cnv[,1] <- as.character(s.cnv[,1])
s.cnv[,4] <- as.character(s.cnv[,4])
s.cnv <- cbind(s.cnv,CNV=1)
colnames(s.cnv) <- c("Chromosome","Start_position","End_position","Aberration_Kind","CNV")

# Draw genomic circos plot

pdf("CircosPlot_ESCA_G2.pdf",width=8,height=8)
par(mar=c(1,1,1,1), cex=1)
circos.initializeWithIdeogram()
# Add CNV results
colors <- c("forestgreen","firebrick")
names(colors)  <- c(0,1)
circos.genomicTrackPlotRegion(s.cnv,  ylim = c(0,1.2),
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region, value, ytop.column = 2, ybottom = 0,
                                                   col = colors[value[[1]]],
                                                   border="white")
                                cell.xlim = get.cell.meta.data("cell.xlim")
                                circos.lines(cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                              })
circos.clear()

legend(-0.2, 0.2, bty="n", y.intersp=1, c("Amplifications","Deletions"), pch=15, col=c("firebrick","forestgreen"), 
       title="ESCA Metastasis Samples\nCopy Number Variations", text.font=3, cex=1.2, title.adj=0)
dev.off()



RecCNV <- RecCNV_h_low

# Prepare selected sample CNV data for circos plot
s.cnv <- as.data.frame(RecCNV[RecCNV[,"q-value"]<=10^-4,c(1:4,6)])
s.cnv <- s.cnv[,c(1,3,4,2)]
xidx <- which(s.cnv$Chromosome==23)
yidx <- which(s.cnv$Chromosome==24)
s.cnv[xidx,"Chromosome"] <- "X"
s.cnv[yidx,"Chromosome"] <- "Y"
Chromosome <- sapply(s.cnv[,1],function(x) paste0("chr",x))
s.cnv <- cbind(Chromosome, s.cnv[,-1])
s.cnv[,1] <- as.character(s.cnv[,1])
s.cnv[,4] <- as.character(s.cnv[,4])
s.cnv <- cbind(s.cnv,CNV=1)
colnames(s.cnv) <- c("Chromosome","Start_position","End_position","Aberration_Kind","CNV")

# Draw genomic circos plot
pdf("CircosPlot_ESCA_G2.pdf",width=15,height=15)
par(mar=c(1,1,1,1), cex=1)
circos.initializeWithIdeogram()
# Add CNV results
colors <- c("forestgreen","firebrick")
names(colors)  <- c(0,1)
circos.genomicTrackPlotRegion(s.cnv,  ylim = c(0,1.2),
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region, value, ytop.column = 2, ybottom = 0,
                                                   col = colors[value[[1]]],
                                                   border="white")
                                cell.xlim = get.cell.meta.data("cell.xlim")
                                circos.lines(cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                              })
circos.clear()

legend(-0.2, 0.2, bty="n", y.intersp=1, c("Amp","Del"), pch=15, col=c("firebrick","forestgreen"), 
       title="CNVs", text.font=3, cex=1.2, title.adj=0)
dev.off()



ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

myups <- up_h_c_low_c4[up_h_c_low_c4 %in% all_genes_met_h]
mydowns <- down_h_c_low_c4[down_h_c_low_c4 %in% all_genes_met_h]

myups[myups %in% hyper_h_n_c_low]
myups[myups %in% hyper_h_n_c1]
myups[myups %in% hyper_h_n_c2]
myups[myups %in% hyper_h_n_c3]
myups[myups %in% hyper_h_n_c4]

myups[myups %in% hypo_h_n_c_low]
myups[myups %in% hypo_h_n_c1]
myups[myups %in% hypo_h_n_c2]
myups[myups %in% hypo_h_n_c3]
myups[myups %in% hypo_h_n_c4]

hyper_probes_c4
ann450k$UCSC_RefGene_Name

dmr3 <- met.dmr.h.n.c3[!met.dmr.h.n.c3$status=="Not Significant",]
dmr3$Position = ann450k[ann450k$Name %in% row.names(dmr3),]$UCSC_RefGene_Group
dmr3$Name = ann450k[ann450k$Name %in% row.names(dmr3),]$Name
dmr3$Gene = ann450k[ann450k$Name %in% row.names(dmr3),]$UCSC_RefGene_Name
dmr3$PositionL = strsplit(ann450k[ann450k$Name %in% row.names(dmr3),]$UCSC_RefGene_Group,";")
dmr3$NameL = strsplit(ann450k[ann450k$Name %in% row.names(dmr3),]$Name,";")
dmr3$GeneL = strsplit(ann450k[ann450k$Name %in% row.names(dmr3),]$UCSC_RefGene_Name,";")

dmrl <- met.dmr.h.n.c.low[!met.dmr.h.n.c.low$status=="Not Significant",]
dmrl$Position = ann450k[ann450k$Name %in% row.names(dmrl),]$UCSC_RefGene_Group
dmrl$Name = ann450k[ann450k$Name %in% row.names(dmrl),]$Name
dmrl$Gene = ann450k[ann450k$Name %in% row.names(dmrl),]$UCSC_RefGene_Name
dmrl$PositionL = strsplit(ann450k[ann450k$Name %in% row.names(dmrl),]$UCSC_RefGene_Group,";")
dmrl$NameL = strsplit(ann450k[ann450k$Name %in% row.names(dmrl),]$Name,";")
dmrl$GeneL = strsplit(ann450k[ann450k$Name %in% row.names(dmrl),]$UCSC_RefGene_Name,";")

dmr4 <- met.dmr.h.n.c4[!met.dmr.h.n.c4$status=="Not Significant",]
dmr4$Position = ann450k[ann450k$Name %in% row.names(dmr4),]$UCSC_RefGene_Group
dmr4$Name = ann450k[ann450k$Name %in% row.names(dmr4),]$Name
dmr4$Gene = ann450k[ann450k$Name %in% row.names(dmr4),]$UCSC_RefGene_Name
dmr4$PositionL = strsplit(ann450k[ann450k$Name %in% row.names(dmr4),]$UCSC_RefGene_Group,";")
dmr4$NameL = strsplit(ann450k[ann450k$Name %in% row.names(dmr4),]$Name,";")
dmr4$GeneL = strsplit(ann450k[ann450k$Name %in% row.names(dmr4),]$UCSC_RefGene_Name,";")

up_h_c_low_c4[up_h_c_low_c4 %in% unlist(dmr4$GeneL)]
up_h_c_low_c4[up_h_c_low_c4 %in% unlist(dmrl$GeneL)]

down_h_c_low_c4[down_h_c_low_c4 %in% unlist(dmr4$GeneL)]
down_h_c_low_c4[down_h_c_low_c4 %in% unlist(dmrl$GeneL)]

up_h_c4[up_h_c4 %in% unlist(dmrl$GeneL)]
up_h_c4[up_h_c4 %in% unlist(dmr4$GeneL)]

uni <- union(union(up_h_c4,up_h_c_low),union(down_h_c_low,down_h_c4))
uni_dmr <- union(unlist(dmrl$GeneL),unlist(dmr4$GeneL))

ups <- union(up_h_c4,down_h_c_low)
downs <- union(down_h_c4,up_h_c_low)
merge <- ups[ups %in% downs]
selected_ups <- ups[!ups %in% merge]
selected_downs <- downs[!downs %in% merge]
ups <- selected_ups[selected_ups %in% uni_dmr]
downs <- selected_downs[selected_downs %in% uni_dmr]

mylist <- up_h_c3[!up_h_c3 %in% merge]
mylist <- mylist[mylist %in% uni_dmr]
mylist <- mylist[mylist %in% up_h_c_low]

dmrl[grep("CDK6",dmr4$GeneL),]

ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes of methylation regulated", RegulonList = ups)
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = all_genes_met_h,nBar = 20,filename = "DEA_ups.pdf")
ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes of methylation regulated", RegulonList = downs)
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = all_genes_met_h,nBar = 20,filename = "DEA_downs.pdf")
ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes of methylation regulated", RegulonList = union(ups,downs))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = all_genes_met_h,nBar = 20,filename = "DEA_union.pdf")

g <- "PCDHGA2"
g %in% ups
g %in% up_h_c4
g %in% down_h_c_low

for(g in finallist[91:91]){
  g <- "SLMO2"
  print(strrep("*",100))
  print(g)
  print(dmrl[grep(g,dmrl$GeneL),c("status","Position","Name","Gene")])
  print(strrep("-",50))
  print(dmr4[grep(g,dmr4$GeneL),c("status","Position","Name","Gene")])
}


for(g in ups[ups %in% unlist(dmr4$GeneL)]){
  print(strrep("*",100))
  print(g)
  print(dmrl[grep(g,dmrl$GeneL),c("status","Position","Name","Gene")])
  print(strrep("-",50))
  print(dmr4[grep(g,dmr4$GeneL),c("status","Position","Name","Gene")])
}

for(g in downs[downs %in% unlist(dmr4$GeneL)]){
  print(strrep("*",100))
  print(g)
  print(dmrl[grep(g,dmrl$GeneL),c("status","Position","Name","Gene")])
  print(strrep("-",50))
  print(dmr4[grep(g,dmr4$GeneL),c("status","Position","Name","Gene")])
}

sort(unique(unlist(strsplit(capture.output(cat(values(hyper)[,"Gene_Symbol"],sep=";")),";"))))

grep(up_h_c_low_c4,dmr3$Gene)

table(met.dmr.h.n.c.low$status)

table(met.dmr.h.n.c1$status)
table(met.dmr.h.n.c1$status)
table(met.dmr.h.n.c2$status)
table(met.dmr.h.n.c3$status)
table(met.dmr.h.n.c4$status)


save(met_h_n,met_h_c,met_h_c_n,met_h_c_1,met_h_c_2,met_h_c_3,met_h_c_4,met_h_c_s,
     met.dmr.h.n.c1,met.dmr.h.n.c2,met.dmr.h.n.c3,met.dmr.h.n.c4,met.dmr.h.low.c4,file = papaste0(p_h,"meth.RData"))
save(genelistDEGs_n_c,genelistDEGs_n_c1,genelistDEGs_n_c2,genelistDEGs_n_c3,genelistDEGs_n_c4,genelistDEGs_n_c_low,
     genelistDEGs_c1_c2,genelistDEGs_c1_c3,genelistDEGs_c1_c4,genelistDEGs_c2_c3,genelistDEGs_c2_c4,genelistDEGs_c3_c4,
     genelistDEGs_c_low_c4,file=papaste0(p_h,"genelist.RData"))
save(ansEA_h_n_c,ansEA_h_n_c1,ansEA_h_n_c2,ansEA_h_n_c3,ansEA_h_n_c4,ansEA_h_n_c_low,ansEA_h_c1_c2,ansEA_h_c1_c3,ansEA_h_c1_c4,
     ansEA_h_c2_c3,ansEA_h_c2_c4,ansEA_h_c3_c4,ansEA_h_c_low_c4,file = papaste0(p_h,"ansEA.RData"))
save(dataDEGs_h_n_c,dataDEGs_h_n_c1,dataDEGs_h_n_c2,dataDEGs_h_n_c3,dataDEGs_h_n_c4,dataDEGs_h_n_c_low,dataDEGs_h_c1_c2,
     dataDEGs_h_c1_c3,dataDEGs_h_c1_c4,dataDEGs_h_c2_c3,dataDEGs_h_c2_c4,dataDEGs_h_c3_c4,dataDEGs_h_c_low_c4, file=papaste0(p_h,"dataDEG.RData"))
save(rse_h_n,rse_h_c,rse_h_c_1,rse_h_c_2,rse_h_c_3,rse_h_c_4,rse_h_c_n,rse_h_c_s,rse_h_c_low, file=papaste0(p_h,"rse.RData"))


save(rse_h_c,met_h_c,rse_h_n,met_h_n, file=paste0(p_h,"_datafile_PCC.RData"))
load("df.RData")


rowData(met_h_c)$Gene_Symbol

up_g1 <- rownames(dataDEGs_n_g1[dataDEGs_n_g1$logFC>0,])
downs_g1 <- rownames(dataDEGs_n_g1[dataDEGs_n_g1$logFC<0,])
up_g2 <- rownames(dataDEGs_n_g2[dataDEGs_n_g2$logFC>0,])
downs_g2 <- rownames(dataDEGs_n_g2[dataDEGs_n_g2$logFC<0,])

mygenes <- union(union(up_g1,downs_g1),union(up_g2,downs_g2))

#my <- rowData(met_h_c)[(rowData(met_h_c)$Gene_Symbol %in% union(ups,downs)) & (rowData(met_h_c)$probeID %in% rownames(met.dmr.h.n.c[met.dmr.h.n.c$status!="Not Significant",])),]

#mygenes <- union(up_g2[!up_g2 %in% up_g1],downs_g2[!downs_g2 %in% downs_g1])
myprobes <- union(union(hyper_probes_g1,hyper_probes_g2),union(hypo_probes_g1,hypo_probes_g2))

my <- rowData(met_h_c)[(rowData(met_h_c)$Gene_Symbol %in% mygenes) & (rowData(met_h_c)$probeID %in% myprobes),]

############# 
############# Calc PCC

print(nrow(my))
c <-0
probe <- c()
gene <- c()
genename <- c()
pvalue <- c()
pcctest <- c()
pos <- c()
chr <- c()
rel <- c()
gbo <- c()

for(i in 0:nrow(my)){
  if(!any(is.na(assay(met_h_c)[my[i,1],]))){
    m <- Winsorize(assay(met_h_c)[my[i,1],])
    r <- Winsorize(assays(rse_h_c)$raw_count[rowData(rse_h_c)$gene_id==my[i,2],])
    
    #scaled_estimate
    #raw_count
    
    names(m) <- substr(names(m),1,12)
    names(r) <- substr(names(r),1,12)
    
    r <- r[names(r) %in% names(m)]
    m <- m[names(m) %in% names(r)]
    
    r <- r[order(names(r))]
    m <- m[order(names(m))]
    
    if(length(m)>0){
      R <- cor(r, m, method = "pearson")
      P <- cor.test(r, m, method= "pearson")
      
      if(!isNA(R) & P$p.value<0.05 & abs(R)>0.5){
        if(!isNA(met.dmr.g1[my[i,1],]$status) & (met.dmr.g1[my[i,1],]$status != met.dmr.g2[my[i,1],]$status))
          if(xor(my[i,2] %in% downs_g1,my[i,2] %in% downs_g2) | xor(my[i,2] %in% up_g1,my[i,2] %in% up_g2)){
            pcctest <- c(pcctest,R)
            pvalue <- c(pvalue,P$p.value)
            gene <- c(gene,my[i,2])
            probe <- c(probe,my[i,1])
            rel <- c(rel,ann450k[my[i,1],]$Relation_to_Island)
            pos <- c(pos,ann450k[my[i,1],]$pos)
            chr <- c(chr,ann450k[my[i,1],]$chr)
            gbo <- c(gbo,ann450k[my[i,1],]$UCSC_RefGene_Group)
            genename <- c(genename,ann450k[my[i,1],]$UCSC_RefGene_Name)
            c <-  c+1
          }
        #      print(i)
        #      print(my[i,])
        #      l <- c(l,my[i,2])
        #      print(P$p.value)
        #      print(P$estimate)
        #      print('--------------------------')
      }
    }
    
    #  print(my[i,1],my[i,2],R,P)
  }
}

print(c)

View(data.frame(probe,gbo,genename,rel,pos,gene,chr,pvalue,pcctest,dataDEGs_n_g1[gene,],dataDEGs_n_g2[gene,],met.dmr.g1[probe,]$status,met.dmr.g2[probe,]$status))

View(data.frame(chr,gene,probe,gbo,rel,met.dmr.g1[probe,]$status,met.dmr.g2[probe,]$status,dataDEGs_n_g1[gene,]$logFC,dataDEGs_n_g2[gene,]$logFC,dataDEGs_n_g1[gene,]$PValue,dataDEGs_n_g2[gene,]$PValue,met.dmr.g1[probe,]$p.value.Solid.Tissue.Normal.Primary.Tumor,met.dmr.g2[probe,]$p.value.Solid.Tissue.Normal.Primary.Tumor,pcctest,pvalue))

############# 
############# Report

mygene <- "IRX4"
myprobe <- "cg00409356"

dataDEGs_n_g1[gene,]
dataDEGs_n_g2[gene,]
met.dmr.g1[probe,]$status
met.dmr.g2[probe,]$status

data.frame(dataDEGs_n_g1[gene,],dataDEGs_n_g2[gene,],met.dmr.g1[probe,]$status,met.dmr.g2[probe,]$status)
############# 

calcPCC <- function(probe,gene){
  m <- Winsorize(assay(met_h_c_4)[probe,])
  r <- Winsorize(assays(rse_h_c_4)$raw_count[rowData(rse_h_c_4)$gene_id==gene,])
  
  names(m) <- substr(names(m),1,12)
  names(r) <- substr(names(r),1,12)
  
  r <- r[names(r) %in% names(m)]
  m <- m[names(m) %in% names(r)]
  
  r <- r[order(names(r))]
  m <- m[order(names(m))]
  
  R <- cor(r, m, method = "pearson")
  P <- cor.test(r, m, method= "pearson")
  
  print(probe)
  print(gene)
  print(P$p.value)
  print(P$estimate)
}

calcPCC("cg04548096","EYA4")

print(c)

mydf <- data.frame(r,m)

ggscatter(mydf, x = "r", y = "m", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "RSE", ylab = "Meth")

ggscatter(mydf, x = "r", y = "m",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)


ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes Normal vs ESCA"), RegulonList = rownames(dataDEGs_h_n_c))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = rownames(dataDEGs_h_n_c),
                        nBar = 10,filename = paste0("DEA_genes_ESCA.pdf"))

ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes Normal vs ESCA early stages"), RegulonList = rownames(dataDEGs_n_g1))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = rownames(dataDEGs_n_g1),
                        nBar = 10,filename = paste0("DEA_genes_G1.pdf"))

ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes Normal vs ESCA metastasis samples"), RegulonList = rownames(dataDEGs_n_g2))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = rownames(dataDEGs_n_g2),
                        nBar = 10,filename = paste0("DEA_genes_G2.pdf"))



probes <- rownames(met.dmr.g1[met.dmr.g1$status=="Hypomethylated in Solid Tissue Normal",])
probes <- union(probes,rownames(met.dmr.g2[met.dmr.g2$status=="Hypomethylated in Solid Tissue Normal",]))
probes <- union(probes,rownames(met.dmr.h.n.c[met.dmr.h.n.c$status=="Hypomethylated in Solid Tissue Normal",]))
probes <- union(probes,rownames(met.dmr.h.n.c.low[met.dmr.h.n.c.low$status=="Hypomethylated in Solid Tissue Normal",]))
probes <- union(probes,rownames(met.dmr.h.n.c1[met.dmr.h.n.c1$status=="Hypomethylated in Solid Tissue Normal",]))
probes <- union(probes,rownames(met.dmr.h.n.c2[met.dmr.h.n.c2$status=="Hypomethylated in Solid Tissue Normal",]))
probes <- union(probes,rownames(met.dmr.h.n.c3[met.dmr.h.n.c3$status=="Hypomethylated in Solid Tissue Normal",]))
probes <- union(probes,rownames(met.dmr.h.n.c4[met.dmr.h.n.c4$status=="Hypomethylated in Solid Tissue Normal",]))
probes_hypo <- probes
length(probes_hypo)

probes <- rownames(met.dmr.g1[met.dmr.g1$status=="Hypermethylated in Solid Tissue Normal",])
probes <- union(probes,rownames(met.dmr.g2[met.dmr.g2$status=="Hypermethylated in Solid Tissue Normal",]))
probes <- union(probes,rownames(met.dmr.h.n.c[met.dmr.h.n.c$status=="Hypermethylated in Solid Tissue Normal",]))
probes <- union(probes,rownames(met.dmr.h.n.c.low[met.dmr.h.n.c.low$status=="Hypermethylated in Solid Tissue Normal",]))
probes <- union(probes,rownames(met.dmr.h.n.c1[met.dmr.h.n.c1$status=="Hypermethylated in Solid Tissue Normal",]))
probes <- union(probes,rownames(met.dmr.h.n.c2[met.dmr.h.n.c2$status=="Hypermethylated in Solid Tissue Normal",]))
probes <- union(probes,rownames(met.dmr.h.n.c3[met.dmr.h.n.c3$status=="Hypermethylated in Solid Tissue Normal",]))
probes <- union(probes,rownames(met.dmr.h.n.c4[met.dmr.h.n.c4$status=="Hypermethylated in Solid Tissue Normal",]))
probes_hyper <- probes

table(ann450k[probes_hypo,]$Relation_to_Island)
table(ann450k[probes_hyper,]$Relation_to_Island)
table(unlist(strsplit(ann450k[probes_hypo,]$UCSC_RefGene_Group,";")))
table(unlist(strsplit(ann450k[probes_hyper,]$UCSC_RefGene_Group,";")))
