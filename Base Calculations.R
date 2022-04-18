################# Methylation

m_query_h_n <- GDCquery(project = p_h, legacy = TRUE,data.category = "DNA methylation",
                        platform = "Illumina Human Methylation 450",sample.type = t_n)

GDCdownload(m_query_h_n)
met_h_n <- GDCprepare(m_query_h_n,summarizedExperiment = TRUE)

met_h_n <- met_h_n[,substr(colnames(met_h_n),1,12) %in% dataClin_h$bcr_patient_barcode]

m_query_h_c <- GDCquery(project = p_h, legacy = TRUE,data.category = "DNA methylation",
                        platform = "Illumina Human Methylation 450",sample.type = t_c)

GDCdownload(m_query_h_c)
met_h_c <- GDCprepare(m_query_h_c,summarizedExperiment = TRUE)

met_h_c_org <- met_h_c
met_h_n_org <- met_h_n

colData(met_h_n) <- colData(met_h_n)[names(colData(met_h_n)[c("barcode","tumor_stage","sample_type")])]

colData(met_h_c) <- colData(met_h_c)[names(colData(met_h_c)[c("barcode","tumor_stage","sample_type")])]

met_h_c <- met_h_c[,substr(colnames(met_h_c),1,12) %in% dataClin_h$bcr_patient_barcode]

met_h_c_n <- met_h_c[,substr(colnames(met_h_c),1,12) %in% dataClin_h_n$bcr_patient_barcode]
met_h_c_1 <- met_h_c[,substr(colnames(met_h_c),1,12) %in% dataClin_h_1$bcr_patient_barcode]
met_h_c_2 <- met_h_c[,substr(colnames(met_h_c),1,12) %in% dataClin_h_2$bcr_patient_barcode]
met_h_c_3 <- met_h_c[,substr(colnames(met_h_c),1,12) %in% dataClin_h_3$bcr_patient_barcode]
met_h_c_4 <- met_h_c[,substr(colnames(met_h_c),1,12) %in% dataClin_h_4$bcr_patient_barcode]
met_h_c_low <- met_h_c[,substr(colnames(met_h_c),1,12) %in% dataClin_h_low$bcr_patient_barcode]
met_h_c_s <- met_h_c[,substr(colnames(met_h_c),1,12) %in% dataClin_h_s$bcr_patient_barcode]

met_G1 <- met_h_c[,substr(colnames(met_h_c),1,12) %in% dataClin_G1$bcr_patient_barcode]
met_G2 <- met_h_c[,substr(colnames(met_h_c),1,12) %in% dataClin_G2$bcr_patient_barcode]

tic("DMR...")
tic("DMR N C...")
met <- SummarizedExperiment::cbind(met_h_n, met_h_c)
met <- subset(met,subset = (rowSums(is.na(assay(met))) == 0))
met <- subset(met,subset = !as.character(seqnames(met)) %in% c("chrNA","chrX","chrY"))
met.dmr.h.n.c <- TCGAanalyze_DMC(met, groupCol = "sample_type", # a column in the colData matrix
                                 p.cut = 0.01, diffmean.cut = 0.25, legend = "State",
                                 plot.filename = paste0(p_h,"_N_C_metvolcano.png"), cores = 1 # if set to 1 there will be a progress bar
)

table(met.dmr.h.n.c$status)

toc()

met <- SummarizedExperiment::cbind(met_h_n, met_h_c_1)
met <- subset(met,subset = (rowSums(is.na(assay(met))) == 0))
met <- subset(met,subset = !as.character(seqnames(met)) %in% c("chrNA","chrX","chrY"))
met.dmr.h.n.c1 <- TCGAanalyze_DMC(met, groupCol = "sample_type", # a column in the colData matrix
                                  p.cut = 0.01, diffmean.cut = 0.25, legend = "State",
                                  plot.filename = paste0(p_h,"_N_C1_metvolcano.png"), cores = 1 # if set to 1 there will be a progress bar
)
met <- SummarizedExperiment::cbind(met_h_n, met_h_c_2)
met <- subset(met,subset = (rowSums(is.na(assay(met))) == 0))
met <- subset(met,subset = !as.character(seqnames(met)) %in% c("chrNA","chrX","chrY"))
met.dmr.h.n.c2 <- TCGAanalyze_DMC(met, groupCol = "sample_type", # a column in the colData matrix
                                  p.cut = 0.01, diffmean.cut = 0.25, legend = "State",
                                  plot.filename = paste0(p_h,"_N_C2_metvolcano.png"), cores = 1 # if set to 1 there will be a progress bar
)
met <- SummarizedExperiment::cbind(met_h_n, met_h_c_3)
met <- subset(met,subset = (rowSums(is.na(assay(met))) == 0))
met <- subset(met,subset = !as.character(seqnames(met)) %in% c("chrNA","chrX","chrY"))
met.dmr.h.n.c3 <- TCGAanalyze_DMC(met, groupCol = "sample_type", # a column in the colData matrix
                                  p.cut = 0.01, diffmean.cut = 0.25, legend = "State",
                                  plot.filename = paste0(p_h,"_N_C3_metvolcano.png"), cores = 1 # if set to 1 there will be a progress bar
)
met <- SummarizedExperiment::cbind(met_h_n, met_h_c_4)
met <- subset(met,subset = (rowSums(is.na(assay(met))) == 0))
met <- subset(met,subset = !as.character(seqnames(met)) %in% c("chrNA","chrX","chrY"))
met.dmr.h.n.c4 <- TCGAanalyze_DMC(met, groupCol = "sample_type", # a column in the colData matrix
                                  p.cut = 0.01, diffmean.cut = 0.25, legend = "State",
                                  plot.filename = paste0(p_h,"_N_C4_metvolcano.png"), cores = 1 # if set to 1 there will be a progress bar
)
met <- SummarizedExperiment::cbind(met_h_n, met_h_c_low)
met <- subset(met,subset = (rowSums(is.na(assay(met))) == 0))
met <- subset(met,subset = !as.character(seqnames(met)) %in% c("chrNA","chrX","chrY"))
met.dmr.h.n.c.low <- TCGAanalyze_DMC(met, groupCol = "sample_type", # a column in the colData matrix
                                     p.cut = 0.01, diffmean.cut = 0.25, legend = "State",
                                     plot.filename = paste0(p_h,"_N_C_low_metvolcano.png"), cores = 1 # if set to 1 there will be a progress bar
)

met <- SummarizedExperiment::cbind(met_h_n, met_G1)
met <- subset(met,subset = (rowSums(is.na(assay(met))) == 0))
met <- subset(met,subset = !as.character(seqnames(met)) %in% c("chrNA","chrX","chrY"))
met.dmr.g1 <- TCGAanalyze_DMC(met, groupCol = "sample_type", # a column in the colData matrix
                              p.cut = 0.01, diffmean.cut = 0.25, legend = "State",
                              plot.filename = paste0(p_h,"_N_G1_metvolcano.png"), cores = 1 # if set to 1 there will be a progress bar
)

met <- SummarizedExperiment::cbind(met_h_n, met_G2)
met <- subset(met,subset = (rowSums(is.na(assay(met))) == 0))
met <- subset(met,subset = !as.character(seqnames(met)) %in% c("chrNA","chrX","chrY"))
met.dmr.g2 <- TCGAanalyze_DMC(met, groupCol = "sample_type", # a column in the colData matrix
                              p.cut = 0.01, diffmean.cut = 0.25, legend = "State",
                              plot.filename = paste0(p_h,"_N_G2_metvolcano.png"), cores = 1 # if set to 1 there will be a progress bar
)


met <- met_h_c_1
table(met$tumor_stage)
met$tumor_stage <- "stage i"
met_h_1 <- met
met <- met_h_c_2
table(met$tumor_stage)
met$tumor_stage <- "stage ii"
met_h_2 <- met

met <- SummarizedExperiment::cbind(met_h_1, met_h_2)
met <- subset(met,subset = (rowSums(is.na(assay(met))) == 0))
met <- subset(met,subset = !as.character(seqnames(met)) %in% c("chrNA","chrX","chrY"))
met.dmr.h.c1.c2 <- TCGAanalyze_DMC(met, groupCol = "tumor_stage", # a column in the colData matrix
                                   p.cut = 0.01, diffmean.cut = 0.25, legend = "State",
                                   plot.filename = paste0(p_h,"_C1_C2_metvolcano.png"), cores = 1 # if set to 1 there will be a progress bar
)

met <- met_h_c_3
table(met$tumor_stage)
met$tumor_stage <- "stage iii"
met_h_3 <- met
met <- SummarizedExperiment::cbind(met_h_2, met_h_3)
met <- subset(met,subset = (rowSums(is.na(assay(met))) == 0))
met <- subset(met,subset = !as.character(seqnames(met)) %in% c("chrNA","chrX","chrY"))
met.dmr.h.c2.c3 <- TCGAanalyze_DMC(met, groupCol = "tumor_stage", # a column in the colData matrix
                                   p.cut = 0.01, diffmean.cut = 0.25, legend = "State",
                                   plot.filename = paste0(p_h,"_C2_C3_metvolcano.png"), cores = 1 # if set to 1 there will be a progress bar
)
met <- met_h_c_4
table(met$tumor_stage)
met$tumor_stage <- "stage iv"
met_h_4 <- met
met <- SummarizedExperiment::cbind(met_h_3, met_h_4)
met <- subset(met,subset = (rowSums(is.na(assay(met))) == 0))
met <- subset(met,subset = !as.character(seqnames(met)) %in% c("chrNA","chrX","chrY"))
met.dmr.h.c3.c4 <- TCGAanalyze_DMC(met, groupCol = "tumor_stage", # a column in the colData matrix
                                   p.cut = 0.01, diffmean.cut = 0.25, legend = "State",
                                   plot.filename = paste0(p_h,"_C3_C4_metvolcano.png"), cores = 1 # if set to 1 there will be a progress bar
)

met <- met_h_c_low
table(met$tumor_stage)
met$tumor_stage <- "stage low"
met_h_low <- met
met <- SummarizedExperiment::cbind(met_h_low, met_h_4)
met <- subset(met,subset = (rowSums(is.na(assay(met))) == 0))
met <- subset(met,subset = !as.character(seqnames(met)) %in% c("chrNA","chrX","chrY"))
met.dmr.h.low.c4 <- TCGAanalyze_DMC(met, groupCol = "tumor_stage", # a column in the colData matrix
                                    p.cut = 0.01, diffmean.cut = 0.25, legend = "State",
                                    plot.filename = paste0(p_h,"_C_low_C4_metvolcano.png"), cores = 1 # if set to 1 there will be a progress bar
)


table(met.dmr.h.n.c1[,"status"])
table(met.dmr.h.n.c2[,"status"])
table(met.dmr.h.n.c3[,"status"])
table(met.dmr.h.n.c4[,"status"])
table(met.dmr.h.n.c[,"status"])
table(met.dmr.h.n.c.low[,"status"])

row.names(met.dmr.h.n.c.low[met.dmr.h.n.c.low$status=="Hypermethylated in Solid Tissue Normal",])
row.names(met.dmr.h.n.c.low[met.dmr.h.n.c.low$status=="Hypomethylated in Solid Tissue Normal",])

row.names(met.dmr.h.n.c4[met.dmr.h.n.c4$status=="Hypermethylated in Solid Tissue Normal",])
row.names(met.dmr.h.n.c4[met.dmr.h.n.c4$status=="Hypomethylated in Solid Tissue Normal",])


met <- met.dmr.h.n.c1
rows <- rownames(met[met$status=="Hypermethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_n)[rowData(met_h_n)$probeID %in% rows,]$Gene_Symbol
hyper_h_n_c1 <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))
rows <- rownames(met[met$status=="Hypomethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_n)[rowData(met_h_n)$probeID %in% rows,]$Gene_Symbol
hypo_h_n_c1 <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))

met <- met.dmr.h.n.c2
rows <- rownames(met[met$status=="Hypermethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_n)[rowData(met_h_n)$probeID %in% rows,]$Gene_Symbol
hyper_h_n_c2 <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))
rows <- rownames(met[met$status=="Hypomethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_n)[rowData(met_h_n)$probeID %in% rows,]$Gene_Symbol
hypo_h_n_c2 <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))

met <- met.dmr.h.n.c3
rows <- rownames(met[met$status=="Hypermethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_n)[rowData(met_h_n)$probeID %in% rows,]$Gene_Symbol
hyper_h_n_c3 <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))
rows <- rownames(met[met$status=="Hypomethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_n)[rowData(met_h_n)$probeID %in% rows,]$Gene_Symbol
hypo_h_n_c3 <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))

met <- met.dmr.h.n.c4
rows <- rownames(met[met$status=="Hypermethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_n)[rowData(met_h_n)$probeID %in% rows,]$Gene_Symbol
hyper_h_n_c4 <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))
rows <- rownames(met[met$status=="Hypomethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_n)[rowData(met_h_n)$probeID %in% rows,]$Gene_Symbol
hypo_h_n_c4 <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))

met <- met.dmr.h.n.c.low
rows <- rownames(met[met$status=="Hypermethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_n)[rowData(met_h_n)$probeID %in% rows,]$Gene_Symbol
hyper_h_n_c_low <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))
rows <- rownames(met[met$status=="Hypomethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_n)[rowData(met_h_n)$probeID %in% rows,]$Gene_Symbol
hypo_h_n_c_low <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))

met <- met.dmr.g1
rows <- rownames(met[met$status=="Hypermethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_n)[rowData(met_h_n)$probeID %in% rows,]$Gene_Symbol
hyper_g1 <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))
rows <- rownames(met[met$status=="Hypomethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_n)[rowData(met_h_n)$probeID %in% rows,]$Gene_Symbol
hypo_g1 <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))

met <- met.dmr.g2
rows <- rownames(met[met$status=="Hypermethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_n)[rowData(met_h_n)$probeID %in% rows,]$Gene_Symbol
hyper_g2 <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))
rows <- rownames(met[met$status=="Hypomethylated in Solid Tissue Normal",])
mylist <- rowData(met_h_n)[rowData(met_h_n)$probeID %in% rows,]$Gene_Symbol
hypo_g2 <- sort(unique(unlist(strsplit(capture.output(cat(mylist,sep=";")),";"))))

hyper_probes_g1 <- rownames(met.dmr.g1[met.dmr.g1$status=="Hypermethylated in Solid Tissue Normal",])
hypo_probes_g1 <- rownames(met.dmr.g1[met.dmr.g1$status=="Hypomethylated in Solid Tissue Normal",])
hyper_probes_g2 <- rownames(met.dmr.g2[met.dmr.g2$status=="Hypermethylated in Solid Tissue Normal",])
hypo_probes_g2 <- rownames(met.dmr.g2[met.dmr.g2$status=="Hypomethylated in Solid Tissue Normal",])


hyper_probes_c4 <- rownames(met.dmr.h.n.c4[met.dmr.h.n.c4$status=="Hypermethylated in Solid Tissue Normal",])
hypo_probes_c4 <- rownames(met.dmr.h.n.c4[met.dmr.h.n.c4$status=="Hypomethylated in Solid Tissue Normal",])
hyper_probes_low <- rownames(met.dmr.h.n.c.low[met.dmr.h.n.c.low$status=="Hypermethylated in Solid Tissue Normal",])
hypo_probes_low <- rownames(met.dmr.h.n.c.low[met.dmr.h.n.c.low$status=="Hypomethylated in Solid Tissue Normal",])

all_genes_hyper_h <- union(union(hyper_h_n_c1,hyper_h_n_c2),union(hyper_h_n_c3,hyper_h_n_c4))
all_genes_hypo_h <- union(union(hypo_h_n_c1,hypo_h_n_c2),union(hypo_h_n_c3,hypo_h_n_c4))
all_genes_met_h <- union(all_genes_hyper_h,all_genes_hypo_h)

hypo_h_n_c_low[!hypo_h_n_c_low %in% hyper_h_n_c4]
hyper_h_n_c_low[!hyper_h_n_c_low %in% all_genes_hyper_h]
hypo_h_n_c3[!hypo_h_n_c3 %in% hypo_h_n_c_low]
hyper_h_n_c3[!hyper_h_n_c3 %in% hyper_h_n_c_low]

hyper_nm <- union(union(hyper_h_n_c1,hyper_h_n_c2),union(hyper_h_n_c3,hyper_h_n_c_low))
hypo_nm <- union(union(hypo_h_n_c1,hypo_h_n_c2),union(hypo_h_n_c3,hypo_h_n_c_low))

#hyper <- met.dmr.h.n.c1[values(met.dmr.h.n.c1)[,"status.Primary.Tumor.Solid.Tissue.Normal"]=="Hypermethylated"]
#hypo <- met.dmr.h.n.c1[values(met.dmr.h.n.c1)[,"status.Primary.Tumor.Solid.Tissue.Normal"]=="Hypomethylated"]

#hyper_h_n_c1 <- sort(unique(unlist(strsplit(capture.output(cat(values(met.dmr.h.n.c1[values(met.dmr.h.n.c1)[,"status.Primary.Tumor.Solid.Tissue.Normal"]=="Hypermethylated"])[,"Gene_Symbol"],sep=";")),";"))))
#hyper_h_n_c2 <- sort(unique(unlist(strsplit(capture.output(cat(values(met.dmr.h.n.c2[values(met.dmr.h.n.c2)[,"status.Primary.Tumor.Solid.Tissue.Normal"]=="Hypermethylated"])[,"Gene_Symbol"],sep=";")),";"))))
#hyper_h_n_c3 <- sort(unique(unlist(strsplit(capture.output(cat(values(met.dmr.h.n.c3[values(met.dmr.h.n.c3)[,"status.Primary.Tumor.Solid.Tissue.Normal"]=="Hypermethylated"])[,"Gene_Symbol"],sep=";")),";"))))
#hyper_h_n_c4 <- sort(unique(unlist(strsplit(capture.output(cat(values(met.dmr.h.n.c4[values(met.dmr.h.n.c4)[,"status.Primary.Tumor.Solid.Tissue.Normal"]=="Hypermethylated"])[,"Gene_Symbol"],sep=";")),";"))))
#hyper_h_n_c_low <- sort(unique(unlist(strsplit(capture.output(cat(values(met.dmr.h.n.c.low[values(met.dmr.h.n.c.low)[,"status.Primary.Tumor.Solid.Tissue.Normal"]=="Hypermethylated"])[,"Gene_Symbol"],sep=";")),";"))))

#hypo_h_n_c1 <- sort(unique(unlist(strsplit(capture.output(cat(values(met.dmr.h.n.c1[values(met.dmr.h.n.c1)[,"status.Primary.Tumor.Solid.Tissue.Normal"]=="Hypomethylated"])[,"Gene_Symbol"],sep=";")),";"))))
#hypo_h_n_c2 <- sort(unique(unlist(strsplit(capture.output(cat(values(met.dmr.h.n.c2[values(met.dmr.h.n.c2)[,"status.Primary.Tumor.Solid.Tissue.Normal"]=="Hypomethylated"])[,"Gene_Symbol"],sep=";")),";"))))
#hypo_h_n_c3 <- sort(unique(unlist(strsplit(capture.output(cat(values(met.dmr.h.n.c3[values(met.dmr.h.n.c3)[,"status.Primary.Tumor.Solid.Tissue.Normal"]=="Hypomethylated"])[,"Gene_Symbol"],sep=";")),";"))))
#hypo_h_n_c4 <- sort(unique(unlist(strsplit(capture.output(cat(values(met.dmr.h.n.c4[values(met.dmr.h.n.c4)[,"status.Primary.Tumor.Solid.Tissue.Normal"]=="Hypomethylated"])[,"Gene_Symbol"],sep=";")),";"))))
#hypo_h_n_c.low <- sort(unique(unlist(strsplit(capture.output(cat(values(met.dmr.h.n.c.low[values(met.dmr.h.n.c.low)[,"status.Primary.Tumor.Solid.Tissue.Normal"]=="Hypomethylated"])[,"Gene_Symbol"],sep=";")),";"))))

#all_genes_hyper_h <- union(union(hyper_h_n_c1,hyper_h_n_c2),union(hyper_h_n_c3,hyper_h_n_c4))
#all_genes_hypo_h <- union(union(hypo_h_n_c1,hypo_h_n_c2),union(hypo_h_n_c3,hypo_h_n_c4))
#all_genes_met_h <- union(all_genes_hyper_h,all_genes_hypo_h)

sig.met <- met.dmr[values(met.dmr)[,"status.stage.ii.stage.i"] %in% c("Hypermethylated","Hypomethylated"),]
dmr_genes <- unique(values(sig.met)["Gene_Symbol"])
toc()

save(met_h_c_org,met_h_n_org,met_h_c,met_h_n,
     met.dmr.h.n.c,met.dmr.h.n.c1,met.dmr.h.n.c2,met.dmr.h.n.c3,met.dmr.h.n.c4,met.dmr.h.n.c.low,file = paste0(p_h,"_meth.RData"))


################



############################# RSE

query_h_n <- GDCquery(project = p_h, 
                      legacy = TRUE,data.category = "Gene expression",data.type = "Gene expression quantification",
                      file.type = "results",platform = "Illumina HiSeq",sample.type = t_n)

query_h_c <- GDCquery(project = p_h, 
                      legacy = TRUE,data.category = "Gene expression",data.type = "Gene expression quantification",
                      file.type = "results",platform = "Illumina HiSeq",sample.type = t_c)

GDCdownload(query_h_n)
rse_h_n <- GDCprepare(query_h_n,summarizedExperiment = TRUE)
#rse_h_n <- rse_h_n[rowSums(assay(rse_h_n))>5,]

GDCdownload(query_h_c)
rse_h_c <- GDCprepare(query_h_c,summarizedExperiment = TRUE)
#rse_h_c <- rse_h_c[rowSums(assay(rse_h_c))>5,]

rse_h_c_org <- rse_h_c
rse_h_n_org <- rse_h_n


newrse <- rse_h_c
rc <- apply(assays(newrse)$raw_count,1,Winsorize)
assays(rse_h_c)$raw_count <- t(rc)
newrse <- rse_h_n
rc <- apply(assays(newrse)$raw_count,1,Winsorize)
assays(rse_h_n)$raw_count <- t(rc)

rse_h_n <- rse_h_n[,substr(colnames(rse_h_n),1,12) %in% dataClin_h$bcr_patient_barcode]

rse_h_c_n <- rse_h_c[,substr(colnames(rse_h_c),1,12) %in% dataClin_h_n$bcr_patient_barcode]
rse_h_c_1 <- rse_h_c[,substr(colnames(rse_h_c),1,12) %in% dataClin_h_1$bcr_patient_barcode]
rse_h_c_2 <- rse_h_c[,substr(colnames(rse_h_c),1,12) %in% dataClin_h_2$bcr_patient_barcode]
rse_h_c_3 <- rse_h_c[,substr(colnames(rse_h_c),1,12) %in% dataClin_h_3$bcr_patient_barcode]
rse_h_c_4 <- rse_h_c[,substr(colnames(rse_h_c),1,12) %in% dataClin_h_4$bcr_patient_barcode]
rse_h_c_low <- rse_h_c[,substr(colnames(rse_h_c),1,12) %in% dataClin_h_low$bcr_patient_barcode]
rse_h_c_s <- rse_h_c[,!(substr(colnames(rse_h_c),1,12) %in% dataClin_h_n$bcr_patient_barcode)]

rse_G1 <- rse_h_c[,substr(colnames(rse_h_c),1,12) %in% dataClin_G1$bcr_patient_barcode]
rse_G2 <- rse_h_c[,substr(colnames(rse_h_c),1,12) %in% dataClin_G2$bcr_patient_barcode]

rse_h_c_n <- rse_h_c_n[,substr(colnames(rse_h_c_n),1,12) %in% substr(colnames(met_h_c_n),1,12)]
rse_h_c_1 <- rse_h_c_1[,substr(colnames(rse_h_c_1),1,12) %in% substr(colnames(met_h_c_1),1,12)]
rse_h_c_2 <- rse_h_c_2[,substr(colnames(rse_h_c_2),1,12) %in% substr(colnames(met_h_c_2),1,12)]
rse_h_c_3 <- rse_h_c_3[,substr(colnames(rse_h_c_3),1,12) %in% substr(colnames(met_h_c_3),1,12)]
rse_h_c_4 <- rse_h_c_4[,substr(colnames(rse_h_c_4),1,12) %in% substr(colnames(met_h_c_4),1,12)]
rse_h_c_low <- rse_h_c_low[,substr(colnames(rse_h_c_low),1,12) %in% substr(colnames(met_h_c_low),1,12)]
rse_h_c_s <- rse_h_c_s[,substr(colnames(rse_h_c_s),1,12) %in% substr(colnames(met_h_c_s),1,12)]

rse_G1 <- rse_G1[,substr(colnames(rse_G1),1,12) %in% substr(colnames(met_G1),1,12)]
rse_G1 <- rse_G1[,substr(colnames(rse_G2),1,12) %in% substr(colnames(met_G2),1,12)]

################# CNV
############
tic("CNV ...")
tic("CNV Prepare...")
query.ESCA.nocnv <- GDCquery(project = p_h,data.category = "Copy number variation",
                             legacy = TRUE, file.type = "nocnv_hg19.seg")

GDCdownload(query.ESCA.nocnv)

cnvMatrix <- GDCprepare(query.ESCA.nocnv)#, save = TRUE, save.filename = "ESCAnocnvhg19.rda")

# Add label (0 for loss, 1 for gain)
cnvMatrix <- cbind(cnvMatrix,Label=NA)
cnvMatrix[cnvMatrix[, "Segment_Mean"] < -0.3, "Label"] <- 0
cnvMatrix[cnvMatrix[, "Segment_Mean"] > 0.3, "Label"] <- 1
cnvMatrix <- cnvMatrix [!is.na(cnvMatrix$Label),]

cnvMatrix <- cnvMatrix[,-6]
colnames(cnvMatrix) <- c("Sample.Name", "Chromosome", "Start", "End", "Num.of.Markers", "Aberration")

# Substitute Chromosomes "X" and "Y" with "23" and "24"
xidx <- which(cnvMatrix$Chromosome=="X")
yidx <- which(cnvMatrix$Chromosome=="Y")
cnvMatrix[xidx,"Chromosome"] <- 23
cnvMatrix[yidx,"Chromosome"] <- 24
cnvMatrix$Chromosome <- sapply(cnvMatrix$Chromosome,as.integer)

# Recurrent CNV identification with GAIA
gdac.root <- "ftp://ftp.broadinstitute.org/pub/GISTIC2.0/hg19_support/"
file <- paste0(gdac.root, "genome.info.6.0_hg19.na31_minus_frequent_nan_probes_sorted_2.1.txt")

if(!file.exists(basename(file))) download(file, basename(file))
markersMatrix <- readr::read_tsv(basename(file), col_names = TRUE, col_types = "ccn", progress = TRUE)
colnames(markersMatrix) <- c("Probe.Name", "Chromosome", "Start")

unique(markersMatrix$Chromosome)
xidx <- which(markersMatrix$Chromosome=="X")
yidx <- which(markersMatrix$Chromosome=="Y")
markersMatrix[xidx,"Chromosome"] <- 23
markersMatrix[yidx,"Chromosome"] <- 24
markersMatrix$Chromosome <- sapply(markersMatrix$Chromosome,as.integer)
markerID <- apply(markersMatrix,1,function(x) paste0(x[2],":",x[3]))
print(table(duplicated(markerID)))

print(table(duplicated(markersMatrix$Probe.Name)))

markersMatrix <- markersMatrix[-which(duplicated(markerID)),]

markerID <- apply(markersMatrix,1,function(x) paste0(x[2],":",x[3]))

file <- paste0(gdac.root, "CNV.hg19.bypos.111213.txt")
if(!file.exists(basename(file))) download(file, basename(file))
commonCNV <- readr::read_tsv(basename(file), progress = TRUE)
commonID <- apply(commonCNV,1,function(x) paste0(x[2],":",x[3]))
print(table(commonID %in% markerID))
print(table(markerID %in% commonID))
markersMatrix_fil <- markersMatrix[!markerID %in% commonID,]
toc()

cnvMatrix_h <- cnvMatrix[substr(cnvMatrix$Sample,1,12) %in% dataClin_h$bcr_patient_barcode,]
cnvMatrix_h_n <- cnvMatrix[substr(cnvMatrix$Sample,1,12) %in% dataClin_h_n$bcr_patient_barcode,]
cnvMatrix_h_1 <- cnvMatrix[substr(cnvMatrix$Sample,1,12) %in% dataClin_h_1$bcr_patient_barcode,]
cnvMatrix_h_2 <- cnvMatrix[substr(cnvMatrix$Sample,1,12) %in% dataClin_h_2$bcr_patient_barcode,]
cnvMatrix_h_3 <- cnvMatrix[substr(cnvMatrix$Sample,1,12) %in% dataClin_h_3$bcr_patient_barcode,]
cnvMatrix_h_4 <- cnvMatrix[substr(cnvMatrix$Sample,1,12) %in% dataClin_h_4$bcr_patient_barcode,]
cnvMatrix_h_low <- cnvMatrix[substr(cnvMatrix$Sample,1,12) %in% dataClin_h_low$bcr_patient_barcode,]
cnvMatrix_G1 <- cnvMatrix[substr(cnvMatrix$Sample,1,12) %in% dataClin_G1$bcr_patient_barcode,]
cnvMatrix_G2 <- cnvMatrix[substr(cnvMatrix$Sample,1,12) %in% dataClin_G2$bcr_patient_barcode,]

cnvMatrix_ESCA <- cnvMatrix

gaiaCNVplot <- function (calls, cancer=NULL, threshold=0.01)
{
  Calls <- calls[order(calls[,"Region Start [bp]"]),]
  Calls <- Calls[order(Calls[,"Chromosome"]),]
  rownames(Calls) <- NULL
  Chromo <- Calls[,"Chromosome"]
  Gains <- apply(Calls,1,function(x) ifelse(x["Aberration Kind"]==1, x["score"], 0))
  Losses <- apply(Calls, 1,function(x) ifelse(x["Aberration Kind"]==0, x["score"], 0))
  plot(Gains, ylim = c(-max(Calls [,"score"]+2), max(Calls[,"score"]+2)), type = "h",
       col = "red", xlab = "Chromosome", ylab = "Score",
       #main = paste("Recurrent Copy Number Variations",cancer, sep =" - "),
       xaxt = "n")
  points(-(Losses), type = "h", col = "blue")
  abline(h = 0, cex = 4)
  abline(h = -log10(threshold), col = "orange", cex = 4, main="test")
  abline(h = log10(threshold), col = "orange", cex = 4, main="test")
  uni.chr <- unique(Chromo)
  temp <- rep(0, length(uni.chr))
  for (i in 1:length(uni.chr)) {
    temp[i] <- max(which(uni.chr[i] == Chromo))
  }
  for (i in 1:length(temp)) {
    abline(v = temp[i], col = "black", lty = "dashed" )
  }
  nChroms <- length(uni.chr)
  begin <- c()
  for (d in 1:nChroms) {
    chrom <- sum(Chromo == uni.chr[d])
    begin <- append(begin, chrom)
  }
  temp2 <- rep(0, nChroms)
  for (i in 1:nChroms) {
    if (i == 1) {
      temp2[1] <- (begin[1] * 0.5)
    }
    else if (i > 1) {
      temp2[i] <- temp[i - 1] + (begin[i] * 0.5)
    }
  }
  uni.chr[uni.chr==23] <- "X"
  uni.chr[uni.chr==24] <- "Y"
  for (i in 1:length(temp)) {
    axis(1, at = temp2[i], labels = uni.chr[i], cex.axis = 1)
  }
  legend(x=1,y=max(Calls[,"score"]+2),y.intersp=0.8, c("Amp"), pch =15, col=c("red"), text.font=3)
  legend(x=1,y=-max(Calls[,"score"]+0.5),y.intersp=0.8, c("Del"), pch =15, col=c("blue"), text.font=3)
}

quadVenn <- function(v1,v2,v3,v4,cats = c("T 1", "T 2", "T 3", "T 4")){
  #dev.off()
  venn.plot <- draw.quad.venn(
    area1 = length(v1),
    area2 = length(v2),
    area3 = length(v3),
    area4 = length(v4),
    n12 = length(intersect(v1,v2)),
    n13 = length(intersect(v1,v3)),
    n14 = length(intersect(v1,v4)),
    n23 = length(intersect(v2,v3)),
    n24 = length(intersect(v2,v4)),
    n34 = length(intersect(v3,v4)),
    n123 = length(intersect(intersect(v1,v2),v3)),
    n124 = length(intersect(intersect(v1,v2),v4)),
    n134 = length(intersect(intersect(v1,v3),v4)),
    n234 = length(intersect(intersect(v2,v3),v4)),
    n1234 = length(intersect(intersect(v1,v2),intersect(v3,v4))),
    category = cats,
    fill = c("orange", "red", "green", "blue"),
    lty = "dashed",
    cex = 2,
    cat.cex = 2,
    cat.col = c("orange", "red", "green", "blue")
  )
  
  grid.draw(venn.plot)
}

tripleVenn <- function(v1,v2,v3,cats = c("T 1", "T 2", "T 3")){
  dev.off()
  venn.plot <- draw.triple.venn(
    area1 = length(v1),
    area2 = length(v2),
    area3 = length(v3),
    n12 = length(intersect(v1,v2)),
    n13 = length(intersect(v1,v3)),
    n23 = length(intersect(v2,v3)),
    n123 = length(intersect(intersect(v1,v2),v3)),
    category = cats,
    fill = c("orange", "red", "green"),
    lty = "dashed",
    cex = 2,
    cat.cex = 2,
    cat.col = c("orange", "red", "green"
    )
  )
  
  grid.draw(venn.plot)
}

markers_obj <- load_markers(as.data.frame(markersMatrix_fil))
# Set q-value threshold
threshold <- 0.0001

tic("CNV ESCA...")
cancer <- "ESCA_G1"
cnvMatrix <- cnvMatrix_G1
nbsamples <- length(unique(cnvMatrix$Sample))
cnv_obj <- load_cnv(cnvMatrix, markers_obj, nbsamples)
results <- runGAIA(cnv_obj, markers_obj, output_file_name=paste0("GAIA_",cancer,"_flt.txt"),
                   aberrations = -1, chromosomes = -1, num_iterations = 100, threshold = 0.25)
RecCNV <- t(apply(results,1,as.numeric))
colnames(RecCNV) <- colnames(results)
RecCNV <- cbind(RecCNV, score=0)
minval <- format(min(RecCNV[RecCNV[,"q-value"]!=0,"q-value"]),scientific=FALSE)
minval <- substring(minval,1, nchar(minval)-1)
RecCNV[RecCNV[,"q-value"]==0,"q-value"] <- as.numeric(minval)
RecCNV[,"score"] <- sapply(RecCNV[,"q-value"],function(x) -log10(as.numeric(x)))
RecCNV[RecCNV[,"q-value"]==as.numeric(minval),]
RecCNV_G1 <- RecCNV

png(filename = paste0("gaia_CNV_",cancer,"_plot.png"),width = 800)
gaiaCNVplot(RecCNV_h,cancer,threshold)
dev.off()
toc()

tic("CNV ESCA...")
cancer <- "ESCA_G2"
cnvMatrix <- cnvMatrix_G2
nbsamples <- length(unique(cnvMatrix$Sample))
cnv_obj <- load_cnv(cnvMatrix, markers_obj, nbsamples)
results <- runGAIA(cnv_obj, markers_obj, output_file_name=paste0("GAIA_",cancer,"_flt.txt"),
                   aberrations = -1, chromosomes = -1, num_iterations = 100, threshold = 0.25)
RecCNV <- t(apply(results,1,as.numeric))
colnames(RecCNV) <- colnames(results)
RecCNV <- cbind(RecCNV, score=0)
minval <- format(min(RecCNV[RecCNV[,"q-value"]!=0,"q-value"]),scientific=FALSE)
minval <- substring(minval,1, nchar(minval)-1)
RecCNV[RecCNV[,"q-value"]==0,"q-value"] <- as.numeric(minval)
RecCNV[,"score"] <- sapply(RecCNV[,"q-value"],function(x) -log10(as.numeric(x)))
RecCNV[RecCNV[,"q-value"]==as.numeric(minval),]
RecCNV_G2 <- RecCNV

png(filename = paste0("gaia_CNV_",cancer,"_plot.png"),width = 800)
gaiaCNVplot(RecCNV_h,cancer,threshold)
dev.off()
toc()

tic("CNV ESCA...")
cancer <- "ESCA"
cnvMatrix <- cnvMatrix_h
nbsamples <- length(unique(cnvMatrix$Sample))
cnv_obj <- load_cnv(cnvMatrix, markers_obj, nbsamples)
results <- runGAIA(cnv_obj, markers_obj, output_file_name=paste0("GAIA_",cancer,"_flt.txt"),
                   aberrations = -1, chromosomes = -1, num_iterations = 100, threshold = 0.25)
RecCNV <- t(apply(results,1,as.numeric))
colnames(RecCNV) <- colnames(results)
RecCNV <- cbind(RecCNV, score=0)
minval <- format(min(RecCNV[RecCNV[,"q-value"]!=0,"q-value"]),scientific=FALSE)
minval <- substring(minval,1, nchar(minval)-1)
RecCNV[RecCNV[,"q-value"]==0,"q-value"] <- as.numeric(minval)
RecCNV[,"score"] <- sapply(RecCNV[,"q-value"],function(x) -log10(as.numeric(x)))
RecCNV[RecCNV[,"q-value"]==as.numeric(minval),]
RecCNV_h <- RecCNV

png(filename = paste0("gaia_CNV_",cancer,"_plot.png"),width = 800)
gaiaCNVplot(RecCNV_h,cancer,threshold)
dev.off()
toc()

cancer <- "ESCA_N"
cnvMatrix <- cnvMatrix_h_n
nbsamples <- length(unique(cnvMatrix$Sample))
cnv_obj <- load_cnv(cnvMatrix, markers_obj, nbsamples)
results <- runGAIA(cnv_obj, markers_obj, output_file_name=paste0("GAIA_",cancer,"_flt.txt"),
                   aberrations = -1, chromosomes = -1, num_iterations = 10, threshold = 0.25)
RecCNV <- t(apply(results,1,as.numeric))
colnames(RecCNV) <- colnames(results)
RecCNV <- cbind(RecCNV, score=0)
minval <- format(min(RecCNV[RecCNV[,"q-value"]!=0,"q-value"]),scientific=FALSE)
minval <- substring(minval,1, nchar(minval)-1)
RecCNV[RecCNV[,"q-value"]==0,"q-value"] <- as.numeric(minval)
RecCNV[,"score"] <- sapply(RecCNV[,"q-value"],function(x) -log10(as.numeric(x)))
RecCNV[RecCNV[,"q-value"]==as.numeric(minval),]
RecCNV_h_n <- RecCNV
png(filename = paste0("gaia_CNV_",cancer,"_plot.png"),width = 800)
gaiaCNVplot(RecCNV_h_n,cancer,threshold)
dev.off()

cancer <- "ESCA_T1"
cnvMatrix <- cnvMatrix_h_1
nbsamples <- length(unique(cnvMatrix$Sample))
cnv_obj <- load_cnv(cnvMatrix, markers_obj, nbsamples)
results <- runGAIA(cnv_obj, markers_obj, output_file_name=paste0("GAIA_",cancer,"_flt.txt"),
                   aberrations = -1, chromosomes = -1, num_iterations = 10, threshold = 0.25)
RecCNV <- t(apply(results,1,as.numeric))
colnames(RecCNV) <- colnames(results)
RecCNV <- cbind(RecCNV, score=0)
minval <- format(min(RecCNV[RecCNV[,"q-value"]!=0,"q-value"]),scientific=FALSE)
minval <- substring(minval,1, nchar(minval)-1)
RecCNV[RecCNV[,"q-value"]==0,"q-value"] <- as.numeric(minval)
RecCNV[,"score"] <- sapply(RecCNV[,"q-value"],function(x) -log10(as.numeric(x)))
RecCNV[RecCNV[,"q-value"]==as.numeric(minval),]
RecCNV_h_1 <- RecCNV
png(filename = paste0("gaia_CNV_",cancer,"_plot.png"),width = 800)
gaiaCNVplot(RecCNV_h_1,cancer,threshold)
dev.off()

cancer <- "ESCA_T2"
cnvMatrix <- cnvMatrix_h_2
nbsamples <- length(unique(cnvMatrix$Sample))
cnv_obj <- load_cnv(cnvMatrix, markers_obj, nbsamples)
results <- runGAIA(cnv_obj, markers_obj, output_file_name=paste0("GAIA_",cancer,"_flt.txt"),
                   aberrations = -1, chromosomes = -1, num_iterations = 10, threshold = 0.25)
RecCNV <- t(apply(results,1,as.numeric))
colnames(RecCNV) <- colnames(results)
RecCNV <- cbind(RecCNV, score=0)
minval <- format(min(RecCNV[RecCNV[,"q-value"]!=0,"q-value"]),scientific=FALSE)
minval <- substring(minval,1, nchar(minval)-1)
RecCNV[RecCNV[,"q-value"]==0,"q-value"] <- as.numeric(minval)
RecCNV[,"score"] <- sapply(RecCNV[,"q-value"],function(x) -log10(as.numeric(x)))
RecCNV[RecCNV[,"q-value"]==as.numeric(minval),]
RecCNV_h_2 <- RecCNV
png(filename = paste0("gaia_CNV_",cancer,"_plot.png"),width = 800)
gaiaCNVplot(RecCNV_h_2,cancer,threshold)
dev.off()

cancer <- "ESCA_T3"
cnvMatrix <- cnvMatrix_h_3
nbsamples <- length(unique(cnvMatrix$Sample))
cnv_obj <- load_cnv(cnvMatrix, markers_obj, nbsamples)
results <- runGAIA(cnv_obj, markers_obj, output_file_name=paste0("GAIA_",cancer,"_flt.txt"),
                   aberrations = -1, chromosomes = -1, num_iterations = 10, threshold = 0.25)
RecCNV <- t(apply(results,1,as.numeric))
colnames(RecCNV) <- colnames(results)
RecCNV <- cbind(RecCNV, score=0)
minval <- format(min(RecCNV[RecCNV[,"q-value"]!=0,"q-value"]),scientific=FALSE)
minval <- substring(minval,1, nchar(minval)-1)
RecCNV[RecCNV[,"q-value"]==0,"q-value"] <- as.numeric(minval)
RecCNV[,"score"] <- sapply(RecCNV[,"q-value"],function(x) -log10(as.numeric(x)))
RecCNV[RecCNV[,"q-value"]==as.numeric(minval),]
RecCNV_h_3 <- RecCNV
png(filename = paste0("gaia_CNV_",cancer,"_plot.png"),width = 800)
gaiaCNVplot(RecCNV_h_3,cancer,threshold)
dev.off()

cancer <- "ESCA_T4"
cnvMatrix <- cnvMatrix_h_4
nbsamples <- length(unique(cnvMatrix$Sample))
cnv_obj <- load_cnv(cnvMatrix, markers_obj, nbsamples)
results <- runGAIA(cnv_obj, markers_obj, output_file_name=paste0("GAIA_",cancer,"_flt.txt"),
                   aberrations = -1, chromosomes = -1, num_iterations = 10, threshold = 0.25)
RecCNV <- t(apply(results,1,as.numeric))
colnames(RecCNV) <- colnames(results)
RecCNV <- cbind(RecCNV, score=0)
minval <- format(min(RecCNV[RecCNV[,"q-value"]!=0,"q-value"]),scientific=FALSE)
minval <- substring(minval,1, nchar(minval)-1)
RecCNV[RecCNV[,"q-value"]==0,"q-value"] <- as.numeric(minval)
RecCNV[,"score"] <- sapply(RecCNV[,"q-value"],function(x) -log10(as.numeric(x)))
RecCNV[RecCNV[,"q-value"]==as.numeric(minval),]
RecCNV_h_4 <- RecCNV
png(filename = paste0("gaia_CNV_",cancer,"_plot.png"),width = 800)
gaiaCNVplot(RecCNV_h_4,cancer,threshold)
dev.off()


cancer <- "ESCA_T_Low"
cnvMatrix <- cnvMatrix_h_low
nbsamples <- length(unique(cnvMatrix$Sample))
cnv_obj <- load_cnv(cnvMatrix, markers_obj, nbsamples)
results <- runGAIA(cnv_obj, markers_obj, output_file_name=paste0("GAIA_",cancer,"_flt.txt"),
                   aberrations = -1, chromosomes = -1, num_iterations = 10, threshold = 0.25)
RecCNV <- t(apply(results,1,as.numeric))
colnames(RecCNV) <- colnames(results)
RecCNV <- cbind(RecCNV, score=0)
minval <- format(min(RecCNV[RecCNV[,"q-value"]!=0,"q-value"]),scientific=FALSE)
minval <- substring(minval,1, nchar(minval)-1)
RecCNV[RecCNV[,"q-value"]==0,"q-value"] <- as.numeric(minval)
RecCNV[,"score"] <- sapply(RecCNV[,"q-value"],function(x) -log10(as.numeric(x)))
RecCNV[RecCNV[,"q-value"]==as.numeric(minval),]
RecCNV_h_low <- RecCNV
png(filename = paste0("gaia_CNV_",cancer,"_plot.png"),width = 800)
gaiaCNVplot(RecCNV,cancer,threshold)
dev.off()

toc()


##################################


############## DEG
##############

################### DEG Pathway
tic("DEG...")
tic("Prepare...")
dp <- TCGAanalyze_Preprocessing(object = rse_h_n,cor.cut = 0.6,datatype = "raw_count",filename = paste0(p_h,"_Norm.png"))
#dn <- TCGAanalyze_Normalization(tabDF = dp,geneInfo = geneInfo,method = "gcContent")
dn <- TCGAanalyze_Normalization(tabDF = dp,geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df_h_n <- subset(df, select = substr(colnames(df),1,12) %in% dataClin_h$bcr_patient_barcode)
dp_h_n <- dp

dp <- TCGAanalyze_Preprocessing(object = rse_h_c_n,cor.cut = 0.6,datatype = "raw_count",filename = paste0(p_h,"_Cancer_Not.png"))
#dn <- TCGAanalyze_Normalization(tabDF = dp,geneInfo = geneInfo,method = "gcContent")
dn <- TCGAanalyze_Normalization(tabDF = dp,geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df_h_c_n <- subset(df, select = substr(colnames(df),1,12) %in% dataClin_h_n$bcr_patient_barcode)
dp_h_c_n <- dp

dp <- TCGAanalyze_Preprocessing(object = rse_h_c_1,cor.cut = 0.6,datatype = "raw_count",filename = paste0(p_h,"_Cancer_T1.png"))
#dn <- TCGAanalyze_Normalization(tabDF = dp,geneInfo = geneInfo,method = "gcContent")
dn <- TCGAanalyze_Normalization(tabDF = dp,geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df_h_c_1 <- subset(df, select = substr(colnames(df),1,12) %in% dataClin_h_1$bcr_patient_barcode)
dp_h_c_1 <- dp

dp <- TCGAanalyze_Preprocessing(object = rse_h_c_2,cor.cut = 0.6,datatype = "raw_count",filename = paste0(p_h,"_Cancer_T2.png"))
#dn <- TCGAanalyze_Normalization(tabDF = dp,geneInfo = geneInfo,method = "gcContent")
dn <- TCGAanalyze_Normalization(tabDF = dp,geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df_h_c_2 <- subset(df, select = substr(colnames(df),1,12) %in% dataClin_h_2$bcr_patient_barcode)
dp_h_c_2 <- dp

dp <- TCGAanalyze_Preprocessing(object = rse_h_c_3,cor.cut = 0.6,datatype = "raw_count",filename = paste0(p_h,"_Cancer_T3.png"))
#dn <- TCGAanalyze_Normalization(tabDF = dp,geneInfo = geneInfo,method = "gcContent")
dn <- TCGAanalyze_Normalization(tabDF = dp,geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df_h_c_3 <- subset(df, select = substr(colnames(df),1,12) %in% dataClin_h_3$bcr_patient_barcode)
dp_h_c_3 <- dp

dp <- TCGAanalyze_Preprocessing(object = rse_h_c_4,cor.cut = 0.6,datatype = "raw_count",filename = paste0(p_h,"_Cancer_T4.png"))
#dn <- TCGAanalyze_Normalization(tabDF = dp,geneInfo = geneInfo,method = "gcContent")
dn <- TCGAanalyze_Normalization(tabDF = dp,geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df_h_c_4 <- subset(df, select = substr(colnames(df),1,12) %in% dataClin_h_4$bcr_patient_barcode)
dp_h_c_4 <- dp

dp <- TCGAanalyze_Preprocessing(object = rse_h_c_low,cor.cut = 0.6,datatype = "raw_count",filename = paste0(p_h,"_Cancer_Low.png"))
#dn <- TCGAanalyze_Normalization(tabDF = dp,geneInfo = geneInfo,method = "gcContent")
dn <- TCGAanalyze_Normalization(tabDF = dp,geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df_h_c_low <- subset(df, select = substr(colnames(df),1,12) %in% dataClin_h_low$bcr_patient_barcode)
dp_h_c_low <- dp

dp <- TCGAanalyze_Preprocessing(object = rse_G1,cor.cut = 0.6,datatype = "raw_count",filename = paste0(p_h,"_G1.png"))
#dn <- TCGAanalyze_Normalization(tabDF = dp,geneInfo = geneInfo,method = "gcContent")
dn <- TCGAanalyze_Normalization(tabDF = dp,geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df_G1 <- subset(df, select = substr(colnames(df),1,12) %in% dataClin_G1$bcr_patient_barcode)
dp_G1 <- dp

dp <- TCGAanalyze_Preprocessing(object = rse_G2,cor.cut = 0.6,datatype = "raw_count",filename = paste0(p_h,"_G2.png"))
#dn <- TCGAanalyze_Normalization(tabDF = dp,geneInfo = geneInfo,method = "gcContent")
dn <- TCGAanalyze_Normalization(tabDF = dp,geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df_G2 <- subset(df, select = substr(colnames(df),1,12) %in% dataClin_G2$bcr_patient_barcode)
dp_G2 <- dp

dp <- TCGAanalyze_Preprocessing(object = rse_h_c,cor.cut = 0.6,datatype = "raw_count",filename = paste0(p_h,"_Cancer.png"))
#dn <- TCGAanalyze_Normalization(tabDF = dp,geneInfo = geneInfo,method = "gcContent")
dn <- TCGAanalyze_Normalization(tabDF = dp,geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df_h_c <- subset(df, select = substr(colnames(df),1,12) %in% dataClin_h$bcr_patient_barcode)
dp_h_c <- dp
toc()
################ DEG 

GenelistComplete <- rownames(assay(rse_h_c,1))
p1 <- p_h
p2 <- p_h

# N vs C

tic("DEG N C...")
dn <- TCGAanalyze_Normalization(tabDF = cbind(dp_h_n, dp_h_c),geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df1 <- subset(df, select = colnames(df) %in% colnames(df_h_n))
S1 <- TCGAquery_SampleTypes(colnames(df1),"NT")
C1 <- "Cancer_N"
df2 <- subset(df, select = colnames(df) %in% colnames(df_h_c))
S2 <- TCGAquery_SampleTypes(colnames(df2),"TP")
C2 <- "Cancer_C"
dataDEGs <- TCGAanalyze_DEA(mat1 = df1[,S1],mat2 = df2[,S2],Cond1type = paste0(p1,"_",C1),Cond2type = paste0(p2,"_",C2),fdr.cut = 0.01 ,logFC.cut = 1)
ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes ",p1,"_",C1," Vs ",p2,"_",C2), RegulonList = rownames(dataDEGs))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = rownames(dataDEGs),nBar = 10,filename = paste0("DEA_genes_",p1,"_",C1,"_vs_",p2,"_",C2,"_20.pdf"))
dataDEGs_h_n_c <- dataDEGs
ansEA_h_n_c <- ansEA

dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,paste0(p1,"_",C1),paste0(p2,"_",C2),df1[,S1],df2[,S2])
dataDEGsFiltLevel$GeneID <- 0
eg = as.data.frame ( bitr(dataDEGsFiltLevel$mRNA,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db"))
eg <- eg[!duplicated(eg$SYMBOL),]
dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL,]
dataDEGsFiltLevel <- dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA,decreasing=FALSE),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]
dataDEGsFiltLevel$GeneID <- eg$ENTREZID
dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID
hsa05200 <-  pathview(gene.data  = genelistDEGs,pathway.id = "hsa05200",species = "hsa",limit = list(gene=as.integer(max(abs(genelistDEGs)))))
file.rename("hsa05200.pathview.png",paste0("hsa05200_",p_h,"_",C1,"_",C2,".png"))
genelistDEGs_n_c <- genelistDEGs

toc()
#N vs C1
dn <- TCGAanalyze_Normalization(tabDF = cbind(dp_h_n, dp_h_c_1),geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df1 <- subset(df, select = colnames(df) %in% colnames(df_h_n))
S1 <- TCGAquery_SampleTypes(colnames(df1),"NT")
C1 <- "Cancer_N"
df2 <- subset(df, select = colnames(df) %in% colnames(df_h_c_1))
S2 <- TCGAquery_SampleTypes(colnames(df2),"TP")
C2 <- "Cancer_C_1"
dataDEGs <- TCGAanalyze_DEA(mat1 = df1[,S1],mat2 = df2[,S2],Cond1type = paste0(p1,"_",C1),Cond2type = paste0(p2,"_",C2),fdr.cut = 0.01 ,logFC.cut = 1)
ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes ",p1,"_",C1," Vs ",p2,"_",C2), RegulonList = rownames(dataDEGs))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = rownames(dataDEGs),nBar = 20,filename = paste0("DEA_genes_",p1,"_",C1,"_vs_",p2,"_",C2,"_20.pdf"))
dataDEGs_h_n_c1 <- dataDEGs
ansEA_h_n_c1 <- ansEA

dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,paste0(p1,"_",C1),paste0(p2,"_",C2),df1[,S1],df2[,S2])
dataDEGsFiltLevel$GeneID <- 0
eg = as.data.frame ( bitr(dataDEGsFiltLevel$mRNA,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db"))
eg <- eg[!duplicated(eg$SYMBOL),]
dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL,]
dataDEGsFiltLevel <- dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA,decreasing=FALSE),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]
dataDEGsFiltLevel$GeneID <- eg$ENTREZID
dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID
hsa05200 <-  pathview(gene.data  = genelistDEGs,pathway.id = "hsa05200",species = "hsa",limit = list(gene=as.integer(max(abs(genelistDEGs)))))
file.rename("hsa05200.pathview.png",paste0("hsa05200_",p_h,"_",C1,"_",C2,".png"))
genelistDEGs_n_c1 <- genelistDEGs

#N vs C2
dn <- TCGAanalyze_Normalization(tabDF = cbind(dp_h_n, dp_h_c_2),geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df1 <- subset(df, select = colnames(df) %in% colnames(df_h_n))
S1 <- TCGAquery_SampleTypes(colnames(df1),"NT")
C1 <- "Cancer_N"
df2 <- subset(df, select = colnames(df) %in% colnames(df_h_c_2))
S2 <- TCGAquery_SampleTypes(colnames(df2),"TP")
C2 <- "Cancer_C_2"
dataDEGs <- TCGAanalyze_DEA(mat1 = df1[,S1],mat2 = df2[,S2],Cond1type = paste0(p1,"_",C1),Cond2type = paste0(p2,"_",C2),fdr.cut = 0.01 ,logFC.cut = 1)
ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes ",p1,"_",C1," Vs ",p2,"_",C2), RegulonList = rownames(dataDEGs))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = rownames(dataDEGs),nBar = 20,filename = paste0("DEA_genes_",p1,"_",C1,"_vs_",p2,"_",C2,"_20.pdf"))
dataDEGs_h_n_c2 <- dataDEGs
ansEA_h_n_c2 <- ansEA

dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,paste0(p1,"_",C1),paste0(p2,"_",C2),df1[,S1],df2[,S2])
dataDEGsFiltLevel$GeneID <- 0
eg = as.data.frame ( bitr(dataDEGsFiltLevel$mRNA,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db"))
eg <- eg[!duplicated(eg$SYMBOL),]
dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL,]
dataDEGsFiltLevel <- dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA,decreasing=FALSE),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]
dataDEGsFiltLevel$GeneID <- eg$ENTREZID
dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID
hsa05200 <-  pathview(gene.data  = genelistDEGs,pathway.id = "hsa05200",species = "hsa",limit = list(gene=as.integer(max(abs(genelistDEGs)))))
file.rename("hsa05200.pathview.png",paste0("hsa05200_",p_h,"_",C1,"_",C2,".png"))
genelistDEGs_n_c2 <- genelistDEGs

#N vs C3
dn <- TCGAanalyze_Normalization(tabDF = cbind(dp_h_n, dp_h_c_3),geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df1 <- subset(df, select = colnames(df) %in% colnames(df_h_n))
S1 <- TCGAquery_SampleTypes(colnames(df1),"NT")
C1 <- "Cancer_N"
df2 <- subset(df, select = colnames(df) %in% colnames(df_h_c_3))
S2 <- TCGAquery_SampleTypes(colnames(df2),"TP")
C2 <- "Cancer_C_3"
dataDEGs <- TCGAanalyze_DEA(mat1 = df1[,S1],mat2 = df2[,S2],Cond1type = paste0(p1,"_",C1),Cond2type = paste0(p2,"_",C2),fdr.cut = 0.01 ,logFC.cut = 1)
ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes ",p1,"_",C1," Vs ",p2,"_",C2), RegulonList = rownames(dataDEGs))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = rownames(dataDEGs),nBar = 20,filename = paste0("DEA_genes_",p1,"_",C1,"_vs_",p2,"_",C2,"_20.pdf"))
dataDEGs_h_n_c3 <- dataDEGs
ansEA_h_n_c3 <- ansEA
genelistDEGs_n_c3 <- genelistDEGs

dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,paste0(p1,"_",C1),paste0(p2,"_",C2),df1[,S1],df2[,S2])
dataDEGsFiltLevel$GeneID <- 0
eg = as.data.frame ( bitr(dataDEGsFiltLevel$mRNA,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db"))
eg <- eg[!duplicated(eg$SYMBOL),]
dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL,]
dataDEGsFiltLevel <- dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA,decreasing=FALSE),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]
dataDEGsFiltLevel$GeneID <- eg$ENTREZID
dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID
hsa05200 <-  pathview(gene.data  = genelistDEGs,pathway.id = "hsa05200",species = "hsa",limit = list(gene=as.integer(max(abs(genelistDEGs)))))
file.rename("hsa05200.pathview.png",paste0("hsa05200_",p_h,"_",C1,"_",C2,".png"))
genelistDEGs_n_c3 <- genelistDEGs

#N vs C4
dn <- TCGAanalyze_Normalization(tabDF = cbind(dp_h_n, dp_h_c_4),geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df1 <- subset(df, select = colnames(df) %in% colnames(df_h_n))
S1 <- TCGAquery_SampleTypes(colnames(df1),"NT")
C1 <- "Cancer_N"
df2 <- subset(df, select = colnames(df) %in% colnames(df_h_c_4))
S2 <- TCGAquery_SampleTypes(colnames(df2),"TP")
C2 <- "Cancer_C_4"
dataDEGs <- TCGAanalyze_DEA(mat1 = df1[,S1],mat2 = df2[,S2],Cond1type = paste0(p1,"_",C1),Cond2type = paste0(p2,"_",C2),fdr.cut = 0.01 ,logFC.cut = 1)
ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes ",p1,"_",C1," Vs ",p2,"_",C2), RegulonList = rownames(dataDEGs))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = rownames(dataDEGs),nBar = 20,filename = paste0("DEA_genes_",p1,"_",C1,"_vs_",p2,"_",C2,"_20.pdf"))
dataDEGs_h_n_c4 <- dataDEGs
ansEA_h_n_c4 <- ansEA

dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,paste0(p1,"_",C1),paste0(p2,"_",C2),df1[,S1],df2[,S2])
dataDEGsFiltLevel$GeneID <- 0
eg = as.data.frame ( bitr(dataDEGsFiltLevel$mRNA,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db"))
eg <- eg[!duplicated(eg$SYMBOL),]
dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL,]
dataDEGsFiltLevel <- dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA,decreasing=FALSE),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]
dataDEGsFiltLevel$GeneID <- eg$ENTREZID
dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID
hsa05200 <-  pathview(gene.data  = genelistDEGs,pathway.id = "hsa05200",species = "hsa",limit = list(gene=as.integer(max(abs(genelistDEGs)))))
file.rename("hsa05200.pathview.png",paste0("hsa05200_",p_h,"_",C1,"_",C2,".png"))
genelistDEGs_n_c4 <- genelistDEGs

#N vs low
dn <- TCGAanalyze_Normalization(tabDF = cbind(dp_h_n, dp_h_c_low),geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df1 <- subset(df, select = colnames(df) %in% colnames(df_h_n))
S1 <- TCGAquery_SampleTypes(colnames(df1),"NT")
C1 <- "Cancer_N"
df2 <- subset(df, select = colnames(df) %in% colnames(df_h_c_low))
S2 <- TCGAquery_SampleTypes(colnames(df2),"TP")
C2 <- "Cancer_C_Low"
dataDEGs <- TCGAanalyze_DEA(mat1 = df1[,S1],mat2 = df2[,S2],Cond1type = paste0(p1,"_",C1),Cond2type = paste0(p2,"_",C2),fdr.cut = 0.01 ,logFC.cut = 1)
ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes ",p1,"_",C1," Vs ",p2,"_",C2), RegulonList = rownames(dataDEGs))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = rownames(dataDEGs),nBar = 20,filename = paste0("DEA_genes_",p1,"_",C1,"_vs_",p2,"_",C2,"_20.pdf"))
dataDEGs_h_n_c_low <- dataDEGs
ansEA_h_n_c_low <- ansEA

dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,paste0(p1,"_",C1),paste0(p2,"_",C2),df1[,S1],df2[,S2])
dataDEGsFiltLevel$GeneID <- 0
eg = as.data.frame ( bitr(dataDEGsFiltLevel$mRNA,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db"))
eg <- eg[!duplicated(eg$SYMBOL),]
dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL,]
dataDEGsFiltLevel <- dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA,decreasing=FALSE),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]
dataDEGsFiltLevel$GeneID <- eg$ENTREZID
dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID
hsa05200 <-  pathview(gene.data  = genelistDEGs,pathway.id = "hsa05200",species = "hsa",limit = list(gene=as.integer(max(abs(genelistDEGs)))))
file.rename("hsa05200.pathview.png",paste0("hsa05200_",p_h,"_",C1,"_",C2,".png"))
genelistDEGs_n_c_low <- genelistDEGs

#N vs G1
dn <- TCGAanalyze_Normalization(tabDF = cbind(dp_h_n, dp_G1),geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df1 <- subset(df, select = colnames(df) %in% colnames(df_h_n))
S1 <- TCGAquery_SampleTypes(colnames(df1),"NT")
C1 <- "Cancer_N"
df2 <- subset(df, select = colnames(df) %in% colnames(df_G1))
S2 <- TCGAquery_SampleTypes(colnames(df2),"TP")
C2 <- "Cancer_G1"
dataDEGs <- TCGAanalyze_DEA(mat1 = df1[,S1],mat2 = df2[,S2],Cond1type = paste0(p1,"_",C1),Cond2type = paste0(p2,"_",C2),fdr.cut = 0.01 ,logFC.cut = 1)
ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes ",p1,"_",C1," Vs ",p2,"_",C2), RegulonList = rownames(dataDEGs))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = rownames(dataDEGs),nBar = 20,filename = paste0("DEA_genes_",p1,"_",C1,"_vs_",p2,"_",C2,"_20.pdf"))
dataDEGs_n_g1 <- dataDEGs
ansEA_n_g1 <- ansEA

dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,paste0(p1,"_",C1),paste0(p2,"_",C2),df1[,S1],df2[,S2])
dataDEGsFiltLevel$GeneID <- 0
eg = as.data.frame ( bitr(dataDEGsFiltLevel$mRNA,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db"))
eg <- eg[!duplicated(eg$SYMBOL),]
dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL,]
dataDEGsFiltLevel <- dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA,decreasing=FALSE),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]
dataDEGsFiltLevel$GeneID <- eg$ENTREZID
dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID
hsa05200 <-  pathview(gene.data  = genelistDEGs,pathway.id = "hsa05200",species = "hsa",limit = list(gene=as.integer(max(abs(genelistDEGs)))))
file.rename("hsa05200.pathview.png",paste0("hsa05200_",p_h,"_",C1,"_",C2,".png"))
genelistDEGs_n_G1 <- genelistDEGs


#N vs G2
dn <- TCGAanalyze_Normalization(tabDF = cbind(dp_h_n, dp_G2),geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df1 <- subset(df, select = colnames(df) %in% colnames(df_h_n))
S1 <- TCGAquery_SampleTypes(colnames(df1),"NT")
C1 <- "Cancer_N"
df2 <- subset(df, select = colnames(df) %in% colnames(df_G2))
S2 <- TCGAquery_SampleTypes(colnames(df2),"TP")
C2 <- "Cancer_G2"
dataDEGs <- TCGAanalyze_DEA(mat1 = df1[,S1],mat2 = df2[,S2],Cond1type = paste0(p1,"_",C1),Cond2type = paste0(p2,"_",C2),fdr.cut = 0.01 ,logFC.cut = 1)
ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes ",p1,"_",C1," Vs ",p2,"_",C2), RegulonList = rownames(dataDEGs))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = rownames(dataDEGs),nBar = 20,filename = paste0("DEA_genes_",p1,"_",C1,"_vs_",p2,"_",C2,"_20.pdf"))
dataDEGs_n_g2 <- dataDEGs
ansEA_n_g2 <- ansEA

dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,paste0(p1,"_",C1),paste0(p2,"_",C2),df1[,S1],df2[,S2])
dataDEGsFiltLevel$GeneID <- 0
eg = as.data.frame ( bitr(dataDEGsFiltLevel$mRNA,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db"))
eg <- eg[!duplicated(eg$SYMBOL),]
dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL,]
dataDEGsFiltLevel <- dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA,decreasing=FALSE),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]
dataDEGsFiltLevel$GeneID <- eg$ENTREZID
dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID
hsa05200 <-  pathview(gene.data  = genelistDEGs,pathway.id = "hsa05200",species = "hsa",limit = list(gene=as.integer(max(abs(genelistDEGs)))))
file.rename("hsa05200.pathview.png",paste0("hsa05200_",p_h,"_",C1,"_",C2,".png"))
genelistDEGs_n_G2 <- genelistDEGs

#G1 vs G2
dn <- TCGAanalyze_Normalization(tabDF = cbind(dp_G1, dp_G2),geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df1 <- subset(df, select = colnames(df) %in% colnames(df_G1))
S1 <- TCGAquery_SampleTypes(colnames(df1),"TP")
C1 <- "Cancer_G_1"
df2 <- subset(df, select = colnames(df) %in% colnames(df_G2))
S2 <- TCGAquery_SampleTypes(colnames(df2),"TP")
C2 <- "Cancer_G_2"
dataDEGs <- TCGAanalyze_DEA(mat1 = df1[,S1],mat2 = df2[,S2],Cond1type = paste0(p1,"_",C1),Cond2type = paste0(p2,"_",C2),fdr.cut = 0.01 ,logFC.cut = 1)
ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes ",p1,"_",C1," Vs ",p2,"_",C2), RegulonList = rownames(dataDEGs))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = rownames(dataDEGs),nBar = 20,filename = paste0("DEA_genes_",p1,"_",C1,"_vs_",p2,"_",C2,"_20.pdf"))
dataDEGs_g1_g2 <- dataDEGs
ansEA_g1_g2 <- ansEA

dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,paste0(p1,"_",C1),paste0(p2,"_",C2),df1[,S1],df2[,S2])
dataDEGsFiltLevel$GeneID <- 0
eg = as.data.frame ( bitr(dataDEGsFiltLevel$mRNA,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db"))
eg <- eg[!duplicated(eg$SYMBOL),]
dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL,]
dataDEGsFiltLevel <- dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA,decreasing=FALSE),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]
dataDEGsFiltLevel$GeneID <- eg$ENTREZID
dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID
hsa05200 <-  pathview(gene.data  = genelistDEGs,pathway.id = "hsa05200",species = "hsa",limit = list(gene=as.integer(max(abs(genelistDEGs)))))
file.rename("hsa05200.pathview.png",paste0("hsa05200_",p_h,"_",C1,"_",C2,".png"))
genelistDEGs_g1_g2 <- genelistDEGs

#C1 vs C2
dn <- TCGAanalyze_Normalization(tabDF = cbind(dp_h_c_1, dp_h_c_2),geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df1 <- subset(df, select = colnames(df) %in% colnames(df_h_c_1))
S1 <- TCGAquery_SampleTypes(colnames(df1),"TP")
C1 <- "Cancer_C_1"
df2 <- subset(df, select = colnames(df) %in% colnames(df_h_c_2))
S2 <- TCGAquery_SampleTypes(colnames(df2),"TP")
C2 <- "Cancer_C_2"
dataDEGs <- TCGAanalyze_DEA(mat1 = df1[,S1],mat2 = df2[,S2],Cond1type = paste0(p1,"_",C1),Cond2type = paste0(p2,"_",C2),fdr.cut = 0.01 ,logFC.cut = 1)
ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes ",p1,"_",C1," Vs ",p2,"_",C2), RegulonList = rownames(dataDEGs))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = rownames(dataDEGs),nBar = 20,filename = paste0("DEA_genes_",p1,"_",C1,"_vs_",p2,"_",C2,"_20.pdf"))

dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,paste0(p1,"_",C1),paste0(p2,"_",C2),df1[,S1],df2[,S2])
dataDEGsFiltLevel$GeneID <- 0
eg = as.data.frame ( bitr(dataDEGsFiltLevel$mRNA,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db"))
eg <- eg[!duplicated(eg$SYMBOL),]
dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL,]
dataDEGsFiltLevel <- dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA,decreasing=FALSE),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]
dataDEGsFiltLevel$GeneID <- eg$ENTREZID
dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID
hsa05200 <-  pathview(gene.data  = genelistDEGs,pathway.id = "hsa05200",species = "hsa",limit = list(gene=as.integer(max(abs(genelistDEGs)))))
file.rename("hsa05200.pathview.png",paste0("hsa05200_",p_h,"_",C1,"_",C2,".png"))
genelistDEGs_c1_c2 <- genelistDEGs

dataDEGs_h_c1_c2 <- dataDEGs
ansEA_h_c1_c2 <- ansEA

#C1 vs C3
dn <- TCGAanalyze_Normalization(tabDF = cbind(dp_h_c_1, dp_h_c_3),geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df1 <- subset(df, select = colnames(df) %in% colnames(df_h_c_1))
S1 <- TCGAquery_SampleTypes(colnames(df1),"TP")
C1 <- "Cancer_C_1"
df2 <- subset(df, select = colnames(df) %in% colnames(df_h_c_3))
S2 <- TCGAquery_SampleTypes(colnames(df2),"TP")
C2 <- "Cancer_C_3"
dataDEGs <- TCGAanalyze_DEA(mat1 = df1[,S1],mat2 = df2[,S2],Cond1type = paste0(p1,"_",C1),Cond2type = paste0(p2,"_",C2),fdr.cut = 0.01 ,logFC.cut = 1)
ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes ",p1,"_",C1," Vs ",p2,"_",C2), RegulonList = rownames(dataDEGs))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = rownames(dataDEGs),nBar = 20,filename = paste0("DEA_genes_",p1,"_",C1,"_vs_",p2,"_",C2,"_20.pdf"))
dataDEGs_h_c1_c3 <- dataDEGs
ansEA_h_c1_c3 <- ansEA

dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,paste0(p1,"_",C1),paste0(p2,"_",C2),df1[,S1],df2[,S2])
dataDEGsFiltLevel$GeneID <- 0
eg = as.data.frame ( bitr(dataDEGsFiltLevel$mRNA,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db"))
eg <- eg[!duplicated(eg$SYMBOL),]
dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL,]
dataDEGsFiltLevel <- dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA,decreasing=FALSE),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]
dataDEGsFiltLevel$GeneID <- eg$ENTREZID
dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID
hsa05200 <-  pathview(gene.data  = genelistDEGs,pathway.id = "hsa05200",species = "hsa",limit = list(gene=as.integer(max(abs(genelistDEGs)))))
file.rename("hsa05200.pathview.png",paste0("hsa05200_",p_h,"_",C1,"_",C2,".png"))
genelistDEGs_c1_c3 <- genelistDEGs

#C1 vs C4
dn <- TCGAanalyze_Normalization(tabDF = cbind(dp_h_c_1, dp_h_c_4),geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df1 <- subset(df, select = colnames(df) %in% colnames(df_h_c_1))
S1 <- TCGAquery_SampleTypes(colnames(df1),"TP")
C1 <- "Cancer_C_1"
df2 <- subset(df, select = colnames(df) %in% colnames(df_h_c_4))
S2 <- TCGAquery_SampleTypes(colnames(df2),"TP")
C2 <- "Cancer_C_4"
dataDEGs <- TCGAanalyze_DEA(mat1 = df1[,S1],mat2 = df2[,S2],Cond1type = paste0(p1,"_",C1),Cond2type = paste0(p2,"_",C2),fdr.cut = 0.01 ,logFC.cut = 1)
ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes ",p1,"_",C1," Vs ",p2,"_",C2), RegulonList = rownames(dataDEGs))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = rownames(dataDEGs),nBar = 20,filename = paste0("DEA_genes_",p1,"_",C1,"_vs_",p2,"_",C2,"_20.pdf"))
dataDEGs_h_c1_c4 <- dataDEGs
ansEA_h_c1_c4 <- ansEA

dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,paste0(p1,"_",C1),paste0(p2,"_",C2),df1[,S1],df2[,S2])
dataDEGsFiltLevel$GeneID <- 0
eg = as.data.frame ( bitr(dataDEGsFiltLevel$mRNA,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db"))
eg <- eg[!duplicated(eg$SYMBOL),]
dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL,]
dataDEGsFiltLevel <- dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA,decreasing=FALSE),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]
dataDEGsFiltLevel$GeneID <- eg$ENTREZID
dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID
hsa05200 <-  pathview(gene.data  = genelistDEGs,pathway.id = "hsa05200",species = "hsa",limit = list(gene=as.integer(max(abs(genelistDEGs)))))
file.rename("hsa05200.pathview.png",paste0("hsa05200_",p_h,"_",C1,"_",C2,".png"))
genelistDEGs_c1_c4 <- genelistDEGs

#C2 vs C3
dn <- TCGAanalyze_Normalization(tabDF = cbind(dp_h_c_2, dp_h_c_3),geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df1 <- subset(df, select = colnames(df) %in% colnames(df_h_c_2))
S1 <- TCGAquery_SampleTypes(colnames(df1),"TP")
C1 <- "Cancer_C_2"
df2 <- subset(df, select = colnames(df) %in% colnames(df_h_c_3))
S2 <- TCGAquery_SampleTypes(colnames(df2),"TP")
C2 <- "Cancer_C_3"
dataDEGs <- TCGAanalyze_DEA(mat1 = df1[,S1],mat2 = df2[,S2],Cond1type = paste0(p1,"_",C1),Cond2type = paste0(p2,"_",C2),fdr.cut = 0.01 ,logFC.cut = 1)
ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes ",p1,"_",C1," Vs ",p2,"_",C2), RegulonList = rownames(dataDEGs))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = rownames(dataDEGs),nBar = 20,filename = paste0("DEA_genes_",p1,"_",C1,"_vs_",p2,"_",C2,"_20.pdf"))
dataDEGs_h_c2_c3 <- dataDEGs
ansEA_h_c2_c3 <- ansEA

dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,paste0(p1,"_",C1),paste0(p2,"_",C2),df1[,S1],df2[,S2])
dataDEGsFiltLevel$GeneID <- 0
eg = as.data.frame ( bitr(dataDEGsFiltLevel$mRNA,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db"))
eg <- eg[!duplicated(eg$SYMBOL),]
dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL,]
dataDEGsFiltLevel <- dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA,decreasing=FALSE),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]
dataDEGsFiltLevel$GeneID <- eg$ENTREZID
dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID
hsa05200 <-  pathview(gene.data  = genelistDEGs,pathway.id = "hsa05200",species = "hsa",limit = list(gene=as.integer(max(abs(genelistDEGs)))))
file.rename("hsa05200.pathview.png",paste0("hsa05200_",p_h,"_",C1,"_",C2,".png"))
genelistDEGs_c2_c3 <- genelistDEGs

#C2 vs C4
dn <- TCGAanalyze_Normalization(tabDF = cbind(dp_h_c_2, dp_h_c_4),geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df1 <- subset(df, select = colnames(df) %in% colnames(df_h_c_2))
S1 <- TCGAquery_SampleTypes(colnames(df1),"TP")
C1 <- "Cancer_C_2"
df2 <- subset(df, select = colnames(df) %in% colnames(df_h_c_4))
S2 <- TCGAquery_SampleTypes(colnames(df2),"TP")
C2 <- "Cancer_C_4"
dataDEGs <- TCGAanalyze_DEA(mat1 = df1[,S1],mat2 = df2[,S2],Cond1type = paste0(p1,"_",C1),Cond2type = paste0(p2,"_",C2),fdr.cut = 0.01 ,logFC.cut = 1)
ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes ",p1,"_",C1," Vs ",p2,"_",C2), RegulonList = rownames(dataDEGs))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = rownames(dataDEGs),nBar = 20,filename = paste0("DEA_genes_",p1,"_",C1,"_vs_",p2,"_",C2,"_20.pdf"))
dataDEGs_h_c2_c4 <- dataDEGs
ansEA_h_c2_c4 <- ansEA

dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,paste0(p1,"_",C1),paste0(p2,"_",C2),df1[,S1],df2[,S2])
dataDEGsFiltLevel$GeneID <- 0
eg = as.data.frame ( bitr(dataDEGsFiltLevel$mRNA,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db"))
eg <- eg[!duplicated(eg$SYMBOL),]
dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL,]
dataDEGsFiltLevel <- dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA,decreasing=FALSE),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]
dataDEGsFiltLevel$GeneID <- eg$ENTREZID
dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID
hsa05200 <-  pathview(gene.data  = genelistDEGs,pathway.id = "hsa05200",species = "hsa",limit = list(gene=as.integer(max(abs(genelistDEGs)))))
file.rename("hsa05200.pathview.png",paste0("hsa05200_",p_h,"_",C1,"_",C2,".png"))
genelistDEGs_c2_c4 <- genelistDEGs

#C3 vs C4
dn <- TCGAanalyze_Normalization(tabDF = cbind(dp_h_c_3, dp_h_c_4),geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df1 <- subset(df, select = colnames(df) %in% colnames(df_h_c_3))
S1 <- TCGAquery_SampleTypes(colnames(df1),"TP")
C1 <- "Cancer_C_3"
df2 <- subset(df, select = colnames(df) %in% colnames(df_h_c_4))
S2 <- TCGAquery_SampleTypes(colnames(df2),"TP")
C2 <- "Cancer_C_4"
dataDEGs <- TCGAanalyze_DEA(mat1 = df1[,S1],mat2 = df2[,S2],Cond1type = paste0(p1,"_",C1),Cond2type = paste0(p2,"_",C2),fdr.cut = 0.01 ,logFC.cut = 1)
ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes ",p1,"_",C1," Vs ",p2,"_",C2), RegulonList = rownames(dataDEGs))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = rownames(dataDEGs),nBar = 20,filename = paste0("DEA_genes_",p1,"_",C1,"_vs_",p2,"_",C2,"_20.pdf"))
dataDEGs_h_c3_c4 <- dataDEGs
ansEA_h_c3_c4 <- ansEA

dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,paste0(p1,"_",C1),paste0(p2,"_",C2),df1[,S1],df2[,S2])
dataDEGsFiltLevel$GeneID <- 0
eg = as.data.frame ( bitr(dataDEGsFiltLevel$mRNA,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db"))
eg <- eg[!duplicated(eg$SYMBOL),]
dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL,]
dataDEGsFiltLevel <- dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA,decreasing=FALSE),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]
dataDEGsFiltLevel$GeneID <- eg$ENTREZID
dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID
hsa05200 <-  pathview(gene.data  = genelistDEGs,pathway.id = "hsa05200",species = "hsa",limit = list(gene=as.integer(max(abs(genelistDEGs)))))
file.rename("hsa05200.pathview.png",paste0("hsa05200_",p_h,"_",C1,"_",C2,".png"))
genelistDEGs_c3_c4 <- genelistDEGs

#low vs C4
dn <- TCGAanalyze_Normalization(tabDF = cbind(dp_h_c_low, dp_h_c_4),geneInfo = geneInfo)
df <- TCGAanalyze_Filtering(tabDF = dn,method = "quantile",qnt.cut =  0.25)
df1 <- subset(df, select = colnames(df) %in% colnames(df_h_c_low))
S1 <- TCGAquery_SampleTypes(colnames(df1),"TP")
C1 <- "Cancer_C_Low"
df2 <- subset(df, select = colnames(df) %in% colnames(df_h_c_4))
S2 <- TCGAquery_SampleTypes(colnames(df2),"TP")
C2 <- "Cancer_C_4"
dataDEGs <- TCGAanalyze_DEA(mat1 = df1[,S1],mat2 = df2[,S2],Cond1type = paste0(p1,"_",C1),Cond2type = paste0(p2,"_",C2),fdr.cut = 0.01 ,logFC.cut = 1)
ansEA <- TCGAanalyze_EAcomplete(TFname=paste0("DEA genes ",p1,"_",C1," Vs ",p2,"_",C2), RegulonList = rownames(dataDEGs))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = rownames(dataDEGs),nBar = 20,filename = paste0("DEA_genes_",p1,"_",C1,"_vs_",p2,"_",C2,"_20.pdf"))
dataDEGs_h_c_low_c4 <- dataDEGs
ansEA_h_c_low_c4 <- ansEA

dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,paste0(p1,"_",C1),paste0(p2,"_",C2),df1[,S1],df2[,S2])
dataDEGsFiltLevel$GeneID <- 0
eg = as.data.frame ( bitr(dataDEGsFiltLevel$mRNA,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db"))
eg <- eg[!duplicated(eg$SYMBOL),]
dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL,]
dataDEGsFiltLevel <- dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA,decreasing=FALSE),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]
dataDEGsFiltLevel$GeneID <- eg$ENTREZID
dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID
hsa05200 <-  pathview(gene.data  = genelistDEGs,pathway.id = "hsa05200",species = "hsa",limit = list(gene=as.integer(max(abs(genelistDEGs)))))
file.rename("hsa05200.pathview.png",paste0("hsa05200_",p_h,"_",C1,"_",C2,".png"))
genelistDEGs_c_low_c4 <- genelistDEGs


ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes of methylation regulated", RegulonList = all_genes_met_h)
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF, PathTab = ansEA$ResPat,nRGTab = all_genes_met_h,nBar = 20,filename = "DEA_genes_meth_reg.pdf")


save(genelistDEGs_n_c,genelistDEGs_n_c1,genelistDEGs_n_c2,genelistDEGs_n_c3,genelistDEGs_n_c4,genelistDEGs_n_c_low,
     genelistDEGs_c1_c2,genelistDEGs_c1_c3,genelistDEGs_c1_c4,genelistDEGs_c2_c3,genelistDEGs_c2_c4,genelistDEGs_c3_c4,
     genelistDEGs_c_low_c4,file=paste0(p_h,"_genelist.RData"))

save(ansEA_h_n_c,ansEA_h_n_c1,ansEA_h_n_c2,ansEA_h_n_c3,ansEA_h_n_c4,ansEA_h_n_c_low,ansEA_h_c1_c2,ansEA_h_c1_c3,ansEA_h_c1_c4,
     ansEA_h_c2_c3,ansEA_h_c2_c4,ansEA_h_c3_c4,ansEA_h_c_low_c4,file = paste0(p_h,"_ansEA.RData"))

save(dataDEGs_h_n_c,dataDEGs_h_n_c1,dataDEGs_h_n_c2,dataDEGs_h_n_c3,dataDEGs_h_n_c4,dataDEGs_h_n_c_low,dataDEGs_h_c1_c2,
     dataDEGs_h_c1_c3,dataDEGs_h_c1_c4,dataDEGs_h_c2_c3,dataDEGs_h_c2_c4,dataDEGs_h_c3_c4,dataDEGs_h_c_low_c4, file=paste0(p_h,"dataDEG.RData"))

save(rse_h_n,rse_h_c,rse_h_n_org,rse_h_c_org,rse_h_c_n,rse_h_c_1,rse_h_c_2,rse_h_c_3,rse_h_c_4,rse_h_c_low,rse_h_c_s,file = paste0(p_h,"_rse.RData"))

toc()


query_count <- GDCquery(project = p_h, 
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "HTSeq - Counts")
GDCdownload(query_count)
data <- GDCprepare(query_count, save=T, save.filename="TCGA.ESCA.Counts.raw.RData")


data$tumor_stage[data$tumor_stage %in% stage_1] <- "stage_i"
data$tumor_stage[data$tumor_stage %in% stage_2] <- "stage_ii"
data$tumor_stage[data$tumor_stage %in% stage_3] <- "stage_iii"
data$tumor_stage[data$tumor_stage %in% stage_4] <- "stage_iv"
data$tumor_stage[data$tumor_stage %in% stage_n] <- "stage_NA"
data$tumor_stage[data$barcode %in% rse_h_n$barcode] <- "Normal"
data$sample_type[data$barcode %in% rse_h_n$barcode] <- "Normal"
data$sample_type[data$barcode %in% rse_h_c$barcode] <- "Cancer"
data$sample_type <- as.factor(data$sample_type)
data$stage_f <- as.factor(data$tumor_stage)

dds <- DESeqDataSet(data[,!is.na(data$sample_type)],design = ~ sample_type)
dds <- DESeqDataSet(data[,!is.na(data$sample_type) & data$tumor_stage %in% c("stage_i","stage_ii")],design = ~ stage_f)

tic()
dds_1_2 <- DESeq(dds)
toc()


#-------------------------------------- 
# DNA methylation data
#-------------------------------------- 
# DNA methylation aligned to hg38
query_met.hg38 <- GDCquery(project= p_h, 
                           data.category = "DNA Methylation", 
                           platform = "Illumina Human Methylation 450")
GDCdownload(query_met.hg38)
data.hg38 <- GDCprepare(query_met.hg38)