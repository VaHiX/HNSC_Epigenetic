
# Query platform Illumina HiSeq with a list of barcode 
query <- GDCquery(project = p_h, 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  legacy = TRUE)

# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(query)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
RnaseqSE <- GDCprepare(query)

Matrix <- assay(RnaseqSE,"raw_count") # or BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")

# For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
Rnaseq_CorOutliers <- TCGAanalyze_Preprocessing(RnaseqSE)


query <- GDCquery(project = p_h,
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  platform = "Illumina HiSeq", 
                  file.type  = "normalized_results",
                  experimental.strategy = "RNA-Seq",
                  legacy = TRUE)
GDCdownload(query, method = "api")
data <- GDCprepare(query)

# Gene expression aligned against hg38
query <- GDCquery(project = p_h,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM-UQ")
GDCdownload(query)
data <- GDCprepare(query)

ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)


query <- GDCquery(project = p_h,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")

samplesDown <- getResults(query,cols = c("cases"))

dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP")

dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "NT")
dataSmTP_short <- dataSmTP[1:10]
dataSmNT_short <- dataSmNT[1:10]

queryDown <- GDCquery(project = p_h, 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - Counts", 
                      barcode = c(dataSmTP_short, dataSmNT_short))

GDCdownload(query = queryDown)


dataPrep1 <- GDCprepare(query = queryDown, 
                        save = TRUE, 
                        save.filename = "TCGA_ESCA_HTSeq_Countds.rda")

dataPrep <- TCGAanalyze_Preprocessing(object = dataPrep1, 
                                      cor.cut = 0.6,
                                      datatype = "HTSeq - Counts")


dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent")
dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "geneLength")

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)
