#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
#BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.knownGene")
#BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")

library("rGADEM") 
library("MotIV") 
library("seqLogo") 
library("Biostrings") 
library("motifStack") 
library("BSgenome.Dmelanogaster.UCSC.dm6") 
library("GenomicRanges") 
library("muscle") 
library ("ChIPpeakAnno")
library(GenomicRanges)
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library(rtracklayer)
library("GenomicFeatures")
library("AnnotationDbi")
library("org.Dm.eg.db")
#library("TxDb.Dmelanogaster.UCSC.dm6.knownGene")


extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "integer")

bed_ecr_6aph=read.table("ecr_6aph.bed",header=FALSE,sep="\t")
bed_ecr_6aph=import("ecr_6aph.bed",extraCols=extraCols_narrowPeak)
bed_ecr_6aph=import("ecr_6aph.bed")

#annoData <- annoGR(TxDb.Dmelanogaster.UCSC.dm6.knownGene, feature="gene")
annoData <- annoGR(TxDb.Dmelanogaster.UCSC.dm6.ensGene, feature="gene")
#annotation type, could be "gene", "exon", "transcript", "CDS", "fiveUTR", "threeUTR", "microRNA", "tRNAs", "geneModel" 
#for TxDb object, or "gene", "exon" "transcript" for EnsDb object


#annotatedPeak_ecr_6aph <-annotatePeakInBatch(bed_ecr_6aph, AnnotationData=annoData)
#annotatedPeak_ecr_6aph <- annotatePeakInBatch(bed_ecr_6aph, AnnotationData=annoData,output="overlapping")
annotatedPeak_ecr_6aph <- annotatePeakInBatch(bed_ecr_6aph, AnnotationData=annoData,output="both")
            
#try with gz file

annotatedPeak_ecr_6aph$symbol <- mapIds(org.Dm.eg.db, 
                     keys=annotatedPeak_ecr_6aph$"feature", 
                     column="SYMBOL", 
                     keytype="FLYBASE",
                     multiVals="first")


annotatedPeak_ecr_6aph$entrez <- mapIds(org.Dm.eg.db, 
                            keys=annotatedPeak_ecr_6aph$"feature", 
                            column="ENTREZID", 
                            keytype="FLYBASE",
                            multiVals="first")

annotatedPeak_ecr_6aph$name2 =   mapIds(org.Dm.eg.db,
                           keys=annotatedPeak_ecr_6aph$"feature", 
                           column="GENENAME",
                           keytype="FLYBASE",
                           multiVals="first")


annotatedPeak_ecr_6aph<-annotatedPeak_ecr_6aph[order(annotatedPeak_ecr_6aph$"symbol",decreasing=F), ]
write.table(annotatedPeak_ecr_6aph, "annotatedPeak_ecr_6aph.csv",sep=";",row.names=F)

df_ecr_6aph<-data.frame(table(unlist(annotatedPeak_ecr_6aph$symbol)))
df_ecr_6aph<-df_ecr_6aph[order(df_ecr_6aph[,2],decreasing=T), ]
write.table(df_ecr_6aph, "df_ecr_6aph.csv")

