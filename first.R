a = list.files("D:/ivan/Harm/NewTestFour/my/R/") 
b=paste("D:/ivan/Harm/NewTestFour/my/R/",a,sep ="")
for (i in b )source(i)


a = list.files("D:/ivan/Harm/NewTestFour/my/RR/") 
b=paste("D:/ivan/Harm/NewTestFour/my/RR/",a,sep ="")
for (i in b )source(i)

library(GenomicFeatures)
library(INSPEcT)
library(rtracklayer)
library(INSPEcT)
library(biomaRt)
library(GenomicFeatures)
library(rtracklayer)
library(RMySQL)
library(RSQLite)
library(AnnotationDbi)
library(compiler)
library(GenomicAlignments)
library(Rsamtools)
#source("https://bioconductor.org/biocLite.R")
#biocLite("INSPEcT")
#dir="D:/ivan/Harm/NewTestFour/my/R/"
#a = list.files(dir) 
#b=paste("D:/ivan/Harm/NewTestFour/my/R/",a,sep ="")
#for (i in b )source(i)
#for(i in 1:23)
#  print(paste("source(b[",i,"])",sep = ""))

#chrominfo <- data.frame(chrom = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrx','chry'),
#                       length=c(249250621,),
#                       is_circular=c(FALSE))

file1="D:/ivan/Harm/NewTestOne/NewTest/gencode.v23.PR81.gtf"

txdb<- makeTxDbFromBiomart(dataset="hsapiens_gene_ensembl")
txdb3 <- makeTxDbFromGFF1(file="D:/ivan/Harm/NewTestOne/NewTest/gencode.v23.PR8.gtf",format="gtf",taxonomyId=9606)

transcripts=load("d:/transcripts.Rdata")
splicings=load("d:/transcripts.Rdata")
genes=load("d:/genes.Rdata")
metadata=load("d:/metadata.Rdata")
txdb <- makeTxDb(transcripts, splicings,genes=genes)
#exon_txdb=exons(txdb)
#tmp=as.data.frame(exon_txdb)
#genes_txdb=genes(txdb)
#tmp=as.data.frame(genes_txdb)
#seqlengths(exon_txdb)
#tr1=as.data.frame(transcripts(txdb))
path1 <- "D:/ProgramFiles/R/R-3.3.1/library/INSPEcT/extdata/4sURNA_0h.bam"
path2<- "D:/ProgramFiles/R/R-3.3.1/library/INSPEcT/extdata/totalRNA_0h.bam"
paths_4su<-c(path1,path1,path1)
paths_foursu<-c(path2,path2,path2) 
by='gene'
countMultiMappingReads=FALSE
allowMultiOverlap=FALSE
strandSpecific=FALSE
isPairedEnd=FALSE

s


makeRPKMsOut <- makeRPKMs1(txdb, paths_4su, paths_total,by='tx')
rpkms <- makeRPKMsOut$rpkms
counts <- makeRPKMsOut$counts
annotation <- makeRPKMsOut$annotation

tpts <- c(0,4,8)
tL <- 4
mycerIds <- newINSPEcT(tpts, tL, rpkms$foursu_exons, rpkms$total_exons,
                       rpkms$foursu_introns, rpkms$total_introns, BPPARAM=SerialParam())

