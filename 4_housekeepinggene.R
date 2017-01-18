rm(list=ls())
load("D:/ivan/ivan6-cell/Results/结果数据/LabAGene.RData")

# Delete the no use varies 
vari=objects()
delva=as.vector(grep("^La",vari,invert=TRUE))  
vari=vari[delva]
rm(list=vari)

HK_genes <- read.delim("D:/ivan/ivan6-cell/Results/downloadedData/HK_genes.txt", header=FALSE)
HK_genes=as.data.frame(HK_genes[,1])
colnames(HK_genes)=c("genesymbol")

# LabA-Gene
# rep1
LabAREP1Express1=LabAREP1Express
synthsisLabAREP1=LabAREP1Express1[,7:9]
degradation