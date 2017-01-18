rm(list=ls())
gentype <- read.delim("D:/ivan/ivan6-cell/Results/downloadedData/gentype.txt")
colnames(gentype)=c("gene","type")
protein_coding=gentype[gentype$type=="protein_coding",]
lincRNA=gentype[gentype$type=="lincRNA",]
pRNA=gentype[gentype$type=="pRNA",]
