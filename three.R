library(sqldf)
require(downloader)  
library(dplyr)
library(sqldf)
library(data.table)
library(ggplot2)
library(compare)
library(plotrix)
rm(list=ls())
path="D:/ivan/Harm/NewTestFour/my/"
Gene_Symbol=read.table(paste(path,"Gene_Tr.csv",sep=""), header = TRUE, sep = ",",row.names=NULL)
colnames(Gene_Symbol) <- c("ID","HGNC")

RNA_4sU_LabA=read.table(paste(path,"RNA_4sU_LabA.txt",sep=""), header = TRUE, sep = " ",row.names=NULL)
colnames(RNA_4sU_LabA)

newdf <- sqldf("select Gene_Symbol.HGNC, from Gene_Symbol join RNA_4sU_LabA where Gene_Symbol$ID=RNA_4sU_LabA$GENE")
write.table(newdf, file=paste(path,"DESeq2_SIUBW2I_sictrl_UT_compare_foldchange2.txt",sep=""),row.names=FALSE, sep="\t", quote=FALSE)




path="D:/ivan/Harm/NewTestFour/my/"
Gene_Symbol=read.csv("D:/ivan/Harm/NewTestFour/my/Gene_Tr.csv", header = TRUE, sep = ",", quote = "\"",dec = ".", fill = TRUE, comment.char = "")
RNA_4sU_LabA <- read.csv("D:/ivan/Harm/NewTestFour/my/RNA_4sU_LabA.txt", sep="")
RNA_4sU_LabAGeneName=substr(RNA_4sU_LabA[,1],1,15);RNA_4sU_LabA[,1]=RNA_4sU_LabAGeneName;RNA_4sU_LabA=as.data.frame(RNA_4sU_LabA)
colnames(RNA_4sU_LabA)

newdf <- sqldf("select RNA_4sU_LabA.GENE,Gene_Symbol.HGNC,RNA_4sU_LabA.Onehour,RNA_4sU_LabA.Fourhour,RNA_4sU_LabA.Eighthour from RNA_4sU_LabA join Gene_Symbol where select RNA_4sU_LabA.GENE=Gene_Symbol.ID")


RNA_total_UnlA <- read.csv("D:/ivan/Harm/NewTestFour/my/RNA_total_UnlA.txt", sep="")
RNA_total_UnlAGeneName=substr(RNA_total_UnlA[,1],1,15);RNA_total_UnlA[,1]=RNA_total_UnlAGeneName;RNA_total_UnlA=as.data.frame(RNA_total_UnlA)
