#Change the ENSG id to GeneSymbol
#
rm(list=ls())
library(sqldf)
path="D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep2LabA/"
path1="D:/ivan/Harm/NewTestFour/my/"
file=list.files("D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep2LabA/")

genesymbol=read.table(paste(path1,"Gene_symbol.txt",sep=""), header = TRUE, sep = "\t",row.names=NULL)
genesymbol=unique(genesymbol)
colnames(genesymbol)<-c("Number","ID","Symbol")



for(i in file){
  Results=read.table(paste(path,i,sep=""), header = TRUE, sep = " ",row.names=NULL)
  Results[,1]=substr(Results[,1],1,15)
  colnames(Results)<-c("ID","Hour0","Hour4","Hour8")
  print(paste(path,i,sep=""))
  print(dim(Results))
  genesymbol$ID=as.character(genesymbol$ID)
  newdf <- sqldf("select genesymbol.Symbol,Results.Hour0,Results.Hour4,Results.Hour8 from Results join genesymbol where genesymbol.ID=Results.ID")
  write.table(newdf, file=paste(path,"GeneSymbol_",i,sep=""),row.names=FALSE, sep="\t", quote=FALSE)
  
}

sy=read.table("D:/result.txt", header = TRUE, sep = "\t",row.names=NULL)

re <- sqldf("select * from sy join genesymbol where genesymbol.ID=sy.gene_id")
write.table(re, file="D:/resultsre.txt",row.names=FALSE, sep="\t", quote=FALSE)

#Merge the data of expreesion and rates into a intergated files
#
#The result file is all.txt
rm(list=ls())
path2="D:/ivan/Harm/NewTestFour/my/TestForDegrate_2016.6.19/REP1/"
file=list.files(path2)
for(i in file){
  assign(paste(substr(i,1,nchar(i)-4),sep=""),read.table(paste(path2,i,sep=""), header = TRUE, sep = "\t",row.names=NULL))    
}
pre_RNA_4sU_LabA=GeneSymbol_REP2_pre_RNA_4sU_LabA
pre_RNA_total_UnlA=GeneSymbol_REP2_pre_RNA_total_UnlA
RNA_4sU_LabA=GeneSymbol_REP2_Total_RNA_4sU_LabA
RNA_total_UnlA=GeneSymbol_REP2_Total_RNA_total_UnlA
LabADegradation=GeneSymbol_Rep2LabADegradation
LabASynthesis=GeneSymbol_Rep2LabASynthesis
LabAProcessing=GeneSymbol_Rep2LabAProcessing

#pre_RNA_4sU_LabA=REP1_pre_RNA_4sU_LabA2
#pre_RNA_total_UnlA=GeneSymbol_REP1_pre_RNA_total_UnlC
#RNA_4sU_LabA=GeneSymbol_REP1_Total_RNA_4sU_LabC
#RNA_total_UnlA=GeneSymbol_REP1_Total_RNA_total_UnlC
#LabADegradation=GeneSymbol_Rep1LabCDegradation
#LabASynthesis=GeneSymbol_Rep1LabCSynthesis
#LabAProcessing=GeneSymbol_REP1_LabCProcessing


colnames(RNA_4sU_LabA)<-c("S","A0","A4","A8")
colnames(RNA_total_UnlA)<-c("St","UA0","UA4","UA8")
colnames(pre_RNA_4sU_LabA)<-c("prS","prA0","prA4","prA8")
colnames(pre_RNA_total_UnlA)<-c("ppS","rpUA0","prUA4","prUA8")

En=cbind(RNA_4sU_LabA,RNA_total_UnlA)
En=En[,c(1,2,3,4,6,7,8)]
Pn=cbind(pre_RNA_4sU_LabA,pre_RNA_total_UnlA)
Pn=Pn[,c(1,2,3,4,6,7,8)]

DE=LabADegradation
colnames(DE)<-c("dS","dA0","dA4","dA8")
PR=LabAProcessing
colnames(PR)<-c("pS","pA0","pA4","pA8")
Sy=LabASynthesis
colnames(Sy)<-c("Sy","sA0","sA4","sA8")

AEN<- sqldf("select * from En join Pn where En.S= Pn.prS")
ADE<- sqldf("select * from AEN join DE where AEN.S= DE.dS")
APR<-sqldf("select * from ADE join PR where PR.pS= ADE.dS")
Ssy<-sqldf("select * from APR join Sy where APR.S= Sy.Sy")

write.table(Ssy, file="D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep2LabA/Gene/Rep2_LabA.txt",row.names=FALSE, sep="\t", quote=FALSE)



#Merge the data of expreesion and rates into a intergated files Labc_Rep1
#
#The result file is all.txt
rm(list=ls())
path2="D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep1LabC/Gene/"
file=list.files(path2)
for(i in file){
  assign(paste(substr(i,1,nchar(i)-4),sep=""),read.table(paste(path2,i,sep=""), header = TRUE, sep = "\t",row.names=NULL))    
}

#pre_RNA_4sU_LabA=GeneSymbol_REP1_pre_RNA_4sU_LabC
#pre_RNA_total_UnlA=GeneSymbol_REP1_pre_RNA_total_UnlC
RNA_4sU_LabA=GeneSymbol_REP1RNA_4sU_LabC
RNA_total_UnlA=GeneSymbol_REP1RNA_total_UnlC
LabADegradation=GeneSymbol_REP1_LabCDegradation
LabASynthesis=GeneSymbol_REP1_LabCSynthesis
LabAProcessing=GeneSymbol_REP1_LabCProcessing


colnames(RNA_4sU_LabA)<-c("S","A0","A4","A8")
colnames(RNA_total_UnlA)<-c("St","UA0","UA4","UA8")
colnames(pre_RNA_4sU_LabA)<-c("prS","prA0","prA4","prA8")
colnames(pre_RNA_total_UnlA)<-c("ppS","rpUA0","prUA4","prUA8")

En=cbind(RNA_4sU_LabA,RNA_total_UnlA)
En=En[,c(1,2,3,4,6,7,8)]
Pn=cbind(pre_RNA_4sU_LabA,pre_RNA_total_UnlA)
Pn=Pn[,c(1,2,3,4,6,7,8)]

DE=LabADegradation
colnames(DE)<-c("dS","dA0","dA4","dA8")
PR=LabAProcessing
colnames(PR)<-c("pS","pA0","pA4","pA8")
Sy=LabASynthesis
colnames(Sy)<-c("Sy","sA0","sA4","sA8")

AEN<- sqldf("select * from En join Pn where En.S= Pn.prS")
ADE<- sqldf("select * from AEN join DE where AEN.S= DE.dS")
APR<-sqldf("select * from ADE join PR where PR.pS= ADE.dS")
Ssy<-sqldf("select * from APR join Sy where APR.S= Sy.Sy")

write.table(Ssy, file="D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep1LabC/Gene/Rep1_LabC.txt",row.names=FALSE, sep="\t", quote=FALSE)




#####Finally Rerults Intergrate
rm(list=ls())
library(sqldf)
path="D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Integer_Results/"

Rep1_LabCResults=read.table(paste(path,"Rep1_LabA.txt",sep=""), header = TRUE, sep = "\t",row.names=NULL)

Rep2_LabCResults=read.table(paste(path,"Rep2_LabA.txt",sep=""), header = TRUE, sep = "\t",row.names=NULL) 
Ssy<-sqldf("select * from Rep1_LabCResults join Rep2_LabCResults where Rep1_LabCResults.S= Rep2_LabCResults.S")
write.table(Ssy, file="D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Integer_Results/LabA.txt",row.names=FALSE, sep="\t", quote=FALSE)


LabA_P=read.table(paste(path,"LabA-P.txt",sep=""), header = TRUE, sep = "\t",row.names=NULL)
LabA_P<-na.omit(LabA_P)
write.table(LabA_P, file="D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Integer_Results/LabA_P_NoNA.txt",row.names=FALSE, sep="\t", quote=FALSE)


LabA_D=read.table(paste(path,"LabA-d.txt",sep=""), header = TRUE, sep = "\t",row.names=NULL)
LabA_D<-na.omit(LabA_D)
write.table(LabA_D, file="D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Integer_Results/LabA_D_NoNA.txt",row.names=FALSE, sep="\t", quote=FALSE)

LabC_P=read.table(paste(path,"LabC-P.txt",sep=""), header = TRUE, sep = "\t",row.names=NULL)
LabC_P<-na.omit(LabC_P)
write.table(LabC_P, file="D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Integer_Results/LabC_P_NoNA.txt",row.names=FALSE, sep="\t", quote=FALSE)


LabC_D=read.table(paste(path,"LabC-d.txt",sep=""), header = TRUE, sep = "\t",row.names=NULL)
LabC_D<-na.omit(LabC_D)
write.table(LabC_D, file="D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Integer_Results/LabC_D_NoNA.txt",row.names=FALSE, sep="\t", quote=FALSE)

rm(list=ls())
Rep1_LabA=read.table("D:/ivan/Harm/NewTestFour/my/TestForDegrate_2016.6.19/Rep1LabARate2.txt", header = TRUE, sep = " ",row.names=NULL)
colnames(Rep1_LabA)<-c("GENE","synthesis_0","synthesis_4","synthesis_8","degradation_0","degradation_4","degradation_8","processing_0","processing_4","processing_8")
Rep2_LabA=read.table("D:/ivan/Harm/NewTestFour/my/TestForDegrate_2016.6.19/Rep2LabARate.txt", header = TRUE, sep = " ",row.names=NULL)
colnames(Rep2_LabA)<-c("GENE","synthesis_0","synthesis_4","synthesis_8","degradation_0","degradation_4","degradation_8","processing_0","processing_4","processing_8")
library(sqldf)
AEN<- sqldf("select * from Rep1_LabA join Rep2_LabA where Rep2_LabA.GENE= Rep1_LabA.GENE")
AEN<-na.omit(AEN)
write.table(AEN, file="D:/ivan/Harm/NewTestFour/my/TestForDegrate_2016.6.19/all.txt",row.names=FALSE, sep="\t", quote=FALSE)
