rm(list=ls())
library(sqldf)
path="D:/ivan/Harm/NewTestFour/my/Results/Fianllyresults/virusforLabAsndLabc/"
file=list.files(path)
for(i in file){
  assign(paste(substr(i,12,nchar(i)-4),sep=""),read.table(paste(path,i,sep=""), header = TRUE, sep = "\t",row.names=NULL))    
}


Virtal_4sU_LabA<- sqldf("select * from RNA_4sU_LabA join RNA_virtalgne where RNA_4sU_LabA.Symbol= RNA_virtalgne.Symbol")
Virtal_4sU_LabC<- sqldf("select * from RNA_4sU_LabC join RNA_virtalgne where RNA_4sU_LabC.Symbol= RNA_virtalgne.Symbol")

Virtal_total_UnlA<- sqldf("select * from RNA_total_UnlA join RNA_virtalgne where RNA_total_UnlA.Symbol= RNA_virtalgne.Symbol")
Virtal_total_UnlC<- sqldf("select * from RNA_total_UnlC join RNA_virtalgne where RNA_total_UnlC.Symbol= RNA_virtalgne.Symbol")

virtal_4su_AC=cbind(Virtal_4sU_LabA,Virtal_4sU_LabC)[,c(1,2,3,4,7,8,9)]
colnames(virtal_4su_AC)<-c("Symbol","LabAHour0","LabAHour4","LabAHour8","LabCHour0","LabCHour4","LabCHour8")
write.table(virtal_4su_AC, file="D:/ivan/Harm/NewTestFour/my/Results/Fianllyresults/virusforLabAsndLabc/virtual_4sU_expression.txt",row.names=FALSE, sep="\t", quote=FALSE)

virtal_total_AC=cbind(Virtal_total_UnlA,Virtal_total_UnlC)[,c(1,2,3,4,7,8,9)]
colnames(virtal_total_AC)<-c("Symbol","TotalAHour0","TotalAHour4","TotalAHour8","TotalAHour0","TotalAHour4","TotalAHour8")
write.table(virtal_total_AC, file="D:/ivan/Harm/NewTestFour/my/Results/Fianllyresults/virusforLabAsndLabc/virtual_total_expression.txt",row.names=FALSE, sep="\t", quote=FALSE)


#利用edgeR来计算差异表达，数据来自于harm D:\ivan\Harm\NewTestFour\my\Results\Fianllyresults\virusforLabAsndLabc，包含了
#病毒基因的表达和病毒mRNA的表达
rm(list=ls())
path="D:/ivan/Harm/NewTestFour/my/Results/Fianllyresults/virusforLabAsndLabc/"
library(edgeR)
library(splines)
#load(url("http://qiubio.com:8080/bioconductor/RNA-seq/ds1.Rdata"))
#write.table(counts, file="D:/ivan/Harm/NewTestFour/my/Results/Fianllyresults/virusforLabAsndLabc/COUNT.txt",row.names=FALSE, sep="\t", quote=FALSE)

rawdata=read.delim(paste(path,"EXOSC3_RNA_labeling_all_counts.txt",sep=""), check.names=FALSE, stringsAsFactors=FALSE)
GeneLength=read.delim(paste(path,"FgeneLenth.txt",sep=""), check.names=FALSE, stringsAsFactors=FALSE)
GeneLength=GeneLength[,3]
group <- c("LA0","LA4","LA8","LC0","LC4","LC8","UA0","UA4","UA8","UC0","UC4","UC8","LA0","LA4","LA8","LC0","LC4","LC8","UA0","UA4","UA8","UC0","UC4","UC8")
d<- DGEList(counts=rawdata[,2:25], genes=data.frame(rawdata[,1],Length=GeneLength),group=group)
d <- calcNormFactors(d)
RPKM1 <- rpkm(d)
write.table(RPKM1, file="D:/ivan/Harm/NewTestFour/my/Results/Fianllyresults/virusforLabAsndLabc/RPKM1.txt",row.names=FALSE, sep="\t", quote=FALSE)


