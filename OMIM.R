rm(list=ls())
library(sqldf)
path2="D:/wzz/OMIM/DATA/"
file=list.files(path2)
for(i in file){
  assign(paste(substr(i,1,nchar(i)-4),sep=""),read.table(paste(path2,i,sep=""), header = TRUE, sep = "\t",row.names=NULL))    
}

#mim2gene <- read.delim("D:/wzz/OMIM/DATA/mim2gene.txt")
#colnames(mim2gene)<-c("DieaseNumber","Type","EntrezID","Symbol","Ensembl")
#mim2gene=unique(mim2gene[mim2gene$Type=="phenotype"|mim2gene$Type=="predominantly phenotypes",c(1,2)])

DieaseName <- read.delim2("D:/wzz/OMIM/DATA/DieaseName.txt")
colnames(DieaseName)<-c("DieaseName","GeneNmuber")

znf <- read.delim2("D:/wzz/OMIM/DATA/znf.txt")
colnames(znf)<-c("GeneNumber","g1","g2","g3")
result<- sqldf("select * from DieaseName join znf where DieaseName.GeneNmuber= znf.GeneNumber")
result<-result[,c(1,4,5,6)]

write.table(result, file="D:/wzz/OMIM/DATA/result.txt",row.names=FALSE, sep="\t", quote=FALSE)
