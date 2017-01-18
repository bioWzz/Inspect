rm(list=ls())
library(ggplot2)
library(sqldf)
library(dplyr)

# 读入labA的结果数据
load("D:/ivan/ivan6-cell/Results/Result_getfromthePROGRAMME/LabAandCresults/LabAGene.RData")

# Delete the no use varies 
vari=objects()
delva=as.vector(grep("^La",vari,invert=TRUE))  
vari=vari[delva]
rm(list=vari)

# 读入管家基因数据,并将管家基因id转化为ENSG
Ense_genesymbol <- read.delim("D:/ivan/ivan6-cell/Results/downloadedData/Ense_genesymbol.txt", stringsAsFactors=FALSE)
colnames(Ense_genesymbol )=c("ENSG","id")

HK_gene <- read.delim("D:/ivan/ivan6-cell/Results/downloadedData/HK_genes.txt", header=FALSE, stringsAsFactors=FALSE)
colnames(HK_gene)=c("id","nebi")
HK_genes<- left_join(HK_gene,Ense_genesymbol,by="id")

HK_genes=as.data.frame(unique(HK_genes[,3]));
HK_genes=na.omit(HK_genes)  
colnames(HK_genes)=c("id")

#从结果数据中筛选出管家基因的数据
# 由于inspect 的ENSG基因的id与 管家基因的不同需要进行转化。例如：ENSE12 EBSEG12.1
name=as.data.frame(rownames(LabAREP1Express))
name1=as.data.frame(strsplit(as.character(name[,1]),'\\.'))
name1<- t(name1)
name1=as.data.frame(name1[,1])
colnames(name1)[1]="id"
LabAREP1Express=cbind(name1,LabAREP1Express)

# 筛选出管家基因的数据
HKgeneResult<-left_join(HK_genes,LabAREP1Express,by="id")
rownames(HKgeneResult)=HKgeneResult[,1]
HKgeneResult=HKgeneResult[,2:16]
HKgeneResult[apply(HKgeneResult<=0, FUN = any, 1), ] = NA
HKgeneResult=na.omit(HKgeneResult)  


# read the viral and viral mRNA data
LabARep1 <- read.table("D:/ivan/ivan6-cell/Results/Result_getfromthePROGRAMME/virsua-Result/LabARep1.txt", quote="\"", comment.char="")
virlabARep1=LabARep1[1:8,]
virmRNlabARep1=LabARep1[9:16,]

# 画出在个时间点最大合成率的分布图
# LabA-Gene
# rep1
LabAREP1Express1=HKgeneResult
# synthsisLabAREP1=apply(LabAREP1Express1[,7:9],1,max)
synthsisLabAREP1=LabAREP1Express1[,7:9]
synthsisLabAREP1=as.data.frame(as.numeric(unlist(synthsisLabAREP1)))
synthsisLabAREP1=log(synthsisLabAREP1)

virlabARep1=LabARep1[1:8,]
# virlabARep1synthsis=apply(virlabARep1[,1:3],1,max)
virlabARep1synthsis=virlabARep1[,1:3]
virlabARep1synthsis=as.data.frame(as.numeric(unlist(virlabARep1synthsis)))
virlabARep1synthsis=log(virlabARep1synthsis)

graphics.off()
data <- data.frame(x = virlabARep1synthsis)
x=as.numeric(unlist(data)) 
hist(x)


# data1<- data.frame(x = pre_RNA_4SU)
# 生成底层和直方图,概率线的图层
p <- ggplot(data, aes(x = synthsisLabAREP1, y = ..density..))
# p <- p + geom_histogram(fill = "navy")
layer1 <- geom_density(data=data,colour = "blue",size = 1)
(p1 <- p + layer1)
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
layer2 <- geom_density(data=virlabARep1synthsis,colour = "red",linetype = 2,size = 1)
(p1 <- p1 + layer2 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

# 画出在个时间点最大降解率的分布图
degradationLabAREP1=apply(LabAREP1Express1[,10:12],1,max)
degradationLabAREP1=as.data.frame(as.numeric(unlist(degradationLabAREP1)))
degradationLabAREP1=log(degradationLabAREP1)
graphics.off()
data <- data.frame(x = degradationLabAREP1)
# data1<- data.frame(x = pre_RNA_4SU)
# 生成底层和直方图,概率线的图层
p <- ggplot(data, aes(x =degradationLabAREP1, y = ..density..))
# p <- p + geom_histogram(fill = "navy")
layer1 <- geom_density(data=data,colour = "blue",size = 1)
(p1 <- p + layer1)
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
# layer2 <- geom_density(data=data1,colour = "blue",linetype = 2,size = 1)
# (p1 <- p1 + layer2 + xlab("log(Expression)") + ylab("fraction of junctions"))
# p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
