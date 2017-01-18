# D:\ivan\ivan6-cell\Results\6-Figure 3-wzz\D-E
# save.image("D:/ivan/ivan6-cell/Results/Result_getfromthePROGRAMME/varbetweenGandT.RData")

rm(list=ls())
library(ggplot2)

# Delete the no use varies 
# vari=objects()
# delva=as.vector(grep("^Sy|^pr|^de",vari,invert=TRUE))  
# vari=vari[delva]
# rm(list=vari)

load("D:/ivan/ivan6-cell/Results/Result_getfromthePROGRAMME/LabAandCresults/LabAandC.RData")

varname<- c("GeneLabAREP1","GeneLabAREP12","GeneLabAREP2","GeneLabCREP1","GeneLabCREP12",
            "GeneLabCREP2","TXLabAREP1","TXLabAREP12","TXLabAREP2","TXLabCREP1","TXLabCREP12",
            "TXLabCREP2")
#将ENSG的id进行转换
for (var in varname)
  assign(var,ConvertName(get(var)))

# 将ENSG的id及数据放置在ENST的后面
varname<- list(c("TXLabAREP1","GeneLabAREP1"),c("TXLabAREP2","GeneLabAREP2"),c("TXLabAREP12","GeneLabAREP12"),
              c("TXLabCREP1","GeneLabCREP1"),c("TXLabCREP2","GeneLabCREP2"),c("TXLabCREP12","GeneLabCREP12"))

for (var in varname)
{ 
  TXExpression=var[1]
  GeneExpression=var[2]
  assign(TXExpression,AddEnstToEnsg(get(TXExpression),get(GeneExpression)))
}

###################################################################################################
# 转录率处理
SynthsisTXLabAREP1=as.data.frame(varcal(DeleteSynthsisNAdata(TXLabAREP1[,c(1,2,9,10,11,24,25,26)])))
SynthsisTXLabAREP2=as.data.frame(varcal(DeleteSynthsisNAdata(TXLabAREP2[,c(1,2,9,10,11,24,25,26)])))
SynthsisTXLabAREP12=as.data.frame(varcal(DeleteSynthsisNAdata(TXLabAREP12[,c(1,2,9,10,11,24,25,26)])))
SynthsisTXLabCREP1=as.data.frame(varcal(DeleteSynthsisNAdata(TXLabCREP1[,c(1,2,9,10,11,24,25,26)])))
SynthsisTXLabCREP2=as.data.frame(varcal(DeleteSynthsisNAdata(TXLabCREP2[,c(1,2,9,10,11,24,25,26)])))
SynthsisTXLabCREP12=as.data.frame(varcal(DeleteSynthsisNAdata(TXLabCREP12[,c(1,2,9,10,11,24,25,26)])))
# 降解率处理
degrateTXLabAREP1=as.data.frame(varcal(DeleteSynthsisNAdata(TXLabAREP1[,c(1,2,12,13,14,27,28,29)])))
degrateTXLabAREP2=as.data.frame(varcal(DeleteSynthsisNAdata(TXLabAREP2[,c(1,2,12,13,14,27,28,29)])))
degrateTXLabAREP12=as.data.frame(varcal(DeleteSynthsisNAdata(TXLabAREP12[,c(1,2,12,13,14,27,28,29)])))
degrateTXLabCREP1=as.data.frame(varcal(DeleteSynthsisNAdata(TXLabCREP1[,c(1,2,12,13,14,27,28,29)])))
degrateTXLabCREP2=as.data.frame(varcal(DeleteSynthsisNAdata(TXLabCREP2[,c(1,2,12,13,14,27,28,29)])))
degrateTXLabCREP12=as.data.frame(varcal(DeleteSynthsisNAdata(TXLabCREP12[,c(1,2,12,13,14,27,28,29)])))

# 处理率处理
processingTXLabAREP1=as.data.frame(varcal(DeleteSynthsisNAdata(TXLabAREP1[,c(1,2,15,16,17,30,31,32)])))
# processingTXLabAREP2=as.data.frame(varcal(DeleteSynthsisNAdata(TXLabAREP2[,c(1,2,15,16,17,30,31,32)])))
# processingTXLabAREP12=as.data.frame(varcal(DeleteSynthsisNAdata(TXLabAREP12[,c(1,2,15,16,17,30,31,32)])))
# processingTXLabCREP1=as.data.frame(varcal(DeleteSynthsisNAdata(TXLabCREP1[,c(1,2,15,16,17,30,31,32)])))
# processingTXLabCREP2=as.data.frame(varcal(DeleteSynthsisNAdata(TXLabCREP2[,c(1,2,15,16,17,30,31,32)])))
# processingTXLabCREP12=as.data.frame(varcal(DeleteSynthsisNAdata(TXLabCREP12[,c(1,2,15,16,17,30,31,32)])))

# 画GENE-TR变异系数的概率分布图

varbetweenGandT(SynthsisTXLabAREP1,degrateTXLabAREP1,1)
varbetweenGandT(SynthsisTXLabAREP2,degrateTXLabAREP2,2)
varbetweenGandT(SynthsisTXLabAREP12,degrateTXLabAREP12,3)

varbetweenGandT(SynthsisTXLabCREP1,degrateTXLabCREP1,4)
varbetweenGandT(SynthsisTXLabCREP2,degrateTXLabCREP2,5)
varbetweenGandT(SynthsisTXLabCREP12,degrateTXLabCREP12,6)

# 画TR-TR变异系数的概率分布图
varbetweenGandT(SynthsisTXLabAREP1,degrateTXLabAREP1,1)




varbetweenGandT<-function(sy,de,i){

  # synthsis
  synthsis=sy[,3]
  synthsis=as.data.frame(as.numeric(unlist(synthsis)))

  # processing
  # processing=Express[,13:15]
  # processing=as.data.frame(as.numeric(unlist(processing)))

  # degration
  degration=de[,3]
  degration=as.data.frame(as.numeric(unlist(degration)))
  

  
  pre_RNA_total=as.data.frame(as.numeric(unlist(synthsis))) 
  pre_RNA_total=pre_RNA_total[pre_RNA_total>0]
  # pre_RNA_total=log(pre_RNA_total)
  
  pre_RNA_4SU=as.data.frame(as.numeric(unlist(degration)))
  pre_RNA_4SU=pre_RNA_4SU[pre_RNA_4SU>0]
  # pre_RNA_4SU=log(pre_RNA_4SU)
  
  # 画pre_RNA_4SU和pre_RNA_total的概率分布图
  graphics.off()
  data <- data.frame(x = pre_RNA_total)
  data1<- data.frame(x = pre_RNA_4SU)
  # 生成底层和直方图,概率线的图层
  p <- ggplot( data=data,aes(x = x, y = ..density..))
  # p <- p + geom_histogram(fill = "navy")
  layer1 <- geom_density(data=data,colour = "black",size = 1)
  (p1 <- p + layer1+ylim(0,4))
  p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
  layer2 <- geom_density(data=data1,colour = "green",size = 1)
  (p1 <- p1 + layer2 + xlab("Fraction of variance explained)") + ylab("fraction of genes")+ scale_colour_manual("PID",values = c("black" = "black"," gre" = "blue")))
  p1+ theme( panel.background = element_blank(),axis.line = element_line(colour = "black"))
 
  ggsave(paste("D:/ivan/ivan6-cell/Results/6-Figure3-wzz-A-Cfinish_Ddoing/",i, '.png', sep=""),width=5,height=5,dpi=300)
}




#convert ENSG id to ENSG id eg. convert ensg12.1 to ensg12  
ConvertName<-function(Express){
  name=as.data.frame(rownames(Express))
  name1=as.data.frame(strsplit(as.character(name[,1]),'\\.'))
  name1<- t(name1)
  name1=as.data.frame(name1[,1])
  colnames(name1)[1]="id"
  Express=cbind(name1,Express)
  rownames(Express)=Express[,1]
  Express=Express[,2:16]
  return(Express)
}

###################################################################################################
#add ENSGid to ENST tables 

AddEnstToEnsg<-function(TXExpress,GeneExpression){
# ENSG.ENST <- read.delim("D:/ivan/ivan6-cell/Results/6-Figure 3-wzz/D-E/ENSG-ENST.txt")
# colnames(ENSG.ENST)[2]="txid"

  TXExpress=cbind(rownames(TXExpress),TXExpress)
  colnames(TXExpress)[1]="txid"
  colnames(TXExpress)=c("txid","TX_total_0","TX_total_4","TX_total_8","TX_preMRNA_0","TX_preMRNA_4","TX_preMRNA_8","TX_synthesis_0", 
                        "TX_synthesis_4","TX_synthesis_8","TX_degradation_0","TX_degradation_4","TX_degradation_8","TX_processing_0","TX_processing_4","TX_processing_8")
  TXExpress=merge(TXExpress,ENSG.ENST, by = "txid")
  colnames(TXExpress)[17]="geneid"

  GeneExpression=cbind(rownames(GeneExpression),GeneExpression)
  colnames(GeneExpression)[1]="geneid"

  TXExpress=merge(TXExpress,GeneExpression, by ="geneid")
  return(TXExpress)
}

Express=TXLabAREP1[,c(1,2,12,13,14,27,28,29)]
# 删除NA和0的数据
  DeleteSynthsisNAdata<-function(Express){
  Express=na.omit(Express)
  Express[apply(Express[,3:8]<=0, FUN = any, 1), ] = NA
  Express=na.omit(Express)
  print(length(Express[,1]))
  return(Express)
}



# calculate the variance of one ensg's ensts
varcal<-function(data)             #calculate the variance of one ensg's ensts
{                                  #the data must be a data.frame start with two culolumns of ensg 

  ensg.name<-unique(data[duplicated(data$geneid),1])      #enst names in character and several columns of expression value in numeric
  
  result<-c()
  for(i in ensg.name)
  {
   
    index.enst<-data[,1]==i    #the index of enst for one ensg

    result<-rbind(result,data[index.enst,])
    
  } 
 
  varexplain<-c()
  for(i in 1:length(result[,1])){
  lm = lm(unlist(result[i,3:5]) ~unlist(result[i,6:8]))
  summary(lm)$r.squared 
  tx=as.data.frame(result[i,1]) 
  gene=as.data.frame(result[i,2])
  varevalue=as.data.frame(cbind(tx,gene,summary(lm)$r.squared ))
  varexplain<- rbind(varexplain,varevalue )
  }
  return(varexplain)
}



# 编写的基因转录本间的变异
Express=dATA
varexplain<-c()

VarianceBetweenGene<-function(Express){
  for(tx in 1:(length(Express[,1])-1))
    for(atx in (tx+1):(length(Express[,1]))){
      lm = lm(unlist(Express[tx,3:4]) ~unlist(Express[atx,3:4]))
      summary(lm)$r.squared 
      tx1=as.data.frame(Express[tx,2]) 
      tx2=as.data.frame(Express[atx,2])
      varevalue=as.data.frame(cbind(tx1,tx2,summary(lm)$r.squared ))
      varexplain<- rbind(varexplain,varevalue )
    }
}




for(tx in 1:(length(Express[,1])-1))
  for(atx in (tx+1):(length(Express[,1]))){
    print(tx)
    print(atx)
  }


