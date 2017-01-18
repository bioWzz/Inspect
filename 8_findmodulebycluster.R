# 从聚类结果中得到转录模式
library(plyr)
# 读入聚类结果
label.A.1.high <- read.delim("D:/label-A-1-high.txt")
ClusterResult=label.A.1.high
# 读入lab结果
load("D:/ivan/ivan6-cell/Results/Result_getfromthePROGRAMME/LabAandCresults/labAandC.RData")
# load("D:/ivan/Harm/NewTestFour/my/LabAGene.RData")
# load("D:/ivan/Harm/NewTestFour/my/Results/NewResults/Rep2LabADegradation.Rdata")
# load("D:/ivan/Harm/NewTestFour/my/Results/NewResults/Rep2LabAProcessing.Rdata")
# load("D:/ivan/Harm/NewTestFour/my/Results/NewResults/Rep2LabASynthesis.Rdata")


#将聚类结果中的ENSG的id进行转换
ClusterResult=ConvertName(ClusterResult)

ClusterResult=ClusterResult[,c(1,9)]
colnames(ClusterResult)[1]="id"

# 将ENSG的id及数据放置在ENSG的后面
DATA=GeneLabAREP1
DATA=ConvertRowNametoFirstCol(DATA)
data=merge(DATA,ClusterResult, by ="id")



# 提取各类
LabAREP1c1=folodchange(extractClusterResult(data,"1"))
LabAREP1c2=folodchange(extractClusterResult(data,"2"))
LabAREP1c3=folodchange(extractClusterResult(data,"3"))
LabAREP1c4=folodchange(extractClusterResult(data,"4"))
LabAREP1c5=folodchange(extractClusterResult(data,"5"))
LabAREP1c6=folodchange(extractClusterResult(data,"6"))
LabAREP1c7=folodchange(extractClusterResult(data,"7"))
LabAREP1c8=folodchange(extractClusterResult(data,"8"))
LabAREP1c9=folodchange(extractClusterResult(data,"9"))

# 提取各类,不转log
LabAREP1c1=extractClusterResult(data,"1")
LabAREP1c2=extractClusterResult(data,"2")
LabAREP1c3=extractClusterResult(data,"3")
LabAREP1c4=extractClusterResult(data,"4")
LabAREP1c5=extractClusterResult(data,"5")
LabAREP1c6=extractClusterResult(data,"6")
LabAREP1c7=extractClusterResult(data,"7")
LabAREP1c8=extractClusterResult(data,"8")
LabAREP1c9=extractClusterResult(data,"9")

# 画箱式图
LabAREP1c1ep=as.data.frame(log2(unlist(LabAREP1c1[,2])))
LabAREP1c2ep=as.data.frame(log2(unlist(LabAREP1c2[,2]))) 
LabAREP1c3ep=as.data.frame(log2(unlist(LabAREP1c3[,2]))) 
LabAREP1c4ep=as.data.frame(log2(unlist(LabAREP1c4[,2])) )
LabAREP1c5ep=as.data.frame(log2(unlist(LabAREP1c5[,2]))) 
LabAREP1c6ep=as.data.frame(log2(unlist(LabAREP1c6[,2]))) 
LabAREP1c7ep=as.data.frame(log2(unlist(LabAREP1c7[,2]))) 
LabAREP1c8ep=as.data.frame(log2(unlist(LabAREP1c8[,2]))) 
LabAREP1c9ep=as.data.frame(log2(unlist(LabAREP1c9[,2])))

LabAREP1c1ep=as.data.frame(log2(unlist(LabAREP1c1[,3])))
LabAREP1c2ep=as.data.frame(log2(unlist(LabAREP1c2[,3]))) 
LabAREP1c3ep=as.data.frame(log2(unlist(LabAREP1c3[,3]))) 
LabAREP1c4ep=as.data.frame(log2(unlist(LabAREP1c4[,3])) )
LabAREP1c5ep=as.data.frame(log2(unlist(LabAREP1c5[,3]))) 
LabAREP1c6ep=as.data.frame(log2(unlist(LabAREP1c6[,3]))) 
LabAREP1c7ep=as.data.frame(log2(unlist(LabAREP1c7[,3]))) 
LabAREP1c8ep=as.data.frame(log2(unlist(LabAREP1c8[,3]))) 
LabAREP1c9ep=as.data.frame(log2(unlist(LabAREP1c9[,3]))) 

LabAREP1c1ep=as.data.frame(log2(unlist(LabAREP1c1[,4])))
LabAREP1c2ep=as.data.frame(log2(unlist(LabAREP1c2[,4]))) 
LabAREP1c3ep=as.data.frame(log2(unlist(LabAREP1c3[,4]))) 
LabAREP1c4ep=as.data.frame(log2(unlist(LabAREP1c4[,4])) )
LabAREP1c5ep=as.data.frame(log2(unlist(LabAREP1c5[,4]))) 
LabAREP1c6ep=as.data.frame(log2(unlist(LabAREP1c6[,4]))) 
LabAREP1c7ep=as.data.frame(log2(unlist(LabAREP1c7[,4]))) 
LabAREP1c8ep=as.data.frame(log2(unlist(LabAREP1c8[,4]))) 
LabAREP1c9ep=as.data.frame(log2(unlist(LabAREP1c9[,4])))

LabAREP1c1ep=as.data.frame(log2(unlist(LabAREP1c1[,2:4])))
LabAREP1c2ep=as.data.frame(log2(unlist(LabAREP1c2[,2:4]))) 
LabAREP1c3ep=as.data.frame(log2(unlist(LabAREP1c3[,2:4]))) 
LabAREP1c4ep=as.data.frame(log2(unlist(LabAREP1c4[,2:4])) )
LabAREP1c5ep=as.data.frame(log2(unlist(LabAREP1c5[,2:4]))) 
LabAREP1c6ep=as.data.frame(log2(unlist(LabAREP1c6[,2:4]))) 
LabAREP1c7ep=as.data.frame(log2(unlist(LabAREP1c7[,2:4]))) 
LabAREP1c8ep=as.data.frame(log2(unlist(LabAREP1c8[,2:4]))) 
LabAREP1c9ep=as.data.frame(log2(unlist(LabAREP1c9[,2:4])))

LabAREP1c1ep=as.data.frame(log2(unlist(LabAREP1c1[,8:10])))
LabAREP1c2ep=as.data.frame(log2(unlist(LabAREP1c2[,8:10]))) 
LabAREP1c3ep=as.data.frame(log2(unlist(LabAREP1c3[,8:10]))) 
LabAREP1c4ep=as.data.frame(log2(unlist(LabAREP1c4[,8:10])) )
LabAREP1c5ep=as.data.frame(log2(unlist(LabAREP1c5[,8:10]))) 
LabAREP1c6ep=as.data.frame(log2(unlist(LabAREP1c6[,8:10]))) 
LabAREP1c7ep=as.data.frame(log2(unlist(LabAREP1c7[,8:10]))) 
LabAREP1c8ep=as.data.frame(log2(unlist(LabAREP1c8[,8:10]))) 
LabAREP1c9ep=as.data.frame(log2(unlist(LabAREP1c9[,8:10])))

LabAREP1c1ep=as.data.frame(log2(unlist(LabAREP1c1[,11:13])))
LabAREP1c2ep=as.data.frame(log2(unlist(LabAREP1c2[,11:13]))) 
LabAREP1c3ep=as.data.frame(log2(unlist(LabAREP1c3[,11:13]))) 
LabAREP1c4ep=as.data.frame(log2(unlist(LabAREP1c4[,11:13])) )
LabAREP1c5ep=as.data.frame(log2(unlist(LabAREP1c5[,11:13]))) 
LabAREP1c6ep=as.data.frame(log2(unlist(LabAREP1c6[,11:13]))) 
LabAREP1c7ep=as.data.frame(log2(unlist(LabAREP1c7[,11:13]))) 
LabAREP1c8ep=as.data.frame(log2(unlist(LabAREP1c8[,11:13]))) 
LabAREP1c9ep=as.data.frame(log2(unlist(LabAREP1c9[,11:13])))

LabAREP1c1ep=as.data.frame(log2(unlist(LabAREP1c1[,14])))
LabAREP1c2ep=as.data.frame(log2(unlist(LabAREP1c2[,14]))) 
LabAREP1c3ep=as.data.frame(log2(unlist(LabAREP1c3[,14]))) 
LabAREP1c4ep=as.data.frame(log2(unlist(LabAREP1c4[,14])) )
LabAREP1c5ep=as.data.frame(log2(unlist(LabAREP1c5[,14]))) 
LabAREP1c6ep=as.data.frame(log2(unlist(LabAREP1c6[,14]))) 
LabAREP1c7ep=as.data.frame(log2(unlist(LabAREP1c7[,14]))) 
LabAREP1c8ep=as.data.frame(log2(unlist(LabAREP1c8[,14]))) 
LabAREP1c9ep=as.data.frame(log2(unlist(LabAREP1c9[,14])))


LabAREP1c1ep=as.data.frame(log2(unlist(LabAREP1c1[,15])))
LabAREP1c2ep=as.data.frame(log2(unlist(LabAREP1c2[,15]))) 
LabAREP1c3ep=as.data.frame(log2(unlist(LabAREP1c3[,15]))) 
LabAREP1c4ep=as.data.frame(log2(unlist(LabAREP1c4[,15])) )
LabAREP1c5ep=as.data.frame(log2(unlist(LabAREP1c5[,15]))) 
LabAREP1c6ep=as.data.frame(log2(unlist(LabAREP1c6[,15]))) 
LabAREP1c7ep=as.data.frame(log2(unlist(LabAREP1c7[,15]))) 
LabAREP1c8ep=as.data.frame(log2(unlist(LabAREP1c8[,15]))) 
LabAREP1c9ep=as.data.frame(log2(unlist(LabAREP1c9[,15])))

LabAREP1c1ep=as.data.frame(log2(unlist(LabAREP1c1[,16])))
LabAREP1c2ep=as.data.frame(log2(unlist(LabAREP1c2[,16]))) 
LabAREP1c3ep=as.data.frame(log2(unlist(LabAREP1c3[,16]))) 
LabAREP1c4ep=as.data.frame(log2(unlist(LabAREP1c4[,16])) )
LabAREP1c5ep=as.data.frame(log2(unlist(LabAREP1c5[,16]))) 
LabAREP1c6ep=as.data.frame(log2(unlist(LabAREP1c6[,16]))) 
LabAREP1c7ep=as.data.frame(log2(unlist(LabAREP1c7[,16]))) 
LabAREP1c8ep=as.data.frame(log2(unlist(LabAREP1c8[,16]))) 
LabAREP1c9ep=as.data.frame(log2(unlist(LabAREP1c9[,16])))

LabAREP1c1ep=as.data.frame(log2(unlist(LabAREP1c1[,14:16])))
LabAREP1c2ep=as.data.frame(log2(unlist(LabAREP1c2[,14:16]))) 
LabAREP1c3ep=as.data.frame(log2(unlist(LabAREP1c3[,14:16]))) 
LabAREP1c4ep=as.data.frame(log2(unlist(LabAREP1c4[,14:16])) )
LabAREP1c5ep=as.data.frame(log2(unlist(LabAREP1c5[,14:16]))) 
LabAREP1c6ep=as.data.frame(log2(unlist(LabAREP1c6[,14:16]))) 
LabAREP1c7ep=as.data.frame(log2(unlist(LabAREP1c7[,14:16]))) 
LabAREP1c8ep=as.data.frame(log2(unlist(LabAREP1c8[,14:16]))) 
LabAREP1c9ep=as.data.frame(log2(unlist(LabAREP1c9[,14:16])))

# cluster为单位
LabAREP1c1ep0=as.data.frame(log2(unlist(LabAREP1c1[,2])))
LabAREP1c1ep4=as.data.frame(log2(unlist(LabAREP1c1[,3])))
LabAREP1c1ep8=as.data.frame(log2(unlist(LabAREP1c1[,4])))

LabAREP1c1s0=as.data.frame(log2(unlist(LabAREP1c1[,8])))
LabAREP1c1s4=as.data.frame(log2(unlist(LabAREP1c1[,9])))
LabAREP1c1s8=as.data.frame(log2(unlist(LabAREP1c1[,10])))

LabAREP1c1d0=as.data.frame(log2(unlist(LabAREP1c1[,11])))
LabAREP1c1d4=as.data.frame(log2(unlist(LabAREP1c1[,12])))
LabAREP1c1d8=as.data.frame(log2(unlist(LabAREP1c1[,13])))

LabAREP1c1p0=as.data.frame(log2(unlist(LabAREP1c1[,14])))
LabAREP1c1p4=as.data.frame(log2(unlist(LabAREP1c1[,15])))
LabAREP1c1p8=as.data.frame(log2(unlist(LabAREP1c1[,16])))

LabAREP1c2ep0=as.data.frame(log2(unlist(LabAREP1c2[,2])))
LabAREP1c2ep4=as.data.frame(log2(unlist(LabAREP1c2[,3])))
LabAREP1c2ep8=as.data.frame(log2(unlist(LabAREP1c2[,4])))

LabAREP1c2s0=as.data.frame(log2(unlist(LabAREP1c2[,8])))
LabAREP1c2s4=as.data.frame(log2(unlist(LabAREP1c2[,9])))
LabAREP1c2s8=as.data.frame(log2(unlist(LabAREP1c2[,10])))

LabAREP1c2d0=as.data.frame(log2(unlist(LabAREP1c2[,11])))
LabAREP1c2d4=as.data.frame(log2(unlist(LabAREP1c2[,12])))
LabAREP1c2d8=as.data.frame(log2(unlist(LabAREP1c2[,13])))

LabAREP1c2p0=as.data.frame(log2(unlist(LabAREP1c2[,14])))
LabAREP1c2p4=as.data.frame(log2(unlist(LabAREP1c2[,15])))
LabAREP1c2p8=as.data.frame(log2(unlist(LabAREP1c2[,16])))

LabAREP1c3ep0=as.data.frame(log2(unlist(LabAREP1c3[,2])))
LabAREP1c3ep4=as.data.frame(log2(unlist(LabAREP1c3[,3])))
LabAREP1c3ep8=as.data.frame(log2(unlist(LabAREP1c3[,4])))

LabAREP1c3s0=as.data.frame(log2(unlist(LabAREP1c3[,8])))
LabAREP1c3s4=as.data.frame(log2(unlist(LabAREP1c3[,9])))
LabAREP1c3s8=as.data.frame(log2(unlist(LabAREP1c3[,10])))

LabAREP1c3d0=as.data.frame(log2(unlist(LabAREP1c3[,11])))
LabAREP1c3d4=as.data.frame(log2(unlist(LabAREP1c3[,12])))
LabAREP1c3d8=as.data.frame(log2(unlist(LabAREP1c3[,13])))

LabAREP1c3p0=as.data.frame(log2(unlist(LabAREP1c3[,14])))
LabAREP1c3p4=as.data.frame(log2(unlist(LabAREP1c3[,15])))
LabAREP1c3p8=as.data.frame(log2(unlist(LabAREP1c3[,16])))


LabAREP1c4ep0=as.data.frame(log2(unlist(LabAREP1c4[,2])))
LabAREP1c4ep4=as.data.frame(log2(unlist(LabAREP1c4[,3])))
LabAREP1c4ep8=as.data.frame(log2(unlist(LabAREP1c4[,4])))

LabAREP1c4s0=as.data.frame(log2(unlist(LabAREP1c4[,8])))
LabAREP1c4s4=as.data.frame(log2(unlist(LabAREP1c4[,9])))
LabAREP1c4s8=as.data.frame(log2(unlist(LabAREP1c4[,10])))

LabAREP1c4d0=as.data.frame(log2(unlist(LabAREP1c4[,11])))
LabAREP1c4d4=as.data.frame(log2(unlist(LabAREP1c4[,12])))
LabAREP1c4d8=as.data.frame(log2(unlist(LabAREP1c4[,13])))

LabAREP1c4p0=as.data.frame(log2(unlist(LabAREP1c4[,14])))
LabAREP1c4p4=as.data.frame(log2(unlist(LabAREP1c4[,15])))
LabAREP1c4p8=as.data.frame(log2(unlist(LabAREP1c4[,16])))

LabAREP1c5ep0=as.data.frame(log2(unlist(LabAREP1c5[,2])))
LabAREP1c5ep4=as.data.frame(log2(unlist(LabAREP1c5[,3])))
LabAREP1c5ep8=as.data.frame(log2(unlist(LabAREP1c5[,4])))

LabAREP1c5s0=as.data.frame(log2(unlist(LabAREP1c5[,8])))
LabAREP1c5s4=as.data.frame(log2(unlist(LabAREP1c5[,9])))
LabAREP1c5s8=as.data.frame(log2(unlist(LabAREP1c5[,10])))

LabAREP1c5d0=as.data.frame(log2(unlist(LabAREP1c5[,11])))
LabAREP1c5d4=as.data.frame(log2(unlist(LabAREP1c5[,12])))
LabAREP1c5d8=as.data.frame(log2(unlist(LabAREP1c5[,13])))

LabAREP1c5p0=as.data.frame(log2(unlist(LabAREP1c5[,14])))
LabAREP1c5p4=as.data.frame(log2(unlist(LabAREP1c5[,15])))
LabAREP1c5p8=as.data.frame(log2(unlist(LabAREP1c5[,16])))

LabAREP1c6ep0=as.data.frame(log2(unlist(LabAREP1c6[,2])))
LabAREP1c6ep4=as.data.frame(log2(unlist(LabAREP1c6[,3])))
LabAREP1c6ep8=as.data.frame(log2(unlist(LabAREP1c6[,4])))

LabAREP1c6s0=as.data.frame(log2(unlist(LabAREP1c6[,8])))
LabAREP1c6s4=as.data.frame(log2(unlist(LabAREP1c6[,9])))
LabAREP1c6s8=as.data.frame(log2(unlist(LabAREP1c6[,10])))

LabAREP1c6d0=as.data.frame(log2(unlist(LabAREP1c6[,11])))
LabAREP1c6d4=as.data.frame(log2(unlist(LabAREP1c6[,12])))
LabAREP1c6d8=as.data.frame(log2(unlist(LabAREP1c6[,13])))

LabAREP1c6p0=as.data.frame(log2(unlist(LabAREP1c6[,14])))
LabAREP1c6p4=as.data.frame(log2(unlist(LabAREP1c6[,15])))
LabAREP1c6p8=as.data.frame(log2(unlist(LabAREP1c6[,16])))

LabAREP1c7ep0=as.data.frame(log2(unlist(LabAREP1c7[,2])))
LabAREP1c7ep4=as.data.frame(log2(unlist(LabAREP1c7[,3])))
LabAREP1c7ep8=as.data.frame(log2(unlist(LabAREP1c7[,4])))

LabAREP1c7s0=as.data.frame(log2(unlist(LabAREP1c7[,8])))
LabAREP1c7s4=as.data.frame(log2(unlist(LabAREP1c7[,9])))
LabAREP1c7s8=as.data.frame(log2(unlist(LabAREP1c7[,10])))

LabAREP1c7d0=as.data.frame(log2(unlist(LabAREP1c7[,11])))
LabAREP1c7d4=as.data.frame(log2(unlist(LabAREP1c7[,12])))
LabAREP1c7d8=as.data.frame(log2(unlist(LabAREP1c7[,13])))

LabAREP1c7p0=as.data.frame(log2(unlist(LabAREP1c7[,14])))
LabAREP1c7p4=as.data.frame(log2(unlist(LabAREP1c7[,15])))
LabAREP1c7p8=as.data.frame(log2(unlist(LabAREP1c7[,16])))

LabAREP1c8ep0=as.data.frame(log2(unlist(LabAREP1c8[,2])))
LabAREP1c8ep4=as.data.frame(log2(unlist(LabAREP1c8[,3])))
LabAREP1c8ep8=as.data.frame(log2(unlist(LabAREP1c8[,4])))

LabAREP1c8s0=as.data.frame(log2(unlist(LabAREP1c8[,8])))
LabAREP1c8s4=as.data.frame(log2(unlist(LabAREP1c8[,9])))
LabAREP1c8s8=as.data.frame(log2(unlist(LabAREP1c8[,10])))

LabAREP1c8d0=as.data.frame(log2(unlist(LabAREP1c8[,11])))
LabAREP1c8d4=as.data.frame(log2(unlist(LabAREP1c8[,12])))
LabAREP1c8d8=as.data.frame(log2(unlist(LabAREP1c8[,13])))

LabAREP1c8p0=as.data.frame(log2(unlist(LabAREP1c8[,14])))
LabAREP1c8p4=as.data.frame(log2(unlist(LabAREP1c8[,15])))
LabAREP1c8p8=as.data.frame(log2(unlist(LabAREP1c8[,16])))

LabAREP1c9ep0=as.data.frame(log2(unlist(LabAREP1c9[,2])))
LabAREP1c9ep4=as.data.frame(log2(unlist(LabAREP1c9[,3])))
LabAREP1c9ep8=as.data.frame(log2(unlist(LabAREP1c9[,4])))

LabAREP1c9s0=as.data.frame(log2(unlist(LabAREP1c9[,8])))
LabAREP1c9s4=as.data.frame(log2(unlist(LabAREP1c9[,9])))
LabAREP1c9s8=as.data.frame(log2(unlist(LabAREP1c9[,10])))

LabAREP1c9d0=as.data.frame(log2(unlist(LabAREP1c9[,11])))
LabAREP1c9d4=as.data.frame(log2(unlist(LabAREP1c9[,12])))
LabAREP1c9d8=as.data.frame(log2(unlist(LabAREP1c9[,13])))

LabAREP1c9p0=as.data.frame(log2(unlist(LabAREP1c9[,14])))
LabAREP1c9p4=as.data.frame(log2(unlist(LabAREP1c9[,15])))
LabAREP1c9p8=as.data.frame(log2(unlist(LabAREP1c9[,16])))

u=rbind.fill(LabAREP1c1ep0,LabAREP1c1ep4,LabAREP1c1ep8
             ,LabAREP1c1s0,LabAREP1c1s4,LabAREP1c1s8
             ,LabAREP1c1d0,LabAREP1c1d4,LabAREP1c1d8
             ,LabAREP1c1p0,LabAREP1c1p4,LabAREP1c1p8)

u=rbind.fill(LabAREP1c2ep0,LabAREP1c2ep4,LabAREP1c2ep8
             ,LabAREP1c2s0,LabAREP1c2s4,LabAREP1c2s8
             ,LabAREP1c2d0,LabAREP1c2d4,LabAREP1c2d8
             ,LabAREP1c2p0,LabAREP1c2p4,LabAREP1c2p8)

u=rbind.fill(LabAREP1c3ep0,LabAREP1c3ep4,LabAREP1c3ep8
             ,LabAREP1c3s0,LabAREP1c3s4,LabAREP1c3s8
             ,LabAREP1c3d0,LabAREP1c3d4,LabAREP1c3d8
             ,LabAREP1c3p0,LabAREP1c3p4,LabAREP1c3p8)

u=rbind.fill(LabAREP1c4ep0,LabAREP1c4ep4,LabAREP1c4ep8
             ,LabAREP1c4s0,LabAREP1c4s4,LabAREP1c4s8
             ,LabAREP1c4d0,LabAREP1c4d4,LabAREP1c4d8
             ,LabAREP1c4p0,LabAREP1c4p4,LabAREP1c4p8)

u=rbind.fill(LabAREP1c5ep0,LabAREP1c5ep4,LabAREP1c5ep8
             ,LabAREP1c5s0,LabAREP1c5s4,LabAREP1c5s8
             ,LabAREP1c5d0,LabAREP1c5d4,LabAREP1c5d8
             ,LabAREP1c5p0,LabAREP1c5p4,LabAREP1c5p8)

u=rbind.fill(LabAREP1c6ep0,LabAREP1c6ep4,LabAREP1c6ep8
             ,LabAREP1c6s0,LabAREP1c6s4,LabAREP1c6s8
             ,LabAREP1c6d0,LabAREP1c6d4,LabAREP1c6d8
             ,LabAREP1c6p0,LabAREP1c6p4,LabAREP1c6p8)

u=rbind.fill(LabAREP1c7ep0,LabAREP1c7ep4,LabAREP1c7ep8
             ,LabAREP1c7s0,LabAREP1c7s4,LabAREP1c7s8
             ,LabAREP1c7d0,LabAREP1c7d4,LabAREP1c7d8
             ,LabAREP1c7p0,LabAREP1c7p4,LabAREP1c7p8)

u=rbind.fill(LabAREP1c8ep0,LabAREP1c8ep4,LabAREP1c8ep8
             ,LabAREP1c8s0,LabAREP1c8s4,LabAREP1c8s8
             ,LabAREP1c8d0,LabAREP1c8d4,LabAREP1c8d8
             ,LabAREP1c8p0,LabAREP1c8p4,LabAREP1c8p8)

u=rbind.fill(LabAREP1c9ep0,LabAREP1c9ep4,LabAREP1c9ep8
             ,LabAREP1c9s0,LabAREP1c9s4,LabAREP1c9s8
             ,LabAREP1c9d0,LabAREP1c9d4,LabAREP1c9d8
             ,LabAREP1c9p0,LabAREP1c9p4,LabAREP1c9p8)

boxplot(u,col=c("steelblue"),ylim = c(-10, 10),ylab="Length of AA",las=1, font.lab=2)


# folod expression

folodchange<-function(Express){
  Express[,3]=log2(Express[,3]/Express[,2])
  Express[,4]=log2(Express[,4]/Express[,2])
  Express[,2]=0
  
  Express[,6]=log2(Express[,6]/Express[,5])
  Express[,7]=log2(Express[,7]/Express[,5])
  Express[,5]=0
  
  Express[,9]=log2(Express[,9]/Express[,8])
  Express[,10]=log2(Express[,10]/Express[,8])
  Express[,8]=0
 
  Express[,12]=log2(Express[,12]/Express[,11])
  Express[,13]=log2(Express[,13]/Express[,11])
  Express[,11]=0

  Express[,15]=log2(Express[,15]/Express[,14])
  Express[,16]=log2(Express[,16]/Express[,14])
  Express[,14]=0
  
  return(Express)
}
  
  
  

# 提取子类，i表示类的名称
extractClusterResult<-function(Express,i){
  index.enst<-Express[,17]==i    #the index of enst for one ensg
  z=Express[index.enst,]
  return(z)
}


















#convert rowname to the first colnames
ConvertRowNametoFirstCol<-function(Express){

  name=as.data.frame(rownames(Express))
  Express=cbind(name,Express)
  colnames(Express)[1]="id"
  return(Express)
}

#convert ENSG id to ENSG id eg. convert ensg12.1 to ensg12  
Express=ClusterResult
ConvertName<-function(Express){
  rownames(Express)=Express[,1]
  Express=Express[,2:length(Express[1,])]
  name=as.data.frame(rownames(Express))
  name1=as.data.frame(strsplit(as.character(name[,1]),'\\.'))
  name1<- t(name1)
  name1=as.data.frame(name1[,1])
  colnames(name1)[1]="id"
  Express=cbind(name1,Express)
  rownames(Express)=Express[,1]
  Express=Express[,1:length(Express[1,])]
  return(Express)
}
