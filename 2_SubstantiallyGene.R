rm(list=ls())
load("D:/ivan/ivan6-cell/Results/结果数据/LabAGene.RData")
# load("D:/ivan/ivan6-cell/Results/结果数据/LabATX.RData")


# 
# Delete the no use varies 
vari=objects()
delva=as.vector(grep("^La",vari,invert=TRUE))  
vari=vari[delva]
rm(list=vari)

# Gene
# 
#对于substantially gene 的提取，保留在所有时刻的内含子，外显子，及计算出来的率的表达都是大于0的
# LabA
# rep1
LabAREP1Express1=LabAREP1Express
LabAREP1Express1[apply(LabAREP1Express1<=0, FUN = any, 1), ] = NA
LabAREP1Express1=na.omit(LabAREP1Express1)

# 对于一个基因计算所有时刻所有表达值的和
allExpression=as.double(LabAREP1Express1[,1])+as.double(LabAREP1Express1[,2])+as.double(LabAREP1Express1[,3])
              +as.double(LabAREP1Express1[,4])+as.double(LabAREP1Express1[,5])+as.double(LabAREP1Express1[,6])

LabAREP1Express1=cbind(as.data.frame(LabAREP1Express1),as.data.frame(allExpression))
LabAREP1Express1=LabAREP1Express1[order(LabAREP1Express1[,16],decreasing=T),]

S10Gene_LabArep1=LabAREP1Express1[1:as.integer(length(LabAREP1Express1[,1])*0.1),]
S70Gene_LabArep1=LabAREP1Express1[1:as.integer(length(LabAREP1Express1[,1])*0.7),]
S17Gene_LabArep1=LabAREP1Express1[1:as.integer(length(LabAREP1Express1[,1])*0.17),]

# S10TX_LabArep1=LabAREP1Express1[1:as.integer(length(LabAREP1Express1[,1])*0.1),]
# S70TX_LabArep1=LabAREP1Express1[1:as.integer(length(LabAREP1Express1[,1])*0.7),]
# S17TX_LabArep1=LabAREP1Express1[1:as.integer(length(LabAREP1Express1[,1])*0.17),]

# rep2
LabAREP2Express1=LabAREP2Express
LabAREP2Express1[apply(LabAREP2Express1<=0, FUN = any, 1), ] = NA
LabAREP2Express1=na.omit(LabAREP2Express1)

# 对于一个基因计算所有时刻所有表达值的和
allExpression=as.double(LabAREP2Express1[,1])+as.double(LabAREP2Express1[,2])+as.double(LabAREP2Express1[,3])
              +as.double(LabAREP2Express1[,4])+as.double(LabAREP2Express1[,5])+as.double(LabAREP2Express1[,6])

LabAREP2Express1=cbind(as.data.frame(LabAREP2Express1),as.data.frame(allExpression))
LabAREP2Express1=LabAREP2Express1[order(LabAREP2Express1[,16],decreasing=T),]

S10Gene_LabArep2=LabAREP2Express1[1:as.integer(length(LabAREP2Express1[,1])*0.1),]
S70Gene_LabArep2=LabAREP2Express1[1:as.integer(length(LabAREP2Express1[,1])*0.7),]
S17Gene_LabArep2=LabAREP2Express1[1:as.integer(length(LabAREP2Express1[,1])*0.17),]

# S10TX_LabArep2=LabAREP2Express1[1:as.integer(length(LabAREP2Express1[,1])*0.1),]
# S70TX_LabArep2=LabAREP2Express1[1:as.integer(length(LabAREP2Express1[,1])*0.7),]
# S17TX_LabArep2=LabAREP2Express1[1:as.integer(length(LabAREP2Express1[,1])*0.17),]

# rep1 and rep2
LabAREP12Express1=LabAREP12Express
LabAREP12Express1[apply(LabAREP12Express1<=0, FUN = any, 1), ] = NA
LabAREP12Express1=na.omit(LabAREP12Express1)

# 对于一个基因计算所有时刻所有表达值的和
allExpression=as.double(LabAREP12Express1[,1])+as.double(LabAREP12Express1[,2])+as.double(LabAREP12Express1[,3])
              +as.double(LabAREP12Express1[,4])+as.double(LabAREP12Express1[,5])+as.double(LabAREP12Express1[,6])

LabAREP12Express1=cbind(as.data.frame(LabAREP12Express1),as.data.frame(allExpression))
LabAREP12Express1=LabAREP12Express1[order(LabAREP12Express1[,16],decreasing=T),]

S10Gene_LabArep12=LabAREP12Express1[1:as.integer(length(LabAREP12Express1[,1])*0.1),]
S70Gene_LabArep12=LabAREP12Express1[1:as.integer(length(LabAREP12Express1[,1])*0.7),]
S17Gene_LabArep12=LabAREP12Express1[1:as.integer(length(LabAREP12Express1[,1])*0.17),]

# S10TX_LabArep12=LabAREP12Express1[1:as.integer(length(LabAREP12Express1[,1])*0.1),]
# S70TX_LabArep12=LabAREP12Express1[1:as.integer(length(LabAREP12Express1[,1])*0.7),]
# S17TX_LabArep12=LabAREP12Express1[1:as.integer(length(LabAREP12Express1[,1])*0.17),]

# Delete the no use varies 
vari=objects()
delva=as.vector(grep("^La",vari))  
vari=vari[delva]
rm(list=vari)

# 
# 保存变量
# 

# save.image(file="D:/ivan/ivan6-cell/Results/结果数据/LabATXSubstantiallyGene.RData")
save.image(file="D:/ivan/ivan6-cell/Results/结果数据/LabAGeneSubstantiallyGene.RData")



rm(list=ls())

# load("D:/ivan/ivan6-cell/Results/结果数据/LabCGene.RData")
 load("D:/ivan/ivan6-cell/Results/结果数据/LabCTX.RData")

# 
# Delete the no use varies 
vari=objects()
delva=as.vector(grep("^La",vari,invert=TRUE))  
vari=vari[delva]
rm(list=vari)

# Gene
# 
#对于substantially gene 的提取，保留在所有时刻的内含子，外显子，及计算出来的率的表达都是大于0的
# LabC
# rep1
LabAREP1Express1=LabCREP1Express
LabAREP1Express1[apply(LabAREP1Express1<=0, FUN = any, 1), ] = NA
LabAREP1Express1=na.omit(LabAREP1Express1)

# 对于一个基因计算所有时刻所有表达值的和
allExpression=as.double(LabAREP1Express1[,1])+as.double(LabAREP1Express1[,2])+as.double(LabAREP1Express1[,3])
              +as.double(LabAREP1Express1[,4])+as.double(LabAREP1Express1[,5])+as.double(LabAREP1Express1[,6])

LabAREP1Express1=cbind(as.data.frame(LabAREP1Express1),as.data.frame(allExpression))
LabAREP1Express1=LabAREP1Express1[order(LabAREP1Express1[,16],decreasing=T),]

# S10Gene_LabCrep1=LabAREP1Express1[1:as.integer(length(LabAREP1Express1[,1])*0.1),]
# S70Gene_LabCrep1=LabAREP1Express1[1:as.integer(length(LabAREP1Express1[,1])*0.7),]
# S17Gene_LabCrep1=LabAREP1Express1[1:as.integer(length(LabAREP1Express1[,1])*0.17),]

S10TX_LabCrep1=LabAREP1Express1[1:as.integer(length(LabAREP1Express1[,1])*0.1),]
S70TX_LabCrep1=LabAREP1Express1[1:as.integer(length(LabAREP1Express1[,1])*0.7),]
S17TX_LabCrep1=LabAREP1Express1[1:as.integer(length(LabAREP1Express1[,1])*0.17),]

# rep2
LabAREP2Express1=LabCREP2Express
LabAREP2Express1[apply(LabAREP2Express1<=0, FUN = any, 1), ] = NA
LabAREP2Express1=na.omit(LabAREP2Express1)

# 对于一个基因计算所有时刻所有表达值的和
allExpression=as.double(LabAREP2Express1[,1])+as.double(LabAREP2Express1[,2])+as.double(LabAREP2Express1[,3])
              +as.double(LabAREP2Express1[,4])+as.double(LabAREP2Express1[,5])+as.double(LabAREP2Express1[,6])

LabAREP2Express1=cbind(as.data.frame(LabAREP2Express1),as.data.frame(allExpression))
LabAREP2Express1=LabAREP2Express1[order(LabAREP2Express1[,16],decreasing=T),]

# S10Gene_LabCrep2=LabAREP2Express1[1:as.integer(length(LabAREP2Express1[,1])*0.1),]
# S70Gene_LabCrep2=LabAREP2Express1[1:as.integer(length(LabAREP2Express1[,1])*0.7),]
# S17Gene_LabCrep2=LabAREP2Express1[1:as.integer(length(LabAREP2Express1[,1])*0.17),]

S10TX_LabCrep2=LabAREP2Express1[1:as.integer(length(LabAREP2Express1[,1])*0.1),]
S70TX_LabCrep2=LabAREP2Express1[1:as.integer(length(LabAREP2Express1[,1])*0.7),]
S17TX_LabCrep2=LabAREP2Express1[1:as.integer(length(LabAREP2Express1[,1])*0.17),]

# rep1 and rep2
LabAREP12Express1=LabCREP12Express
LabAREP12Express1[apply(LabAREP12Express1<=0, FUN = any, 1), ] = NA
LabAREP12Express1=na.omit(LabAREP12Express1)

# 对于一个基因计算所有时刻所有表达值的和
allExpression=as.double(LabAREP12Express1[,1])+as.double(LabAREP12Express1[,2])+as.double(LabAREP12Express1[,3])
              +as.double(LabAREP12Express1[,4])+as.double(LabAREP12Express1[,5])+as.double(LabAREP12Express1[,6])

LabAREP12Express1=cbind(as.data.frame(LabAREP12Express1),as.data.frame(allExpression))
LabAREP12Express1=LabAREP12Express1[order(LabAREP12Express1[,16],decreasing=T),]
# S10Gene_LabCrep12=LabAREP12Express1[1:as.integer(length(LabAREP12Express1[,1])*0.1),]
# S70Gene_LabCrep12=LabAREP12Express1[1:as.integer(length(LabAREP12Express1[,1])*0.7),]
# S17Gene_LabCrep12=LabAREP12Express1[1:as.integer(length(LabAREP12Express1[,1])*0.17),]

S10TX_LabCrep12=LabAREP12Express1[1:as.integer(length(LabAREP12Express1[,1])*0.1),]
S70TX_LabCrep12=LabAREP12Express1[1:as.integer(length(LabAREP12Express1[,1])*0.7),]
S17TX_LabCrep12=LabAREP12Express1[1:as.integer(length(LabAREP12Express1[,1])*0.17),]

# Delete the no use varies 
vari=objects()
delva=as.vector(grep("^La",vari))  
vari=vari[delva]
rm(list=vari)

# 
# 保存变量
# 


# save.image(file="D:/ivan/ivan6-cell/Results/结果数据/LabCGeneSubstantiallyGene.RData")
save.image(file="D:/ivan/ivan6-cell/Results/结果数据/LabCTXSubstantiallyGene.RData")
