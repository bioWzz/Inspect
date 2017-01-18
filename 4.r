rm(list=ls())
library(GenomicFeatures)
library(INSPEcT)
library(Rsamtools)
library(deSolve)
library(BiocParallel)
library(compiler)
library(sqldf)
library(edgeR)
library(splines)
###数据预处理：
# 利用edgeR 软件分别计算出EXON和intron的表达，对数据进行预处理
# (1)
# 
# 
# 
txdb<- makeTxDbFromGFF(file="D:/ivan/Harm/NewTestOne/NewTest/Newgencode.v23.PR8.gtf",format="gtf",taxonomyId=9606)
su4<-c("D:/ProgramFilesRRR/R-3.3.1/library/INSPEcT/extdata/MyRep1/LabA0hpi.bam",
             "D:/ProgramFilesRRR/R-3.3.1/library/INSPEcT/extdata/MyRep1/LabA4hpi.bam",
             "D:/ProgramFilesRRR/R-3.3.1/library/INSPEcT/extdata/MyRep1/LabA8hpi.bam")

total<-c("D:/ProgramFilesRRR/R-3.3.1/library/INSPEcT/extdata/MyRep1/UnlA0hpi.bam",
               "D:/ProgramFilesRRR/R-3.3.1/library/INSPEcT/extdata/MyRep1/UnlA4hpi.bam",
               "D:/ProgramFilesRRR/R-3.3.1/library/INSPEcT/extdata/MyRep1/UnlA8hpi.bam")

REP1RPKM <- makeRPKMs(txdb,su4,total)
# RPKM
REP14SU_Exon<-REP1RPKM$rpkms$foursu_exons
REP14SU_Intron<-REP1RPKM$rpkms$foursu_introns
REP1Total_Exon<-REP1RPKM$rpkms$total_exons
REP1Total_Intron<-REP1RPKM$rpkms$total_introns

# count
CountREP14SU_Exon<-REP1RPKM$counts$foursu$exonCounts
CountREP14SU_Intron<-REP1RPKM$counts$foursu$intronCounts
CountREP1Total_Exon<-REP1RPKM$counts$total$exonCounts
CountREP1Total_Intron<-REP1RPKM$counts$total$intronCounts

# CPM=C/N*1000000

#caulate the EXON_RPKM value with the package of edge
# write.table(Exonwidth, "D:/Exonwidth.txt", sep="\t", row.names=T, col.names=T)
# write.table(Intronwidth, "D:/Intronwidth.txt", sep="\t", row.names=T, col.names=T)

REP1RPKM_Exon=as.data.frame(cbind(rownames(CountREP14SU_Exon),CountREP14SU_Exon,CountREP1Total_Exon))
colnames(REP1RPKM_Exon)<-c("GeneSymbol","REP14SU_ExonHour0","REP14SU_ExonHour4","REP14SU_ExonHour8","REP1Total_ExonHour0","REP1Total_ExonHour4","REP1Total_ExonHour8")

Exonwidth=as.data.frame(Exonwidth);Exonwidths=cbind(rownames(Exonwidth),Exonwidth);colnames(Exonwidths)<-c("GeneSymbol","Length")

ALL<- sqldf("select * from REP1RPKM_Exon join Exonwidths where REP1RPKM_Exon.GeneSymbol= Exonwidths.GeneSymbol")

rawdataSU_Exon=ALL[,1:7]
write.table(rawdataSU_Exon, file = "d:/REP1RPKM_Exon.txt", append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
GeneLength=ALL[,8:9]
write.table(GeneLength, file = "d:/FgeneLenth.txt", append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

rawdataSU_Exon=read.delim("d:/REP1RPKM_Exon.txt", check.names=FALSE, stringsAsFactors=FALSE)
GeneLength=read.delim("d:/FgeneLenth.txt", check.names=FALSE, stringsAsFactors=FALSE)
GeneLength=GeneLength[,2]


group <- c("pre0","pre4","pre8","total0","total4","total8")
d<- DGEList(counts=rawdataSU_Exon[,2:7], genes=data.frame(rawdataSU_Exon[,1],Length=GeneLength),group=group)
d <- calcNormFactors(d)
RPKM1 <- rpkm(d)
REP1EXONRPKM<-cbind(rawdataSU_Exon[,1],RPKM1)
write.table(REP1EXONRPKM, "D:/REP1EXONRPKM.txt", sep="\t", row.names=F, col.names=T)



#caulate the Intron_RPKM value with the package of edge
REP1RPKM_Intron=as.data.frame(cbind(rownames(CountREP14SU_Intron),CountREP14SU_Intron,CountREP1Total_Intron))
colnames(REP1RPKM_Intron)<-c("GeneSymbol","REP14SU_IntronHour0","REP14SU_IntronHour4","REP14SU_IntronHour8","REP1Total_IntronHour0","REP1Total_IntronHour4","REP1Total_IntronHour8")

Intronwidth=as.data.frame(Intronwidth)
Intronwidths=cbind(rownames(Intronwidth),Intronwidth) 
colnames(Intronwidths)<-c("GeneSymbol","Length")

IntronALL<- sqldf("select * from REP1RPKM_Intron join Intronwidths where REP1RPKM_Intron.GeneSymbol= Intronwidths.GeneSymbol")
rawdataSU_Intron=IntronALL[,1:7]
write.table(rawdataSU_Intron, file = "d:/REP1RPKM_Intron.txt", append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

GeneLength=IntronALL[,8:9]
write.table(GeneLength, file = "d:/IntronFgeneLenth.txt", append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

rawdataSU_Intron=read.delim("d:/REP1RPKM_Intron.txt", check.names=FALSE, stringsAsFactors=FALSE)
GeneLength=read.delim("d:/IntronFgeneLenth.txt", check.names=FALSE, stringsAsFactors=FALSE)
GeneLength=GeneLength[,2]
group <- c("pre0","pre4","pre8","total0","total4","total8")
d<- DGEList(counts=rawdataSU_Intron[,2:7], genes=data.frame(rawdataSU_Intron[,1],Length=GeneLength),group=group)
d <- calcNormFactors(d)
REP1RPKM1 <- rpkm(d)
REP1INTRONRPKM<-cbind(rawdataSU_Intron[,1],REP1RPKM1)
write.table(REP1INTRONRPKM, "D:/REP1IntronRPKM.txt", sep="\t", row.names=F, col.names=T)











paths_4su<-c("D:/ProgramFilesRRR/R-3.3.1/library/INSPEcT/extdata/MyRep2/LabA0.bam"
             ,"D:/ProgramFilesRRR/R-3.3.1/library/INSPEcT/extdata/MyRep2/LabA4.bam",
             "D:/ProgramFilesRRR/R-3.3.1/library/INSPEcT/extdata/MyRep2/LabA8.bam")

paths_total<-c("D:/ProgramFilesRRR/R-3.3.1/library/INSPEcT/extdata/MyRep2/ULabA0.bam"
,"D:/ProgramFilesRRR/R-3.3.1/library/INSPEcT/extdata/MyRep2/ULabA4.bam",
               "D:/ProgramFilesRRR/R-3.3.1/library/INSPEcT/extdata/MyRep2/ULabA8.bam")
REP2RPKM <- makeRPKMs(txdb, paths_4su, paths_total)
# RPKM
REP24SU_Exon<-REP2RPKM$rpkms$foursu_exons
REP24SU_Intron<-REP2RPKM$rpkms$foursu_introns
REP2Total_Exon<-REP2RPKM$rpkms$total_exons
REP2Total_Intron<-REP2RPKM$rpkms$total_introns

# count
CountREP24SU_Exon<-REP2RPKM$counts$foursu$exonCounts
CountREP24SU_Intron<-REP2RPKM$counts$foursu$intronCounts
CountREP2Total_Exon<-REP2RPKM$counts$total$exonCounts
CountREP2Total_Intron<-REP2RPKM$counts$total$intronCounts

#caulate the RPKM value with the package of edge

exonsDB <- reduce(exonsBy(txdb ,'gene'))
exonsDB <- exonsDB[elementNROWS(range(exonsDB))==1]
intronsDB <- psetdiff(unlist(range(exonsDB)),exonsDB)
intronsDB <- intronsDB[elementNROWS(intronsDB)>0]
Exonwidth=sapply(width(exonsDB),sum)
Intronwidth=sapply(width(intronsDB),sum)

REP2_Exon=as.data.frame(cbind(rownames(CountREP24SU_Exon),CountREP24SU_Exon,CountREP2Total_Exon))
colnames(REP2_Exon)<-c("GeneSymbol","REP24SU_ExonHour0","REP24SU_ExonHour4","REP24SU_ExonHour8","REP2Total_ExonHour0","REP2Total_ExonHour4","REP2Total_ExonHour8")

Exonwidth=as.data.frame(Exonwidth)
Exonwidths=cbind(rownames(Exonwidth),Exonwidth)
colnames(Exonwidths)<-c("GeneSymbol","Length")

REP2ALL<- sqldf("select * from REP2_Exon join Exonwidths where REP2_Exon.GeneSymbol= Exonwidths.GeneSymbol")
rawdata=REP2ALL[,1:7]
write.table(rawdata, file = "d:/REP2_Exon.txt", append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

GeneLength=REP2ALL[,8:9]
write.table(GeneLength, file = "d:/FgeneLenth.txt", append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

rawdata=read.delim("d:/REP2_Exon.txt", check.names=FALSE, stringsAsFactors=FALSE)
GeneLength=read.delim("d:/FgeneLenth.txt", check.names=FALSE, stringsAsFactors=FALSE)
GeneLength=GeneLength[,2]
group <- c("pre0","pre4","pre8","total0","total4","total8")
d<- DGEList(counts=rawdata[,2:7], genes=data.frame(rawdata[,1],Length=GeneLength),group=group)
d <- calcNormFactors(d)
RPKM1 <- rpkm(d)
REP2EXONRPKM<-cbind(rawdata[,1],RPKM1)
write.table(REP2EXONRPKM, "D:/REP2EXONRPKM.txt", sep="\t", row.names=F, col.names=T)


#caulate the Intron_RPKM value with the package of edge
REP2RPKM_Intron=as.data.frame(cbind(rownames(CountREP24SU_Intron),CountREP24SU_Intron,CountREP2Total_Intron))
colnames(REP2RPKM_Intron)<-c("GeneSymbol","REP24SU_IntronHour0","REP24SU_IntronHour4","REP24SU_IntronHour8","REP2Total_IntronHour0","REP2Total_IntronHour4","REP2Total_IntronHour8")

Intronwidth=as.data.frame(Intronwidth)
Intronwidths=cbind(rownames(Intronwidth),Intronwidth) 
colnames(Intronwidths)<-c("GeneSymbol","Length")

IntronALL<- sqldf("select * from REP2RPKM_Intron join Intronwidths where REP2RPKM_Intron.GeneSymbol= Intronwidths.GeneSymbol")
rawdata=IntronALL[,1:7]
write.table(rawdata, file = "d:/REP2RPKM_Intron.txt", append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

GeneLength=IntronALL[,8:9]
write.table(GeneLength, file = "d:/IntronFgeneLenth.txt", append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

rawdata=read.delim("d:/REP2RPKM_Intron.txt", check.names=FALSE, stringsAsFactors=FALSE)
GeneLength=read.delim("d:/IntronFgeneLenth.txt", check.names=FALSE, stringsAsFactors=FALSE)
GeneLength=GeneLength[,2]
group <- c("pre0","pre4","pre8","total0","total4","total8")
d<- DGEList(counts=rawdata[,2:7], genes=data.frame(rawdata[,1],Length=GeneLength),group=group)
d <- calcNormFactors(d)
REP2RPKM1 <- rpkm(d)
REP2INTRONRPKM<-cbind(rawdata[,1],REP2RPKM1)
write.table(REP2INTRONRPKM, "D:/REP2IntronRPKM.txt", sep="\t", row.names=F, col.names=T)

CountREP14SU_Exon<-REP1RPKM$counts$foursu$exonCounts
CountREP14SU_Intron<-REP1RPKM$counts$foursu$intronCounts
CountREP1Total_Exon<-REP1RPKM$counts$total$exonCounts
CountREP1Total_Intron<-REP1RPKM$counts$total$intronCounts

CountREP24SU_Exon<-REP2RPKM$counts$foursu$exonCounts
CountREP24SU_Intron<-REP2RPKM$counts$foursu$intronCounts
CountREP2Total_Exon<-REP2RPKM$counts$total$exonCounts
CountREP2Total_Intron<-REP2RPKM$counts$total$intronCounts

countEXON=as.data.frame(cbind(CountREP14SU_Exon,CountREP1Total_Exon,CountREP24SU_Exon,CountREP2Total_Exon))
colnames(countEXON)<-c("REP14SU_ExonHour0","REP14SU_ExonHour4","REP14SU_ExonHour8","REP1Total_ExonHour0","REP1Total_ExonHour4","REP1Total_ExonHour8"
                       ,"REP24SU_ExonHour0","REP24SU_ExonHour4","REP24SU_ExonHour8","REP2Total_ExonHour0","REP2Total_ExonHour4","REP2Total_ExonHour8")

EXONRPKM=as.matrix(countEXON)
ALLexons=EXONRPKM
#EXON___REP1和REP2合起来分析
EXONRPKM=cbind(REP1EXONRPKM,REP2EXONRPKM[,2:7])
rownames(EXONRPKM)=EXONRPKM[,1];
EXONRPKM=EXONRPKM[,2:13]

EXONRPKM[apply(EXONRPKM==0, FUN = all, 1), ] = NA
EXONRPKM[apply(EXONRPKM<10, FUN = all, 1), ] = NA
# EXONRPKM[apply(!is.finite(EXONRPKM), FUN = any, 1), ] = NA
ALLexons=na.omit(EXONRPKM)
# ALLexons=EXONRPKM

# 
# 对count的个数根据测序的深度进行标准化
# 

# model<-lm(ALLexons[,1]~ALLexons[,7])
# plot(ALLexons[,7],ALLexons[,1],xlim=c(0,10000),ylim=c(0,10000))
# abline(model,col="blue")
# summary(model)
# model$coefficients[[2]]
graphics.off()
r14E0=((as.double(ALLexons[,1])/7148281)*100000000*1.071)/as.double(unlist(Exonwidth))
r24E0=((as.double(ALLexons[,7])/6042502)*100000000)/as.double(unlist(Exonwidth))

model<-lm(r14E0~r24E0)
plot(r24E0,r14E0,xlim=c(0,300),ylim=c(0,300))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]

# model<-lm(ALLexons[,2]~ALLexons[,8])
# plot(ALLexons[,8],ALLexons[,2],xlim=c(0,10000),ylim=c(0,10000))
# abline(model,col="blue")
# summary(model)
# model$coefficients[[2]]

graphics.off()
r14E4=((as.double(ALLexons[,2])/7863883)*100000000)/as.double(unlist(Exonwidth))
r24E4=((as.double(ALLexons[,8])/9014666)*100000000)/as.double(unlist(Exonwidth))

model<-lm(r14E4~r24E4)
plot(r24E4,r14E4,xlim=c(0,100),ylim=c(0,100))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]



# model<-lm(ALLexons[,3]~ALLexons[,9])
# plot(ALLexons[,9],ALLexons[,3],xlim=c(0,10000),ylim=c(0,10000))
# abline(model,col="blue")
# summary(model)
# model$coefficients[[2]]
graphics.off()
r14E8=(as.double(ALLexons[,3])/6361178)*10000000*0.97
r24E8=(as.double(ALLexons[,9])/7297282)*10000000


model<-lm(r14E8~r24E8)
plot(r24E8,r14E8,xlim=c(0,10000),ylim=c(0,10000))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]



# model<-lm(ALLexons[,4]~ALLexons[,10])
# plot(ALLexons[,10],ALLexons[,4],xlim=c(0,10000),ylim=c(0,10000))
# abline(model,col="blue")
# summary(model)
# model$coefficients[[2]]

graphics.off()
r1TE0=(as.double(ALLexons[,4])/16433062)*100000000*0.98
r2TE0=(as.double(ALLexons[,10])/29329621)*100000000


model<-lm(r1TE0~r2TE0)
plot(r2TE0,r1TE0,xlim=c(0,10000),ylim=c(0,10000))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]



# model<-lm(ALLexons[,5]~ALLexons[,11])
# plot(ALLexons[,11],ALLexons[,5],xlim=c(0,10000),ylim=c(0,10000))
# abline(model,col="blue")
# summary(model)
# model$coefficients[[2]]

r1TE4=(as.double(ALLexons[,5])/24313241)*100000000*0.97
r2TE4=(as.double(ALLexons[,11])/29523912)*100000000


model<-lm(r1TE4~r2TE4)
plot(r2TE4,r1TE4,xlim=c(0,10000),ylim=c(0,10000))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]


# model<-lm(ALLexons[,6]~ALLexons[,12])
# plot(ALLexons[,12],ALLexons[,6],xlim=c(0,10000),ylim=c(0,10000))
# abline(model,col="blue")
# summary(model)
# model$coefficients[[2]]
graphics.off()
r1TE8=(as.double(ALLexons[,6])/11760276)*100000000
r2TE8=(as.double(ALLexons[,12])/27588260)*100000000*1.005


model<-lm(r1TE8~r2TE8)
plot(r2TE8,r1TE8,xlim=c(0,10000),ylim=c(0,10000))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]

write.table(ALLexons, "D:/EXONRPKM.txt", sep="\t", row.names=F, col.names=T)

#计算rep1和rep2两组的偏差
ALLexonsp=cbind(r14E0,r14E4,r14E8,r1TE0,r1TE4,r1TE8,r24E0,r24E4,r24E8,r2TE0,r2TE4,r2TE8,as.data.frame(r14E0/r24E0),
                as.data.frame(r14E4/r24E4),
                as.data.frame(r14E8/r24E8),
                as.data.frame(r1TE0/r2TE0),
                as.data.frame(r1TE4/r2TE4),
                as.data.frame(r1TE8/r2TE8)
                )
rownames(ALLexonsp)=rownames(ALLexons)
#删除在两个重复样本中变异差的基因
va=as.matrix(ALLexonsp[,13:18])
va[apply(!is.finite(va), FUN = any, 1), ] = NA
va[apply(va>1.5, FUN = any, 1), ] = NA
va[apply(va<0.75, FUN = any, 1), ] = NA

exonsp=cbind(ALLexonsp[,1:12],va)
exonsp=na.omit(exonsp)
write.table(exonsp, "D:/exonsp.txt", sep="\t", row.names=T, col.names=T)

model<-lm(CountREP2Total_Exon[,i]~CountREP1Total_Exon[,i])
plot(exonsp[,7],exonsp[,1],xlim=c(0,10000),ylim=c(0,10000))
plot(exonsp[,8],exonsp[,2],xlim=c(0,10000),ylim=c(0,10000))
plot(exonsp[,9],exonsp[,3],xlim=c(0,10000),ylim=c(0,10000))
plot(exonsp[,10],exonsp[,4],xlim=c(0,20000),ylim=c(0,20000))
plot(exonsp[,11],exonsp[,5],xlim=c(0,20000),ylim=c(0,20000))
plot(exonsp[,12],exonsp[,6],xlim=c(0,20000),ylim=c(0,20000))

exonsp=as.data.frame(exonsp)
#Intron___REP1和REP2合起来分析
INTRONRPKM=cbind(REP1INTRONRPKM,REP2INTRONRPKM[,2:7])
rownames(INTRONRPKM)=INTRONRPKM[,1];
INTRONRPKM=INTRONRPKM[,2:13]

INTRONRPKM[apply(INTRONRPKM==0, FUN = all, 1), ] = NA
INTRONRPKM[apply(INTRONRPKM<0.1, FUN = all, 1), ] = NA
# EXONRPKM[apply(!is.finite(EXONRPKM), FUN = any, 1), ] = NA
ALLIntron=na.omit(INTRONRPKM)
# ALLexons=EXONRPKM
# model<-lm(ALLexons[1:5000,1]~ALLexons[1:5000,7])
plot(ALLexons[,7],ALLexons[,1],xlim=c(0,2000),ylim=c(0,2000))
plot(ALLexons[,8],ALLexons[,2],xlim=c(0,2000),ylim=c(0,2000))
plot(ALLexons[,9],ALLexons[,3],xlim=c(0,2000),ylim=c(0,2000))
plot(ALLexons[,10],ALLexons[,4],xlim=c(0,2000),ylim=c(0,2000))
plot(ALLexons[,11],ALLexons[,5],xlim=c(0,2000),ylim=c(0,2000))
plot(ALLexons[,12],ALLexons[,6],xlim=c(0,2000),ylim=c(0,2000))
abline(model,col="blue")

write.table(ALLexons, "D:/EXONRPKM.txt", sep="\t", row.names=F, col.names=T)

#计算rep1和rep2两组的偏差
ALLintronsp=cbind(ALLIntron,as.data.frame(as.numeric(ALLIntron[,1]) /as.numeric(ALLIntron[,7])),
                as.data.frame(as.numeric(ALLIntron[,2]) /as.numeric(ALLIntron[,8])),
                as.data.frame(as.numeric(ALLIntron[,3]) /as.numeric(ALLIntron[,9])),
                as.data.frame(as.numeric(ALLIntron[,4]) /as.numeric(ALLIntron[,10])),
                as.data.frame(as.numeric(ALLIntron[,5]) /as.numeric(ALLIntron[,11])),
                as.data.frame(as.numeric(ALLIntron[,6]) /as.numeric(ALLIntron[,12]))
)
#删除在两个重复样本中变异差的基因
va=as.matrix(ALLintronsp[,13:18])
va[apply(!is.finite(va), FUN = any, 1), ] = NA
va[apply(va>1.2, FUN = any, 1), ] = NA
va[apply(va<0.85, FUN = any, 1), ] = NA

intronsp=cbind(ALLintronsp,va)
intronsp=na.omit(intronsp)
write.table(exonsp, "D:/intronsp.txt", sep="\t", row.names=T, col.names=T)
intronsp=as.matrix(intronsp)
model<-lm(CountREP2Total_Exon[,i]~CountREP1Total_Exon[,i])
plot(intronsp[,7],intronsp[,1],xlim=c(0,200),ylim=c(0,200))
plot(intronsp[,8],intronsp[,2],xlim=c(0,200),ylim=c(0,200))
plot(intronsp[,9],intronsp[,3],xlim=c(0,200),ylim=c(0,200))
plot(intronsp[,10],intronsp[,4],xlim=c(0,200),ylim=c(0,200))
plot(intronsp[,11],intronsp[,5],xlim=c(0,200),ylim=c(0,200))
plot(intronsp[,12],intronsp[,6],xlim=c(0,200),ylim=c(0,200))


##
exonsp=cbind(rownames(exonsp),exonsp)
colnames(exonsp)[1]=c("GeneSymbol")
exonsp=as.data.frame(exonsp)
colnames(exonsp)[14:19]=c("suexonhour0","suexonhour4","suexonhour8","totalexonhour0","totalexonhour4","totalexonhour8")

intronsp=cbind(rownames(intronsp),intronsp)
colnames(intronsp)[1]=c("GeneSymbol")
intronsp=as.data.frame(intronsp)
colnames(intronsp)[14:19]=c("suintronhour0","suintronhour4","suintronhour8","totalintronhour0","totalintronhour4","totalintronhour8")

ALL<- sqldf("select * from exonsp join intronsp where exonsp.GeneSymbol= intronsp .GeneSymbol")

write.table(ALL, "D:/ALLRPKM85_12.txt", sep="\t", row.names=F, col.names=T)










# for(i in 1:3){
# png(file=paste("D:/TotalExon",i,".png"))
# par(mfrow=c(1,2))
# model<-lm(CountREP2Total_Exon[,i]~CountREP1Total_Exon[,i])
# plot(CountREP1Total_Exon[,i],CountREP2Total_Exon[,i])
# abline(model,col="blue")
# 
# model<-lm(REP2Total_Exon[,i]~REP1Total_Exon[,i])
# plot(REP1Total_Exon[,i],REP2Total_Exon[,i])
# abline(model,col="blue")
# 
# dev.off()
# }
# 
# graphics.off()
# for(i in 1:3){
# png(file=paste("D:/TotalIntron",i,".png"))
# par(mfrow=c(1,2))
# model<-lm(CountREP2Total_Intron[,i]~CountREP1Total_Intron[,i])
# plot(CountREP1Total_Intron[,i],CountREP2Total_Intron[,i])
# abline(model,col="blue")
# 
# model<-lm(REP2Total_Intron[,i]~REP1Total_Intron[,i])
# plot(REP1Total_Intron[,i],REP2Total_Intron[,i])
# abline(model,col="blue")
# dev.off()
# }
