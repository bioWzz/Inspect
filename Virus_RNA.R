rm(list=ls())
library(limma)
library(edgeR)
library(sqldf)
##Produce the RPKM value and GENE EXPRESSION ANALYSIS
path="D:/ivan/Harm/NewTestFour/my/Results/VirusInput/limma/"
#Targets <- read.delim(paste(path,"EXOSC3_RNA_labeling_design.txt",sep=""), sep = "\t")
Counts=read.delim(paste(path,"EXOSC3_RNA_labeling_all_counts.txt",sep=""), sep = "\t")
Counts=cbind(row.names(Counts),Counts)
Counts=as.data.frame(Counts[order(Counts[,1]),])
Counts=Counts[,2:25]
dim(Counts)
y <- DGEList(counts = Counts)
y <- calcNormFactors(y)
logcpm <- cpm(y, prior.count=2, log=TRUE)

Geneanno<-read.table(paste(path,"gencode.v23.PR8.anno",sep=""), header = TRUE, sep = "\t")
Geneanno<- as.vector(Geneanno[,c(1,3)])
Geneanno=Geneanno[order(Geneanno[,1]),]
Geneanno=as.vector(Geneanno[,2])

RPKM=as.data.frame(rpkm(y, gene.length=Geneanno))
write.table(RPKM,file=paste(path,"VirusExpression.txt",sep=""),append = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
RPKM=RPKM[c(60499:60514),]
Vrna=RPKM[c(1,3,5,7,9,11,13,15),]
Vmrna=RPKM[c(2,4,6,8,10,12,14,16),]
v=rbind(Vrna,Vmrna)
## Expression analysis of compare LabA and LabC
Hour4=v[,c("LabA.4hpi","LabA.4_2","LabC.4hpi","LabC.4_2","UnlA.4hpi","UnlabA.4_2","UnlC.4hpi","UnlabC.4_2")]
Hour8=v[,c("LabA.8hpi","LabA.8_2","LabC.8hpi","LabC.8_2","UnlA.8hpi","UnlabA.8_2","UnlC.8hpi","UnlabC.8_2")]
write.table(Hour4, file="D:/ivan/Harm/NewTestFour/my/Results/VirusInput/ExpressionAnalysis/Hour4.txt",row.names=T, sep="\t", quote=FALSE)
write.table(Hour8, file="D:/ivan/Harm/NewTestFour/my/Results/VirusInput/ExpressionAnalysis/Hour8.txt",row.names=T, sep="\t", quote=FALSE)

## Expression analysis of compare mRNA and RNA
LabA4=cbind(Vrna$LabA.4hpi,Vrna$LabA.4_2,Vmrna$LabA.4hpi,Vmrna$LabA.4_2,Vrna$UnlA.4hpi,Vrna$UnlabA.4_2,Vmrna$UnlA.4hpi,Vmrna$UnlabA.4_2)
colnames(LabA4)<-c("R1_RNA_LabA4","R2_RNA_LabA4","R1_mRNA_LabA4","R2_mRNA_LabA4","R1_RNA_ULabA4","R2_RNA_ULabA4","R1_mRNA_ULabA4","R2_mRNA_ULabA4")
row.names(LabA4)<-row.names(Vrna)

LabA8=cbind(Vrna$LabA.8hpi,Vrna$LabA.8_2,Vmrna$LabA.8hpi,Vmrna$LabA.8_2,Vrna$UnlA.8hpi,Vrna$UnlabA.8_2,Vmrna$UnlA.8hpi,Vmrna$UnlabA.8_2)
colnames(LabA8)<-c("R1_RNA_LabA8","R2_RNA_LabA8","R1_mRNA_LabA8","R2_mRNA_LabA8","R1_RNA_ULabA8","R2_RNA_ULabA8","R1_mRNA_ULabA8","R2_mRNA_ULabA8")
row.names(LabA8)<-row.names(Vrna)

LabC4=cbind(Vrna$LabC.4hpi,Vrna$LabC.4_2,Vmrna$LabC.4hpi,Vmrna$LabC.4_2,Vrna$UnlC.4hpi,Vrna$UnlabC.4_2,Vmrna$UnlC.4hpi,Vmrna$UnlabC.4_2)
colnames(LabC4)<-c("R1_RNA_LabC4","R2_RNA_LabC4","R1_mRNA_LabC4","R2_mRNA_LabC4","R1_RNA_ULabC4","R2_RNA_ULabC4","R1_mRNA_ULabC4","R2_mRNA_ULabC4")
row.names(LabC4)<-row.names(Vrna)

LabC8=cbind(Vrna$LabC.8hpi,Vrna$LabC.8_2,Vmrna$LabC.8hpi,Vmrna$LabC.8_2,Vrna$UnlC.8hpi,Vrna$UnlabC.8_2,Vmrna$UnlC.8hpi,Vmrna$UnlabC.8_2)
colnames(LabC8)<-c("R1_RNA_LabC8","R2_RNA_LabC8","R1_mRNA_LabC8","R2_mRNA_LabC8","R1_RNA_ULabC8","R2_RNA_ULabC8","R1_mRNA_ULabC8","R2_mRNA_ULabC8")
row.names(LabC8)<-row.names(Vrna)

write.table(LabA4, file="D:/ivan/Harm/NewTestFour/my/Results/VirusInput/ExpressionAnalysis/LabA4.txt",row.names=T, sep="\t", quote=FALSE)
write.table(LabA8, file="D:/ivan/Harm/NewTestFour/my/Results/VirusInput/ExpressionAnalysis/LabA8.txt",row.names=T, sep="\t", quote=FALSE)
write.table(LabC4, file="D:/ivan/Harm/NewTestFour/my/Results/VirusInput/ExpressionAnalysis/LabC4.txt",row.names=T, sep="\t", quote=FALSE)
write.table(LabC8, file="D:/ivan/Harm/NewTestFour/my/Results/VirusInput/ExpressionAnalysis/LabC8.txt",row.names=T, sep="\t", quote=FALSE)


##caculate the rate
RPKM=as.data.frame(rpkm(y, gene.length=Geneanno))
RPKM=RPKM[c(60499:60514),]
Vrna=RPKM[c(1,3,5,7,9,11,13,15),]
Vmrna=RPKM[c(2,4,6,8,10,12,14,16),]
v=rbind(Vrna,Vmrna)
write.table(v,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/virsua-Result/RPKM.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)
library(GenomicFeatures)
library(INSPEcT)
library(Rsamtools)
colnames(v)
v=RPKM
#LabA-Rep1
Rep1foursu_exons=v[,c("LabA.0hpi","LabA.4hpi","LabA.8hpi")]
colnames(Rep1foursu_exons) <-c("RA1Onehour","RA1Fourhour","RA1Eighthour")


Rep1foursu_introns=Rep1foursu_exons
Rep1foursu_introns[,]=0
colnames(Rep1foursu_introns) <-c("R1A0","R1A4","R1A8")

Rep1total_exons=v[,c("UnlA.0hpi","UnlA.4hpi","UnlA.8hpi")]
colnames(Rep1total_exons) <-c("RA1Onehour","RA1Fourhour","RA1Eighthour")

Rep1total_introns=Rep1total_exons
Rep1total_introns[,]=0
colnames(Rep1total_introns) <-c("RA1Onehour","RA1Fourhour","RA1Eighthour")

tpts <- c(0,4,8)
tL <- 4
VirusalRNA <- newINSPEcT(tpts, tL,Rep1foursu_exons, Rep1total_exons,
                         Rep1foursu_introns, Rep1total_introns, BPPARAM=SerialParam())
VirusalRNALabASynthesis=(ratesFirstGuess(VirusalRNA, 'synthesis'))
VirusalRNALabADegradation=(ratesFirstGuess(VirusalRNA, 'degradation'))
VirusalRNALabAProcessing=(ratesFirstGuess(VirusalRNA, 'processing'))

LabARep1=cbind(VirusalRNALabASynthesis,VirusalRNALabADegradation,VirusalRNALabAProcessing)
write.table(LabARep1,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/virsua-Result/LabARep1.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

re=cbind()

#LabA-Rep2
Rep2foursu_exons=v[,c("LabA.0_2","LabA.4_2","LabA.8_2")]
colnames(Rep2foursu_exons) <-c("Onehour","Fourhour","Eighthour")


Rep2foursu_introns=Rep2foursu_exons
Rep2foursu_introns[,]=0
colnames(Rep2foursu_introns) <-c("Onehour","Fourhour","Eighthour")

Rep2total_exons=v[,c("UnlabA.0_2","UnlabA.4_2","UnlabA.8_2")]
colnames(Rep2total_exons) <-c("Onehour","Fourhour","Eighthour")

Rep2total_introns=Rep2total_exons
Rep2total_introns[,]=0
colnames(Rep2total_introns) <-c("Onehour","Fourhour","Eighthour")

tpts <- c(0,4,8)
tL <- 4
VirusalRNA <- newINSPEcT(tpts, tL, Rep2foursu_exons, Rep2total_exons,
                         Rep2foursu_introns, Rep2total_introns, BPPARAM=SerialParam())
VirusalRNALabASynthesis=(ratesFirstGuess(VirusalRNA, 'synthesis'))
VirusalRNALabADegradation=(ratesFirstGuess(VirusalRNA, 'degradation'))
VirusalRNALabAProcessing=(ratesFirstGuess(VirusalRNA, 'processing'))

LabARep2=cbind(VirusalRNALabASynthesis,VirusalRNALabADegradation,VirusalRNALabAProcessing)
write.table(LabARep2,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/virsua-Result/LabARep2.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)



tpts <- c(0,4,8,0,4,8)
tL <- 4
foursu_exons=cbind(Rep1foursu_exons,Rep2foursu_exons)
colnames(foursu_exons) <-c("foursu_exonsREP1Onehour","foursu_exonsREP1Fourhour","foursu_exonsREP1Eighthour","foursu_exonsREP2Onehour","foursu_exonsREP2Fourhour","foursu_exonsREP2Eighthour")

total_exons=cbind(Rep1total_exons,Rep2total_exons)
colnames(total_exons) <-c("total_exonsREP1Onehour","total_exonsREP1Fourhour","total_exonsREP1Eighthour","total_exonsREP2Onehour","total_exonsREP2Fourhour","total_exonsREP2Eighthour")

foursu_introns=cbind(Rep1foursu_introns,Rep2foursu_introns)
colnames(foursu_introns) <-c("foursu_intronsREP1Onehour","foursu_intronsREP1Fourhour","foursu_intronsREP1Eighthour","foursu_intronsREP2Onehour","foursu_intronsREP2Fourhour","foursu_intronsREP2Eighthour")

total_introns=cbind(Rep1total_introns,Rep2total_introns)
colnames(total_introns) <-c("total_intronsREP1Onehour","total_intronsREP1Fourhour","total_intronsREP1Eighthour","total_intronsREP2Onehour","total_intronsREP2Fourhour","total_intronsREP2Eighthour")

virusLabA12<-newINSPEcT(tpts,tL,foursu_exons, total_exons,foursu_introns,total_introns,BPPARAM=SerialParam())
LabAREP12Express=LabAREP12@ratesFirstGuess@assayData$exprs


#Labc-Rep1
foursu_exons=v[,c("LabC.0hpi","LabC.4hpi","LabC.8hpi")]
colnames(foursu_exons) <-c("Onehour","Fourhour","Eighthour")


foursu_introns=foursu_exons
foursu_introns[,]=0
colnames(foursu_introns) <-c("Onehour","Fourhour","Eighthour")

total_exons=v[,c("UnlC.0hpi","UnlC.4hpi","UnlC.8hpi")]
colnames(total_exons) <-c("Onehour","Fourhour","Eighthour")

total_introns=total_exons
total_introns[,]=0
colnames(total_introns) <-c("Onehour","Fourhour","Eighthour")

tpts <- c(0,4,8)
tL <- 4
VirusalRNA <- newINSPEcT(tpts, tL, foursu_exons, total_exons,
                         foursu_introns, total_introns, BPPARAM=SerialParam())
VirusalRNALabASynthesis=(ratesFirstGuess(VirusalRNA, 'synthesis'))
VirusalRNALabADegradation=(ratesFirstGuess(VirusalRNA, 'degradation'))
VirusalRNALabAProcessing=(ratesFirstGuess(VirusalRNA, 'processing'))

LabARep2=cbind(VirusalRNALabASynthesis,VirusalRNALabADegradation,VirusalRNALabAProcessing)
write.table(LabARep2,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/virsua-Result/LabCRep1.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

#Labc-Rep2
foursu_exons=v[,c("LabC.0_2","LabC.4_2","LabC.8_2")]
colnames(foursu_exons) <-c("Onehour","Fourhour","Eighthour")


foursu_introns=foursu_exons
foursu_introns[,]=0
colnames(foursu_introns) <-c("Onehour","Fourhour","Eighthour")

total_exons=v[,c("UnlabC.0_2","UnlabC.4_2","UnlabC.8_2")]
colnames(total_exons) <-c("Onehour","Fourhour","Eighthour")

total_introns=total_exons
total_introns[,]=0
colnames(total_introns) <-c("Onehour","Fourhour","Eighthour")

tpts <- c(0,4,8)
tL <- 4
VirusalRNA <- newINSPEcT(tpts, tL, foursu_exons, total_exons,
                         foursu_introns, total_introns, BPPARAM=SerialParam())
VirusalRNALabASynthesis=(ratesFirstGuess(VirusalRNA, 'synthesis'))
VirusalRNALabADegradation=(ratesFirstGuess(VirusalRNA, 'degradation'))
VirusalRNALabAProcessing=(ratesFirstGuess(VirusalRNA, 'processing'))

LabARep2=cbind(VirusalRNALabASynthesis,VirusalRNALabADegradation,VirusalRNALabAProcessing)
write.table(LabARep2,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/virsua-Result/LabCRep2.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)


# tpts <- c(0,4,8,0,4,8)
# tL <- 4
# foursu_exons=cbind(REP1_foursu_exons,Rep2_foursu_exons)
# colnames(foursu_exons) <-c("foursu_exonsREP1Onehour","foursu_exonsREP1Fourhour","foursu_exonsREP1Eighthour","foursu_exonsREP2Onehour","foursu_exonsREP2Fourhour","foursu_exonsREP2Eighthour")
# 
# total_exons=cbind(REP1_total_exons, Rep2_total_exons)
# colnames(total_exons) <-c("total_exonsREP1Onehour","total_exonsREP1Fourhour","total_exonsREP1Eighthour","total_exonsREP2Onehour","total_exonsREP2Fourhour","total_exonsREP2Eighthour")
# 
# foursu_introns=cbind(REP1_foursu_introns,Rep2_foursu_introns)
# colnames(foursu_introns) <-c("foursu_intronsREP1Onehour","foursu_intronsREP1Fourhour","foursu_intronsREP1Eighthour","foursu_intronsREP2Onehour","foursu_intronsREP2Fourhour","foursu_intronsREP2Eighthour")
# 
# total_introns=cbind(REP1_total_introns,Rep2_total_introns)
# colnames(total_introns) <-c("total_intronsREP1Onehour","total_intronsREP1Fourhour","total_intronsREP1Eighthour","total_intronsREP2Onehour","total_intronsREP2Fourhour","total_intronsREP2Eighthour")
# 
# LabAREP12<-newINSPEcT(tpts,tL,foursu_exons, total_exons,foursu_introns,total_introns,BPPARAM=SerialParam())
# LabAREP12Express=LabAREP12@ratesFirstGuess@assayData$exprs











