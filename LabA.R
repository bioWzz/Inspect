rm(list=ls())
library(GenomicFeatures)
library(INSPEcT)
library(Rsamtools)
library(deSolve)
library(BiocParallel)
library(compiler)
library(sqldf)
load("D:/ivan/Harm/NewTestFour/my/Result/txdb.Rdata")
# Gene
# 
 ResultDir="D:/ivan/Harm/NewTestFour/my/Result/"

txdb3 <- makeTxDbFromGFF(file="D:/ivan/Harm/NewTestOne/NewTest/Newgencode.v23.PR8.gtf",format="gtf",taxonomyId=9606)
 # save(txdb3, file = paste(ResultDir,"txdb.Rdata",sep=""))
 load(paste(ResultDir,"txdb.Rdata",sep=""))
#### replication one tarnscprition.bam #######
# path=D:ProgramFilesRRR\R-3.3.1\library\INSPEcT\extdata  Rep1\Genome
###################################################################################################################
#####REP1---LabA---------------------------------------------------------------------------------
###calcutation of Expression of LabA
# LabA-0
paths_4su <- system.file('extdata/MyRep1', 'LabA0hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep1', 'UnlA0hpi.bam', package="INSPEcT")
Rep1LabA0hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)


# LabA-4
paths_4su <- system.file('extdata/MyRep1', 'LabA4hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep1', 'UnlA4hpi.bam', package="INSPEcT")
Rep1LabA4hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)

# LabA-8
paths_4su <- system.file('extdata/MyRep1', 'LabA8hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep1', 'UnlA8hpi.bam', package="INSPEcT")
Rep1LabA8hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)

######################
##LabA rate
LabA0hpifoursu_exons=Rep1LabA0hpiRPKMsOut$rpkms$foursu_exons
LabA4hpifoursu_exons=Rep1LabA4hpiRPKMsOut$rpkms$foursu_exons
LabA8hpifoursu_exons=Rep1LabA8hpiRPKMsOut$rpkms$foursu_exons

REP1_foursu_exons=cbind(LabA0hpifoursu_exons,LabA4hpifoursu_exons,LabA8hpifoursu_exons)
colnames(REP1_foursu_exons) <-c("Onehour","Fourhour","Eighthour")

LabA0hpifoursu_introns=Rep1LabA0hpiRPKMsOut$rpkms$foursu_introns
LabA4hpifoursu_introns=Rep1LabA4hpiRPKMsOut$rpkms$foursu_introns
LabA8hpifoursu_introns=Rep1LabA8hpiRPKMsOut$rpkms$foursu_introns

REP1_foursu_introns=cbind( LabA0hpifoursu_introns,LabA4hpifoursu_introns,LabA8hpifoursu_introns)
colnames(REP1_foursu_introns) <-c("Onehour","Fourhour","Eighthour")


# UnlA0hpitotal_exons=Rep1LabA0hpiRPKMsOut$rpkms$total_exons
UnlA0hpitotal_exons=Rep1LabA0hpiRPKMsOut$rpkms$total_exons
UnlA4hpitotal_exons=Rep1LabA4hpiRPKMsOut$rpkms$total_exons
UnlA8hpitotal_exons=Rep1LabA8hpiRPKMsOut$rpkms$total_exons
REP1_total_exons=cbind(UnlA0hpitotal_exons,UnlA4hpitotal_exons,UnlA8hpitotal_exons)
colnames(REP1_total_exons) <-c("Onehour","Fourhour","Eighthour")

UnlA0hpitotal_introns=Rep1LabA0hpiRPKMsOut$rpkms$total_introns
UnlA4hpitotal_introns=Rep1LabA4hpiRPKMsOut$rpkms$total_introns
UnlA8hpitotal_introns=Rep1LabA8hpiRPKMsOut$rpkms$total_introns
REP1_total_introns=cbind(UnlA0hpitotal_introns,UnlA4hpitotal_introns,UnlA8hpitotal_introns)
colnames(REP1_total_introns) <-c("Onehour","Fourhour","Eighthour")

tpts <- c(0,4,8)
tL <- 4
GeneLabAREP1 <- newINSPEcT(tpts, tL, REP1_foursu_exons, REP1_total_exons,
                   REP1_foursu_introns, REP1_total_introns, BPPARAM=SerialParam())
GeneLabAREP1Express=GeneLabAREP1@ratesFirstGuess@assayData$exprs





# # ###################################################################################################################
# #####REP2---LabA---------------------------------------------------------------------------------
ResultDir="D:/ivan/Harm/NewTestFour/my/Result/"
# ###calcutation of Expression of LabA
# LabA-0
paths_4su <- system.file('extdata/MyRep2', 'LabA0.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep2', 'ULabA0.bam', package="INSPEcT")
Rep2LabA0hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)


# LabA-4
paths_4su <- system.file('extdata/MyRep2', 'LabA4.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep2', 'ULabA4.bam', package="INSPEcT")
Rep2LabA4hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)


# LabA-8
paths_4su <- system.file('extdata/MyRep2', 'LabA8.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep2', 'ULabA8.bam', package="INSPEcT")
Rep2LabA8hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)


######################
##LabA rate
LabA0hpifoursu_exons=Rep2LabA0hpiRPKMsOut$rpkms$foursu_exons
LabA4hpifoursu_exons=Rep2LabA4hpiRPKMsOut$rpkms$foursu_exons
LabA8hpifoursu_exons=Rep2LabA8hpiRPKMsOut$rpkms$foursu_exons
Rep2_foursu_exons=cbind(LabA0hpifoursu_exons,LabA4hpifoursu_exons,LabA8hpifoursu_exons)
colnames(Rep2_foursu_exons) <-c("Onehour","Fourhour","Eighthour")

LabA0hpifoursu_introns=Rep2LabA0hpiRPKMsOut$rpkms$foursu_introns
LabA4hpifoursu_introns=Rep2LabA4hpiRPKMsOut$rpkms$foursu_introns
LabA8hpifoursu_introns=Rep2LabA8hpiRPKMsOut$rpkms$foursu_introns
Rep2_foursu_introns=cbind(LabA0hpifoursu_introns,LabA4hpifoursu_introns,LabA8hpifoursu_introns)
colnames(Rep2_foursu_introns) <-c("Onehour","Fourhour","Eighthour")

UnlA0hpitotal_exons=Rep2LabA0hpiRPKMsOut$rpkms$total_exons
UnlA4hpitotal_exons=Rep2LabA4hpiRPKMsOut$rpkms$total_exons
UnlA8hpitotal_exons=Rep2LabA8hpiRPKMsOut$rpkms$total_exons
Rep2_total_exons=cbind(UnlA0hpitotal_exons,UnlA4hpitotal_exons,UnlA8hpitotal_exons)
colnames(Rep2_total_exons) <-c("Onehour","Fourhour","Eighthour")

UnlA0hpitotal_introns=Rep2LabA0hpiRPKMsOut$rpkms$total_introns
UnlA4hpitotal_introns=Rep2LabA4hpiRPKMsOut$rpkms$total_introns
UnlA8hpitotal_introns=Rep2LabA8hpiRPKMsOut$rpkms$total_introns
Rep2_total_introns=cbind(UnlA0hpitotal_introns,UnlA4hpitotal_introns,UnlA8hpitotal_introns)
colnames(Rep2_total_introns) <-c("Onehour","Fourhour","Eighthour")

GeneLabAREP2 <- newINSPEcT(tpts, tL, Rep2_foursu_exons, Rep2_total_exons,
                   Rep2_foursu_introns, Rep2_total_introns, BPPARAM=SerialParam())
GeneLabAREP2Express=GeneLabAREP2@ratesFirstGuess@assayData$exprs

#Rep1+Rep2
tpts <- c(0,4,8,0,4,8)
tL <- 4
foursu_exons=cbind(REP1_foursu_exons,Rep2_foursu_exons)
colnames(foursu_exons) <-c("foursu_exonsREP1Onehour","foursu_exonsREP1Fourhour","foursu_exonsREP1Eighthour","foursu_exonsREP2Onehour","foursu_exonsREP2Fourhour","foursu_exonsREP2Eighthour")

total_exons=cbind(REP1_total_exons, Rep2_total_exons)
colnames(total_exons) <-c("total_exonsREP1Onehour","total_exonsREP1Fourhour","total_exonsREP1Eighthour","total_exonsREP2Onehour","total_exonsREP2Fourhour","total_exonsREP2Eighthour")

foursu_introns=cbind(REP1_foursu_introns,Rep2_foursu_introns)
colnames(foursu_introns) <-c("foursu_intronsREP1Onehour","foursu_intronsREP1Fourhour","foursu_intronsREP1Eighthour","foursu_intronsREP2Onehour","foursu_intronsREP2Fourhour","foursu_intronsREP2Eighthour")

total_introns=cbind(REP1_total_introns,Rep2_total_introns)
colnames(total_introns) <-c("total_intronsREP1Onehour","total_intronsREP1Fourhour","total_intronsREP1Eighthour","total_intronsREP2Onehour","total_intronsREP2Fourhour","total_intronsREP2Eighthour")

GeneLabAREP12<-newINSPEcT(tpts,tL,foursu_exons, total_exons,foursu_introns,total_introns,BPPARAM=SerialParam())
GeneLabAREP12Express=GeneLabAREP12@ratesFirstGuess@assayData$exprs

# 
# TX
# 
# ResultDir="D:/ivan/Harm/NewTestFour/my/Result/"

txdb3 <- makeTxDbFromGFF(file="D:/ivan/Harm/NewTestOne/NewTest/Newgencode.v23.PR8.gtf",format="gtf",taxonomyId=9606)
# save(txdb3, file = paste(ResultDir,"txdb.Rdata",sep=""))
# load(paste(ResultDir,"txdb.Rdata",sep=""))
#### replication one tarnscprition.bam #######
# path=D:ProgramFilesRRR\R-3.3.1\library\INSPEcT\extdata  Rep1\Genome
###################################################################################################################
#####REP1---LabA---------------------------------------------------------------------------------
###calcutation of Expression of LabA
# LabA-0
paths_4su <- system.file('extdata/MyRep1', 'LabA0hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep1', 'UnlA0hpi.bam', package="INSPEcT")
Rep1LabA0hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total,by ="tx")


# LabA-4
paths_4su <- system.file('extdata/MyRep1', 'LabA4hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep1', 'UnlA4hpi.bam', package="INSPEcT")
Rep1LabA4hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total,by ="tx")

# LabA-8
paths_4su <- system.file('extdata/MyRep1', 'LabA8hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep1', 'UnlA8hpi.bam', package="INSPEcT")
Rep1LabA8hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total,by ="tx")

######################
##LabA rate
LabA0hpifoursu_exons=Rep1LabA0hpiRPKMsOut$rpkms$foursu_exons
LabA4hpifoursu_exons=Rep1LabA4hpiRPKMsOut$rpkms$foursu_exons
LabA8hpifoursu_exons=Rep1LabA8hpiRPKMsOut$rpkms$foursu_exons

REP1_foursu_exons=cbind(LabA0hpifoursu_exons,LabA4hpifoursu_exons,LabA8hpifoursu_exons)
colnames(REP1_foursu_exons) <-c("Onehour","Fourhour","Eighthour")

LabA0hpifoursu_introns=Rep1LabA0hpiRPKMsOut$rpkms$foursu_introns
LabA4hpifoursu_introns=Rep1LabA4hpiRPKMsOut$rpkms$foursu_introns
LabA8hpifoursu_introns=Rep1LabA8hpiRPKMsOut$rpkms$foursu_introns

REP1_foursu_introns=cbind( LabA0hpifoursu_introns,LabA4hpifoursu_introns,LabA8hpifoursu_introns)
colnames(REP1_foursu_introns) <-c("Onehour","Fourhour","Eighthour")


# UnlA0hpitotal_exons=Rep1LabA0hpiRPKMsOut$rpkms$total_exons
UnlA0hpitotal_exons=Rep1LabA0hpiRPKMsOut$rpkms$total_exons
UnlA4hpitotal_exons=Rep1LabA4hpiRPKMsOut$rpkms$total_exons
UnlA8hpitotal_exons=Rep1LabA8hpiRPKMsOut$rpkms$total_exons
REP1_total_exons=cbind(UnlA0hpitotal_exons,UnlA4hpitotal_exons,UnlA8hpitotal_exons)
colnames(REP1_total_exons) <-c("Onehour","Fourhour","Eighthour")

UnlA0hpitotal_introns=Rep1LabA0hpiRPKMsOut$rpkms$total_introns
UnlA4hpitotal_introns=Rep1LabA4hpiRPKMsOut$rpkms$total_introns
UnlA8hpitotal_introns=Rep1LabA8hpiRPKMsOut$rpkms$total_introns
REP1_total_introns=cbind(UnlA0hpitotal_introns,UnlA4hpitotal_introns,UnlA8hpitotal_introns)
colnames(REP1_total_introns) <-c("Onehour","Fourhour","Eighthour")

tpts <- c(0,4,8)
tL <- 4
TXLabAREP1 <- newINSPEcT(tpts, tL, REP1_foursu_exons, REP1_total_exons,
                   REP1_foursu_introns, REP1_total_introns, BPPARAM=SerialParam())
TXLabAREP1Express=TXLabAREP1@ratesFirstGuess@assayData$exprs





# # ###################################################################################################################
# #####REP2---LabA---------------------------------------------------------------------------------
ResultDir="D:/ivan/Harm/NewTestFour/my/Result/"
# ###calcutation of Expression of LabA
# LabA-0
paths_4su <- system.file('extdata/MyRep2', 'LabA0.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep2', 'ULabA0.bam', package="INSPEcT")
Rep2LabA0hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total,by ="tx")


# LabA-4
paths_4su <- system.file('extdata/MyRep2', 'LabA4.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep2', 'ULabA4.bam', package="INSPEcT")
Rep2LabA4hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total,by ="tx")


# LabA-8
paths_4su <- system.file('extdata/MyRep2', 'LabA8.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep2', 'ULabA8.bam', package="INSPEcT")
Rep2LabA8hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total,by ="tx")


######################
##LabA rate
LabA0hpifoursu_exons=Rep2LabA0hpiRPKMsOut$rpkms$foursu_exons
LabA4hpifoursu_exons=Rep2LabA4hpiRPKMsOut$rpkms$foursu_exons
LabA8hpifoursu_exons=Rep2LabA8hpiRPKMsOut$rpkms$foursu_exons
Rep2_foursu_exons=cbind(LabA0hpifoursu_exons,LabA4hpifoursu_exons,LabA8hpifoursu_exons)
colnames(Rep2_foursu_exons) <-c("Onehour","Fourhour","Eighthour")

LabA0hpifoursu_introns=Rep2LabA0hpiRPKMsOut$rpkms$foursu_introns
LabA4hpifoursu_introns=Rep2LabA4hpiRPKMsOut$rpkms$foursu_introns
LabA8hpifoursu_introns=Rep2LabA8hpiRPKMsOut$rpkms$foursu_introns
Rep2_foursu_introns=cbind(LabA0hpifoursu_introns,LabA4hpifoursu_introns,LabA8hpifoursu_introns)
colnames(Rep2_foursu_introns) <-c("Onehour","Fourhour","Eighthour")

UnlA0hpitotal_exons=Rep2LabA0hpiRPKMsOut$rpkms$total_exons
UnlA4hpitotal_exons=Rep2LabA4hpiRPKMsOut$rpkms$total_exons
UnlA8hpitotal_exons=Rep2LabA8hpiRPKMsOut$rpkms$total_exons
Rep2_total_exons=cbind(UnlA0hpitotal_exons,UnlA4hpitotal_exons,UnlA8hpitotal_exons)
colnames(Rep2_total_exons) <-c("Onehour","Fourhour","Eighthour")

UnlA0hpitotal_introns=Rep2LabA0hpiRPKMsOut$rpkms$total_introns
UnlA4hpitotal_introns=Rep2LabA4hpiRPKMsOut$rpkms$total_introns
UnlA8hpitotal_introns=Rep2LabA8hpiRPKMsOut$rpkms$total_introns
Rep2_total_introns=cbind(UnlA0hpitotal_introns,UnlA4hpitotal_introns,UnlA8hpitotal_introns)
colnames(Rep2_total_introns) <-c("Onehour","Fourhour","Eighthour")

TXLabAREP2 <- newINSPEcT(tpts, tL, Rep2_foursu_exons, Rep2_total_exons,
                   Rep2_foursu_introns, Rep2_total_introns, BPPARAM=SerialParam())
TXLabAREP2Express=TXLabAREP2@ratesFirstGuess@assayData$exprs

#Rep1+Rep2
tpts <- c(0,4,8,0,4,8)
tL <- 4
foursu_exons=cbind(REP1_foursu_exons,Rep2_foursu_exons)
colnames(foursu_exons) <-c("foursu_exonsREP1Onehour","foursu_exonsREP1Fourhour","foursu_exonsREP1Eighthour","foursu_exonsREP2Onehour","foursu_exonsREP2Fourhour","foursu_exonsREP2Eighthour")

total_exons=cbind(REP1_total_exons, Rep2_total_exons)
colnames(total_exons) <-c("total_exonsREP1Onehour","total_exonsREP1Fourhour","total_exonsREP1Eighthour","total_exonsREP2Onehour","total_exonsREP2Fourhour","total_exonsREP2Eighthour")

foursu_introns=cbind(REP1_foursu_introns,Rep2_foursu_introns)
colnames(foursu_introns) <-c("foursu_intronsREP1Onehour","foursu_intronsREP1Fourhour","foursu_intronsREP1Eighthour","foursu_intronsREP2Onehour","foursu_intronsREP2Fourhour","foursu_intronsREP2Eighthour")

total_introns=cbind(REP1_total_introns,Rep2_total_introns)
colnames(total_introns) <-c("total_intronsREP1Onehour","total_intronsREP1Fourhour","total_intronsREP1Eighthour","total_intronsREP2Onehour","total_intronsREP2Fourhour","total_intronsREP2Eighthour")

TXLabAREP12<-newINSPEcT(tpts,tL,foursu_exons, total_exons,foursu_introns,total_introns,BPPARAM=SerialParam())
TXLabAREP12Express=TXLabAREP12@ratesFirstGuess@assayData$exprs
