rm(list=ls())
library(GenomicFeatures)
library(INSPEcT)
library(Rsamtools)
library(deSolve)
library(BiocParallel)
library(compiler)
library(sqldf)

# Gene
# 
# ResultDir="D:/ivan/Harm/NewTestFour/my/Result/"

txdb3 <- makeTxDbFromGFF(file="D:/ivan/Harm/NewTestOne/NewTest/Newgencode.v23.PR8.gtf",format="gtf",taxonomyId=9606)
# save(txdb3, file = paste(ResultDir,"txdb.Rdata",sep=""))
# load(paste(ResultDir,"txdb.Rdata",sep=""))
#### replication one tarnscprition.bam #######
# path=D:ProgramFilesRRR\R-3.3.1\library\INSPEcT\extdata  Rep1\Genome
###################################################################################################################
#####REP1---LabC---------------------------------------------------------------------------------
###calcutation of Expression of LabC
# LabC-0
paths_4su <- system.file('extdata/MyRep1', 'LabC0hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep1', 'UnlC0hpi.bam', package="INSPEcT")
Rep1LabC0hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)


# LabC-4
paths_4su <- system.file('extdata/MyRep1', 'LabC4hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep1', 'UnlC4hpi.bam', package="INSPEcT")
Rep1LabC4hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)

# LabC-8
paths_4su <- system.file('extdata/MyRep1', 'LabC8hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep1', 'UnlC8hpi.bam', package="INSPEcT")
Rep1LabC8hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)

######################
##LabC rate
LabC0hpifoursu_exons=Rep1LabC0hpiRPKMsOut$rpkms$foursu_exons
LabC4hpifoursu_exons=Rep1LabC4hpiRPKMsOut$rpkms$foursu_exons
LabC8hpifoursu_exons=Rep1LabC8hpiRPKMsOut$rpkms$foursu_exons

REP1_foursu_exons=cbind(LabC0hpifoursu_exons,LabC4hpifoursu_exons,LabC8hpifoursu_exons)
colnames(REP1_foursu_exons) <-c("Onehour","Fourhour","Eighthour")

LabC0hpifoursu_introns=Rep1LabC0hpiRPKMsOut$rpkms$foursu_introns
LabC4hpifoursu_introns=Rep1LabC4hpiRPKMsOut$rpkms$foursu_introns
LabC8hpifoursu_introns=Rep1LabC8hpiRPKMsOut$rpkms$foursu_introns

REP1_foursu_introns=cbind( LabC0hpifoursu_introns,LabC4hpifoursu_introns,LabC8hpifoursu_introns)
colnames(REP1_foursu_introns) <-c("Onehour","Fourhour","Eighthour")


# UnlA0hpitotal_exons=Rep1LabC0hpiRPKMsOut$rpkms$total_exons
UnlA0hpitotal_exons=Rep1LabC0hpiRPKMsOut$rpkms$total_exons
UnlA4hpitotal_exons=Rep1LabC4hpiRPKMsOut$rpkms$total_exons
UnlA8hpitotal_exons=Rep1LabC8hpiRPKMsOut$rpkms$total_exons
REP1_total_exons=cbind(UnlA0hpitotal_exons,UnlA4hpitotal_exons,UnlA8hpitotal_exons)
colnames(REP1_total_exons) <-c("Onehour","Fourhour","Eighthour")

UnlA0hpitotal_introns=Rep1LabC0hpiRPKMsOut$rpkms$total_introns
UnlA4hpitotal_introns=Rep1LabC4hpiRPKMsOut$rpkms$total_introns
UnlA8hpitotal_introns=Rep1LabC8hpiRPKMsOut$rpkms$total_introns
REP1_total_introns=cbind(UnlA0hpitotal_introns,UnlA4hpitotal_introns,UnlA8hpitotal_introns)
colnames(REP1_total_introns) <-c("Onehour","Fourhour","Eighthour")

tpts <- c(0,4,8)
tL <- 4
GeneLabCREP1 <- newINSPEcT(tpts, tL, REP1_foursu_exons, REP1_total_exons,
                           REP1_foursu_introns, REP1_total_introns, BPPARAM=SerialParam())
GeneLabCREP1Express=GeneLabCREP1@ratesFirstGuess@assayData$exprs





# # ###################################################################################################################
# #####REP2---LabC---------------------------------------------------------------------------------
ResultDir="D:/ivan/Harm/NewTestFour/my/Result/"
# ###calcutation of Expression of LabC
# LabC-0
paths_4su <- system.file('extdata/MyRep2', 'LabC0.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep2', 'ULabC0.bam', package="INSPEcT")
Rep2LabC0hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)


# LabC-4
paths_4su <- system.file('extdata/MyRep2', 'LabC4.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep2', 'ULabC4.bam', package="INSPEcT")
Rep2LabC4hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)


# LabC-8
paths_4su <- system.file('extdata/MyRep2', 'LabC8.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep2', 'ULabC8.bam', package="INSPEcT")
Rep2LabC8hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)


######################
##LabC rate
LabC0hpifoursu_exons=Rep2LabC0hpiRPKMsOut$rpkms$foursu_exons
LabC4hpifoursu_exons=Rep2LabC4hpiRPKMsOut$rpkms$foursu_exons
LabC8hpifoursu_exons=Rep2LabC8hpiRPKMsOut$rpkms$foursu_exons
Rep2_foursu_exons=cbind(LabC0hpifoursu_exons,LabC4hpifoursu_exons,LabC8hpifoursu_exons)
colnames(Rep2_foursu_exons) <-c("Onehour","Fourhour","Eighthour")

LabC0hpifoursu_introns=Rep2LabC0hpiRPKMsOut$rpkms$foursu_introns
LabC4hpifoursu_introns=Rep2LabC4hpiRPKMsOut$rpkms$foursu_introns
LabC8hpifoursu_introns=Rep2LabC8hpiRPKMsOut$rpkms$foursu_introns
Rep2_foursu_introns=cbind(LabC0hpifoursu_introns,LabC4hpifoursu_introns,LabC8hpifoursu_introns)
colnames(Rep2_foursu_introns) <-c("Onehour","Fourhour","Eighthour")

UnlA0hpitotal_exons=Rep2LabC0hpiRPKMsOut$rpkms$total_exons
UnlA4hpitotal_exons=Rep2LabC4hpiRPKMsOut$rpkms$total_exons
UnlA8hpitotal_exons=Rep2LabC8hpiRPKMsOut$rpkms$total_exons
Rep2_total_exons=cbind(UnlA0hpitotal_exons,UnlA4hpitotal_exons,UnlA8hpitotal_exons)
colnames(Rep2_total_exons) <-c("Onehour","Fourhour","Eighthour")

UnlA0hpitotal_introns=Rep2LabC0hpiRPKMsOut$rpkms$total_introns
UnlA4hpitotal_introns=Rep2LabC4hpiRPKMsOut$rpkms$total_introns
UnlA8hpitotal_introns=Rep2LabC8hpiRPKMsOut$rpkms$total_introns
Rep2_total_introns=cbind(UnlA0hpitotal_introns,UnlA4hpitotal_introns,UnlA8hpitotal_introns)
colnames(Rep2_total_introns) <-c("Onehour","Fourhour","Eighthour")

GeneLabCREP2 <- newINSPEcT(tpts, tL, Rep2_foursu_exons, Rep2_total_exons,
                           Rep2_foursu_introns, Rep2_total_introns, BPPARAM=SerialParam())
GeneLabCREP2Express=GeneLabCREP2@ratesFirstGuess@assayData$exprs

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

GeneLabCREP12<-newINSPEcT(tpts,tL,foursu_exons, total_exons,foursu_introns,total_introns,BPPARAM=SerialParam())
GeneLabCREP12Express=GeneLabCREP12@ratesFirstGuess@assayData$exprs

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
#####REP1---LabC---------------------------------------------------------------------------------
###calcutation of Expression of LabC
# LabC-0
paths_4su <- system.file('extdata/MyRep1', 'LabC0hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep1', 'UnlC0hpi.bam', package="INSPEcT")
Rep1LabC0hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total,by ="tx")


# LabC-4
paths_4su <- system.file('extdata/MyRep1', 'LabC4hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep1', 'UnlC4hpi.bam', package="INSPEcT")
Rep1LabC4hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total,by ="tx")

# LabC-8
paths_4su <- system.file('extdata/MyRep1', 'LabC8hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep1', 'UnlC8hpi.bam', package="INSPEcT")
Rep1LabC8hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total,by ="tx")

######################
##LabC rate
LabC0hpifoursu_exons=Rep1LabC0hpiRPKMsOut$rpkms$foursu_exons
LabC4hpifoursu_exons=Rep1LabC4hpiRPKMsOut$rpkms$foursu_exons
LabC8hpifoursu_exons=Rep1LabC8hpiRPKMsOut$rpkms$foursu_exons

REP1_foursu_exons=cbind(LabC0hpifoursu_exons,LabC4hpifoursu_exons,LabC8hpifoursu_exons)
colnames(REP1_foursu_exons) <-c("Onehour","Fourhour","Eighthour")

LabC0hpifoursu_introns=Rep1LabC0hpiRPKMsOut$rpkms$foursu_introns
LabC4hpifoursu_introns=Rep1LabC4hpiRPKMsOut$rpkms$foursu_introns
LabC8hpifoursu_introns=Rep1LabC8hpiRPKMsOut$rpkms$foursu_introns

REP1_foursu_introns=cbind( LabC0hpifoursu_introns,LabC4hpifoursu_introns,LabC8hpifoursu_introns)
colnames(REP1_foursu_introns) <-c("Onehour","Fourhour","Eighthour")


# UnlA0hpitotal_exons=Rep1LabC0hpiRPKMsOut$rpkms$total_exons
UnlA0hpitotal_exons=Rep1LabC0hpiRPKMsOut$rpkms$total_exons
UnlA4hpitotal_exons=Rep1LabC4hpiRPKMsOut$rpkms$total_exons
UnlA8hpitotal_exons=Rep1LabC8hpiRPKMsOut$rpkms$total_exons
REP1_total_exons=cbind(UnlA0hpitotal_exons,UnlA4hpitotal_exons,UnlA8hpitotal_exons)
colnames(REP1_total_exons) <-c("Onehour","Fourhour","Eighthour")

UnlA0hpitotal_introns=Rep1LabC0hpiRPKMsOut$rpkms$total_introns
UnlA4hpitotal_introns=Rep1LabC4hpiRPKMsOut$rpkms$total_introns
UnlA8hpitotal_introns=Rep1LabC8hpiRPKMsOut$rpkms$total_introns
REP1_total_introns=cbind(UnlA0hpitotal_introns,UnlA4hpitotal_introns,UnlA8hpitotal_introns)
colnames(REP1_total_introns) <-c("Onehour","Fourhour","Eighthour")

tpts <- c(0,4,8)
tL <- 4
TXLabCREP1 <- newINSPEcT(tpts, tL, REP1_foursu_exons, REP1_total_exons,
                         REP1_foursu_introns, REP1_total_introns, BPPARAM=SerialParam())
TXLabCREP1Express=TXLabCREP1@ratesFirstGuess@assayData$exprs





# # ###################################################################################################################
# #####REP2---LabC---------------------------------------------------------------------------------
ResultDir="D:/ivan/Harm/NewTestFour/my/Result/"
# ###calcutation of Expression of LabC
# LabC-0
paths_4su <- system.file('extdata/MyRep2', 'LabC0.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep2', 'ULabC0.bam', package="INSPEcT")
Rep2LabC0hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total,by ="tx")


# LabC-4
paths_4su <- system.file('extdata/MyRep2', 'LabC4.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep2', 'ULabC4.bam', package="INSPEcT")
Rep2LabC4hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total,by ="tx")


# LabC-8
paths_4su <- system.file('extdata/MyRep2', 'LabC8.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep2', 'ULabC8.bam', package="INSPEcT")
Rep2LabC8hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total,by ="tx")


######################
##LabC rate
LabC0hpifoursu_exons=Rep2LabC0hpiRPKMsOut$rpkms$foursu_exons
LabC4hpifoursu_exons=Rep2LabC4hpiRPKMsOut$rpkms$foursu_exons
LabC8hpifoursu_exons=Rep2LabC8hpiRPKMsOut$rpkms$foursu_exons
Rep2_foursu_exons=cbind(LabC0hpifoursu_exons,LabC4hpifoursu_exons,LabC8hpifoursu_exons)
colnames(Rep2_foursu_exons) <-c("Onehour","Fourhour","Eighthour")

LabC0hpifoursu_introns=Rep2LabC0hpiRPKMsOut$rpkms$foursu_introns
LabC4hpifoursu_introns=Rep2LabC4hpiRPKMsOut$rpkms$foursu_introns
LabC8hpifoursu_introns=Rep2LabC8hpiRPKMsOut$rpkms$foursu_introns
Rep2_foursu_introns=cbind(LabC0hpifoursu_introns,LabC4hpifoursu_introns,LabC8hpifoursu_introns)
colnames(Rep2_foursu_introns) <-c("Onehour","Fourhour","Eighthour")

UnlA0hpitotal_exons=Rep2LabC0hpiRPKMsOut$rpkms$total_exons
UnlA4hpitotal_exons=Rep2LabC4hpiRPKMsOut$rpkms$total_exons
UnlA8hpitotal_exons=Rep2LabC8hpiRPKMsOut$rpkms$total_exons
Rep2_total_exons=cbind(UnlA0hpitotal_exons,UnlA4hpitotal_exons,UnlA8hpitotal_exons)
colnames(Rep2_total_exons) <-c("Onehour","Fourhour","Eighthour")

UnlA0hpitotal_introns=Rep2LabC0hpiRPKMsOut$rpkms$total_introns
UnlA4hpitotal_introns=Rep2LabC4hpiRPKMsOut$rpkms$total_introns
UnlA8hpitotal_introns=Rep2LabC8hpiRPKMsOut$rpkms$total_introns
Rep2_total_introns=cbind(UnlA0hpitotal_introns,UnlA4hpitotal_introns,UnlA8hpitotal_introns)
colnames(Rep2_total_introns) <-c("Onehour","Fourhour","Eighthour")

TXLabCREP2 <- newINSPEcT(tpts, tL, Rep2_foursu_exons, Rep2_total_exons,
                         Rep2_foursu_introns, Rep2_total_introns, BPPARAM=SerialParam())
TXLabCREP2Express=TXLabCREP2@ratesFirstGuess@assayData$exprs

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

TXLabCREP12<-newINSPEcT(tpts,tL,foursu_exons, total_exons,foursu_introns,total_introns,BPPARAM=SerialParam())
TXLabCREP12Express=TXLabCREP12@ratesFirstGuess@assayData$exprs
