rm(list=ls())
library(GenomicFeatures)
library(INSPEcT)
library(Rsamtools)
library(deSolve)
library(BiocParallel)
library(compiler)
library(sqldf)

txdb3 <- makeTxDbFromGFF(file="D:/ivan/Harm/NewTestOne/NewTest/Newgencode.v23.PR8.gtf",format="gtf",taxonomyId=9606)

###################################################################################################################
#####REP1---LabC---------------------------------------------------------------------------------
###Expression of LabC
# LabC-0
paths_4su <- system.file('extdata/MyRep1', 'LabC0hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep1', 'UnlC0hpi.bam', package="INSPEcT")
Rep1LabC0hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total,by ="tx")
# save(Rep1LabC0hpiRPKMsOut, file = paste(ResultDir,"Rep1LabC0hpi.Rdata",sep=""))

# LabC-4
paths_4su <- system.file('extdata/MyRep1', 'LabC4hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep1', 'UnlC4hpi.bam', package="INSPEcT")
Rep1LabC4hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total,by ="tx")
# save(Rep1LabC4hpiRPKMsOut, file = paste(ResultDir,"Rep1LabC4hpi.Rdata",sep=""))

# LabC-8
paths_4su <- system.file('extdata/MyRep1', 'LabC8hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep1', 'UnlC8hpi.bam', package="INSPEcT")
Rep1LabC8hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total,by ="tx")
# save(Rep1LabC8hpiRPKMsOut, file = paste(ResultDir,"Rep1LabC8hpi.Rdata",sep=""))

##LabC rate
LabC0hpifoursu_exons=Rep1LabC0hpiRPKMsOut$rpkms$foursu_exons
LabC4hpifoursu_exons=Rep1LabC4hpiRPKMsOut$rpkms$foursu_exons
LabC8hpifoursu_exons=Rep1LabC8hpiRPKMsOut$rpkms$foursu_exons
REP1_foursu_exons=cbind(LabC0hpifoursu_exons,LabC4hpifoursu_exons,LabC8hpifoursu_exons)
colnames(REP1_foursu_exons) <-c("Onehour","Fourhour","Eighthour")

LabC0hpifoursu_introns=Rep1LabC0hpiRPKMsOut$rpkms$foursu_introns
LabC4hpifoursu_introns=Rep1LabC4hpiRPKMsOut$rpkms$foursu_introns
LabC8hpifoursu_introns=Rep1LabC8hpiRPKMsOut$rpkms$foursu_introns
REP1_foursu_introns=cbind(LabC0hpifoursu_introns,LabC4hpifoursu_introns,LabC8hpifoursu_introns)
colnames(REP1_foursu_introns) <-c("Onehour","Fourhour","Eighthour")

UnlC0hpitotal_exons=Rep1LabC0hpiRPKMsOut$rpkms$total_exons
UnlC4hpitotal_exons=Rep1LabC4hpiRPKMsOut$rpkms$total_exons
UnlC8hpitotal_exons=Rep1LabC8hpiRPKMsOut$rpkms$total_exons
REP1_total_exons=cbind(UnlC0hpitotal_exons,UnlC4hpitotal_exons,UnlC8hpitotal_exons)
colnames(REP1_total_exons) <-c("Onehour","Fourhour","Eighthour")

UnlC0hpitotal_introns=Rep1LabC0hpiRPKMsOut$rpkms$total_introns
UnlC4hpitotal_introns=Rep1LabC4hpiRPKMsOut$rpkms$total_introns
UnlC8hpitotal_introns=Rep1LabC8hpiRPKMsOut$rpkms$total_introns
REP1_total_introns=cbind(UnlC0hpitotal_introns,UnlC4hpitotal_introns,UnlC8hpitotal_introns)
colnames(REP1_total_introns) <-c("Onehour","Fourhour","Eighthour")


tpts <- c(0,4,8)
tL <- 4
LabCREP1 <- newINSPEcT(tpts, tL, REP1_foursu_exons, REP1_total_exons,
                   REP1_foursu_introns, REP1_total_introns, BPPARAM=SerialParam())
LabCREP1Express=LabCREP1@ratesFirstGuess@assayData$exprs


###################################################################################################################
#####REP2---LabC---------------------------------------------------------------------------------------------------

###Expression of LabC
# LabC-0
paths_4su <- system.file('extdata/MyRep2', 'LabC0.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep2', 'ULabC0.bam', package="INSPEcT")
Rep2LabC0hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total,by ="tx")
# save(Rep2LabC0hpiRPKMsOut, file = paste(ResultDir,"Rep2LabC0hpi.Rdata",sep=""))

# LabC-4
paths_4su <- system.file('extdata/MyRep2', 'LabC4.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep2', 'ULabC4.bam', package="INSPEcT")
Rep2LabC4hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total,by ="tx")
# save(Rep2LabC4hpiRPKMsOut, file = paste(ResultDir,"Rep2LabC4hpi.Rdata",sep=""))

# LabC-8
paths_4su <- system.file('extdata/MyRep2', 'LabC8.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep2', 'ULabC8.bam', package="INSPEcT")
Rep2LabC8hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total,by ="tx")
# save(Rep2LabC8hpiRPKMsOut, file = paste(ResultDir,"Rep2LabC8hpi.Rdata",sep=""))

##LabC rate
LabC0hpifoursu_exons=Rep2LabC0hpiRPKMsOut$rpkms$foursu_exons
LabC4hpifoursu_exons=Rep2LabC4hpiRPKMsOut$rpkms$foursu_exons
LabC8hpifoursu_exons=Rep2LabC8hpiRPKMsOut$rpkms$foursu_exons
REP2_foursu_exons=cbind(LabC0hpifoursu_exons,LabC4hpifoursu_exons,LabC8hpifoursu_exons)
colnames(REP2_foursu_exons) <-c("Onehour","Fourhour","Eighthour")

LabC0hpifoursu_introns=Rep2LabC0hpiRPKMsOut$rpkms$foursu_introns
LabC4hpifoursu_introns=Rep2LabC4hpiRPKMsOut$rpkms$foursu_introns
LabC8hpifoursu_introns=Rep2LabC8hpiRPKMsOut$rpkms$foursu_introns
REP2_foursu_introns=cbind(LabC0hpifoursu_introns,LabC4hpifoursu_introns,LabC8hpifoursu_introns)
colnames(REP2_foursu_introns) <-c("Onehour","Fourhour","Eighthour")

UnlC0hpitotal_exons=Rep2LabC0hpiRPKMsOut$rpkms$total_exons
UnlC4hpitotal_exons=Rep2LabC4hpiRPKMsOut$rpkms$total_exons
UnlC8hpitotal_exons=Rep2LabC8hpiRPKMsOut$rpkms$total_exons
REP2_total_exons=cbind(UnlC0hpitotal_exons,UnlC4hpitotal_exons,UnlC8hpitotal_exons)
colnames(REP2_total_exons) <-c("Onehour","Fourhour","Eighthour")

UnlC0hpitotal_introns=Rep2LabC0hpiRPKMsOut$rpkms$total_introns
UnlC4hpitotal_introns=Rep2LabC4hpiRPKMsOut$rpkms$total_introns
UnlC8hpitotal_introns=Rep2LabC8hpiRPKMsOut$rpkms$total_introns
REP2_total_introns=cbind(UnlC0hpitotal_introns,UnlC4hpitotal_introns,UnlC8hpitotal_introns)
colnames(REP2_total_introns) <-c("Onehour","Fourhour","Eighthour")

tpts <- c(0,4,8)
tL <- 4
LabCREP2 <- newINSPEcT(tpts, tL,REP2_foursu_exons,REP2_total_exons,
                       REP2_foursu_introns,REP2_total_introns, BPPARAM=SerialParam())
LabCREP2Express=LabCREP2@ratesFirstGuess@assayData$exprs


#Rep1+Rep2
tpts <- c(0,4,8,0,4,8)
tL <- 4
foursu_exons=cbind(REP1_foursu_exons,REP2_foursu_exons)
colnames(foursu_exons) <-c("foursu_exonsREP1Onehour","foursu_exonsREP1Fourhour","foursu_exonsREP1Eighthour","foursu_exonsREP2Onehour","foursu_exonsREP2Fourhour","foursu_exonsREP2Eighthour")

total_exons=cbind(REP1_total_exons,REP2_total_exons)
colnames(total_exons) <-c("total_exonsREP1Onehour","total_exonsREP1Fourhour","total_exonsREP1Eighthour","total_exonsREP2Onehour","total_exonsREP2Fourhour","total_exonsREP2Eighthour")

foursu_introns=cbind(REP1_foursu_introns, REP2_foursu_introns)
colnames(foursu_introns) <-c("foursu_intronsREP1Onehour","foursu_intronsREP1Fourhour","foursu_intronsREP1Eighthour","foursu_intronsREP2Onehour","foursu_intronsREP2Fourhour","foursu_intronsREP2Eighthour")

total_introns=cbind(REP1_total_introns,REP2_total_introns)
colnames(total_introns) <-c("total_intronsREP1Onehour","total_intronsREP1Fourhour","total_intronsREP1Eighthour","total_intronsREP2Onehour","total_intronsREP2Fourhour","total_intronsREP2Eighthour")

LabCREP12<-newINSPEcT(tpts,tL,foursu_exons, total_exons,foursu_introns,total_introns,BPPARAM=SerialParam())
LabCREP12Express=LabCREP12@ratesFirstGuess@assayData$exprs


