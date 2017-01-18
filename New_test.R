rm(list=ls())
library(GenomicFeatures)
library(INSPEcT)
library(Rsamtools)
ResultDir="D:/ivan/Harm/NewTestFour/my/Result/"

txdb3 <- makeTxDbFromGFF(file="D:/ivan/Harm/NewTestOne/NewTest/Newgencode.v23.PR8.gtf",format="gtf",taxonomyId=9606)
save(txdb3, file = paste(ResultDir,"txdb.Rdata",sep=""))
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
colnames(foursu_exons) <-c("Onehour","Fourhour","Eighthour")

LabA0hpifoursu_introns=Rep1LabA0hpiRPKMsOut$rpkms$foursu_introns
LabA4hpifoursu_introns=Rep1LabA4hpiRPKMsOut$rpkms$foursu_introns
LabA8hpifoursu_introns=Rep1LabA8hpiRPKMsOut$rpkms$foursu_introns
REP1_foursu_introns=cbind(LabA0hpifoursu_introns,LabA4hpifoursu_introns,LabA8hpifoursu_introns)
colnames(foursu_introns) <-c("Onehour","Fourhour","Eighthour","Onehour","Fourhour","Eighthour")

UnlA0hpitotal_exons=Rep1LabA0hpiRPKMsOut$rpkms$total_exons
UnlA4hpitotal_exons=Rep1LabA4hpiRPKMsOut$rpkms$total_exons
UnlA8hpitotal_exons=Rep1LabA8hpiRPKMsOut$rpkms$total_exons
REP1_total_exons=cbind(UnlA0hpitotal_exons,UnlA4hpitotal_exons,UnlA8hpitotal_exons)
colnames(total_exons) <-c("Onehour","Fourhour","Eighthour")

UnlA0hpitotal_introns=Rep1LabA0hpiRPKMsOut$rpkms$total_introns
UnlA4hpitotal_introns=Rep1LabA4hpiRPKMsOut$rpkms$total_introns
UnlA8hpitotal_introns=Rep1LabA8hpiRPKMsOut$rpkms$total_introns
REP1_total_introns=cbind(UnlA0hpitotal_introns,UnlA4hpitotal_introns,UnlA8hpitotal_introns)
colnames(total_introns) <-c("Onehour","Fourhour","Eighthour")


###################################################################################################################
#####REP2---LabA---------------------------------------------------------------------------------
ResultDir="D:/ivan/Harm/NewTestFour/my/Result/"
###calcutation of Expression of LabA
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
colnames(foursu_exons) <-c("Onehour","Fourhour","Eighthour")

LabA0hpifoursu_introns=Rep2LabA0hpiRPKMsOut$rpkms$foursu_introns
LabA4hpifoursu_introns=Rep2LabA4hpiRPKMsOut$rpkms$foursu_introns
LabA8hpifoursu_introns=Rep2LabA8hpiRPKMsOut$rpkms$foursu_introns
Rep2_foursu_introns=cbind(LabA0hpifoursu_introns,LabA4hpifoursu_introns,LabA8hpifoursu_introns)
colnames(foursu_introns) <-c("Onehour","Fourhour","Eighthour")

UnlA0hpitotal_exons=Rep2LabA0hpiRPKMsOut$rpkms$total_exons
UnlA4hpitotal_exons=Rep2LabA4hpiRPKMsOut$rpkms$total_exons
UnlA8hpitotal_exons=Rep2LabA8hpiRPKMsOut$rpkms$total_exons
Rep2_total_exons=cbind(UnlA0hpitotal_exons,UnlA4hpitotal_exons,UnlA8hpitotal_exons)
colnames(total_exons) <-c("Onehour","Fourhour","Eighthour")

UnlA0hpitotal_introns=Rep2LabA0hpiRPKMsOut$rpkms$total_introns
UnlA4hpitotal_introns=Rep2LabA4hpiRPKMsOut$rpkms$total_introns
UnlA8hpitotal_introns=Rep2LabA8hpiRPKMsOut$rpkms$total_introns
Rep2_total_introns=cbind(UnlA0hpitotal_introns,UnlA4hpitotal_introns,UnlA8hpitotal_introns)
colnames(total_introns) <-c("Onehour","Fourhour","Eighthour")


tpts <- c(0,4,8)
tL <- 4
REP1 <- newINSPEcT(tpts, tL, REP1_foursu_exons, REP1_total_exons,
                   REP1_foursu_introns, REP1_total_introns, BPPARAM=SerialParam(),degDuringPulse=TRUE)

REP2 <- newINSPEcT(tpts, tL, REP2_foursu_exons, REP2_total_exons,
                   REP2_foursu_introns, REP2_total_introns, BPPARAM=SerialParam(),degDuringPulse=TRUE)

tpts <- c(0,0,4,4,8,8)
tL <- 4

REP12<-newINSPEcT(tpts, tL, REP1_foursu_exons,REP2_foursu_exons, REP1_total_exons,REP2_total_exons,
                  REP1_foursu_introns, REP2_foursu_introns,REP1_total_introns,REP2_total_introns,BPPARAM=SerialParam(),degDuringPulse=TRUE)

tpts <- c(0,4,8,0,4,8)
tL <- 4

REP1_2<- newINSPEcT(tpts, tL, REP1_foursu_exons, REP1_total_exons,
                    REP1_foursu_introns, REP1_total_introns,  REP2_foursu_exons, REP2_total_exons,
                    REP2_foursu_introns, REP2_total_introns,BPPARAM=SerialParam(),degDuringPulse=TRUE)


#NO degDuringPulse
tpts <- c(0,4,8)
tL <- 4
REP1 <- newINSPEcT(tpts, tL, REP1_foursu_exons, REP1_total_exons,
                   REP1_foursu_introns, REP1_total_introns, BPPARAM=SerialParam())

REP2 <- newINSPEcT(tpts, tL, REP2_foursu_exons, REP2_total_exons,
                   REP2_foursu_introns, REP2_total_introns, BPPARAM=SerialParam())

tpts <- c(0,0,4,4,8,8)
tL <- 4

REP12<-newINSPEcT(tpts, tL, REP1_foursu_exons,REP2_foursu_exons, REP1_total_exons,REP2_total_exons,
                  REP1_foursu_introns, REP2_foursu_introns,REP1_total_introns,REP2_total_introns,BPPARAM=SerialParam())

tpts <- c(0,4,8,0,4,8)
tL <- 4

REP1_2<- newINSPEcT(tpts, tL, REP1_foursu_exons, REP1_total_exons,
                    REP1_foursu_introns, REP1_total_introns,  REP2_foursu_exons, REP2_total_exons,
                    REP2_foursu_introns, REP2_total_introns,BPPARAM=SerialParam())

#mycerIds11<-mycerIds11[1:10]
#geneClass(mycerIds11)
#modelingParams(mycerIds11)$useSigmoidFun<-FALSE
#modelingParams(mycerIds11)$nInit<-20
#modelingParams(mycerIds11)$nIter<-1000
#head(ratesFirstGuess(mycerIds11,'synthesis'))
#mycerIds12<-mycerIds11[1:3]
#mycerIds13<-modelRates(mycerIds12,seed=1,BPPARAM=SerialParam())
#head(ratesFirstGuess(mycerIds13,'synthesis'))
save(mycerIds, file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep1_LabAmycerIds.Rdata")









###################################################################################################################
#####REP1---LabC---------------------------------------------------------------------------------
###Expression of LabC
# LabC-0
paths_4su <- system.file('extdata/MyRep1', 'LabC0hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep1', 'UnlC0hpi.bam', package="INSPEcT")
Rep1LabC0hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)
save(Rep1LabC0hpiRPKMsOut, file = paste(ResultDir,"Rep1LabC0hpi.Rdata",sep=""))

# LabC-4
paths_4su <- system.file('extdata/MyRep1', 'LabC4hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep1', 'UnlC4hpi.bam', package="INSPEcT")
Rep1LabC4hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)
save(Rep1LabC4hpiRPKMsOut, file = paste(ResultDir,"Rep1LabC4hpi.Rdata",sep=""))

# LabC-8
paths_4su <- system.file('extdata/MyRep1', 'LabC8hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep1', 'UnlC8hpi.bam', package="INSPEcT")
Rep1LabC8hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)
save(Rep1LabC8hpiRPKMsOut, file = paste(ResultDir,"Rep1LabC8hpi.Rdata",sep=""))

##LabC rate
LabC0hpifoursu_exons=Rep1LabC0hpiRPKMsOut$rpkms$foursu_exons
LabC4hpifoursu_exons=Rep1LabC4hpiRPKMsOut$rpkms$foursu_exons
LabC8hpifoursu_exons=Rep1LabC8hpiRPKMsOut$rpkms$foursu_exons
REP1_foursu_exons=cbind(LabC0hpifoursu_exons,LabC4hpifoursu_exons,LabC8hpifoursu_exons)
colnames(foursu_exons) <-c("Onehour","Fourhour","Eighthour")

LabC0hpifoursu_introns=Rep1LabC0hpiRPKMsOut$rpkms$foursu_introns
LabC4hpifoursu_introns=Rep1LabC4hpiRPKMsOut$rpkms$foursu_introns
LabC8hpifoursu_introns=Rep1LabC8hpiRPKMsOut$rpkms$foursu_introns
REP1_foursu_introns=cbind(LabC0hpifoursu_introns,LabC4hpifoursu_introns,LabC8hpifoursu_introns)
colnames(foursu_introns) <-c("Onehour","Fourhour","Eighthour")

UnlC0hpitotal_exons=Rep1LabC0hpiRPKMsOut$rpkms$total_exons
UnlC4hpitotal_exons=Rep1LabC4hpiRPKMsOut$rpkms$total_exons
UnlC8hpitotal_exons=Rep1LabC8hpiRPKMsOut$rpkms$total_exons
REP1_total_exons=cbind(UnlC0hpitotal_exons,UnlC4hpitotal_exons,UnlC8hpitotal_exons)
colnames(total_exons) <-c("Onehour","Fourhour","Eighthour")

UnlC0hpitotal_introns=Rep1LabC0hpiRPKMsOut$rpkms$total_introns
UnlC4hpitotal_introns=Rep1LabC4hpiRPKMsOut$rpkms$total_introns
UnlC8hpitotal_introns=Rep1LabC8hpiRPKMsOut$rpkms$total_introns
REP1_total_introns=cbind(UnlC0hpitotal_introns,UnlC4hpitotal_introns,UnlC8hpitotal_introns)
colnames(total_introns) <-c("Onehour","Fourhour","Eighthour")


###################################################################################################################
#####REP2---LabC---------------------------------------------------------------------------------------------------

###Expression of LabC
# LabC-0
paths_4su <- system.file('extdata/MyRep2', 'LabC0.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep2', 'ULabC0.bam', package="INSPEcT")
Rep2LabC0hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)
save(Rep2LabC0hpiRPKMsOut, file = paste(ResultDir,"Rep2LabC0hpi.Rdata",sep=""))

# LabC-4
paths_4su <- system.file('extdata/MyRep2', 'LabC4.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep2', 'ULabC4.bam', package="INSPEcT")
Rep2LabC4hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)
save(Rep2LabC4hpiRPKMsOut, file = paste(ResultDir,"Rep2LabC4hpi.Rdata",sep=""))

# LabC-8
paths_4su <- system.file('extdata/MyRep2', 'LabC8.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep2', 'ULabC8.bam', package="INSPEcT")
Rep2LabC8hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)
save(Rep2LabC8hpiRPKMsOut, file = paste(ResultDir,"Rep2LabC8hpi.Rdata",sep=""))

##LabC rate
LabC0hpifoursu_exons=Rep2LabC0hpiRPKMsOut$rpkms$foursu_exons
LabC4hpifoursu_exons=Rep2LabC4hpiRPKMsOut$rpkms$foursu_exons
LabC8hpifoursu_exons=Rep2LabC8hpiRPKMsOut$rpkms$foursu_exons
REP2_foursu_exons=cbind(LabC0hpifoursu_exons,LabC4hpifoursu_exons,LabC8hpifoursu_exons)
colnames(foursu_exons) <-c("Onehour","Fourhour","Eighthour")

LabC0hpifoursu_introns=Rep2LabC0hpiRPKMsOut$rpkms$foursu_introns
LabC4hpifoursu_introns=Rep2LabC4hpiRPKMsOut$rpkms$foursu_introns
LabC8hpifoursu_introns=Rep2LabC8hpiRPKMsOut$rpkms$foursu_introns
REP2_foursu_introns=cbind(LabC0hpifoursu_introns,LabC4hpifoursu_introns,LabC8hpifoursu_introns)
colnames(foursu_introns) <-c("Onehour","Fourhour","Eighthour")

UnlC0hpitotal_exons=Rep2LabC0hpiRPKMsOut$rpkms$total_exons
UnlC4hpitotal_exons=Rep2LabC4hpiRPKMsOut$rpkms$total_exons
UnlC8hpitotal_exons=Rep2LabC8hpiRPKMsOut$rpkms$total_exons
REP2_total_exons=cbind(UnlC0hpitotal_exons,UnlC4hpitotal_exons,UnlC8hpitotal_exons)
colnames(total_exons) <-c("Onehour","Fourhour","Eighthour")

UnlC0hpitotal_introns=Rep2LabC0hpiRPKMsOut$rpkms$total_introns
UnlC4hpitotal_introns=Rep2LabC4hpiRPKMsOut$rpkms$total_introns
UnlC8hpitotal_introns=Rep2LabC8hpiRPKMsOut$rpkms$total_introns
REP2_total_introns=cbind(UnlC0hpitotal_introns,UnlC4hpitotal_introns,UnlC8hpitotal_introns)
colnames(total_introns) <-c("Onehour","Fourhour","Eighthour")


#Saving LabC the rating
LabCSynthesis=(ratesFirstGuess(mycerIds, 'synthesis'))
save(LabCSynthesis,  file = paste(ResultDir,"Rep2LabCSynthesis.Rdata",sep=""))
write.table(LabCSynthesis,file = paste(ResultDir,"REP2_LabCSynthesis.txt",sep=""),append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

LabCDegradation=(ratesFirstGuess(mycerIds, 'degradation'))
save(LabCDegradation, file = paste(ResultDir,"Rep2LabCDegradation.Rdata",sep=""))
write.table(LabCDegradation,file = paste(ResultDir,"REP2_LabCDegradation.txt",sep=""),append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

LabCProcessing=(ratesFirstGuess(mycerIds, 'processing'))
save(LabCProcessing, file = paste(ResultDir,"Rep2LabCProcessing.Rdata",sep=""))
write.table(LabCProcessing,file = paste(ResultDir,"REP2_LabCProcessing.txt",sep=""),append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

Rep2C_AllRate=cbind(LabCSynthesis,LabCDegradation,LabCProcessing)
write.table(Rep2C_AllRate,file = paste(ResultDir,"Rep2LabCRate.txt",sep=""),append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)






