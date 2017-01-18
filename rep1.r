rm(list=ls())
library(GenomicFeatures)
library(INSPEcT)
library(Rsamtools)
library(deSolve)
library(BiocParallel)
library(compiler)
library(sqldf)

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
REP1_foursu_introns=cbind(LabA0hpifoursu_introns,LabA4hpifoursu_introns,LabA8hpifoursu_introns)
colnames(REP1_foursu_introns) <-c("Onehour","Fourhour","Eighthour")

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
REP1 <- newINSPEcT(tpts, tL, REP1_foursu_exons, REP1_total_exons,
                   REP1_foursu_introns, REP1_total_introns, BPPARAM=SerialParam())
REP1Express=REP1@ratesFirstGuess@assayData$exprs





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

REP2 <- newINSPEcT(tpts, tL, Rep2_foursu_exons, Rep2_total_exons,
                   Rep2_foursu_introns, Rep2_total_introns, BPPARAM=SerialParam())
REP2Express=REP2@ratesFirstGuess@assayData$exprs




# Rep2Totalcount=Rep2LabA8hpiRPKMsOut$counts$total
# Rep2LabA8hpiRPKMsOut$counts$foursu$exonCounts
REP2Express1=REP2Express
R2=as.table(cbind(rownames(REP2Express1),REP2Express1))
colnames(R2)[1] <- c("R2GeneSymbol")
write.table(R2, file="D:/R2.txt",row.names=FALSE, sep="\t", quote=FALSE)
R2=read.table("D:/R2.txt", header = TRUE, sep = "\t",row.names=NULL)

REP1Express1=REP1Express
R1=as.table(cbind(rownames(REP1Express1),REP1Express1))
colnames(R1)[1] <- c("R1GeneSymbol")
write.table(R1, file="D:/R1.txt",row.names=FALSE, sep="\t", quote=FALSE)
R1=read.table("D:/R1.txt", header = TRUE, sep = "\t",row.names=NULL)

ALL<- sqldf("select * from R1 join R2 where R2.R2GeneSymbol= R1.R1GeneSymbol")
ALL<-na.omit(ALL)
R1=ALL[,1:16]
R2=ALL[,17:32]





#Rep1+Rep2
tpts <- c(0,4,8,0,4,8)
tL <- 4
foursu_exons=cbind(REP1_foursu_exons,Rep2_foursu_exons)
colnames(foursu_exons) <-c("foursu_exonsREP1Onehour","foursu_exonsREP1Fourhour","foursu_exonsREP1Eighthour","foursu_exonsREP2Onehour","foursu_exonsREP2Fourhour","foursu_exonsREP2Eighthour")

total_exons=cbind(REP1_total_exons,Rep2_total_exons)
colnames(total_exons) <-c("total_exonsREP1Onehour","total_exonsREP1Fourhour","total_exonsREP1Eighthour","total_exonsREP2Onehour","total_exonsREP2Fourhour","total_exonsREP2Eighthour")

foursu_introns=cbind(REP1_foursu_introns, Rep2_foursu_introns)
colnames(foursu_introns) <-c("foursu_intronsREP1Onehour","foursu_intronsREP1Fourhour","foursu_intronsREP1Eighthour","foursu_intronsREP2Onehour","foursu_intronsREP2Fourhour","foursu_intronsREP2Eighthour")

total_introns=cbind(REP1_total_introns,Rep2_total_introns)
colnames(total_introns) <-c("total_intronsREP1Onehour","total_intronsREP1Fourhour","total_intronsREP1Eighthour","total_intronsREP2Onehour","total_intronsREP2Fourhour","total_intronsREP2Eighthour")

ALLexons=cbind.fill(foursu_exons,total_exons)
ALLexons[apply(ALLexons==0, FUN = any, 1), ] = NA
ALLexons[apply(ALLexons<2, FUN = any, 1), ] = NA
ALLexons[apply(!is.finite(ALLexons), FUN = any, 1), ] = NA
ALLexons=na.omit(ALLexons)

ALLintron=cbind.fill(foursu_introns,total_introns) 
ALLintron[apply(ALLintron==0, FUN = any, 1), ] = NA
ALLintron[apply(ALLintron<0.3, FUN = any, 1), ] = NA
ALLintron[apply(!is.finite(ALLintron), FUN = any, 1), ] = NA
ALLintron=na.omit(ALLintron)

ALL=cbind.fill(ALLexons,ALLintron)
ALL[apply(ALL==0, FUN = any, 1), ] = NA
ALL[apply(!is.finite(ALL), FUN = any, 1), ] = NA
ALL=na.omit(ALL)

foursu_exons=ALL[,1:6]
total_exons=ALL[,7:12]
foursu_introns=ALL[,13:18]
total_introns=ALL[,19:24]
 



REP12<-newINSPEcT(tpts,tL,foursu_exons, total_exons,foursu_introns,total_introns,BPPARAM=SerialParam())
save(REP12, file = "D:/ivan/Harm/NewTestFour/my/Result/LabA_REP1211_degDuringPulse.Rdata")
modelingParams(REP12)
mycerIdsOneGene <- modelRates(REP12, seed=1, BPPARAM=SerialParam())






REPtest<-REP12[1:900]
object<-REP1[1:1000]

simRates <- makeSimModel(REPtest,50, newTpts=NULL, seed=1)
simData1rep <- makeSimDataset(simRates, tpts, 2, seed=NULL)
simData1rep1 <- modelRates(simData1rep, seed=NULL,BPPARAM = bpparam(),
                           verbose = NULL)

rocCurve(simRates, simData1rep1); title("1 replicate - 3 time points", line=3)


REP12<-newINSPEcT(tpts,tL,foursu_exons, total_exons,foursu_introns,total_introns,BPPARAM=SerialParam())
REP12Express=REP12@ratesFirstGuess@assayData$exprs




















save(REP2, file = "D:/ivan/Harm/NewTestFour/my/Result/LabA_REP2_degDuringPulse.Rdata")
rocCurve(simRates, simData1rep); title("1 replicate - 9 time points", line=3)
rocCurve(simRates, simData3rep); title("3 replicates - 12 time points", line=3)













REPtest<-REP1[1:1700]
object<-REP1[1:1000]

simRates <- makeSimModel(REPtest,1500, newTpts=NULL, seed=1)
simData1rep <- makeSimDataset(simRates, tpts, 1, seed=NULL)
simData1rep1 <- modelRates(simData1rep, seed=1)

rocCurve(simRates, simData1rep1); title("1 replicate - 3 time points", line=3)

simData2rep<-makeSimDataset(simRates,tpts,2,seed=1)
simData2rep2<-modelRates(simData2rep,seed=1)

rocCurve(simRates, simData2rep2); title("2 replicate - 3 time points", line=3)



simData3rep<-makeSimDataset(simRates,tpts,3,seed=1)
simData3rep2<-modelRates(simData3rep,seed=1)

rocCurve(simRates, simData3rep2); title("3 replicate - 3 time points", line=3)

tpts1 <- c(0,4,8,10,11,12,13,14)
simData4rep<-makeSimDataset(simRates,tpts1,2,seed=NULL)
simData4rep1<-modelRates(simData4rep,seed=NULL)
rocCurve(simRates, simData4rep1); title("2 replicate - 8 time points", line=3)















z=REPtest@ratesFirstGuess@assayData$exprs
p=z
p[apply(p<0.01, FUN = any, 1), ] = NA
p[apply(p==0, FUN = any, 1), ] = NA
p[apply(!is.finite(p), FUN = any, 1), ] = NA
p=na.omit(p)




simData2rep<-makeSimDataset(simRates,tpts,2,seed=NULL)
simData2rep<-modelRates(simData2rep,seed=NULL)

simData3rep<-makeSimDataset(simRates,tpts,5,seed=NULL)
simData3rep<-modelRates(simData3rep,seed=NULL)

tpts1 <- c(0,4,8,10,11,12,13,14)
simData4rep<-makeSimDataset(simRates,tpts1,10,seed=NULL)
simData4rep<-modelRates(simData4rep,seed=NULL)


par(mfrow=c(1,2))
rocCurve(simRates,simData1rep)
object=simRates
object2=simData1rep



rocCurve(simRates,simData2rep)
rocCurve(simRates,simData3rep)
rocCurve(simRates,simData4rep)

rocThresholds(simRates,simData1rep,bTsh=c(.01,.01,.05),cTsh=.1)
rocThresholds(simRates,simData2rep,bTsh=c(.01,.01,.05),cTsh=.1)
rocThresholds(simRates,simData3rep,bTsh=c(.01,.01,.05),cTsh=.1)
rocThresholds(simRates,simData4rep,bTsh=c(.01,.01,.05),cTsh=.1)







tpts <- c(0,4,8,0,4,8)
tL <- 4
foursu_exons=cbind(REP1_foursu_exons,Rep2_foursu_exons)
colnames(foursu_exons) <-c("REP1Onehour","REP1Fourhour","REP1Eighthour","REP2Onehour","REP2Fourhour","REP2Eighthour")
total_exons=cbind(REP1_total_exons,Rep2_total_exons)
colnames(total_exons) <-c("REP1Onehour","REP1Fourhour","REP1Eighthour","REP2Onehour","REP2Fourhour","REP2Eighthour")
foursu_introns=cbind(REP1_foursu_introns, Rep2_foursu_introns)
colnames(foursu_introns) <-c("REP1Onehour","REP1Fourhour","REP1Eighthour","REP2Onehour","REP2Fourhour","REP2Eighthour")
total_introns=cbind(REP1_total_introns,Rep2_total_introns)
colnames(total_introns) <-c("REP1Onehour","REP1Fourhour","REP1Eighthour","REP2Onehour","REP2Fourhour","REP2Eighthour")

REP12<-newINSPEcT(tpts,tL,foursu_exons, total_exons,foursu_introns,total_introns,BPPARAM=SerialParam())
save(REP12, file = "D:/ivan/Harm/NewTestFour/my/Result/LabA_REP1211_degDuringPulse.Rdata")

tpts <- c(0,4,8,0,4,8)
tL <- 4




REP1_2<- newINSPEcT(tpts, tL, REP1_foursu_exons, REP1_total_exons,
                    REP1_foursu_introns, REP1_total_introns,  Rep2_foursu_exons, Rep2_total_exons,
                    Rep2_foursu_introns, Rep2_total_introns,BPPARAM=SerialParam())


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
                       REP2_total_exons,REP2_total_introns, BPPARAM=SerialParam())
LabCREP2Express=LabCREP2@ratesFirstGuess@assayData$exprs



