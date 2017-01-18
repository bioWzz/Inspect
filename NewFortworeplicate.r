rm(list=ls())
library(GenomicFeatures)
library(INSPEcT)
library(Rsamtools)
txdb3 <- makeTxDbFromGFF(file="D:/ivan/Harm/NewTestOne/NewTest/Newgencode.v23.PR8.gtf",format="gtf",taxonomyId=9606)
save(txdb3, file = "D:/ivan/Harm/NewTestFour/my/NewResults/txdb.Rdata")
txdb3=load("D:/ivan/Harm/NewTestFour/my/NewResults/txdb.Rdata")
#### replication one tarnscprition.bam #######
# path=D:ProgramFilesRRR\R-3.3.1\library\INSPEcT\extdata  Rep1\Genome
ResultDir="D:/ivan/Harm/NewTestFour/my/TestForDegrate_2016.6.19/"

###################################################################################################################
#####REP1---LabA---------------------------------------------------------------------------------
###calcutation of Expression of LabA
# LabA-0
paths_4su <- system.file('extdata/MyRep1', 'LabA0hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep1', 'UnlA0hpi.bam', package="INSPEcT")
Rep1LabA0hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)
save(Rep1LabA0hpiRPKMsOut, file = paste(ResultDir,"Rep1LabA0hpi.Rdata",sep=""))


# LabA-4
paths_4su <- system.file('extdata/MyRep1', 'LabA4hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep1', 'UnlA4hpi.bam', package="INSPEcT")
Rep1LabA4hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)
save(Rep1LabA4hpiRPKMsOut, file = paste(ResultDir,"Rep1LabA4hpi.Rdata",sep=""))

# LabA-8
paths_4su <- system.file('extdata/MyRep1', 'LabA8hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep1', 'UnlA8hpi.bam', package="INSPEcT")
Rep1LabA8hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)
save(Rep1LabA8hpiRPKMsOut, file = paste(ResultDir,"Rep1LabA8hpi.Rdata",sep=""))

######################
##LabA rate
LabA0hpifoursu_exons=Rep1LabA0hpiRPKMsOut$rpkms$foursu_exons
LabA4hpifoursu_exons=Rep1LabA4hpiRPKMsOut$rpkms$foursu_exons
LabA8hpifoursu_exons=Rep1LabA8hpiRPKMsOut$rpkms$foursu_exons
foursu_exons=cbind(LabA0hpifoursu_exons,LabA4hpifoursu_exons,LabA8hpifoursu_exons)
colnames(foursu_exons) <-c("Onehour","Fourhour","Eighthour")

LabA0hpifoursu_introns=Rep1LabA0hpiRPKMsOut$rpkms$foursu_introns
LabA4hpifoursu_introns=Rep1LabA4hpiRPKMsOut$rpkms$foursu_introns
LabA8hpifoursu_introns=Rep1LabA8hpiRPKMsOut$rpkms$foursu_introns
foursu_introns=cbind(LabA0hpifoursu_introns,LabA4hpifoursu_introns,LabA8hpifoursu_introns)
colnames(foursu_introns) <-c("Onehour","Fourhour","Eighthour")

UnlA0hpitotal_exons=Rep1LabA0hpiRPKMsOut$rpkms$total_exons
UnlA4hpitotal_exons=Rep1LabA4hpiRPKMsOut$rpkms$total_exons
UnlA8hpitotal_exons=Rep1LabA8hpiRPKMsOut$rpkms$total_exons
total_exons=cbind(UnlA0hpitotal_exons,UnlA4hpitotal_exons,UnlA8hpitotal_exons)
colnames(total_exons) <-c("Onehour","Fourhour","Eighthour")

UnlA0hpitotal_introns=Rep1LabA0hpiRPKMsOut$rpkms$total_introns
UnlA4hpitotal_introns=Rep1LabA4hpiRPKMsOut$rpkms$total_introns
UnlA8hpitotal_introns=Rep1LabA8hpiRPKMsOut$rpkms$total_introns
total_introns=cbind(UnlA0hpitotal_introns,UnlA4hpitotal_introns,UnlA8hpitotal_introns)
colnames(total_introns) <-c("Onehour","Fourhour","Eighthour")


write.table(foursu_exons,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP1_Total_RNA_4sU_LabA.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)
write.table(total_exons,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP1_Total_RNA_total_UnlA.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)
write.table(foursu_introns,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP1_pre_RNA_4sU_LabA.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)
write.table(total_introns,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP1_pre_RNA_total_UnlA.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)



tpts <- c(0,4,8)
tL <- 4
mycerIds <- newINSPEcT(tpts, tL, foursu_exons, total_exons,
                       foursu_introns, total_introns, BPPARAM=SerialParam())
save(mycerIds, file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep1_LabAmycerIds.Rdata")

#Saving LabA the rating
LabASynthesis=(ratesFirstGuess(mycerIds, 'synthesis'))
save(LabASynthesis, file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep1LabASynthesis.Rdata")
write.table(LabASynthesis,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep1LabASynthesis.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

LabADegradation=(ratesFirstGuess(mycerIds, 'degradation'))
save(LabADegradation, file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep1LabADegradation.Rdata")
write.table(LabADegradation,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep1LabADegradation.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

LabAProcessing=(ratesFirstGuess(mycerIds, 'processing'))
save(LabAProcessing, file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep1LabAProcessing.Rdata")
write.table(LabAProcessing,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep1LabAProcessing.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

Rep1A_AllRate=cbind(LabASynthesis,LabADegradation,LabAProcessing)
write.table(Rep1A_AllRate,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep1LabARate.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)



###################################################################################################################
#####REP2---LabA---------------------------------------------------------------------------------
#### replication two tarnscprition.bam#######

ResultDir="D:/ivan/Harm/NewTestFour/my/Results/NewResultsNewResults9_6/"
###calcutation of Expression of LabA
# LabA-0
paths_4su <- system.file('extdata/MyRep2', 'LabA0.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep2', 'ULabA0.bam', package="INSPEcT")
Rep2LabA0hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)
save(Rep2LabA0hpiRPKMsOut, file = paste(ResultDir,"Rep2LabA0hpi.Rdata",sep=""))

# LabA-4
paths_4su <- system.file('extdata/MyRep2', 'LabA4.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep2', 'ULabA4.bam', package="INSPEcT")
Rep2LabA4hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)
save(Rep2LabA4hpiRPKMsOut, file = paste(ResultDir,"Rep2LabA4hpi.Rdata",sep=""))

# LabA-8
paths_4su <- system.file('extdata/MyRep2', 'LabA8.bam', package="INSPEcT")
paths_total <- system.file('extdata/MyRep2', 'ULabA8.bam', package="INSPEcT")
Rep2LabA8hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)
save(Rep2LabA8hpiRPKMsOut, file = paste(ResultDir,"Rep2LabA8hpi.Rdata",sep=""))

######################
##LabA rate
LabA0hpifoursu_exons=Rep2LabA0hpiRPKMsOut$rpkms$foursu_exons
LabA4hpifoursu_exons=Rep2LabA4hpiRPKMsOut$rpkms$foursu_exons
LabA8hpifoursu_exons=Rep2LabA8hpiRPKMsOut$rpkms$foursu_exons
foursu_exons=cbind(LabA0hpifoursu_exons,LabA4hpifoursu_exons,LabA8hpifoursu_exons)
colnames(foursu_exons) <-c("Onehour","Fourhour","Eighthour")

LabA0hpifoursu_introns=Rep2LabA0hpiRPKMsOut$rpkms$foursu_introns
LabA4hpifoursu_introns=Rep2LabA4hpiRPKMsOut$rpkms$foursu_introns
LabA8hpifoursu_introns=Rep2LabA8hpiRPKMsOut$rpkms$foursu_introns
foursu_introns=cbind(LabA0hpifoursu_introns,LabA4hpifoursu_introns,LabA8hpifoursu_introns)
colnames(foursu_introns) <-c("Onehour","Fourhour","Eighthour")

UnlA0hpitotal_exons=Rep2LabA0hpiRPKMsOut$rpkms$total_exons
UnlA4hpitotal_exons=Rep2LabA4hpiRPKMsOut$rpkms$total_exons
UnlA8hpitotal_exons=Rep2LabA8hpiRPKMsOut$rpkms$total_exons
total_exons=cbind(UnlA0hpitotal_exons,UnlA4hpitotal_exons,UnlA8hpitotal_exons)
colnames(total_exons) <-c("Onehour","Fourhour","Eighthour")

UnlA0hpitotal_introns=Rep2LabA0hpiRPKMsOut$rpkms$total_introns
UnlA4hpitotal_introns=Rep2LabA4hpiRPKMsOut$rpkms$total_introns
UnlA8hpitotal_introns=Rep2LabA8hpiRPKMsOut$rpkms$total_introns
total_introns=cbind(UnlA0hpitotal_introns,UnlA4hpitotal_introns,UnlA8hpitotal_introns)
colnames(total_introns) <-c("Onehour","Fourhour","Eighthour")


write.table(foursu_exons,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP2_Total_RNA_4sU_LabA.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)
write.table(total_exons,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP2_Total_RNA_total_UnlA.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)
write.table(foursu_introns,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP2_pre_RNA_4sU_LabA.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)
write.table(total_introns,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP2_pre_RNA_total_UnlA.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

tpts <- c(0,4,8)
tL <- 4
mycerIds <- newINSPEcT(tpts, tL, foursu_exons, total_exons,
                       foursu_introns, total_introns, BPPARAM=SerialParam())
save(mycerIds, file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep2_LabAmycerIds.Rdata")

#Saving LabA the rating
LabASynthesis=(ratesFirstGuess(mycerIds, 'synthesis'))
save(LabASynthesis, file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep2LabASynthesis.Rdata")
write.table(LabASynthesis,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep2LabASynthesis.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

LabADegradation=(ratesFirstGuess(mycerIds, 'degradation'))
save(LabADegradation, file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep2LabADegradation.Rdata")
write.table(LabADegradation,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep2LabADegradation.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

LabAProcessing=(ratesFirstGuess(mycerIds, 'processing'))
save(LabAProcessing, file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep2LabAProcessing.Rdata")
write.table(LabAProcessing,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep2LabAProcessing.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

AllRate=cbind(LabASynthesis,LabADegradation,LabAProcessing)
write.table(AllRate,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults/Rep2LabARate.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)




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
foursu_exons=cbind(LabC0hpifoursu_exons,LabC4hpifoursu_exons,LabC8hpifoursu_exons)
colnames(foursu_exons) <-c("Onehour","Fourhour","Eighthour")

LabC0hpifoursu_introns=Rep1LabC0hpiRPKMsOut$rpkms$foursu_introns
LabC4hpifoursu_introns=Rep1LabC4hpiRPKMsOut$rpkms$foursu_introns
LabC8hpifoursu_introns=Rep1LabC8hpiRPKMsOut$rpkms$foursu_introns
foursu_introns=cbind(LabC0hpifoursu_introns,LabC4hpifoursu_introns,LabC8hpifoursu_introns)
colnames(foursu_introns) <-c("Onehour","Fourhour","Eighthour")

UnlC0hpitotal_exons=Rep1LabC0hpiRPKMsOut$rpkms$total_exons
UnlC4hpitotal_exons=Rep1LabC4hpiRPKMsOut$rpkms$total_exons
UnlC8hpitotal_exons=Rep1LabC8hpiRPKMsOut$rpkms$total_exons
total_exons=cbind(UnlC0hpitotal_exons,UnlC4hpitotal_exons,UnlC8hpitotal_exons)
colnames(total_exons) <-c("Onehour","Fourhour","Eighthour")

UnlC0hpitotal_introns=Rep1LabC0hpiRPKMsOut$rpkms$total_introns
UnlC4hpitotal_introns=Rep1LabC4hpiRPKMsOut$rpkms$total_introns
UnlC8hpitotal_introns=Rep1LabC8hpiRPKMsOut$rpkms$total_introns
total_introns=cbind(UnlC0hpitotal_introns,UnlC4hpitotal_introns,UnlC8hpitotal_introns)
colnames(total_introns) <-c("Onehour","Fourhour","Eighthour")


write.table(foursu_exons,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP1RNA_4sU_LabC.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)
write.table(total_exons,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP1RNA_total_UnlC.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)
write.table(foursu_introns,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP1_pre_RNA_4sU_LabC.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)
write.table(total_introns,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP1_pre_RNA_total_UnlC.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)


tpts <- c(0,4,8)
tL <- 4
mycerIds <- newINSPEcT(tpts, tL, foursu_exons, total_exons,
                       foursu_introns, total_introns, BPPARAM=SerialParam())
save(mycerIds, file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP1_LabCmycerIds.Rdata")

#Saving LabC the rating
LabCSynthesis=(ratesFirstGuess(mycerIds, 'synthesis'))
save(LabCSynthesis, file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep1LabCSynthesis.Rdata")
write.table(LabCSynthesis,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP1_LabCSynthesis.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

LabCDegradation=(ratesFirstGuess(mycerIds, 'degradation'))
save(LabCDegradation, file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep1LabCDegradation.Rdata")
write.table(LabCDegradation,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP1_LabCDegradation.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

LabCProcessing=(ratesFirstGuess(mycerIds, 'processing'))
save(LabCProcessing, file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep1LabCProcessing.Rdata")
write.table(LabCProcessing,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP1_LabCProcessing.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

Rep1C_AllRate=cbind(LabCSynthesis,LabCDegradation,LabCProcessing)
write.table(Rep1C_AllRate,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep1LabCRate.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)



###################################################################################################################
#####REP2---LabC---------------------------------------------------------------------------------------------------

###Expression of LabC
# LabC-0
ResultDir="D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/"
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
foursu_exons=cbind(LabC0hpifoursu_exons,LabC4hpifoursu_exons,LabC8hpifoursu_exons)
colnames(foursu_exons) <-c("Onehour","Fourhour","Eighthour")

LabC0hpifoursu_introns=Rep2LabC0hpiRPKMsOut$rpkms$foursu_introns
LabC4hpifoursu_introns=Rep2LabC4hpiRPKMsOut$rpkms$foursu_introns
LabC8hpifoursu_introns=Rep2LabC8hpiRPKMsOut$rpkms$foursu_introns
foursu_introns=cbind(LabC0hpifoursu_introns,LabC4hpifoursu_introns,LabC8hpifoursu_introns)
colnames(foursu_introns) <-c("Onehour","Fourhour","Eighthour")

UnlC0hpitotal_exons=Rep2LabC0hpiRPKMsOut$rpkms$total_exons
UnlC4hpitotal_exons=Rep2LabC4hpiRPKMsOut$rpkms$total_exons
UnlC8hpitotal_exons=Rep2LabC8hpiRPKMsOut$rpkms$total_exons
total_exons=cbind(UnlC0hpitotal_exons,UnlC4hpitotal_exons,UnlC8hpitotal_exons)
colnames(total_exons) <-c("Onehour","Fourhour","Eighthour")

UnlC0hpitotal_introns=Rep2LabC0hpiRPKMsOut$rpkms$total_introns
UnlC4hpitotal_introns=Rep2LabC4hpiRPKMsOut$rpkms$total_introns
UnlC8hpitotal_introns=Rep2LabC8hpiRPKMsOut$rpkms$total_introns
total_introns=cbind(UnlC0hpitotal_introns,UnlC4hpitotal_introns,UnlC8hpitotal_introns)
colnames(total_introns) <-c("Onehour","Fourhour","Eighthour")


write.table(foursu_exons,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP2RNA_4sU_LabC.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)
write.table(total_exons,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP2RNA_total_UnlC.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)
write.table(foursu_introns,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP2_pre_RNA_4sU_LabC.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)
write.table(total_introns,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP2_pre_RNA_total_UnlC.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

###Rate of LabC
tpts <- c(0,4,8)
tL <- 4
mycerIds <- newINSPEcT(tpts, tL, foursu_exons, total_exons,
                       foursu_introns, total_introns, BPPARAM=SerialParam())
save(mycerIds, file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP2_LabCmycerIds.Rdata")

#Saving LabC the rating
LabCSynthesis=(ratesFirstGuess(mycerIds, 'synthesis'))
save(LabCSynthesis, file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep2LabCSynthesis.Rdata")
write.table(LabCSynthesis,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP2_LabCSynthesis.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

LabCDegradation=(ratesFirstGuess(mycerIds, 'degradation'))
save(LabCDegradation, file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep2LabCDegradation.Rdata")
write.table(LabCDegradation,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP2_LabCDegradation.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

LabCProcessing=(ratesFirstGuess(mycerIds, 'processing'))
save(LabCProcessing, file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep2LabCProcessing.Rdata")
write.table(LabCProcessing,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/REP2_LabCProcessing.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

Rep2C_AllRate=cbind(LabCSynthesis,LabCDegradation,LabCProcessing)
write.table(Rep2C_AllRate,file = "D:/ivan/Harm/NewTestFour/my/Results/NewResults9_6/Rep2LabCRate.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)






