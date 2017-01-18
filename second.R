#txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
## Getting the results of expression and rate 

library(GenomicFeatures)
library(INSPEcT)
library(Rsamtools)
txdb3 <- makeTxDbFromGFF(file="D:/ivan/Harm/NewTestOne/NewTest/gencode.v23.PR8.gtf",format="gtf",taxonomyId=9606)
save(txdb3, file = "D:/ivan/Harm/txdb.Rdata")
#txdb3 <- makeTxDbFromGFF(file="/sc/orga/projects/sumo/harm/RNA_labing/hg19/gencode.v23.PR8.gtf",format="gtf",taxonomyId=9606)

paths_4su <- system.file('extdata', 'LabA0hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata', 'UnlA0hpi.bam', package="INSPEcT")
LabA0hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)
save(LabA0hpiRPKMsOut, file = "D:/ivan/Harm/LabA0hpi.Rdata")

paths_4su <- system.file('extdata', 'LabA4hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata', 'UnlA4hpi.bam', package="INSPEcT")
LabA4hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)
save(LabA4hpiRPKMsOut, file = "D:/ivan/Harm/LabA4hpi.Rdata")

paths_4su <- system.file('extdata', 'LabA8hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata', 'UnlA8hpi.bam', package="INSPEcT")
LabA8hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)
save(LabA8hpiRPKMsOut, file = "D:/ivan/Harm/LabA8hpi.Rdata")

paths_4su <- system.file('extdata', 'LabC0hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata', 'UnlC0hpi.bam', package="INSPEcT")
LabC0hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)
save(LabC0hpiRPKMsOut, file = "D:/ivan/Harm/LabC0hpi.Rdata")

paths_4su <- system.file('extdata', 'LabC4hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata', 'UnlC4hpi.bam', package="INSPEcT")
LabC4hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)
save(LabC4hpiRPKMsOut, file = "D:/ivan/Harm/LabC4hpi.Rdata")

paths_4su <- system.file('extdata', 'LabC8hpi.bam', package="INSPEcT")
paths_total <- system.file('extdata', 'UnlC8hpi.bam', package="INSPEcT")
LabC8hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)
save(LabC8hpiRPKMsOut, file = "D:/ivan/Harm/LabC8hpi.Rdata")
Rowname=rownames(LabA0hpifoursu_exons)

LabA0hpifoursu_exons=LabA0hpiRPKMsOut$rpkms$foursu_exons
LabA4hpifoursu_exons=LabA4hpiRPKMsOut$rpkms$foursu_exons
LabA8hpifoursu_exons=LabA8hpiRPKMsOut$rpkms$foursu_exons
foursu_exons=cbind(LabA0hpifoursu_exons,LabA4hpifoursu_exons,LabA8hpifoursu_exons)
colnames(foursu_exons) <-c("Onehour","Fourhour","Eighthour")

LabA0hpifoursu_introns=LabA0hpiRPKMsOut$rpkms$foursu_introns
LabA4hpifoursu_introns=LabA4hpiRPKMsOut$rpkms$foursu_introns
LabA8hpifoursu_introns=LabA8hpiRPKMsOut$rpkms$foursu_introns
foursu_introns=cbind(LabA0hpifoursu_introns,LabA4hpifoursu_introns,LabA8hpifoursu_introns)
colnames(foursu_introns) <-c("Onehour","Fourhour","Eighthour")

UnlA0hpitotal_exons=LabA0hpiRPKMsOut$rpkms$total_exons
UnlA4hpitotal_exons=LabA4hpiRPKMsOut$rpkms$total_exons
UnlA8hpitotal_exons=LabA8hpiRPKMsOut$rpkms$total_exons
total_exons=cbind(UnlA0hpifoursu_exons,UnlA4hpifoursu_exons,UnlA8hpifoursu_exons)
colnames(total_exons) <-c("Onehour","Fourhour","Eighthour")

UnlA0hpitotal_introns=LabA0hpiRPKMsOut$rpkms$total_introns
UnlA4hpitotal_introns=LabA4hpiRPKMsOut$rpkms$total_introns
UnlA8hpitotal_introns=LabA8hpiRPKMsOut$rpkms$total_introns
total_introns=cbind(UnlA0hpitotal_introns,UnlA4hpitotal_introns,UnlA8hpitotal_introns)
colnames(total_introns) <-c("Onehour","Fourhour","Eighthour")

write.table(foursu_exons,file = "D:/ivan/Harm/RNA_4sU_LabA.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)
write.table(total_exons,file = "D:/ivan/Harm/RNA_total_UnlA.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

tpts <- c(0,4,8)
tL <- 4
mycerIds <- newINSPEcT(tpts, tL, foursu_exons, total_exons,
                       foursu_introns, total_introns, BPPARAM=SerialParam())
save(mycerIds, file = "D:/ivan/Harm/LabAmycerIds.Rdata")
mycerIds10 <- mycerIds
inHeatmap(mycerIds10, clustering=FALSE)
mycerIds10 <- modelRates(mycerIds10, seed=1)

write.table(foursu_exons,file = "D:/ivan/Harm/RNA_4sU_LabC.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)
write.table(total_exons,file = "D:/ivan/Harm/RNA_total_UnlC.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)




#Saving the rating
LabASynthesis=(ratesFirstGuess(mycerIds, 'synthesis'))
save(LabASynthesis, file = "D:/ivan/Harm/LabASynthesis.Rdata")
write.table(LabASynthesis,file = "D:/ivan/Harm/LabASynthesis.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

LabADegradation=(ratesFirstGuess(mycerIds, 'degradation'))
save(LabADegradation, file = "D:/ivan/Harm/LabADegradation.Rdata")
write.table(LabADegradation,file = "D:/ivan/Harm/LabADegradation.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

LabAProcessing=(ratesFirstGuess(mycerIds, 'processing'))
save(LabAProcessing, file = "D:/ivan/Harm/LabAProcessing.Rdata")
write.table(LabAProcessing,file = "D:/ivan/Harm/LabAProcessing.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

mycerIds10 <- modelRates(mycerIds10, seed=1)
head(viewModelRates(mycerIds10, 'synthesis'))
inHeatmap(mycerIds10, type='model', clustering=FALSE)
geneClass(mycerIds10)
plotGene(mycerIds10, 2, fix.yaxis=TRUE)
chisq <- chisqmodel(mycerIds10)
hist(log10(chisq), main='', xlab='log10 chi-squared p-value')
discard <- which(chisq>.05)
featureNames(mycerIds10)[discard]
## [1] "56372" "74194"
mycerIds10new <- mycerIds10[-discard]





load("D:/ivan/Harm/NewTestFour/my/LabC0hpi.Rdata")
load("D:/ivan/Harm/NewTestFour/my/LabC4hpi.Rdata")
load("D:/ivan/Harm/NewTestFour/my/LabC8hpi.Rdata")

LabC0hpifoursu_exons=LabC0hpiRPKMsOut$rpkms$foursu_exons
LabC4hpifoursu_exons=LabC4hpiRPKMsOut$rpkms$foursu_exons
LabC8hpifoursu_exons=LabC8hpiRPKMsOut$rpkms$foursu_exons
LabCfoursu_exons=cbind(LabC0hpifoursu_exons,LabC4hpifoursu_exons,LabC8hpifoursu_exons)
colnames(LabCfoursu_exons) <-c("Onehour","Fourhour","Eighthour")

LabC0hpifoursu_introns=LabC0hpiRPKMsOut$rpkms$foursu_introns
LabC4hpifoursu_introns=LabC4hpiRPKMsOut$rpkms$foursu_introns
LabC8hpifoursu_introns=LabC8hpiRPKMsOut$rpkms$foursu_introns
LabCfoursu_introns=cbind(LabC0hpifoursu_introns,LabC4hpifoursu_introns,LabC8hpifoursu_introns)
colnames(LabCfoursu_introns) <-c("Onehour","Fourhour","Eighthour")

UnlC0hpitotal_exons=LabC0hpiRPKMsOut$rpkms$total_exons
UnlC4hpitotal_exons=LabC4hpiRPKMsOut$rpkms$total_exons
UnlC8hpitotal_exons=LabC8hpiRPKMsOut$rpkms$total_exons
LabCtotal_exons=cbind(UnlC0hpitotal_exons,UnlC4hpitotal_exons,UnlC8hpitotal_exons)
colnames(LabCtotal_exons) <-c("Onehour","Fourhour","Eighthour")

UnlC0hpitotal_introns=LabC0hpiRPKMsOut$rpkms$total_introns
UnlC4hpitotal_introns=LabC4hpiRPKMsOut$rpkms$total_introns
UnlC8hpitotal_introns=LabC8hpiRPKMsOut$rpkms$total_introns
LabCtotal_introns=cbind(UnlC0hpitotal_introns,UnlC4hpitotal_introns,UnlC8hpitotal_introns)
colnames(LabCtotal_introns) <-c("Onehour","Fourhour","Eighthour")


write.table(LabCfoursu_exons,file = "D:/ivan/Harm/NewTestFour/my/RNA_4sU_LabC.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)
write.table(LabCtotal_exons,file = "D:/ivan/Harm/NewTestFour/my/RNA_total_UnlC.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)






LabCmycerIds <- newINSPEcT(tpts, tL, LabCfoursu_exons, LabCtotal_exons,
                           LabCfoursu_introns, LabCtotal_introns, BPPARAM=SerialParam())
save(LabCmycerIds, file = "D:/ivan/Harm/LabCmycerIds.Rdata")
load("D:/ivan/Harm/NewTestFour/my/LabCmycerIds.Rdata")

LabCSynthesis=(ratesFirstGuess(LabCmycerIds, 'synthesis'))
save(LabCSynthesis, file = "D:/ivan/Harm/NewTestFour/my/LabCSynthesis.Rdata")
write.table(LabCSynthesis,file = "D:/ivan/Harm/NewTestFour/my/LabCSynthesis.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

LabCDegradation=(ratesFirstGuess(LabCmycerIds, 'degradation'))
save(LabCDegradation, file = "D:/ivan/Harm/NewTestFour/my/LabCDegradation.Rdata")
write.table(LabCDegradation,file = "D:/ivan/Harm/NewTestFour/my/LabCDegradation.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)

LabCProcessing=(ratesFirstGuess(LabCmycerIds, 'processing'))
save(LabCProcessing, file = "D:/ivan/Harm/NewTestFour/my/LabCProcessing.Rdata")
write.table(LabCProcessing,file = "D:/ivan/Harm/NewTestFour/my/LabCProcessing.txt",append = FALSE, quote = F,row.names = TRUE,col.names = TRUE)


