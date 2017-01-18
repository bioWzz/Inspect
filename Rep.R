rm(list=ls())
library(GenomicFeatures)
library(INSPEcT)
library(Rsamtools)
library(deSolve)
library(BiocParallel)
library(compiler)
library(sqldf)
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




exons=cbind(REP1_foursu_exons, REP1_total_exons)
exons[apply(exons==0, FUN = any, 1), ] = NA
exons[apply(exons<0.1, FUN = any, 1), ] = NA
exons[apply(!is.finite(exons), FUN = any, 1), ] = NA
exons=na.omit(exons)
write.table(exons, file="D:/exon.txt",row.names=TRUE, sep="\t", quote=FALSE)
exons=read.table("D:/exon.txt", header = TRUE, sep = "\t",row.names=NULL)


intron=cbind(REP1_foursu_introns, REP1_total_introns)
intron[apply(intron==0, FUN = any, 1), ] = NA
intron[apply(intron<0.1, FUN = any, 1), ] = NA
intron[apply(!is.finite(intron), FUN = any, 1), ] = NA
intron=na.omit(intron)
write.table(intron, file="D:/intron.txt",row.names=TRUE, sep="\t", quote=FALSE)
intron=read.table("D:/intron.txt", header = TRUE, sep = "\t",row.names=NULL)

exons<- sqldf("select * from exons join intron where exons.Gensymbol= intron.Gensymbol")
intron=exons[,8:14]
exons=exons[,1:7]

write.table(exons, file="D:/exon1.txt",row.names=F, sep="\t", quote=FALSE)
exons=read.table("D:/exon1.txt", header = TRUE, sep = "\t",row.names=1)

write.table(intron, file="D:/intron1.txt",row.names=TRUE, sep="\t", quote=FALSE)
intron=read.table("D:/intron1.txt", header = TRUE, sep = "\t",row.names=1)


REP1_foursu_exons=exons[,1:3]
REP1_total_exons=exons[,4:6]

REP1_foursu_introns=intron[,1:3]
REP1_total_introns=intron[,4:6]


tpts <- c(0,4,8)
tL <- 4
REP1 <- newINSPEcT(tpts, tL, REP1_foursu_exons, REP1_total_exons,
                   REP1_foursu_introns, REP1_total_introns, BPPARAM=SerialParam())
save(REP1, file = "D:/ivan/Harm/NewTestFour/my/Result/LabA_REP1_degDuringPulse.Rdata")

z=REP1@ratesFirstGuess@assayData$exprs

# ###################################################################################################################
# #####REP2---LabA---------------------------------------------------------------------------------
# ResultDir="D:/ivan/Harm/NewTestFour/my/Result/"
# ###calcutation of Expression of LabA
# # LabA-0
# paths_4su <- system.file('extdata/MyRep2', 'LabA0.bam', package="INSPEcT")
# paths_total <- system.file('extdata/MyRep2', 'ULabA0.bam', package="INSPEcT")
# Rep2LabA0hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)
# 
# 
# # LabA-4
# paths_4su <- system.file('extdata/MyRep2', 'LabA4.bam', package="INSPEcT")
# paths_total <- system.file('extdata/MyRep2', 'ULabA4.bam', package="INSPEcT")
# Rep2LabA4hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)
# 
# 
# # LabA-8
# paths_4su <- system.file('extdata/MyRep2', 'LabA8.bam', package="INSPEcT")
# paths_total <- system.file('extdata/MyRep2', 'ULabA8.bam', package="INSPEcT")
# Rep2LabA8hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su, paths_total)
# 
# 
# ######################
# ##LabA rate
# LabA0hpifoursu_exons=Rep2LabA0hpiRPKMsOut$rpkms$foursu_exons
# LabA4hpifoursu_exons=Rep2LabA4hpiRPKMsOut$rpkms$foursu_exons
# LabA8hpifoursu_exons=Rep2LabA8hpiRPKMsOut$rpkms$foursu_exons
# Rep2_foursu_exons=cbind(LabA0hpifoursu_exons,LabA4hpifoursu_exons,LabA8hpifoursu_exons)
# colnames(Rep2_foursu_exons) <-c("Onehour","Fourhour","Eighthour")
# 
# LabA0hpifoursu_introns=Rep2LabA0hpiRPKMsOut$rpkms$foursu_introns
# LabA4hpifoursu_introns=Rep2LabA4hpiRPKMsOut$rpkms$foursu_introns
# LabA8hpifoursu_introns=Rep2LabA8hpiRPKMsOut$rpkms$foursu_introns
# Rep2_foursu_introns=cbind(LabA0hpifoursu_introns,LabA4hpifoursu_introns,LabA8hpifoursu_introns)
# colnames(Rep2_foursu_introns) <-c("Onehour","Fourhour","Eighthour")
# 
# UnlA0hpitotal_exons=Rep2LabA0hpiRPKMsOut$rpkms$total_exons
# UnlA4hpitotal_exons=Rep2LabA4hpiRPKMsOut$rpkms$total_exons
# UnlA8hpitotal_exons=Rep2LabA8hpiRPKMsOut$rpkms$total_exons
# Rep2_total_exons=cbind(UnlA0hpitotal_exons,UnlA4hpitotal_exons,UnlA8hpitotal_exons)
# colnames(Rep2_total_exons) <-c("Onehour","Fourhour","Eighthour")
# 
# UnlA0hpitotal_introns=Rep2LabA0hpiRPKMsOut$rpkms$total_introns
# UnlA4hpitotal_introns=Rep2LabA4hpiRPKMsOut$rpkms$total_introns
# UnlA8hpitotal_introns=Rep2LabA8hpiRPKMsOut$rpkms$total_introns
# Rep2_total_introns=cbind(UnlA0hpitotal_introns,UnlA4hpitotal_introns,UnlA8hpitotal_introns)
# colnames(Rep2_total_introns) <-c("Onehour","Fourhour","Eighthour")
# 
# all=cbind(REP1_foursu_exons, REP1_total_exons,
#           REP1_foursu_introns, REP1_total_introns)



#rm(simData1rep)
#rm(simRates)
#rm(REPtest)



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



#Make the makeSimModel by myself
nGenes=1000
newTpts=NULL
probs=c(constant=.5,sigmoid=.3,impulse=.2)
na.rm=TRUE
seed=NULL

# tpts <- object@tpts
# 
# concentrations <- list(
#   total=p[,1:3]
#   , total_var=as.double(total)
#   , preMRNA=p[,4:6]
#   , preMRNA_var=as.double(preMRNA)
# )
# 
# rates <- list(
#    alpha=p[,7:9]
#   , alpha_var=as.double(alpha)
#   , beta=p[,10:12]
#   , gamma=p[,13:15]
# )
tpts <- object@tpts
concentrations <- list(
  total=ratesFirstGuess(object, 'total')
  , total_var=ratesFirstGuessVar(object, 'total')
  , preMRNA=ratesFirstGuess(object, 'preMRNA')
  , preMRNA_var=ratesFirstGuessVar(object, 'preMRNA')
)
rates <- list(
  alpha=ratesFirstGuess(object, 'synthesis')
  , alpha_var=ratesFirstGuessVar(object, 'synthesis')
  , beta=ratesFirstGuess(object, 'degradation')
  , gamma=ratesFirstGuess(object, 'processing')
)

out <-makeSimData(nGenes, tpts,concentrations, rates
                    , newTpts=newTpts, probs=probs, na.rm=FALSE, seed=seed)

simdataSpecs <- out$simdataSpecs
simdataSpecs <- lapply(simdataSpecs, function(x) list(x))
newObject <- new('INSPEcT_model')
newObject@ratesSpecs <- simdataSpecs
newObject@params$sim$flag <- TRUE
newObject@params$sim$foldchange <- out$simulatedFC
newObject@params$sim$noiseVar <- out$noiseVar

simRates<- newObject
#Make the makeSimModel is finished

simData1rep <- makeSimDataset(simRates, tpts, 1, seed=1)
simData1rep<-modelRates(simData1rep,seed=NULL)
#simData1rep<-modelRates(simData1rep,seed=NULL)

#Make the modelRates by myself
object<-simData1rep

ll=object@ratesFirstGuess@assayData$exprs
ppp=ll
ppp[apply(ppp==0, FUN = any, 1), ] = 0.015
ppp[apply(!is.finite(ppp), FUN = any, 1), ] =0.01
ppp=na.omit(ppp)

seed=NULL
# , nCores=NULL
#BPPARAM=bpparam()
verbose=NULL

if( length(object@model@ratesSpecs) > 0 )
  stop('Remove the model before running the model again. (See "?removeModel")')
if( !is.null(seed) && !is.numeric(seed) )
  stop('Seed argument must be either NULL or numeric.')
# if( !is.null(nCores) && !is.numeric(nCores) )
# 	stop('nCores argument must be either NULL or numeric.')
if( !is.null(verbose) && !is.logical(verbose) )
  stop('verbose argument must be either NULL or logical.')

tpts <- object@tpts
log_shift <- .find_tt_par(tpts)

concentrations <- list(
  total=ppp[,1:3]
  , total_var=as.double(total)
  , preMRNA=ppp[,4:6]
  , preMRNA_var=as.double(preMRNA)
)

rates <- list(
  alpha=ppp[,7:9]
  , alpha_var=as.double(alpha)
  , beta=ppp[,10:12]
  , gamma=ppp[,13:15]
)

ratesSpecs <- inspect.engine(tpts, log_shift, concentrations, rates
                              , nInit=object@params$nInit
                              , nIter=object@params$nIter
                              , na.rm=object@params$na.rm
                              , verbose=if(is.null(verbose)) object@params$verbose else verbose
                              , estimateRatesWith=object@params$estimateRatesWith
                              , sigmoidDegradation=object@params$useSigmoidFun
                              , sigmoidSynthesis=object@params$useSigmoidFun
                              , sigmoidTotal=object@params$useSigmoidFun
                              , sigmoidProcessing=object@params$useSigmoidFun
                              , sigmoidPre=object@params$useSigmoidFun
                              , testOnSmooth=object@params$testOnSmooth
                              , seed=seed
)
ratesSpecs<-paramSpecs
#names(ratesSpecs) <- featureNames(object@ratesFirstGuess)
names(ratesSpecs) <- featureNames(object@ratesFirstGuess)
object@model@ratesSpecs <- ratesSpecs
object <- makeModelRates(object)
simData1rep<-object 
# tpts
# log_shift
# concentrations
# rates
# nInit=object@params$nInit
#  nIter=object@params$nIter
#  na.rm=object@params$na.rm
#  verbose=if(is.null(verbose)) object@params$verbose else verbose
#  estimateRatesWith=object@params$estimateRatesWith
#  sigmoidDegradation=object@params$useSigmoidFun
#  sigmoidSynthesis=object@params$useSigmoidFun
#  sigmoidTotal=object@params$useSigmoidFun
#  sigmoidProcessing=object@params$useSigmoidFun
#  sigmoidPre=object@params$useSigmoidFun
#  testOnSmooth=object@params$testOnSmooth
#  seed=seed




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



data('simRates', package='INSPEcT')
data('simData1rep', package='INSPEcT')
data('simData3rep', package='INSPEcT')
par(mfrow=c(1,2))
data('')

REP2 <- newINSPEcT(tpts, tL, Rep2_foursu_exons, Rep2_total_exons,
                   Rep2_foursu_introns, Rep2_total_introns, BPPARAM=SerialParam())
save(REP2, file = "D:/ivan/Harm/NewTestFour/my/Result/LabA_REP2_degDuringPulse.Rdata")
rocCurve(simRates, simData1rep); title("1 replicate - 9 time points", line=3)
rocCurve(simRates, simData3rep); title("3 replicates - 12 time points", line=3)



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






