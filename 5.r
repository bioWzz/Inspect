
EXON=ALLRPKM[,1:12]
INTRON=ALLRPKM[,13:24]
REP1_foursu_exons=EXON[,1:3]
# colnames(REP1_foursu_exons) <-c("Onehour","Fourhour","Eighthour")

REP1_foursu_introns=INTRON[,1:3]
# colnames(REP1_foursu_introns) <-c("Onehour","Fourhour","Eighthour")

REP1_total_exons=EXON[,4:6]
# colnames(REP1_total_exons) <-c("Onehour","Fourhour","Eighthour")

REP1_total_introns=INTRON[,4:6]

# colnames(REP1_total_introns) <-c("Onehour","Fourhour","Eighthour")


tpts <- c(0,4,8)
tL <- 4
REP1 <- newINSPEcT(tpts, tL, REP1_foursu_exons, REP1_total_exons,
                   REP1_foursu_introns, REP1_total_introns, BPPARAM=SerialParam())
REP1Express=REP1@ratesFirstGuess@assayData$exprs







Rep2_foursu_exons=EXON[,7:9]
# colnames(Rep2_foursu_exons) <-c("Onehour","Fourhour","Eighthour")
Rep2_foursu_introns=INTRON[,7:9]
# colnames(Rep2_foursu_introns) <-c("Onehour","Fourhour","Eighthour")

Rep2_total_exons=EXON[,10:12]
# colnames(Rep2_total_exons) <-c("Onehour","Fourhour","Eighthour")

Rep2_total_introns=INTRON[,10:12]
# colnames(Rep2_total_introns) <-c("Onehour","Fourhour","Eighthour")


REP2 <- newINSPEcT(tpts, tL, Rep2_foursu_exons, Rep2_total_exons,
                   Rep2_foursu_introns, Rep2_total_introns, BPPARAM=SerialParam())
REP2Express=REP2@ratesFirstGuess@assayData$exprs


graphics.off()
for(i in 1:15){
png(file=paste("D:/ALL",i,".png"))
par(mfrow=c(1,1))
model<-lm(REP2Express[,i]~REP1Express[,i])
plot(REP1Express[,i],REP2Express[,i],xlab=colnames(REP1Express)[i],ylab=colnames(REP2Express)[i],xlim=c(0,max(na.omit(REP1Express[,i]))), ylim=c(0, max(na.omit(REP1Express[,i]))))
abline(model,col="blue")

dev.off()
}

graphics.off()
par(mfrow=c(4,4))
model<-lm(REP2Express[,1]~REP1Express[,1])
plot(REP1Express[,1],REP2Express[,1],xlim=c(0,1500),ylim=c(0,1500))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]

model<-lm(REP2Express[,2]~REP1Express[,2])
plot(REP1Express[,2],REP2Express[,2])
abline(model,col="blue")
summary(model)
model$coefficients[[2]]

model<-lm(REP2Express[,3]~REP1Express[,3])
plot(REP1Express[,3],REP2Express[,3])
abline(model,col="blue")
summary(model)
model$coefficients[[2]]

model<-lm(REP2Express[,4]~REP1Express[,4])
plot(REP1Express[,4],REP2Express[,4])
abline(model,col="blue")
summary(model)
model$coefficients[[2]]

model<-lm(REP2Express[,5]~REP1Express[,5])
plot(REP1Express[,5],REP2Express[,5])
abline(model,col="blue")
summary(model)
model$coefficients[[2]]

model<-lm(REP2Express[,6]~REP1Express[,6])
plot(REP1Express[,6],REP2Express[,6])
abline(model,col="blue")
summary(model)
model$coefficients[[2]]

model<-lm(REP2Express[,7]~REP1Express[,7])
plot(REP1Express[,7],REP2Express[,7])
abline(model,col="blue")
summary(model)
model$coefficients[[2]]

model<-lm(REP2Express[,8]~REP1Express[,8])
plot(REP1Express[,8],REP2Express[,8],xlim=c(0,350),ylim=c(0,350))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]

model<-lm(REP2Express[,9]~REP1Express[,9])
plot(REP1Express[,9],REP2Express[,9],xlim=c(0,250),ylim=c(0,250))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]


model<-lm(REP2Express[,10]~REP1Express[,10])
plot(REP1Express[,10],REP2Express[,10],xlim=c(0,10),ylim=c(0,10))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]

model<-lm(REP2Express[,11]~REP1Express[,11])
plot(REP1Express[,11],REP2Express[,11],xlim=c(0,10),ylim=c(0,10))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]

model<-lm(REP2Express[,12]~REP1Express[,12])
plot(REP1Express[,12],REP2Express[,12],xlim=c(0,6),ylim=c(0,6))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]

model<-lm(REP2Express[,13]~REP1Express[,13])
plot(REP1Express[,13],REP2Express[,13],xlim=c(0,6),ylim=c(0,6))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]

model<-lm(REP2Express[,14]~REP1Express[,14])
plot(REP1Express[,14],REP2Express[,14],xlim=c(0,6),ylim=c(0,6))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]

model<-lm(REP2Express[,15]~REP1Express[,15])
plot(REP1Express[,15],REP2Express[,15],xlim=c(0,6),ylim=c(0,6))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]
