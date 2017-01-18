#inspect RPKM 
# 4su-exon
# inspect
graphics.off()
par(mfrow=c(2,2))

graphics.off()
par(mfrow=c(4,4))
model<-lm(REP2Express[,1]~REP1Express[,1])
plot(REP1Express[,1],REP2Express[,1],xlim=c(0,1500),ylim=c(0,1500))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]


model<-lm(REP24SU_Exon[,1]~REP14SU_Exon[,1])
plot(REP14SU_Exon[,1],REP24SU_Exon[,1],xlim=c(0,max(REP14SU_Exon[,1])), ylim=c(0, max(na.omit(REP14SU_Exon[,1]))))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]

plot(REP14SU_Exon[,2],REP24SU_Exon[,2],xlim=c(0,500),ylim=c(0,500))
plot(REP14SU_Exon[,3],REP24SU_Exon[,3],xlim=c(0,500),ylim=c(0,500))

# count
model<-lm(CountREP24SU_Exon[,1]~CountREP14SU_Exon[,1])
plot(CountREP14SU_Exon[,1],CountREP24SU_Exon[,1],xlim=c(0,max(CountREP14SU_Exon[,1])), ylim=c(0, max(na.omit(CountREP14SU_Exon[,1]))))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]

plot(CountREP14SU_Exon[,2],CountREP24SU_Exon[,2],xlim=c(0,10000),ylim=c(0,10000))
plot(CountREP14SU_Exon[,3],CountREP24SU_Exon[,3],xlim=c(0,10000),ylim=c(0,10000))
k=as.data.frame(CountREP14SU_Exon[,1]/CountREP24SU_Exon[,1])

# edge RPKM
model<-lm(REP2EXONRPKM[1:15000,2]~REP1EXONRPKM[1:15000,2])
plot(REP1EXONRPKM[,2],REP2EXONRPKM[,2])
abline(model,col="blue")
summary(model)
model$coefficients[[2]]




plot(REP1EXONRPKM[,3],REP2EXONRPKM[,3],xlim=c(0,500),ylim=c(0,500))
plot(REP1EXONRPKM[,4],REP2EXONRPKM[,4],xlim=c(0,500),ylim=c(0,500))


# 4su-inron
graphics.off()
par(mfrow=c(1,3))
plot(REP14SU_Intron[,1],REP24SU_Intron[,1],xlim=c(0,20),ylim=c(0,20))
plot(REP14SU_Intron[,2],REP24SU_Intron[,2],xlim=c(0,20),ylim=c(0,20))
plot(REP14SU_Intron[,3],REP24SU_Intron[,3],xlim=c(0,20),ylim=c(0,20))


graphics.off()
par(mfrow=c(1,3))
plot(REP1Total_Exon[,1],REP2Total_Exon[,1],xlim=c(0,500),ylim=c(0,500))
plot(REP1Total_Exon[,2],REP2Total_Exon[,2],xlim=c(0,500),ylim=c(0,500))
plot(REP1Total_Exon[,3],REP2Total_Exon[,3],xlim=c(0,500),ylim=c(0,500))



graphics.off()
par(mfrow=c(1,3))
plot(REP1Total_Intron[,1],REP2Total_Intron[,1],xlim=c(0,20),ylim=c(0,20))
plot(REP1Total_Intron[,2],REP2Total_Intron[,2],xlim=c(0,20),ylim=c(0,20))
plot(REP1Total_Intron[,3],REP2Total_Intron[,3],xlim=c(0,20),ylim=c(0,20))
