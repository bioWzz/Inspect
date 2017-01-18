CountREP14SU_Exon<-REP1RPKM$counts$foursu$exonCounts
CountREP14SU_Intron<-REP1RPKM$counts$foursu$intronCounts
CountREP1Total_Exon<-REP1RPKM$counts$total$exonCounts
CountREP1Total_Intron<-REP1RPKM$counts$total$intronCounts

CountREP24SU_Exon<-REP2RPKM$counts$foursu$exonCounts
CountREP24SU_Intron<-REP2RPKM$counts$foursu$intronCounts
CountREP2Total_Exon<-REP2RPKM$counts$total$exonCounts
CountREP2Total_Intron<-REP2RPKM$counts$total$intronCounts

EXON=as.data.frame(cbind(rownames(CountREP14SU_Exon),CountREP14SU_Exon,CountREP1Total_Exon,CountREP24SU_Exon,CountREP2Total_Exon))
colnames(EXON)=c("GeneSymbol","r1E0","r1E4","r1E8","r1T0","r1T4","r1T8","r2E0","r2E4","r2E8","r2T0","r2T4","r2T8")
EXON1<- sqldf("select * from EXON join Exonwidths where EXON.GeneSymbol= Exonwidths.GeneSymbol")
rownames(EXON1)=EXON1[,1]
EXON1=EXON1[,c(2:13,15)]
EXON1=as.matrix(EXON1)



  

# calculate the expression of exons using count numbers
# 
graphics.off()
par(mfrow=c(2,3))

r14E0=((as.double(EXON1[,1])/7148281)*9014666*1.079)/as.double(EXON1[,13])
r24E0=((as.double(EXON1[,7])/6042502)*9014666)/as.double(EXON1[,13])
model<-lm(r14E0~r24E0)
plot(r24E0,r14E0,xlim=c(0,max(as.double(r14E0))),ylim=c(0,max(as.double(r14E0))))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]
# EXON1[,1]=r14E0
# EXON1[,7]=r24E0


r14E4=((as.double(EXON1[,2])/7863883)*9014666*1.007)/as.double(EXON1[,13])
r24E4=((as.double(EXON1[,8])/9014666)*9014666)/as.double(EXON1[,13])
model<-lm(r14E4~r24E4)
plot(r24E4,r14E4,xlim=c(0,max(as.double(r14E4))),ylim=c(0,max(as.double(r14E4))))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]
# EXON1[,2]=r14E4
# EXON1[,8]=r24E4


r14E8=((as.double(EXON1[,3])/6361178)*9014666*0.9)/as.double(EXON1[,13])
r24E8=((as.double(EXON1[,9])/7297282)*9014666)/as.double(EXON1[,13])
model<-lm(r14E8~r24E8)
plot(r24E8,r14E8,xlim=c(0,max(as.double(r14E8))),ylim=c(0,max(as.double(r14E8))))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]
# EXON1[,3]=r14E8
# EXON1[,9]=r24E8


r1TE0=((as.double(EXON1[,4])/16433062)*29523912*1.17)/as.double(EXON1[,13])
r2TE0=((as.double(EXON1[,10])/29329621)*29523912)/as.double(EXON1[,13])
model<-lm(r1TE0~r2TE0)
plot(r2TE0,r1TE0,xlim=c(0,max(as.double(r1TE0))),ylim=c(0,max(as.double(r1TE0))))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]
# EXON1[,4]=r1TE0
# EXON1[,10]=r2TE0


r1TE4=((as.double(EXON1[,5])/24313241)*29523912*1.17)/as.double(EXON1[,13])
r2TE4=((as.double(EXON1[,11])/29523912)*29523912)/as.double(EXON1[,13])
model<-lm(r1TE4~r2TE4)
plot(r2TE4,r1TE4,xlim=c(0,max(as.double(r1TE4))),ylim=c(0,max(as.double(r1TE4))))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]
# EXON1[,5]=r1TE4
# EXON1[,11]=r2TE4

r1TE8=((as.double(EXON1[,6])/11760276)*29523912*1.36)/as.double(EXON1[,13])
r2TE8=((as.double(EXON1[,12])/27588260)*29523912)/as.double(EXON1[,13])
model<-lm(r1TE8~r2TE8)
plot(r2TE8,r1TE8,xlim=c(0,max(as.double(r1TE8))),ylim=c(0,max(as.double(r1TE8))))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]
# EXON1[,6]=r1TE8
# EXON1[,12]=r2TE8


INTRON=as.data.frame(cbind(rownames(CountREP14SU_Intron),CountREP14SU_Intron,CountREP1Total_Intron,CountREP24SU_Intron,CountREP2Total_Intron))
colnames(INTRON)=c("GeneSymbol","r1I0","r1I4","r1I8","r1IT0","r1IT4","r1IT8","r2I0","r2I4","r2I8","r2IT0","r2IT4","r2IT8")
INTRON1<-sqldf("select * from INTRON join Intronwidths where INTRON.GeneSymbol= Intronwidths.GeneSymbol")
rownames(INTRON1)=INTRON1[,1]
INTRON1=INTRON1[,c(2:13,15)]
INTRON1=as.matrix(INTRON1)

# calculate the expression of introns using count numbers
# REP1RPKM$counts$foursu$stat

graphics.off()
par(mfrow=c(2,3))

r1I0=((as.double(INTRON1[,1])/7381141)*7381141*1.05)/as.double(INTRON1[,13]) 
r2I0=((as.double(INTRON1[,7])/1753339)*7381141)/as.double(INTRON1[,13])
model<-lm(r1I0~r2I0)
plot(r2I0,r1I0,xlim=c(0,max(as.double(r1I0))),ylim=c(0,max(as.double(r1I0))))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]
INTRON1[,1]=r1I0
INTRON1[,7]=r2I0

r1I4=((as.double(INTRON1[,2])/5542195)*7381141*1.08)/as.double(INTRON1[,13])    
r2I4=((as.double(INTRON1[,8])/6329159)*7381141)/as.double(INTRON1[,13])
model<-lm(r1I4~r2I4)
plot(r2I4,r1I4,xlim=c(0,max(as.double(r1I0))),ylim=c(0,max(as.double(r1I0))))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]
INTRON1[,2]=r1I4
INTRON1[,8]=r2I4

r1I8=((as.double(INTRON1[,3])/5430455)*7381141*1.029)/as.double(INTRON1[,13])
r2I8=((as.double(INTRON1[,9])/6357131)*7381141)/as.double(INTRON1[,13])
model<-lm(r1I8~r2I8)
plot(r2I8,r1I8,xlim=c(0,max(as.double(r1I0))),ylim=c(0,max(as.double(r1I0))))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]
INTRON1[,3]=r1I8
INTRON1[,9]=r2I8

r1TI0=((as.double(INTRON1[,4])/2561666)*7019734*1.08)/as.double(INTRON1[,13])
r2TI0=((as.double(INTRON1[,10])/4720736)*7019734)/as.double(INTRON1[,13])
model<-lm(r1TI0~r2TI0)
plot(r2TI0,r1TI0,xlim=c(0,max(as.double(r1I0))),ylim=c(0,max(as.double(r1I0))))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]
INTRON1[,4]=r1TI0
INTRON1[,10]=r2TI0

r1TI4=((as.double(INTRON1[,5])/4585186)*7019734*1.055)/as.double(INTRON1[,13])
r2TI4=((as.double(INTRON1[,11])/5685720)*7019734)/as.double(INTRON1[,13])
model<-lm(r1TI4~r2TI4)
plot(r2TI4,r1TI4,xlim=c(0,max(as.double(r1I0))),ylim=c(0,max(as.double(r1I0))))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]
INTRON1[,5]=r1TI4
INTRON1[,11]=r2TI4

r1TI8=((as.double(INTRON1[,6])/3906025)*7019734*1.07)/as.double(INTRON1[,13])
r2TI8=((as.double(INTRON1[,12])/7019734)*7019734)/as.double(INTRON1[,13])
model<-lm(r1TI8~r2TI8)
plot(r2TI8,r1TI8,xlim=c(0,max(as.double(r1I0))),ylim=c(0,max(as.double(r1I0))))
abline(model,col="blue")
summary(model)
model$coefficients[[2]]
INTRON1[,6]=r1TI8
INTRON1[,12]=r2TI8


