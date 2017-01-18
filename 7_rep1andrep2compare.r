# 比较labA 样本间的重复性

load("D:/ivan/ivan6-cell/Results/Result_getfromthePROGRAMME/LabAandCresults/labAandC.RData")
library(ggplot2)
relatedfigure(GeneLabCREP1,GeneLabCREP2,1)
relatedfigure(GeneLabCREP1,GeneLabCREP12,1)
relatedfigure(GeneLabCREP2,GeneLabCREP12,1)

relatedfigure(GeneLabAREP1,GeneLabAREP2,1)
relatedfigure(GeneLabAREP1,GeneLabAREP12,1)
relatedfigure(GeneLabAREP2,GeneLabAREP12,1)

relatedfigure<-function(rep1,rep2,i){
 
 rep1=cbind(rownames(rep1),rep1)
 colnames(rep1)[1]="geneid"
 
 rep2=cbind(rownames(rep2),rep2)
 colnames(rep2)[1]="geneid"
 
 t0rep1=rep1[,c(1,2)]
 t4rep1=rep1[,c(1,3)]
 t8rep1=rep1[,c(1,4)]
 
 p0rep1=rep1[,c(1,5)]
 p4rep1=rep1[,c(1,6)]
 p8rep1=rep1[,c(1,7)]
 
 t0rep2=rep2[,c(1,2)]
 t4rep2=rep2[,c(1,3)]
 t8rep2=rep2[,c(1,4)]
 
 p0rep2=rep2[,c(1,5)]
 p4rep2=rep2[,c(1,6)]
 p8rep2=rep2[,c(1,7)]
 
t0=merge(t0rep1,t0rep2, by ="geneid")
t4=merge(t4rep1,t4rep2, by ="geneid")
t8=merge(t8rep1,t8rep2, by ="geneid")
p0=merge(p0rep1,p0rep2, by ="geneid")
p4=merge(p4rep1,p4rep2, by ="geneid")
p8=merge(p8rep1,p8rep2, by ="geneid")
 
 graphics.off()
 par(mfrow=c(2,3))
 
 model<-lm(t0[,2]~t0[,3])
 summary(model)
 model$coefficients[[2]]
 plot(t0[,3],t0[,2],xlim=c(0,1500),ylim=c(0,1500),main=as.character( model$coefficients[[2]]))
 abline(model,col="blue")

 model<-lm(t4[,2]~t4[,3])
 summary(model)
 model$coefficients[[2]]
 plot(t4[,3],t4[,2],xlim=c(0,1500),ylim=c(0,1500),main=as.character( model$coefficients[[2]]))
 abline(model,col="blue")
 
 model<-lm(t8[,2]~t8[,3])
 summary(model)
 model$coefficients[[2]]
 plot(t8[,3],t8[,2],xlim=c(0,1500),ylim=c(0,1500),main=as.character( model$coefficients[[2]]))
 abline(model,col="blue")
 
 
 model<-lm(p0[,2]~p0[,3])
 summary(model)
 model$coefficients[[2]]
 plot(p0[,3],p0[,2],xlim=c(0,50),ylim=c(0,50),main=as.character( model$coefficients[[2]]))
 abline(model,col="blue")
 
 model<-lm(p4[,2]~p4[,3])
 summary(model)
 model$coefficients[[2]]
 plot(p4[,3],p4[,2],xlim=c(0,50),ylim=c(0,50),main=as.character( model$coefficients[[2]]))
 abline(model,col="blue")
 
 model<-lm(p8[,2]~p8[,3])
 summary(model)
 model$coefficients[[2]]
 plot(p8[,3],p8[,2],xlim=c(0,50),ylim=c(0,50),main=as.character( model$coefficients[[2]]))
 abline(model,col="blue")
 
}


threerelatedfigure(GeneLabCREP1,GeneLabCREP2,1)
threerelatedfigure(GeneLabCREP1,GeneLabCREP12,1)
threerelatedfigure(GeneLabCREP2,GeneLabCREP12,1)

threerelatedfigure(GeneLabAREP1,GeneLabAREP2,1)
threerelatedfigure(GeneLabAREP1,GeneLabAREP12,1)
threerelatedfigure(GeneLabAREP2,GeneLabAREP12,1)

threerelatedfigure<-function(rep1,rep2,i){
  
  rep1=cbind(rownames(rep1),rep1)
  colnames(rep1)[1]="geneid"
  
  rep2=cbind(rownames(rep2),rep2)
  colnames(rep2)[1]="geneid"
  
  s0rep1=rep1[,c(1,8)]
  s4rep1=rep1[,c(1,9)]
  s8rep1=rep1[,c(1,10)]
  
  p0rep1=rep1[,c(1,11)]
  p4rep1=rep1[,c(1,12)]
  p8rep1=rep1[,c(1,13)]
  
  d0rep1=rep1[,c(1,14)]
  d4rep1=rep1[,c(1,15)]
  d8rep1=rep1[,c(1,16)]
  
  s0rep2=rep2[,c(1,8)]
  s4rep2=rep2[,c(1,9)]
  s8rep2=rep2[,c(1,10)]
  
  p0rep2=rep2[,c(1,11)]
  p4rep2=rep2[,c(1,12)]
  p8rep2=rep2[,c(1,13)]
  
  d0rep2=rep2[,c(1,14)]
  d4rep2=rep2[,c(1,15)]
  d8rep2=rep2[,c(1,16)]
  
  
  
  
  
  s0=merge(s0rep1,s0rep2, by ="geneid")
  s4=merge(s4rep1,s0rep2, by ="geneid")
  s8=merge(s8rep1,s0rep2, by ="geneid")
  p0=merge(p0rep1,p0rep2, by ="geneid")
  p4=merge(p4rep1,p4rep2, by ="geneid")
  p8=merge(p8rep1,p8rep2, by ="geneid")
  d0=merge(d0rep1,d0rep2, by ="geneid")
  d4=merge(d4rep1,d4rep2, by ="geneid")
  d8=merge(d8rep1,d8rep2, by ="geneid")
  
  graphics.off()
  par(mfrow=c(3,3))
  
  model<-lm(s0[,2]~s0[,3])
  summary(model)
  model$coefficients[[2]]
  plot(s0[,3],s0[,2],xlim=c(0,1000),ylim=c(0,1000),main=as.character( model$coefficients[[2]]))
  abline(model,col="blue")
  
  model<-lm(s4[,2]~s4[,3])
  summary(model)
  model$coefficients[[2]]
  plot(s4[,3],s4[,2],xlim=c(0,1000),ylim=c(0,1000),main=as.character( model$coefficients[[2]]))
  abline(model,col="blue")
  
  model<-lm(s8[,2]~s8[,3])
  summary(model)
  model$coefficients[[2]]
  plot(s8[,3],s8[,2],xlim=c(0,1000),ylim=c(0,1000),main=as.character( model$coefficients[[2]]))
  abline(model,col="blue")
  
  
  model<-lm(p0[,2]~p0[,3])
  summary(model)
  model$coefficients[[2]]
  plot(p0[,3],p0[,2],xlim=c(0,1000),ylim=c(0,1000),main=as.character( model$coefficients[[2]]))
  abline(model,col="blue")
  
  model<-lm(p4[,2]~p4[,3])
  summary(model)
  model$coefficients[[2]]
  plot(p4[,3],p4[,2],xlim=c(0,1000),ylim=c(0,1000),main=as.character( model$coefficients[[2]]))
  abline(model,col="blue")
  
  model<-lm(p8[,2]~p8[,3])
  summary(model)
  model$coefficients[[2]]
  plot(p8[,3],p8[,2],xlim=c(0,1000),ylim=c(0,1000),main=as.character( model$coefficients[[2]]))
  abline(model,col="blue")
  
  model<-lm(d0[,2]~d0[,3])
  summary(model)
  model$coefficients[[2]]
  plot(d0[,3],d0[,2],xlim=c(0,1000),ylim=c(0,1000),main=as.character( model$coefficients[[2]]))
  abline(model,col="blue")
  
  model<-lm(d4[,2]~d4[,3])
  summary(model)
  model$coefficients[[2]]
  plot(d4[,3],d4[,2],xlim=c(0,1000),ylim=c(0,1000),main=as.character( model$coefficients[[2]]))
  abline(model,col="blue")
  
  model<-lm(d8[,2]~d8[,3])
  summary(model)
  model$coefficients[[2]]
  plot(d8[,3],d8[,2],xlim=c(0,1000),ylim=c(0,1000),main=as.character( model$coefficients[[2]]))
  abline(model,col="blue")
  
}


