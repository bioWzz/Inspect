#Change the ENSG id to GeneSymbol

rm(list=ls())
library(sqldf)
path="D:/ivan/Harm/NewTestFour/my/Results/CompareOneRepAndTwoRep2/CompareOneAndTwo/"

GeneSymbol_LabARate=read.table(paste(path,"GeneSymbol_LabARate.txt",sep=""), header = TRUE, sep = "\t",row.names=NULL)
GeneSymbol_LabCRate=read.table(paste(path,"GeneSymbol_LabCRate.txt",sep=""), header = TRUE, sep = "\t",row.names=NULL)
LabARp1_2=read.table(paste(path,"LabARp1_2.txt",sep=""), header = TRUE, sep = "\t",row.names=NULL)
LabCRp1_2=read.table(paste(path,"LabCRp1_2.txt",sep=""), header = TRUE, sep = "\t",row.names=NULL)

LabA <- sqldf("select * from GeneSymbol_LabARate join LabARp1_2 where GeneSymbol_LabARate.Symbol=LabARp1_2.GeneSymbol")
LabC <- sqldf("select * from GeneSymbol_LabCRate join LabCRp1_2 where GeneSymbol_LabCRate.Symbol=LabCRp1_2.GeneSymbol")

write.table(LabA, file="D:/ivan/Harm/NewTestFour/my/Results/CompareOneRepAndTwoRep2/CompareOneAndTwo/LabA.txt",row.names=FALSE, sep="\t", quote=FALSE)
write.table(LabC, file="D:/ivan/Harm/NewTestFour/my/Results/CompareOneRepAndTwoRep2/CompareOneAndTwo/LabC.txt",row.names=FALSE, sep="\t", quote=FALSE)

colnames(LabA)
par(mfrow=c(3,3))
with(LabA,{plot(SA0,Rp1_sA0)})
with(LabA,{cor.test(SA0,Rp1_sA0)})

with(LabA,{plot(SA0,Rp2_sA0)})
with(LabA,{cor.test(SA0,Rp2_sA0)})


with(LabA,{plot(SA4,Rp1_sA4)})
with(LabA,{cor.test(SA0,Rp2_sA0)})

with(LabA,{plot(SA4,Rp2_sA4)})
with(LabA,{cor.test(SA4,Rp2_sA4)})

with(LabA,{plot(SA8,Rp1_sA8)})
with(LabA,{cor.test(SA8,Rp1_sA8)})

with(LabA,{plot(SA8,Rp2_sA8)})
with(LabA,{cor.test(SA8,Rp2_sA8)})


with(LabA,{plot(Rp1_sA0,Rp2_sA0)})
with(LabA,{cor.test(Rp1_sA0,Rp2_sA0)})

with(LabA,{plot(Rp1_sA4,Rp2_sA4)})
with(LabA,{cor.test(Rp1_sA4,Rp2_sA4)})

with(LabA,{plot(Rp1_sA8,Rp2_sA8)})
with(LabA,{cor.test(Rp1_sA8,Rp2_sA8)})

par(mfrow=c(3,3))
with(LabA,{plot(pA0,Rp1_pA0)})
with(LabA,{plot(pA0,Rp2_pA0)})
with(LabA,{plot(pA4,Rp1_pA4)})
with(LabA,{plot(pA4,Rp2_pA4)})
with(LabA,{plot(pA8,Rp1_pA8)})
with(LabA,{plot(pA8,Rp2_pA8)})

with(LabA,{cor.test(pA0,Rp1_pA0)})
with(LabA,{cor.test(pA0,Rp2_pA0)})
with(LabA,{cor.test(pA4,Rp1_pA4)})
with(LabA,{cor.test(pA4,Rp2_pA4)})
with(LabA,{cor.test(pA8,Rp1_pA8)})
with(LabA,{cor.test(pA8,Rp2_pA8)})

with(LabA,{plot(Rp1_pA0,Rp2_pA0)})
with(LabA,{plot(Rp1_pA4,Rp2_pA4)})
with(LabA,{plot(Rp1_pA8,Rp2_pA8)})

with(LabA,{cor.test(Rp1_pA0,Rp2_pA0)})
with(LabA,{cor.test(Rp1_pA4,Rp2_pA4)})
with(LabA,{cor.test(Rp1_pA8,Rp2_pA8)})


par(mfrow=c(3,3))
with(LabA,{plot(dA0,Rp1_dA0)})
with(LabA,{plot(dA0,Rp2_dA0)})
with(LabA,{plot(dA4,Rp1_dA4)})
with(LabA,{plot(dA4,Rp2_dA4)})
with(LabA,{plot(dA8,Rp1_dA8)})
with(LabA,{plot(dA8,Rp2_dA8)})

with(LabA,{cor.test(dA0,Rp1_dA0)})
with(LabA,{cor.test(dA0,Rp2_dA0)})
with(LabA,{cor.test(dA4,Rp1_dA4)})
with(LabA,{cor.test(dA4,Rp2_dA4)})
with(LabA,{cor.test(dA8,Rp1_dA8)})
with(LabA,{cor.test(dA8,Rp2_dA8)})

with(LabA,{plot(Rp1_dA0,Rp2_dA0)})
with(LabA,{plot(Rp1_dA4,Rp2_dA4)})
with(LabA,{plot(Rp1_dA8,Rp2_dA8)})

with(LabA,{cor.test(Rp1_dA0,Rp2_dA0)})
with(LabA,{cor.test(Rp1_dA4,Rp2_dA4)})
with(LabA,{cor.test(Rp1_dA8,Rp2_dA8)})
