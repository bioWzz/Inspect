#edge--d diff gene between wt_PR8 dNS1/sictrl_PR8 dNS1
library(edgeR)
library(splines)
rm(list=ls())
path="D:/ivan/ivan2_RNASeq/DiffGenes/TR_C_INF_SiUBE2IAffectdgenes/"
counts=read.table(paste(path,"RPKM1.txt",sep=""), header = TRUE, sep = "", quote = "\t",numerals = c("allow.loss"),row.names=1,stringsAsFactors = default.stringsAsFactors(),nrows =60508)
ls()
class(counts)
dim(counts)
colnames(counts)
counts=counts[,c(2,3,4,5,10,11,12,13,26,27,28,29,34,35,36,37)]
#head(counts[, 1:7], 3)
grp <- as.factor(substr(colnames(counts), 1, 3)) ##substr截取列名的前两个字母
table(grp)
o <- order(grp)
pairs(log2(1+counts[2:60508,o[1:8]]), pch=".",lower.panel=NULL) ##这个图给我们样品间关系的一揽印象
d <- DGEList(counts=counts, group=grp)
d <- calcNormFactors(d)
d$samples
#cps <- cpm(d) ## count per million
#k <- rowSums(cps>=1) > 2
#d <- d[k,]
dim(d)
cols <- as.numeric(d$samples$group)
plotMDS(d,col=cols)
mm <- model.matrix(~-1+grp)
d <- estimateGLMCommonDisp(d,mm)
d <- estimateGLMTrendedDisp(d,mm)
d <- estimateGLMTagwiseDisp(d,mm)
plotBCV(d)
d$common.dispersion
sqrt(d$common.dispersion)
plotMeanVar(d, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE)
par(mfrow=c(1,2))
plotSmear(d, pair=c("UBE","Con"), ylim=c(-5,5))
f <- glmFit(d,mm)
con <- makeContrasts("UBE-Con"=grpsiC-grpWT_,levels=colnames(mm))
con
lrt <- glmLRT(f,contrast=con)
topTags(lrt,20)
cps <- cpm(d)
o <- order(colnames(counts))
tt <- topTags(lrt, n=Inf)$table
write.table(tt, file=paste(path,"edge_INF_WT_SIC.txt",sep=""), sep="\t", quote=FALSE)

#DESeq2--diff gene between wt_PR8 dNS1/sictrl_PR8 dNS1
rm(list=ls())
library(DESeq2)
path="D:/ivan/ivan2_RNASeq/DiffGenes/TR_C_INF_SiUBE2IAffectdgenes/"
counts=read.table(paste(path,"RPKM1.txt",sep=""), header = TRUE, sep = "", quote = "\t",numerals = c("allow.loss"),row.names=1,stringsAsFactors = default.stringsAsFactors(),nrows =60508)
colnames(counts)
counts=counts[,c(2,3,4,5,10,11,12,13,26,27,28,29,34,35,36,37)]
ls()
grp <- as.factor(substr(colnames(counts), 1, 3)) ##substr截取列名的前两个字母
dds <- DESeqDataSetFromMatrix(counts, colData=data.frame(grp), design=formula(~-1+grp)) ##这里的colData必须是一个DataFrame或者data.frame。每一行都对应着counts中的一列。design中的公式和limma中的方法一致。
design(dds)
dds <- DESeq(dds,betaPrior=FALSE)
res <- results(dds) ##得到结果
write.table(res, file=paste(path,"DESeq2_INF_WT_SIC.txt",sep=""), sep="\t", quote=FALSE)

#edge--diff gene between siBUE2I_PR8 dNS1/sictrl_PR8 dNS1
library(edgeR)
library(splines)
rm(list=ls())
path="D:/ivan/ivan2_RNASeq/DiffGenes/TR_C_INF_SiUBE2IAffectdgenes/"
counts=read.table(paste(path,"RPKM1.txt",sep=""), header = TRUE, sep = "", quote = "\t",numerals = c("allow.loss"),row.names=1,stringsAsFactors = default.stringsAsFactors(),nrows =60508)
ls()
class(counts)
dim(counts)
colnames(counts)
counts=counts[,c(2,3,4,5,10,11,12,13,14,15,16,17,22,23,24,25)]
grp <- as.factor(substr(colnames(counts), 1, 3)) ##substr截取列名的前两个字母
table(grp)
o <- order(grp)
pairs(log2(1+counts[2:60508,o[1:8]]), pch=".",lower.panel=NULL) ##这个图给我们样品间关系的一揽印象
d <- DGEList(counts=counts, group=grp)
d <- calcNormFactors(d)
d$samples
#cps <- cpm(d) ## count per million
#k <- rowSums(cps>=1) > 2
#d <- d[k,]
dim(d)
cols <- as.numeric(d$samples$group)
plotMDS(d,col=cols)
mm <- model.matrix(~-1+grp)
d <- estimateGLMCommonDisp(d,mm)
d <- estimateGLMTrendedDisp(d,mm)
d <- estimateGLMTagwiseDisp(d,mm)
plotBCV(d)
d$common.dispersion
sqrt(d$common.dispersion)
plotMeanVar(d, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE)
par(mfrow=c(1,2))
plotSmear(d, pair=c("UBE","Con"), ylim=c(-5,5))
f <- glmFit(d,mm)
con <- makeContrasts("UBE-Con"=grpsiC-grpsiU,levels=colnames(mm))
con
lrt <- glmLRT(f,contrast=con)
topTags(lrt,20)
cps <- cpm(d)
o <- order(colnames(counts))
tt <- topTags(lrt, n=Inf)$table
write.table(tt, file=paste(path,"edge_INF_SI_SIC.txt",sep=""), sep="\t", quote=FALSE)

#DESeq2--diff gene between siBUE2I_PR8 dNS1/sictrl_PR8 dNS1
library(DESeq2)
path="D:/ivan/ivan2_RNASeq/DiffGenes/TR_C_INF_SiUBE2IAffectdgenes/"
counts=read.table(paste(path,"RPKM1.txt",sep=""), header = TRUE, sep = "", quote = "\t",numerals = c("allow.loss"),row.names=1,stringsAsFactors = default.stringsAsFactors(),nrows =60508)
counts=counts[,c(2,3,4,5,10,11,12,13,14,15,16,17,22,23,24,25)]
ls()
grp <- as.factor(substr(colnames(counts), 1, 3)) ##substr截取列名的前两个字母
dds <- DESeqDataSetFromMatrix(counts, colData=data.frame(grp), design=formula(~-1+grp)) ##这里的colData必须是一个DataFrame或者data.frame。每一行都对应着counts中的一列。design中的公式和limma中的方法一致。
design(dds)
dds <- DESeq(dds,betaPrior=FALSE)
res <- results(dds) ##得到结果
write.table(res, file=paste(path,"DESeq2_INF_SI_SIC.txt",sep=""), sep="\t", quote=FALSE)






#Differnet gene of edge between siBUE2I_PR8 dNS1/sictrl_PR8 dNS1 for fold_change=2 or 1.5
path="D:/ivan/ivan2_RNASeq/DiffGenes/TR_C_INF_SiUBE2IAffectdgenes/"
edgeDiffGeneRaw=read.table(paste(path,"edge_INF_SI_SIC.txt",sep=""), header = TRUE, sep = "\t")
colnames( edgeDiffGeneRaw )
attach(edgeDiffGeneRaw)
edgeDiffGeneRaw1<-edgeDiffGeneRaw[(logFC>1|logFC<(-1))&FDR<0.01,]
write.table(edgeDiffGeneRaw1, file=paste(path,"Diffgene_edge_INF_SI_SIC2_foldchange2.txt",sep=""), sep="\t", quote=FALSE)
edgeDiffGeneRaw2<-edgeDiffGeneRaw[(logFC>0.59|logFC<(-0.59))&FDR<0.01,]
write.table(edgeDiffGeneRaw2, file=paste(path,"Diffgene_edge_UI_siBEU2_foldchange1.5.txt",sep=""), sep="\t", quote=FALSE)


#Differnet gene of DESeq2 between siBUE2I_PR8 dNS1/sictrl_PR8 dNS1 for fold_change=2 or 1.5
path="D:/ivan/ivan2_RNASeq/DiffGenes/TR_C_INF_SiUBE2IAffectdgenes/"
DESeq2DiffGeneRaw=read.table(paste(path,"DESeq2_INF_SI_SIC.txt",sep=""), header = TRUE, sep = "\t")
colnames( DESeq2DiffGeneRaw )
attach(DESeq2DiffGeneRaw)
DESeq2DiffGeneRaw1<-DESeq2DiffGeneRaw[(log2FoldChange>1|log2FoldChange<(-1))&padj<0.01&! is.na(padj),]
write.table(DESeq2DiffGeneRaw1, file=paste(path,"Diffgene_DESeq2_INF_SI_SIC_foldchange2.txt",sep=""), sep="\t", quote=FALSE)

DESeq2DiffGeneRaw2<-DESeq2DiffGeneRaw[(log2FoldChange>0.59|log2FoldChange<(-0.59))&padj<0.01&! is.na(padj),]
write.table(DESeq2DiffGeneRaw2, file=paste(path,"Diffgene_DESeq2_UI_siBEU2_sictrl_foldchange1.5.txt",sep=""), sep="\t", quote=FALSE)

#check the same diffgenes of three methods by venn plot
#http://www.cnblogs.com/xianghang123/archive/2013/03/25/2980623.html
#install.packages("gplots")

#foldchang_2
library(gplots)
DESeq2name=read.table(paste(path,"Diffgene_DESeq2_INF_SI_SIC_foldchange2.txt",sep=""), header = TRUE, sep = "\t",row.names=NULL)
DESeq2name=DESeq2name[1]
edgename=read.table(paste(path,"Diffgene_edge_INF_SI_SIC2_foldchange2.txt",sep=""), header = TRUE, sep = "\t",row.names=NULL)
edgename=edgename[1]
input<-list(DESeq2name,edgename)
venn(input)


##DESeq2:produce the fold change values of wt_PR8 dNS1/sictrl_PR8 AND siBUE2I_PR8 dNS1/sictrl_PR8 dNS1
## fold change 2
library(sqldf)
rm(list=ls())
path="D:/ivan/ivan2_RNASeq/DiffGenes/TR_C_INF_SiUBE2IAffectdgenes/"
DESeq2_siUBE_sicontrol=read.table(paste(path,"Diffgene_DESeq2_INF_SI_SIC_foldchange1.5.txt",sep=""), header = TRUE, sep = "\t",row.names=NULL)
colnames( DESeq2_siUBE_sicontrol)
DESeq2_siUBE_sicontrol=DESeq2_siUBE_sicontrol[,c(1,3)]
colnames(DESeq2_siUBE_sicontrol) <- c("gene","log2FoldChange") 

UTDESeq2DiffGeneRaw=read.table(paste(path,"DESeq2_INF_WT_SIC.txt",sep=""), header = TRUE, sep = "\t",row.names=NULL)
colnames( UTDESeq2DiffGeneRaw )
UTDESeq2DiffGeneRaw=UTDESeq2DiffGeneRaw[,c(1,3)]
colnames(UTDESeq2DiffGeneRaw) <- c("gene","log2FoldChange") 
newdf <- sqldf("select DESeq2_siUBE_sicontrol.gene,DESeq2_siUBE_sicontrol.log2FoldChange,UTDESeq2DiffGeneRaw.gene,UTDESeq2DiffGeneRaw.log2FoldChange from DESeq2_siUBE_sicontrol join UTDESeq2DiffGeneRaw where DESeq2_siUBE_sicontrol.gene=UTDESeq2DiffGeneRaw.gene")
write.table(newdf, file=paste(path,"DESeq2_SIUBW2I_sictrl_UT_compare_foldchange2.txt",sep=""),row.names=FALSE, sep="\t", quote=FALSE)

##画热图
library(gplots)
path="D:/Program Files/R/R-3.3.1/library/DEGseq/extdata/C_INF_SiUBE2IAffectdgenes/"
heatmapdata=read.table(paste(path,"DESeq2_SIUBW2I_sictrl_UT_compare_foldchange2.txt",sep=""), header = TRUE, sep = "", quote = "\t",numerals = c("allow.loss"),row.names=1,stringsAsFactors = default.stringsAsFactors(),nrows =430)
heatmapdata=na.omit(heatmapdata) 
heatmapdata<- as.matrix(heatmapdata) 
heatmap.2(heatmapdata)




##edge:produce the fold change values of SIUB/CONTRL AND UT/CONTRL
## fold change 2
rm(list=ls())
library(sqldf)
path="D:/ivan/ivan2_RNASeq/DiffGenes/TR_C_INF_SiUBE2IAffectdgenes/"

edgeDiffGeneRaw1=read.table(paste(path,"Diffgene_edge_INF_SI_SIC2_foldchange2.txt",sep=""), header = TRUE, sep = "\t",row.names=NULL)
colnames(edgeDiffGeneRaw1)
edgeDiffGeneRaw1=edgeDiffGeneRaw1[,c(1,2)]
colnames( edgeDiffGeneRaw1)
colnames(edgeDiffGeneRaw1) <- c("gene","logFC") 


UTedgeDiffGeneRaw=read.table(paste(path,"edge_INF_WT_SIC.txt",sep=""), header = TRUE, sep = "\t",row.names=NULL)
colnames(UTedgeDiffGeneRaw)
UTedgeDiffGeneRaw=UTedgeDiffGeneRaw[,c(1,2)]
colnames(UTedgeDiffGeneRaw )
colnames(UTedgeDiffGeneRaw) <- c("gene","logFC") 
newdf <- sqldf("select edgeDiffGeneRaw1.gene,edgeDiffGeneRaw1.logFC,UTedgeDiffGeneRaw.gene,UTedgeDiffGeneRaw.logFC from edgeDiffGeneRaw1 join UTedgeDiffGeneRaw where edgeDiffGeneRaw1.gene=UTedgeDiffGeneRaw.gene")
write.table(newdf, file=paste(path,"edge_SIUBW2I_sictrl_UT_compare_foldchange2.txt",sep=""), sep="\t", quote=FALSE)

##画热图
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(made4)
path="D:/Program Files/R/R-3.3.1/library/DEGseq/extdata/C_INF_SiUBE2IAffectdgenes/"
heatmapdata=read.table(paste(path,"edge_SIUBW2I_sictrl_UT_compare_foldchange2.txt",sep=""), header = TRUE, sep = "", quote = "\t",numerals = c("allow.loss"),row.names=NULL,stringsAsFactors = default.stringsAsFactors(),nrows =430)
heatmapdata=na.omit(heatmapdata) 
heatmapdata<- as.matrix(heatmapdata) 
heatmap.2(heatmapdata)

heatmapdata=as.data.frame(heatmapdata)
colnames(heatmapdata )
heatmapdata <- heatmapdata[order(heatmapdata$siU_siC),]
row.names(heatmapdata) <- heatmapdata$gene
heatmapdata <- heatmapdata[,2:3]
heatmapdata_matrix <- data.matrix(heatmapdata)
mycolors <- colorRampPalette(c("red","white","darkblue" ))
nba_heatmap <- heatmap(heatmapdata_matrix, Rowv=NA, Colv=NA, col =mycolors(599), scale="none",na.rm = TRUE, ColSideColors,RowSideColors, margins=c(5,10))



