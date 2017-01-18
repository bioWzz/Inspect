# rm(list=ls())
# rep1
library(ggplot2)
pre_RNA_total=as.data.frame(as.numeric(Rep2_total_introns)) 
pre_RNA_total=pre_RNA_total[pre_RNA_total>0]
pre_RNA_total=log(pre_RNA_total)

pre_RNA_4SU=as.data.frame(as.numeric(Rep2_foursu_introns))
pre_RNA_4SU=pre_RNA_4SU[pre_RNA_4SU>0]
pre_RNA_4SU=log(pre_RNA_4SU)



# 画pre_RNA_4SU和pre_RNA_total的概率分布图
graphics.off()
data <- data.frame(x = pre_RNA_total)
data1<- data.frame(x = pre_RNA_4SU)
# 生成底层和直方图,概率线的图层
p <- ggplot( data=data,aes(x = x, y = ..density..))
# p <- p + geom_histogram(fill = "navy")
layer1 <- geom_density(data=data,colour = "blue",size = 1)
(p1 <- p + layer1 +xlim(-20, 20))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
layer2 <- geom_density(data=data1,colour = "blue",linetype = 2,size = 1)
(p1 <- p1 + layer2 + xlab("log(Expression)") + ylab("fraction of genes"))
p1+ theme( panel.background = element_blank(),axis.line = element_line(colour = "black"))

p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))




total_RNA_total=as.data.frame(as.numeric(Rep2_foursu_exons))
total_RNA_total=total_RNA_total[total_RNA_total>0]
total_RNA_total=log(total_RNA_total)

total_RNA_4SU=as.data.frame(as.numeric(Rep2_total_exons))
total_RNA_4SU=total_RNA_4SU[total_RNA_4SU>0]
total_RNA_4SU=log(total_RNA_4SU)

# 画total_RNA_4SU和total_RNA_total的概率分布图
graphics.off()
data3 <- data.frame(x = total_RNA_total)
data4<- data.frame(x = total_RNA_4SU)
# 生成底层和直方图,概率线的图层
p <- ggplot(data=data, aes(x = x, y = ..density..))
# p <- p + geom_histogram(fill = "navy")
layer1 <- geom_density(data=data3,colour = "red",size = 1)
(p1 <- p + layer1 +xlim(-20, 20))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
layer2 <- geom_density(data=data4,colour = "red",linetype = 2,size = 1)
(p1 <- p1 + layer2 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))


# 为了坐标统一
graphics.off()
data <- data.frame(x = pre_RNA_total)
data1<- data.frame(x = pre_RNA_4SU)
data3 <- data.frame(x = total_RNA_total)
data4<- data.frame(x = total_RNA_4SU)
# 生成底层和直方图,概率线的图层
p <- ggplot( data=data,aes(x = x, y = ..density..))
# p <- p + geom_histogram(fill = "navy")
layer1 <- geom_density(data=data,colour = "white",size = 1)
(p1 <- p + layer1)
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
layer2 <- geom_density(data=data1,colour = "white",linetype = 2,size = 1)
(p1 <- p1 + layer2 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer3 <- geom_density(data=data3,colour = "red",size = 1)
(p1 <- p1 + layer3 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer4 <- geom_density(data=data4,colour = "red",linetype = 2,size = 1)
(p1 <- p1 + layer4 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))





graphics.off()
data <- data.frame(x = pre_RNA_total)
data1<- data.frame(x = pre_RNA_4SU)
data3 <- data.frame(x = total_RNA_total)
data4<- data.frame(x = total_RNA_4SU)
# 生成底层和直方图,概率线的图层
p <- ggplot( aes(x = x, y = ..density..))
# p <- p + geom_histogram(fill = "navy")
layer1 <- geom_density(data=data,colour = "blue",size = 1)
(p1 <- p + layer1)
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
layer2 <- geom_density(data=data1,colour = "blue",linetype = 2,size = 1)
(p1 <- p1 + layer2 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer3 <- geom_density(data=data3,colour = "white",size = 1)
(p1 <- p1 + layer3 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer4 <- geom_density(data=data4,colour = "white",linetype = 2,size = 1)
(p1 <- p1 + layer4 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))



# rep2
library(ggplot2)
pre_RNA_total=as.data.frame(as.numeric(Rep2_total_introns)) 
pre_RNA_total=pre_RNA_total[pre_RNA_total>0]
pre_RNA_total=log(pre_RNA_total)

pre_RNA_4SU=as.data.frame(as.numeric(Rep2_foursu_introns))
pre_RNA_4SU=pre_RNA_4SU[pre_RNA_4SU>0]
pre_RNA_4SU=log(pre_RNA_4SU)

graphics.off()
data <- data.frame(x = pre_RNA_total)
data1<- data.frame(x = pre_RNA_4SU)
# 生成底层和直方图,概率线的图层
p <- ggplot(data, aes(x = x, y = ..density..))
# p <- p + geom_histogram(fill = "navy")
layer1 <- geom_density(data=data,colour = "blue",size = 1)
(p1 <- p + layer1)
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
layer2 <- geom_density(data=data1,colour = "blue",linetype = 2,size = 1)
(p1 <- p1 + layer2 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))


total_RNA_total=as.data.frame(as.numeric(Rep2_foursu_exons))
total_RNA_total=total_RNA_total[total_RNA_total>0]
total_RNA_total=log(total_RNA_total)

total_RNA_4SU=as.data.frame(as.numeric(Rep2_total_exons))
total_RNA_4SU=total_RNA_4SU[total_RNA_4SU>0]
total_RNA_4SU=log(total_RNA_4SU)

graphics.off()
data3 <- data.frame(x = total_RNA_total)
data4<- data.frame(x = total_RNA_4SU)
# 生成底层和直方图,概率线的图层
p <- ggplot( aes(x = x, y = ..density..))
# p <- p + geom_histogram(fill = "navy")
layer1 <- geom_density(data=data3,colour = "red",size = 1)
(p1 <- p + layer1)
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
layer2 <- geom_density(data=data4,colour = "red",linetype = 2,size = 1)
(p1 <- p1 + layer2 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))





graphics.off()
data <- data.frame(x = pre_RNA_total)
data1<- data.frame(x = pre_RNA_4SU)
data3 <- data.frame(x = total_RNA_total)
data4<- data.frame(x = total_RNA_4SU)
# 生成底层和直方图,概率线的图层
p <- ggplot(aes(x = x, y = ..density..))
# p <- p + geom_histogram(fill = "navy")
layer1 <- geom_density(data=data,colour = "white",size = 1)
(p1 <- p + layer1)
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
layer2 <- geom_density(data=data1,colour = "white",linetype = 2,size = 1)
(p1 <- p1 + layer2 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer3 <- geom_density(data=data3,colour = "red",size = 1)
(p1 <- p1 + layer3 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer4 <- geom_density(data=data4,colour = "red",linetype = 2,size = 1)
(p1 <- p1 + layer4 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))





graphics.off()
data <- data.frame(x = pre_RNA_total)
data1<- data.frame(x = pre_RNA_4SU)
data3 <- data.frame(x = total_RNA_total)
data4<- data.frame(x = total_RNA_4SU)
# 生成底层和直方图,概率线的图层
p <- ggplot(data,aes(x = x, y = ..density..))
# p <- p + geom_histogram(fill = "navy")
layer1 <- geom_density(data=data,colour = "blue",size = 1)
(p1 <- p + layer1)
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
layer2 <- geom_density(data=data1,colour = "blue",linetype = 2,size = 1)
(p1 <- p1 + layer2 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer3 <- geom_density(data=data3,colour = "white",size = 1)
(p1 <- p1 + layer3 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer4 <- geom_density(data=data4,colour = "white",linetype = 2,size = 1)
(p1 <- p1 + layer4 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))


# Combing of rep1 and rep2


library(ggplot2)

pre_RNA_total=as.data.frame(as.numeric(cbind(REP1_total_introns,Rep2_total_introns))) 
pre_RNA_total=pre_RNA_total[pre_RNA_total>0]
pre_RNA_total=log(pre_RNA_total)

pre_RNA_4SU=as.data.frame(as.numeric(cbind(REP1_foursu_introns,Rep2_foursu_introns)))
pre_RNA_4SU=pre_RNA_4SU[pre_RNA_4SU>0]
pre_RNA_4SU=log(pre_RNA_4SU)




graphics.off()
data <- data.frame(x = pre_RNA_total)
data1<- data.frame(x = pre_RNA_4SU)
# 生成底层和直方图,概率线的图层
p <- ggplot(data, aes(x = x, y = ..density..))
# p <- p + geom_histogram(fill = "navy")
layer1 <- geom_density(data=data,colour = "blue",size = 1)
(p1 <- p + layer1)
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
layer2 <- geom_density(data=data1,colour = "blue",linetype = 2,size = 1)
(p1 <- p1 + layer2 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))




total_RNA_total=as.data.frame(as.numeric(cbind(REP1_foursu_exons,Rep2_foursu_exons)))
total_RNA_total=total_RNA_total[total_RNA_total>0]
total_RNA_total=log(total_RNA_total)

total_RNA_4SU=as.data.frame(as.numeric(cbind(REP1_total_exons,Rep2_total_exons)))
total_RNA_4SU=total_RNA_4SU[total_RNA_4SU>0]
total_RNA_4SU=log(total_RNA_4SU)

graphics.off()
data <- data.frame(x = total_RNA_total)
data1<- data.frame(x = total_RNA_4SU)
# 生成底层和直方图,概率线的图层
p <- ggplot(data, aes(x = x, y = ..density..))
# p <- p + geom_histogram(fill = "navy")
layer1 <- geom_density(data=data,colour = "red",size = 1)
(p1 <- p + layer1)
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
layer2 <- geom_density(data=data1,colour = "red",linetype = 2,size = 1)
(p1 <- p1 + layer2 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))




graphics.off()
data <- data.frame(x = pre_RNA_total)
data1<- data.frame(x = pre_RNA_4SU)
data3 <- data.frame(x = total_RNA_total)
data4<- data.frame(x = total_RNA_4SU)
# 生成底层和直方图,概率线的图层
p <- ggplot(data,aes(x = x, y = ..density..))
# p <- p + geom_histogram(fill = "navy")
layer1 <- geom_density(data=data,colour = "white",size = 1)
(p1 <- p + layer1)
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
layer2 <- geom_density(data=data1,colour = "white",linetype = 2,size = 1)
(p1 <- p1 + layer2 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer3 <- geom_density(data=data3,colour = "red",size = 1)
(p1 <- p1 + layer3 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer4 <- geom_density(data=data4,colour = "red",linetype = 2,size = 1)
(p1 <- p1 + layer4 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))





graphics.off()
data <- data.frame(x = pre_RNA_total)
data1<- data.frame(x = pre_RNA_4SU)
data3 <- data.frame(x = total_RNA_total)
data4<- data.frame(x = total_RNA_4SU)
# 生成底层和直方图,概率线的图层
p <- ggplot(data,aes(x = x, y = ..density..))
# p <- p + geom_histogram(fill = "navy")
layer1 <- geom_density(data=data,colour = "blue",size = 1)
(p1 <- p + layer1)
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
layer2 <- geom_density(data=data1,colour = "blue",linetype = 2,size = 1)
(p1 <- p1 + layer2 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer3 <- geom_density(data=data3,colour = "white",size = 1)
(p1 <- p1 + layer3 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer4 <- geom_density(data=data4,colour = "white",linetype = 2,size = 1)
(p1 <- p1 + layer4 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
