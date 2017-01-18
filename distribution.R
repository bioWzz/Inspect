e1=REP1_foursu_exons
rm(list=ls())
distributionfigure<-function(e1,e2,color1,color2,linetype1,linetype2,size)
{
  library(ggplot2)
  e1=as.data.frame(as.numeric(unlist(e1))) 
  e1=e1[e1>0]
  e1=log(e1)
  
  e2=as.data.frame(as.numeric(unlist(e2)))
  e2=e2[e2>0]
  e2=log(e2)
  # 画e1和e2的概率分布图
  graphics.off()
  data <- data.frame(x = e1)
  data1<- data.frame(x = e2)
  
  p <- ggplot( aes(x = x, y = ..density..))
  layer1 <- geom_density(data=data,colour = "blue",size = 1)
  (p1 <- p + layer1)
  p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
  layer2 <- geom_density(data=data1,colour = "blue",linetype = 2,size = 1)
  (p1 <- p1 + layer2 + xlab("log(Expression)") + ylab("fraction of junctions"))
  p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
  
  
}
