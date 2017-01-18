rm(list=ls())
library(ggplot2)


load("D:/ivan/ivan6-cell/Results/Result_getfromthePROGRAMME/LabAandCresults/LabAGene.RData")
z=distributionofrate(LabAREP1Express,1,0.1,0.01)
z1=distributionofrate(LabAREP2Express,1,0.1,0.01)
z2=distributionofrate(LabAREP12Express,1,0.1,0.01)

load("D:/ivan/ivan6-cell/Results/Result_getfromthePROGRAMME/LabAandCresults/LabCGene.RData")
z=distributionofrate(LabCREP1Express,1,0.1,0.01)
z1=distributionofrate(LabCREP2Express,1,0.1,0.01)
z2=distributionofrate(LabCREP12Express,1,0.1,0.01)



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# 
# 
# 定义生成转录，合成降解的分布函数级中位数线图的函数
distributionofrate<-function(Express,totalThreshold,preRNAThreshold,rateThreshold){
  
  # 数据进行初筛
  Express[apply(Express[,1:3]<=totalThreshold, FUN = any, 1), ] = NA
  Express[apply(Express[,4:6]<=preRNAThreshold, FUN = any, 1), ] = NA
  Express[apply(Express[,7:15]<=rateThreshold, FUN = any, 1), ] = NA
  Express=na.omit(Express)
  print(length(Express[,1]))
  
  # synthsis
  synthsis=Express[,7:9]
  synthsis=as.data.frame(as.numeric(unlist(synthsis)))
  s=median(as.numeric(unlist(synthsis)))
  print(s)
  synthsis=log2(synthsis)
  
  # processing
  processing=Express[,13:15]
  processing=as.data.frame(as.numeric(unlist(processing)))
  p=median(as.numeric(unlist(processing)))
  print(p)
  processing=log2(processing)
  
  # degration
  degration=Express[,10:12]
  degration=as.data.frame(as.numeric(unlist(degration)))
  d=median(as.numeric(unlist(degration)))
  print(d)
  degration=log2(degration)
  
  # #画概率分布图 
  graphics.off()
  par(mfrow=c(1,3))
  
  p1=ggplot(synthsis, aes(x = synthsis),y = ..density..)+
    # geom_histogram(bins = 20,fill = "steelblue", colour = "black")+
    geom_density(data=synthsis,colour = "black",size = 1)+
    geom_vline(xintercept =log2(s) ,colour="red")+
    xlim(-5,5)+
    labs(x = "Gene's transcription rate ", y = "fraction of genes", title = "")+
    theme(panel.grid =element_blank(),panel.background = element_blank())  ## 删去网格线
  
  p2=ggplot(processing, aes(x = processing),y = ..density..)+
    geom_density(data=processing,colour = "black",size = 1)+
    geom_vline(xintercept =log2(p) ,colour="red")+
    xlim(-5,5)+
    labs(x = "Gene's processing rate ", y = "fraction of genes", title = "")+
    theme(panel.grid =element_blank(),panel.background = element_blank())  ## 删去网格线
  
  p3=ggplot(degration, aes(x = degration),y = ..density..)+
    geom_density(data=degration,colour = "black",size = 1)+
    geom_vline(xintercept =log2(d) ,colour="red")+
    xlim(-5,5)+
    labs(x = "Gene's degration rate ", y = "fraction of genes", title = "")+
    theme(panel.grid =element_blank(),panel.background = element_blank())  ## 删去网格线
  
  mp=multiplot(p1, p2, p3, cols=1)
  return(mp)
}

# 生成的结果存储在
# D:\ivan\ivan6-cell\Results\6-Figure 3-wzz\A-C
