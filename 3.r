P1=REP1@ratesFirstGuess@assayData$exprs
P2=REP2@ratesFirstGuess@assayData$exprs
write.table(P1, file = "D:/Rep1.txt", append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
write.table(P2, file = "D:\REP2.txt", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")




P=cbind.fill(P1,P2)
P[apply(P==0, FUN = any, 1), ] = NA
P[apply(!is.finite(P), FUN = any, 1), ] = NA
P=as.data.frame(na.omit(P))
P1=P[,1:15]
P2=P[,16:30]
plot(P1[,1],P2[,1])


length(P1[,1])
length(P2[,1])
