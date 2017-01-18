
All=left_join(pseudogene,result_gai)
Z=na.omit(All)


write.table(Z,file = "D:/ivan/ivan-jusual/pseudogene.txt", append = FALSE, quote =F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")