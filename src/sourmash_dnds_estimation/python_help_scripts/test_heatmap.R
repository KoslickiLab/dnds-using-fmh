library(RColorBrewer) 
library(gplots)
data_matrix <- cbind(c(rnorm(30,-2.5,sd= 0.85),rnorm(30,25,sd= 8),rnorm(30,6,sd= 3)),
                     c(rnorm(30,-2.5,sd= 0.85),rnorm(30,25,sd= 8),rnorm(30,6,sd= 3)),
                     c(rnorm(30,-2.5,sd= 0.85),rnorm(30,25,sd= 8),rnorm(30,6,sd= 3)),
                     c(rnorm(30,-2.5,sd= 0.85),rnorm(30,25,sd= 8),rnorm(30,6,sd= 3)))

breaks = c(5/6* -5, 5/6* -4, 5/6* -3, 5/6* -2, 5/6* -1, 0 ,5/6* 1, 50/6*1, 50/6*2, 50/6*3, 50/6*4, 50/6*5)

setwd("/Users/jzr5814/Desktop")

data_matrix <- as.matrix(read.csv("heatmap_pairwise_mean.csv",
                        row.names = 1))

#hv <- heatmap.2(data_matrix,
#                scale = "none")

#hv
hv <- heatmap(data_matrix, 
                scale="none",
                Rowv=NA,
                Colv=NA,
                col = rev(brewer.pal(11,"RdBu")),
                margins=c(5,5),
                cexRow=0.5, cexCol=1.0, 
                breaks=breaks,
                ylab= "Mutations",
                main = "heatmap",
                key=TRUE,keysize=1.5, trace="none")
