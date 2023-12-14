library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)
library(data.table)

############################## Files #################################

folder_fmhdnds="/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/marine_bacteria/"
file_fmhdnds='all_dnds_constant.csv'
folder_labels = 'compare_protein/'
file_labels = "compare.prot.15.mat.labels.txt"
file_containments = "compare.prot.15.csv"

############################ Read in files and create matrix for containments ############################

#containments = paste(folder_fmhdnds,folder_labels,file_containments,sep="")
containments_mat <- as.matrix(fread('containments_matrix.csv', sep = ",", header = T)[,-1])
containments_mat
hclust_result <- hclust(dist(containments_mat), method = "average")
x <- colnames(containments_mat)

############################ Read in files and create matrix for dnds ############################

#fmh_dnds_file=paste(folder_fmhdnds,file_fmhdnds,sep="")
dnds_mat <- as.matrix(fread('dNdS_matrix.csv', sep = ",", header = T)[,-1])
dnds_mat
#clean_dnds_data <- dnds_data %>%
#  mutate_all(~ifelse(is.na(.), 1, .))
#dnds_data
heatmap(dnds_mat, Rowv = hclust_result, Colv = FALSE, col = heat.colors(256), scale = "row", main = "Heatmap with Dendrogram")
