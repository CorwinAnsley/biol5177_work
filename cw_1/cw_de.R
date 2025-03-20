

# Install Dependencies (uncomment if not installed)
#install.packages("BiocManager") 
#BiocManager::install("DESeq2")
#BiocManager::install("limma")

# Load required packages
library(DESeq2)
library(ggplot2)
library("limma")

gene_count_table = read.csv("./data/gene_count_matrix.csv",row.names=1)
tran_count_table = read.csv("./data/transcript_count_matrix.csv",row.names=1)

gene_count_table = subset(gene_count_table, select=-c(s1t.c2,s2t.c2))
tran_count_table = subset(tran_count_table, select=-c(s1t.c2,s2t.c2))

sample_names = c('s1','s10','s11','s12','s2','s3','s4','s5','s6','s7','s8','s9')

colnames(gene_count_table) = sample_names
colnames(tran_count_table) = sample_names
