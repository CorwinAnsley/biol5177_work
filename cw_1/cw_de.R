

# Install dependencies (uncomment if not installed)
install.packages("BiocManager") 
BiocManager::install("DESeq2")
BiocManager::install("vsn")
BiocManager::install("limma")

install.packages("ggrepel")
install.packages("ashr")

# Load required packages
library(DESeq2)
library(ggplot2)
library(ggrepel)
library("limma")
library(vsn)
library(ashr)

# Load in and format data
gene_count_table = read.csv("./data/gene_count_matrix.csv",row.names=1)
tran_count_table = read.csv("./data/transcript_count_matrix.csv",row.names=1)

gene_count_table = subset(gene_count_table, select=-c(s1t.c2,s2t.c2))
tran_count_table = subset(tran_count_table, select=-c(s1t.c2,s2t.c2))

sample_names = c('s1','s10','s11','s12','s2','s3','s4','s5','s6','s7','s8','s9')

colnames(gene_count_table) = sample_names
colnames(tran_count_table) = sample_names

col_table = read.csv("./data/Design-Exp.csv",row.names=1)

# Create dds object
dds_genes = DESeqDataSetFromMatrix(countData=gene_count_table,colData=col_table,design= ~ group)
dds_genes

dds_tran = DESeqDataSetFromMatrix(countData=tran_count_table,colData=col_table,design= ~ group)
dds_tran

# Carry out pre-filtering

min_group_size = 3
keep = rowSums(counts(dds_genes) >= 10) >= min_group_size
dds_genes = dds_genes[keep,]

keep = rowSums(counts(dds_tran) >= 10) >= min_group_size
dds_tran = dds_tran[keep,]

# Run the DE analysis 
dds_genes = DESeq(dds_genes)
dds_tran = DESeq(dds_tran)

# Make dispersion plots
plotDispEsts(dds_genes)
plotDispEsts(dds_tran)

# Perform rlog on both objects
rld_genes = rlog(dds_genes, blind=TRUE)
rld_tran = rlog(dds_tran, blind=TRUE)

# Plot PCAs
ggp = plotPCA(rld_genes,intgroup=c("group")) +
  geom_text_repel(aes(label=sample_names),show.legend = FALSE)
ggp

ggp = plotPCA(rld_tran,intgroup=c("group")) +
  geom_text_repel(aes(label=sample_names),show.legend = FALSE)
ggp

# Following steps only for genes

# matrix of log2(raw-counts)
#lgc.raw = log2(counts(dds_genes,normalized=FALSE)+1)

# matrix of log2(normalized-counts) for genes
lgc_norm_genes = log2(counts(dds_genes,normalized=TRUE)+1)

meanSdPlot(assay(rld_genes))

meanSdPlot(lgc_norm_genes)

# Get unshrunken results for LFC = 0
res_lfc0_BvsA = results(dds_genes,contrast=c("group","B","A"))
res_lfc0_CvsA = results(dds_genes,contrast=c("group","C","A"))

# LFC = 1 
res_lfc1_BvsA = results(dds_genes,contrast=c("group","B","A"),lfcThreshold=1)
res_lfc1_CvsA = results(dds_genes,contrast=c("group","C","A"),lfcThreshold=1)

DESeq2::plotMA(res_lfc0,alpha=0.001,main='Gene level B vs A LFC 0')
#plotMA(res_lfc0 ,alpha=0.5,ylim=c(-6,6),main='Gene level B vs A LFC 0')

res_lfcshrink = lfcShrink(dds_genes,contrast=list(c("group_B_vs_A")),type="ashr")
DESeq2::plotMA(res_lfcshrink ,alpha=0.001,ylim=c(-6,6),main='Gene level B vs A LFC 0')

res_lfc1 = results(dds_genes,contrast=c("group","B","A"),lfcThreshold=1)
DESeq2::plotMA(res_lfc1 ,alpha=0.001,main='Gene level B vs A')
plotM