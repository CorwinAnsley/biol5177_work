

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

library(reshape2)

# Function which gets sorted sig genes and writes to csv
write_sig_results = function(result, p_threshold,filepath){
  # Get subset of sig genes
  result = subset(result,result$padj<p_threshold)
  
  # Sort the results
  result = result[order(result$padj),]
  
  # Write the sorted sig genes to csv
  write.csv(result,file=filepath,quote=FALSE)
}

# Function to create pca plot
pca_graph = function(genes_frame,col_table){
  pca = prcomp(t(as.matrix(sapply(genes_frame,as.numeric))))
  pca_coordinates = data.frame(pca$x)
  
  vars = apply(pca$x, 2, var)
  prop_x = round(vars["PC1"] / sum(vars),4) * 100
  prop_y = round(vars["PC2"] / sum(vars),4) * 100
  x_axis_label = paste("PC1 ", " ",prop_x, "% variance",sep="")
  y_axis_label = paste("PC2 ", " ",prop_y, "% variance",sep="")
  
  ggp = ggplot(pca_coordinates, aes(x=PC1, y=PC2, colour=col_table$group)) +
    geom_point() +
    scale_color_manual(values=as.vector(c("darkcyan", "black", "darkred"))) +
    labs(color = "group") +
    xlab(x_axis_label) +
    ylab(y_axis_label) +
    geom_text_repel(aes(label=sample_names),show.legend = FALSE)
  
  return(ggp)
}

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
plotDispEsts(dds_genes,main='Genes')
dev.print(pdf, './plots/disp_genes.pdf')

plotDispEsts(dds_tran,main='Transcripts')
dev.print(pdf, './plots/disp_tran.pdf')

# Perform rlog on both objects
rld_genes = rlog(dds_genes, blind=TRUE)
rld_tran = rlog(dds_tran, blind=TRUE)

# Plot PCAs
ggp = plotPCA(rld_genes,intgroup=c("group")) +
  geom_text_repel(aes(label=sample_names),show.legend = FALSE) +
  scale_color_manual(values=as.vector(c("darkcyan", "black", "darkred"))) +
  ggtitle("Genes") + 
  theme(aspect.ratio = 1)
ggsave("./plots/pca_genes.pdf", width = 9, height = 9)

ggp = plotPCA(rld_tran,intgroup=c("group")) +
  geom_text_repel(aes(label=sample_names),show.legend = FALSE) +
  scale_color_manual(values=as.vector(c("darkcyan", "black", "darkred"))) +
  ggtitle("Transcripts") + 
  theme(aspect.ratio = 1)
ggsave("./plots/pca_tran.pdf", width = 9, height = 9)

# Following steps only for genes

# matrix of log2(normalized-counts) for genes
lgc_norm_genes = log2(counts(dds_genes,normalized=TRUE)+1)

ggp = meanSdPlot(assay(rld_genes)) +
  theme(aspect.ratio = 1)
ggsave("./plots/meanSd_rlog.pdf", width = 9, height = 9)

ggp = meanSdPlot(lgc_norm_genes) +
  theme(aspect.ratio = 1)
ggsave("./plots/meanSd_log2_norm.pdf", width = 9, height = 9)

# Get unshrunken results for LFC = 0
res_lfc0_BvsA = results(dds_genes,contrast=c("group","B","A"),lfcThreshold=0)
res_lfc0_CvsA = results(dds_genes,contrast=c("group","C","A"),lfcThreshold=0)

# LFC = 1 
res_lfc1_BvsA = results(dds_genes,contrast=c("group","B","A"),lfcThreshold=1)
res_lfc1_CvsA = results(dds_genes,contrast=c("group","C","A"),lfcThreshold=1)

# Get shrunken results
res_lfc0_shrink_BvsA = lfcShrink(dds_genes,contrast=list(c("group_B_vs_A")),type="ashr",lfcThreshold=0)
res_lfc0_shrink_CvsA = lfcShrink(dds_genes,contrast=list(c("group_C_vs_A")),type="ashr",lfcThreshold=0)
res_lfc1_shrink_BvsA = lfcShrink(dds_genes,contrast=list(c("group_B_vs_A")),type="ashr",lfcThreshold=1)
res_lfc1_shrink_CvsA = lfcShrink(dds_genes,contrast=list(c("group_C_vs_A")),type="ashr",lfcThreshold=1)

results = list(res_lfc0_BvsA,
               res_lfc0_CvsA,
               res_lfc1_BvsA,
               res_lfc1_CvsA,
               res_lfc0_shrink_BvsA,
               res_lfc0_shrink_CvsA,
               res_lfc1_shrink_BvsA,
               res_lfc1_shrink_CvsA)

ma_filenames = list('lfc0_BvsA',
                    'lfc0_CvsA',
                    'lfc1_BvsA',
                    'lfc1_CvsA',
                    'lfc0_shrink_BvsA',
                    'lfc0_shrink_CvsA',
                    'lfc1_shrink_BvsA',
                    'lfc1_shrink_CvsA')

for (i in 1:length(results)){
  ggp = DESeq2::plotMA(results[[i]],alpha=0.001,ylim=c(-6,6))
  filepath = paste("./plots/MAplot_",ma_filenames[i], sep = "")
  filepath = paste(filepath,".pdf", sep = "")
  ggsave(filepath, width = 9, height = 9)
}

# Write significant results
write_sig_results(res_lfc1_shrink_BvsA,0.05,'./data/Significant.Genes.BvsA.csv')
write_sig_results(res_lfc1_shrink_CvsA,0.05,'./data/Significant.Genes.CvsA.csv')

# Perform batch correction
ex_design = model.matrix(design(dds_genes),colData(dds_genes))
batch_corrected_genes = limma::removeBatchEffect(assay(rld_genes), batch=colData(dds_genes)$batch, design=ex_design)
# Write to csv
write.csv(batch_corrected_genes,file="./data/BatchCorrected.Rlog.csv",quote=FALSE)

# Convert to data frame for plotting
batch_corrected_genes = data.frame(batch_corrected_genes)

# Create PCA of batch corrected results
ggp = pca_graph(batch_corrected_genes,col_table)
ggsave("./plots/pca_genes_batch_corrected.pdf", width = 9, height = 9)





