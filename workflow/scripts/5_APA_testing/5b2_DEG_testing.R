#https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html

log_smk <- function() {
  if (exists("snakemake") & length(snakemake@log) != 0) {
    log <- file(snakemake@log[1][[1]], open = "wt")
    sink(log, append = TRUE)
    sink(log, append = TRUE, type = "message")
  }
}

log_smk()

# Load libraries
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(apeglm)
library(tibble)

setwd(snakemake@input[["work_dir"]])
script_dir <- getwd()
counts_dir <- snakemake@input[["counts_dir"]]
meta_file <- snakemake@input[["meta_file"]]
outdir <- paste0(script_dir,'/',snakemake@output[[1]],'/')
fprefix <- snakemake@params[["file_prefix"]]
group_one <- snakemake@params[["group_one"]]
group_two <- snakemake@params[["group_two"]]
testing_column <- snakemake@params[["meta_testing_column"]]
cutoff_pct <- snakemake@params[["cutoff_pct"]]
min_cell_per_group <- snakemake@params[["min_cell_per_group"]]


source(paste0(script_dir,'/scripts/5_APA_testing/scPASU_functions_for_differential_testing.R'))

if(!dir.exists(outdir))
{dir.create(outdir, recursive = TRUE)}

# Read merged counts
merged_counts<-read.delim(paste0(script_dir,'/',counts_dir,'/',fprefix,'_counts.txt'))
colnames(merged_counts) <- gsub('\\.','-',colnames(merged_counts))

# Read meta data [save meta data from Seurat object - ensure the set of cells match that of the count matrix]
meta<-read.csv(paste0(script_dir,'/',meta_file),row.names=1)
meta[,testing_column] <- sub(" ", "_", meta[,testing_column])
meta[,testing_column] <- sub("-", "_", meta[,testing_column])

counts_by_cellstate <- stratify_matrix(merged_counts = merged_counts, meta = meta, vars=c(testing_column),
                                     cutoff_pct = cutoff_pct, min_cell_per_group = min_cell_per_group)

###### group 1 v. group 2 testing ######

inputs <- create_test_inputs(test = 'DEG', all_groups = counts_by_cellstate, ident1=c(group_one),ident2=c(group_two), APA.feature.filter = FALSE, 
                             replicate = 'random', nrep = 3 ,p = 0.7)

#colnames(inputs) <- c(paste0('basal_rep',1:3),paste0('intermediate_rep',1:3))

sig_genes <- DEG_DESeq2_pseudobulk(inputs = inputs, comp = c(group_one, group_two), test.used = 'LRT',
                                   padj_thres=0.01,l2fc_frac_thres=log2(1.5),outdir)

# Print session info
sessionInfo()
cat('\n')

cat('End of script ','\n')
