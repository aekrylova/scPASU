#!/usr/bin/env Rscript

log_smk <- function() {
  if (exists("snakemake") & length(snakemake@log) != 0) {
    log <- file(snakemake@log[1][[1]], open = "wt")
    sink(log, append = TRUE)
    sink(log, append = TRUE, type = "message")
  }
}

log_smk()

library(DEXSeq)
library(dplyr)
library(ggplot2)
library(reshape)
library(stringr)
library(data.table)
library(edgeR)
library(Seurat)
library(readr)

setwd(snakemake@input[["work_dir"]])
script_dir <- getwd()
counts_dir <- snakemake@input[["counts_dir"]]
meta_file <- snakemake@input[["meta_file"]]
peak_ref_dir <- snakemake@input[["peak_ref_dir"]]
outdir <- paste0(script_dir,'/',snakemake@output[[1]],'/')
fprefix <- snakemake@params[["file_prefix"]]
group_one <- snakemake@params[["group_one"]]
group_two <- snakemake@params[["group_two"]]
testing_column <- snakemake@params[["meta_testing_column"]]
min_cell_expr_pct <- snakemake@params[["min_cell_expr_pct"]]

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

# Peak ref is jtu$join and final_annotation is the peak column 
peak_ref<-read.table(paste0(script_dir,'/',peak_ref_dir,'/',fprefix,'_final_peak_universe_updated.txt'),header=TRUE,sep='\t')

# Delete original peak column and rename final annotation column as peak
cat('Use final annotations as the peak names \n')
col<-which(colnames(peak_ref)=='final_annotation')
peak_ref$peak<-peak_ref[,col]
peak_ref<-peak_ref[,-col]

## Remove P0 peaks ##
cat('Remove P0 peaks \n')
plist<-strsplit(peak_ref$peak,split=':') %>% sapply(.,'[[',3)
rem<-which(plist=='P0')
peak_ref<-peak_ref[-rem,]

### DEXSeq test ###

apa <- stratify_matrix(merged_counts=merged_counts, meta=meta, vars=c(testing_column),cutoff_pct = 0, min_cell_per_group = 0)
if (!dir.exists(outdir)){dir.create(outdir, recursive = TRUE)}

###### group 1 v. group 2 testing ######

inputs <- create_test_inputs(test='APA',apa,ident1=c(group_one),ident2=c(group_two),min_cell_expr_pct=min_cell_expr_pct,
                             expr.thres=1,pseudo.count=1,replicate='random',nrep=3,p=0.7)

#colnames(inputs) <- c(paste0('im_rep',1:3),paste0('pu_rep',1:3))

### DEXSeq test ###
DEXseq_res <- APA_DEXseq_test(ident1=group_one,ident2=group_two,inputs,min_peak=2,ncpu=4,dispersion.plot.save=TRUE,peak_ref=peak_ref,outdir=outdir)

### T Test ###
DEXseq_res <- APA_TTest(DEXseq_res,ident1=group_one,ident2=group_two)

## Check significance ##
DEXseq_res <- callsig(DEXseq_res,ident1=group_one,ident2=group_two,delta_frac_thres=0.1, padj_thres=0.01,
                      l2fc_frac_thres=log2(1.5),mean_frac_thres=0.05)

comp=paste0(group_one,'_v_',group_two)
saveRDS(DEXseq_res,paste0(outdir,comp,'_APAtest.res.rds'))


### Generate raw counts files and UCSC bedGraph tracks for significant APA genes

ident1 <- sapply(str_split(comp,pattern='_v_',n = 2),`[`,1)
ident2 <- sapply(str_split(comp,pattern='_v_',n = 2),`[`,2)
bedtracks_list <- list()
bedtracks_list[[1]] <- as.data.frame(DEXseq_res$res)[,c('chr','start','end',paste0('mean_',ident1,'_frac'),'strand','int_sig')]
names(bedtracks_list)[1] <- paste0(comp,'-',ident1)
bedtracks_list[[2]] <- as.data.frame(DEXseq_res$res)[,c('chr','start','end',paste0('mean_',ident2,'_frac'),'strand','int_sig')]
names(bedtracks_list)[2] <- paste0(comp,'-',ident2)

# Create bedGraph tracks
for (i in 1:length(bedtracks_list)) {
  bedGraph <- bedtracks_list[[i]] 
  ident <- names(bedtracks_list[i])
  plus <- bedGraph[which(bedGraph$strand == '+'),]
  plus <- plus[,!(names(plus) %in% c('strand','int_sig'))]
  minus <- bedGraph[which(bedGraph$strand == '-'),]
  minus <- minus[,!(names(minus) %in% c('strand','int_sig'))]
  write.table(plus,paste0(outdir,ident,'_plus.bedGraph'),col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
  write.table(minus,paste0(outdir,ident,'_minus.bedGraph'),col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
}

raw.count <- DEXseq_res$counts.raw
write.table(raw.count,paste0(outdir,comp,'_raw.count.txt'),col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')
res <- DEXseq_res$res
res <- res[,-c('dexseq_sig','ttest_sig','all_tests_sig', 'delta_sig', 'perfc_sig')]
write.table(res,paste0(outdir,comp,'_res.txt'),col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
