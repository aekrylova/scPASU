log_smk <- function() {
  if (exists("snakemake") & length(snakemake@log) != 0) {
    log <- file(snakemake@log[1][[1]], open = "wt")
    sink(log, append = TRUE)
    sink(log, append = TRUE, type = "message")
  }
}

log_smk()

library(data.table)
library(dplyr)
library(Seurat)
library(stringr)

# Local run
counts_dir <- snakemake@input[["counts_dir"]]
work_dir <- snakemake@input[["work_dir"]]
fprefix <- snakemake@params[["file_prefix"]]
seurat_obj_path <- snakemake@params[["seurat_obj"]]
outdir <- snakemake@output[[1]]
samples <- snakemake@params[["samples"]]

if(!dir.exists(outdir))
{dir.create(outdir, recursive = TRUE)}

f <- paste0(work_dir,'/',counts_dir,'/',samples,'/outs/raw_feature_bc_matrix/')

counts<-lapply(f,function(x) {
  mtx <- Read10X(data.dir = x)
  mtx <- as.data.frame(mtx)
  return(mtx)
}
)
cat("Converted to dataframe \n")

names(counts) <- samples

# Modify cell barcodes in accordance with previous scRNA-seq analyses
if (seurat_obj_path != '') {
  cat('Modifying cell barcodes to match Seurat when integrated \n')
  z <- readRDS(seurat_obj_path)
  meta <- z@meta.data
  write.table(meta,paste0(outdir,'/',fprefix,'_meta.txt'),col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')
  meta$sample_idx <-as.numeric(sub(".*_([0-9]+)$", "\\1", row.names(meta)))
  sample_idx <- unique(meta$sample_idx)
  names(sample_idx) <- unlist(meta[match(sample_idx,meta$sample_idx),'sample'])
  
  for (i in 1:length(counts)) {
    for (s in names(sample_idx)) {
      if (names(counts)[i]==s) {
        colnames(counts[[i]]) <- sub(".*_(.*)", "\\1",colnames(counts[[i]]))
        colnames(counts[[i]]) <- paste0(colnames(counts[[i]]),'_',sample_idx[names(sample_idx)==s])
      }
    }
  }
  names(counts) <- NULL
  merged_counts<-do.call(cbind,counts)
  stopifnot(identical(sort(colnames(merged_counts)),sort(row.names(z@meta.data))))
} else {
  names(counts) <- NULL
  merged_counts<-do.call(cbind,counts)
}

write.table(merged_counts,paste0(work_dir,'/',outdir,'/',fprefix,'_counts.txt'),row.names = TRUE, col.names = TRUE, quote = FALSE, sep='\t')
