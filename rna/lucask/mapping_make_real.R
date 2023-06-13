suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(future))

source("/home/lucas/Documents/git/mouse_organogenesis_10x_multiome_publication/utils.R")

args <- list()
args$transcriptome_samples <- c("Day_4_R","Day_4_L")
args$transcriptome_dir <- "/home/lucas/Documents/Work/agrAnalysis/data/original/"
args$atac_samples <- c("E7.5_rep1","E7.5_rep2","E7.75_rep1")
args$atac_dir <- "/home/lucas/Documents/Work/mouse_ATAC/ftp1.babraham.ac.uk/data/original"
args$cell2metacell <- "/home/lucas/Documents/Work/agrAnalysis/results/rna/agr_mapping" #agr_mapping_mnn_Day_4_L.txt.gz


#loading atac data
atac_data <- list()
for (sample in args$atac_samples){
  atac_data[[sample]] <- Read10X(sprintf("%s/%s/filtered_feature_bc_matrix",args$atac_dir,sample))  
}

atac_granges <- list()
for (sample in args$atac_samples){
  atac_peaks <- fread(file = sprintf("%s/%s/atac_fragments.tsv.gz",args$atac_dir,sample), skip = 52) %>% 
    .[,c("V1","V2","V3")] %>% setnames(old=c("V1","V2","V3"), new = c("chr", "start", "end"))
  atac_granges[[sample]] <- makeGRangesFromDataFrame(atac_peaks) 
}
rm(atac_peaks)
combined_peaks <- reduce(x = c(atac_granges[[args$atac_samples[1]]],atac_granges[[args$atac_samples[2]]],atac_granges[[args$atac_samples[3]]])
rm(atac_granges)



sample_data <- Read10X(sprintf("%s%s/outs/filtered_feature_bc_matrix",args$transcriptome_dir,args$transcriptome_samples[2]))
cell2metacell <- fread(file=file.path(args$cell2metacell,sprintf("/agr_mapping_mnn_%s.txt.gz",args$transcriptome_samples[2])))

atac_mtx <- sparseMatrix()

for (sample in args$samples){
  sample_data <- Read10X(sprintf("%s/%s/outs/filtered_feature_bc_matrix",args$transcriptome_dir,sample))
  cell2metacell <- fread(file=file.path(args$cell2metacell,sprintf("/agr_mapping_mnn_%s.txt.gz",sample)))
}
