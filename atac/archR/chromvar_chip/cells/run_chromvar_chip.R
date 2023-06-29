here::i_am("atac/archR/chromvar_chip/cells/run_chromvar_chip.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(GenomicRanges))
# suppressPackageStartupMessages(library(chromVAR))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--metadata',  type="character",              help='Cell metadata') 
p$add_argument('--motif_annotation',  type="character",              help='Motif annotation') 
p$add_argument('--atac_peak_matrix',  type="character",              help='ATAC Peak matrix (pseudobulk)') 
p$add_argument('--motifmatcher',  type="character",              help='Motif annotation') 
p$add_argument('--peak_metadata',  type="character",              help='') 
p$add_argument('--background_peaks',  type="character",              help='') 
p$add_argument('--min_number_peaks',     type="integer",    default=30,    help='Minimum number of peaks per TF')
p$add_argument('--min_chip_score',     type="double",    default=0.15,    help='Minimum ChIP score')
p$add_argument('--ignore_negative_values',  action="store_true",  help='Ignote negative ChIP-seq scores')
p$add_argument('--outdir',          type="character",                help='Output directory')
p$add_argument('--test',  action="store_true",  help='Test mode?')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$metadata <- file.path(io$basedir,"results/atac/archR/celltype_assignment/sample_metadata_after_celltype_assignment.txt.gz")
# args$motif_annotation <- "JASPAR"
# args$atac_peak_matrix <- file.path(io$basedir,"processed/atac/archR/Matrices/PeakMatrix_summarized_experiment.rds")
# args$motifmatcher <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/pseudobulk/%s/motifmatchr_virtual_chip.rds",args$motif_annotation))
# args$background_peaks <- file.path(io$basedir,"/processed/atac/archR/Background-Peaks.rds")
# args$peak_metadata <- file.path(io$basedir,"/processed/atac/archR/PeakCalls/peak_metadata.tsv.gz")
# args$min_number_peaks <- 30
# args$min_chip_score <- 0.10
# args$outdir <- file.path(io$basedir,"results/atac/archR/chromvar_chip/cells")
# args$ignore_negative_values <- TRUE
# args$test <- TRUE
## END TEST ##

dir.create(args$outdir, showWarnings=F, recursive=T)

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE]

if (args$test) {
  sample_metadata <- sample_metadata %>% head(n=250)
}

##########################
## Load ATAC PeakMatrix ##
##########################

atac_peakMatrix.se <- readRDS(args$atac_peak_matrix)[,sample_metadata$cell]

# Load peak metadata
peak_metadata.dt <- fread(args$peak_metadata) %>% 
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[,c("idx","score","GC")] %>% 
  setkey(idx) %>% .[rownames(atac_peakMatrix.se)]
stopifnot(peak_metadata.dt$idx==rownames(atac_peakMatrix.se))

# temporary
if (any(!c("start","strand")%in%colnames(rowData(atac_peakMatrix.se)))) {
  tmp <- rownames(atac_peakMatrix.se) %>% strsplit(":") %>% map_chr(2)
  rowData(atac_peakMatrix.se)$start <- tmp %>% strsplit("-") %>% map_chr(1)
  rowData(atac_peakMatrix.se)$end <- tmp %>% strsplit("-") %>% map_chr(2)
}

###############################
## Load motifmatcher results ##
###############################

print("Loading motifmatcher results...")

motifmatcher.se <- readRDS(args$motifmatcher)

stopifnot(c("motifMatches","VirtualChipScores")%in%names(assays(motifmatcher.se)))

###################################################################
## Update motifmatcher results using the virtual ChIP-seq scores ##
###################################################################

if (args$ignore_negative_values) {
  print(sprintf("Number of matches before filtering negative TF binding values: %d",sum(assay(motifmatcher.se,"motifMatches"))))
  # assays(motifmatcher.se) <- assays(motifmatcher.se)[""]
  assay(motifmatcher.se,"motifMatches")[assay(motifmatcher.se,"VirtualChipScores")<0] <- F
  print(sprintf("Number of matches after filtering negative TF binding values: %d",sum(assay(motifmatcher.se,"motifMatches"))))
}

print(sprintf("Number of matches before filtering based on minimum ChIP-seq score: %d",sum(assay(motifmatcher.se,"motifMatches"))))
# assays(motifmatcher.se) <- assays(motifmatcher.se)[""]
assay(motifmatcher.se,"motifMatches")[abs(assay(motifmatcher.se,"VirtualChipScores"))<=args$min_chip_score] <- F
print(sprintf("Number of matches after filtering based on minimum ChIP-seq score: %d",sum(assay(motifmatcher.se,"motifMatches"))))

assays(motifmatcher.se) <- assays(motifmatcher.se)["motifMatches"]

################
## Filter TFs ##
################

# Filter TFs with too few peaks
TFs <- which(colSums(assay(motifmatcher.se,"motifMatches"))>=args$min_number_peaks) %>% names
TFs.removed <- which(colSums(assay(motifmatcher.se,"motifMatches"))<args$min_number_peaks) %>% names

cat(sprintf("%s TFs removed because they don't have enough binding sites: %s", length(TFs.removed), paste(TFs.removed, collapse=" ")))

if (args$test) {
  TFs <- c("FOXA2","MIXL1","GATA1","EOMES","FOXC1")
}
motifmatcher.se <- motifmatcher.se[,TFs]

###########################
## Load background peaks ##
###########################

print("Loading background peaks...")

bgdPeaks.se <- readRDS(args$background_peaks)
tmp <- rowRanges(bgdPeaks.se)
rownames(bgdPeaks.se) <- sprintf("%s:%s-%s",seqnames(tmp), start(tmp), end(tmp))
stopifnot(sort(rownames(bgdPeaks.se))==sort(rownames(atac_peakMatrix.se)))
bgdPeaks.se <- bgdPeaks.se[rownames(atac_peakMatrix.se),]

###################
## Sanity checks ##
###################

#LJK - added - 230619
#motifmatcher.se and atac_peakMatrix.se rownames were not in same order.
#written crude function to get them into order.

sortSummarizedExperiment <- function(se, matrix_name) {
	tmp_assay <- assay(se)
	tmp_colDat <- colData(se)
  #copying rowData is specific per summarized matrix, one column (idx) are the indexes of the array.
  tmp_rowDat <- rowData(se)
  
	#sorting the matrix
	tmp_assay <- tmp_assay[sort(rownames(tmp_assay)),]
	#put assay in list
	tmp_assay_list <- list()
	if (matrix_name == "PeakMatrix"){
		tmp_assay_list$PeakMatrix <- tmp_assay
	} else if (matrix_name == "motifMatches"){
    tmp_assay_list$motifMatches <- tmp_assay
  }
	#make new gRanges
	tmp_gRanges.df <- data.frame(seqnames=sapply(str_split(rownames(tmp_assay),":"),"[[",1),
								 ranges=sapply(str_split(rownames(tmp_assay),":"),"[[",2),
								 strand=rep("*",length(rownames(tmp_assay))))
	tmp_gRanges.df$start <- sapply(str_split(tmp_gRanges.df$ranges,"-"),"[[",1)
	tmp_gRanges.df$end <- sapply(str_split(tmp_gRanges.df$ranges,"-"),"[[",2)
	tmp_gRanges <- makeGRangesFromDataFrame(tmp_gRanges.df)
 
  #sorting the rowData
  tmp_rowDat <- tmp_rowDat[sort(rownames(tmp_rowDat)),]
  tmp_rowDat$idx <- seq(1, length(rownames(tmp_rowDat)))
  

	new_se <- SummarizedExperiment(assays=tmp_assay_list,rowRanges=tmp_gRanges,colData=tmp_colDat)
  rowData(new_se) <- tmp_rowDat
	return(new_se)
}

args$matrix <- "motifMatches"
motifmatcher.se <- sortSummarizedExperiment(motifmatcher.se, args$matrix)
args$matrix <- "PeakMatrix"
atac_peakMatrix.se <- sortSummarizedExperiment(atac_peakMatrix.se, args$matrix)


stopifnot(rownames(atac_peakMatrix.se)==rownames(motifmatcher.se))

#########################################
## Use chromVAR default implementation ##
#########################################

# TAKES TOO MUCH MEMORY

# print("Running chromVAR's default implementation...")

# # prepare data for chromvar
# assayNames(atac_peakMatrix.se) <- "counts"

# # Compute deviations
# chromvar_deviations_chromvar.se <- computeDeviations(
#   object = atac_peakMatrix.se,
#   annotations = motifmatcher.se,
#   background_peaks = assay(bgdPeaks.se)
# )

# # Save
# saveRDS(chromvar_deviations_chromvar.se, file.path(args$outdir,sprintf("chromVAR_chip_%s_chromvar.rds",args$motif_annotation)))

#########################################
## Use ArchR's chromVAR implementation ##
#########################################

print("Running ArchR's chromVAR implementation...")

source(here::here("atac/archR/chromvar/utils.R"))

featureDF <- data.frame(
  rowSums = rowSums(assay(atac_peakMatrix.se)),
  start = rowData(atac_peakMatrix.se)$start,
  end = rowData(atac_peakMatrix.se)$end,
  GC = peak_metadata.dt$GC
)

chromvar_deviations_archr.se <- .customDeviations(
  countsMatrix = assay(atac_peakMatrix.se),
  annotationsMatrix = as(assay(motifmatcher.se),"dgCMatrix"),
  backgroudPeaks = assay(bgdPeaks.se),
  expectation = featureDF$rowSums/sum(featureDF$rowSums),
  prefix = "",
  out = c("deviations", "z"),
  threads = 1,
  verbose = TRUE
)

# Save
saveRDS(chromvar_deviations_archr.se, file.path(args$outdir,sprintf("chromVAR_chip_%s_archr.rds",args$motif_annotation)))
