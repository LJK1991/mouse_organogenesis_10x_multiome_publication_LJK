here::i_am("rna_atac/gene_regulatory_networks/metacells/trajectories/build_GRN_metacells_trajectory.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(furrr))

################################
## Initialize argument parser ##
################################

#LJK 230626
# Made it so that the entire script can be run with input options

#LJK - note
# is --trajectory even needed? is never called anywhere

p <- ArgumentParser(description='')
p$add_argument('--trajectory', type="character", help="trajectory file (rna)")
p$add_argument('--sce',                type="character",              help='RNA SingleCellExperiment of trajectory(metacells)')
p$add_argument('--tf2gene_virtual_chip', type="character",help="Virtual CHiP-seq")
p$add_argument('--trajectory_name', type="character",help="name of trajectory")
p$add_argument('--min_chip_score', type="numeric", default=0.15, help='minimum CHiP-score')
p$add_argument('--max_distance', type="numeric", default=5e4, help="maximum distance") #to what? lets figure out
p$add_argument('--threads', type="numeric", default=4, help="threads") 
p$add_argument('--outdir',          type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

##########
## TEST ##
##########

# I/O
#io$basedir <- file.path(io$basedir,"test")
#io$trajectory <- file.path("results/rna/trajectories/nmp/nmp_trajectory.txt.gz")
#io$sce <- file.path(io$basedir, 'results/rna/metacells/trajectories/nmp/SingleCellExperiment_metacells.rds')
#io$tf2gene_virtual_chip <- file.path(io$basedir,"results/rna_atac/virtual_chipseq/metacells/trajectories/nmp/JASPAR/TF2gene_after_virtual_chip.txt.gz")
#io$outdir <-  file.path(io$basedir,"results/rna_atac/gene_regulatory_networks/metacells/trajectories/nmp"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
# opts$celltypes <- setdiff(opts$celltypes, c("ExE_endoderm","ExE_ectoderm","Parietal_endoderm"))
#opts$trajectory_name <- "nmp"
#opts$min_chip_score <- 0.15
#opts$max_distance <- 5e4
#opts$ncores <- 4

dir.create(args$outdir, showWarnings=F)

####################################################
## Load TF2gene links based on in silico ChIP-seq ##
####################################################

#tf2gene_chip.dt <- fread(io$tf2gene_virtual_chip) %>%
#  .[chip_score>=opts$min_chip_score & dist<=opts$max_distance] %>% 
#  .[,c("tf","gene")] %>% unique # Only keep TF-gene links
  
tf2gene_chip.dt <- fread(args$tf2gene_virtual_chip) %>%
  .[chip_score>=args$min_chip_score & dist<=args$max_distance] %>% 
  .[,c("tf","gene")] %>% unique # Only keep TF-gene links

# tf2gene_chip.dt[tf=="T"] %>% View

##############################
## Load RNA expression data ##
##############################

#sce <- readRDS(io$sce)
sce <- readRDS(args$sce)

# (Optional) restrict to marker genes
# marker_genes.dt <- fread(io$rna.atlas.marker_genes)
# sce <- sce[rownames(sce)%in%unique(marker_genes.dt$gene),]

##########################
## Filter TFs and genes ##
##########################

TFs <- intersect(unique(tf2gene_chip.dt$tf),toupper(rownames(sce)))
genes <- intersect(unique(tf2gene_chip.dt$gene),rownames(sce))

tf2gene_chip.dt <- tf2gene_chip.dt[tf%in%TFs & gene%in%genes,]

# Fetch RNA expression matrices
rna_tf.mtx <- logcounts(sce)[str_to_title(unique(tf2gene_chip.dt$tf)),]; rownames(rna_tf.mtx) <- toupper(rownames(rna_tf.mtx))
rna_targets.mtx <- logcounts(sce)[unique(tf2gene_chip.dt$gene),]

# Filter out lowly variable genes and TFs
rna_tf.mtx <- rna_tf.mtx[apply(rna_tf.mtx,1,var)>=1,]
rna_targets.mtx <- rna_targets.mtx[apply(rna_targets.mtx,1,var)>=0.1,]

TFs <- intersect(unique(tf2gene_chip.dt$tf),rownames(rna_tf.mtx))
genes <- intersect(unique(tf2gene_chip.dt$gene),rownames(rna_targets.mtx))
tf2gene_chip.dt <- tf2gene_chip.dt[tf%in%TFs & gene%in%genes,]

print(sprintf("Number of TFs: %s",length(TFs)))
print(sprintf("Number of genes: %s",length(genes)))

####################
## run regression ##
####################

# cvfit <- cv.glmnet(x, y, alpha=0)
# opt_ridge <- glmnet(x, y, alpha = 0, lambda  = cvfit$lambda.min)
# df <- data.frame(tf=rownames(opt_ridge$beta), gene=i, beta=as.matrix(opt_ridge$beta)[,1])

if (args$threads>1){
  plan(multicore, workers=args$threads)
  genes_split <- split(genes, cut(seq_along(genes), args$threads, labels = FALSE)) 
} else {
  plan(sequential)
  genes_split <- list(genes)
}

# length(genes_split); sapply(genes_split,length)

GRN_coef.dt <- genes_split %>% future_map(function(genes) {
  tmp <- tf2gene_chip.dt[gene%in%genes]
  genes %>% map(function(i) {
    # print(sprintf("%s (%d/%d)",i,match(i,genes),length(genes)))
    tfs <- tmp[gene==i,tf]
    tfs %>% map(function(j) {
      x <- rna_tf.mtx[j,]
      y <- rna_targets.mtx[i,]
      lm.fit <- lm(y~x)
      data.frame(tf=j, gene=i, beta=round(coef(lm.fit)[[2]],3), pvalue=format(summary(lm.fit)$coefficients[2,4], digits=3))
    }) %>% rbindlist
  }) %>% rbindlist %>% return
}) %>% rbindlist

# save
#fwrite(GRN_coef.dt, file.path(io$outdir,'global_chip_GRN_coef.txt.gz'), sep="\t")
fwrite(GRN_coef.dt, file.path(args$outdir,'global_chip_GRN_coef.txt.gz'), sep="\t")

