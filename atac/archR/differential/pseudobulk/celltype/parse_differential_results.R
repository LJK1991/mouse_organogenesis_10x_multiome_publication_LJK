here::i_am("atac/archR/differential/pseudobulk/celltype/parse_differential_results.R")

# Load default settings
source(here::here("settings.R"))

#LJK
# for some reason this script seems to crashes in snakemake, unable to either make the outdir ater which it complains that the outfile cannot be made.
# when running it manually with the input given by snakemake it runs perfectly fine. not sure what is going on.

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--diff_results_dir',   type="character",     help='File')
p$add_argument('--outdir',             type="character",     help='File')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$diff_results_dir <- file.path(io$basedir,"results/atac/archR/differential/metacells/celltype/PeakMatrix")
# args$outdir <- file.path(io$basedir,"results/atac/archR/differential/metacells/celltype/PeakMatrix/parsed")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings=F, recursive=T)

################################################
## Load differential expression and fetch TFs ##
################################################

diff_results_list <- list()

# i <- "Visceral_endoderm"; j <- "Surface_ectoderm"
for (i in 1:length(opts$celltypes)) {
  for (j in i:length(opts$celltypes)) {
    
    if (i!=j) {
      file <- file.path(args$diff_results_dir,sprintf("%s_vs_%s.txt.gz",opts$celltypes[[i]],opts$celltypes[[j]]))
      if (file.exists(file)) {
        tmp <- fread(file) %>% .[,c("celltypeA","celltypeB"):=list(opts$celltypes[[i]],opts$celltypes[[j]])]
        diff_results_list[[sprintf("%s_vs_%s",opts$celltypes[[i]],opts$celltypes[[j]])]] <- tmp
      } else {
        print(sprintf("%s not found...",file))
      }
    }
  }
}
 
##########
## Save ##
##########

fwrite(rbindlist(diff_results_list), file.path(args$outdir,"diff_results.txt.gz"), sep="\t", quote=F, na="NA")

##########
## TEST ##
##########
#LJK - modify - 230613
#why is this done with a static path and not an argument like args$ / opts$
#"/bi/group/reik/ricard/data/gastrulation_multiome_10x/test/results/atac/archR/differential/pseudobulk/celltype/PeakMatrix/parsed/diff_results.txt.gz"
#additionally if it is "TEST" why is it not commented out?, will leave it in for now, just in case it does something important, which i dont think it does
tmp <- fread(file.path(args$diff_results_dir, "parsed/diff_results.txt.gz")) %>% .[abs(logFC)>=2 & padj_fdr<=0.01 & (mean_groupA>=2.5 | mean_groupB>=2.5)] %>% .[,.N,by="feature"]

sum(tmp$N>=2)