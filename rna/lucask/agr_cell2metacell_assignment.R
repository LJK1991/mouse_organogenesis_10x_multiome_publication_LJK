here::i_am("rna/lucask/agr_cell2metacell_assignment.R")

source(here::here("settings.R"))

p <- ArgumentParser(description='')
p$add_argument('--mapping_metadata',type="character",help='Cell metadata file')
p$add_argument('--metacell_assignment',type="character",help='cell2metacell file')
p$add_argument('--output',type="character",help='output file')
p$add_argument('--trajectories', type="character",help='trajectories')

args <- p$parse_args(commandArgs(TRUE))

args$output <- "/home/lucas/Documents/Work/agrAnalysis/results/rna/metacells/trajectories/Mesen/agr_cell2metacell_assignment.txt.gz"
args$mapping_metadata <- "/home/lucas/Documents/Work/agrAnalysis/results/rna/agr_mapping/sample_metadata_after_mapping.txt.gz"
args$metacell_assignment <- "/home/lucas/Documents/Work/agrAnalysis/results/rna/metacells/trajectories/Mesen/cell2metacell_assignment.txt.gz"


####################################
## reading in metacell assignment ##
####################################

cell2metacell <- fread(args$metacell_assignment)

#########################
## reading in metadata ##
#########################
sample_metadata <- fread(args$mapping_metadata) %>% .[cell%in%cell2metacell$cell]

#sanity check
stopifnot(sample_metadata$cell == cell2metacell$cell)

######################################
## assigning agr cells to metacells ##
######################################
agr_cell2metacell <- merge(sample_metadata, cell2metacell, by = "cell") %>% .[,c("closest.cell","metacell")] %>% setnames("closest.cell","cell")

stopifnot(cell2metacell$metacell == agr_cell2metacell$metacell)

fwrite(agr_cell2metacell, file = file.path(args$output))

