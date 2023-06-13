######################
## Import libraries ##
######################

import os
from re import search
import SEACells

###########################
## Load default settings ##
###########################

#LJK-modify-230410
#removed recognition and just added file.paths
#if search("BI2404M", os.uname()[1]):
#    exec(open('/Users/argelagr/gastrulation_multiome_10x/settings.py').read())
#    exec(open('/Users/argelagr/gastrulation_multiome_10x/utils.py').read())
#elif search("pebble|headstone", os.uname()[1]):
#    exec(open('/bi/group/reik/ricard/scripts/gastrulation_multiome_10x/settings.py').read())
#    exec(open('/bi/group/reik/ricard/scripts/gastrulation_multiome_10x/utils.py').read())
#else:
#    exit("Computer not recognised")
    
exec(open('/home/lucas/Documents/git/mouse_organogenesis_10x_multiome_publication/settings.py').read())
exec(open('/home/lucas/Documents/git/mouse_organogenesis_10x_multiome_publication/utils.py').read())

################################
## Initialise argument parser ##
################################

p = argparse.ArgumentParser( description='' )
p.add_argument( '--anndata',               type=str,                required=True,           help='Anndata file')
p.add_argument( '--metadata',               type=str,                required=True,           help='Cell metadata file')
p.add_argument( '--outdir',               type=str,                required=True,           help='Output directory')
p.add_argument( '--trajectory',               type=str,                required=True,           help='Trajectory')
p.add_argument( '--samples',            type=str, nargs="+",             default="all",             help='samples to use')
p.add_argument( '--percent_metacells',            type=float,              default=0.05,             help='Number of metacells (as a fraction of the total number of cells)')
p.add_argument( '--n_features',            type=int,              default=1500,             help='Number of features')
p.add_argument( '--n_pcs',            type=int,              default=15,             help='Number of PCs')
# p.add_argument( '--test',               action="store_true",                           help='Test mode?' )
# p.add_argument( '--seed',                  type=int,                default=42,               help='Random seed')
# p.add_argument( '--n_iter',       type=int,              default=50,              help='Number of iterations')
args = p.parse_args()

# convert args to dictionary
args = vars(args)

## START TEST ##
# args = {}
# args["anndata"] = io["basedir"] + "/processed/rna/anndata.h5ad"
# args["metadata"] = io["basedir"] + "/results/rna/mapping/sample_metadata_after_mapping.txt.gz"
# args["samples"] = ["all"]
# args["trajectory"] = "nmp"
# args["percent_metacells"] = 0.05
# args["n_features"] = 1500
# args["n_pcs"] = 15
# args["outdir"] = io["basedir"] + "/results/rna/metacells/trajectories/nmp"
## END TEST ##

#####################
## Parse arguments ##
#####################

# I/O
if not os.path.isdir(args["outdir"]): os.makedirs(args["outdir"])
args["outdir"] = Path(args["outdir"])

sc.settings.figdir = args["outdir"] / "pdf"

# Options
#LJK-modify-230410
#changed trajectory specifications.
if args["trajectory"]=="nmp":
  # args["samples"] = ["E8.5_CRISPR_T_KO", "E8.5_CRISPR_T_WT"]
  opts["samples"] = ["E8.0_rep1", "E8.0_rep2", "E8.5_rep1", "E8.5_rep2", "E8.75_rep1", "E8.75_rep2", "E8.5_CRISPR_T_KO", "E8.5_CRISPR_T_WT"]
  opts["celltypes"] = ["Caudal_Mesoderm", "Somitic_mesoderm", "NMP", "Spinal_cord"]
elif args["trajectory"]=="Haemo":
  opts["samples"] = ["Day_4_L","Day_4_R"]
  opts["celltypes"] = ["Epiblast","PGC","Primitive_Streak","Nascent_mesoderm","Mixed_mesoderm","Haematoendothelial_progenitors","Blood_progenitors_1","Blood_progenitors_2"]
elif args["trajectory"]=="Mesen":
  opts["samples"] = ["Day_4_L","Day_4_R"]
  opts["celltypes"] = ["Epiblast","PGC","Primitive_Streak","Nascent_mesoderm","Mixed_mesoderm","Mesenchyme"]
elif args["trajectory"]=="HeaMen":
  opts["samples"] = ["Day_4_L","Day_4_R"]
  opts["celltypes"] = ["Nascent_mesoderm","Mixed_mesoderm","Mesenchyme","Haematoendothelial_progenitors"]
else:
  print("Trajectory not known")
  exit()

print("Infering metacells for %s trajectory..." % args["trajectory"])

print(args["samples"])
print(opts["samples"])
if isinstance(args["samples"],list): 
  if args["samples"][0]=="all":
    args["samples"] = opts["samples"]
  else:
    assert set(args["samples"]).issubset(opts["samples"])
else:
  print('args["samples"] has to be a list')

print(args)

###################
## Load metadata ##
###################

#LJK-modify-230406
#added a piece to change celltype.mapped to celltype. (starting to think that during the mapping there should be a 'celltype' column instead.
metadata = pd.read_table(args["metadata"])
print(metadata.columns)
new_columns = {}
for name in metadata.columns:
	if name == 'celltype.mapped':
		new_columns[name] = 'celltype'
	else:
		new_columns[name] = name
metadata.rename(columns=new_columns, inplace=True)
print(metadata.columns)

metadata = (metadata >>
    # mask(X["pass_rnaQC"]==True, X["pass_atacQC"]==True, X["doublet_call"]==False, X["celltype"].isin(opts["celltypes"])) >>
    mask(X["pass_rnaQC"]==True, X["doublet_call"]==False, X["celltype"].isin(opts["celltypes"])) >>
    mask(X["sample"].isin(args["samples"]))
).set_index("cell", drop=False)

print(metadata.shape)
print(metadata.head())

##################
## Load AnnData ##
##################

adata = load_adata(
	adata_file = args["anndata"], 
	metadata_file = args["metadata"], 
	cells = metadata.index.values, 
	normalise = True, 
  keep_counts = True,
	filter_lowly_expressed_genes = True, 
	set_colors = False
)

# Set colors
adata.obs = adata.obs.rename(columns={"celltype.mapped":"celltype"})
colPalette_celltypes = [opts["celltype_colors"][i.replace(" ","_")] for i in sorted(np.unique(adata.obs['celltype']))]
adata.uns['celltype_colors'] = colPalette_celltypes
#colPalette_stages = [opts["stages_colors"][i.replace(" ","_")] for i in sorted(np.unique(adata.obs['stage']))]
#adata.uns['stage_colors'] = colPalette_stages

#LJK-added-230419
# testing out regression

sc.pp.regress_out(adata, ['nFeature_RNA'])

#######################
## Feature selection ##
#######################

sc.pp.highly_variable_genes(adata, n_top_genes=args["n_features"])

##############################
## Dimensionality reduction ##
##############################

# Load precomputed PCA coordinates
# pca_mtx = pd.read_csv(io["pca_rna"]).set_index("cell", drop=True).loc[adata.obs.index].to_numpy()
# pca_mtx = pd.read_csv(io["pca_atac"]).set_index("cell", drop=True).loc[adata.obs.index].to_numpy()
# adata.obsm["X_pca"] = pca_mtx

# Run PCA
sc.tl.pca(adata, n_comps=args["n_pcs"], svd_solver='arpack')

# Plot PCA
# sc.pl.pca(adata, components=[1,2], color=["celltype","stage"], size=25, legend_loc=None)

# Batch effect correction 
#LJK-modify-230414
#perform batch correction on sample instead of stage.
sc.external.pp.harmony_integrate(adata,"sample", basis='X_pca', adjusted_basis='X_pca_harmony')

# Build kNN graph
sc.pp.neighbors(adata, n_neighbors=25, use_rep='X_pca_harmony')

# Run UMAP
# sc.tl.umap(adata, min_dist=0.5, n_components=2)

# Plot UMAP
# sc.pl.umap(adata, color=["celltype"], size=25, legend_loc=None, save=args["outdir"] / "pdf/umap.pdf")

# Run Force Atlas
#LJK-modify-230411
#removed genotype from the draw_graph color options
sc.tl.draw_graph(adata, layout='fa', init_pos=None)
sc.pl.draw_graph(adata, color=["celltype","stage","sample"], size=20, save="_trajectory.pdf")

########################
## Fit metacell model ##
########################

n_metacells = round(args["percent_metacells"] * adata.shape[0])

print("Fitting SEACells with %d metacells..." % (n_metacells))

#LJK-modify-230406
#https://github.com/dpeerlab/SEACells/blob/main/SEACells/core.py
#according to this it should be 'waypoint_proportion' instead of 'waypt_proportion'
#somehow that doesnt work, but default is 1 anyway so will remove it instead.
##looks like i have to run .construct_kernel matrix() before calling .fit, will have to ask later if this alters anything
model = SEACells.core.SEACells(adata, 
                  build_kernel_on = 'X_draw_graph_fa', 
                  n_SEACells = n_metacells, 
                  n_waypoint_eigs=10,
                  convergence_epsilon = 1e-6)
model.construct_kernel_matrix()
model.initialize_archetypes()
model.fit(min_iter=10, max_iter=150)

adata.obs[['SEACell']].head()

#######################
## Plot model output ##
#######################

model.plot_convergence(save_as=args["outdir"] / "pdf/model_convergence.pdf")

SEACells.plot.plot_2D(adata, key='X_draw_graph_fa', colour_metacells=False, save_as=args["outdir"] / "pdf/forceatlas_highlight_metacells.pdf")

################################################################
## Aggregate counts and plot trajectory at the metacell level ##
################################################################

#test (added list capability)
adata_metacells = SEACells.core.summarize_by_SEACell(adata, SEACells_label='SEACell',celltype_label=["celltype","sample"], summarize_layer='raw')
#LJK-modify-230411
#removed genotype from the select list
#adata_metacells.obs = (adata.obs.loc[adata_metacells.obs.index] >> 
#    select(["sample","celltype"])
#)

sc.pp.normalize_total(adata_metacells)
sc.pp.log1p(adata_metacells)
sc.pp.highly_variable_genes(adata_metacells, n_top_genes=1500)
sc.tl.pca(adata_metacells, n_comps=10)
sc.pp.neighbors(adata_metacells, n_neighbors=15, use_rep='X_pca')
sc.tl.draw_graph(adata_metacells, layout="fa", init_pos=None)
#LJK-modify-230411
#removed genotype from the draw_graph color options
sc.pl.draw_graph(adata_metacells, color=["celltype","sample"], size=150, save="_metacell_trajectory.pdf")

##########
## Save ##
##########

# Save trajectory coordinates
to_save = pd.DataFrame(adata_metacells.obsm["X_draw_graph_fa"], index=adata_metacells.obs_names, columns=["FA1","FA2"])
to_save.to_csv(args["outdir"] / "metacell_trajectory.txt.gz", sep='\t')

# Save cell2metacell assignemnt
to_save = adata.obs[['SEACell']].reset_index()
to_save.columns = ["cell","metacell"]
to_save.to_csv(args["outdir"] / "cell2metacell_assignment.txt.gz", sep="\t", header=True, index=False)

# Save anndata
adata_metacells.write_h5ad(args["outdir"] / "anndata_metacells.h5ad")
