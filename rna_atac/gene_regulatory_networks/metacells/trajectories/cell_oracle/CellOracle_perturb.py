#pasted everything in here, and removed the 'out:'
#attempt to read and make it ito a py script which can be put into the snakemake pipeline?

###############
## Importing ##
###############
import os
from re import search
from dfply import *
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import importlib
import argparse
from itertools import compress

from pathlib import Path

import celloracle as co
co.__version__

%load_ext autoreload
%autoreload 2

#exec(open('/Users/argelagr/gastrulation_multiome_10x/settings.py').read())
#exec(open('/Users/argelagr/gastrulation_multiome_10x/utils.py').read())
exec(open('/media/draco/lucask/chapter_three/mouse_organogenesis_10x_multiome_publication/settings.py').read())
exec(open("/media/draco/lucask/chapter_three/mouse_organogenesis_10x_multiome_publication/utils.py").read())


###################
## Get arguments ##
###################
p = argparse.ArgumentParser()
p.add_argument('--metadata',type=str,help='RNA (metacells) trajectory metadata')
p.add_argument('--anndata',type=str,help='RNA (metacells) trajectory anndata')
p.add_argument('--tf2gene',type=str,help='TF2gene after virtual ChIP')
#p.add_argument('--tf_markers_up', type=str,help='atlas markers tf')
p.add_argument('--trajectory_file',type=str,help='trajectory file')
p.add_argument('--trajectory_name',type=str,help='name of trajectory')
p.add_argument('--samples',type=str,help='samples to use')
p.add_argument('--outdir',type=str,help='output directory')
p.add_argument('--chip_score',type=float, default=0.20, help='ChIP score')
p.add_argument('--max_distance',type=int, default=50000, help='Max genomic distance')
#p.add_argument('--scale_simulation',type=float,default=0.55, help='')
args = parser.parse_args()

os.mkdir(args["outdir"])

#io['anndata'] = Path(io["basedir"]) / "results/rna_atac/rna_vs_acc/trajectories/blood_trajectory/anndata.h5ad"
#io['metadata'] = Path(io["basedir"]) / "results/rna/metacells/trajectories/nmp/metacells_metadata.txt.gz"
#io['anndata'] = Path(io["basedir"]) / "results/rna/metacells/trajectories/nmp/anndata_metacells.h5ad"
#io['tf2gene'] = Path(io["basedir"]) / 'results/rna_atac/virtual_chipseq/metacells/trajectories/nmp/JASPAR/TF2gene_after_virtual_chip.txt.gz'
#io["outdir"] = Path(io["basedir"]) / "results/rna_atac/gene_regulatory_networks/metacells/trajectories/nmp/celloracle"

#########################
## Additional settings ##
#########################
#visual setting
%config InlineBackend.figure_format = 'retina'
%matplotlib inline
plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300
# %%capture
# sc.settings.verbosity = 3
# sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(8, 7), facecolor='white')

if args["samples"] == "all":
	args["samples"] = opts["samples"]
else:
	assert all(x in opts["samples"] for x in args["samples"]) in , "Unkown samples selected"
	
#opts["samples"] = [
#	"E7.5_rep1",
#	"E7.5_rep2",
#	"E7.75_rep1",
#	"E8.0_rep1",
#	"E8.0_rep2",
#	"E8.5_rep1",
#	"E8.5_rep2",
#	"E8.75_rep1",
#	"E8.75_rep2",
#	"E8.5_CRISPR_T_KO",
#	"E8.5_CRISPR_T_WT"
#]
# opts["samples"] = ["E8.5_rep1"]


#Define your celltypes per trajectory
if args["trajectory_name"] == 'nmp':
	args["celltypes"] = ["Caudal_Mesoderm", "Somitic_mesoderm", "NMP", "Spinal_cord"]
	args["TOI"] = ["Sox2","T"]
elif args["trajectory_name"] == "Haemo":
	args["celltypes"] = ["Epiblast","PGC","Primitive_Streak","Nascent_mesoderm","Mixed_mesoderm","Haematoendothelial_progenitors","Blood_progenitors_1","Blood_progenitors_2"]
	args["TOI"] = ["Runx1","Tal1","Nfe2","Fli1","Gata2","Gata1"]
elif args["trajectory_name"] == "Mesen":
	args["celltypes"] = ["Epiblast","PGC","Primitive_Streak","Nascent_mesoderm","Mixed_mesoderm","Mesenchyme"]
    args["TOI"] = ["Hand2","Otx2","Elf3","Mesp1","Lhx1"]
elif args["trajectory_name"] == "HeaMen":
	args["celltypes"] = ["Nascent_mesoderm","Mixed_mesoderm","Mesenchyme","Haematoendothelial_progenitors"]
    args["TOI"] = ["Etv2","Mixl1","Ebf1","Eomes","Snai2"]

assert all(x in args["celltypes"] for x in opts["celltypes"]), "Detected unknown celltypes"


##################
## Loading data ##
##################

metadata = (pd.read_table(args["metadata"]) >>
    mask(X["sample"].isin(args["samples"]), X["celltype"].isin(args["celltypes"]))
).set_index("metacell", drop=False)
metadata.shape

# adata = load_adata(
#     adata_file = io["anndata"], 
#     metadata_file = io["metadata"], 
#     cells = metadata.index.values, 
#     normalise = True, 
#     keep_counts = True,
# 	filter_lowly_expressed_genes = True, 
# 	set_colors = False
# )

adata = load_adata(
    adata_file = args['anndata'], 
    normalise = True, 
    keep_counts = True,
	set_colors = False
)
adata

colPalette_celltypes = [opts["celltype_colors"][i.replace(" ","_")] for i in sorted(np.unique(adata.obs['celltype']))]
adata.uns['celltype_colors'] = colPalette_celltypes

colPalette_celltypes = [opts["stage_colors"][i.replace(" ","_")] for i in sorted(np.unique(adata.obs['stage']))]
adata.uns['stage_colors'] = colPalette_celltypes

#LJK - note
#why this step? checking dataframe?
adata.X[:5,:5].todense()

#######################
## Feature selection ##
#######################
sc.pp.highly_variable_genes(adata, n_top_genes=25000)

#Subset anndata with highly variable genes
adata = adata[:,adata.var["highly_variable"]]
adata

####################
## TF of interest ##
####################

#print like this so we can see which factors are in or not in the 'log'
np.array(args["TOI"])
np.isin(np.array(args["TOI"]),adata.var_names)

#remove TOI if not in list, otherwise bound to give errors
args["TOI"] = list(compress(args["TOI"],np.isin(np.array(args["TOI"]),adata.var_names))

assert all(np.isin(np.array(args["TOI"]),adata.var_names)), "not all TF are in adata.var_names"

#########
## PCA ##
#########
sc.tl.pca(adata, n_comps=15, svd_solver='arpack')

#############################
## Batch effect correction ##
#############################

#LJK - note
#why no batch correction?, will put it in to be sure (but on sample)
# sc.external.pp.harmony_integrate(adata, "stage", basis='X_pca', adjusted_basis='X_pca_harmony')
sc.external.pp.harmony_integrate(adata, "sample", basis='X_pca', adjusted_basis='X_pca_harmony')

################
## k-NN graph ##
################
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=15, use_rep="X_pca")

###########################
## Force-directed layout ##
###########################
sc.tl.draw_graph(adata, layout='fa', init_pos=None)
#adata.obsm["X_draw_graph_fa"] = pd.read_csv("/Users/argelagr/data/gastrulation_multiome_10x/test/results/rna/metacells/trajectories/nmp/metacell_trajectory.txt.gz", sep="\t", index_col=0).loc[adata.obs.index.values].values
adata.obsm["X_draw_graph_fa"] = pd.read_csv(args["trajectory_file"], sep="\t", index_col=0).loc[adata.obs.index.values].values
adata.obsm

sc.pl.draw_graph(adata, color=["celltype"], size=200, legend_loc="on data")

#################################
## Prepare data for celloracle ##
#################################
adata.X = adata.layers["raw"].copy()
adata.X[1:5,1:5].todense()

TF_markers_pandas = pd.read_csv(io["rna_atlas_marker_TFs_all"], sep="\t") >> (
    mask(X["score"]>=0.75, X["celltype"].isin(args["celltypes"])) >>
    mutate(gene=X["gene"].str.upper())
)

TFs = np.unique(TF_markers_pandas.gene.values)
TFs

########################
## Load TF2gene links ##
########################
tf2gene_df = pd.read_csv(args['tf2gene'],header=0, sep="\t")
tf2gene_df

TFs

########################################
## Filter target genes to also be TFs ##
########################################
tf2gene_df = pd.read_csv(args['tf2gene'],header=0, sep="\t") >> (
    mask(X["chip_score"]>=args['chip_score'], X["dist"]<=args['max_distance']) >> 
    mask(X["gene"].str.upper().isin(TFs),X["tf"].str.upper().isin(TFs)) >> 
    select([1,3]) >> 
    mutate(gene=X["gene"].str.upper())
)
tf2gene_df.shape
tf2gene_df.head()

########################################
## Filter genes in the anndata object ##
########################################
tf2gene_df["gene"].values

all_genes = np.unique(np.concatenate((tf2gene_df["tf"].values,tf2gene_df["gene"].values)))
len(all_genes)
# adata.var.index.str.upper().isin()
adata.var.index = adata.var.index.str.upper()
adata = adata[:,all_genes]
tmp = adata.var["gene"].str.upper().isin(tf2gene_df["tf"].values).values
adata.var["gene"][tmp] = adata.var["gene"][tmp].str.upper()
adata.var["gene"].str.upper().isin(tf2gene_df["tf"].values).sum()

################
## Filter TFs ##
################
tf2gene_df = tf2gene_df[tf2gene_df["tf"].isin(adata.var_names)]
tf2gene_df.shape

##################
## Filter genes ##
##################
tf2gene_df = tf2gene_df[tf2gene_df["gene"].isin(adata.var_names)]
tf2gene_df.shape

#Convert to a dictionary
tf2gene_dic = tf2gene_df.groupby('gene')['tf'].apply(lambda x: x.tolist()).to_dict()
# np.isin(np.array(list(tf2gene_dic.keys())),adata.var_names).sum()
len(tf2gene_dic.keys())

############################
## Initiate Oracle object ##
############################
oracle = co.Oracle()

##Load gene expression data into oracle object.
adata.var["symbol"] = adata.var.index.values
#adata.var["isin_top1000_var_mean_genes"] = True

#LJK - note
#what does the dummy column do? also why is it being assigned a color?
##create dummy column
adata.obs["foo"] = "foo"
adata.obs["foo"] = adata.obs["foo"].astype('category')
adata.uns["foo_colors"] = np.array(["#9e6762"])

adata.obsm
oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name="foo", embedding_name="X_draw_graph_fa")

####################################
## Load TFinfo into oracle object ##
####################################
oracle.import_TF_data(TFdict=tf2gene_dic)

#oracle.adata.var.loc[:,["isin_TFdict_targets","isin_TFdict_regulators"]]
oracle.adata.var.loc[:,"isin_TFdict_targets"].sum()
oracle.adata.var.loc[:,"isin_TFdict_regulators"].sum()


##Knn imputation > ??
#########
## PCA ##
#########
# oracle.perform_PCA()
# Select important PCs
# plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
# n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
# plt.axvline(n_comps, c="k")
#KNN imputation. It creates the attributes:
# - knn (scipy.sparse.csr_matrix): knn contiguity matrix
# - knn_smoothing_w (scipy.sparse.lil_matrix): the weights used for the smoothing
# k = 25
# oracle.knn_imputation(n_pca_dims=15, k=k, balanced=True, b_sight=k*8, b_maxl=k*4, n_jobs=4)
# oracle.knn.shape
# oracle.adata.layers

oracle.adata.layers["normalized_count"].todense()[:5,:5]

#####################
## GRN calculation ##
#####################
#The next step is constructing a cluster-specific GRN for all clusters.
#The "get_links" function returns GRNs as a Links object which stores inferred GRNs and the corresponding metadata. 
#You can do network analysis with the Links object.

importlib.reload(co)
# co.__file__

links = oracle.get_links(cluster_name_for_GRN_unit="foo", alpha=10, verbose_level=10, test_mode=False)
links

links.links_dict["foo"].shape
links.links_dict["foo"]

#################
## Exploration ##
#################
links.filter_links(p=0.0001, weight="coef_abs")
links.links_dict["foo"].shape

##########
## Save ##
##########
#LJK - note 
#why is this not being saved?
#io['links_outfile'] = io['outdir'] / 'test.celloracle.links'
# links.to_hdf5(str(io['links_outfile']))
#io['oracle_outfile'] = io['outdir'] / "test.celloracle.oracle"
# oracle.to_hdf5(str(io['oracle_outfile']))
links.to_hdf5(str(args['outdir'] + "/test.celloracle.links"))
oracle.to_hdf5(str(args['outdir'] + "/test.celloracle.oracle"))

##Explore
for TF in args["TOI"]:
	links.links_dict["foo"] >> mask(X.source==TF) >> arrange("coef_mean")
	
#links.links_dict["foo"] >> mask(X.source=="T") >> arrange("coef_mean")
#links.links_dict["foo"] >> mask(X.source=="SOX2") >> arrange("coef_mean")

oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)

#LJK - note
#not sure if everything below should be in a loop?
#will do it and then see later
for TF in args["TOI"]:
	TF_dir = str(args["outdir"] + "/" + TF)
	os.mkdir(TF_dir)
	sc.pl.draw_graph(oracle.adata, color=[TF,'celltype'], size=200, layer="normalized_count", cmap="viridis", save=str(TF_dir + "expression_across_trajectory.pdf"))
	oracle.simulate_shift(perturb_condition={TF: 0.0}, n_propagation=3)
	oracle.estimate_transition_prob(n_neighbors=15, knn_random=True, sampled_fraction=1)
	oracle.calculate_embedding_shift(sigma_corr = 0.05)
	oracle.calculate_p_mass(smooth=0.8, n_grid=40, n_neighbors=15)
	thresh_suggestions = oracle.suggest_mass_thresholds(n_suggestion=8)
	plt.savefig(fname=str(TF_dir + "/threshold_suggestions.pdf"),thresh_suggestions,format="pdf")
	
	mass_filter  = oracle.calculate_mass_filter(min_mass=0.01, plot=True)
	plt.savefig(fname=str(TF_dir + "/mass_filter.pdf"),mass_filter,format="pdf")
	
	#is this a number that should be an argment option?
	#scale_simulation = 0.55
	fig, ax = plt.subplots(2, 5,  figsize=[13, 6])
	ax_pos = 0
	for scale in np.arange(0.1,1,0.2):
		oracle.plot_simulation_flow_on_grid(scale=scale, ax=ax[0,ax_pos])
		oracle.plot_simulation_flow_random_on_grid(scale=scale, ax=ax[1,ax_pos])
		ax_pos += 1
	
	ax[0,0].set_title(f"Perturbation simulation results: {TF} KO")
	ax[1,0].set_title(f"Perturbation simulation with randomized GRNs")
	
	plt.savefig(fname=str(TF_dir + "/scale_simulations.pdf"),fig,format="pdf")

	cell_colors = [opts["celltype_colors"][i] for i in oracle.adata.obs["celltype"].values]
	for scale in np.arange(0.1,1,0.2):
		fig, ax = plt.subplots(figsize=[10, 8])
		ax.scatter(oracle.embedding[:,0], oracle.embedding[:,1], c=cell_colors, alpha=0.45)
		oracle.plot_simulation_flow_on_grid(scale=scale, ax=ax, show_background=False)
		plt.savefig(str(args["outdir"] + "/TF_force_scale_" + scale +".pdf"))
	
#gene = "T"
#sc.pl.draw_graph(oracle.adata, color=[gene,'celltype'], size=200, layer="normalized_count", cmap="viridis")
#oracle.simulate_shift(perturb_condition={gene: 0.0}, n_propagation=3)
#oracle.estimate_transition_prob(n_neighbors=15, knn_random=True, sampled_fraction=1)
#oracle.calculate_embedding_shift(sigma_corr = 0.05)
#oracle.calculate_p_mass(smooth=0.8, n_grid=40, n_neighbors=15)
#oracle.suggest_mass_thresholds(n_suggestion=8)

#oracle.calculate_mass_filter(min_mass=0.01, plot=True)

#scale_simulation = 0.55
#oracle.plot_simulation_flow_on_grid(scale=scale_simulation)

# fig, ax = plt.subplots(1, 2,  figsize=[15, 7])
# oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
# ax[0].set_title(f"Perturbation simulation results: {gene} KO")
# oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
# ax[1].set_title(f"Perturbation simulation with randomized GRNs")
# plt.show()

#cell_colors = [opts["celltype_colors"][i] for i in oracle.adata.obs["celltype"].values]
#fig, ax = plt.subplots(figsize=[10, 8])
#ax.scatter(oracle.embedding[:,0], oracle.embedding[:,1], c=cell_colors, alpha=0.45)
#oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax, show_background=False)
#plt.savefig(io['outdir'] / "test.pdf")

##LJK - note
##There is more on the celloracle website https://morris-lab.github.io/CellOracle.documentation/notebooks/05_simulation/Gata1_KO_simulation_with_Paul_etal_2015_data.html
##but this might/may be unecesarry or even just too much.

##make completion token for snakemake
with open(str(args["outdir"] + "/CellOracle_" + args["trajectory_name"] + "_completion-token.txt"),'w') as f:
    pass
