#!/usr/bin/env python
# coding: utf-8

# Import the necessary packages.
import diffxpy.api as de
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import seaborn as sns
from gff3 import Gff3
import gffutils
import scanpy as sc
sc.settings.set_figure_params(dpi=2000, facecolor='white',)# Set resolution, styling of figures.
results_file = 'AMP_scanpy.h5ad'

# Read the matrix.
AMP_0h = pd.read_csv('AMP-0h.selected.counts.tsv', sep='\t')
AMP_06275_0h.set_index('gene', inplace=True)
AMP_06275_1h = pd.read_csv('AMP-1h.selected.counts.tsv', sep='\t')
AMP_06275_1h.set_index('gene', inplace=True)
AMP_06275_2h = pd.read_csv('AMP-2h.selected.counts.tsv', sep='\t')
AMP_06275_2h.set_index('gene', inplace=True)
AMP_4h = pd.read_csv('AMP-4h.selected.counts.tsv', sep='\t')
AMP_4h.set_index('gene', inplace=True)

# Read the E.coli bw25113 genome database file.
filename = "Escherichia_coli_bw25113.ASM75055v1.46.gff3"
database_filename = "Escherichia_coli_bw25113-AMP"
db = gffutils.create_db(filename, database_filename, merge_strategy="warning")

# Create a dataframe from the gene list.
genes = [i for i in [f.id for f in db.all_features()] if 'gene' in i]# Parse all the genes out of the database created from the .gff3 file.
gene_names = []
for gene in genes:
    gene_names.append(db[gene][8]['Name'][0]) # Extract the gene name that correspomnd to a gene ID number and append to the list of genes.
gene_dict_rev = dict(zip(gene_names, genes))# Create dictionaries to quickly go from gene ID to gene name and vice versa.
gene_dict = dict(zip(genes, gene_names)) 
all_df = pd.DataFrame()
all_df['gene'] = genes

# Annotate the different samples (0h, 1h, 2h, and 4h after antibiotic exposure).
times = []
t = ['AMP_0h','AMP_1h','AMP_2h','AMP_4h']
for i,frame in enumerate(frames):
    frames[i].columns = [name +'_'+ t[i] for name in frame.columns]
    times = times+ (len(frame.T)*[t[i]])

# Create dataframe with all the samples joined on gene name.
all_df = pd.merge(all_df, AMP_0h, on='gene', how='outer').fillna(0, downcast='infer')
all_df = pd.merge(all_df, AMP_1h, on='gene', how='outer').fillna(0, downcast='infer')
all_df = pd.merge(all_df, AMP_2h, on='gene', how='outer').fillna(0, downcast='infer')
all_df = pd.merge(all_df, AMP_4h, on='gene', how='outer').fillna(0, downcast='infer')

# Map the gene ID from the raw data to the gene names from the .gff3 file stored in gene_dict.
all_df.gene = all_df.gene.apply(lambda x: gene_dict[x])
all_df.set_index('gene', inplace=True)

# Basic quality control
all_df = all_df.T
all_df['time_point'] = times
all_df.loc[:,'Total UMIs'] = all_df.drop(['time_point'],axis = 1).sum(axis=1)# Create a column with the detected UMI count of each cell.
sns.violinplot(data=all_df, x='time_point', y='Total UMIs', cut=0)# Visualize the distrubution of UMI count of each cell.
all_df = all_df[(all_df['Total UMIs'] <3000) & (all_df['Total UMIs'] >300)]# Filter the most extreme outliers.
sns.violinplot(data=all_df, x='time_point', y='Total UMIs', cut=0)# Visualize the distrubution of UMI count of each filtered cell.
all_df['gene_count'] = all_df.drop(['time_point','Total UMIs'],axis = 1).astype(bool).sum(axis=1)# Create a column with the detected gene count of each cell.
sns.violinplot(data=all_df, x='time_point', y='gene_count', cut=0)# Visualize the distrubution of gene count of each filtered cell.
all_df.iloc[:, :-3].to_csv('AMP.csv')# Save the raw data in csv file.

# Read the csv file with scanpy.
data = sc.read('AMP.csv')
data.obs['time'] = all_df['time_point'] # Label samples with the time point.
data.var_names_make_unique() # Make duplicate gene names unique.

# rRNA genes filtering.
rRNA = pd.read_table('rRNA22.csv',sep=',')# Read the list of rRNA genes.
rRNA.head()
pattern = '|'.join(rRNA['gene'].to_list())
rRNAbool = data.var_names.str.contains(pattern)# Annotate the group of rRNA genes
data.obs['percent_rRNA'] = np.sum(data[:,rRNAbool].X,axis=1) / np.sum(data.X,axis=1)# Create a column with the percent of rRNA gene of each cell.
sc.pl.violin(data, keys='percent_rRNA',groupby = 'time')# Visualize the distrubution of the percent of rRNA gene of each cell.
data = data[:,~rRNAbool]# Remove the rRNA genes.

# Calculate and visualize those genes that yield the highest fraction of counts in each single cell, across all cells.
sc.pl.highest_expr_genes(data, n_top=20)

# Low quality cells and genes filtering.
sc.pp.filter_cells(data, min_genes=10)# Remove all cells with fewer than 10 genes.
sc.pp.filter_genes(data, min_cells=3)# Remove all genes expressed in fewer than 3 cells.

# Normalization and dimensionality reduction.
sc.pp.normalize_total(data, target_sum=1e4)# Total-count normalize (library-size correct) the data matrix 𝐗 to 10,000 reads per cell.
sc.pp.log1p(data)# The data are log-transformed with a pseudocount of 1.
adata= data
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)# Identify highly-variable genes to prioritize in further analysis.
sc.pl.highly_variable_genes(adata)# Plot the highly variable genes.
sc.pp.scale(adata, max_value=10)# Scale each gene to unit variance. Clip values exceeding standard deviation 10.
sc.tl.pca(adata, svd_solver='arpack')# Reduce the dimensionality of the data by running principal component analysis (PCA).
sc.pl.pca_variance_ratio(adata, log=True)# Inspect the contribution of single PCs to the total variance in the data. 
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=17)# Compute the neighborhood graph of cells using the PCA representation of the data matrix. 
sc.tl.umap(adata)# Dimensionality reduction using UMAP.
sc.pl.umap(adata, color='time',size=15)# Plot the data with UMAP, with coloring according to time groups.
sc.tl.leiden(adata, resolution =1)# Cluster the neighborhood graph of cells using Leiden graph-clustering method (community detection based on optimizing modularity).
sc.pl.umap(adata, color='leiden',size=15,legend_loc='on data',legend_fontsize =10)# Visualize the data with UMAP, with coloring according to clustering.

# Marker genes and Differential expression.
sc.pl.correlation_matrix(data, groupby='leiden')# Computes a hierarchical clustering for the subclustets ('leiden').

sc.tl.rank_genes_groups(data, 'time')# Make gene markers among time groups.
sc.pl.rank_genes_groups_dotplot(data, n_genes=5,standard_scale ='var')# Plot ranking of top 5 DEGs among time groups using dotplot plot.

sc.tl.dendrogram(adata, 'leiden') # Computes a hierarchical clustering for the subclustets ('leiden').
sc.tl.rank_genes_groups(adata, 'leiden', n_genes=20)# Make gene markers among subclusters.
result = data.uns['rank_genes_groups']
groups = result['names'].dtype.names
degs=pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals','logfoldchanges']}).head(20)# Export DEG file containing 'names', 'logfoldchanges','pvals'.
degs.to_csv('AMP-degs.csv') # Save the top20 DEGs subclusters in csv file.
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5,standard_scale = 'var',groups=('6','9','10'))# Plot ranking of top 5 DEGs in subcluster 6, 9, and 10 using dotplot plot.

# Visualisation of known marker genes.
marker_genes=['glpK','aceA','osmY','fusA','ssrA','fhuA','secY','recA','sulA','ftsN']# Define a set of markers that are known from the literature.
sc.pl.umap(data, color= marker_genes,size=30,legend_loc='on data',legend_fontsize =10,color_map='coolwarm')# Visualize the marker genes using UMAP.
SOS_ROS_genes={'SOS response':[ 'uvrB','recF', 'uvrA', 'dinG', 'uvrD', 'lexA', 'tisB', 'cho', 'dinI', 'dinB', 'yebG', 'umuD','umuC', 'sulA', 'ruvB', 'ruvA', 'recN','recA'],
               'ROS degradation':['sodB', 'sodA', 'katG', 'katE']}# Define a set of markers involved in SOS-response and ROS degradation.
sc.pl.dotplot(data, SOS_ROS_genes, groupby='leiden', standard_scale='var',dendrogram=True)# Plot mean expression levels of marker_genes in subclusters.
Metabolic_genes={'TCA':['mdh', 'gltA', 'acnB', 'acnA','icd','sucA','sucB','lpd','sucD','sucC','sdhA','sdhB','sdhC'],
                   'glycine metabolism':['glyA'],
                   'glycolysis':['gapA'],
                   'glyoxylate cycle':['mdh','gltA','acnB','acnA','aceA','aceB']}# Define a set of markers involved in different metabolic pathways.
sc.pl.dotplot(data, Metabolic_genes, groupby='leiden', standard_scale='var',dendrogram=True)# Plot mean expression levels of marker_genes in subclusters.

