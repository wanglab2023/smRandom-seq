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
results_file = 'technical repeatability.h5ad'

# Read the matrix.
Repeat1 = pd.read_csv('Repeat1.selected.counts.tsv', sep='\t')
Repeat1.columns = ["gene"] + Repeat1.columns[1:].tolist()
Repeat1.set_index('gene', inplace=True)
Repeat2 = pd.read_csv('CIP-T1.selected.counts.tsv', sep='\t')
Repeat2.columns = ["gene"] + Repeat2.columns[1:].tolist()
Repeat2.set_index('gene', inplace=True)

#Read the E.coli bw25113 genome database file.
filename = "Escherichia_coli_bw25113.ASM75055v1.46.gff3"
database_filename = "Escherichia_coli_bw25113-technical repeatability"
db = gffutils.create_db(filename, database_filename, merge_strategy="warning")

# Create a dataframe from the gene list.
genes = [i for i in [f.id for f in db.all_features()] if 'gene' in i]# Parse all the genes out of the database created from the .gff3 file.
gene_names = []
for gene in genes:
    gene_names.append(db[gene][8]['Name'][0])# Extract the gene name that correspomnd to a gene ID number and append to the list of genes.
gene_dict_rev = dict(zip(gene_names, genes))# Create dictionaries to quickly go from gene ID to gene name and vice versa.
gene_dict = dict(zip(genes, gene_names))
all_df = pd.DataFrame()
all_df['gene'] = genes

samples = []
s = ['Repeat1','Repeat2']
for i,frame in enumerate(frames):
    frames[i].columns = [name +'_'+ s[i] for name in frame.columns]
   samples  = samples + (len(frame.T)*[s[i]])

# Create dataframe with all the samples joined on gene name.
all_df = pd.merge(all_df, Repeat1, on='gene', how='outer').fillna(0, downcast='infer')
all_df = pd.merge(all_df, Repeat2, on='gene', how='outer').fillna(0, downcast='infer')

# Map the gene ID from the raw data to the gene names from the .gff3 file stored in gene_dict.
all_df.gene = all_df.gene.apply(lambda x: gene_dict[x])
all_df.set_index('gene', inplace=True)

# Save the raw data in csv file.
all_df = all_df.T
all_df['samples'] = samples 
all_df.iloc[:, :-3].to_csv('technical repeatability.csv')

# Read the csv file with scanpy.
data = sc.read('technical repeatability.csv')
data.obs['sample'] = all_df['samples']# Label samples with the time point.
data.var_names_make_unique()# Make duplicate gene names unique.

# Normalization and dimensionality reduction.
sc.pp.normalize_total(data, target_sum=1e4)# Total-count normalize (library-size correct) the data matrix 𝐗 to 10,000 reads per cell.
sc.pp.log1p(data)# The data are log-transformed with a pseudocount of 1.
adata= data
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)# Identify highly-variable genes to prioritize in further analysis.
sc.pl.highly_variable_genes(adata)# Plot the highly variable genes.
sc.pp.scale(adata, max_value=10)# Scale each gene to unit variance. Clip values exceeding standard deviation 10.
sc.tl.pca(adata, svd_solver='arpack')# Reduce the dimensionality of the data by running principal component analysis (PCA).
sc.pl.pca_variance_ratio(adata, log=True)# Inspect the contribution of single PCs to the total variance in the data. 
sc.pp.neighbors(adata_normed, n_neighbors=10, random_state=1, n_pcs=5) #  Compute the neighborhood graph of cells using the PCA representation of the data matrix. 
sc.tl.umap(adata_normed, min_dist=0.7, spread=2.0, random_state=1, n_components=2,)# Cluster the neighborhood graph of cells using Leiden graph-clustering method (community detection based on optimizing modularity).
sc.pl.umap(adata_normed, color='sample',size=20)# Plot the data with UMAP, with coloring according to samples.
