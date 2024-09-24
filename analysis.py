# %%
import pandas as pd
from scipy.io import mmread
import numpy as np
# %%
# Paper original: LIM-domain-only 4 (LMO4) enhances CD8+ T-cell stemness and tumor rejection by boosting IL-21-STAT3 signaling
# https://www.nature.com/articles/s41392-024-01915-z
# Refiere a este otro paper: CRISPR activation and interference screens decode stimulation responses in primary human T cells
# https://www.science.org/doi/10.1126/science.abj4008#bibliography
# De donde sale la data (Gene Expression Omnibus)
# La data la obtuve de aca: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=xxxxxxx
# Donde xxxxxx es GSE174292 por ejemplo
# 
# %%
df = pd.read_excel('data/GSE174255_sgRNA-Read-Counts.xlsx')
# %%
# %%
df = pd.read_csv('data/GSE174284_gene_counts_raw.txt', sep='\t')
# %%
df.columns
# %%
# %% 
# Gene identifiers
features = pd.read_csv('data/GSE190604_features.tsv', sep='\t', header=None)
features.columns = ['Ensembl Gene ID', 'Gene Name', 'Method']

if len(features['Gene Name'].unique()) != len(features):
    print('Gene Name repetidos, por que?')

if len(features['Ensembl Gene ID'].unique()) != len(features):
    print('Ensembl Gene ID repetidos, por que?')
features
# %%
features[features['Gene Name'] == 'CYB561D2']
# %% 
# individual cells in the single-cell RNA sequencing (scRNA-seq) data.
# scRNA-seq is a high-throughput technique used to analyze the transcriptomes of individual cells, allowing researchers to measure gene expression at the single-cell level.
bar_codes_df = pd.read_csv('data/GSE190604_barcodes.tsv', sep='\t', header=None)
bar_codes_df.columns = ['scRNA-seq']
bar_codes_df
# %%
# Matrix
# These values are typically counts of RNA transcripts detected in each cell, 
# indicating how much of each gene is expressed

mat = mmread('data/GSE190604_matrix.mtx')
# %%
mat.shape
# %%
# Get all nonzero elements of sparse matrix
nonzero_elements = mat.data
nonzero_indices = mat.nonzero()

# Create a DataFrame with nonzero elements and their indices
nonzero_df = pd.DataFrame({
    'row': nonzero_indices[0],
    'col': nonzero_indices[1],
    'value': nonzero_elements
})
# %%
nonzero_df
# %% Sparcity
1 - len(nonzero_df)/np.prod(mat.shape)
# %%
# Todos los elementos son positivos?
print(len(nonzero_df) == (nonzero_df['value'] > 0).sum())
# %%
outliers_min_value = 50
nonzero_df[nonzero_df['value'] < outliers_min_value]['value'].hist(bins=outliers_min_value)
# %% 
# Histogram y log axis
outliers_min_value = 3000
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
nonzero_df[nonzero_df['value'] < outliers_min_value]['value'].hist(bins=outliers_min_value-1)
plt.yscale('log')
plt.xlabel('Value')
plt.ylabel('Count (log scale)')
plt.show()
# %%
# %%
cell_ranger_df = pd.read_csv('data/GSE190604_cellranger-guidecalls-aggregated-unfiltered.txt', sep='\t')
# %%
screens_read_counts_df = pd.read_csv('data/GSE190846_supp_CD4_CRISPR_screens_read_counts.tsv', sep='\t')
screens_read_counts_df
# %%
# Esta es la del paper que mandaron originalmente
# 
raw_read_count_matrix_df = pd.read_csv('data/GSE248171_Raw_Read_Count_Matrix.txt', sep='\t').set_index('gene_name')
raw_read_count_matrix_df
# %%
raw_read_count_matrix_df
# %%
values_series = raw_read_count_matrix_df.stack().reset_index()[0]
outliers_min_value = 5000
values_series[(values_series>0)&(values_series<outliers_min_value)].hist()
# plt.yscale('log')
plt.xlabel('Value')
plt.ylabel('Count (log scale)')
plt.show()
# %%
raw_read_count_matrix_df.index.values
# %%
def is_gene_in_features(gene_names):
    indexes = np.array([gene_name.lower() in set(features['Gene Name'].apply(str.lower)) for gene_name in gene_names])
    
    return gene_names[~indexes]
# %%
is_gene_in_features(raw_read_count_matrix_df.index.values[:50])
# %%
gene_name = 'D830031N03Rik'
[gn for gn in features['Gene Name'].apply(str.lower) if (gn in gene_name.lower()) or (gene_name.lower() in gn)]
# %%

# %%
