import scanpy as sc

# Import adata object
adata = sc.read("adata_leiden.h5ad")

# Define cluster annotations
cluster_annotations = {
    '0': 'Fibroblast',
    '1': 'Fibroblast',
    '2': 'Fibroblast',
    '3': 'Proliferating ciliated progenitor',
    '4': 'Goblet',
    '5': 'Fibroblast',
    '6': 'Proliferating secretory/ciliated progenitor',
    '7': 'Secretory',
    '8': 'Secretory club',
    '9': 'Fetal progenitor',
    '10': 'Basal',
    '11': 'Basal',
    '12': 'Ciliated'
}


# Add a new column with cell type annotations based on Leiden clusters
adata.obs['cell_type'] = adata.obs['leiden_res_0.75'].map(cluster_annotations)

# Visualize UMAP with the newly assigned cell type annotations
sc.pl.umap(adata, color='cell_type', legend_loc='on data', save='_new_coarse_cell_type_annotations.png')

adata.write("bc_adata.h5ad")

