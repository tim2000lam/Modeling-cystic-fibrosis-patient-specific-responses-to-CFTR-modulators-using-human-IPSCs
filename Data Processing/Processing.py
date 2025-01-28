import os
import numpy as np
import scanpy as sc
import scanpy.external as sce
import bbknn

adata = sc.read_h5ad("complete_adata_combined.h5ad")
adata.var_names_make_unique()

# Define parameters
neighbors_within_batch = 3
n_pcs = 30
min_dist = 0.1
spread = 5.0
trim = None
save_folder = "combined_DMSO_BC"
key = "cell_line"

# Ensure the save folder exists
os.makedirs(save_folder, exist_ok=True)

# Set the seed for reproducibility
np.random.seed(42)
sc.settings.figdir = save_folder

# Annotate mitochondrial, ribosomal, and hemoglobin genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

# Calculate QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

# Filter cells with more than 10% mitochondrial gene expression
adata = adata[adata.obs["pct_counts_mt"] < 10]

# Ensure adata is a full object and not a view
adata = adata.copy()

# Filter cells and genes
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)

# Doublet detection
sc.pp.scrublet(adata, batch_key="sample")

# Save raw counts layer
adata.layers["counts"] = adata.X.copy()

# Normalization to median total counts and log1p transformation
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# Feature selection (highly variable genes)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")

# Perform PCA for dimensionality reduction
sc.pp.pca(adata, n_comps=30)

sce.pp.harmony_integrate(adata, key ="cell_line")

# Compute the neighborhood graph using BBKNN with trim
sc.external.pp.bbknn(
    adata,
    batch_key="cell_line",
    neighbors_within_batch=3,
    n_pcs=30,
    trim=500,
)

# Compute UMAP embedding
sc.tl.umap(adata, min_dist=0.1, spread=5.0)

# Perform Leiden clustering at different resolutions
for res in [0.10, 0.25, 0.5, 0.75,1.00,1.25,1.50]:
    sc.tl.leiden(adata, key_added=f"leiden_res_{res:.2f}", resolution=res, flavor="igraph")

# Save the final AnnData object
adata.write(f"{save_folder}/adata.h5ad")
