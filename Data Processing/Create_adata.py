import scanpy as sc
import anndata as ad#
import numpy as np

# Set seed for reproducibility 
np.random.seed(42)

# Load samples
samples = {
    "CF3_DMSO": "Filtered_bc_matrixes/CF003_control.h5",
    "CF3_ETI": "Filtered_bc_matrixes/CF003_ETI.h5",
    "CF3_LI": "Filtered_bc_matrixes/CF003_Ork.h5",
    "CF9_DMSO": "Filtered_bc_matrixes/CF009_DMSO.h5",
    "CF9_ETI": "Filtered_bc_matrixes/CF009_ETI.h5",
    "CF9_LI": "Filtered_bc_matrixes/CF009_LI.h5",
    "CF9_GE": "Filtered_bc_matrixes/CF009_GE.h5",
    "CF83_ETI": "Filtered_bc_matrixes/CF083_ETI.h5",
    "CF83_LI": "Filtered_bc_matrixes/CF083_LI.h5",
    "CF83_DMSO": "Filtered_bc_matrixes/CF083_DMSO.h5"
}

# Read in data for .h5 samples
adatas = {}
for sample_id, path in samples.items():
    sample_adata = sc.read_10x_h5(path)
    sample_adata.var_names_make_unique()

    # Make observation names unique by adding the sample ID as a suffix
    sample_adata.obs_names = [f"{barcode}_{sample_id}" for barcode in sample_adata.obs_names]

    adatas[sample_id] = sample_adata

# Save the adata
adata.write("adata.h5ad")
