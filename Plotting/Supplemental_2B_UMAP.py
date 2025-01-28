import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd

def create_and_save_combined_umap(adata_list, genes, titles, save_path, color_map="viridis", vmax=None):
    """
    Creates UMAP plots for multiple AnnData objects and saves them vertically stacked in a single image.

    Parameters:
    - adata_list (list): List of AnnData objects.
    - genes (list): List of genes to plot UMAPs for.
    - titles (list): List of titles for each UMAP plot.
    - save_path (str): Path to save the combined UMAP image.
    - color_map (str): The colormap for the UMAP plot. Default is 'viridis'.

    Returns:
    - None: Saves the UMAP plots to the specified path.
    """
    # Calculate total figure height based on the number of UMAP plots
    total_figsize = (10, 8 * len(adata_list))

    # Create a combined figure
    fig, axes = plt.subplots(len(adata_list), 1, figsize=total_figsize)

    if len(adata_list) == 1:  # Handle case with one AnnData object
        axes = [axes]

    for i, (adata, gene, title) in enumerate(zip(adata_list, genes, titles)):
        # Check if the gene exists in `adata.var_names`
        if gene not in adata.var_names:
            raise ValueError(f"The gene '{gene}' is not found in adata.var_names.")

        # Temporarily add the gene expression as a column in adata.obs
        temp_col = f"{gene}_temp"
        gene_expression = adata[:, gene].X
        
        # Convert sparse matrices to dense, handle NaNs, and ensure flattening
        if hasattr(gene_expression, "toarray"):
            gene_expression = gene_expression.toarray()

        gene_expression = pd.Series(gene_expression.flatten()).fillna(0).values
        adata.obs[temp_col] = gene_expression

        # Plot the UMAP
        sc.pl.umap(
            adata,
            color=temp_col,
            ax=axes[i],
            color_map=color_map,
            vmax=vmax,
            size=100,
            show=False
        )
        axes[i].set_title(title, fontsize=16)
        
        # Remove the temporary column
        del adata.obs[temp_col]


    # Adjust layout and save the combined figure
    plt.tight_layout()
    plt.savefig(save_path, dpi=900, bbox_inches='tight')
    print(f"UMAP plot saved to {save_path}")

def update_obs_sample(adata, mapping):
    """
    Updates the values in adata.obs['sample'] based on a given mapping.

    Parameters:
        adata (AnnData): An AnnData object.
        mapping (dict): A dictionary with keys representing original values 
                        and values representing the new mapped values.

    Returns:
        None: Modifies adata in place.
    """
    if 'sample' not in adata.obs:
        raise KeyError("The .obs attribute does not contain 'sample'.")
    
    # Apply the mapping to update the values
    adata.obs['sample'] = adata.obs['sample'].replace(mapping)


# Adjust text sizes
plt.rcParams['font.size'] = 18  # Adjust base font size
plt.rcParams['axes.titlesize'] = 18  # Title font size
plt.rcParams['axes.labelsize'] = 18  # Axis label font size
plt.rcParams['legend.fontsize'] = 18  # Legend font size

# Read in adata
adata = sc.read("path_to_adata.h5ad")

genes = ["CFTR"]
titles = ["CFTR Expression"]
save_path = "UMAP.png"

# Generate and save the combined UMAP plots
create_and_save_combined_umap([adata], genes, titles, save_path,vmax = 5)
