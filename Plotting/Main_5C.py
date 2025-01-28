import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def create_cftr_boxplot_with_shading_and_visual_short_names(
    adata,
    gene="CFTR",
    sample_col="sample",
    cell_line_col="cell_line",
    threshold=0.1,
    colors=["#fff2cc", "#d9f0a3", "#a6bddb"],
    output_file=None
):
    """
    Filters CFTR-positive cells and creates a boxplot combined with a scatter plot of gene expression for each sample.
    Adds background shading corresponding to cell_line groups and shortens sample names visually.

    Parameters:
    - adata (AnnData): The input AnnData object.
    - gene (str): The gene to filter and plot (default is "CFTR").
    - sample_col (str): The column in adata.obs representing samples.
    - cell_line_col (str): The column in adata.obs representing cell lines for shading.
    - threshold (float): The minimum expression value to consider a cell CFTR-positive.
    - colors (list): Colors for shading different cell_line values.
    - output_file (str): Path to save the plot. If None, the plot will just be shown.
    """
    # Check if the gene exists in adata
    if gene not in adata.var_names:
        raise ValueError(f"The gene '{gene}' is not found in adata.var_names.")

    # Extract gene expression from raw counts
    if "counts" not in adata.layers:
        raise ValueError("The layer 'counts' is not found in adata.layers.")
        
    gene_expression = adata.layers["counts"][:, adata.var_names.get_loc(gene)]
    if hasattr(gene_expression, "toarray"):  # Handle sparse matrix
        gene_expression = gene_expression.toarray()
    gene_expression = gene_expression.flatten()  # Ensure 1D array

    # Add gene expression to adata.obs
    adata.obs[gene] = gene_expression

    # Filter to only CFTR-positive cells
    adata_cftr_positive = adata[adata.obs[gene] > threshold]

    # Create a DataFrame for plotting
    df = pd.DataFrame({
        "Sample": adata_cftr_positive.obs[sample_col],
        "Cell Line": adata_cftr_positive.obs[cell_line_col],
        gene: adata_cftr_positive.obs[gene]
    })

    # Convert the 'Cell Line' column to a categorical type with a specified order
    df["Cell Line"] = pd.Categorical(df["Cell Line"], categories=["CF1", "CF4", "CF5"], ordered=True)

    # Define the desired order for samples
    df['Sample'] = pd.Categorical(df['Sample'], ordered=True)

    # Sort the DataFrame by cell line and sample
    df.sort_values(["Cell Line", "Sample"], inplace=True)

    # Create the plot
    plt.figure(figsize=(10, 6))
    ax = plt.gca()

    # Add background shading for each unique cell line
    unique_cell_lines = df["Cell Line"].unique()
    for i, cell_line in enumerate(unique_cell_lines):
        samples_in_cell_line = df[df["Cell Line"] == cell_line]["Sample"].unique()
        sample_indices = [df["Sample"].unique().tolist().index(s) for s in samples_in_cell_line]
        if sample_indices:
            ax.axvspan(
                min(sample_indices) - 0.5,
                max(sample_indices) + 0.5,
                color=colors[i % len(colors)],
                alpha=0.3,
                label=cell_line if i < len(colors) else None
            )

    # Plot boxplot and scatterplot
    sns.boxplot(
        data=df, x="Sample", y=gene, color="white", width=0.6, linewidth=1.5, showfliers=False
    )
    sns.stripplot(
        data=df, x="Sample", y=gene, color="black", alpha=0.6, jitter=True, size=3
    )

    # Customize plot appearance
    plt.title(f"{gene} Expression", fontsize=20)
    plt.xlabel("Sample", fontsize=18)
    plt.ylabel(f"{gene} Expression", fontsize=18)
    plt.xticks(rotation=45, fontsize=16)
    plt.legend(title="Cell Line", loc="upper center", fontsize=16, frameon=True)
    plt.tight_layout()

    # Save or show the plot
    if output_file:
        plt.savefig(output_file, dpi=900, bbox_inches="tight")
        print(f"Plot saved to {output_file}")
    else:
        plt.show()
    return(df)

# Example usage:
adata_full = sc.read_h5ad("path_to_fulladata.h5ad")

# Create the plot
df = create_cftr_boxplot_with_shading_and_visual_short_names(
    adata_full,
    gene="CFTR",
    sample_col="sample",
    cell_line_col="cell_line",
    threshold=0.1,
    output_file="5C.png"
)
