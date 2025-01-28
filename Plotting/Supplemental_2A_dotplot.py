import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
import os

def scale_array_to_01(array):
    array = np.asarray(array)  # Ensure the input is a NumPy array
    scaled_array = (array - np.min(array)) / (np.max(array) - np.min(array))
    return scaled_array

def round_list_values(input_list):
    result = []
    for item in input_list:
        if item in ['Mean Expression', 'Percent Expressing']:
            result.append(item)  # Keep these strings as-is
        else:
            try:
                rounded_value = f"{float(item):.2g}"
                result.append(rounded_value)
            except ValueError:
                result.append(item)
    return result

def create_stacked_dotplot(adata_list, genes, titles, save_path, figsize=(14, 8), cmap="viridis"):
    os.makedirs(os.path.dirname(save_path), exist_ok=True)

    results = []
    for adata, gene, title in zip(adata_list, genes, titles):
        if gene not in adata.var_names:
            raise ValueError(f"Gene '{gene}' not found in adata.var_names.")

        gene_expression = adata[:, gene].X
        if hasattr(gene_expression, "toarray"):  # Handle sparse matrix
            gene_expression = gene_expression.toarray().flatten()
        scaled_expression = scale_array_to_01(gene_expression)

        adata.obs[gene] = scaled_expression

        grouped = adata.obs.groupby(['cell_type', 'treatment'])

        mean_expr = adata.obs.groupby(['cell_type', 'treatment'])[gene].apply(lambda x: x.mean())
        pct_expr = grouped[gene].apply(lambda x: (x > 0.0).mean() * 100)

        df = pd.DataFrame({
            'Mean Expression': mean_expr,
            'Percent Expressing': pct_expr,
            'Cell Type': grouped.size().index.get_level_values(0),
            'Treatment': grouped.size().index.get_level_values(1),
            'Title': title
        })
        results.append(df)

    combined_df = pd.concat(results).reset_index(drop=True)

    desired_order = ["DMSO", "LI", "ETI"]
    combined_df["Treatment"] = pd.Categorical(combined_df["Treatment"], categories=desired_order, ordered=True)
    combined_df = combined_df.sort_values(by=["Cell Type", "Treatment"])

    num_subplots = len(combined_df['Title'].unique())
    fig, axes = plt.subplots(num_subplots, 1, figsize=(figsize[0], figsize[1] * num_subplots))

    if num_subplots == 1:
        axes = [axes]

    for ax, title in zip(axes, combined_df['Title'].unique()):
        subset = combined_df[combined_df['Title'] == title]
        scatter = sns.scatterplot(
            data=subset,
            x='Treatment',
            y='Cell Type',
            size='Percent Expressing',
            hue='Mean Expression',
            palette=cmap,
            sizes=(60, 600),
            edgecolor="k",
            alpha=0.8,
            ax=ax,
            legend=True
        )
        ax.set_title(title, fontsize=16)
        ax.set_xlabel("", fontsize=12)
        ax.set_ylabel("Cell Type", fontsize=12)
        ax.set_xlim(-0.5, len(subset['Treatment'].unique()) - 0.5)

        handles, labels = ax.get_legend_handles_labels()
        for i, handle in enumerate(handles[:len(handles) // 2]):
            handle.set_markersize(handle.get_markersize() * 3)

        labels = round_list_values(labels)
        ax.legend(handles=handles, labels=labels, loc="upper left", bbox_to_anchor=(1, 1), title="Legend")

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.2)
    plt.savefig(save_path, dpi=900, bbox_inches="tight")
    print(f"Stacked dot plot saved to {save_path}")

plt.rcParams['font.size'] = 18
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18

# Load the data
adata = sc.read("path_to_adata.h5ad")

# Split the data based on cell line
adata_cf1 = adata[adata.obs['cell_line'] == 'CF1'].copy()
adata_cf4 = adata[adata.obs['cell_line'] == 'CF4'].copy()
adata_cf5 = adata[adata.obs['cell_line'] == 'CF5'].copy()

# Specify the genes to plot and titles
genes = ["CFTR","CFTR","CFTR"]
titles = ["CF1","CF4","CF5"]

# Define save path for the plot
save_path = "figures/stacked_dotplot.png"

# Generate and save the plot
create_stacked_dotplot([adata_cf1,adata_cf4, adata_cf5], genes, titles, save_path)
