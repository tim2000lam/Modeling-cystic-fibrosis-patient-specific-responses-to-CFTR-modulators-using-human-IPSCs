import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex

# Adjust text sizes
plt.rcParams['font.size'] = 18  # Adjust base font size
plt.rcParams['axes.titlesize'] = 18  # Title font size
plt.rcParams['axes.labelsize'] = 18  # Axis label font size
plt.rcParams['legend.fontsize'] = 18  # Legend font size

# Read in a single adata object
adata = sc.read("path_to_adata.h5ad")


# Split the adata into datasets based on `cell_line` values
adata_cf1 = adata[adata.obs['cell_line'] == 'CF1'].copy()
adata_cf4 = adata[adata.obs['cell_line'] == 'CF4'].copy()
adata_cf5 = adata[adata.obs['cell_line'] == 'CF5'].copy()

# Create list of unique cell types for colour assignment
adatas = [adata_cf1, adata_cf4, adata_cf5]
unique_cell_types_full = []
for adata in adatas:
    unique_cell_types = adata.obs['cell_type'].unique()
    for cell_type in unique_cell_types:
        if cell_type not in unique_cell_types_full:
            unique_cell_types_full.append(cell_type)

# Assign colour palette
colors = sns.color_palette("tab20", len(unique_cell_types_full))
cell_type_to_color = dict(zip(unique_cell_types_full, colors))

# Convert colors to a format compatible with Scanpy (hex values)
cell_type_to_hex_color = {ctype: to_hex(color) for ctype, color in cell_type_to_color.items()}

# Update colour assignment to align with cluster annotations
colour_assignment = {ctype: cell_type_to_hex_color[ctype] for ctype in unique_cell_types_full}

# Initialize a figure with subplots
fig, axes = plt.subplots(len(adatas), 1, figsize=(12, 18))  # Adjust size for visibility; one row per dataset

# Generate stacked bar plots for each adata
for i, (adata, ax) in enumerate(zip(adatas, axes)):
    # Group by 'treatment' and 'cell_type' to get the cell counts for each combination
    group_counts = adata.obs.groupby(['treatment', 'cell_type']).size().unstack(fill_value=0)

    # Calculate proportions by dividing each count by the total cells per treatment
    proportions = group_counts.div(group_counts.sum(axis=1), axis=0)

    # Keep only columns (cell types) that exist in `colour_assignment`
    valid_cell_types = [ctype for ctype in proportions.columns if ctype in colour_assignment]
    proportions = proportions[valid_cell_types]

    proportions = proportions.reindex(["DMSO", "LI", "ETI"])

    # Plotting with custom colors and ordered stacking
    proportions.plot(
        kind='bar',
        stacked=True,
        color=[colour_assignment[cell_type] for cell_type in proportions.columns],
        ax=ax,  # Plot to the specific subplot axis
    )
    ax.set_title(f'{["CF1", "CF4", "CF5"][i]}')
    ax.set_xlabel('Treatment')
    ax.set_ylabel('Proportion of Cells')
    ax.legend(
        title='Cell Type',
        bbox_to_anchor=(1.05, 1),  # Position the legend outside the plot
        loc='upper left',
        frameon=False,
        reverse=True
    )
    ax.tick_params(axis='x', labelrotation=0, labelsize=10)

# Adjust layout to ensure subplots fit nicely
plt.tight_layout()

# Save the combined plot
output_path = "combined_stacked_bar_plots_individual_legends.png"
plt.savefig(output_path, dpi=300)  # Save with high resolution
plt.show()

print(f"Stacked bar plots with individual legends saved as {output_path}")
