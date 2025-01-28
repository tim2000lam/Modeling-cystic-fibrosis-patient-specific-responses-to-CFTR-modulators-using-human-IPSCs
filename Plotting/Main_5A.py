import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex

# Adjust text sizes
plt.rcParams['font.size'] = 18  # Adjust base font size
plt.rcParams['axes.titlesize'] = 18  # Title font size
plt.rcParams['axes.labelsize'] = 18  # Axis label font size
plt.rcParams['legend.fontsize'] = 18  # Legend font size

# Read in adata
adata = sc.read("path_to_bc_adata.h5ad")

# Create list of unique cell types for colour assignment
unique_cell_types = adata.obs['cell_type'].unique()

# Assign colour palette
colors = sns.color_palette("tab20", len(unique_cell_types))
cell_type_to_color = dict(zip(unique_cell_types, colors))

# Convert colors to a format compatible with Scanpy (hex values)
cell_type_to_hex_color = {ctype: to_hex(color) for ctype, color in cell_type_to_color.items()}

# Update colour assignment to align with cluster annotations
colour_assignment = {ctype: cell_type_to_hex_color[ctype] for ctype in unique_cell_types}

# Set the color map in Scanpy for the cell_type annotations
adata.uns['cell_type_colors'] = [colour_assignment[ctype] for ctype in unique_cell_types]

# Define plot title and file name for saving
plot_title = "UMAP by Cell Type"
file_name = "cell_type_annotations_umap.png"

# Generate UMAP plot without showing it immediately
sc.pl.umap(
    adata,
    color='cell_type',
    palette=cell_type_to_color,
    legend_loc='none',  # Disable default legend placement
    title=plot_title,
    save=False,  # Do not save yet
    show=False  # Suppress immediate display
)

# Access the current figure
fig, ax = plt.gcf(), plt.gca()

# Adjust layout to make space for the legend
plt.subplots_adjust(bottom=0.3)  # Adjust bottom margin to fit the legend

# Create a custom legend below the plot
legend_handles = [
    plt.Line2D([0], [0], marker='o', color=color, linestyle='', markersize=8, label=ctype)
    for ctype, color in cell_type_to_color.items()
]
fig.legend(
    handles=legend_handles,
    loc='lower center',  # Position legend below the plot
    bbox_to_anchor=(0.5, -0.1),  # Center it horizontally below the plot
    ncol=5,  # Adjust the number of columns as needed
    frameon=False,  # Remove the box around the legend
)

# Save the updated figure
fig.savefig(file_name, bbox_inches='tight')

# Show the plot with the updated legend
plt.show()
