#!/usr/bin/env python3
import argparse
import os
import torch
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D # For custom legend

# Define mappings for Channel 0 (Bases) visualization
# These should match what was used for generating the tensor if possible
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4, '*': 5, ' ': 6, '-': 6}
INDEX_TO_BASE_VIS = {v: k for k, v in BASE_TO_INDEX.items()}
# Ensure all indices up to the max are covered for colormap generation
MAX_BASE_INDEX = max(BASE_TO_INDEX.values())

def visualize_channel(tensor_filepath, channel_index, output_image_path=None, custom_title=None, cmap_name=None):
    """
    Loads a 3-channel tensor from a .pth file and visualizes the specified channel as an image.

    Args:
        tensor_filepath (str): Path to the .pth tensor file.
        channel_index (int): Index of the channel to visualize (0, 1, or 2).
        output_image_path (str, optional): Path to save the output image. If None, displays the image.
        custom_title (str, optional): Custom title for the plot.
        cmap_name (str, optional): Specific matplotlib colormap to use.
    """
    if not os.path.exists(tensor_filepath):
        print(f"❌ Error: Tensor file not found at {tensor_filepath}")
        return

    try:
        # Load the tensor (assuming it's saved on CPU or can be moved to CPU)
        tensor = torch.load(tensor_filepath, map_location=torch.device('cpu'))
    except Exception as e:
        print(f"❌ Error loading tensor from {tensor_filepath}: {e}")
        return

    if not isinstance(tensor, torch.Tensor):
        print(f"❌ Error: Loaded file {tensor_filepath} does not contain a PyTorch tensor.")
        return

    if tensor.ndim != 3 or tensor.shape[0] != 3:
        print(f"❌ Error: Tensor is not 3-channel. Shape found: {tensor.shape}")
        return

    if not (0 <= channel_index <= 2):
        print(f"❌ Error: Channel index must be 0, 1, or 2. Got {channel_index}.")
        return

    # Select the specified channel and convert to NumPy array
    channel_data = tensor[channel_index].numpy().astype(float) # imshow prefers float

    num_rows, window_size = channel_data.shape

    fig, ax = plt.subplots(figsize=(max(10, window_size / 10), max(6, num_rows / 20)))

    # --- Colormap and normalization setup ---
    cmap = cmap_name
    norm = None
    cbar = None
    legend_elements = None

    if channel_index == 0:  # Bases
        channel_name = "Bases (Channel 0)"
        # Define a discrete colormap for bases
        # Ensure colors cover all possible indices from 0 to MAX_BASE_INDEX
        base_colors_list = ['lightcoral', 'lightskyblue', 'lightgreen', 'gold', 'plum', 'grey', 'whitesmoke']
        # Extend or wrap colors if MAX_BASE_INDEX + 1 > len(base_colors_list)
        if MAX_BASE_INDEX + 1 > len(base_colors_list):
             # Simple repetition if more indices than colors
            extended_colors = (base_colors_list * ((MAX_BASE_INDEX + 1) // len(base_colors_list) + 1))[:MAX_BASE_INDEX+1]
        else:
            extended_colors = base_colors_list[:MAX_BASE_INDEX+1]

        cmap_bases = mcolors.ListedColormap(extended_colors)
        norm = mcolors.BoundaryNorm(np.arange(-0.5, MAX_BASE_INDEX + 1.5, 1), cmap_bases.N)
        cmap = cmap_bases

        # Prepare legend handles for bases
        legend_elements = []
        unique_indices_present = np.unique(channel_data)
        for i in range(MAX_BASE_INDEX + 1):
            if i in unique_indices_present: # Only add legend items for present bases
                base_char = INDEX_TO_BASE_VIS.get(i, f"Idx {i}")
                legend_elements.append(Line2D([0], [0], marker='s', color='w',
                                              label=f'{base_char} ({i})',
                                              markerfacecolor=cmap(norm(i)), markersize=10))

    elif channel_index == 1:  # Qualities
        channel_name = "Qualities (Channel 1)"
        if cmap is None:
            cmap = 'viridis' # Good for sequential data like qualities
        # Qualities usually range from 0 to ~40 (or up to 93 for Phred).
        # You might want to adjust vmin/vmax if you know the typical range.
        im = ax.imshow(channel_data, cmap=cmap, interpolation='nearest', aspect='auto')
        cbar = fig.colorbar(im, ax=ax, orientation='vertical', label='Quality Score (Phred-like)')

    elif channel_index == 2:  # Mismatches
        channel_name = "Mismatches (Channel 2)"
        # Binary data: 0 for match/padding, 1 for mismatch
        if cmap is None:
            cmap = mcolors.ListedColormap(['lightgray', 'red']) # 0: lightgray (match/pad), 1: red (mismatch)
        norm = mcolors.BoundaryNorm([-0.5, 0.5, 1.5], cmap.N)
        im = ax.imshow(channel_data, cmap=cmap, norm=norm, interpolation='nearest', aspect='auto')
        # Custom colorbar for mismatches
        cbar = fig.colorbar(im, ax=ax, ticks=[0, 1], orientation='vertical')
        cbar.ax.set_yticklabels(['Match/Pad (0)', 'Mismatch (1)'])

    if channel_index == 0: # Bases need specific imshow for discrete colormap
        im = ax.imshow(channel_data, cmap=cmap, norm=norm, interpolation='nearest', aspect='auto')
        if legend_elements:
             # Place legend outside the plot to avoid overlap
            ax.legend(handles=legend_elements, title="Base Legend", bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.subplots_adjust(right=0.75) # Adjust layout to make space for legend
    elif channel_index != 1 and channel_index != 2: # If not already handled by specific cbar logic
        im = ax.imshow(channel_data, cmap=cmap, interpolation='nearest', aspect='auto')
        fig.colorbar(im, ax=ax, orientation='vertical')


    # --- Titles and labels ---
    plot_title = custom_title if custom_title else f"Tensor Visualization: {os.path.basename(tensor_filepath)}\n{channel_name}"
    ax.set_title(plot_title, fontsize=14)
    ax.set_xlabel("Window Position", fontsize=12)
    ax.set_ylabel(f"Read / Ref Index (0 = Ref, Total Rows: {num_rows})", fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=10)

    # Adjust y-ticks to be more readable if many rows
    if num_rows > 20:
        step = num_rows // 10
        ax.set_yticks(np.arange(0, num_rows, step))
    else:
        ax.set_yticks(np.arange(0, num_rows, 1))


    plt.tight_layout() # Adjust layout to prevent labels from overlapping (if legend is not too wide)
    if legend_elements: # If legend was created, tight_layout might need to be called after legend placement
        pass # Already called plt.subplots_adjust for legend


    if output_image_path:
        try:
            plt.savefig(output_image_path, dpi=150, bbox_inches='tight')
            print(f"🖼️ Image saved to {output_image_path}")
        except Exception as e:
            print(f"❌ Error saving image to {output_image_path}: {e}")
    else:
        plt.show()

def main():
    parser = argparse.ArgumentParser(description="Visualize a specific channel of a 3-channel tensor from a .pth file.")
    parser.add_argument("pth_file", help="Path to the .pth tensor file.")
    parser.add_argument("channel", type=int, choices=[0, 1, 2],
                        help="Channel index to visualize (0: Bases, 1: Qualities, 2: Mismatches).")
    parser.add_argument("-o", "--output_image", help="Optional: Path to save the visualization as an image file.")
    parser.add_argument("-t", "--title", help="Optional: Custom title for the plot.")
    parser.add_argument("-c", "--cmap", help="Optional: Matplotlib colormap name to use (e.g., 'viridis', 'gray'). Overrides defaults.")

    args = parser.parse_args()

    visualize_channel(args.pth_file, args.channel, args.output_image, args.title, args.cmap)

if __name__ == "__main__":
    main()