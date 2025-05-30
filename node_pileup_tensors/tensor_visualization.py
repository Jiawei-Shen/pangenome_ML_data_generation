#!/usr/bin/env python3
import argparse
import os
import torch
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D  # For custom legend

# Define mappings for Channel 0 (Bases) visualization
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4, '*': 5, ' ': 6, '-': 6}
INDEX_TO_BASE_VIS = {v: k for k, v in BASE_TO_INDEX.items()}
MAX_BASE_INDEX = max(BASE_TO_INDEX.values())

def plot_on_ax(ax, channel_data, channel_index, tensor_basename, num_rows, window_size, custom_subplot_title=None,
               user_cmap_name=None):
    cmap = user_cmap_name
    norm = None
    legend_elements = None
    current_channel_name = f"Channel {channel_index}"  # Default name

    if channel_index == 0:  # Bases
        current_channel_name = "Bases"
        base_colors_list = ['lightcoral', 'lightskyblue', 'lightgreen', 'gold', 'plum', 'grey', 'whitesmoke']
        extended_colors = (base_colors_list * ((MAX_BASE_INDEX + 1) // len(base_colors_list) + 1))[:MAX_BASE_INDEX + 1]
        cmap_bases = mcolors.ListedColormap(extended_colors)
        norm = mcolors.BoundaryNorm(np.arange(-0.5, MAX_BASE_INDEX + 1.5, 1), cmap_bases.N)
        if cmap is None:
            cmap = cmap_bases

        legend_elements = []
        unique_indices = np.unique(channel_data)
        legend_cmap = plt.get_cmap(cmap) if isinstance(cmap, str) else cmap
        for i in range(MAX_BASE_INDEX + 1):
            if i in unique_indices:
                base_char = INDEX_TO_BASE_VIS.get(i, f"Idx {i}")
                legend_elements.append(
                    Line2D([0], [0], marker='s', color='w', label=f'{base_char} ({i})',
                           markerfacecolor=legend_cmap(norm(i)), markersize=6)
                )
        im = ax.imshow(channel_data, cmap=cmap, norm=norm, interpolation='nearest', aspect='auto')
        if legend_elements:
            ax.legend(handles=legend_elements, title="Base Legend", fontsize='small', title_fontsize='small',
                      bbox_to_anchor=(1.02, 1), loc='upper left')

    elif channel_index == 1:  # Qualities
        current_channel_name = "Qualities"
        if cmap is None:
            cmap = 'viridis'
        im = ax.imshow(channel_data, cmap=cmap, interpolation='nearest', aspect='auto')
        ax.figure.colorbar(im, ax=ax, orientation='vertical', label='Quality Score', fraction=0.046, pad=0.04)

    elif channel_index == 2:  # Mismatches
        current_channel_name = "Mismatches"
        default_cmap = mcolors.ListedColormap(['lightgray', 'red'])
        default_norm = mcolors.BoundaryNorm([-0.5, 0.5, 1.5], default_cmap.N)
        effective_cmap = cmap if cmap is not None else default_cmap
        effective_norm = norm if norm is not None else default_norm
        im = ax.imshow(channel_data, cmap=effective_cmap, norm=effective_norm,
                       interpolation='nearest', aspect='auto')
        cbar = ax.figure.colorbar(im, ax=ax, ticks=[0, 1], orientation='vertical', fraction=0.046, pad=0.04)
        cbar.ax.set_yticklabels(['Match/Pad (0)', 'Mismatch (1)'])

    else:
        im = ax.imshow(channel_data, cmap='gray', interpolation='nearest', aspect='auto')
        ax.figure.colorbar(im, ax=ax, orientation='vertical')

    title = custom_subplot_title if custom_subplot_title else f"{current_channel_name} from {tensor_basename}"
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Window Position", fontsize=9)
    ax.set_ylabel(f"Read/Ref (Rows: {num_rows})", fontsize=9)
    ax.tick_params(axis='both', which='major', labelsize=8)
    if num_rows > 10:
        step = max(1, num_rows // (5 if num_rows <= 40 else 10))
        ax.set_yticks(np.arange(0, num_rows, step))


def process_file(pth_path, output_dir, args):
    tensor = torch.load(pth_path, map_location='cpu')
    if not isinstance(tensor, torch.Tensor) or tensor.ndim != 3 or tensor.shape[0] != 3:
        print(f"Skipping {pth_path}: not a 3-channel tensor.")
        return

    basename = os.path.splitext(os.path.basename(pth_path))[0]
    out_name = f"{basename}.png"
    out_path = os.path.join(output_dir, out_name)

    if args.channel is None:
        first = tensor[0].numpy()
        rows, cols = first.shape
        fig, axes = plt.subplots(3, 1, figsize=(max(10, cols/8), 3*max(4, rows/10)+2), sharex=True)
        titles = ["Bases (0)", "Qualities (1)", "Mismatches (2)"]
        for i, ax in enumerate(axes):
            data = tensor[i].numpy().astype(float)
            plot_on_ax(ax, data, i, basename, data.shape[0], data.shape[1],
                       custom_subplot_title=titles[i], user_cmap_name=args.cmap)
        fig.suptitle(args.title or f"Tensor Viz: {basename}", fontsize=16)
        fig.tight_layout(rect=[0, 0.03, 0.90, 0.95])
    else:
        data = tensor[args.channel].numpy().astype(float)
        rows, cols = data.shape
        fig, ax = plt.subplots(figsize=(max(10, cols/8), max(6, rows/10)+1))
        plot_on_ax(ax, data, args.channel, basename, rows, cols,
                   custom_subplot_title=args.title, user_cmap_name=args.cmap)
        fig.tight_layout(rect=[0, 0, 0.90, 1])

    if args.output_image:
        out_path = os.path.join(output_dir, args.output_image)
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Batch visualize 3-channel .pth tensors in a directory.")
    parser.add_argument("input_dir", help="Directory containing .pth files.")
    parser.add_argument("output_dir", help="Directory to save generated images.")
    parser.add_argument("-ch", "--channel", type=int, choices=[0,1,2], default=None,
                        help="Single channel to plot (0-2). Default: all channels.")
    parser.add_argument("-o", "--output_image",
                        help="Override output image filename for each tensor.")
    parser.add_argument("-t", "--title", help="Custom title or suptitle for plots.")
    parser.add_argument("-c", "--cmap", help="Matplotlib colormap name to override defaults.")
    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        print(f"Error: Input directory not found: {args.input_dir}")
        return
    os.makedirs(args.output_dir, exist_ok=True)

    for fname in os.listdir(args.input_dir):
        if fname.endswith('.pth'):
            pth_path = os.path.join(args.input_dir, fname)
            process_file(pth_path, args.output_dir, args)

if __name__ == "__main__":
    main()
