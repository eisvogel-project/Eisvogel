import argparse, json, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def visualize_tree(treepath, vizpath):

    with open(treepath, 'r') as treefile:
        tree = json.load(treefile)

    fig, ax = plt.subplots()

    plot_start_coords = None
    plot_end_coords = None

    colormap = plt.cm.nipy_spectral
    colors = colormap(np.linspace(0, 1, 10))
    
    for entry in tree:

        entry_type = entry["type"]
        
        start_coords = np.array(entry["start_coords"])
        end_coords = np.array(entry["end_coords"])
        shape = np.array(entry["shape"])
            
        assert all(end_coords - start_coords == shape)

        if entry_type == "Node":
            tree_level = entry["level"]
            margin = (tree_level + 1) * 0.3
            edgecolor = colors[tree_level]

            # if tree_level != 1 and tree_level != 0:
            #     continue

            # if tree_level != 2 and tree_level != 1 and tree_:
            #     continue
            
        else:            
            margin = 0
            edgecolor = "gray"

            continue;

        start_coords_plot = start_coords - [margin, margin]
        end_coords_plot = end_coords + [margin, margin]
        shape_plot = end_coords_plot - start_coords_plot

        if plot_start_coords is None:
            plot_start_coords = np.copy(start_coords_plot)

        if plot_end_coords is None:
            plot_end_coords = np.copy(end_coords_plot)

        plot_start_coords = np.minimum(plot_start_coords, start_coords_plot)
        plot_end_coords = np.maximum(plot_end_coords, end_coords_plot)
            
        rect = patches.Rectangle(start_coords_plot, shape_plot[0], shape_plot[1],
                                 linewidth = tree_level / 2, edgecolor = edgecolor, facecolor = 'none',
                                 zorder = tree_level)
        ax.add_patch(rect)    

    plot_margin = 1
    ax.set_xlim(plot_start_coords[0] - plot_margin, plot_end_coords[0] + plot_margin)
    ax.set_ylim(plot_start_coords[1] - plot_margin, plot_end_coords[1] + plot_margin)

    # ax.set_xlim(2000, 3000)
    # ax.set_ylim(2000, 3000)
    
    fig.savefig(vizpath)
        
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", action = "store", dest = "treepath")
    parser.add_argument("--viz", action = "store", dest = "vizpath")
    args = vars(parser.parse_args())

    visualize_tree(**args)
