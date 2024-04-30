import argparse, math, string, yaml
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable

import plotting_utils

def pos_to_index(pos_xy, shape_xy, range_x, range_y):

    if not isinstance(pos_xy, np.ndarray):
        pos_xy = np.array(pos_xy)

    if not isinstance(shape_xy, np.ndarray):
        shape_xy = np.array(shape_xy)

    start_pos_xy = np.array([range_x[0], range_y[0]])
    end_pos_xy = np.array([range_x[1], range_y[1]])

    return (pos_xy - start_pos_xy) / (end_pos_xy - start_pos_xy) * shape_xy

def greens_function_absmax_annotated(exported_green_path, outpath, config_path, fs = 15):

    with open(config_path) as configfile:
        config = yaml.safe_load(configfile)
    
    # Make sure we have a list of 2-dim coordinates
    annotations = config["annotations"]
    number_annotations = len(annotations)

    data = np.load(exported_green_path) # data comes in [r, z, t, field_dim], where field_dim = {E_r, E_z}
    data *= float(config["field_scale_fact"])
    E_r_data = data[:,:,:,0]
    E_z_data = data[:,:,:,1]
    E_abs_data = np.linalg.norm(data, axis = 3)
    E_abs_max_data = np.max(E_abs_data, axis = -1)
    
    fig = plt.figure(figsize = (12, 4), layout = "constrained")
    gs = GridSpec(1, 2, figure = fig)

    # GridSpec for time=domain graphs
    gs_detail = GridSpecFromSubplotSpec(config["annotation_cols"], config["annotation_rows"], subplot_spec = gs[0])

    annotation_labels = string.ascii_uppercase
    def field_annotation_epilog(ax):
        for cur_annotation, cur_label in zip(annotations, annotation_labels):
            ax.scatter(*cur_annotation["pos"], color = "white", marker = "*")
            ax.text(cur_annotation["pos"][0] + 4, cur_annotation["pos"][1], cur_label, color = "white", fontsize = fs, ha = "left", va = "center")
            
    # Field plot
    ax_field = fig.add_subplot(gs[1])
    fieldplot = plotting_utils.plot_field(ax_field, E_abs_max_data, range_x = config["range_x"], range_y = config["range_y"], epilog = field_annotation_epilog, fs = fs)

    divider = make_axes_locatable(ax_field)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(fieldplot, cax = cax)
    cbar.set_label("$|\mathbf{K}|_\mathrm{max}$ [a.u.]", fontsize = fs)
    cbar.ax.tick_params(labelsize = fs)

    # Annotations
    tvals = np.linspace(*config["range_t"], E_abs_data.shape[2])
    for ind, (cur_annotation, cur_title) in enumerate(zip(annotations, annotation_labels)):
        ax_annotation = fig.add_subplot(gs_detail[ind])
        cur_detail_pos_ind = pos_to_index(cur_annotation["pos"], E_abs_max_data.shape, config["range_x"], config["range_y"])
        ax_annotation.plot(tvals, E_z_data[int(cur_detail_pos_ind[0]), int(cur_detail_pos_ind[1]), :], color = "black", label = r"$K_z$")
        ax_annotation.plot(tvals, E_r_data[int(cur_detail_pos_ind[0]), int(cur_detail_pos_ind[1]), :], color = "gray", label = r"$K_r$")
        ax_annotation.set_xlim(*cur_annotation["range_t"])
        ax_annotation.set_xlabel("Time [ns]", fontsize = fs)
        # ax_annotation.set_title(cur_title, fontsize = fs)
        ax_annotation.text(cur_title, fontsize = fs)
        ax_annotation.tick_params(axis = "y", direction = "in", left = True, right = True, labelsize = fs)
        ax_annotation.tick_params(axis = "x", direction = "in", bottom = True, top = True, labelsize = fs)    
        ax_annotation.legend(frameon = False)

    ax.set_xlim(*config["lim_x"])
    ax.set_ylim(*config["lim_y"])
        
    fig.savefig(outpath)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--green_exported", action = "store", dest = "exported_green_path")
    parser.add_argument("--out", action = "store", dest = "outpath")
    parser.add_argument("--config", action = "store", dest = "config_path")
    args = vars(parser.parse_args())
    
    greens_function_absmax_annotated(**args)
