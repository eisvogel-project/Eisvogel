import argparse, yaml
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import plotting_utils

def greens_function_absmax(exported_green_path, outpath, config_path, fs = 15, normalized = True):

    with open(config_path) as configfile:
        config = yaml.safe_load(configfile)
    
    data = np.load(exported_green_path) # data comes in [r, z, t, field_dim], where field_dim = {E_r, E_z}
    E_abs_data = np.linalg.norm(data, axis = 3)
    E_abs_max_data = np.max(E_abs_data, axis = -1)

    if normalized:
        normval = np.max(E_abs_max_data)
        E_abs_max_data /= normval

    figsize_y = 5
    figsize_x = (config["range_x"][1] - config["range_x"][0]) / (config["range_y"][1] - config["range_y"][0]) * figsize_y * 1.2
    fig = plt.figure(figsize = (figsize_x, figsize_y), layout = "constrained")
    ax = fig.add_subplot(111)    

    def show_ice_surface(ax):
        ax.axhline(0.0, color = "gray", ls = "dashed")
    
    fieldplot = plotting_utils.plot_field(ax, E_abs_max_data, range_x = config["range_x"], range_y = config["range_y"], fs = fs, epilog = show_ice_surface, cmap = "coolwarm")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(fieldplot, cax = cax)
    cbar.set_label("$|\mathbf{K}|_\mathrm{max}$ [a.u.]", fontsize = fs)
    cbar.ax.tick_params(labelsize = fs)

    ax.set_xlim(*config["lim_x"])
    ax.set_ylim(*config["lim_y"])
    
    fig.savefig(outpath, dpi = 300)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--green_exported", action = "store", dest = "exported_green_path")
    parser.add_argument("--out", action = "store", dest = "outpath")
    parser.add_argument("--config", action = "store", dest = "config_path")
    args = vars(parser.parse_args())

    greens_function_absmax(**args)

