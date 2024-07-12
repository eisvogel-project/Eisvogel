import argparse, yaml, os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import plotting_utils

def make_field_plot(outpath, field_data, epilog, range_x, range_y, xlim, ylim, zlabel, fs = 13):

    figsize_y = 5
    figsize_x = (xlim[1] - xlim[0]) / (ylim[1] - ylim[0]) * figsize_y * 1.2
    fig = plt.figure(figsize = (figsize_x, figsize_y), layout = "constrained")
    ax = fig.add_subplot(111)    

    fieldplot = plotting_utils.plot_field(ax, field_data, range_x = range_x, range_y = range_y, fs = fs, epilog = epilog,
                                          vmin = -2e-3, vmax = 2e-3, linthresh = 1e-8)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(fieldplot, cax = cax)
    cbar.set_label(zlabel, fontsize = fs)
    cbar.ax.tick_params(labelsize = fs)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    
    fig.savefig(outpath, dpi = 300)

def greens_function_still(exported_green_path, tval, outdir, config_path):

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    with open(config_path) as configfile:
        config = yaml.safe_load(configfile)
    
    data = np.load(exported_green_path, mmap_mode = 'r') # data comes in [r, z, t, field_dim], where field_dim = {E_r, E_z}
    E_r_data = data[:,:,:,0]
    E_z_data = data[:,:,:,1]

    t_ind = int(data.shape[2] / (config["range_t"][1] - config["range_t"][0]) * tval)
    
    def show_ice_surface(ax):
        ax.axhline(0.0, color = "gray", ls = "dashed")
    
    outpath = os.path.join(outdir, "E_r.png")
    make_field_plot(outpath, E_r_data[:, :, t_ind], show_ice_surface, config["range_x"], config["range_y"], config["lim_x"], config["lim_y"], zlabel = r"$K_r$ [a.u.]")

    outpath = os.path.join(outdir, "E_z.png")
    make_field_plot(outpath, E_z_data[:, :, t_ind], show_ice_surface, config["range_x"], config["range_y"], config["lim_x"], config["lim_y"], zlabel = r"$K_z$ [a.u.]")
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--green_exported", action = "store", dest = "exported_green_path")
    parser.add_argument("--tval", action = "store", dest = "tval", type = float)
    parser.add_argument("--outdir", action = "store", dest = "outdir")
    parser.add_argument("--config", action = "store", dest = "config_path")
    args = vars(parser.parse_args())

    greens_function_still(**args)
