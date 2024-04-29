import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import plotting_utils

def greens_function_absmax(exported_green_path, outpath, range_x, range_y, fs = 15):
    data = np.load(exported_green_path) # data comes in [r, z, t, field_dim], where field_dim = {E_r, E_z}

    E_abs_data = np.linalg.norm(data, axis = 3)
    E_abs_max_data = np.max(E_abs_data, axis = -1)

    figsize_y = 5
    figsize_x = (range_x[1] - range_x[0]) / (range_y[1] - range_y[0]) * figsize_y * 1.35
    fig = plt.figure(figsize = (figsize_x, figsize_y), layout = "constrained")
    ax = fig.add_subplot(111)
    
    fieldplot = plotting_utils.plot_field(ax, E_abs_max_data, range_x = range_x, range_y = range_y, fs = fs)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(fieldplot, cax = cax)
    cbar.set_label("$K_\mathrm{max}$ [a.u.]", fontsize = fs)
    cbar.ax.tick_params(labelsize = fs)
    
    fig.savefig(outpath)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--green_exported", action = "store", dest = "exported_green_path")
    parser.add_argument("--out", action = "store", dest = "outpath")
    parser.add_argument("--range_x", action = "store", nargs = "+", type = float, default = [0, 500.0 / 3.0])
    parser.add_argument("--range_y", action = "store", nargs = "+", type = float, default = [-300.0 / 3.0, 200.0 / 3.0])
    args = vars(parser.parse_args())

    greens_function_absmax(**args)

