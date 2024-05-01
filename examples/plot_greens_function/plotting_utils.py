import numpy as np
import matplotlib.colors as colors

def plot_field(ax, data, range_x, range_y, xlabel = "r [m]", ylabel = "z [m]", epilog = None, fs = 13, vmin = 0.0, vmax = 2e-3, linthresh = 1e-5, cmap = "bwr"):
    vals_x = np.linspace(range_x[0], range_x[1], data.shape[0])
    vals_y = np.linspace(range_y[0], range_y[1], data.shape[1])
    vals_mesh_x, vals_mesh_y = np.meshgrid(vals_x, vals_y)

    fieldplot = ax.pcolormesh(vals_mesh_x, vals_mesh_y, np.transpose(data), 
                              cmap = cmap, shading = "gouraud",
                              norm = colors.SymLogNorm(linthresh = linthresh, linscale = 1, base = 10, vmin = vmin, vmax = vmax))
    
    ax.set_aspect("equal")
    ax.set_xlabel(xlabel, fontsize = fs)
    ax.set_ylabel(ylabel, fontsize = fs)

    if epilog:
        epilog(ax)
    
    ax.tick_params(axis = "y", direction = "in", left = True, right = True, labelsize = fs)
    ax.tick_params(axis = "x", direction = "in", bottom = True, top = True, labelsize = fs)

    return fieldplot

def autoscale_y(ax, margin = 0.1):

    def get_visible_ylim(line):
        xvals, yvals = line.get_xdata(), line.get_ydata()
        x_low, x_high = ax.get_xlim()
        y_vis = yvals[np.logical_and(xvals > x_low, xvals < x_high)]
        y_min_vis, y_max_vis = np.min(y_vis), np.max(y_vis)
        y_range = y_max_vis - y_min_vis

        ylim = [y_min_vis - y_range * margin, y_max_vis + y_range * margin]
        return ylim

    lines = ax.get_lines()
    ylim = [np.inf, -np.inf]

    for line in lines:
        cur_ylim = get_visible_ylim(line)

        if cur_ylim[0] < ylim[0]:
            ylim[0] = cur_ylim[0]

        if cur_ylim[1] > ylim[1]:
            ylim[1] = cur_ylim[1]

    ax.set_ylim(*ylim)
