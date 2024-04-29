import argparse
import numpy as np
import matplotlib.pyplot as plt

def plot_field(ax, data, range_x, range_y, xlabel = "r [m]", ylabel = "z [m]", epilog = None, fs = 13):
    vals_x = np.linspace(range_x[0], range_x[1], data.shape[0])
    vals_y = np.linspace(range_y[0], range_y[1], data.shape[1])
    vals_mesh_x, vals_mesh_y = np.meshgrid(vals_x, vals_y)
    
    linthresh = 1e-4
    vmin = 0.0
    vmax = 2e-3
    
    fieldplot = ax.pcolormesh(vals_mesh_x, vals_mesh_y, np.transpose(data), 
                              cmap = "coolwarm", shading = "gouraud",
                              norm = colors.SymLogNorm(linthresh = linthresh, linscale = 1, base = 10, vmin = vmin, vmax = vmax))
    
    ax.set_aspect("equal")
    ax.set_xlabel(xlabel, fontsize = fs)
    ax.set_ylabel(ylabel, fontsize = fs)

    if epilog:
        epilog(ax)
    
    ax.tick_params(axis = "y", direction = "in", left = True, right = True, labelsize = fs)
    ax.tick_params(axis = "x", direction = "in", bottom = True, top = True, labelsize = fs)

    return fieldplot

def greens_function_absmax(exported_green_path, outpath):
    data = np.load(exported_green_path) # data comes in [r, z, t, field_dim], where field_dim = {E_r, E_z}

    E_abs_data = np.linalg.norm(data, axis = 3)
    E_abs_max_data = np.max(E_abs_data, axis = -1)

    

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--green_exported", action = "store", dest = "exported_green_path")
    parser.add_argument("--out", action = "store", dest = "outpath")
    args = vars(parser.parse_args())

    greens_function_absmax(**args)

