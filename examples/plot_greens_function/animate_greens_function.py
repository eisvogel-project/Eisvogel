import argparse, yaml, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors

def generate_animation(outpath, data, range_x, range_y, range_t, xlabel = "r [m]", ylabel = "z [m]", zlabel = "", epilog = None, fs = 13, downsample_t = 1):
    num_frames = data.shape[2]
    vals_x = np.linspace(range_x[0], range_x[1], data.shape[0])
    vals_y = np.linspace(range_y[0], range_y[1], data.shape[1])

    vals_mesh_x, vals_mesh_y = np.meshgrid(vals_x, vals_y)

    figsize_y = 5
    figsize_x = (range_x[1] - range_x[0]) / (range_y[1] - range_y[0]) * figsize_y * 1.35
    
    fig = plt.figure(figsize = (figsize_x, figsize_y), dpi = 400)
    ax = fig.add_subplot(111)

    fig.subplots_adjust(left = 0.15, right = 0.8, top = 0.95, bottom = 0.1)

    linthresh = 1e-8
    vmin = -2e-3
    vmax = 2e-3

    fieldplot = ax.pcolormesh(vals_mesh_x, vals_mesh_y, np.transpose(data[:, :, 0]), 
                            cmap = "bwr", shading = "gouraud",
                            norm = colors.SymLogNorm(linthresh = linthresh, linscale = 1, base = 10, vmin = vmin, vmax = vmax))
    ax.set_aspect("equal")
    ax.set_xlabel(xlabel, fontsize = fs)
    ax.set_ylabel(ylabel, fontsize = fs)
    
    ax.tick_params(axis = "y", direction = "in", left = True, right = True, labelsize = fs)
    ax.tick_params(axis = "x", direction = "in", bottom = True, top = True, labelsize = fs)

    timestep_text = ax.text(0.6, 0.85, r"$t = {:.2f}\,\mathrm{{ns}}$".format(range_t[0]),
                            transform = ax.transAxes, fontsize = fs)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(fieldplot, cax = cax)
    cbar.set_label(zlabel, fontsize = fs)
    cbar.ax.tick_params(labelsize = fs)

    if epilog:
        epilog(ax)        
        
    def animate(it):
        ind_t = downsample_t * it
        fieldplot.set_array(np.transpose(data[:, :, ind_t]).ravel())
        timestep_text.set_text(r"$t = {:.2f}\,\mathrm{{ns}}$".format(range_t[0] + ind_t * (range_t[1] - range_t[0]) / num_frames))
        print(f"Rendered frame {ind_t} / {num_frames}")
        return fieldplot
    
    anim = animation.FuncAnimation(fig, animate, frames = int(num_frames / downsample_t), interval = 5, repeat = False)

    writer = animation.PillowWriter(fps = 15,
                                    metadata = dict(artist='Me'),
                                    bitrate = 3000)
    anim.save(outpath, writer = writer)

    plt.close()

def animate_greens_function(exported_green_path, outdir, config_path):

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    with open(config_path) as configfile:
        config = yaml.safe_load(configfile)
    
    data = np.load(exported_green_path)
    
    E_r_data = data[:,:,:,0]
    E_z_data = data[:,:,:,1]
    E_abs_data = np.linalg.norm(data, axis = 3)

    def show_interface(ax):
        ax.axhline(0.0, color = "gray", ls = "dashed")

    outpath = os.path.join(outdir, "E_z.gif")
    generate_animation(outpath, E_z_data, range_x = config["range_x"], range_y = config["range_y"], range_t = config["range_t"],
                       zlabel = r"$E_z$ [a.u.]", epilog = show_interface)

    outpath = os.path.join(outdir, "E_r.gif")
    generate_animation(outpath, E_r_data, range_x = config["range_x"], range_y = config["range_y"], range_t = config["range_t"],
                       zlabel = r"$E_r$ [a.u.]", epilog = show_interface)
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--green_exported", action = "store", dest = "exported_green_path")
    parser.add_argument("--outdir", action = "store", dest = "outdir")
    parser.add_argument("--config", action = "store", dest = "config_path")
    args = vars(parser.parse_args())

    animate_greens_function(**args)
