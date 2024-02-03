import argparse, os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from eisvogel import CylindricalWeightingField, CoordVector

def plot_2d(outpath, vals_xy, vals_z, xlabel = "", ylabel = "", zlabel = "", title = "", fs = 13,
            num_pts_interp = (1000, 1000), cmap = "bwr"):

    vals_x = [cur[0] for cur in vals_xy]
    vals_y = [cur[1] for cur in vals_xy]
    
    vals_fine_x = np.linspace(min(vals_x), max(vals_x), num_pts_interp[0])
    vals_fine_y = np.linspace(min(vals_y), max(vals_y), num_pts_interp[1])
    
    triang = tri.Triangulation(vals_x, vals_y)
    interpolator = tri.LinearTriInterpolator(triang, vals_z)
    vals_mesh_x, vals_mesh_y = np.meshgrid(vals_fine_x, vals_fine_y)
    vals_interp_z = interpolator(vals_mesh_x, vals_mesh_y)

    fig, ax = plt.subplots(1, 1, figsize = (5, 5))
    surveyplot = ax.pcolormesh(vals_mesh_x, vals_mesh_y, vals_interp_z, cmap = cmap,
                               norm = matplotlib.colors.SymLogNorm(linthresh = 1e-6))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(surveyplot, cax = cax)    
    cbar.set_label(zlabel, fontsize = fs)
    
    cbar.ax.tick_params(labelsize = fs)
    ax.tick_params(axis = "y", direction = "in", left = True, right = True, labelsize = fs)
    ax.tick_params(axis = "x", direction = "in", bottom = True, top = True, labelsize = fs)

    ax.set_xlabel(xlabel, fontsize = fs)
    ax.set_ylabel(ylabel, fontsize = fs)
    ax.set_title(title, fontsize = fs)

    plt.tight_layout()
    
    fig.savefig(outpath, dpi = 300)
    plt.close()

def plot_fields(outdir, cwf, tval, num_pts = 100):

    start_coords_txyz = cwf.GetStartCoords()
    end_coords_txyz = cwf.GetEndCoords()
    
    vals_xz = []
    vals_E_r = []
    vals_E_z = []
    vals_E_abs = []
    for cur_xval in np.linspace(start_coords_txyz[1] + 2, end_coords_txyz[1] - 2, num_pts):
        for cur_zval in np.linspace(start_coords_txyz[3] + 2, end_coords_txyz[3] - 2, num_pts):
            
            vals_xz.append([cur_xval, cur_zval])
            cur_yval = 0.0
            
            (E_r, E_z) = cwf.E_rz([tval, cur_xval, cur_yval, cur_zval])
            vals_E_r.append(E_r)
            vals_E_z.append(E_z)
            vals_E_abs.append(np.sqrt(E_r**2 + E_z**2))

    outpath = os.path.join(outdir, f"E_r_{tval}.png")
    plot_2d(outpath, vals_xz, vals_E_r, xlabel = "x", ylabel = "z", zlabel = r"$E_r$", title = rf"t = {tval}", cmap = "bwr")

    outpath = os.path.join(outdir, f"E_z_{tval}.png")
    plot_2d(outpath, vals_xz, vals_E_z, xlabel = "x", ylabel = "z", zlabel = r"$E_z$", title = rf"t = {tval}", cmap = "bwr")

    outpath = os.path.join(outdir, f"E_abs_{tval}.png")
    plot_2d(outpath, vals_xz, vals_E_abs, xlabel = "x", ylabel = "z", zlabel = "|E|", title = rf"t = {tval}", cmap = "Blues")
    
def plot_weighting_field(wf_path, outdir):
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    cwf = CylindricalWeightingField(wf_path)
    start_coords_txyz = cwf.GetStartCoords()
    end_coords_txyz = cwf.GetEndCoords()
    
    num_time_slices = 5
    
    for tval in np.linspace(start_coords_txyz[0], end_coords_txyz[0] - 2, num_time_slices):
        print(f"Plotting for t = {tval}")
        plot_fields(outdir, cwf, tval)
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--wf_path", action = "store", dest = "wf_path")
    parser.add_argument("--outdir", action = "store", dest = "outdir")    
    args = vars(parser.parse_args())

    plot_weighting_field(**args)
