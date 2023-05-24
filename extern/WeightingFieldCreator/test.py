import meep as mp
import numpy as np

cell = mp.Vector3(20, 20, 20)

geometry = [
    mp.Block(
        mp.Vector3(mp.inf, mp.inf, 20),
        center=mp.Vector3(0, 0, -10),
        material=mp.Medium(epsilon=15),
    )
]

delta_func = lambda t: np.exp(-50 * (t-10)**2)

sources = [
    mp.Source(
        mp.CustomSource(src_func = delta_func,
                        start_time = 1.0,
                        end_time = 20.0,
                        is_integrated = False),
        component = mp.Ez,
        center = mp.Vector3(0, 0, -1)
    )
]

pml_layers = [mp.PML(1.0)]

resolution = 10

sim = mp.Simulation(
    cell_size=cell,
    boundary_layers=pml_layers,
    geometry=geometry,
    sources=sources,
    resolution=resolution,
)

sim.run(until=20)

import matplotlib.pyplot as plt
import numpy as np

ex_data = np.flip(sim.get_array(center = mp.Vector3(), size = mp.Vector3(20, 0, 20), component = mp.Ex), axis = 1)
ey_data = np.flip(sim.get_array(center = mp.Vector3(), size = mp.Vector3(20, 0, 20), component = mp.Ey), axis = 1)
ez_data = np.flip(sim.get_array(center = mp.Vector3(), size = mp.Vector3(20, 0, 20), component = mp.Ez), axis = 1)

fig, axes = plt.subplots(nrows = 1, ncols = 3, figsize = (15, 3))
axes[0].imshow(ex_data.transpose(), interpolation="spline36", cmap="RdBu", alpha=0.9)
axes[1].imshow(ey_data.transpose(), interpolation="spline36", cmap="RdBu", alpha=0.9)
axes[2].imshow(ez_data.transpose(), interpolation="spline36", cmap="RdBu", alpha=0.9)

fig.savefig("Ex_Ey_Ez.pdf")

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (3, 3))
ax.imshow(ez_data.transpose(), interpolation="spline36", cmap="RdBu", alpha=0.9)

fig.savefig("Ez.pdf")
