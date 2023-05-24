import meep as mp

cell = mp.Vector3(16, 8, 0)

geometry = [
    mp.Block(
        mp.Vector3(mp.inf, 1, mp.inf),
        center=mp.Vector3(),
        material=mp.Medium(epsilon=12),
    )
]

sources = [
    mp.Source(
        mp.ContinuousSource(frequency=0.15), component=mp.Ez, center=mp.Vector3(-7, 0)
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

sim.run(until=200)


import matplotlib.pyplot as plt
import numpy as np

eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize = (11, 3))
axes[0].imshow(eps_data.transpose(), interpolation="spline36", cmap="binary")

ez_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)
axes[1].imshow(eps_data.transpose(), interpolation="spline36", cmap="binary")
axes[1].imshow(ez_data.transpose(), interpolation="spline36", cmap="RdBu", alpha=0.9)

fig.savefig("test.pdf")

print("done")
