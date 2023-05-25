import argparse
import numpy as np
from eisvogel import DeltaVector, CoordVector, FieldVector, SignalCalculator, SparseCurrentDensity3D

def compute_signal(wf_path):

    b = 10
    tstart = -250.0
    tend = 250.0
    charge = 1
    beta = 0.9

    voxel_size = DeltaVector.FromDeltaTXYZ(0.1, 0.1, 0.1, 0.1)
    curr_dist = SparseCurrentDensity3D(voxel_size)

    for cur_t in np.linspace(tstart, tend, 5000):
        
        cur_x = beta * cur_t
        cur_y = 0
        cur_z = b

        cur_jx = charge * beta / (0.1 ** 3)
        cur_jy = 0
        cur_jz = 0

        curr_dist.addCurrentElement(CoordVector.FromTXYZ(cur_t, cur_x, cur_y, cur_z),
                                    FieldVector.FromXYZ(cur_jx, cur_jy, cur_jz))

    calc = SignalCalculator(wf_path)

    for cur_t in range(-10, 20):
        cur_sig = calc.ComputeSignal(curr_dist, cur_t)
        print(f"{cur_t}: {cur_sig}")

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--wf_path", action = "store", dest = "wf_path")
    args = vars(parser.parse_args())
    
    compute_signal(**args)
