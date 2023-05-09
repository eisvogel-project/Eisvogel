import argparse
from eisvogel import CoordVector, SignalCalculator, Current0D

def compute_signal(wf_path):
    
    b = 10
    tstart = -250.0
    tend = 250.0
    charge = 1
    beta = 0.9

    points = [CoordVector.FromTXYZ(tstart, beta * tstart, 0, b),
              CoordVector.FromTXYZ(tend, beta * tend, 0, b)]
    charges = [charge]
    track = Current0D.FromSegments(points, charges)

    calc = SignalCalculator(wf_path)

    for cur_t in range(-10, 20):
        cur_sig = calc.ComputeSignal(track, cur_t)
        print(f"{cur_t}: {cur_sig}")

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--wf_path", action = "store", dest = "wf_path")
    args = vars(parser.parse_args())
    
    compute_signal(**args)
