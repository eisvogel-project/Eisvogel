import argparse
import pyeisvogel

def compute_signal(wf_path):

    calc = pyeisvogel.SignalCalculator(wf_path)
    print(calc)

    vec = pyeisvogel.CoordVector.FromTXYZ(1, 2, 3, 4)
    print(vec)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--wf_path", action = "store", dest = "wf_path")
    args = vars(parser.parse_args())
    
    compute_signal(**args)
