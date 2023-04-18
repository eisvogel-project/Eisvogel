import argparse
import pyeisvogel

def compute_signal(wf_path):

    wf = pyeisvogel.WeightingField.from_path(wf_path)
    print(wf)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--wf_path", action = "store", dest = "wf_path")
    args = vars(parser.parse_args())
    
    compute_signal(**args)
