from eisvogel import CoordVector, SignalCalculator, Current0D, CreateElectricDipoleWeightingField
import unittest
import numpy as np

def impulse_response(t, tp, N):
    return 1.0 / (tp * np.math.factorial(N - 1)) * np.exp(-N * t / tp) * (N * t / tp) ** N

def filter(sig_raw, delta_t, tp, N, tmax = None):
    if tmax is None:
        tmax = 5 * tp

    tkerns = np.linspace(0, tmax, int(tmax / delta_t))
    kernel = impulse_response(tkerns, tp, N)
    sig_filt = np.convolve(sig_raw, kernel) * delta_t

    return sig_filt[:len(sig_raw)]

def compute_filtered_signal_analytic(vvec, bvec, q, tp, N, t_min, t_max, delta_t):
    
    # velocity and impact parameter must be orthogonal
    assert np.inner(vvec, bvec) == 0

    ts = np.linspace(t_min, t_max, 1 + int((t_max - t_min) / delta_t))
    
    # compute raw signal
    v = np.linalg.norm(vvec)
    b = np.linalg.norm(bvec)
    gamma = 1.0 / np.sqrt(1 - v ** 2)
    
    Eperp = 1.0 / (4 * np.pi) * gamma * q * b / np.power(b ** 2 + (gamma * v * ts) ** 2, 1.5)
    Epara = -1.0 / (4 * np.pi) * gamma * q * v * ts / np.power(b ** 2 + (gamma * v * ts) ** 2, 1.5)

    vunitvec = vvec / v
    bunitvec = bvec / b
    antennaunitvec = np.array([0.0, 0.0, 1.0])

    sig_raw = Epara * np.inner(vunitvec, antennaunitvec) + Eperp * np.inner(-bunitvec, antennaunitvec)
    sig_filt = filter(sig_raw, delta_t, tp, N)

    return ts, sig_filt

class TestPointChargeCoulomb(unittest.TestCase):
    
    def test_field(self):

        threshold = 4e-2
        
        test_interval = [-2, 20]
        beta = 0.9
        b = 10
        q = 1
        
        tp = 5
        N = 6
        
        vvec = np.array([beta, 0.0, 0.0])
        bvec = np.array([0.0, 0.0, b])

        delta_t = tp / 100
        t_min_analytic = test_interval[0] - 5 * tp
        t_max_analytic = test_interval[1]
        
        tvals, sigvals = compute_filtered_signal_analytic(vvec, bvec, q, tp, N, t_min_analytic, t_max_analytic, delta_t)
        for_test = np.logical_and(tvals > test_interval[0], tvals < test_interval[1])

        tvals_test = tvals[for_test]
        sigvals_test = sigvals[for_test]

        v = np.linalg.norm(vvec)
        b = np.linalg.norm(bvec)
        t_valid_sig = test_interval[0] - 5 * tp
        t_min_ev = 1.0 / (1 - v ** 2) * (t_valid_sig - np.sqrt(b**2 + (t_valid_sig ** 2 - b ** 2) * v ** 2))
        t_max_ev = test_interval[1]
            
        start_coords = CoordVector.FromTRZ(-300.0, -10.0, -20.0)
        end_coords = CoordVector.FromTRZ(300.0, 300.0, 20.0)
        os_factor = 10
        r_min = 0.1
        wf_path = "./dipole.bin"
        
        CreateElectricDipoleWeightingField(wf_path, start_coords, end_coords, tp, N, r_min, os_factor)
        
        points = [CoordVector.FromTXYZ(t_min_ev, beta * t_min_ev, 0, b),
                  CoordVector.FromTXYZ(t_max_ev, beta * t_max_ev, 0, b)]
        charges = [q]
        track = Current0D.FromSegments(points, charges)
        
        calc = SignalCalculator(wf_path)
        
        for cur_t, cur_sig in zip(tvals_test, sigvals_test):
            cur_sig_ev = calc.ComputeSignal(track, cur_t)
            assert abs(cur_sig_ev / cur_sig - 1) < threshold

if __name__ == "__main__":
    unittest.main()
