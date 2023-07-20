from eisvogel import CoordVector, SignalCalculator, Current0D, CreateElectricDipoleWeightingField
import unittest, itertools
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

def compare_eisvogel_to_analytic(bvec, vvec, charge, tp, N, time_test_interval, threshold):

    # compute filtered signal analytically    
    delta_t = tp / 100
    t_min_analytic = time_test_interval[0] - 5 * tp
    t_max_analytic = time_test_interval[1]
        
    tvals, sigvals = compute_filtered_signal_analytic(vvec, bvec, charge, tp, N, t_min_analytic, t_max_analytic, delta_t)
    for_test = np.logical_and(tvals > time_test_interval[0], tvals < time_test_interval[1])
        
    tvals_test = tvals[for_test]
    sigvals_test = sigvals[for_test]

    # compute filtered signal with Eisvogel
    v = np.linalg.norm(vvec)
    b = np.linalg.norm(bvec)
    t_valid_sig = time_test_interval[0] - 5 * tp
    t_min_ev = 1.0 / (1 - v ** 2) * (t_valid_sig - np.sqrt(b**2 + (t_valid_sig ** 2 - b ** 2) * v ** 2))
    t_max_ev = time_test_interval[1]

    pos_start_vec = bvec + vvec * t_min_ev
    pos_end_vec = bvec + vvec * t_max_ev        

    def r_xy(vec):
        return np.linalg.norm(vec[0:2])

    def z(vec):
        return vec[-1]
        
    padding = 25
    start_coords = CoordVector.FromTRZ(tvals_test[0] - t_max_ev - padding,
                                       min(r_xy(pos_start_vec), r_xy(pos_end_vec)) - padding,
                                       min(z(pos_start_vec), z(pos_end_vec)) - padding)
    end_coords = CoordVector.FromTRZ(tvals_test[1] - t_min_ev + padding,
                                     max(r_xy(pos_start_vec), r_xy(pos_end_vec)) + padding,
                                     max(z(pos_start_vec), z(pos_end_vec)) + padding)
    os_factor = 10
    r_min = 0.1
    wf_path = "./dipole.bin"
    
    CreateElectricDipoleWeightingField(wf_path, start_coords, end_coords, tp, N, r_min, os_factor)        
    
    points = [CoordVector.FromTXYZ(t_min_ev, *pos_start_vec),
              CoordVector.FromTXYZ(t_max_ev, *pos_end_vec)]
    charges = [charge]
    track = Current0D.FromSegments(points, charges)
    
    calc = SignalCalculator(wf_path)
    
    for cur_t, cur_sig in zip(tvals_test, sigvals_test):
        cur_sig_ev = calc.ComputeSignal(track, cur_t)
        if abs(cur_sig_ev / cur_sig - 1) > threshold:
            return False

    return True

class TestPointChargeCoulomb(unittest.TestCase):
    
    def test_field(self):

        threshold = 4e-2
        
        test_interval = [-2, 20]
        q = 1
        
        beta_vals = [0.8, 0.9]
        b_vals = [10, 20]        
        tp_vals = [5]
        N_vals = [6]

        for ind, (beta, b, tp, N) in enumerate(itertools.product(beta_vals, b_vals, tp_vals, N_vals)):
            with self.subTest(i = ind):        
                vvec = np.array([beta, 0.0, 0.0])
                bvec = np.array([0.0, 0.0, b])
                self.assertTrue(compare_eisvogel_to_analytic(bvec, vvec, q, tp, N, test_interval, threshold))
        
if __name__ == "__main__":
    unittest.main()
