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

        vvec = np.array([0.9, 0.0, 0.0])
        bvec = np.array([0.0, 0.0, 10.0])
        q = 1

        delta_t = 0.01
        t_min = -20
        t_max = 20
        tp = 5
        N = 6
        
        tvals, sigvals = compute_filtered_signal_analytic(vvec, bvec, q, tp, N, t_min, t_max, delta_t)

        for cur_t, cur_sig in zip(tvals, sigvals):
            print(f"{cur_t}: {cur_sig}")
            
        start_coords = CoordVector.FromTRZ(-300.0, -10.0, -20.0)
        end_coords = CoordVector.FromTRZ(300.0, 300.0, 20.0)
        os_factor = 10
        r_min = 0.1
        wf_path = "./dipole.bin"
        
        CreateElectricDipoleWeightingField(wf_path, start_coords, end_coords, tp, N, r_min, os_factor)

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
        
        for cur_t, cur_sig in zip(tvals[::100], sigvals[::100]):
            cur_sig_ev = calc.ComputeSignal(track, cur_t)
            print(f"{cur_t}: {cur_sig_ev} <-> {cur_sig}")

        self.assertTrue(True)

if __name__ == "__main__":
    unittest.main()
