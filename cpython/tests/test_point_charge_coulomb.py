from eisvogel import CoordVector, SignalCalculator, Current0D, CreateElectricDipoleWeightingField
import unittest

class TestPointChargeCoulomb(unittest.TestCase):

    def test_field(self):
        
        start_coords = CoordVector.FromTRZ(-300.0, -10.0, -20.0)
        end_coords = CoordVector.FromTRZ(300.0, 300.0, 20.0)
        tp = 5.0
        N = 4
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
        
        for cur_t in range(-10, 20):
            cur_sig = calc.ComputeSignal(track, cur_t)
            print(f"{cur_t}: {cur_sig}")

        self.assertTrue(False)

if __name__ == "__main__":
    unittest.main()
