import numpy as np
import matplotlib.pyplot as plt
from NuRadioReco.utilities import units, fft
import NuRadioMC.SignalGen.askaryan
import argparse
import scipy.signal

parser = argparse.ArgumentParser()
parser.add_argument('--eisvogel_output', type=str, default='shower_radio_emission.csv', help='name of the file containing the Eisvogel result')
parser.add_argument('--shower_energy', type=float, default=1.e18, help='energy of the shower, in eV')
parser.add_argument('--shower_type', type=str, default='EM', help='Options are "EM" to simulate an electromagnetic and "HAD" to simulate a hadronic shower.')
parser.add_argument('--i_shower', type=int, default=9, help='Index of the shower profile in the library. Can be between 0 and 9 for the NuRadioMC library.')
args = parser.parse_args()

shower_position = [-112., -162.]
index_of_refraction = 1.78

cherenkov_angle = np.arccos(1. / index_of_refraction)
viewing_angle = np.arctan(shower_position[1] / shower_position[0])

def filter_func(frequencies, tp=1., N=6):
    return 1. / (1 + 2.j * np.pi  * frequencies * tp / N) ** (N + 1)

data_eisvogel = np.genfromtxt(
    args.eisvogel_output,
    delimiter=','
)
times = data_eisvogel[:, 0]
sampling_rate = 1. / (times[1] - times[0])
arz_sampling_rate = 20.
arz_trace = NuRadioMC.SignalGen.askaryan.get_time_trace(
    args.shower_energy,
    viewing_angle,
    data_eisvogel.shape[0] * (arz_sampling_rate / sampling_rate),
    1. / arz_sampling_rate,
    args.shower_type,
    index_of_refraction,
    np.sqrt(shower_position[0]**2 + shower_position[1]**2),
    "ARZ2020",
    iN=args.i_shower
)
arz_trace = scipy.signal.resample(arz_trace, times.shape[0])
freqs = np.fft.rfftfreq(times.shape[0], 1. / sampling_rate)
arz_spec = fft.time2freq(arz_trace, sampling_rate) * filter_func(freqs)
arz_trace = fft.freq2time(arz_spec, sampling_rate)
spec_eisvogel = fft.time2freq(data_eisvogel[:, 1], 1)

fig1 = plt.figure(figsize=(8, 6))
ax1_1 = fig1.add_subplot(211)
ax1_2 = fig1.add_subplot(212)
ax1_1.plot(
    times,
    data_eisvogel[:, 1] / np.max(data_eisvogel[:, 1]),
    label='Eisvogel'
)
ax1_1.plot(
    times[:arz_trace.shape[0]],
    arz_trace / np.max(arz_trace),
    label='ARZ',
    alpha=.7
)
ax1_2.plot(
    freqs,
    np.abs(spec_eisvogel) / np.max(np.abs(spec_eisvogel)),
    label='Eisvogel'
)
ax1_2.plot(
    freqs,
    np.abs(arz_spec) / np.max(np.abs(arz_spec)),
    label='ARZ'
)
ax1_1.set_xlim([1150, 1300])
ax1_1.set_xlabel('t [ns]')
ax1_1.set_ylabel('U [a.u.]')
ax1_1.grid()
ax1_1.legend()
ax1_2.grid()
ax1_2.legend()
ax1_2.set_xlim([0, 1])
ax1_2.set_xlabel('f [GHz]')
ax1_2.set_ylabel('U [a.u.]')
plt.show()