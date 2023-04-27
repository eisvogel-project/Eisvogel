import pickle
import json
import os
import requests
import h5py
import numpy as np
from NuRadioReco.utilities import units

if not os.path.exists('./ARZ_library_v1.2.pkl'):
	URL = 'https://rnog-data.zeuthen.desy.de/shower_library/library_v1.2.pkl'
	r = requests.get(URL)
	if r.status_code != requests.codes.ok:
		print("error in download of antenna model")
		raise IOError("error in download of antenna model")
	with open('ARZ_library_v1.2.pkl', "wb") as code:
		code.write(r.content)

import array

def write_shower(f, E, type, grammage, q):
	np.array([len(grammage)], dtype=np.int32).tofile(f)
	np.array([type == "HAD"], dtype=np.int32).tofile(f)
	np.array([E], dtype='d').tofile(f)
	(np.array(grammage, dtype='d') / (units.kg / units.m**2)).tofile(f) #assuming it's double precision float, works with both array or ndarrays
	np.array(q, dtype='d').tofile(f)

shower_library = pickle.load(open('ARZ_library_v1.2.pkl', 'rb'), encoding='latin1')
output_dict = {}
hf = h5py.File('ARZ_library_v1.2.hdf5', 'w')
f = open('shower_file', 'wb')
for key_1 in shower_library.keys():
	output_dict[key_1] = {}
	group = hf.create_group(key_1)
	energies = []
	for key_2 in shower_library[key_1].keys():
		energies.append(float(key_2))
		profiles = []
		for i in range(len(shower_library[key_1][key_2]['charge_excess'])):
			profiles.append(list(shower_library[key_1][key_2]['charge_excess'][i]))
			write_shower(
				f,
				key_2,
				key_1,
				shower_library[key_1][key_2]['depth'],
				list(shower_library[key_1][key_2]['charge_excess'][i])
			)
		output_dict[key_1][key_2] = {
			'depth': list(shower_library[key_1][key_2]['depth']),
			'charge_excess': profiles
		}
		group_2 = group.create_group(str(key_2))
		group_2.create_dataset('charge_excess', data=profiles)
		group_2.create_dataset('depth', data=list(shower_library[key_1][key_2]['depth']))

	group.create_dataset('energies', data=list(np.log10(energies)))
hf.close()
json.dump(output_dict, open('ARZ_library_v1.2.json', 'w'))
