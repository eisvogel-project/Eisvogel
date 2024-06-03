#!/bin/bash

export EISVOGELDIR=/home/windischhofer/Eisvogel/
source ${EISVOGELDIR}/setup_mdwy.sh

cd ${EISVOGELDIR}/build/
mpirun --mca orte_base_help_aggregate 0 -np 128 /home/windischhofer/Eisvogel/build/examples/dipole_ice/dipole_ice /home/windischhofer/data/windischhofer/eisvogel/gf_dipole_homogeneous_n_1.78_meep/ /home/windischhofer/Eisvogel/examples/dipole_ice/complicated_ice.csv
