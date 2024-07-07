#!/bin/bash

export EISVOGELDIR=/home/windischhofer/Eisvogel/
source ${EISVOGELDIR}/setup_mdwy.sh

cd ${EISVOGELDIR}/build/
source ./setup.sh
mpirun --mca orte_base_help_aggregate 0 -np 64 applications/dipole/dipole
