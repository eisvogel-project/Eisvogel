#!/bin/bash

export EISVOGELDIR=/home/windischhofer/Eisvogel/
source ${EISVOGELDIR}/setup_mdwy.sh

cd ${EISVOGELDIR}/build/
mpirun --mca orte_base_help_aggregate 0 -np 128 applications/dipole/dipole
