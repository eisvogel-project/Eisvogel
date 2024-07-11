#!/bin/bash

export EISVOGELDIR=/home/weipow/Eisvogel

source /project/avieregg/software/midway3-setup.sh 
module load python
source /project/avieregg/eisvogel/setup_mdwy.sh

cd ${EISVOGELDIR}/build/
mpirun --mca orte_base_help_aggregate 0 -np 128 applications/dipole/dipole
