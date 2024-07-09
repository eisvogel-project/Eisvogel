#!/bin/bash

source /project/avieregg/software/midway3-setup.sh 
module load python
source /project/avieregg/eisvogel/setup_mdwy.sh

mpirun --mca orte_base_help_aggregate 0 -np 128 /home/weipow/Eisvogel/build/applications/dipole/dipole
