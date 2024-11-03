#!/bin/bash
#SBATCH --job-name=dipole
#SBATCH --output=dipole.out
#SBATCH --error=dipole.err
#SBATCH --account=pi-avieregg
#SBATCH --time=24:00:00
#SBATCH --partition=avieregg
#SBATCH --nodes=5
#SBATCH --exclusive
#SBATCH --ntasks-per-node=32
#SBATCH --mem=257660
#SBATCH --mail-type=NONE
#SBATCH --mail-user=None

export EISVOGELDIR=/home/windischhofer/Eisvogel/
source ${EISVOGELDIR}/setup_mdwy.sh

cd ${EISVOGELDIR}/build/
mpirun --mca orte_base_help_aggregate 0 -np 160 applications/dipole/dipole
