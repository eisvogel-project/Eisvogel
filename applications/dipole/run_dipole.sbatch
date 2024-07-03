#!/bin/bash
#SBATCH --job-name=dipole
#SBATCH --output=dipole.out
#SBATCH --error=dipole.err
#SBATCH --account=pi-avieregg
#SBATCH --time=24:00:00
#SBATCH --partition=avieregg
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks-per-node=64
#SBATCH --mem=257660
#SBATCH --mail-type=NONE
#SBATCH --mail-user=None

sh /home/windischhofer/Eisvogel/examples/dipole_ice/run_dipole.sh
