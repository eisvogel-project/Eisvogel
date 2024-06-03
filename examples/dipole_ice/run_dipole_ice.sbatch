#!/bin/bash
#SBATCH --job-name=dipole_ice
#SBATCH --output=/home/windischhofer/data/windischhofer/eisvogel/dipole_ice.out
#SBATCH --error=/home/windischhofer/data/windischhofer/eisvogel/dipole_ice.err
#SBATCH --account=pi-avieregg
#SBATCH --time=24:00:00
#SBATCH --partition=avieregg
#SBATCH --nodes=4
#SBATCH --exclusive
#SBATCH --ntasks-per-node=32
#SBATCH --mem=257660
#SBATCH --mail-type=NONE
#SBATCH --mail-user=None

sh /home/windischhofer/Eisvogel/examples/dipole_ice/run_dipole_ice.sh
