#!/bin/bash

#----------------------------------------------------

#SBATCH -J m12i_find_clusters_b4n5     # Job name Give a name for the job here good convention would include the simulation name
#SBATCH -o find_clusters_b4n5_m12i.o%j # Name of stdout output file
#SBATCH -e find_clusters_b4n5_m12i.e%j # Name of stderr error file
#SBATCH -p normal           # Queue (partition) name skx-normal more memory if you use normal only it has less memory
#SBATCH -N 1                # Total # of nodes (must be 1 for serial) Not required for now
#SBATCH -n 1                # Total # of mpi tasks (should be 1 for serial) Not required for now
#SBATCH -t 06:59:00         # Run time (hh:mm:ss)
#SBATCH --mail-user=binod.bhattarai@sxc.edu.np
#SBATCH --mail-type=all    # Send email at begin and end of job

# Other commands must follow all #SBATCH directives...

# Launch serial code...

/home1/07428/binod/anaconda3/bin/python identify_clusters_585_to_693.py
# ---------------------------------------------------