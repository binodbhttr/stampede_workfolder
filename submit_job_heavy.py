#!/bin/bash

#----------------------------------------------------

#SBATCH -J tracking_star_clusters     # Job name Give a name for the job here good convention would include the simulation name
#SBATCH -o tracking_star_clusters.o%j # Name of stdout output file
#SBATCH -e tracking_star_clusters.e%j # Name of stderr error file
#SBATCH -p skx-normal           # Queue (partition) name skx-normal more memory if you use normal only it has less memory
#SBATCH -N 1                # Total # of nodes (must be 1 for serial) Not required for now
#SBATCH -n 1                # Total # of mpi tasks (should be 1 for serial) Not required for now
#SBATCH -t 04:00:00         # Run time (hh:mm:ss)
#SBATCH --mail-user=binod.bhattarai@sxc.edu.np
#SBATCH --mail-type=all    # Send email at begin and end of job

# Other commands must follow all #SBATCH directives...

# Launch serial code...

/home1/07428/binod/anaconda3/bin/python load_export_all_snapshots_all_clusters.py


# ---------------------------------------------------

