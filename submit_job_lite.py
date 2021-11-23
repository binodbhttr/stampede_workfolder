#!/bin/bash

#----------------------------------------------------

#SBATCH -J fire2_b4n4_clusters     # Job name Give a name for the job here good convention would include the simulation name
#SBATCH -o fire2_b4n4_clusters.o%j # Name of stdout output file
#SBATCH -e fire2_b4n4_clusters.e%j # Name of stderr error file
#SBATCH -p normal           # Queue (partition) name skx-normal more memory if you use normal only it has less memory
#SBATCH -N 1                # Total # of nodes (must be 1 for serial) Not required for now
#SBATCH -n 1                # Total # of mpi tasks (should be 1 for serial) Not required for now
#SBATCH -t 6:00:00         # Run time (hh:mm:ss)
#SBATCH --mail-user=binod.bhattarai@sxc.edu.np
#SBATCH --mail-type=all    # Send email at begin and end of job

# Other commands must follow all #SBATCH directives...

# Launch serial code...

ipython /home1/07428/binod/stampede_workfolder/fire2/clusters_at_multiple_snapshots/fire2_chemistry_analysis_multiple_snapshots/identify_clusters_at_multiple_snapshots.py

# ---------------------------------------------------

