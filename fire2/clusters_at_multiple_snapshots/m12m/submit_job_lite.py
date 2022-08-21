#!/bin/bash

#----------------------------------------------------

#SBATCH -J m12m_identify_clusters_all_at_once     # Job name Give a name for the job here good convention would include the simulation name
#SBATCH -o m12m_identify_clusters_all_at_once.o%j # Name of stdout output file
#SBATCH -e m12m_identify_clusters_all_at_once.e%j # Name of stderr error file
#SBATCH -p normal           # Queue (partition) name skx-normal more memory if you use normal only it has less memory
#SBATCH -N 1                # Total # of nodes (must be 1 for serial) Not required for now
#SBATCH -n 1                # Total # of mpi tasks (should be 1 for serial) Not required for now
#SBATCH -t 07:00:00         # Run time (hh:mm:ss)
#SBATCH --mail-user=binodbhattarai.info.np@gmail.com
#SBATCH --mail-type=all    # Send email at begin and end of job

# Other commands must follow all #SBATCH directives...

# Launch serial code...

/home1/07428/binod/anaconda3/bin/python m12m_identify_clusters_all_at_once.py


# ---------------------------------------------------

