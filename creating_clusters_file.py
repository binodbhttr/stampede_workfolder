import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors
import numpy as np
import gizmo_analysis as gizmo
import utilities as ut
from fof_analysis import fof
from matplotlib import rc #to use Latex math symbols like 'phi'

############################################################################
#read in sim files and find relevant particles
############################################################################
#STAMPEDE
simname = 'm12i_res7100_mhdcv'
simdir = '/scratch/projects/xsede/GalaxiesOnFIRE/mhdcv/m12i_res7100_mhdcv/1Myr/fire2/'

#reading in snapshot_times.txt file to get snapshot numbers and times
# columns are: snapshot scale-factor redshift time[Gyr] time_width[Myr]
snapnumber = 596 #the particular snapshot we're looking at right now

#for troubleshooting on local machine
#LOCAL
#simname='m12i.res57000'
#simdir='/Users/sloebman/Dropbox/RESEARCH/ANDREW/SIMS/m12i.res57000/'
#snapnumber = 600

part = gizmo.io.Read.read_snapshots(['all'],'snapshot_index', snapnumber, simulation_name=simname, simulation_directory=simdir, assign_hosts_rotation=True, assign_hosts=True)  

srxyz     = part['star'].prop('host.distance.total')
sRxy      = part['star'].prop('host.distance.principal.cylindrical')[:,0]
sx        = part['star'].prop('host.distance.principal')[:,0]
sy        = part['star'].prop('host.distance.principal')[:,1]
sz        = part['star'].prop('host.distance.principal')[:,2] 
smass     = part['star'].prop('mass')
sfeh      = part['star'].prop('metallicity.fe')
sids      = part['star'].prop('id')
sid_child = part['star'].prop('id.child')
sage      = part['star'].prop('age')


############################################################################
#identifying star clusters/associations
###########################################################################

#select young stars within the disk for star clusters
si = np.where((sage <= .003) & (srxyz < 20) & (abs(sz) < 1.5))

linking_length = 0.01 #10 parsec
ncut           = 10 #10 star particles

#Running fof
ind, xsp, ysp, zsp, msp, grpid, r90, r50, rmax =fof.find(sx[si],sy[si],sz[si], b=linking_length, mass=smass[si], ncut=ncut)
ngroup = len(msp)

