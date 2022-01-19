import matplotlib as plt2 #need this for patches for shaded circles
import numpy as np
import matplotlib.pyplot as plt
import gizmo_analysis as gizmo 
import utilities as ut
import matplotlib.colors as colors
from matplotlib import rc #to use Latex math symbols like 'phi'
import astropy
from astropy.io import ascii
import matplotlib
import pdb
from importlib import reload
from sl_utilities import distance_functions
#pdb.set_trace()  #<--in case need to troubleshoot
import pickle
import os



simname = 'm12i_res7100_mhdcv'
simdir = '/scratch/projects/xsede/GalaxiesOnFIRE/cr_suite/m12i_res7100/mhdcv/1Myr/fire2/'
start= 596
datapath="/home1/07428/binod/stampede_workfolder/fire2/young_stars_data/"
for snapnumber in range(start,630,1):
    export_stars={}
    part = gizmo.io.Read.read_snapshots(['all'],'snapshot_index', snapnumber, simulation_name=simname, simulation_directory=simdir, assign_hosts_rotation=True, assign_hosts=True)  
    x        = part['star'].prop('host.distance.principal')[:,0]
    y        = part['star'].prop('host.distance.principal')[:,1]
    age      = part['star'].prop('age')
    keep = np.where(age <= .001)
    x_young=x[keep]
    y_young=y[keep]
    export_stars={"x_young":x_young,"y_young":y_young}
    file_name="fire2_young_clusters__snapshot_"+str(snapnumber)+".pkl"

    with open(datapath+file_name, 'wb') as output:
      # Pickle dictionary using protocol 0.
      pickle.dump(export_stars, output)


