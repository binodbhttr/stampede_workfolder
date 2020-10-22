#import sys
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
#import fof_v
import pdb
import gizmo_analysis as gizmo
import utilities as ut
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import pickle
import seaborn as sns
import pandas as pd
from astropy.io import ascii
from astropy.table import Table
#import alpha_functions

params = {"font.family":"serif","mathtext.fontset":"stix"}
matplotlib.rcParams.update(params)
matplotlib.rcParams['pdf.fonttype'] = 42

#These are constansts. Do not change !!!!
MsunToGm = 1.99e33
KpcToCm = 3.086e21
mp = 1.67e-24
#bin_edge = 10. #changing the bin_edge
bin_edge = 30.

snap = 596 #snapshot to begin with
simname = ['/scratch/projects/xsede/GalaxiesOnFIRE/mhdcv/m12i_res7100_mhdcv/1Myr/fire2']
#catname = ['cloud_props_m12m_600.txt','cloud_props_m12i_rl_600.txt','cloud_props_m12f_600.txt']
#picklefile = ['m12m_599_clusters_indices.pkl','m12i_clusters_599.pkl','m12f_clusters_599.pkl']

simdir = "/scratch/projects/xsede/GalaxiesOnFIRE/mhdcv/m12i_res7100_mhdcv/1Myr/fire2/"
snapnumber = 596
part = gizmo.io.Read.read_snapshots(["all"],"snapshot_index", snapnumber, simulation_name=simname, simulation_directory=simdir, assign_hosts=True,assign_hosts_rotation=True)
#part = gizmo.io.Read.read_snapshots('all', 'snapshot_index', snap,simname[s],assign_host_principal_axes=True)