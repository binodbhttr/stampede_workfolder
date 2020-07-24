#!/usr/bin/env python
# coding: utf-8

#this file will load snapshots and store details of the tracked star cluster thoroughout the snapshots into individual files available at ./data

import gizmo_analysis as gizmo
import utilities as ut
from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

from sl_utilities import distinct_colours as dc
from sl_utilities import distance_functions


simname = 'm12i_res7100_mhdcv'
simdir = '/scratch/projects/xsede/GalaxiesOnFIRE/mhdcv/m12i_res7100_mhdcv/1Myr/1Myr_fire2'

# In[ ]:
snapshot_start=596
snapshot_end=598
part={} #part is a dictionary 
id={}
id_child={}
id_generation={}
age={}
x={}
y={}
z={}
n={}
mass={}
vx={}
vy={}
vz={}
total_snaps=snapshot_end-snapshot_start+1
snap=snapshot_start
for i in range(total_snaps): 
  part[snap]=gizmo.io.Read.read_snapshots(['star'],'snapshot_index', snap, simulation_name=simname, simulation_directory=simdir, assign_hosts=True, assign_hosts_rotation=True) #snap is the snapshot number here that changes everytime the loop iterates
  id[snap]=part[snap]['star'].prop('id')
  id_child[snap]=part[snap]['star'].prop('id.child')
  id_generation[snap]=part[snap]['star'].prop('id.generation')
  age[snap]=part[snap]['star'].prop('age')
  x[snap]=part[snap]['star'].prop('host.distance.principal')[:,0] #x component of the position of all stars 
  y[snap]=part[snap]['star'].prop('host.distance.principal')[:,1] #y component of the position of all stars
  z[snap]=part[snap]['star'].prop('host.distance.principal')[:,2] #z component of the position of all stars
  n[snap]=len(x[snap]) # counting the total no. of star particles
  vx[snap]=part[snap]['star'].prop('host.velocity.principal')[:,0]
  vy[snap]=part[snap]['star'].prop('host.velocity.principal')[:,1]
  vz[snap]=part[snap]['star'].prop('host.velocity.principal')[:,2]
  mass[snap]=part[snap]['star']['mass'] #mass of all stars in snapshot 691
  print("\n#######################\n#######################\nLoaded id,id_child,age,x,y,z,vx,vy,vz,mass and number of particles for snapshot no.",snap)
  print("Total no of particles in this shapshot no.",snap,"is",n[snap])
  snap=snap+1
  
  
