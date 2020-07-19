#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import gizmo_analysis as gizmo
import utilities as ut
from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D



# In[ ]:


simname = 'm12i_res7100_mhdcv'
simdir = '/scratch/projects/xsede/GalaxiesOnFIRE/mhdcv/m12i_res7100_mhdcv/1Myr/'
snapnumber = 691


# In[ ]:
snapshot_start=671
snapshot_end=672
part=[] #part is a list of arrarys
id=[]
id_child=[]
age=[]
x=[]
y=[]
z=[]
n=[]
mass=[]
for i in range(snapshot_end+1):
  if i<snapshot_start:
    part.append(0) 
    id.append(0)
    id_child.append(0) 
    age.append(0)
    x.append(0)
    y.append(0)
    z.append(0)
    n.append(0)
    mass.append(0)
  else:  
    part.append(gizmo.io.Read.read_snapshots(['star'],'snapshot_index', i, simulation_name=simname, simulation_directory=simdir, assign_hosts=True, assign_hosts_rotation=True)) #691 is the snapnumber here
    id.append(part[i]['star'].prop('id'))
    id_child.append(part[i]['star'].prop('id.child'))
    age.append(part[i]['star'].prop('age'))
    x.append(part[i]['star'].prop('host.distance.principal')[:,0]) #x component of the position of all stars 
    y.append(part[i]['star'].prop('host.distance.principal')[:,1]) #y component of the position of all stars
    z.append(part[i]['star'].prop('host.distance.principal')[:,2]) #z component of the position of all stars
    n.append(len(x[i])) # counting the total no. of star particles
    mass.append(part[i]['star']['mass']) #mass of all stars in snapshot 691
    print("\n#######################\n#######################\nLoaded id,id_child,age,x,y,z, mass and number of particles for snapshot no.",i)
    print("Total no of particles in this shapshot no.",i,"is",n[i])



print("All the ages: of the last snapshot is",age[snapshot_end])