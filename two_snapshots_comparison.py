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
snapnumber = 696


# In[ ]:


part = gizmo.io.Read.read_snapshots(['star'],'snapshot_index', snapnumber, simulation_name=simname, simulation_directory=simdir, assign_hosts=True, assign_hosts_rotation=True)

part_691 = gizmo.io.Read.read_snapshots(['star'],'snapshot_index',691 , simulation_name=simname, simulation_directory=simdir, assign_hosts=True, assign_hosts_rotation=True)
# In[ ]:


part['star'].keys()
id=part['star'].prop('id')
print("Loaded the ids of the latest snapshot (696)",id)

part_691['star'].keys()
id_691=part_691['star'].prop('id')
print("Loaded the ids of the early snapshot (691)",id_691)

# In[ ]:


# 3-D position of star particle (particle number x dimension number) [kpc comoving]

starposition=part['star']['position'] # starposition is the array of position of all stars
ageall=part['star'].prop('age') #stored all ages in a

starposition_691=part_691['star']['position']
ageall_691=part_691['star'].prop('age')


 
xall = part['star'].prop('host.distance.principal')[:,0] #x component of the position of all stars 
yall = part['star'].prop('host.distance.principal')[:,1] #y component of the position of all stars
zall = part['star'].prop('host.distance.principal')[:,2] #z component of the position of all stars
nall=len(xall) # counting the total no. of star particles


xall_691 = part_691['star'].prop('host.distance.principal')[:,0] #x component of the position of all stars 
yall_691 = part_691['star'].prop('host.distance.principal')[:,1] #y component of the position of all stars
zall_691 = part_691['star'].prop('host.distance.principal')[:,2] #z component of the position of all stars
nall_691=len(xall_691) # counting the total no. of star particles



n=nall # This is the no. of stars we are going to use. n=nall means choose all, n=100 means choose 100
n_691=nall_691
print("Total particles in 696 snapshots is",n)
print("Total particles in 691 snapshots is",n_691)                
smallpart=starposition[0:n] # This is selecting a small part from the collection of all star positions

x=xall[0:n] # select x components from 0 to n
y=yall[0:n] # select y components from 0 to n
z=zall[0:n] # select z components from 0 to n
age=ageall[0:n]*1e9 #age converted to years Note: NOT IN GYR ANYMORE !!!!!!!!!!!

x_691=xall_691[0:n_691] # select x components from 0 to n
y_691=yall_691[0:n_691] # select y components from 0 to n
z_691=zall_691[0:n_691] # select z components from 0 to n
age_691=ageall_691[0:n_691]*1e9 #age converted to years Note: NOT IN GYR ANYMORE !!!!!!!!!!!


R=np.sqrt(np.square(x)+np.square(y)) #calculate the radius in the xy plane
print("Calculated the cylindrical radius R in xy plane for the snapshot 696",R)
r=np.sqrt(x**2+y**2+z**2) #calculate the spherical radius in the xyz plane
print("Calculated the spherical radius r in the xyz plane")

R_691=np.sqrt(np.square(x_691)+np.square(y_691)) #calculate the radius in the xy plane
print("Calculated the cylindrical radius R in xy plane for the snapshot 691",R_691)


# In[ ]:





# In[ ]:




