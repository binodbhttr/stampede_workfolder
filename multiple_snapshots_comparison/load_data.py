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


part_691 = gizmo.io.Read.read_snapshots(['star'],'snapshot_index', snapnumber, simulation_name=simname, simulation_directory=simdir, assign_hosts=True, assign_hosts_rotation=True) #691 is the snapnumber here

part_692 = gizmo.io.Read.read_snapshots(['star'],'snapshot_index',692 , simulation_name=simname, simulation_directory=simdir, assign_hosts=True, assign_hosts_rotation=True)

part_693 = gizmo.io.Read.read_snapshots(['star'],'snapshot_index',693 , simulation_name=simname, simulation_directory=simdir, assign_hosts=True, assign_hosts_rotation=True)

part_694 = gizmo.io.Read.read_snapshots(['star'],'snapshot_index',694 , simulation_name=simname, simulation_directory=simdir, assign_hosts=True, assign_hosts_rotation=True)

part_695 = gizmo.io.Read.read_snapshots(['star'],'snapshot_index',695 , simulation_name=simname, simulation_directory=simdir, assign_hosts=True, assign_hosts_rotation=True)

part_696 = gizmo.io.Read.read_snapshots(['star'],'snapshot_index',696 , simulation_name=simname, simulation_directory=simdir, assign_hosts=True, assign_hosts_rotation=True)

# Loading age and id data from snapshot 691
part_691['star'].keys()
id_691=part_691['star'].prop('id')
id_child_691=part_691['star'].prop('id.child')
print(" \n Loaded the ids and id.child of the early snapshot (691)",id_691,id_child_691)
ageall_691=part_691['star'].prop('age')






# Loading age and id data from snapshot 692
part_692['star'].keys()
id_692=part_692['star'].prop('id')
id_child_692=part_692['star'].prop('id.child')
print(" \n Loaded the ids and id.child of the latest snapshot (692)",id_692,id_child_692)
ageall_692=part_692['star'].prop('age')

# Loading age and id data from snapshot 693
part_693['star'].keys()
id_693=part_693['star'].prop('id')
id_child_693=part_693['star'].prop('id.child')
print(" \n Loaded the ids and id.child of the latest snapshot (693)",id_693,id_child_693)
ageall_693=part_693['star'].prop('age')


# Loading age and id data from snapshot 694
part_694['star'].keys()
id_694=part_694['star'].prop('id')
id_child_694=part_694['star'].prop('id.child')
print(" \n Loaded the ids and id.child of the latest snapshot (694)",id_694,id_child_694)
ageall_694=part_694['star'].prop('age')


# Loading data from snapshot 695
part_695['star'].keys()
id_695=part_695['star'].prop('id')
id_child_695=part_695['star'].prop('id.child')
print(" \n Loaded the ids and id.child of the latest snapshot (695)",id_695,id_child_695)
ageall_695=part_695['star'].prop('age')


# Loading age and id data from snapshot 696
part_696['star'].keys()
id_696=part_696['star'].prop('id')
id_child_696=part_696['star'].prop('id.child')
print(" \n Loaded the ids and id.child of the latest snapshot (696)",id_696,id_child_696)
ageall_696=part_696['star'].prop('age') 



###################################################
xall_691 = part_691['star'].prop('host.distance.principal')[:,0] #x component of the position of all stars 
yall_691 = part_691['star'].prop('host.distance.principal')[:,1] #y component of the position of all stars
zall_691 = part_691['star'].prop('host.distance.principal')[:,2] #z component of the position of all stars
nall_691=len(xall_691) # counting the total no. of star particles
n_691=nall_691
print("\n Total particles in 691 snapshots is",n_691) 
x_691=xall_691[0:n_691] # select x components from 0 to n
y_691=yall_691[0:n_691] # select y components from 0 to n
z_691=zall_691[0:n_691] # select z components from 0 to n
age_691=ageall_691[0:n_691]*1e9 #age converted to years Note: NOT IN GYR ANYMORE !!!!!!!!!!!
             
############################
###################################################
xall_692 = part_692['star'].prop('host.distance.principal')[:,0] #x component of the position of all stars 
yall_692 = part_692['star'].prop('host.distance.principal')[:,1] #y component of the position of all stars
zall_692 = part_692['star'].prop('host.distance.principal')[:,2] #z component of the position of all stars
nall_692=len(xall_692) # counting the total no. of star particles
n_692=nall_692
print("\n Total particles in 692 snapshots is",n_692) 
x_692=xall_692[0:n_692] # select x components from 0 to n
y_692=yall_692[0:n_692] # select y components from 0 to n
z_692=zall_692[0:n_692] # select z components from 0 to n
age_692=ageall_692[0:n_692]*1e9 #age converted to years Note: NOT IN GYR ANYMORE !!!!!!!!!!!
             
############################
###################################################
xall_693 = part_693['star'].prop('host.distance.principal')[:,0] #x component of the position of all stars 
yall_693 = part_693['star'].prop('host.distance.principal')[:,1] #y component of the position of all stars
zall_693 = part_693['star'].prop('host.distance.principal')[:,2] #z component of the position of all stars
nall_693=len(xall_693) # counting the total no. of star particles
n_693=nall_693
print("\n Total particles in 693 snapshots is",n_693) 
x_693=xall_693[0:n_693] # select x components from 0 to n
y_693=yall_693[0:n_693] # select y components from 0 to n
z_693=zall_693[0:n_693] # select z components from 0 to n
age_693=ageall_693[0:n_693]*1e9 #age converted to years Note: NOT IN GYR ANYMORE !!!!!!!!!!!
             
############################
###################################################
xall_694 = part_694['star'].prop('host.distance.principal')[:,0] #x component of the position of all stars 
yall_694 = part_694['star'].prop('host.distance.principal')[:,1] #y component of the position of all stars
zall_694 = part_694['star'].prop('host.distance.principal')[:,2] #z component of the position of all stars
nall_694=len(xall_694) # counting the total no. of star particles
n_694=nall_694
print("\n Total particles in 694 snapshots is",n_694) 
x_694=xall_694[0:n_694] # select x components from 0 to n
y_694=yall_694[0:n_694] # select y components from 0 to n
z_694=zall_694[0:n_694] # select z components from 0 to n
age_694=ageall_694[0:n_694]*1e9 #age converted to years Note: NOT IN GYR ANYMORE !!!!!!!!!!!
             
############################
###################################################
xall_695 = part_695['star'].prop('host.distance.principal')[:,0] #x component of the position of all stars 
yall_695 = part_695['star'].prop('host.distance.principal')[:,1] #y component of the position of all stars
zall_695 = part_695['star'].prop('host.distance.principal')[:,2] #z component of the position of all stars
nall_695=len(xall_695) # counting the total no. of star particles
n_695=nall_695
print("\n Total particles in 695 snapshots is",n_695) 
x_695=xall_695[0:n_695] # select x components from 0 to n
y_695=yall_695[0:n_695] # select y components from 0 to n
z_695=zall_695[0:n_695] # select z components from 0 to n
age_695=ageall_695[0:n_695]*1e9 #age converted to years Note: NOT IN GYR ANYMORE !!!!!!!!!!!
             
############################
###################################################
xall_696 = part_696['star'].prop('host.distance.principal')[:,0] #x component of the position of all stars 
yall_696 = part_696['star'].prop('host.distance.principal')[:,1] #y component of the position of all stars
zall_696 = part_696['star'].prop('host.distance.principal')[:,2] #z component of the position of all stars
nall_696=len(xall_696) # counting the total no. of star particles
n_696=nall_696
print("\n Total particles in 696 snapshots is",n_696) 
x_696=xall_696[0:n_696] # select x components from 0 to n
y_696=yall_696[0:n_696] # select y components from 0 to n
z_696=zall_696[0:n_696] # select z components from 0 to n
age_696=ageall_696[0:n_696]*1e9 #age converted to years Note: NOT IN GYR ANYMORE !!!!!!!!!!!
             
############################



R_691=np.sqrt(np.square(x_691)+np.square(y_691)) #calculate the radius in the xy plane
print("\n ################## \n Calculated the cylindrical radius R in xy plane for the snapshot 691",R_691)




