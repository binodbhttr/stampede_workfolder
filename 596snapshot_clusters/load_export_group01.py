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
import os




simname = 'm12i_res7100_mhdcv'
simdir = '/scratch/projects/xsede/GalaxiesOnFIRE/mhdcv/m12i_res7100_mhdcv/1Myr/1Myr_fire2'

# In[ ]:
snapshot_start=596
snapshot_end=696 #ran out of memory after 646


#Loading the sample cluster to be tracked and sorting its id and id_child
cluster_groupid="snapshot596_cluster_group01" #Remember to change it if you  are changing the star culuster you are tracking in given snapshot

path="./data/"+cluster_groupid+"/" #creating a path to store the data only if it does not exist
if not os.path.exists(path):
  os.makedirs(path)

id_test_cluster=np.array([21964074, 13285005, 30978536, 68630449, 45857179, 47786639, 57702149, 19618367, 15329324, 34194763, 47378181, 13826769, 59570031, 40986229, 8908953, 19549744, 47837702, 22810807, 29187485, 15137299, 37115631, 61430465, 16546558, 36814748, 34745906, 14417153, 55390648, 16401221, 26360923, 7530384, 19930456, 18125528, 24685347, 27644515, 24071454, 22664528, 35943583, 26608080, 20640400, 13409766, 13890784, 40596381, 60695413]) #ids of the cluster to begin with 
id_child_test_cluster=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
sortind=np.argsort(id_test_cluster)
id_test_cluster_sorted=id_test_cluster[sortind]
id_child_test_cluster_sorted=id_child_test_cluster[sortind]
print("The total no. of stars in this cluster is",len(id_test_cluster_sorted))
print("Sorted ids of this cluster is",id_test_cluster_sorted)
################################################################
################################################################
#colors=['cyan','blue','green','magenta','yellow','orange','purple','tan','lime','brown','grey','pink','navy','teal'] #as the ids are sorted and the matching would also give sorted ids, the arrays created from tracking would also be sorted in term 
#colors=['cyan','blue','green','magenta','yellow','orange','purple','tan','lime','brown','grey']

################################################################
################################################################
###Now let's create a function for find the matching ids in the next snapshot using a function
def matchids(id_current,id_child_current,id_next,id_child_next,id_generation_next): #this function returns the index of the ids in the next snapshot that match with the current
  ind=np.array(0)
  for i in range(len(id_current)):
    match=np.where((id_next==id_current[i])&(id_child_next==id_child_current[i])&(id_generation_next<30)) #Also add id_generations<30
    print("\nFound the id",id_next[match],"at the index",match[0],"in this snapshot")
    ind=np.append(ind,match)
  ind_tracked_id_next=ind[1:len(ind)]  #The extra element in the beginning is removed by this process
  print("\nThese are the indices of the ids that matched in current snapshot\n",ind_tracked_id_next)
  return ind_tracked_id_next  
################################################################
################################################################



total_snaps=snapshot_end-snapshot_start+1
snap=snapshot_start
ind_tracked={} #finding the indices of the tracked stars in each snapshot.
for i in range(total_snaps): 
  part=gizmo.io.Read.read_snapshots(['star'],'snapshot_index', snap, simulation_name=simname, simulation_directory=simdir, assign_hosts=True, assign_hosts_rotation=True) #snap is the snapshot number here that changes everytime the loop iterates
  
  id=part['star'].prop('id')
  id_child=part['star'].prop('id.child')
  id_generation=part['star'].prop('id.generation')
  age=part['star'].prop('age')
  x=part['star'].prop('host.distance.principal')[:,0] #x component of the position of all stars 
  y=part['star'].prop('host.distance.principal')[:,1] #y component of the position of all stars
  z=part['star'].prop('host.distance.principal')[:,2] #z component of the position of all stars
  n=len(x) # counting the total no. of star particles
  vx=part['star'].prop('host.velocity.principal')[:,0]
  vy=part['star'].prop('host.velocity.principal')[:,1]
  vz=part['star'].prop('host.velocity.principal')[:,2]
  mass=part['star']['mass'] #mass of all stars in snapshot 691
  print("\n#######################\n#######################\nLoaded id,id_child,age,x,y,z,vx,vy,vz,mass and number of particles for snapshot no.",snap)
  print("Total no of particles in this shapshot no.",snap,"is",n)
  
  ################################################################
  ################################################################
  ### Now matching the IDs for all our snapshots that are loaded within the same loop
  print("\n\nNow Matching the ids for the snapshot",snap,"########\n")
  ind_tracked=matchids(id_test_cluster_sorted,id_child_test_cluster_sorted,id,id_child,id_generation)
  
  ################################################################
  ################################################################
  #Now finding the postion, velocities, age and mass of the tracked stars in each snapshots
  age_tracked=age[ind_tracked]
  x_tracked=x[ind_tracked]
  y_tracked=y[ind_tracked]
  z_tracked=z[ind_tracked]
  vx_tracked=vx[ind_tracked]
  vy_tracked=vy[ind_tracked]
  vz_tracked=vz[ind_tracked]
  mass_tracked=mass[ind_tracked]
  
  ################################################################
  ################################################################
  #Now calculating center of mass and realted properties  
  xcm=distance_functions.cm(x_tracked,mass_tracked)
  ycm=distance_functions.cm(y_tracked,mass_tracked)
  zcm=distance_functions.cm(z_tracked,mass_tracked)
  delta_rxyz=distance_functions.dr(x_tracked,y_tracked,z_tracked,mass_tracked)
  rmax=distance_functions.drmax(x_tracked,y_tracked,z_tracked,mass_tracked)
  ymax=(ycm+1.1*rmax)
  ymin=(ycm-1.1*rmax)
  xmax=(xcm+1.1*rmax)
  xmin=(xcm-1.1*rmax)
  avg_delta_rxyz=np.mean(np.absolute(delta_rxyz))
  
  
  ################################################################
  ################################################################
  #Now lets write our tracked data from each snapshots to a file
  
  
    
  dict_exportdata={"ind_tracked":ind_tracked,
"age_tracked":age_tracked,"x_tracked":x_tracked,"y_tracked":y_tracked,"z_tracked":z_tracked,"vx_tracked":vx_tracked,"vy_tracked":vy_tracked,"vz_tracked":vz_tracked,"mass_tracked":mass_tracked,"xcm":xcm,"ycm":ycm,"zcm":zcm,
"delta_rxyz":delta_rxyz,"rmax":rmax,"ymax":ymax,"ymin":ymin,"xmax":xmax,"xmin":xmin,"avg_delta_rxyz":avg_delta_rxyz}
  file_name=cluster_groupid+"_tracked_data_snapshot"+str(snap)
  ut.io.file_hdf5(path+file_name, dict_exportdata)
  print("\n Stored data from the snapshot no.",snap,"to filename:",file_name,".hdf5\n#####\n")
  snap=snap+1
  




