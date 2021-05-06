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
import pickle



simname = 'm12i_res7100_mhdcv'
simdir = '/scratch/projects/xsede/GalaxiesOnFIRE/cr_suite/m12i_res7100/mhdcv/1Myr/fire2/'

snapshot_start=596
snapshot_end=696 #ran out of memory after 646


#Loading the sample cluster to be tracked and sorting its id and id_child
cluster_group="snapshot596" #Remember to change it if you  are changing the star cluster you are tracking in given snapshot

data_path="./fire2_data_pkl/" #creating a path to store the data only if it does not exist
gas_datapath="./fire2_gas_data_pkl/"

###############################################
###############################################
#loading data of all clusters (id, id_children and all)

cluster_file_name="clusters_"+simname+"_snapshot_"+str(snapshot_start) 
with open(data_path+cluster_file_name, "rb") as fp:
    import_cluster = pickle.load(fp)




total_clusters=len(import_cluster)         #count the total no. of clusters we have in clusters file
total_snaps=snapshot_end-snapshot_start+1  #count the total no. of snapshots we are going to track
snap=snapshot_start                        #Mark the beginning of the snapshot
test_cluster=1                             #Mark the beginning of the cluster that we would be tracking
ind_tracked={}                             #finding the indices of the tracked stars in each snapshot for each cluster.

tracked_stars_each_snap={}     #dictionary for storing data from all clusters and in each snapshot

colors=['cyan','blue','green','magenta','yellow','teal','brown','darkslategray','lime','red','orange','purple','rosybrown','pink','navy','olive','cornflowerblue']
##################################

for i in range(total_snaps):               #we run a for loop until the end of all snapshots
  
  cluster_data_name="all_clusters_at_snapshot_"+str(snapnumber)+".pkl" 
  with open(datapath+cluster_data_name, "rb") as input:
    tracked_data= pickle.load(input)
    
  part=gizmo.io.Read.read_snapshots(['star'],'snapshot_index', snap, simulation_name=simname, simulation_directory=simdir, assign_hosts=True, assign_hosts_rotation=True)               #snap is the snapshot number here that changes everytime the loop iterates. It starts with sanpshot_start

  age=part['star'].prop('age')
  x=part['star'].prop('host.distance.principal')[:,0] #x component of the position of all stars 
  y=part['star'].prop('host.distance.principal')[:,1] #y component of the position of all stars
  z=part['star'].prop('host.distance.principal')[:,2] #z component of the position of all stars
  Rxy = part['star'].prop('host.distance.principal.cylindrical')[:,0]
  
  cluster_count=0                       #This is the test cluster to begin with. It resets to 1 at each snapshot 
  tracked_data_all_clusters={}         #This dictionary keeps the tracked information for all star clusters in given snapshot. It also resets with snapshots.
  for j in range(total_clusters):
    new_stars={}                    #This dictionary holds the tracked information for a cluster temporarily until it is pushed to tracked_data_all_clusters
    xcm=tracked_data[cluster_count+1]["xcm"]
    ycm=tracked_data[cluster_count+1]["ycm"]
    #circle_radius=((x-xcm)**2+(y-ycm)**2)**(1/2)
    #region=5
    #keep_young = np.where((age <= .003) & (circle_radius<=region) & (abs(z) < 1.5))
    x_cluster=tracked_data[cluster_count+1]["x_tracked"]
    y_cluster=tracked_data[cluster_count+1]["y_tracked"]
    ax.scatter(x_cluster,y_cluster,c=colors[j])
    
    test_cluster+=1
    
  gas_file_name=simtype+"_gas_data"+str(snapnumber)+".pkl"
  with open(gas_datapath+gas_file_name, "rb") as input:
      import_gasdata = pickle.load(input)
  
  
  
  
  
  keep_young = np.where((age <= .003) & (abs(z) < 1.5))
  x0=x[keep_young]
  y0=y[keep_young]
    
    
  ################################################################
  #Now lets write our tracked data from each snapshots to a file
  new_stars={"ind_tracked":ind_tracked,"age_tracked":age_tracked,"x_tracked":x_tracked,"y_tracked":y_tracked,"z_tracked":z_tracked,
    "vx_tracked":vx_tracked,"vy_tracked":vy_tracked,"vz_tracked":vz_tracked,
    "mass_tracked":mass_tracked,"xcm":xcm,"ycm":ycm,"zcm":zcm,"delta_rxyz":delta_rxyz,"rmax":rmax,
    "ymax":ymax,"ymin":ymin,"xmax":xmax,"xmin":xmin,"avg_delta_rxyz":avg_delta_rxyz,
    "vR_cyl_tracked":vR_cyl_tracked,"vphi_cyl_tracked":vphi_cyl_tracked,"vz_cyl_tracked":vz_cyl_tracked}
    
    
  tracked_data_all_clusters.update({test_cluster:tracked_data}) #access it using tracked_data_all_clusters[clusterid]["key"]
  
  
  #We have collected all tracked information for our test clusters for one snapshot now which is in dictionary tracked_data_all_clusters.
  #We would store that dictionary in a pickle file as a backup data from each snapshot about or all star cluster.
  ###Now moving out of the j loop. We scanned all clusters and stored the tracked data into a dictionary tracked_data_all_clusters  
  file_name="all_clusters_at_snapshot_"+str(snap)+".pkl"
  with open(path+file_name, 'wb') as output:
    pickle.dump(tracked_data_all_clusters, output)
  print("\n Stored tracked data of all stars clusters in the snapshot no.",snap,"to filename:",file_name,"\n#####\n")
  #Note we are still in i loop which is scanning each snapshot !!!!!!!!!!!!!!!!!!
  tracked_data_all_clusters_each_snap.update({snap:tracked_data_all_clusters}) 
  snap=snap+1

#Came out of i loop 
  
#Finally we scanned through all the snapshots and now we come out of the loop to store the final file that contains all the tracked information of all clusters 
#Now we store the tracked data of all clusters into a dictionary with snapshot number as the key. This would be our final dictionary of dictinaries. 
with open(path+"total_data_all_clusters_all_snapshots.pkl", 'wb') as output:
  pickle.dump(tracked_data_all_clusters_each_snap, output) #access the data using tracked_data_all_clusters_each_snap[snapshot][cluster_id]["x_tracked"])
  
#Expected Result: tracked_data_all_clusters_each_snap={596:{1:tracked_data,2:tracked_data..},597:{1:tracked_data,2:tracked_data..upto total clusters},598....     upto total snapshots}
  
  
  
  