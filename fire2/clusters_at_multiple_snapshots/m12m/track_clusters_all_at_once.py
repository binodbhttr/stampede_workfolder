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


simtype="fire2"
simname = 'm12m_res7100_mhdcv'
simdir = '/scratch/projects/xsede/GalaxiesOnFIRE/cr_suite/m12m_res7100/mhdcv/1Myr/'

data_start=596
data_end=668

for s in range(data_start,data_end+3,3): #doing the tracking for clusters starting from data_start until data_end. Note we identified clusters in 3Myrs gap

  snapshot_start=s
  snapshot_end=snapshot_start+30 #ran out of memory after 646
  
  
  #Loading the sample cluster to be tracked and sorting its id and id_child
  cluster_group="_snapshot_"+str(snapshot_start) #Remember to change it if you  are changing the star cluster you are tracking in given snapshot
  
  path="./fire2_data_pkl/" #creating a path to store the data only if it does not exist
  #if not os.path.exists(path):
  #  os.makedirs(path)
  
  
  
  
  ###############################################
  ###############################################
  #loading data of all clusters (id, id_children and all)
  
  cluster_path="./fire2_data_pkl/"
  cluster_file_name=simtype+"_clusters_"+simname+cluster_group+".pkl"
  with open(cluster_path+cluster_file_name, "rb") as fp:
      import_cluster = pickle.load(fp)
  
  '''
  import_cluster is a dictionary which is a collection of dictionaries for each cluster.
  It was exported as follows from the program that tracked the clusters.
  Note keys might have changed. Make a test run by importing the pkl file that contains the cluster info. 
  
  cluster1={"cluster_groupid":grpid1,"no_of_star":nstar1,"id":id1,"id_children":id_children1}
  cluster2={"cluster_groupid":grpid2,"no_of_star":nstar2,"id":id2,"id_children":id_children2}
  cluster3={"cluster_groupid":grpid3,"no_of_star":nstar3,"id":id3,"id_children":id_children3} and so on
  export_cluster={1:cluster1,2:cluster2,3:cluster3,4:cluster4,5:cluster5,6:cluster6,7:cluster7,8:cluster8,9:cluster9,10:cluster10}
  path="./"
  file_name="clusters_"+simname+"_snapshot_"+str(snapshot) 
  with open(path+file_name, 'wb') as output:
      pickle.dump(export_cluster, output)
  
  #####3
  To access id from say cluster 2, use import_cluster[2]["id"]
  '''
  ################################################################
  ################################################################
  
  
  
  
  
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
  
  
  total_clusters=len(import_cluster)         #count the total no. of clusters we have in clusters file
  total_snaps=snapshot_end-snapshot_start+1  #count the total no. of snapshots we are going to track
  snap=snapshot_start                        #Mark the beginning of the snapshot
  test_cluster=1                             #Mark the beginning of the cluster that we would be tracking
  ind_tracked={}                             #finding the indices of the tracked stars in each snapshot for each cluster.
  
  tracked_data_all_clusters_each_snap={}     #dictionary for storing data from all clusters and in each snapshot
  
  
  #Algorithm followed for scanning each snapshot for all given star clusters is as follows
  '''
  Step1: Load the star clusters from the given file
  Step2: snap=snapshot to begin with
  Step3: Load id, id_child, id_genearation, age, mass, positoins, velocities and everything you want form the snapshot snap
  Step4: Now assign test_cluster=1 which is the star cluster to begin with for tracking
  Step5: Track the stars from the test_cluster using the id, id_child and id_generation
  Step6: Get the information of the tracked stars id, id_child, id_genearation, age, mass, positoins, velocities and everything you want form the snapshot snap
          Save the information you gathered into a dictionary tracked_data_all_clusters
          test_cluster+=1
  Step7: Is test_cluster<total_clusters? If yes goto step 5 else goto step 8
  Step8: All clusters were tracked for the given snapshot
          Store the data gathered from all clusters in given snapshot to a dictionary tracked_data_all_clusters_each_snap
          It has the format {596:{1:{"x_tracked":[array],...},2:{"x_tracked":[array],...}},597:{1:{"x_tracked":[array],..},2:{"x_tracked":[array],..}..}}
          snap+=1
          Is snap<total_snaps? If yes goto step 3 else goto Step 9
  Step9: End
  
  '''
  
  
  
  
  for i in range(2):               #we run a for loop for just two snapshots
    part=gizmo.io.Read.read_snapshots(['star'],'snapshot_index', snap, simulation_name=simname, simulation_directory=simdir, assign_hosts=True, assign_hosts_rotation=True)               #snap is the snapshot number here that changes everytime the loop iterates. It starts with sanpshot_start
    
    id=part['star'].prop('id')
    id_child=part['star'].prop('id.child')
    id_generation=part['star'].prop('id.generation')
    age=part['star'].prop('age')
    feh = part['star'].prop('metallicity.fe')
    mgh = part['star'].prop('metallicity.mg')
    x=part['star'].prop('host.distance.principal')[:,0] #x component of the position of all stars 
    y=part['star'].prop('host.distance.principal')[:,1] #y component of the position of all stars
    z=part['star'].prop('host.distance.principal')[:,2] #z component of the position of all stars
    n=len(x) # counting the total no. of star particles
    vx=part['star'].prop('host.velocity.principal')[:,0]
    vy=part['star'].prop('host.velocity.principal')[:,1]
    vz=part['star'].prop('host.velocity.principal')[:,2]
    mass=part['star']['mass'] #mass of all stars in snapshot 691
    # compute 3-D velocity in cylindrical coordinates
    # first value is along the major axes (positive definite)
    # secod value is azimuthal velocity in the plane of the disk (positive definite)
    # third value is vertical velocity wrt the disk (signed)
    vR_cyl=part['star'].prop('host.velocity.principal.cylindrical')[:,0]
    vphi_cyl=part['star'].prop('host.velocity.principal.cylindrical')[:,1]
    vz_cyl=part['star'].prop('host.velocity.principal.cylindrical')[:,2]
    print("\n#########\n############\nLoaded id,id_child,age,x,y,z,vx,vy,vz,mass, cylindrical velocities and number of particles for snapshot no.",snap)
    print("Total no of particles in this shapshot no.",snap,"is",n)
    
    
    
    
    
    ################################################################
    ################################################################
    ### Now matching the IDs of this snapshot with the ids of all clusters
    test_cluster=1                       #This is the test cluster to begin with. It resets to 1 at each snapshot 
    tracked_data_all_clusters={}         #This dictionary keeps the tracked information for all star clusters in given snapshot. It also resets with snapshots.
    for j in range(total_clusters):
      tracked_data={}                    #This dictionary holds the tracked information for a cluster temporarily until it is pushed to tracked_data_all_clusters
      print("\n\nNow Matching the ids of the snapshot",snap," with the cluster group id",test_cluster)
      id_test_cluster=import_cluster[test_cluster]["id"]
      id_child_test_cluster=import_cluster[test_cluster]["id_children"]
      sortind=np.argsort(id_test_cluster)
      id_test_cluster_sorted=id_test_cluster[sortind]
      id_child_test_cluster_sorted=id_child_test_cluster[sortind]
      print("The total no. of stars in cluster id ",test_cluster," is",len(id_test_cluster_sorted))
      print("Sorted ids of this cluster is",id_test_cluster_sorted)
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
      vR_cyl_tracked=vR_cyl[ind_tracked]
      vphi_cyl_tracked=vphi_cyl[ind_tracked]
      vz_cyl_tracked=vz_cyl[ind_tracked]
      feh_tracked=feh[ind_tracked]
      mgh_tracked=mgh[ind_tracked]
      
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
      tracked_data={"ind_tracked":ind_tracked,"age_tracked":age_tracked,"x_tracked":x_tracked,"y_tracked":y_tracked,"z_tracked":z_tracked,
      "vx_tracked":vx_tracked,"vy_tracked":vy_tracked,"vz_tracked":vz_tracked,
      "mass_tracked":mass_tracked,"xcm":xcm,"ycm":ycm,"zcm":zcm,"delta_rxyz":delta_rxyz,"rmax":rmax,
      "ymax":ymax,"ymin":ymin,"xmax":xmax,"xmin":xmin,"avg_delta_rxyz":avg_delta_rxyz,
      "vR_cyl_tracked":vR_cyl_tracked,"vphi_cyl_tracked":vphi_cyl_tracked,"vz_cyl_tracked":vz_cyl_tracked,"feh_tracked":feh_tracked,"mgh_tracked":mgh_tracked}
      
      
      tracked_data_all_clusters.update({test_cluster:tracked_data}) #access it using tracked_data_all_clusters[clusterid]["key"]
      test_cluster+=1
    
    #We have collected all tracked information for our test clusters for one snapshot now which is in dictionary tracked_data_all_clusters.
    #We would store that dictionary in a pickle file as a backup data from each snapshot about or all star cluster.
    ###Now moving out of the j loop. We scanned all clusters and stored the tracked data into a dictionary tracked_data_all_clusters  
    file_name=str(snapshot_start)+"clusters_at_snapshot_"+str(snap)+".pkl"
    with open(path+file_name, 'wb') as output:
      pickle.dump(tracked_data_all_clusters, output)
    print("\n Stored tracked data of all stars clusters in the snapshot no.",snap,"to filename:",file_name,"\n#####\n")
    #Note we are still in i loop which is scanning each snapshot !!!!!!!!!!!!!!!!!!
    tracked_data_all_clusters_each_snap.update({snap:tracked_data_all_clusters}) 
    snap=snap+30
  
  #Came out of i loop 
''' 
#Finally we scanned through all the snapshots and now we come out of the loop to store the final file that contains all the tracked information of all clusters 
#Now we store the tracked data of all clusters into a dictionary with snapshot number as the key. This would be our final dictionary of dictinaries. 
with open(path+"total_data_all_clusters_all_snapshots.pkl", 'wb') as output:
  pickle.dump(tracked_data_all_clusters_each_snap, output) #access the data using tracked_data_all_clusters_each_snap[snapshot][cluster_id]["x_tracked"])
  
#Expected Result: tracked_data_all_clusters_each_snap={596:{1:tracked_data,2:tracked_data..},597:{1:tracked_data,2:tracked_data..upto total clusters},598....     upto total snapshots}
  
'''
  
  