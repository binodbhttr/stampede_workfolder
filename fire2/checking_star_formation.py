#!/usr/bin/env python
# coding: utf-8
import gizmo_analysis as gizmo
import utilities as ut
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from sl_utilities import distinct_colours as dc
from sl_utilities import distance_functions
import os
import pickle

simname = 'm12i_res7100_mhdcv'
simdir = '/scratch/projects/xsede/GalaxiesOnFIRE/cr_suite/m12i_res7100/mhdcv/1Myr/fire2/'
simtype="fire2"
snapshot_start=596
snapshot_end=696 #ran out of memory after 646
#Loading the sample cluster to be tracked and sorting its id and id_child
cluster_group="snapshot596" #Remember to change it if you  are changing the star cluster you are tracking in given snapshot

data_path="./fire2_data_pkl/" #creating a path to store the data only if it does not exist
gas_datapath="./fire2_gas_data_pkl/"
young_star_datapath="./fire2_young_stars/"
young_star_plotpath="./fire2_young_stars_plot/"

if not os.path.exists(young_star_datapath):
  os.makedirs(young_star_datapath)

if not os.path.exists(young_star_plotpath):
  os.makedirs(young_star_plotpath)

#############################################################################
#constants
#############################################################################

MsunToGm = 1.99e33
KpcToCm = 3.086e21
mp = 1.67e-24
#bin_edge = 10.
bin_edge = 30.

bins = np.arange(-25,25,0.1)

cluster_file_name="clusters_"+simname+"_snapshot_"+str(snapshot_start) 
with open(data_path+cluster_file_name, "rb") as fp:
    import_cluster = pickle.load(fp)

total_clusters=len(import_cluster)         #count the total no. of clusters we have in clusters file
total_snaps=snapshot_end-snapshot_start+1  #count the total no. of snapshots we are going to track
snap=snapshot_start                        #Mark the beginning of the snapshot

tracked_stars_each_snap={}     #dictionary for storing data from all clusters and in each snapshot

colors=['cyan','blue','green','magenta','yellow','teal','brown','darkslategray','lime','red','orange','purple','rosybrown','pink','navy','olive','cornflowerblue']
##################################

snap=snapshot_start
for i in range(total_snaps):               #we run a for loop until the end of all snapshots
  new_stars={}
  fig1=plt.figure()
  fig1.set_size_inches(7,7)
  ax=fig1.add_axes([0.17, 0.185, 0.65, 0.65]) #left, bottom, width, height
  cluster_data_name="all_clusters_at_snapshot_"+str(snap)+".pkl" 
  with open(data_path+cluster_data_name, "rb") as input:
    tracked_data= pickle.load(input)
    
  part=gizmo.io.Read.read_snapshots(['star'],'snapshot_index', snap, simulation_name=simname, simulation_directory=simdir, assign_hosts=True, assign_hosts_rotation=True)               #snap is the snapshot number here that changes everytime the loop iterates. It starts with sanpshot_start

  age=part['star'].prop('age')
  x=part['star'].prop('host.distance.principal')[:,0] #x component of the position of all stars 
  y=part['star'].prop('host.distance.principal')[:,1] #y component of the position of all stars
  z=part['star'].prop('host.distance.principal')[:,2] #z component of the position of all stars
  Rxy = part['star'].prop('host.distance.principal.cylindrical')[:,0]
  
  cluster_count=0                       #This is the test cluster to begin with. It resets to 1 at each snapshot 
  for j in range(total_clusters):
    xcm=tracked_data[cluster_count+1]["xcm"]
    ycm=tracked_data[cluster_count+1]["ycm"]
    x_cluster=tracked_data[cluster_count+1]["x_tracked"]
    y_cluster=tracked_data[cluster_count+1]["y_tracked"]
    ax.scatter(x_cluster,y_cluster,c=colors[j])
    cluster_count+=1
    
  gas_file_name=simtype+"_gas_data"+str(snap)+".pkl"
  with open(gas_datapath+gas_file_name, "rb") as input:
      import_gasdata = pickle.load(input)
  
  face=import_gasdata["face"]
  snaptime=import_gasdata["snaptime"]
  norm = matplotlib.colors.LogNorm(vmin=1, vmax=1000) #the color range plotted
  im = ax.imshow(face/(((bins[1]-bins[0])*1000)**2),origin='lower',interpolation='nearest',norm=norm,extent=(-25,25,-25,25),cmap='binary') 
    
  #colorbar for the background gas density
  cmap_gray = matplotlib.cm.get_cmap('binary')
  norm1 = matplotlib.colors.LogNorm(vmin=1,vmax=1000)
  cbar_ax1 = fig1.add_axes([0.04, 0.185, 0.04, 0.64]) # position of gray colorbar (left, bottom, width, height)
  cb1 = fig1.colorbar(im, cax=cbar_ax1, ticklocation='left')
  cb1.set_label('$\Sigma$ (M$_{{\odot}}$/pc$^2$)', labelpad=-5, fontsize=12)
  
  keep_young = np.where((age <= .003) & (abs(z) < 1.5))
  x0=x[keep_young]
  y0=y[keep_young]
  
  ax.scatter(x0,y0,c="red",s=5,label="Young Stars (<3Myr)")
  ax.set_xlim(np.min(x0)-5,np.max(x0)+5)
  ax.set_ylim(np.min(y0)-5,np.max(y0)+5)
  ax.legend(bbox_to_anchor=(1,0.5), loc='center left')
  
  fig1.savefig(str(snap)+"test.png",bbox_inches='tight',dpi=200)  
  fig1.clf()
  plt.close()  
  ################################################################
  #Now lets write our tracked data from each snapshots to a file
  new_stars.update({"x_new_stars":x0,"y_new_stars":y0}) #access it using tracked_data_all_clusters[clusterid]["key"]
  
  file_name="young_stars_at_snapshot_"+str(snap)+".pkl"
  with open(young_star_datapath+file_name, 'wb') as output:
    pickle.dump(new_stars, output)
  snap=snap+1 