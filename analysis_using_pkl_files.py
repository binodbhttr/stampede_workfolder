import utilities as ut
from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from sl_utilities import distinct_colours as dc
from sl_utilities import distance_functions
import os
import pickle


path="./data_pkl/" 
snap=596
simname = 'm12i_res7100_mhdcv'
simdir = '/scratch/projects/xsede/GalaxiesOnFIRE/mhdcv/m12i_res7100_mhdcv/1Myr/1Myr_fire2'
file_name="total_data_all_clusters_all_snapshots.pkl"
snapshot_start=596
snapshot_end=696



with open(path+file_name, "rb") as input:
    importdata = pickle.load(input)

print("Testing if the loading of data was successful !! \n")
print("x_tracked of snapshot 596 for cluster 1",importdata[596][1]["x_tracked"])
print("x_cm of snapshot 596 for cluster 1",importdata[596][1]["xcm"])
print("x_tracked of snapshot 597 for cluster 1",importdata[597][1]["x_tracked"])
print("x_cm of snapshot 597 for cluster 1",importdata[597][1]["xcm"])
print("x_tracked of snapshot 597 for cluster 2",importdata[597][2]["x_tracked"])
print("x_cm of snapshot 597 for cluster 2",importdata[597][2]["xcm"])



cluster_groupid=[]
total_clusters=len(importdata[snapshot_start])
for i in range(total_clusters):
    cluster_groupid.append("snapshot"+str(snapshot_start)+"_cluster_group"+str(i+1))

print("These are the clusters groups we have tracked data of:\n",cluster_groupid)




plot_path="./plots_pkl/allculsters_single_figure/" #creating a path to store the plots only if it does not exist
if not os.path.exists(plot_path):
  os.makedirs(plot_path)


##################################
#Now plotting the average distances from the center of mass for all snapshots we tracked in a single figure for all clusters
fig1=plt.figure(figsize=(8,5))
ax1=fig1.add_subplot(1,1,1)
ax1.minorticks_on()
ax1.set_xlabel('Time in Myr')
ax1.set_ylabel('Average Distance of stars from CM of the cluster in parsec')
snapshot_list=np.arange(snapshot_start,snapshot_end+1) #create a list of snapshots to plot
time=snapshot_list-snapshot_start #time starts from zero here where t=0 is at the snapshot where we start from
n=snapshot_end-snapshot_start+1
cluster_count=0
for i in range(len(cluster_groupid)):
    avg_r_cm_temp=np.array(0)
    snapshot_count=snapshot_start
    for j in range(n):
        avg_r_cm_temp=np.append(avg_r_cm_temp,importdata[snapshot_count][cluster_count+1]["avg_delta_rxyz"])
        snapshot_count+=1
    avg_r_cm=avg_r_cm_temp[1:len(avg_r_cm_temp)]*1000 #converted into parsec
    #print(avg_r_cm)
    slope, intercept = np.polyfit(time,avg_r_cm, 1)
    ax1.plot(time,avg_r_cm,label=cluster_groupid[cluster_count]+",slope="+str(round(slope,4)))
    #ax1.plot(time,slope*time+intercept,label=cluster_groupid[cluster_count]+"fitted")
    ax1.legend(loc='upper left')
    cluster_count+=1
plt.tight_layout()
fig1.savefig(plot_path+"avg_delta_rxyz_all_clusters.png",dpi=150)
#########################################################
#########################################################