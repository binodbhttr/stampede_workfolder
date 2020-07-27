#from __main__ import *
from sl_utilities import distinct_colours as dc
from sl_utilities import distance_functions
import utilities as ut
from matplotlib import pyplot as plt
import numpy as np
import os
#this file is targetted to help avoid running tonns of snapshots each time to analyse the tracked data
# I am working towards loading each data from this file and generating all the kinds of plots here
# It is assumed we know the snapshot_start and snapshot_end value
################################################
#Now loading data from each file and storing it to a dictionary importdata which can be accessed as importdata[snapshotnumber]['key']
#For example to load the x coordinates of the tracked stars in snapshot no. 690, importdata[690]['x_tracked']
cluster_groupid="snapshot596_cluster_group8" #Remember to change it if you  are changing the star culuster you are tracking in given snapshot
data_path="./data/"+cluster_groupid+"/" #this is the path of the directory where our data is present

plot_path="./plots/"+cluster_groupid+"/" #creating a path to store the plots only if it does not exist
if not os.path.exists(plot_path):
  os.makedirs(plot_path)
  

snapshot_start=596
snapshot_end=696
n=snapshot_end-snapshot_start+1
count=snapshot_start
importdata={}
for i in range(n):
  file_name=cluster_groupid+"_tracked_data_snapshot"+str(count)
  importdata[count]=ut.io.file_hdf5(data_path+file_name) #reading data from each file and storing it to a dictionary importdata
  os.system('clear')
  print("\n Loaded data from the snapshot no.",count,"located in filename:",file_name,"to importdata[",count,"]\n#####\n")
  count=count+1
################################################
colors=['cyan','blue','green','magenta','yellow','orange','purple','tan','lime','brown','grey','pink','navy','teal']


##Now plotting all snapshots
#################################################  
#total_subplots=snapshot_end-snapshot_start+1  
#print("\n Now we are going to plot these snapshots !!!!!! \n####\nTotal plots we would need is",total_subplots)
##cols=int(total_subplots**0.5)
#cols=3
#print("Total columns we need in this plot is",cols)
#rows=total_subplots//cols
#rows=rows+total_subplots%cols
#position = range(1,total_subplots + 1)
#
#fig1 = plt.figure(figsize=(9,rows*3))
#snap=snapshot_start
#for i in range(total_subplots): # add every single subplot to the figure with a for loop
#    print("\n\nx, y and z of snapshot",snap,"is",importdata[snap]['x_tracked'],importdata[snap]['y_tracked'],importdata[snap]['z_tracked'])
#    print("xcm of snapshot",snap,"is",importdata[snap]['xcm'],importdata[snap]['ycm'],importdata[snap]['zcm'])
#    os.system('clear')
#    ax = fig1.add_subplot(rows,cols,position[i])
#    ax.scatter(np.absolute(importdata[snap]['x_tracked']),np.absolute(importdata[snap]['y_tracked']),marker=".",s=importdata[snap]['mass_tracked']/100,c=colors)
#    ax.plot(np.absolute(importdata[snap]['xcm']),np.absolute(importdata[snap]['ycm']),color='red',marker=".",markersize=1)
#    ax.set_xlabel('x')
#    ax.set_ylabel('y')
#    ax.minorticks_on()
#    ax.set_xlim(np.absolute(importdata[snap]['xmin']),np.absolute(importdata[snap]['xmax']))
#    ax.set_ylim(np.absolute(importdata[snap]['ymin']),np.absolute(importdata[snap]['ymax']))
#    title="Snapshot " + str(snap)
#    ax.set_title(title)
#    snap=snap+1
#      
#plt.tight_layout()
#plotname=cluster_groupid+"_"+str(snapshot_start)+"_to_"+str(snapshot_end)+"_using_exported_data.png"
#print("###################\nSaving plot of all snapshots as filename:",plotname)
#fig1.savefig(plot_path+plotname,dpi=600)
#print("\nSaving Complete !!!")
#plt.close()
###############################################################################
###############################################################################



##################################
#Now plotting the verage distances from the center of mass for each snapshots we tracked
snapshot_list=np.arange(snapshot_start,snapshot_end+1) #create a list of snapshots to plot with the average distance from the CM
avg_r_cm_temp=np.array(0)
count=snapshot_start
for i in range(len(snapshot_list)):
  avg_r_cm_temp=np.append(avg_r_cm_temp,importdata[count]['avg_delta_rxyz'])
  count=count+1
avg_r_cm=avg_r_cm_temp[1:len(avg_r_cm_temp)]
print("\nThese are the average distances from the center of mass for each snapshots we tracked:\n",avg_r_cm)


fig2=plt.figure()
ax1=fig2.add_subplot(1,1,1)
ax1.plot(snapshot_list-snapshot_start,avg_r_cm*1000, color='b') #snapshot_list and avg_r_cm both are arrays
ax1.set_xlabel('Time in Myr')
ax1.set_ylabel('Average Distance of stars from CM of the cluster in parsec')
ax1.set_title(cluster_groupid)
ax1.minorticks_on()
plt.tight_layout()
plotname="averageDistanceFromCM_"+str(snapshot_start)+"_to_"+str(snapshot_end)+".png"
print("###################\nSaving plot of Average Distance from CM vs Snapshots as:",plotname)
fig2.savefig(plot_path+plotname)
plt.close()
###################################