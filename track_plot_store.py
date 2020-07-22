from __main__ import *
from sl_utilities import distinct_colours as dc
from sl_utilities import distance_functions
import utilities as ut
from matplotlib import pyplot as plt
import numpy as np


#run this file after running the the file "load_snapshots_start_to_end.py"
#this file allows you to give the start cluster in the snapshot to begin with and allow you to plot them across all other snapshots in future 
# and export the data for future use
################################################################
################################################################
#Loading the sample cluster to be tracked and sorting its id and id_child
cluster_groupid="snapshot671_cluster_group3" #Remember to change it if you  are changing the star culuster you are tracking in given snapshot
id_test_cluster=np.array([68937285, 22084940, 62548983, 9584721, 19068644, 15790620, 11621407, 18313194, 64000598, 16844755, 61271023, 26250753, 8928920, 56087355, 5936263]) #ids of the cluster to begin with 
id_child_test_cluster=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
sortind=np.argsort(id_test_cluster)
id_test_cluster_sorted=id_test_cluster[sortind]
id_child_test_cluster_sorted=id_child_test_cluster[sortind]
print("The total no. of stars in this cluster is",len(id_test_cluster_sorted))
print("Sorted ids of this cluster is",id_test_cluster_sorted)
################################################################
################################################################



################################################################
################################################################
###Now let's find the matching ids in the next snapshot using a function
def matchids(id_current,id_child_current,id_next,id_child_next): #this function returns the index of the ids in the next snapshot that match with the current
  ind=np.array(0)
  for i in range(len(id_current)):
    match=np.where((id_next==id_current[i])&(id_child_next==id_child_current[i]))
    print("\nFound the id",id_next[match],"at the index",match[0],"in this snapshot")
    ind=np.append(ind,match)
  ind_tracked_id_next=ind[1:len(ind)]  #The extra element in the beginning is removed by this process
  print("\nThese are the indices of the ids that matched in current snapshot\n",ind_tracked_id_next)
  return ind_tracked_id_next  
################################################################
################################################################





################################################################
################################################################
### Now matching the IDs for all our snapshots that are loaded
ind_tracked={} #finding the indices of the tracked stars in each snapshot. To access indices for each snapshot use id[snapshot_n0][ind_tracked[snapshot_no]]
total_snaps=snapshot_end-snapshot_start+1
snap=snapshot_start
for j in range(total_snaps): 
  id_next=id[snap]
  id_child_next=id_child[snap]
  print("\n\nNow Matching the ids for the snapshot",snap,"########\n")
  ind_tracked[snap]=matchids(id_test_cluster_sorted,id_child_test_cluster_sorted,id_next,id_child_next)
  snap=snap+1
    
 
###Testing if the matching worked
print("\nTesting:(You should get the same IDs as that of cluster you are tracking.)These are the ids that matached in snapshot 671",id[671][ind_tracked[671]]) 
################################################################
################################################################






################################################################
################################################################
#Now finding the postion, velocities, age and mass of the tracked stars in each snapshots
age_tracked={}
x_tracked={}
y_tracked={}
z_tracked={}
vx_tracked={}
vy_tracked={}
vz_tracked={}
mass_tracked={}
total_snaps=snapshot_end-snapshot_start+1
snap=snapshot_start
for i in range(total_snaps): # finding the x, y, z and mass of the tracked stars in each snapshot and storing each value in a list
  age_tracked[snap]=age[snap][ind_tracked[snap]]
  x_tracked[snap]=x[snap][ind_tracked[snap]]
  y_tracked[snap]=y[snap][ind_tracked[snap]]
  z_tracked[snap]=z[snap][ind_tracked[snap]]
  vx_tracked[snap]=vx[snap][ind_tracked[snap]]
  vy_tracked[snap]=vy[snap][ind_tracked[snap]]
  vz_tracked[snap]=vz[snap][ind_tracked[snap]]
  mass_tracked[snap]=mass[snap][ind_tracked[snap]]
  snap=snap+1
################################################################
################################################################




################################################################
################################################################  
#getting distinct color for each star particle
#colors = dc.get_distinct(len(x_tracked[snapshot_start])) #note it works only for upto 12 particles
################################################################
################################################################





################################################################
################################################################
#Now calculating center of mass and realted properties
xcm={}
ycm={}
zcm={}
delta_rxyz={}
rmax={}
ymax={}
ymin={}
xmax={}
xmin={}
avg_delta_rxyz={}
total_snaps=snapshot_end-snapshot_start+1
snap=snapshot_start
for i in range(total_snaps): #Calculating the center of mass and related features using the x y and z and mass values from tracked stars of all snapshots
  xcm[snap]=distance_functions.cm(x_tracked[snap],mass_tracked[snap])
  ycm[snap]=distance_functions.cm(y_tracked[snap],mass_tracked[snap])
  zcm[snap]=distance_functions.cm(z_tracked[snap],mass_tracked[snap])
  delta_rxyz[snap]=distance_functions.dr(x_tracked[snap],y_tracked[snap],z_tracked[snap],mass_tracked[snap])
  rmax[snap]=distance_functions.drmax(x_tracked[snap],y_tracked[snap],z_tracked[snap],mass_tracked[snap])
  ymax[snap]=(ycm[snap]+3*rmax[snap])
  ymin[snap]=(ycm[snap]-3*rmax[snap])
  xmax[snap]=(xcm[snap]+3*rmax[snap])
  xmin[snap]=(xcm[snap]-3*rmax[snap])
  avg_delta_rxyz[snap]=np.mean(np.absolute(delta_rxyz[snap]))
  snap=snap+1
################################################################
################################################################




################################################################
################################################################
#Now going to plot those snapshots with tracked data from the given star cluster
total_subplots=snapshot_end-snapshot_start+1
print("\n Now we are going to plot these snapshots !!!!!! \n####\nTotal plots we would need is",total_subplots)
#cols=int(total_subplots**0.5) #It was a bad idea
cols=3
print("Total columns we need in this plot is",cols)
rows=total_subplots//cols
rows=rows+total_subplots%cols
position = range(1,total_subplots + 1)

fig = plt.figure(figsize=(10,20))
snap=snapshot_start
for i in range(total_subplots): # add every single subplot to the figure with a for loop
    print("\n\n x of snapshot",snap,"is",x_tracked[snap])
    print("xcm of snapshot",snap,"is",xcm[snap])
    ax = fig.add_subplot(rows,cols,position[i])
    ax.scatter(np.absolute(x_tracked[snap]),np.absolute(y_tracked[snap]),marker=".",s=1)
    ax.plot(np.absolute(xcm[snap]),np.absolute(ycm[snap]),color='red',marker=".",markersize=1)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    
    ax.minorticks_on()
    ax.set_xlim(np.absolute(xmin[snap]),np.absolute(xmax[snap]))
    ax.set_ylim(np.absolute(ymin[snap]),np.absolute(ymax[snap]))
    print("test")
    title="Snapshot " + str(snap)
    ax.set_title(title)
    print("test")
    snap=snap+1
      
plt.tight_layout()
plotname="./plots/snapshots_"+str(snapshot_start)+"_to_"+str(snapshot_end)+cluster_groupid+".png"
plt.savefig(plotname,dpi=600)
plt.close()
print("Plot generated and saved as filename:",plotname)
################################################################
################################################################





################################################################
################################################################
#plotting the average distance of paricles from the center of mass of each snapshot
snapshot_list=np.arange(snapshot_start,snapshot_end+1)
avg_r_cm_temp=np.array(0)
count=snapshot_start
for i in range(len(snapshot_list)):
  avg_r_cm_temp=np.append(avg_r_cm_temp,avg_delta_rxyz[count])
  count=count+1

avg_r_cm=avg_r_cm_temp[1:len(avg_r_cm_temp)]
print(avg_r_cm)


fig2=plt.figure()
ax1=fig2.add_subplot(1,1,1)
ax1.plot(snapshot_list,avg_r_cm,marker='*', color='b')
ax1.set_xlabel('Snapshots')
ax1.set_ylabel('Average Distance from CM')
ax1.minorticks_on()
plt.tight_layout()
plotname="./plots/averageDistanceFromCM_"+str(snapshot_start)+"_to_"+str(snapshot_end)
fig2.savefig(plotname)
plt.close()
################################################################
################################################################






################################################################
################################################################
#Now lets write our tracked data from each snapshots to a file
total_snaps=snapshot_end-snapshot_start+1
count=snapshot_start
for i in range(total_snaps):
  dict_exportdata={"ind_tracked":ind_tracked[count],
"age_tracked":age_tracked[count],"x_tracked":x_tracked[count],"y_tracked":y_tracked[count],"z_tracked":z_tracked[count],"vx_tracked":vx_tracked[count],"vy_tracked":vy_tracked[count],"vz_tracked":vz_tracked[count],"mass_tracked":mass_tracked[count],"xcm":xcm[count],"ycm":ycm[count],"zcm":zcm[count],
"delta_rxyz":delta_rxyz[count],"rmax":rmax[count],"ymax":ymax[count],"ymin":ymin[count],"xmax":xmax[count],"xmin":xmin[count],"avg_delta_rxyz":avg_delta_rxyz[count]}
  file_name=cluster_groupid+"_tracked_data_snapshot"+str(count)
  path="./data/"
  ut.io.file_hdf5(path+file_name, dict_exportdata)
  print("\n Stored data from the snapshot no.",count,"to filename:",file_name,".hdf5\n#####\n")
  count=count+1  


################################################################
################################################################  