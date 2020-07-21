from __main__ import *
from sl_utilities import distinct_colours as dc
from sl_utilities import distance_functions
import utilities as ut

###########Loading the sample cluster to be tracked and sorting its id and id_child
id_test_cluster=np.array([68937285, 22084940, 62548983, 9584721, 19068644, 15790620, 11621407, 18313194, 64000598, 16844755, 61271023, 26250753, 8928920, 56087355, 5936263]) #ids of the cluster to begin with 
id_child_test_cluster=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
sortind=np.argsort(id_test_cluster)
id_test_cluster_sorted=id_test_cluster[sortind]
id_child_test_cluster_sorted=id_child_test_cluster[sortind]
print("The total no. of stars in this cluster is",len(id_test_cluster_sorted))
print("Sorted ids of this cluster is",id_test_cluster_sorted)

##########################################################################################

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
###################################################################


### Now matching the IDs for all our snapshots that are loaded
ind_tracked={} #finding the indices of the tracked stars in each snapshot. To access indices for each snapshot use id[snapshot_n0][ind_tracked[snapshot_no]]
for j in range(snapshot_end+1): 
  if j<snapshot_start:
    ind_tracked[j]=0
  else:
    id_next=id[j]
    id_child_next=id_child[j]
    print("\n\nNow Matching the ids for the snapshot",j,"########\n")
    ind_tracked[j]=matchids(id_test_cluster_sorted,id_child_test_cluster_sorted,id_next,id_child_next)
    
 
###Testing if the matching worked
print("\nTesting:(You should get the same IDs as that of cluster your are tracking.)These are the ids that matached in snapshot 671",id[671][ind_tracked[671]]) 
##################################################################

##################################################################
age_tracked={}
x_tracked={}
y_tracked={}
z_tracked={}
mass_tracked={}
for i in range(snapshot_end+1): # finding the x, y, z and mass of the tracked stars in each snapshot and storing each value in a list
  if i<snapshot_start:
    age_tracked[i]=0
    x_tracked[i]=0
    y_tracked[i]=0
    z_tracked[i]=0
    mass_tracked[i]=0
    
  else:
    age_tracked[i]=age[i][ind_tracked[i]]
    x_tracked[i]=x[i][ind_tracked[i]]
    y_tracked[i]=y[i][ind_tracked[i]]
    z_tracked[i]=z[i][ind_tracked[i]]
    mass_tracked[i]=mass[i][ind_tracked[i]]
  
################################################################
#getting distinct color for each star particle
#colors = dc.get_distinct(len(x_tracked[snapshot_start])) #note it works only for upto 12 particles

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
for i in range(snapshot_end+1): #Calculating the center of mass and related features using the x y and z and mass values from tracked stars of all snapshots
  if i<snapshot_start:
    xcm[i]={}
    ycm[i]={}
    zcm[i]={}
    delta_rxyz[i]={}
    rmax[i]={}
    ymax[i]={}
    ymin[i]={}
    xmax[i]={}
    xmin[i]={}
    avg_delta_rxyz[i]={}
  else:
    xcm[i]=distance_functions.cm(x_tracked[i],mass_tracked[i])
    ycm[i]=distance_functions.cm(y_tracked[i],mass_tracked[i])
    zcm[i]=distance_functions.cm(z_tracked[i],mass_tracked[i])
    delta_rxyz[i]=distance_functions.dr(x_tracked[i],y_tracked[i],z_tracked[i],mass_tracked[i])
    rmax[i]=distance_functions.drmax(x_tracked[i],y_tracked[i],z_tracked[i],mass_tracked[i])
    ymax[i]=(ycm[i]+3*rmax[i])
    ymin[i]=(ycm[i]-3*rmax[i])
    xmax[i]=(xcm[i]+3*rmax[i])
    xmin[i]=(xcm[i]-3*rmax[i])
    avg_delta_rxyz[i]=np.mean(np.absolute(delta_rxyz[i]))

total_subplots=snapshot_end-snapshot_start+1
print("\n Now we are going to plot these snapshots !!!!!! \n####\nTotal plots we would need is",total_subplots)
#cols=int(total_subplots**0.5)
cols=3
print("Total columns we need in this plot is",cols)
rows=total_subplots//cols
rows=rows+total_subplots%cols
position = range(1,total_subplots + 1)

fig = plt.figure(figsize=(10,20))
snap=snapshot_start
for i in range(total_subplots): # add every single subplot to the figure with a for loop
    print("\n\nx of snapshot",snap,"is",x_tracked[snap])
    print("xcm of snapshot",snap,"is",xcm[snap])
    ax = fig.add_subplot(rows,cols,position[i])
    ax.scatter(np.absolute(x_tracked[snap]),np.absolute(y_tracked[snap]),marker=".",s=1)
    ax.plot(np.absolute(xcm[snap]),np.absolute(ycm[snap]),color='red',marker=".",markersize=1)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.minorticks_on()
    ax.set_xlim(np.absolute(xmin[snap]),np.absolute(xmax[snap]))
    ax.set_ylim(np.absolute(ymin[snap]),np.absolute(ymax[snap]))
    title="Snapshot " + str(snap)
    ax.set_title(title)
    snap=snap+1
      
plt.tight_layout()
plotname="./plots/snapshots_"+str(snapshot_start)+"_to_"+str(snapshot_end)
plt.savefig(plotname,dpi=600)
plt.close()
print("Plot generated and saved as filename:",plotname)

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
fig2.savefig("./plots/snapshots_vs_averageDistanceFromCM")
plt.close()

#Now lets write our tracked data from each snapshots to a file

count=snapshot_start
for i in range(len(snapshot_list)):
  dict_exportdata={"ind_tracked":ind_tracked[count],
"age_tracked":age_tracked[count],"x_tracked":x_tracked[count],"y_tracked":y_tracked[count],"z_tracked":z_tracked[count],"mass_tracked":mass_tracked[count],"xcm":xcm[count],"ycm":ycm[count],"zcm":zcm[count],
"delta_rxyz":delta_rxyz[count],"rmax":rmax[count],"ymax":ymax[count],"ymin":ymin[count],"xmax":xmax[count],"xmin":xmin[count],"avg_delta_rxyz":avg_delta_rxyz[count]}
  file_name="export_tracked_data_snapshot"+str(count)
  path="./data/"
  ut.io.file_hdf5(path+file_name, dict_exportdata)
  print("\n Stored data from the snapshot no.",count,"to filename:",file_name,"\n#####\n")
  count=count+1




'''
# below, cluster4 is a dictionary i generate that i want to write to an hdf5 file
path = '/home1/04712/tg840119/' # where i want to write the output file (my home directory on stampede)
file_name = 'm12f_cluster4.hdf5' # the name of the file i want to write
ut.io.file_hdf5(path+file_name, cluster4) # call the file_hdf5 function from utilities.io
# at this point i have an hdf5 file containing my dictionary saved to the directory i specified
# to read this file in later, just type: 
cluster4 = ut.io.file_hdf5(path+file_name) # now i have the dictionary i want loaded into cluster4 again.
'''


'''
#Now plotting ind 2D the stars in two snapshots
fig8 = plt.figure()
#fig8.suptitle("Stars aged 7 to 8 dec within 1 from 0,8,0 \n \n")

ax1 = fig8.add_subplot(321)
for i in range(len(x_691_tracked)):
  ax1.scatter([-1*x_691_tracked[i]],[-1*y_691_tracked[i]],color=colors[i],marker=".",s=10)
ax1.plot(xcm_691,ycm_691,color='red',marker=".")
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.minorticks_on()
ax1.set_xlim(xmin_691,xmax_691)
ax1.set_ylim(ymin_691,ymax_691)
ax1.set_title('Snapshot 691')

ax2=fig8.add_subplot(322)
for i in range(len(x_692_tracked)):
  ax2.scatter([x_692_tracked[i]],[y_692_tracked[i]],color=colors[i],marker=".",s=10)
ax2.plot(xcm_692,ycm_692,color='red',marker=".")
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.minorticks_on()
ax2.set_xlim(xmin_692,xmax_692)
ax2.set_ylim(ymin_692,ymax_692)
ax2.set_title('Snapshot 692')


ax3=fig8.add_subplot(323)
for i in range(len(x_693_tracked)):
  ax3.scatter([x_693_tracked[i]],[y_693_tracked[i]],color=colors[i],marker=".",s=10)
ax3.plot(xcm_693,ycm_693,color='red',marker=".")
ax3.set_xlabel('x')
ax3.set_ylabel('y')
ax3.set_xlim(xmin_693,xmax_693)
ax3.set_ylim(ymin_693,ymax_693)
ax3.minorticks_on()
ax3.set_title('Snapshot 693')

#plt.subplots_adjust(hspace=.5)

ax4=fig8.add_subplot(324)
for i in range(len(x_694_tracked)):
  ax4.scatter([x_694_tracked[i]],[y_694_tracked[i]],color=colors[i],marker=".",s=10)
ax4.plot(xcm_694,ycm_694,color='red',marker=".")
ax4.set_xlabel('x')
ax4.set_ylabel('y')
ax4.minorticks_on()
ax4.set_xlim(xmin_694,xmax_694)
ax4.set_ylim(ymin_694,ymax_694)
ax4.set_title('Snapshot 694')


ax5=fig8.add_subplot(325)
for i in range(len(x_695_tracked)):  
  ax5.scatter([x_695_tracked[i]],[y_695_tracked[i]],color=colors[i],marker=".",s=10)
ax5.plot(xcm_695,ycm_695,color='red',marker=".")
ax5.set_xlabel('x')
ax5.set_ylabel('y')
ax5.minorticks_on()
ax5.set_xlim(xmin_695,xmax_695)
ax5.set_ylim(ymin_695,ymax_695)
ax5.set_title('Snapshot 695')

#plt.subplots_adjust(hspace=.5)

ax6=fig8.add_subplot(326)
for i in range(len(x_696_tracked)):  
  ax6.scatter([x_696_tracked[i]],[y_696_tracked[i]],color=colors[i],marker=".",s=10)
ax6.plot(xcm_696,ycm_696,color='red',marker=".")
ax6.set_xlabel('x')
ax6.set_ylabel('y')
ax6.minorticks_on()
ax6.set_xlim(xmin_696,xmax_696)
ax6.set_ylim(ymin_696,ymax_696)
ax6.set_title('Snapshot 696')

plt.tight_layout()
fig8.savefig("./plots/cluster_group_15_snapshots691to696.png")
'''