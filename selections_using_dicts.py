from __main__ import *
from sl_utilities import distinct_colours as dc
from sl_utilities import distance_functions

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
    print("\nMatched the id",id_next[match])
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
    #print(id_array[j])
    #print(id_child_array[j])
    id_next=id[j]
    id_child_next=id_child[j]
    ind_tracked[j]=matchids(id_test_cluster_sorted,id_child_test_cluster_sorted,id_next,id_child_next)
    
 
###Testing if the matching worked
print("These are the ids that were tracked in snapshot 671",id[671][ind_tracked[671]]) 
##################################################################


'''
##################################################################
age_tracked=[]
x_tracked=[]
y_tracked=[]
z_tracked=[]
mass_tracked=[]
for i in range(snapshot_end+1): # finding the x, y, z and mass of the tracked stars in each snapshot and storing each value in a list
  x_tracked.append(x[i][ind_tracked[i]])
  y_tracked.append(y[i][ind_tracked[i]])
  z_tracked.append(z[i][ind_tracked[i]])
  mass_tracked.append(mass_tracked[i][ind_tracked[i]])
  
################################################################

#detting distinct color for each star particle
colors = dc.get_distinct(len(x_tracked[snapshot_start]))


#we had these values from the data 




xcm_691=8.2097
ycm_691=-0.9572
ycm_691=ycm_691*-1
xcm_691=xcm_691*-1
rmax_691=11.4
rmax_kpc_691=0.0114

ymax_691=ycm_691+3*rmax_kpc_691
ymin_691=ycm_691-3*rmax_kpc_691
xmax_691=xcm_691+3*rmax_kpc_691
xmin_691=xcm_691-3*rmax_kpc_691


xcm=[]
ycm=[]
zcm=[]
delta_rxyz=[]
rmax=[]
ymax=[]
ymin=[]
xmax=[]
xmin=[]

for i in range(snapshot_end+1): #Calculating the center of mass and related features using the x y and z and mass values from tracked stars of all snapshots
  xcm.append(distance_functions.cm(x__tracked[i],mass_tracked[i]))
  ycm.append(distance_functions.cm(y__tracked[i],mass_tracked[i]))
  zcm.append(distance_functions.cm(z__tracked[i],mass_tracked[i]))
  delta_rxyz.append(distance_functions.dr(x_tracked[i],y_tracked[i],z_tracked[i],mass_tracked[i]))
  rmax.append(distance_functions.drmax(x_tracked[i],y_tracked[i],z_tracked[i],mass_tracked[i]))
  ymax.append(ycm[i]+3*rmax[i])
  ymin.append(ycm[i]-3*rmax[i])
  xmax.append(xcm[i]+3*rmax[i])
  xmin.append(xcm[i]-3*rmax[i])

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