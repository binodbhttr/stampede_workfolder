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
cluster_groupid="snapshot596_cluster_group15" #Remember to change it if you  are changing the star culuster you are tracking in given snapshot
data_path="./data/"+cluster_groupid+"/" #this is the path of the directory where our data is present

plot_path="./plots/"+cluster_groupid+"/" #creating a path to store the plots only if it does not exist
if not os.path.exists(plot_path):
  os.makedirs(plot_path)
  
#dict_keys(['age_tracked', 'avg_delta_rxyz', 'delta_rxyz', 'ind_tracked', 'mass_tracked', 'rmax', 'vx_tracked', 'vy_tracked', 'vz_tracked', 'x_tracked', 'xcm', 'xmax', 'xmin', 'y_tracked', 'ycm', 'ymax', 'ymin', 'z_tracked', 'zcm']) #these are the keys of dictionary importdata[snapshotnumber]


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




##################################
#Now plotting the average distances from the center of mass for each snapshots we tracked
snapshot_list=np.arange(snapshot_start,snapshot_end+1) #create a list of snapshots to plot with the average distance from the CM
avg_r_cm_temp=np.array(0)
count=snapshot_start
for i in range(len(snapshot_list)):
  avg_r_cm_temp=np.append(avg_r_cm_temp,importdata[count]['avg_delta_rxyz'])
  count=count+1
avg_r_cm=avg_r_cm_temp[1:len(avg_r_cm_temp)]
time=snapshot_list-snapshot_start #time starts from zero here where t=0 is at the snapshot where we start from
print("\nThese are the average distances from the center of mass for each snapshots we tracked:\n",avg_r_cm)


fig2=plt.figure()
ax1=fig2.add_subplot(1,1,1)
ax1.plot(time,avg_r_cm*1000, color='b') #snapshot_list and avg_r_cm both are arrays
ax1.set_xlabel('Time in Myr')
ax1.set_ylabel('Average Distance of stars from CM of the cluster in parsec')
ax1.set_title(cluster_groupid)
ax1.minorticks_on()
plt.tight_layout()
plotname="averageDistanceFromCM_"+str(snapshot_start)+"_to_"+str(snapshot_end)+".png"
print("###################\nSaving plot of Average Distance from CM vs Snapshots as:",plotname)
fig2.savefig(plot_path+plotname)
plt.close()




#########################################################
#########################################################
#plotting the standard deviation of the max r from the center of Mass for each snapshot
std_r_cm_temp=np.array(0)
count=snapshot_start
for i in range(len(snapshot_list)):
  std_r_cm_temp=np.append(std_r_cm_temp,np.nanstd(importdata[count]['delta_rxyz']))
  count=count+1
std_r_cm=std_r_cm_temp[1:len(std_r_cm_temp)]
time=snapshot_list-snapshot_start #time starts from zero here where t=0 is at the snapshot where we start from
print("\nThese are the stdev of distances from the center of mass for each snapshots we tracked:\n",std_r_cm)




fig3=plt.figure()
ax1=fig3.add_subplot(1,1,1)
ax1.plot(time,std_r_cm*1000,color='b') #snapshot_list and avg_r_cm both are arrays
ax1.set_xlabel('Time in Myr')
ax1.set_ylabel('stdev of Distance of stars from CM of the cluster in parsec')
ax1.set_title(cluster_groupid)
ax1.minorticks_on()
plt.tight_layout()
plotname="stdev_DistanceFromCM_"+str(snapshot_start)+"_to_"+str(snapshot_end)+".png"
print("###################\nSaving plot of stdev of Distance from CM vs Snapshots as:",plotname)
fig3.savefig(plot_path+plotname)
plt.close()


#########################################################
#########################################################
#plotting the standard deviation of the magnitude of velocitiy of the stars in each snapshot
std_v_mag_temp=np.array(0)
count=snapshot_start
for i in range(len(snapshot_list)):
  v_mag=(importdata[count]['vx_tracked']**2+importdata[count]['vy_tracked']**2+importdata[count]['vz_tracked']**2)**(1/2)*1.023
  std_v_mag_temp=np.append(std_v_mag_temp,np.nanstd(v_mag))
  count=count+1
std_v_mag=std_v_mag_temp[1:len(std_v_mag_temp)]
time=snapshot_list-snapshot_start #time starts from zero here where t=0 is at the snapshot where we start from
print("\nThese are the stdev of distances from the center of mass for each snapshots we tracked:\n",std_v_mag)


fig4=plt.figure()
ax1=fig4.add_subplot(1,1,1)
ax1.plot(time,std_v_mag,color='b') #time and std_v_mag both are arrays
ax1.set_xlabel('Time in Myr')
ax1.set_ylabel('stdev of Velocity of stars in the cluster in parsec/Myr')
ax1.set_title(cluster_groupid)
ax1.minorticks_on()
plt.tight_layout()
plotname="stdev_Velocity_"+str(snapshot_start)+"_to_"+str(snapshot_end)+".png"
print("###################\nSaving plot of stdev of Velocity of stars in the cluster in each Snapshot as:",plotname)
fig4.savefig(plot_path+plotname)
plt.close()



#########################################################
#########################################################
#plotting the largest distance from the Center of Mass of the stars in each snapshot
rmax_temp=np.array(0)
count=snapshot_start
for i in range(len(snapshot_list)):
  rmax_temp=np.append(rmax_temp,importdata[count]['rmax']*1000)
  count=count+1
rmax=rmax_temp[1:len(rmax_temp)]
time=snapshot_list-snapshot_start #time starts from zero here where t=0 is at the snapshot where we start from
print("\nThese are the largest distances (in pc) from the center of mass for each snapshots we tracked:\n",rmax)


fig4=plt.figure()
ax1=fig4.add_subplot(1,1,1)
ax1.plot(time,rmax,color='b') 
ax1.set_xlabel('Time in Myr')
ax1.set_ylabel('largest distances (in parsec) from the center of mass')
ax1.set_title(cluster_groupid)
ax1.minorticks_on()
plt.tight_layout()
plotname="rmax_fromCM"+str(snapshot_start)+"_to_"+str(snapshot_end)+".png"
print("###################\nSaving plot of largest distances (in pc) from the center of mass in each Snapshot as:",plotname)
fig4.savefig(plot_path+plotname)
plt.close()
