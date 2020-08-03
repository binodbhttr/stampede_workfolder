from sl_utilities import distinct_colours as dc
from sl_utilities import distance_functions
import utilities as ut
from matplotlib import pyplot as plt
import numpy as np
import os
#this file is targetted to help avoid running tonns of snapshots each time to analyse the tracked data
# I am working towards loading each data from all files that are present in data folder and generating all the kinds of plots here
# It is assumed we know the snapshot_start and snapshot_end value
################################################
#Now loading data from each file and storing it to a dictionary importdata which can be accessed as importdata[snapshotnumber]['key']
#For example to load the x coordinates of the tracked stars in snapshot no. 690, importdata[690]['x_tracked']

cluster_groupid=["snapshot596_cluster_group01","snapshot596_cluster_group02","snapshot596_cluster_group03","snapshot596_cluster_group04","snapshot596_cluster_group05","snapshot596_cluster_group8","snapshot596_cluster_group15"] #Remember to change it if you  are changing the star culuster you are tracking in given snapshot

#cluster_groupid=["snapshot596_cluster_group01","snapshot596_cluster_group02"] #Remember to change it if you  are changing the star culuster you are tracking in given snapshot

plot_path="./plots/allculsters_single_figure/" #creating a path to store the plots only if it does not exist
if not os.path.exists(plot_path):
  os.makedirs(plot_path)
  
#dict_keys(['age_tracked', 'avg_delta_rxyz', 'delta_rxyz', 'ind_tracked', 'mass_tracked', 'rmax', 'vx_tracked', 'vy_tracked', 'vz_tracked', 'x_tracked', 'xcm', 'xmax', 'xmin', 'y_tracked', 'ycm', 'ymax', 'ymin', 'z_tracked', 'zcm']) #these are the keys of dictionary importdata[snapshotnumber]


snapshot_start=596
snapshot_end=696
n=snapshot_end-snapshot_start+1
count=snapshot_start
cluster_count=0
importdata={}
for i in range(len(cluster_groupid)):
    data_each_cluster={}
    for j in range(n):
        file_name=cluster_groupid[i]+"_tracked_data_snapshot"+str(count)
        data_path="./data/"+cluster_groupid[cluster_count]+"/" #this is the path of the directory where our data is present
        data_each_cluster[count]=ut.io.file_hdf5(data_path+file_name) #reading data from each file and storing it to a dictionary importdata
        count=count+1
    importdata.update({cluster_groupid[cluster_count]:data_each_cluster}) #storing data from each starcluster and storing it to a dictionary importdata with the cluster id as the key

    os.system('clear')
    #print("\n Loaded data from the snapshot no.",count,"located in filename:",file_name,"to importdata[",count,"]\n#####\n")
    count=snapshot_start
    cluster_count+=1

    
########################################################
########################################################
#Testing if the data was loaded sucessfully       
x=importdata["snapshot596_cluster_group01"][596]["x_tracked"]
a=importdata[cluster_groupid[1]][596]["x_tracked"]
print("Making sure if the data from different clusters are loaded\n",cluster_groupid[0],"\n",x,"\n\n",cluster_groupid[1],"\n",a) 
########################################################
########################################################



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
        avg_r_cm_temp=np.append(avg_r_cm_temp,importdata[cluster_groupid[cluster_count]][snapshot_count]["avg_delta_rxyz"])
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







#########################################################
#########################################################
#plotting the standard deviation of the max r from the center of Mass for each snapshot for all clusters
fig2=plt.figure(figsize=(8,5))
ax1=fig2.add_subplot(1,1,1)
ax1.minorticks_on()
ax1.set_xlabel('Time in Myr')
ax1.set_ylabel('stdev of Distance of stars from CM of the cluster in parsec')
snapshot_list=np.arange(snapshot_start,snapshot_end+1) #create a list of snapshots to plot
time=snapshot_list-snapshot_start #time starts from zero here where t=0 is at the snapshot where we start from
n=snapshot_end-snapshot_start+1
cluster_count=0
for i in range(len(cluster_groupid)):
    std_r_cm_temp=np.array(0)
    snapshot_count=snapshot_start
    for j in range(n):
        std_r_cm_temp=np.append(std_r_cm_temp,np.nanstd(importdata[cluster_groupid[cluster_count]][snapshot_count]['delta_rxyz']))
        snapshot_count+=1
    std_r_cm=std_r_cm_temp[1:len(std_r_cm_temp)]*1000 #converting to parsec as it was in kpc
    slope, intercept = np.polyfit(time,std_r_cm, 1)
    ax1.plot(time,std_r_cm,label=cluster_groupid[cluster_count]+",slope="+str(round(slope,4)))
    #ax1.plot(time,slope*time+intercept,label=cluster_groupid[cluster_count]+"fitted")
    ax1.legend(loc='upper left')
    cluster_count+=1
plt.tight_layout()
fig2.savefig(plot_path+"stdev_of_distance_of_stars_from_CM_all_clusters.png",dpi=150)
#########################################################
#########################################################






#########################################################
#########################################################
#plotting the standard deviation of the magnitude of velocitiy of the stars in each snapshot for all clusters

fig3=plt.figure(figsize=(8,5))
ax1=fig3.add_subplot(1,1,1)
ax1.minorticks_on()
ax1.set_xlabel('Time in Myr')
ax1.set_ylabel('stdev of Magnitude of Velocity of stars in the cluster in parsec/Myr')
snapshot_list=np.arange(snapshot_start,snapshot_end+1) #create a list of snapshots to plot
time=snapshot_list-snapshot_start #time starts from zero here where t=0 is at the snapshot where we start from
n=snapshot_end-snapshot_start+1
cluster_count=0
for i in range(len(cluster_groupid)):
    std_v_mag_temp=np.array(0)
    snapshot_count=snapshot_start
    for j in range(n):
        vx=importdata[cluster_groupid[cluster_count]][snapshot_count]['vx_tracked']
        vy=importdata[cluster_groupid[cluster_count]][snapshot_count]['vy_tracked']
        vz=importdata[cluster_groupid[cluster_count]][snapshot_count]['vz_tracked']
        v_mag=(vx**2+vy**2+vz**2)**(1/2)*1.023 #calculated magnitude of velocity and converted it into parsec/Myr
        std_v_mag_temp=np.append(std_v_mag_temp,np.nanstd(v_mag))
        snapshot_count+=1
    std_v_mag=std_v_mag_temp[1:len(std_v_mag_temp)]
    slope, intercept = np.polyfit(time,std_v_mag, 1)
    ax1.plot(time,std_v_mag,label=cluster_groupid[cluster_count]+",slope="+str(round(slope,4)))
    #ax1.plot(time,slope*time+intercept,label=cluster_groupid[cluster_count]+"fitted")
    ax1.legend(loc='upper right')
    cluster_count+=1
plt.tight_layout()    
fig3.savefig(plot_path+"stdev_of_Vmag_of_stars_all_clusters.png",dpi=150)
#########################################################
#########################################################





#########################################################
#########################################################
fig4=plt.figure(figsize=(8,7))
ax1=fig4.add_subplot(1,1,1)
ax1.minorticks_on()
ax1.set_xlabel('Time in Myr')
ax1.set_ylabel('stdev of Magnitude of Velocity of stars in the cluster in parsec/Myr')
snapshot_list=np.arange(snapshot_start,snapshot_end+1) #create a list of snapshots to plot
time=snapshot_list-snapshot_start #time starts from zero here where t=0 is at the snapshot where we start from
n=snapshot_end-snapshot_start+1
cluster_count=0
for i in range(len(cluster_groupid)):
    std_v_mag_temp=np.array(0)
    snapshot_count=snapshot_start
    for j in range(n):
        vx=importdata[cluster_groupid[cluster_count]][snapshot_count]['vx_tracked']
        vy=importdata[cluster_groupid[cluster_count]][snapshot_count]['vy_tracked']
        vz=importdata[cluster_groupid[cluster_count]][snapshot_count]['vz_tracked']
        v_mag=(vx**2+vy**2+vz**2)**(1/2)*1.023 #calculated magnitude of velocity and converted it into parsec/Myr
        std_v_mag_temp=np.append(std_v_mag_temp,np.nanstd(v_mag))
        snapshot_count+=1
    std_v_mag=std_v_mag_temp[1:len(std_v_mag_temp)]
    slope, intercept = np.polyfit(time,std_v_mag, 1)
    #ax1.plot(time,std_v_mag,label=cluster_groupid[cluster_count]+",slope="+str(round(slope,4)))
    ax1.plot(time,slope*time+intercept,label=cluster_groupid[cluster_count]+"_fitted, "+"slope = "+str(round(slope,4)))
    #ax1.legend(loc='upper right')
    ax1.legend(bbox_to_anchor=(0,-0.15), loc='upper left')
    cluster_count+=1
plt.tight_layout()    
fig4.savefig(plot_path+"linear_fit_stdev_of_Vmag_of_stars_all_clusters.png",dpi=150)
#########################################################
#########################################################





#########################################################
#########################################################
#plotting the largest distance from the Center of Mass of the stars in each snapshot for all clusters

fig5=plt.figure(figsize=(8,5))
ax1=fig5.add_subplot(1,1,1)
ax1.minorticks_on()
ax1.set_xlabel('Time in Myr')
ax1.set_ylabel('largest distances (in parsec) from the center of mass')
snapshot_list=np.arange(snapshot_start,snapshot_end+1) #create a list of snapshots to plot
time=snapshot_list-snapshot_start #time starts from zero here where t=0 is at the snapshot where we start from
n=snapshot_end-snapshot_start+1
cluster_count=0
for i in range(len(cluster_groupid)):
    rmax_temp=np.array(0)
    snapshot_count=snapshot_start
    for i in range(n):
        rmax_temp=np.append(rmax_temp,importdata[cluster_groupid[cluster_count]][snapshot_count]['rmax'])
        snapshot_count+=1
    rmax=rmax_temp[1:len(rmax_temp)]*1000 #converted rmax to parsec
    slope, intercept = np.polyfit(time,rmax, 1)
    ax1.plot(time,rmax,label=cluster_groupid[cluster_count]+",slope="+str(round(slope,4)))
    #ax1.plot(time,slope*time+intercept,label=cluster_groupid[cluster_count]+"fitted")
    ax1.legend(loc='upper left')
    cluster_count+=1
plt.tight_layout()
fig5.savefig(plot_path+"largest_distance_from_CM_all_clusters.png",dpi=150)
#########################################################
#########################################################




#########################################################
#########################################################
#plotting the distance of CM from the galactic center in each snapshot for all clusters (in kpc)

fig6=plt.figure(figsize=(8,7))
ax1=fig6.add_subplot(1,1,1)
ax1.minorticks_on()
ax1.set_xlabel('Time in Myr')
ax1.set_ylabel('Distance of Center of Mass from the galactic center in kpc')
snapshot_list=np.arange(snapshot_start,snapshot_end+1) #create a list of snapshots to plot
time=snapshot_list-snapshot_start #time starts from zero here where t=0 is at the snapshot where we start from
n=snapshot_end-snapshot_start+1
cluster_count=0
for i in range(len(cluster_groupid)):
    dist_cm_temp=np.array(0)
    snapshot_count=snapshot_start
    for j in range(n):
        xcm=importdata[cluster_groupid[cluster_count]][snapshot_count]['xcm']
        ycm=importdata[cluster_groupid[cluster_count]][snapshot_count]['ycm']
        dis=(xcm**2+ycm**2)**(1/2) #left the calulation in kpc
        dist_cm_temp=np.append(dist_cm_temp,dis)
        snapshot_count+=1
    dis_cm=dist_cm_temp[1:len(dist_cm_temp)]
    slope, intercept = np.polyfit(time,dis_cm, 1)
    ax1.plot(time,dis_cm,label=cluster_groupid[cluster_count]+",slope="+str(round(slope,4)))
    #ax1.plot(time,slope*time+intercept,label=cluster_groupid[cluster_count]+"fitted")
    #ax1.legend(loc='upper right')
    ax1.legend(bbox_to_anchor=(0,-0.15), loc='upper left')
    cluster_count+=1
plt.tight_layout()    
fig6.savefig(plot_path+"distance_of_CM_from_galactic_center_all_clusters.png",dpi=150)
#########################################################
#########################################################



