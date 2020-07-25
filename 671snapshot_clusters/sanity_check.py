from sl_utilities import distinct_colours as dc
from sl_utilities import distance_functions
import utilities as ut
from matplotlib import pyplot as plt
import numpy as np


cluster_groupid="snapshot671_cluster_group3" #Remember to change it if you  are changing the star culuster you are tracking in given snapshot
snapshot_start=671
snapshot_end=696
n=snapshot_end-snapshot_start+1
count=snapshot_start
importdata={}

#dict_keys(['age_tracked', 'avg_delta_rxyz', 'delta_rxyz', 'ind_tracked', 'mass_tracked', 'rmax', 'vx_tracked', 'vy_tracked', 'vz_tracked', 'x_tracked', 'xcm', 'xmax', 'xmin', 'y_tracked', 'ycm', 'ymax', 'ymin', 'z_tracked', 'zcm'])

colors=['cyan','blue','green','magenta','yellow','orange','purple','tan','lime','brown','grey','pink','navy','teal','wheat']


for i in range(n):
  file_name=cluster_groupid+"_tracked_data_snapshot"+str(count)
  path="./data/"
  importdata[count]=ut.io.file_hdf5(path+file_name) #reading data from each file and storing it to a dictionary importdata
  #print("\n Loaded data from the snapshot no.",count,"located in filename:",file_name,"to importdata[",count,"]\n#####\n")
  count=count+1
################################################

x1=importdata[671]['x_tracked']
y1=importdata[671]['y_tracked']
z1=importdata[671]['z_tracked']
mass1=importdata[671]['mass_tracked']
xcm1=distance_functions.cm(x1,mass1)
ycm1=distance_functions.cm(y1,mass1)
zcm1=distance_functions.cm(z1,mass1)
rcm1=distance_functions.dr(x1,y1,z1,mass1)

x2=importdata[672]['x_tracked']
y2=importdata[672]['y_tracked']
z2=importdata[672]['z_tracked']
mass2=importdata[672]['mass_tracked']

xcm2=distance_functions.cm(x2,mass2)
ycm2=distance_functions.cm(y2,mass2)
zcm2=distance_functions.cm(z2,mass2)
rcm2=distance_functions.dr(x2,y2,z2,mass2)

fig1=plt.figure(figsize=(8,10))
ax1=fig1.add_subplot(2,1,1)
ax1.scatter(x1,y1,c=colors)
ax1.scatter(xcm1,ycm1,c="red",marker="*",s=10)
ax1.set_xlabel("x")
ax1.set_ylabel("y")

ax2=fig1.add_subplot(2,1,2)
ax2.scatter(x2,y2,c=colors)
ax2.scatter(xcm2,ycm2,c="red",marker="*",s=10)





fig1.savefig("./plots/sanity_check.png")