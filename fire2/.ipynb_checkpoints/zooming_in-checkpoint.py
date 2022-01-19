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

snapnumber=650
snap=snapnumber
data_path="./fire2_data_pkl/" 
gas_datapath="./fire2_gas_data_pkl/"

gas_file_name=simtype+"_gas_data"+str(snapnumber)+".pkl"
cluster_data_name="all_clusters_at_snapshot_"+str(snapnumber)+".pkl" 


#############################################################################
#constants
#############################################################################

MsunToGm = 1.99e33
KpcToCm = 3.086e21
mp = 1.67e-24
#bin_edge = 10.
bin_edge = 30.

bins = np.arange(-5,5,0.1)
#bins = np.arange(-25,25,0.1)

with open(data_path+cluster_data_name, "rb") as input:
    tracked_data= pickle.load(input)
      

with open(gas_datapath+gas_file_name, "rb") as input:
      import_gasdata = pickle.load(input)

cluster_group=17


fig1=plt.figure()
fig1.set_size_inches(7,7)
ax=fig1.add_axes([0.17, 0.185, 0.65, 0.65]) #left, bottom, width, height


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
#plot a scale bar 5 kpc long
#ax.plot([-22.5,-17.5], [22.5,22.5], 'k-', linewidth=5)
label1 = "5 kpc"
#ax.text(-22.5, 20, label1, fontsize=12.5)

label2="m12i mhdcv "+simtype.upper() 
#ax.text(-23.5,-23.5,label2,fontsize=13.5) 

#label the time from the snapshot_times.txt file
label3 = 'time = ' + f'{snaptime:.3f}' + ' Gyr'
#ax.text(3.5,22,label3,fontsize=13.5) #display at the top right  
#ax.text(-23.5,-21,label3,fontsize=13.5) #display time on the bottom left above simtype


x=tracked_data[cluster_group]["x_tracked"]
y=tracked_data[cluster_group]["y_tracked"]

ax.scatter(x,y,c="cornflowerblue")


#################################################Looking for young stars in the region###################
part=gizmo.io.Read.read_snapshots(['star'],'snapshot_index', snap, simulation_name=simname, simulation_directory=simdir, assign_hosts=True, assign_hosts_rotation=True)               #snap is the snapshot number here that changes everytime the loop iterates. It starts with sanpshot_start

age=part['star'].prop('age')
x=part['star'].prop('host.distance.principal')[:,0] #x component of the position of all stars 
y=part['star'].prop('host.distance.principal')[:,1] #y component of the position of all stars
z=part['star'].prop('host.distance.principal')[:,2] #z component of the position of all stars
Rxy = part['star'].prop('host.distance.principal.cylindrical')[:,0]

xcm=tracked_data[cluster_group]["xcm"]
ycm=tracked_data[cluster_group]["ycm"]
radius=((xcm)**2+(ycm)**2)**(1/2)
print("The cluster CM is at radius: ",radius)

circle_radius=((x-xcm)**2+(y-ycm)**2)**(1/2)

#keep = np.where((age <= .003) & ((Rxy < 20) & (Rxy>(radius))) & (abs(z) < 1.5))
region=5
keep = np.where((age <= .003) & (circle_radius<=region) & (abs(z) < 1.5))
x0=x[keep]
y0=y[keep]

ax.scatter(x0,y0,c="red",s=5)

circle=plt.Circle((xcm,ycm),region,fill=False)
ax.add_patch(circle)

fig1.savefig("zoomed_"+str(snapnumber)+"test.png")
