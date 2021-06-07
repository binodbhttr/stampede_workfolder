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

snapnumber=656
snapshot_start=596
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

bins = np.arange(-25,25,0.1)


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
cbar_ax1 = fig1.add_axes([0.09, 0.185, 0.03, 0.65]) # position of gray colorbar (left, bottom, width, height)
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


x_tracked=tracked_data[cluster_group]["x_tracked"]
y_tracked=tracked_data[cluster_group]["y_tracked"]

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

keep_old=np.where((age <= .040) & (age >= 0.035) & (circle_radius<=region) & (abs(z) < 1.5))
#time=(snapnumber-snapshot_start)/1000
#keep_old=np.where((age <= time) & (circle_radius<=region) & (abs(z) < 1.5))
x0_old=x[keep_old]
y0_old=y[keep_old]
age_old=age[keep_old]
ax.scatter(x0_old,y0_old,c="red",s=7,label="Old Stars (<=40 & >=35) Myr")
#ax.scatter(x0,y0,c="blue",label="Young Stars (<3Myr)")

#age_old=1000*age_old #Converting to Myr
#cb2=ax.scatter(x0_old,y0_old,c=age_old,s=7,cmap=plt.cm.get_cmap('bwr'), vmin=np.min(age_old),vmax=np.max(age_old))
#cbar_ax = fig1.add_axes([0.81, 0.185, 0.03, 0.65]) # position of the colorbar (left, bottom, width, height)
#fig1.colorbar(cb2, cax=cbar_ax)
#cbar_ax.set_ylabel('Age(Myr)')
#cbar_ax.yaxis.label.set_size(12)

ax.scatter(x_tracked,y_tracked,c="green",s=100,marker="*",label="clusters 17",edgecolor='black')#cluster 17 is cornflowerblue
circle=plt.Circle((xcm,ycm),region,fill=False)
ax.add_patch(circle)

ax.set_xlim(np.min(x0)-5,np.max(x0)+5)
ax.set_ylim(np.min(y0)-5,np.max(y0)+5)

ax.legend(bbox_to_anchor=(1,0.5), loc='center left')
ax.set_title("time= "+str(snapnumber-snapshot_start)+" Myr")
fig1.savefig(str(snapnumber)+"test_sne.png",bbox_inches='tight',dpi=200)
