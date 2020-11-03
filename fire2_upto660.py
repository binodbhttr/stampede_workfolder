##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
#This file is designed to plot animation of all 17 clusters using the file total_data_all_clusters_all_snapshots.pkl that contains tracked data of all clusters

import utilities as ut
from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from sl_utilities import distinct_colours as dc
from sl_utilities import distance_functions
import os
import pickle
import imageio
import matplotlib
import matplotlib.colors as colors
import gizmo_analysis as gizmo
import utilities as ut
from fof_analysis import fof
from matplotlib import rc #to use Latex math symbols like 'phi'
import astropy
from astropy.io import ascii



matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['xtick.labelsize'] = 16
matplotlib.rcParams['ytick.labelsize'] = 16
matplotlib.rc('text', usetex=False)
#############################################################################
#constants
#############################################################################

MsunToGm = 1.99e33
KpcToCm = 3.086e21
mp = 1.67e-24
#bin_edge = 10.
bin_edge = 30.

bins = np.arange(-25,25,0.1)

############################################################################
#read in sim files and find relevant particles
############################################################################
#STAMPEDE
simname = 'm12i_res7100_mhdcv'
simdir = '/scratch/projects/xsede/GalaxiesOnFIRE/mhdcv/m12i_res7100_mhdcv/1Myr/fire2/'

#reading in snapshot_times.txt file to get snapshot numbers and times
# columns are: snapshot scale-factor redshift time[Gyr] time_width[Myr]
snapshot_times = simdir + 'snapshot_times.txt'
snaptime_data = astropy.io.ascii.read(snapshot_times, guess=False, comment="#")
snaptime_data = np.genfromtxt(snapshot_times, usecols=(0,3), skip_header=4, dtype=float) #the first and fourth columns are the only ones we need 
snaps = np.array(snaptime_data[:,0]) #col1 = first column saved from text file #This is a collection of all snapshot nos.
times = np.array(snaptime_data[:,1]) #col4 = second column saved #This is a collection of times equivalent to those snapshot nos.
#######################################
#######################################
#######################################

datapath="./data_pkl/"  #this is teh path where the data of our tracked clusters is
snap=596                #this is the snapshot at which the clusters were taken from using the fof algorithm
simname = 'm12i_res7100_mhdcv' #this is the name of the simulation
simdir = '/scratch/projects/xsede/GalaxiesOnFIRE/mhdcv/m12i_res7100_mhdcv/1Myr/fire2'
#file_name="total_data_all_clusters_all_snapshots.pkl"  #this data contains the the tracked information of all clusters
snapshot_start=641  #snapshot to begin creating the figure
snapshot_end=660    #snapshot to stop at
#print("Testing if the loading of data was successful !! \n")
#print("x_tracked of snapshot 596 for cluster 1",importdata[596][1]["x_tracked"])
#print("x_cm of snapshot 596 for cluster 1",importdata[596][1]["xcm"])
#print("x_tracked of snapshot 597 for cluster 1",importdata[597][1]["x_tracked"])
#print("x_cm of snapshot 597 for cluster 1",importdata[597][1]["xcm"])
#print("x_tracked of snapshot 597 for cluster 2",importdata[597][2]["x_tracked"])
#print("x_cm of snapshot 597 for cluster 2",importdata[597][2]["xcm"])
################################
################################
#This section is just to create a list of clusters and display a print like this: These are the clusters groups we have tracked data of:
#['snapshot660_cluster_group1', 'snapshot660_cluster_group2', 'snapshot660_cluster_group3', 'snapshot660_cluster_group4', 'snapshot660_cluster_group5', 
#'snapshot660_cluster_group6', 'snapshot660_cluster_group7', 'snapshot660_cluster_group8', 'snapshot660_cluster_group9', 'snapshot660_cluster_group10']

cluster_groupid=[]
#total_clusters=len(importdata[snapshot_start])
total_clusters=17
for i in range(total_clusters):
    cluster_groupid.append("snapshot"+str(snapshot_start)+"_cluster_group"+str(i+1))

print("These are the clusters groups we have tracked data of:\n",cluster_groupid)

##############################
##############################
plot_path="./plots/plots_with_gas/" #creating a path to store the plots only if it does not exist
if not os.path.exists(plot_path):
  os.makedirs(plot_path)

n=snapshot_end-snapshot_start+1 # total no. of snapshots to generate image 

colors=['cyan','blue','green','magenta','yellow','teal','brown','darkslategray','lime','red','orange','purple','rosybrown','pink','navy','olive','cornflowerblue'] #17 colors
#colors=['cyan','blue','green','magenta','yellow','teal','brown','tan','lime','red'] #10 colors

#total_clusters=10 #just in case you want to change the no. of clusters. Setting 10 means select the first 10
plot_count=0       #counts the plot as we move forward with each snapshot
snapshot_list=np.arange(snapshot_start,snapshot_end+1) #create a list of snapshot numbers to plot to plot eg. [596,597, ...]
time=snapshot_list-snapshot_start #time starts from zero here where t=0 is at the snapshot where we start from eg. [0,1,2,3,....]
#time=[[0],[20],[40],[60],[80],[100]] #another way of providing fixed time stamps to look at


############################
############################
#Doing all this trouble to stor the initial mass of the clusters
fn="all_clusters_at_snapshot_"+str(snap)+".pkl"
with open(datapath+fn, "rb") as input:
  initial_cluster_data = pickle.load(input)
clustermass=[]

for c in range(len(initial_cluster_data)):
    m=np.sum(initial_cluster_data[c+1]["mass_tracked"]) #the cluster group name starts with 1 so we are adding extra 1 to the counter c
    clustermass.append(m)
############################
###############################






for i in range(len(time)):                              
    snapnumber=time[i]+snapshot_start      #this is the true snapshot number each time
    ###################
    file_name="all_clusters_at_snapshot_"+str(snapnumber)+".pkl"
    with open(datapath+file_name, "rb") as input:
      importdata = pickle.load(input)
    #loading data from each snapshot into importdata
    ##################
    ##################
    snaptime = times[np.where(snaps == snapnumber)][0] #time of snapshot in Gyr
    fig1=plt.figure()
    fig1.set_size_inches(7,7)
    ax=fig1.add_axes([0.17, 0.185, 0.65, 0.65]) #left, bottom, width, height
    cluster_count=0   #reset cluster_count to 0 as the j loop scans through all clusters to create a plot for given snapshot but i loop is for next snapshot
    for j in range(total_clusters):  #total_clusters tells us how many clusters to look for
        
        x=importdata[cluster_count+1]["x_tracked"]
        y=importdata[cluster_count+1]["y_tracked"]
        xcm=importdata[cluster_count+1]["xcm"]
        ycm=importdata[cluster_count+1]["ycm"]
        z=importdata[cluster_count+1]["z_tracked"]
        #cluster_mass=np.sum(importdata[cluster_count+1]["mass_tracked"])
        #cluster_mass_initial=np.sum(initial_cluster_data[cluster_count+1]["mass_tracked"])
        s1=ax.scatter(x,y,label=f'{clustermass[cluster_count]:.2e}'+' M$_{{\odot}}$',c=colors[cluster_count])
        #ax.scatter(np.abs(xcm),np.abs(ycm),c="black")
        ax.minorticks_on()
        ax.tick_params(labelsize=12)
        ax.set_xlabel("x (kpc)",fontsize=12,labelpad=3)
        ax.set_ylabel("y (kpc)",fontsize=12,labelpad=-5)
        #time_label= 'time = ' + f'{snaptime:.3f}' + ' Gyr'
        #ax.set_title("Clusters at T="+str(time[plot_count])+" in Myr")
        #ax.set_title("Clusters at "+time_label)
        cluster_count+=1
    
    #ax.legend(bbox_to_anchor=(1,0.5), loc='center left')
    handles,labels = ax.get_legend_handles_labels()
    sorted_handles= [x for _,x in sorted(zip(clustermass,handles),reverse=True)] #sort the handles (the colors next to the labels) based on clustermass
    sorted_legends= [x for _,x in sorted(zip(clustermass,labels),reverse=True)] #sort the labels based on the clustermass which is a list
    ax.legend(sorted_handles,sorted_legends,bbox_to_anchor=(1,0.5), loc='center left')
    #plt.tight_layout()
    part = gizmo.io.Read.read_snapshots(['all'],'snapshot_index', time[i]+snapshot_start, simulation_name=simname, simulation_directory=simdir, assign_hosts_rotation=True, assign_hosts=True)  
    t = np.max(part['star'].prop('form.time'))  
    
    rGas = part['gas'].prop('host.distance.principal.cylindrical')[:,0]
    zGas = part['gas'].prop('host.distance.principal.cylindrical')[:,1]
    
    xGas = part['gas'].prop('host.distance.principal.cartesian')[:,0]
    yGas = part['gas'].prop('host.distance.principal.cartesian')[:,1]
    zGas = part['gas'].prop('host.distance.principal.cartesian')[:,2]
    
    vxGas = part['gas'].prop('host.velocity.principal.cartesian')[:,0]
    vyGas = part['gas'].prop('host.velocity.principal.cartesian')[:,1]
    vzGas = part['gas'].prop('host.velocity.principal.cartesian')[:,2]
    
    mGas = part['gas']['mass']
    rhoGas = part['gas']['density']
    tGas = part['gas']['temperature']
    idGas = part['gas']['id']
    
    i_gas = np.where((rGas <= bin_edge) & (np.fabs(zGas) <= 1.5) & (part['gas']['density']*((MsunToGm/KpcToCm**3)/mp) >= 10.) & (tGas <= 1e4))
    
    x = xGas[i_gas]
    y = yGas[i_gas]
    z = zGas[i_gas]
    vx = vxGas[i_gas]
    vy = vyGas[i_gas]
    vz = vzGas[i_gas]
    m = mGas[i_gas]
    rho = part['gas'].prop('number.density')[i_gas]
    id = part['gas']['id'][i_gas]
    
    ###########################################################################
    #gas image (2d histogram)
    ###########################################################################
    #cold (< 10^4 K) gas in the midplane (|z| <= 1.5 kpc within bin_edge
    v =  np.where((rGas <= bin_edge) & (np.fabs(zGas) <= 1.5) & (tGas <= 1e4))
    face, xh, yh = np.histogram2d(part['gas'].prop('host.distance.principal.cartesian')[v,1][0],part['gas'].prop('host.distance.principal.cartesian')[v,0][0],bins=[bins,bins], weights=part['gas']['mass'][v])

    ###########################################################################

    norm = matplotlib.colors.LogNorm(vmin=1, vmax=1000) #the color range plotted
    im = ax.imshow(face/(((bins[1]-bins[0])*1000)**2),origin='lower',interpolation='nearest',norm=norm,extent=(-25,25,-25,25),cmap='binary') 
    
    #colorbar for the background gas density
    cmap_gray = matplotlib.cm.get_cmap('binary')
    norm1 = matplotlib.colors.LogNorm(vmin=1,vmax=1000)
    cbar_ax1 = fig1.add_axes([0.04, 0.185, 0.04, 0.65]) # position of gray colorbar (left, bottom, width, height)
    cb1 = fig1.colorbar(im, cax=cbar_ax1, ticklocation='left')
    cb1.set_label('$\Sigma$ (M$_{{\odot}}$/pc$^2$)', labelpad=-5, fontsize=12)
    
    #plot a scale bar 5 kpc long
    ax.plot([-22.5,-17.5], [22.5,22.5], 'k-', linewidth=5)
    label1 = "5 kpc"
    ax.text(-22.5, 20, label1, fontsize=12.5)
    
    #label the name of the galaxy on plot 
    label2="m12i mhdcv FIRE-2" 
    ax.text(-23.5,-23.5,label2,fontsize=13.5) 
    
    #label the time from the snapshot_times.txt file
    label3 = 'time = ' + f'{snaptime:.3f}' + ' Gyr'
    ax.text(3.5,-23.5,label3,fontsize=13.5)  
    
    fig1.savefig(plot_path+"snap"+str(plot_count+snapshot_start)+".png",bbox_inches='tight',dpi=200)
    plot_count+=1
    fig1.clf()
    plt.close()



     