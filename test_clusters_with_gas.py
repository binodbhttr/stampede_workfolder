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

datapath="./data_pkl/"
snap=690
simname = 'm12i_res7100_mhdcv'
simdir = '/scratch/projects/xsede/GalaxiesOnFIRE/mhdcv/m12i_res7100_mhdcv/1Myr/fire2'
file_name="total_data_all_clusters_all_snapshots.pkl"
snapshot_start=660
snapshot_end=660


with open(datapath+file_name, "rb") as input:
    importdata = pickle.load(input)

print("Testing if the loading of data was successful !! \n")
print("x_tracked of snapshot 596 for cluster 1",importdata[596][1]["x_tracked"])
print("x_cm of snapshot 596 for cluster 1",importdata[596][1]["xcm"])
print("x_tracked of snapshot 597 for cluster 1",importdata[597][1]["x_tracked"])
print("x_cm of snapshot 597 for cluster 1",importdata[597][1]["xcm"])
print("x_tracked of snapshot 597 for cluster 2",importdata[597][2]["x_tracked"])
print("x_cm of snapshot 597 for cluster 2",importdata[597][2]["xcm"])

cluster_groupid=[]
total_clusters=len(importdata[snapshot_start])
for i in range(total_clusters):
    cluster_groupid.append("snapshot"+str(snapshot_start)+"_cluster_group"+str(i+1))

print("These are the clusters groups we have tracked data of:\n",cluster_groupid)


plot_path="./plots/plots_with_gas/" #creating a path to store the plots only if it does not exist
if not os.path.exists(plot_path):
  os.makedirs(plot_path)


n=snapshot_end-snapshot_start+1

colors=['cyan','blue','green','magenta','yellow','teal','brown','tan','lime','red']

total_clusters=10
plot_count=0
snapshot_list=np.arange(snapshot_start,snapshot_end+1) #create a list of snapshots to plot
time=snapshot_list-snapshot_start #time starts from zero here where t=0 is at the snapshot where we start from
#time=[[0],[20],[40],[60],[80],[100]] #another way of providing fixed time stamps to look at
for i in range(len(time)):
    snapnumber=time[i]+snapshot_start
    snaptime = times[np.where(snaps == snapnumber)][0] #time of snapshot in Gyr
    fig1=plt.figure()
    fig1.set_size_inches(7,7)
    ax=fig1.add_axes([0.17, 0.185, 0.65, 0.65]) #left, bottom, width, height
    cluster_count=0
    for j in range(total_clusters):
        
        x=importdata[time[i]+snapshot_start][cluster_count+1]["x_tracked"]
        y=importdata[time[i]+snapshot_start][cluster_count+1]["y_tracked"]
        xcm=importdata[time[i]+snapshot_start][cluster_count+1]["xcm"]
        ycm=importdata[time[i]+snapshot_start][cluster_count+1]["ycm"]
        z=importdata[time[i]+snapshot_start][cluster_count+1]["z_tracked"]
        
        s1=ax.scatter(x,y,label="cluster_group"+str(cluster_count+1),c=colors[cluster_count])
        #ax.scatter(np.abs(xcm),np.abs(ycm),c="black")
        ax.minorticks_on()
        ax.set_xlabel("x position in kpc")
        ax.set_ylabel("y position in kpc")
        time_label= 'time = ' + f'{snaptime:.3f}' + ' Gyr'
        #ax.set_title("Clusters at T="+str(time[plot_count])+" in Myr")
        ax.set_title("Clusters at "+time_label)
        cluster_count+=1
    
    ax.legend(bbox_to_anchor=(1,0.5), loc='center left')
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
    cb1.set_label('$\Sigma$ (M$_{{\odot}}$/pc$^2$)', labelpad=-5, fontsize=14)
    
    #plot a scale bar 5 kpc long
    ax.plot([-22.5,-17.5], [22.5,22.5], 'k-', linewidth=5)
    label1 = "5 kpc"
    ax.text(-22.5, 20, label1, fontsize=12.5)
    
    #label the name of the galaxy on plot 
    label2="m12i mhdcv" 
    ax.text(-23.5,-23.5,label2,fontsize=13.5) 
    
    #label the time from the snapshot_times.txt file
    #label3 = 'time = ' + f'{snaptime:.3f}' + ' Gyr'
    #ax.text(3.5,-23.5,label3,fontsize=13.5)  
    
    fig1.savefig(plot_path+"snap"+str(plot_count+snapshot_start)+".png",bbox_inches='tight',dpi=200)
    plot_count+=1
    fig1.clf()
    plt.close()






images=[]
p=[600,630,660,690]
for img in p:
    images.append(imageio.imread(plot_path+"snap"+str(img) + '.png'))
    #os.remove(plot_path+"snap"+str(k) + ".png")
imageio.mimsave(plot_path+"animation.gif", images, duration = 1/8)
     