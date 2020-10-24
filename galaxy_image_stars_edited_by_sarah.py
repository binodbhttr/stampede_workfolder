import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors
import numpy as np
import gizmo_analysis as gizmo
import utilities as ut
import os
import pickle
#from fof_analysis import fof
from matplotlib import rc #to use Latex math symbols like 'phi'

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
snapnumber = 596

#LOCAL
#simname='m12i.res57000'
#simdir='/Users/sloebman/Dropbox/RESEARCH/ANDREW/SIMS/m12i.res57000/'
#snapnumber = 600

part = gizmo.io.Read.read_snapshots(['all'],'snapshot_index', snapnumber, simulation_name=simname, simulation_directory=simdir, assign_hosts_rotation=True, assign_hosts=True)  
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
    
srxyz     = part['star'].prop('host.distance.total')
sRxy      = part['star'].prop('host.distance.principal.cylindrical')[:,0]
sx        = part['star'].prop('host.distance.principal')[:,0]
sy        = part['star'].prop('host.distance.principal')[:,1]
sz        = part['star'].prop('host.distance.principal')[:,2] 
smass     = part['star'].prop('mass')
sfeh      = part['star'].prop('metallicity.fe')
sids      = part['star'].prop('id')
sid_child = part['star'].prop('id.child')
sage      = part['star'].prop('age')

############################################################################
#identifying star clusters/associations
###########################################################################

#select young stars within the disk for star clusters
si = np.where((sage <= .003) & (srxyz < 20) & (abs(sz) < 1.5))

linking_length = 0.01 #10 parsec
ncut           = 10 #10 star particles

#Running fof




clusterpath="./data_pkl/" 
snap=596
simname = 'm12i_res7100_mhdcv'
simdir = '/scratch/projects/xsede/GalaxiesOnFIRE/mhdcv/m12i_res7100_mhdcv/1Myr/1Myr_fire2'
file_name="total_data_all_clusters_all_snapshots.pkl"
snapshot_start=596
snapshot_end=696
with open(clusterpath+file_name, "rb") as input:
    importdata = pickle.load(input)
xsp=importdata[596][1]["x_tracked"]
ysp=importdata[596][1]["y_tracked"]
zsp=importdata[596][1]["z_tracked"]
ind=importdata[596][1]["ind_tracked"]
rmax=importdata[596][1]["rmax"]
msp=np.sum(importdata[596][1]["mass_tracked"])
#ind, xsp, ysp, zsp, msp, grpid, r90, r50, rmax =fof.find(sx[si],sy[si],sz[si], b=linking_length, mass=smass[si], ncut=ncut)
#ngroup = len(msp)

###########################################################################
#gas image (2d histogram)
###########################################################################
#cold (< 10^4 K) gas in the midplane (|z| <= 1.5 kpc within bin_edge
i =  np.where((rGas <= bin_edge) & (np.fabs(zGas) <= 1.5) & (tGas <= 1e4))
    
face, xh, yh = np.histogram2d(part['gas'].prop('host.distance.principal.cartesian')[i,1][0],part['gas'].prop('host.distance.principal.cartesian')[i,0][0],bins=[bins,bins], weights=part['gas']['mass'][i])

###########################################################################
#plotting
###########################################################################
#plt.switch_backend('Qt5Agg')
fig = plt.figure()
fig.set_size_inches(7,7)

ax = fig.add_axes([0.17, 0.185, 0.65, 0.65]) #left, bottom, width, height
cbar_ax1 = fig.add_axes([0.11, 0.185, 0.04, 0.65]) # position of gray colorbar
cbar_ax2 = fig.add_axes([0.84, 0.185, 0.04, 0.65]) # position of colored colorbar

#params = {"font.family":"serif","mathtext.fontset":"stix"}
#matplotlib.rcParams.update(params)
#Latex (which works on stampede & my laptop)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['xtick.labelsize'] = 16
matplotlib.rcParams['ytick.labelsize'] = 16
matplotlib.rc('text', usetex=False)

ax.set_xlabel('x (kpc)', fontsize=17, labelpad=2)
ax.set_ylabel('y (kpc)', fontsize=17, labelpad=-7) #neg labelpad = closer to axis
ax.set_xlim(-25, 25) 
ax.set_ylim(-25, 25)

norm = matplotlib.colors.LogNorm(vmin=1, vmax=1000)
#background gray image of cold gas surface density
im = ax.imshow(face/(((bins[1]-bins[0])*1000)**2),origin='lower',interpolation='nearest',norm=norm,extent=(-25,25,-25,25),cmap='binary') #min & max color range specified here

#overplotted star clusters scaled by their rmax (larger rmax = larger star)
rmax_parsec = np.array(rmax)*1000.
norm2 = matplotlib.colors.LogNorm(vmin=1e4,vmax=1e6)
im2 = ax.scatter(xsp,ysp,c="red",s=rmax_parsec*2,alpha=0.8,cmap='plasma',norm=norm2,marker='*') #i changed from line 152 to 151
#im2 = ax.scatter(xsp,ysp,c=msp,s=rmax_parsec*2,alpha=0.8,cmap='plasma',norm=norm2,marker='*') #color coded interms of mass
#im2 = ax.scatter(xsp,ysp,c=msp,s=rmax_parsec*2,alpha=0.8,cmap='plasma',norm=matplotlib.colors.LogNorm(),vmin=1e4,vmax=1e6,marker='*') #min & max color range specified here
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)

label='m12i mhdcv'
ax.text(-24,-24,label) #labeling the name of the galaxy on top of plot 

#colorbar for the background gas density
cmap_gray = matplotlib.cm.get_cmap('binary')
norm1 = matplotlib.colors.LogNorm(vmin=1,vmax=1000)
cb1 = fig.colorbar(im, cax=cbar_ax1, ticklocation='left')
cb1.set_label('$\Sigma$ (M$_{{\odot}}$/pc$^2$)', labelpad=-5, fontsize=14)

#colorbar for the star clusters
cmap_plasma = matplotlib.cm.get_cmap('plasma')
norm2 = matplotlib.colors.LogNorm(vmin=1e4,vmax=1e6)
cb2 = fig.colorbar(im2, cax=cbar_ax2, ticklocation='right')
cb2.set_label('$M_{{cluster}}$ (M$_{{\odot}}$)', labelpad=5, fontsize=14)
cb2.ax.yaxis.set_label_position('right')

###########################################################################
#saving plot to pdf
###########################################################################

#odir = '/Users/sloebman/Desktop/'
odir = '/home1/07428/binod/stampede_workfolder/'
ofile = odir + 'm12i_mhdcv_faceon_coldgas_starclusters_596.pdf' 
plt.savefig(ofile, bbox_inches='tight')
plt.close(fig)
plt.clf()
