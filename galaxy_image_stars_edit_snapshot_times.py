import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors
import numpy as np
import gizmo_analysis as gizmo
import utilities as ut
from fof_analysis import fof
from matplotlib import rc #to use Latex math symbols like 'phi'
import astropy
from astropy.io import ascii

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
snapnumber = 596 #the particular snapshot we're looking at right now

snaptime_data = astropy.io.ascii.read(snapshot_times, guess=False, comment="#")
snaptime_data = np.genfromtxt(snapshot_times, usecols=(0,3), skip_header=4, dtype=float) #the first and fourth columns are the only ones we need 
snaps = np.array(snaptime_data[:,0]) #col1 = first column saved from text file
times = np.array(snaptime_data[:,1]) #col4 = second column saved

snaptime = times[np.where(snaps == snapnumber)][0] #time of snapshot in Gyr
#print(f'{snaptime:.3f}') #this format saves 3 places after the decimal

#for troubleshooting on local machine
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
ind, xsp, ysp, zsp, msp, grpid, r90, r50, rmax =fof.find(sx[si],sy[si],sz[si], b=linking_length, mass=smass[si], ncut=ncut)
ngroup = len(msp)


ind=ind[0:10]
xsp=xsp[0:10]
ysp=ysp[0:10]
zsp=zsp[0:10]
msp=msp[0:10]
grpid=grpid[0:10]
rmax=rmax[0:10]
r90=r90[0:10]
r50=r50[0:10]


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

#background gray image of cold gas surface density
norm = matplotlib.colors.LogNorm(vmin=1, vmax=1000) #the color range plotted
im = ax.imshow(face/(((bins[1]-bins[0])*1000)**2),origin='lower',interpolation='nearest',norm=norm,extent=(-25,25,-25,25),cmap='binary') 

#overplotted star clusters scaled by their rmax (larger rmax = larger star)
rmax_parsec = np.array(rmax)*1000.
norm2 = matplotlib.colors.LogNorm(vmin=1e4,vmax=1e6)
im2 = ax.scatter(xsp,ysp,c=msp,s=rmax_parsec*2,alpha=0.8,cmap='plasma',norm=norm2,marker='*')
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)

#plot a scale bar 5 kpc long
ax.plot([-22.5,-17.5], [22.5,22.5], 'k-', linewidth=5)
label = '5 kpc'
ax.text(-22.5, 20, label, fontsize=12.5)

#label the name of the galaxy on plot 
label2='m12i mhdcv' 
ax.text(-23.5,-23.5,label2, fontsize=13.5) 

#label the time from the snapshot_times.txt file
label3 = 'time = ' + f'{snaptime:.3f}' + ' Gyr'
ax.text(3.5,-23.5,label3, fontsize=13.5)  

#colorbar for the background gas density
cb1 = fig.colorbar(im, cax=cbar_ax1, ticklocation='left')
cb1.set_label('$\Sigma$ (M$_{{\odot}}$/pc$^2$)', labelpad=-5, fontsize=14)

#colorbar for the star clusters
cb2 = fig.colorbar(im2, cax=cbar_ax2, ticklocation='right')
cb2.set_label('$M_{{cluster}}$ (M$_{{\odot}}$)', labelpad=5, fontsize=14)
cb2.ax.yaxis.set_label_position('right')

###########################################################################
#saving plot to pdf
###########################################################################

#odir = '/Users/sloebman/Desktop/'
odir = '/home1/07428/binod/stampede_workfolder/plots/plots_with_gas/'
ofile = odir + 'm12i_mhdcv_faceon_coldgas_starclusters_596.png' 
plt.savefig(ofile, bbox_inches='tight')
plt.close(fig)
plt.clf()



############################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
############################################################################
############################################################################
############################################################################
#Now we are going to track the clusters in snapshot 596 and see where they are at snapshot 686
#reading in snapshot_times.txt file to get snapshot numbers and times
# columns are: snapshot scale-factor redshift time[Gyr] time_width[Myr]
snapshot_times = simdir + 'snapshot_times.txt'
snapnumber_later = 685 #the particular snapshot we're looking at right now

snaptime_data = astropy.io.ascii.read(snapshot_times, guess=False, comment="#")
snaptime_data = np.genfromtxt(snapshot_times, usecols=(0,3), skip_header=4, dtype=float) #the first and fourth columns are the only ones we need 
snaps = np.array(snaptime_data[:,0]) #col1 = first column saved from text file
times = np.array(snaptime_data[:,1]) #col4 = second column saved

snaptime = times[np.where(snaps == snapnumber_later)][0] #time of snapshot in Gyr
#print(f'{snaptime:.3f}') #this format saves 3 places after the decimal

#for troubleshooting on local machine
#LOCAL
#simname='m12i.res57000'
#simdir='/Users/sloebman/Dropbox/RESEARCH/ANDREW/SIMS/m12i.res57000/'
#snapnumber = 600

part = gizmo.io.Read.read_snapshots(['all'],'snapshot_index', snapnumber_later, simulation_name=simname, simulation_directory=simdir, assign_hosts_rotation=True, assign_hosts=True)  
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
#tracking clusters/associations
###########################################################################
datapath="/home1/07428/binod/stampede_workfolder/data_pkl/" 
file_name="total_data_all_clusters_all_snapshots.pkl"

with open(datapath+file_name, "rb") as input:
    importdata = pickle.load(input)

print("Testing if the loading of data was successful !! \n")
print("x_cm of snapshot 596 for cluster 1 from the data I have is",importdata[596][1]["xcm"])
print("xcm of the snapshot 596 cluster 1 from the fof algorithm in real time is ",xsp[0])
print("x_cm of snapshot 597 for cluster 2",importdata[597][2]["xcm"])



ind_685=[]
xsp_685=[]
ysp_685=[]
zsp_685=[]
msp_685=[]
group=1
for i in range(10):
  ind_685.append(importdata[snapnumber_later][group]["ind_tracked"])
  xsp_685.append(importdata[snapnumber_later][group]["xcm"])
  ysp_685.append(importdata[snapnumber_later][group]["ycm"])
  zsp_685.append(importdata[snapnumber_later][group]["zcm"])
  msp_685.append(np.sum(importdata[snapnumber_later][group]["mass_tracked"])



###########################################################################
#gas image (2d histogram)
###########################################################################
#cold (< 10^4 K) gas in the midplane (|z| <= 1.5 kpc within bin_edge


###########################################################################
#plotting
###########################################################################
#plt.switch_backend('Qt5Agg')
fig2 = plt.figure()
fig2.set_size_inches(7,7)

ax = fig2.add_axes([0.17, 0.185, 0.65, 0.65]) #left, bottom, width, height
cbar_ax1 = fig2.add_axes([0.11, 0.185, 0.04, 0.65]) # position of gray colorbar
cbar_ax2 = fig2.add_axes([0.84, 0.185, 0.04, 0.65]) # position of colored colorbar

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

#background gray image of cold gas surface density
norm = matplotlib.colors.LogNorm(vmin=1, vmax=1000) #the color range plotted
im = ax.imshow(face/(((bins[1]-bins[0])*1000)**2),origin='lower',interpolation='nearest',norm=norm,extent=(-25,25,-25,25),cmap='binary') 

#overplotted star clusters scaled by their rmax (larger rmax = larger star)
rmax_parsec = np.array(rmax_685)*1000.
norm2 = matplotlib.colors.LogNorm(vmin=1e4,vmax=1e6)
im2 = ax.scatter(xsp_685,ysp_685,c=msp_685,s=rmax_parsec*2,alpha=0.8,cmap='plasma',norm=norm2,marker='*')
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)

#plot a scale bar 5 kpc long
ax.plot([-22.5,-17.5], [22.5,22.5], 'k-', linewidth=5)
label = '5 kpc'
ax.text(-22.5, 20, label, fontsize=12.5)

#label the name of the galaxy on plot 
label2='m12i mhdcv' 
ax.text(-23.5,-23.5,label2, fontsize=13.5) 

#label the time from the snapshot_times.txt file
label3 = 'time = ' + f'{snaptime:.3f}' + ' Gyr'
ax.text(3.5,-23.5,label3, fontsize=13.5)  

#colorbar for the background gas density
cb1 = fig2.colorbar(im, cax=cbar_ax1, ticklocation='left')
cb1.set_label('$\Sigma$ (M$_{{\odot}}$/pc$^2$)', labelpad=-5, fontsize=14)

#colorbar for the star clusters
cb2 = fig2.colorbar(im2, cax=cbar_ax2, ticklocation='right')
cb2.set_label('$M_{{cluster}}$ (M$_{{\odot}}$)', labelpad=5, fontsize=14)
cb2.ax.yaxis.set_label_position('right')

###########################################################################
#saving plot to pdf
###########################################################################

#odir = '/Users/sloebman/Desktop/'
odir = '/home1/07428/binod/stampede_workfolder/plots/plots_with_gas/'
ofile = odir + 'm12i_mhdcv_faceon_coldgas_starclusters_685.png' 
plt.savefig(ofile, bbox_inches='tight')
plt.close(fig2)
plt.clf()




