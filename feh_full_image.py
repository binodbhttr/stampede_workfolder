import numpy as np
import matplotlib.pyplot as plt
import gizmo_analysis as gizmo 
import utilities as ut
from sl_utilities import distinct_colours as dc
from fof_analysis import fof
import matplotlib.colors as colors
from sl_utilities import load_pickle_dict
from sl_utilities import read_pickles
from sl_utilities import cp_dict
from importlib import reload
from scipy import stats
from matplotlib import rc #to use Latex math symbols like 'phi'

#to get metallicity for clusters later
#m12i_cluster_dict, m12i_complex_dict, m12f_cluster_dict, m12f_complex_dict, m12m_cluster_dict, m12m_complex_dict = read_pickles.get_final_dict()

simname = 'm12m.res7100'
simdir = '/Users/sloebman/Dropbox/RESEARCH/ANDREW/SIMS/'

part = gizmo.io.Read.read_snapshots(['star'], 'redshift', 0, assign_host_principal_axes=True, simulation_name=simname, simulation_directory=simdir+simname)

species='star'

x = part[species].prop('host.distance.principal')[:,0]
y = part[species].prop('host.distance.principal')[:,1]
z = part[species].prop('host.distance.principal')[:,2]
vx = part[species].prop('host.velocity.principal')[:,0]
vy = part[species].prop('host.velocity.principal')[:,1]
vz = part[species].prop('host.velocity.principal')[:,2]

mass = part[species].prop('mass')
massi = part[species].prop('form.mass')
feh = part[species].prop('metallicity.fe')
feh = np.array(feh)
    
age = part[species].prop('age')
age = np.array(age)
r = np.sqrt(x**2 + y**2 + z**2)
young = 0.003

keep = np.where((age < .003) & (r < 20) & (abs(z) < 1.5))
m50  = np.where((age <= 0.050) & (r < 20) & (abs(z) < 1.5))
m1000 = np.where((age <= 1) & (r < 20) & (abs(z) < 1.5))

mass_m50 = mass[m50]
x_m50 = x[m50]
y_m50 = y[m50]
z_m50 = z[m50]
r_m50 = r[m50]
age_m50 = age[m50]
feh_m50 = feh[m50]

mass_m1000 = mass[m1000]
x_m1000 = x[m1000]
y_m1000 = y[m1000]
z_m1000 = z[m1000]
r_m1000 = r[m1000]
age_m1000 = age[m1000]
feh_m1000 = feh[m1000]

xkeep   = x[keep]
ykeep   = y[keep]
zkeep   = z[keep]
mkeep   = mass[keep]
mikeep  = massi[keep]
vxkeep  = vx[keep]
vykeep  = vy[keep]
vzkeep  = vz[keep]
rkeep   = r[keep]
agekeep = age[keep]
fehkeep = feh[keep]

ind1, xcm1, ycm1, zcm1, mtot1, grpid1, r901, r501, rmax1 = fof.find(xkeep,ykeep,zkeep,b=.01, mass=mkeep, ncut=5)
ind2, xcm2, ycm2, zcm2, mtot2, grpid2, r902, r502, rmax2 = fof.find(xkeep,ykeep,zkeep,b=.02, mass=mkeep, ncut=5)

dict1 = cp_dict.get_dict(ind1, xcm1, ycm1, zcm1, mtot1, grpid1, r901, r501, rmax1, xkeep, ykeep, zkeep, mkeep, vxkeep, vykeep, vzkeep, rkeep, agekeep, fehkeep, mikeep)
dict2 = cp_dict.get_dict(ind2, xcm2, ycm2, zcm2, mtot2, grpid2, r902, r502, rmax2, xkeep, ykeep, zkeep, mkeep, vxkeep, vykeep, vzkeep, rkeep, agekeep, fehkeep, mikeep)

#----------------------------------------------------------------------------------------------------
#do histogram work outside of plot
bin_edge = 15 #how far out in the disk
width = .2 #sam: 100-200 pc for pretty pics, 750 for analysis 
ksbin = np.arange(-1.*(bin_edge)+width/2.,(bin_edge)+width/2.,width) 

den, xh, yh = np.histogram2d(y_m1000, x_m1000, weights=mass_m1000, bins=[ksbin, ksbin]) #good stars & histogram line up now <--  x & y have to be inverted in histogram
feh2, xedges2, yedges2, bin2 = stats.binned_statistic_2d(x_m1000, y_m1000, feh_m1000, 'mean', bins=[ksbin,ksbin])
m2 = np.ma.masked_where(den < 30000,  feh2) #masks below cut

#density peaks of young (< 50 Myr) stars: shows spiral structure 
cbin_edge = 15 
cwidth = .4 
cksbin = np.arange(-1.*(cbin_edge)+cwidth/2.,(cbin_edge)+cwidth/2.,cwidth) 

cden, cxh, cyh = np.histogram2d(y_m50, x_m50, weights=mass_m50, bins=[cksbin, cksbin]) 
#----------------------------------------------------------------------------------------------------
simname = 'm12i.res7100.reionize-late'
simdir = '/Users/sloebman/Dropbox/RESEARCH/ANDREW/SIMS/'

parti = gizmo.io.Read.read_snapshots(['star'], 'redshift', 0, assign_host_principal_axes=True, simulation_name=simname, simulation_directory=simdir+simname)

species='star'

ix = parti[species].prop('host.distance.principal')[:,0]
iy = parti[species].prop('host.distance.principal')[:,1]
iz = parti[species].prop('host.distance.principal')[:,2]
ivx = parti[species].prop('host.velocity.principal')[:,0]
ivy = parti[species].prop('host.velocity.principal')[:,1]
ivz = parti[species].prop('host.velocity.principal')[:,2]

imass = parti[species].prop('mass')
imassi = parti[species].prop('form.mass')
ifeh = parti[species].prop('metallicity.fe')
ifeh = np.array(ifeh)
    
iage = parti[species].prop('age')
iage = np.array(iage)
ir = np.sqrt(ix**2 + iy**2 + iz**2)

ikeep = np.where((iage < .003) & (ir < 20) & (abs(iz) < 1.5))
im50  = np.where((iage <= 0.050) & (ir < 20) & (abs(iz) < 1.5))
im1000 = np.where((iage <= 1) & (ir < 20) & (abs(iz) < 1.5))

imass_m50 = imass[im50]
ix_m50 = ix[im50]
iy_m50 = iy[im50]
iz_m50 = iz[im50]
ir_m50 = ir[im50]
iage_m50 = iage[im50]
ifeh_m50 = ifeh[im50]

imass_m1000 = imass[im1000]
ix_m1000 = ix[im1000]
iy_m1000 = iy[im1000]
iz_m1000 = iz[im1000]
ir_m1000 = ir[im1000]
iage_m1000 = iage[im1000]
ifeh_m1000 = ifeh[im1000]

ixkeep   = ix[ikeep]
iykeep   = iy[ikeep]
izkeep   = iz[ikeep]
imkeep   = imass[ikeep]
imikeep  = imassi[ikeep]
ivxkeep  = ivx[ikeep]
ivykeep  = ivy[ikeep]
ivzkeep  = ivz[ikeep]
irkeep   = ir[ikeep]
iagekeep = iage[ikeep]
ifehkeep = ifeh[ikeep]

iind1, ixcm1, iycm1, izcm1, imtot1, igrpid1, ir901, ir501, irmax1 = fof.find(ixkeep,iykeep,izkeep,b=.01, mass=imkeep, ncut=5)
iind2, ixcm2, iycm2, izcm2, imtot2, igrpid2, ir902, ir502, irmax2 = fof.find(ixkeep,iykeep,izkeep,b=.02, mass=imkeep, ncut=5)

#idict1 = cp_dict.get_dict(iind1, ixcm1, iycm1, izcm1, imtot1, igrpid1, ir901, ir501, irmax1, ixkeep, iykeep, izkeep, imkeep, ivxkeep, ivykeep, ivzkeep, irkeep, iagekeep, ifehkeep, imikeep)
#idict2 = cp_dict.get_dict(iind2, ixcm2, iycm2, izcm2, imtot2, igrpid2, ir902, ir502, irmax2, ixkeep, iykeep, izkeep, imkeep, ivxkeep, ivykeep, ivzkeep, irkeep, iagekeep, ifehkeep, imikeep)

#----------------------------------------------------------------------------------------------------
#do histogram work outside of plot
ibin_edge = 15 #how far out in the disk
iwidth = .2 #sam: 100-200 pc for pretty pics, 750 for analysis 
iksbin = np.arange(-1.*(ibin_edge)+iwidth/2.,(ibin_edge)+iwidth/2.,iwidth) 

iden, ixh, iyh = np.histogram2d(iy_m1000, ix_m1000, weights=imass_m1000, bins=[iksbin, iksbin]) #good stars & histogram line up now <--  x & y have to be inverted in histogram
ifeh2, ixedges2, iyedges2, ibin2 = stats.binned_statistic_2d(ix_m1000, iy_m1000, ifeh_m1000, 'mean', bins=[iksbin,iksbin])
m2i = np.ma.masked_where(iden < 30000,  ifeh2) #masks below cut

#density peaks of young (< 50 Myr) stars: shows spiral structure 
icbin_edge = 15 
icwidth = .4 
icksbin = np.arange(-1.*(icbin_edge)+icwidth/2.,(icbin_edge)+icwidth/2.,icwidth) 

icden, icxh, icyh = np.histogram2d(iy_m50, ix_m50, weights=imass_m50, bins=[icksbin, icksbin]) 
#----------------------------------------------------------------------------------------------------
#PLOTTING
fig = plt.figure()
fig.set_size_inches(10,4.5)
fig.subplots_adjust(wspace=.25, hspace=.25)
rc('font', **{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams["axes.edgecolor"] = "0.15"
plt.rcParams["axes.linewidth"]  = 1.75
plt.rcParams["axes.grid"] = False

ax1 = plt.subplot2grid((1, 2), (0, 0))
ax2 = plt.subplot2grid((1, 2), (0, 1))

#ax1 = fig.add_axes([0.07, 0.1, 0.85, 0.85]) #left, bottom, width, height 

ax1.set_xlabel('x (kpc)', fontsize=14)
ax1.set_ylabel('y (kpc)', fontsize=14, labelpad=-5)
ax1.set_xlim(-15, 15) 
ax1.set_ylim(-15, 15)

im1 = ax1.imshow(m2, interpolation='nearest', origin='low', extent=[xedges2[0],xedges2[-1],yedges2[0],yedges2[-1]], cmap = plt.cm.get_cmap('RdBu_r'), aspect='auto', vmin=-.3, vmax=0.5)

#im = ax.imshow(den/(width**2), interpolation='nearest',extent=(-1.*bin_edge,bin_edge,-1.*bin_edge,bin_edge), origin='lower', vmin=1e3, vmax=5e7, norm=colors.LogNorm(), cmap='pink')
#cb.set_label(r'$\Sigma_{*}$ (M$_{\odot}$/pc$^2$)')

ax1.contour(cden,extent=[cxh[0],cxh[-1],cyh[0],cyh[-1]], linewidths=[1], colors='gray', levels = [1.1e5])

#colorbar same for image and complexes
mean_feh = fehkeep
xcm = xcm2
ycm = ycm2
rmax = np.array(rmax2)*1000.

im1 = ax1.scatter(xcm, ycm, c=mean_feh, cmap=plt.cm.get_cmap('RdBu_r'), vmin=-.3,vmax=0.5, marker='o', s=rmax*2, zorder=100, edgecolor='black', linewidth=1.2)

ax1.text(10.7,13.7,'m12m', color='white', fontsize=9.5)

fig.subplots_adjust(left=0.18, bottom=0.14)
#cbar_ax = fig.add_axes([0.84, 0.1, 0.05, 0.85]) #left, bottom, width, height
cbar_ax = fig.add_axes([0.074, 0.14, 0.022, 0.74])
cb = fig.colorbar(im1, cax=cbar_ax, ticklocation='left')
cb.set_label('[Fe/H] (dex)', labelpad=-5, fontsize=12)

#--------------------------------------------------------------------------------
ax2.set_xlabel('x (kpc)', fontsize=14)
ax2.set_ylabel('y (kpc)', fontsize=14, labelpad=-5)
ax2.set_xlim(-15, 15) 
ax2.set_ylim(-15, 15)

im2 = ax2.imshow(m2i, interpolation='nearest', origin='low', extent=[ixedges2[0],ixedges2[-1],iyedges2[0],iyedges2[-1]], cmap = plt.cm.get_cmap('RdBu_r'), aspect='auto', vmin=-.3, vmax=0.5)

ax2.contour(icden,extent=[icxh[0],icxh[-1],icyh[0],icyh[-1]], linewidths=[1], colors='gray', levels = [1.1e5])

#colorbar same for image and complexes
imean_feh = ifehkeep
ixcm2, iycm2
ixcm = ixcm2
iycm = iycm2
irmax = np.array(irmax2)*1000.

#imean_feh = idict2['mean_feh']
#ixcm = idict2['xcm']
#iycm = idict2['ycm']
#irmax = np.array(idict2['rmax'])*1000.

im2 = ax2.scatter(ixcm, iycm, c=imean_feh, cmap=plt.cm.get_cmap('RdBu_r'), vmin=-.3,vmax=0.5, marker='o', s=irmax*2, zorder=100, edgecolor='black', linewidth=1.2)

ax2.text(10.7,13.7,'m12i', color='white', fontsize=9.5)

outfile = '/Users/sloebman/Dropbox/RESEARCH/ANDREW/clusters/figures/grants/m12im_hist2d_mean_feh_1Gyr_contour_stars_50Myr_complexes_3Myr_xy.eps'
plt.savefig(outfile)

plt.close()
plt.clf()

#---------------------------------------------------------------------------

#PLOTTING
fig = plt.figure()
fig.set_size_inches(7,7)
#fig.subplots_adjust(wspace=0.25, hspace=0.25)

#ax3 = plt.subplot2grid((1, 1), (0, 0))
#ax3 = plt.subplot(111)
#ax3.get_position() #[left, bottom], [width, height]] #Bbox([[0.125, 0.11], [0.9, 0.88]])
#ax3.set_position([0.15, 0.15, 0.7, 0.7])
ax3 = fig.add_axes([0.31, 0.185, 0.65, 0.65]) #left, bottom, width, height

rc('font', **{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16

ax3.set_xlabel('x (kpc)', fontsize=17, labelpad=2)
ax3.set_ylabel('y (kpc)', fontsize=17, labelpad=-7) #neg labelpad = closer to axis
ax3.set_xlim(-15, 15) 
ax3.set_ylim(-15, 15)

im3 = ax3.imshow(m2, interpolation='nearest', origin='low', extent=[xedges2[0],xedges2[-1],yedges2[0],yedges2[-1]], cmap = plt.cm.get_cmap('RdBu_r'), aspect='auto', vmin=-.3, vmax=0.5)

ax3.contour(cden,extent=[cxh[0],cxh[-1],cyh[0],cyh[-1]], linewidths=[1], colors='gray', levels = [1.1e5])

#colorbar same for image and complexes
mean_feh = feh2
xcm = xcm2
ycm = ycm2
rmax = np.array(rmax2)*1000.

im3 = ax3.scatter(xcm, ycm, c=mean_feh, cmap=plt.cm.get_cmap('RdBu_r'), vmin=-.3,vmax=0.5, marker='o', s=rmax*2, zorder=100, edgecolor='black', linewidth=1.2)

ax3.text(11.25,13.5,'m12m', color='white', fontsize=12)

#fig.subplots_adjust(left=0.2, bottom=0.15)
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14


cbar_ax = fig.add_axes([0.13, 0.185, 0.04, 0.65]) #left, bottom, width, height
cb = fig.colorbar(im1, cax=cbar_ax, ticklocation='left')
cb.set_label('[Fe/H] (dex)', labelpad=-5, fontsize=14)


outfile = '/Users/sloebman/Dropbox/RESEARCH/ANDREW/clusters/figures/grants/m12m_hist2d_mean_feh_1Gyr_contour_stars_50Myr_complexes_3Myr_xy.eps'
plt.savefig(outfile)

plt.close()
plt.clf()





