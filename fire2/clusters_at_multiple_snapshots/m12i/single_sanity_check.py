import matplotlib as plt2 #need this for patches for shaded circles
import numpy as np
import matplotlib.pyplot as plt
import gizmo_analysis as gizmo 
import utilities as ut
import matplotlib.colors as colors
from matplotlib import rc #to use Latex math symbols like 'phi'
import astropy
from astropy.io import ascii
import matplotlib
import pdb
from importlib import reload
from sl_utilities import distance_functions
#pdb.set_trace()  #<--in case need to troubleshoot
import pickle
import os
from fof_analysis import fof


#########################################################################
#########################################################################
#load particle data here as you normally do to get the star information
#cut on Rxy, |z| and age
#run fof
#sind, sxcm, sycm, szcm, smtot, sgrpid, sr90, sr50, srmax =fof.find(x[si],y[si],z[si], b=b_kpc, mass=mass[si], ncut=ncut_min)
#srcm = np.sqrt(sxcm**2. + sycm**2.)

simname = 'm12i_res7100_mhdcv'
simdir = '/scratch/projects/xsede/GalaxiesOnFIRE/cr_suite/m12i_res7100/mhdcv/1Myr/fire2/'
snapnumber = 586

part = gizmo.io.Read.read_snapshots(['all'],'snapshot_index', snapnumber, simulation_name=simname, simulation_directory=simdir, assign_hosts_rotation=True, assign_hosts=True)  

rxyz     = part['star'].prop('host.distance.total')
Rxy      = part['star'].prop('host.distance.principal.cylindrical')[:,0]
x        = part['star'].prop('host.distance.principal')[:,0]
y        = part['star'].prop('host.distance.principal')[:,1]
z        = part['star'].prop('host.distance.principal')[:,2] 
mass     = part['star'].prop('mass')
feh      = part['star'].prop('metallicity.fe')
mgh      = part['star'].prop('metallicity.mg')
ch      = part['star'].prop('metallicity.c')
nh      = part['star'].prop('metallicity.n')
oh      = part['star'].prop('metallicity.o')
neh      = part['star'].prop('metallicity.ne')
sih      = part['star'].prop('metallicity.si')
sh      = part['star'].prop('metallicity.s')
cah      = part['star'].prop('metallicity.ca')
mgfe    = part['star'].prop('metallicity.mg - metallicity.fe')
ofe    = part['star'].prop('metallicity.o - metallicity.fe')
ids      = part['star'].prop('id')
id_child = part['star'].prop('id.child')
age      = part['star'].prop('age')
#select young stars within the disk for star clusters
keep = np.where((age <= .003) & ((Rxy < 20) & (Rxy>2)) & (abs(z) < 1.5))
#to run cluster finding on
rxyz0     = rxyz[keep]
Rxy0      = Rxy[keep]
x0        = x[keep]
y0        = y[keep]
z0        = z[keep]
mass0     = mass[keep]
feh0      = feh[keep]
mgh0=mgh[keep]
ch0=ch[keep]
nh0=nh[keep]
oh0=oh[keep]
neh0=neh[keep]
sih0=sih[keep]
sh0=sh[keep]
cah0=cah[keep]
mgfe0=mgfe[keep]
ofe0=ofe[keep]
id0       = ids[keep]
id_child0 = id_child[keep]
age0      = age[keep]

linking_length = 0.01 #10 parsec (unit here is in kpc)
ncut           = 5 #5 star particles

ind, xcm, ycm, zcm, mtot, grpid, r90, r50, rmax =fof.find(x0,y0,z0, b=linking_length, mass=mass0, ncut=ncut)
ngroup = len(mtot)



grp_index=1


ids_in_cluster = id0[ind[grp_index]]  #these are the star particle ids in each cluster
id_children_in_cluster = id_child0[ind[grp_index]]
age=age0[ind[grp_index]]
nstar = len(ids_in_cluster)
groupid=grpid[grp_index]
print('------------------------------------------------------------------------------------------------------------------')
print('grpid, nstar, xcm (kpc), ycm (kpc), zcm (kpc), mtot (msun), rmax (pc)')
print('%s     %i     %.4f     %.4f    %.4f     %.2e     %.1f ' % (grpid[grp_index], nstar, xcm[grp_index], ycm[grp_index], zcm[grp_index], mtot[grp_index], 1000*rmax[grp_index]))
print("age os stars in the cluster",age)
print('ids')
string = '[' 
for i in ids_in_cluster:
    string = string + str(i) + ', '

#get rid of last extra ,
length = len(string)-2
string = string[0:length] + ']'
print(string)

print('id children')
string = '[' 
for i in id_children_in_cluster:
    string = string + str(i) + ', '

#get rid of last extra ,
length = len(string)-2
string = string[0:length] + ']'
print(string)
print("These are the ids printed",ids_in_cluster)


cluster={"cluster_groupid":groupid,"no_of_star":nstar,"id":ids_in_cluster,"id_children":id_children_in_cluster,"xcm":xcm[grp_index],"ycm":ycm[grp_index],"zcm":zcm[grp_index],"mtot":mtot[grp_index],"r90":r90[grp_index],"r50":r50[grp_index],"rmax":rmax[grp_index],"x":x0[ind[grp_index]],"y":y0[ind[grp_index]],"z":z0[ind[grp_index]],"age":age0[ind[grp_index]],"feh":feh0[ind[grp_index]],"mgh":mgh0[ind[grp_index]],"ch":ch0[ind[grp_index]],"nh":nh0[ind[grp_index]],"oh":oh0[ind[grp_index]],"neh":neh0[ind[grp_index]],"sih":sih0[ind[grp_index]],"sh":sh0[ind[grp_index]],"cah":cah0[ind[grp_index]],"mgfe":mgfe0[ind[grp_index]],"ofe":ofe0[ind[grp_index]]}