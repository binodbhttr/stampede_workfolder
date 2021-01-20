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

#this is just for m12i!  need to change name for m12m cloud catalog and input data for fof

dirname_gmc = '/home1/07428/binod/stampede_workfolder/sarah_paper_plots/m12m/' #change this! 
filename_gmc = 'cloud_props_m12m_mhd_stamp_fire2650.txt'


#load gmc catalog data
data_gmc  = astropy.io.ascii.read(dirname_gmc+filename_gmc) #data_gmc.keys() #to see content
#gind = np.where(data_gmc['mass'] > 1e6)  #leave this commented for now

gxcm  = data_gmc['xcm']
gycm  = data_gmc['ycm']
gzcm  = data_gmc['zcm']
gmtot = data_gmc['mass']
gr90  = data_gmc['r_90']
gn = len(gxcm)


#########################################################################
#########################################################################
#load particle data here as you normally do to get the star information
#cut on Rxy, |z| and age
#run fof
#sind, sxcm, sycm, szcm, smtot, sgrpid, sr90, sr50, srmax =fof.find(x[si],y[si],z[si], b=b_kpc, mass=mass[si], ncut=ncut_min)
#srcm = np.sqrt(sxcm**2. + sycm**2.)

simname = 'm12m_res7100_mhdcv'
simdir = '/scratch/projects/xsede/GalaxiesOnFIRE/cr_suite/m12m_res7100/mhdcv/1Myr/'
snapnumber = 650
part = gizmo.io.Read.read_snapshots(['all'],'snapshot_index', snapnumber, simulation_name=simname, simulation_directory=simdir, assign_hosts_rotation=True, assign_hosts=True)  

rxyz     = part['star'].prop('host.distance.total')
Rxy      = part['star'].prop('host.distance.principal.cylindrical')[:,0]
x        = part['star'].prop('host.distance.principal')[:,0]
y        = part['star'].prop('host.distance.principal')[:,1]
z        = part['star'].prop('host.distance.principal')[:,2] 
mass     = part['star'].prop('mass')
feh      = part['star'].prop('metallicity.fe')
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
id0       = ids[keep]
id_child0 = id_child[keep]
age0      = age[keep]

linking_length = 0.008 #4 parsec
ncut           = 4 #4 star particles

ind, xcm, ycm, zcm, mtot, grpid, r90, r50, rmax =fof.find(x0,y0,z0, b=linking_length, mass=mass0, ncut=ncut)
ngroup = len(mtot)

bool_arr = distance_functions.array_embedded_check(xcm, ycm, zcm, rmax, gxcm, gycm, gzcm, gr90)  #this is just a little function that i wrote (but didn't check) to see if the cluster is inside the radius of gmc.  it will need to be visualized to confirm this

export_cluster={}
for grp_index in range(ngroup):  #iterate over each group
    cluster={}
    ids_in_cluster = id0[ind[grp_index]]  #these are the star particle ids in each cluster
    id_children_in_cluster = id_child0[ind[grp_index]]
    nstar = len(ids_in_cluster)
    groupid=grpid[grp_index]
    print('------------------------------------------------------------------------------------------------------------------')
    print('grpid, nstar, xcm (kpc), ycm (kpc), zcm (kpc), mtot (msun), rmax (pc)')
    print('%s     %i     %.4f     %.4f    %.4f     %.2e     %.1f ' % (grpid[grp_index], nstar, xcm[grp_index], ycm[grp_index], zcm[grp_index], mtot[grp_index], 1000*rmax[grp_index]))
   
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
    #feh_in_cluster=feh[ids_in_cluster]
    cluster={"cluster_groupid":groupid,"no_of_star":nstar,"id":ids_in_cluster,"id_children":id_children_in_cluster,"is_embedded":bool_arr[grp_index],"xcm":xcm[grp_index],"ycm":ycm[grp_index],"zcm":zcm[grp_index],"mtot":mtot[grp_index],"r90":r90[grp_index],"r50":r50[grp_index],"rmax":rmax[grp_index]}
    export_cluster.update({groupid:cluster})



######################################################################
#######################################################################

print('------------------------------------------------------------------------------------------------------------------')

path="./fire2_data_pkl/" #creating a path to store the data only if it does not exist
if not os.path.exists(path):
  os.makedirs(path)

file_name="fire2_associations_"+simname+"_snapshot_"+str(snapnumber) 

with open(path+file_name, 'wb') as output:
    # Pickle dictionary using protocol 0.
    pickle.dump(export_cluster, output)
    
########################
#to test if the information was stored properly
with open(path+file_name, "rb") as fp:
    import_cluster = pickle.load(fp)
        
print("\n\n\nThe ids of the cluster with group id 1 is",import_cluster[1]["id"])
########################
print("Total associations present is",len(import_cluster))


###################################
####################################
#####################################
#######################################
#########################################

#modify this loop as needed
'''
sn = len(xcm)
for i in range(sn):
    if (srcm[i] >= 6) and (srcm[i] <= 7) and (bool_arr[i] == True):
        #print(srcm[i], bool_arr[i], smtot[i], scount[i])

        gind, drmin = distance_functions.closest_cluster_index(sxcm[i], sycm[i], szcm[i], gxcm, gycm, gzcm)
        print(i,sxcm[i], sycm[i], gxcm[gind][0], gycm[gind][0])
#i = 11; gind = 44
        print(gind[0], cloud_number[gind][0])
'''