import numpy as np
import gizmo_analysis as gizmo
import utilities as ut
from fof_analysis import fof
import pickle
import os
simname = 'm12i_res7100_mhdcv'
simdir = '/scratch/projects/xsede/GalaxiesOnFIRE/mhdcv/m12i_res7100_mhdcv/1Myr/fire3'
snapnumber = 596

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
keep = np.where((age <= .003) & (rxyz < 20) & (abs(z) < 1.5))

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

linking_length = 0.01 #10 parsec
ncut           = 10 #10 star particles

ind, xcm, ycm, zcm, mtot, grpid, r90, r50, rmax =fof.find(x0,y0,z0, b=linking_length, mass=mass0, ncut=ncut)
ngroup = len(mtot)

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
    cluster={"cluster_groupid":groupid,"no_of_star":nstar,"id":ids_in_cluster,"id_children":id_children_in_cluster}
    export_cluster.update({groupid:cluster})

print('------------------------------------------------------------------------------------------------------------------')

path="./fire3data_pkl/" #creating a path to store the data only if it does not exist
if not os.path.exists(path):
  os.makedirs(path)

file_name="fire3_clusters_"+simname+"_snapshot_"+str(snapnumber) 

with open(path+file_name, 'wb') as output:
    # Pickle dictionary using protocol 0.
    pickle.dump(export_cluster, output)
    
########################
#to test if the information was stored properly
with open(path+file_name, "rb") as fp:
    import_cluster = pickle.load(fp)
        
print("\n\n\nThe ids of the cluster with group id 1 is",import_cluster[1]["id"])
########################
print("Total clusters present is",len(import_cluster))


