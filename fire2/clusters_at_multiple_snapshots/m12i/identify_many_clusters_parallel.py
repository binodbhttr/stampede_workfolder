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
import sys

#########################################################################
#########################################################################
#load particle data here as you normally do to get the star information
#cut on Rxy, |z| and age
#run fof
#sind, sxcm, sycm, szcm, smtot, sgrpid, sr90, sr50, srmax =fof.find(x[si],y[si],z[si], b=b_kpc, mass=mass[si], ncut=ncut_min)
#srcm = np.sqrt(sxcm**2. + sycm**2.)

simname = 'm12i_res7100_mhdcv'
simdir = '/scratch/projects/xsede/GalaxiesOnFIRE/cr_suite/m12i_res7100/mhdcv/1Myr/fire2/'

argdex = int(sys.argv[1])
start  = 666+int(argdex*6)-6
finish = 666+int(argdex*6)

snapnumber = start

for s in range(snapnumber,finish,3):
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
  sife    = part['star'].prop('metallicity.si - metallicity.fe')
  cafe    = part['star'].prop('metallicity.ca - metallicity.fe')
  ofe    = part['star'].prop('metallicity.o - metallicity.fe')
  sfe    = part['star'].prop('metallicity.s - metallicity.fe')
  nefe    = part['star'].prop('metallicity.ne - metallicity.fe')
  ids      = part['star'].prop('id')
  id_child = part['star'].prop('id.child')
  age      = part['star'].prop('age')
  del part
  #select young stars within the disk for star clusters
  keep = np.where((age <= .003) & ((Rxy < 20) & (Rxy>2)) & (abs(z) < 1.5))
  #to run cluster finding on
  rxyz0     = rxyz[keep]
  del rxyz
  Rxy0      = Rxy[keep]
  del Rxy
  x0        = x[keep]
  del x
  y0        = y[keep]
  del y
  z0        = z[keep]
  del z
  mass0     = mass[keep]
  del mass
  feh0      = feh[keep]
  del feh
  mgh0=mgh[keep]
  del mgh
  ch0=ch[keep]
  del ch
  nh0=nh[keep]
  del nh
  oh0=oh[keep]
  del oh
  neh0=neh[keep]
  del neh
  sih0=sih[keep]
  del sih
  sh0=sh[keep]
  del sh
  cah0=cah[keep]
  del cah
  mgfe0=mgfe[keep]
  del mgfe
  ofe0=ofe[keep]
  del ofe
  sife0=sife[keep]
  del sife
  cafe0=cafe[keep]
  del cafe
  sfe0=sfe[keep]
  del sfe
  nefe0=nefe[keep]
  del nefe
  
  id0       = ids[keep]
  del ids
  id_child0 = id_child[keep]
  del id_child
  age0      = age[keep]
  del age
  
  linking_length = 0.004 #4 parsec (unit here is in kpc)
  ncut           = 6 #6 star particles
  
  ind, xcm, ycm, zcm, mtot, grpid, r90, r50, rmax =fof.find(x0,y0,z0, b=linking_length, mass=mass0, ncut=ncut)
  ngroup = len(mtot)
  

  export_cluster={}
  for grp_index in range(ngroup):  #iterate over each group
      cluster={}
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
      #feh_in_cluster=feh0[ind[grp_index]]
      cluster={"cluster_groupid":groupid,"no_of_star":nstar,"id":ids_in_cluster,"id_children":id_children_in_cluster,"xcm":xcm[grp_index],"ycm":ycm[grp_index],"zcm":zcm[grp_index],"mtot":mtot[grp_index],"r90":r90[grp_index],"r50":r50[grp_index],"rmax":rmax[grp_index],"x":x0[ind[grp_index]],"y":y0[ind[grp_index]],"z":z0[ind[grp_index]],"age":age0[ind[grp_index]],"feh":feh0[ind[grp_index]],"mgh":mgh0[ind[grp_index]],"ch":ch0[ind[grp_index]],"nh":nh0[ind[grp_index]],"oh":oh0[ind[grp_index]],"neh":neh0[ind[grp_index]],"sih":sih0[ind[grp_index]],"sh":sh0[ind[grp_index]],"cah":cah0[ind[grp_index]],"mgfe":mgfe0[ind[grp_index]],"ofe":ofe0[ind[grp_index]],"sife":sife0[ind[grp_index]],"sfe":sfe0[ind[grp_index]],"cafe":cafe0[ind[grp_index]],"nefe":nefe0[ind[grp_index]]}
      export_cluster.update({groupid:cluster})
  
  
  
  ######################################################################
  #######################################################################
  
  print('------------------------------------------------------------------------------------------------------------------')
  
  path="./fire2_clusterdata_pkl_b4n6/" #creating a path to store the data only if it does not exist
  if not os.path.exists(path):
    os.makedirs(path)
  
  file_name="fire2_clusters_"+simname+"_snapshot_"+str(snapnumber)+".pkl" 
  
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
  snapnumber+=3

###################################
####################################
#####################################
#######################################
#########################################
