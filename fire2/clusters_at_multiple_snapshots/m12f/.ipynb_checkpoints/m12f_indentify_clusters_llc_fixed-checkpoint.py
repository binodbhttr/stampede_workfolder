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

simname = 'm12f_res7100_mhdcv'
#simdir = '/scratch/projects/xsede/GalaxiesOnFIRE/cr_suite/m12f_res7100/mhdcv/1Myr/'
simdir = '/scratch/projects/xsede/GalaxiesOnFIRE/cr_suite/m12f_r7100/mhdcv/1Myr/'

snapnumber = 588



## NOTE: Inputs can be either position or velocity in the "Translated Box Coordinate Frame"
## If inputs are positions in the "Translated Box Coordinate Frame", outputs will be positions in the "Loebman Lab Coordinate Frame"
## If inputs are velocities in the "Translated Box Coordinate Frame", outputs will be velocities in the "Loebman Lab Coordinate Frame"

def ll_coordinate_transformation(rotation_matrix, x_translated_box_coordinate, y_translated_box_coordinate, z_translated_box_coordinate):
    
    # identifying the elements of the True Rotation Matrix
    rot_mat11 = rotation_matrix[0][0][0]
    rot_mat12 = rotation_matrix[0][0][1]
    rot_mat13 = rotation_matrix[0][0][2]
    
    rot_mat21 = rotation_matrix[0][1][0]
    rot_mat22 = rotation_matrix[0][1][1]
    rot_mat23 = rotation_matrix[0][1][2]
    
    rot_mat31 = rotation_matrix[0][2][0]
    rot_mat32 = rotation_matrix[0][2][1]
    rot_mat33 = rotation_matrix[0][2][2]
    
    ''' 
        Finding the three "True Rotation Angles": rot_angle_0, rot_angle_1 and rot_angle_2
        
        Note that there are eight possible sets of rotation angles.

        We use the following approach to know which set of rotation angles is the "True Rotation Angles" set.
        For each of the eight possible sets of rotation angles, 
        we compute the rotation matrix (let's call it the Recovered Rotation Matrix), 
        and we check which of these "Recovered Rotation Matrices" matches the "True Rotation Matrix."
        
    '''
    
    rot_angle_0 = [] 
    rot_angle_1 = []
    rot_angle_2 = []

    for k in range (2):
        for j in range (2):
            for i in range (2):
                recovered_rotation_angles = [0, 0, 0]
                recovered_rotation_angles[0] =  np.arctan(-rot_mat32/rot_mat33) + (i*np.pi)
                recovered_rotation_angles[1] =  np.arctan(rot_mat31/(np.sqrt(1- (rot_mat31**2)))) + (j*np.pi)
                recovered_rotation_angles[2] =  np.arctan(-rot_mat21/rot_mat11) + (k*np.pi)

                recovered_rot_mat11 = np.cos(recovered_rotation_angles[1]) * np.cos(recovered_rotation_angles[2])
                recovered_rot_mat12 = np.cos(recovered_rotation_angles[0]) * np.sin(recovered_rotation_angles[2]) + np.sin(recovered_rotation_angles[0]) * np.sin(recovered_rotation_angles[1]) * np.cos(recovered_rotation_angles[2])
                recovered_rot_mat13 = np.sin(recovered_rotation_angles[0]) * np.sin(recovered_rotation_angles[2]) - np.cos(recovered_rotation_angles[0]) * np.sin(recovered_rotation_angles[1]) * np.cos(recovered_rotation_angles[2])
                recovered_rot_mat21 = -np.cos(recovered_rotation_angles[1]) * np.sin(recovered_rotation_angles[2])
                recovered_rot_mat22 = np.cos(recovered_rotation_angles[0]) * np.cos(recovered_rotation_angles[2]) - np.sin(recovered_rotation_angles[0]) * np.sin(recovered_rotation_angles[1]) * np.sin(recovered_rotation_angles[2])
                recovered_rot_mat23 = np.sin(recovered_rotation_angles[0]) * np.cos(recovered_rotation_angles[2]) + np.cos(recovered_rotation_angles[0]) * np.sin(recovered_rotation_angles[1]) * np.sin(recovered_rotation_angles[2])
                recovered_rot_mat31 = np.sin(recovered_rotation_angles[1])
                recovered_rot_mat32 = -np.sin(recovered_rotation_angles[0]) * np.cos(recovered_rotation_angles[1])
                recovered_rot_mat33 = np.cos(recovered_rotation_angles[0]) * np.cos(recovered_rotation_angles[1])

                if ((np.round(recovered_rot_mat11.astype('float64'), 3) == np.round(rot_mat11.astype('float64'), 3)) & 
                    (np.round(recovered_rot_mat12.astype('float64'), 3) == np.round(rot_mat12.astype('float64'), 3)) & 
                    (np.round(recovered_rot_mat13.astype('float64'), 3) == np.round(rot_mat13.astype('float64'), 3)) &
                    (np.round(recovered_rot_mat21.astype('float64'), 3) == np.round(rot_mat21.astype('float64'), 3)) &
                    (np.round(recovered_rot_mat22.astype('float64'), 3) == np.round(rot_mat22.astype('float64'), 3)) &
                    (np.round(recovered_rot_mat23.astype('float64'), 3) == np.round(rot_mat23.astype('float64'), 3)) &
                    (np.round(recovered_rot_mat31.astype('float64'), 3) == np.round(rot_mat31.astype('float64'), 3)) &
                    (np.round(recovered_rot_mat32.astype('float64'), 3) == np.round(rot_mat32.astype('float64'), 3)) &
                    (np.round(recovered_rot_mat33.astype('float64'), 3) == np.round(rot_mat33.astype('float64'), 3))):

                    rot_angle_0.append(recovered_rotation_angles[0])
                    rot_angle_1.append(recovered_rotation_angles[1])
                    rot_angle_2.append(recovered_rotation_angles[2])

                    break
                    del(recovered_rotation_angles)
                    
    ''' 
        Constructing the "Loebman Lab Transformation Matrix": ll_tans_matij   
    '''
    
    ll_trans_mat11 = np.cos(rot_angle_1)
    ll_trans_mat12 = np.sin(rot_angle_0) * np.sin(rot_angle_1)
    ll_trans_mat13 = - np.cos(rot_angle_0) * np.sin(rot_angle_1)
    ll_trans_mat21 = np.zeros(1)
    ll_trans_mat22 = np.cos(rot_angle_0)
    ll_trans_mat23 = np.sin(rot_angle_0) 
    ll_trans_mat31 = np.sin(rot_angle_1)
    ll_trans_mat32 = -np.sin(rot_angle_0) * np.cos(rot_angle_1)
    ll_trans_mat33 = np.cos(rot_angle_0) * np.cos(rot_angle_1)
    
    '''
        Finding coordinate of particles in the "Loebman Lab Coordinate System": x_ll, y_ll and z_ll
    '''
    
    x_ll = (ll_trans_mat11*x_translated_box_coordinate) + (ll_trans_mat12*y_translated_box_coordinate) + (ll_trans_mat13*z_translated_box_coordinate)
    y_ll = (ll_trans_mat21*x_translated_box_coordinate) + (ll_trans_mat22*y_translated_box_coordinate) + (ll_trans_mat23*z_translated_box_coordinate)
    z_ll = (ll_trans_mat31*x_translated_box_coordinate) + (ll_trans_mat32*y_translated_box_coordinate) + (ll_trans_mat33*z_translated_box_coordinate)
    
    position_ll = np.asarray((x_ll, y_ll, z_ll)).T
    
    return position_ll   


for s in range(snapnumber,696,3):
    part = gizmo.io.Read.read_snapshots(['all'],'snapshot_index', snapnumber, simulation_name=simname, simulation_directory=simdir, assign_hosts_rotation=True, assign_hosts=True)  

    rxyz     = part['star'].prop('host.distance.total')
    Rxy      = part['star'].prop('host.distance.principal.cylindrical')[:,0]
    #x        = part['star'].prop('host.distance.principal')[:,0]
    #y        = part['star'].prop('host.distance.principal')[:,1]
    #z        = part['star'].prop('host.distance.principal')[:,2] 
    
    x_box=part['star'].prop('host.distance')[:,0]
    y_box=part['star'].prop('host.distance')[:,1]
    z=part['star'].prop('host.distance.principal')[:,2] 
    z_box=part['star'].prop('host.distance')[:,2]
    
    rm = part.host['rotation']     

    x_llc=ll_coordinate_transformation(rm, x_box, y_box, z_box)[:,0]
    y_llc=ll_coordinate_transformation(rm, x_box, y_box, z_box)[:,1]
    z_llc=ll_coordinate_transformation(rm, x_box, y_box, z_box)[:,2]
    
    mass     = part['star'].prop('mass')

    feh      = part['star'].prop('metallicity.fe')
    mgh      = part['star'].prop('metallicity.mg')
    cah      = part['star'].prop('metallicity.ca')
    sh      = part['star'].prop('metallicity.s')
    sih      = part['star'].prop('metallicity.si')
    neh      = part['star'].prop('metallicity.ne')
    oh      = part['star'].prop('metallicity.o')
    nh      = part['star'].prop('metallicity.n')
    ch      = part['star'].prop('metallicity.c')

    #metal_list=["mgh","cah","sh","sih","neh","oh","nh","ch"]
    #metal_over_fe_list=['mgfe','cafe','sfe','sife','nefe','ofe','nfe','cfe']

    mgfe  = part['star'].prop('metallicity.mg - metallicity.fe')
    cafe  = part['star'].prop('metallicity.ca - metallicity.fe')
    sfe   = part['star'].prop('metallicity.s - metallicity.fe')
    sife  = part['star'].prop('metallicity.si - metallicity.fe')
    nefe  = part['star'].prop('metallicity.ne - metallicity.fe')
    ofe   = part['star'].prop('metallicity.o - metallicity.fe')
    nfe   = part['star'].prop('metallicity.n - metallicity.fe')
    cfe   = part['star'].prop('metallicity.c - metallicity.fe')

    ids      = part['star'].prop('id')
    id_child = part['star'].prop('id.child')
    age      = part['star'].prop('age')
    #select young stars within the disk for star clusters
    #keep = np.where((age <= .003) & ((Rxy < 20) & (Rxy>2)) & (abs(z) < 1.5))
    keep = np.where((age <= .003) & (rxyz < 20) & (abs(z) < 1.5))

    #to run cluster finding on
    rxyz0     = rxyz[keep]
    Rxy0      = Rxy[keep]
    
    x0_llc    = x_llc[keep]
    y0_llc    = y_llc[keep]
    z0_llc    = z_llc[keep]
    
    #x0       = x[keep]
    #y0       = y[keep]
    #z0       = z[keep]
    mass0     = mass[keep]
    
    feh0=feh[keep]
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
    sife0=sife[keep]
    cafe0=cafe[keep]
    sfe0=sfe[keep]
    nefe0=nefe[keep]
    nfe0=nfe[keep]
    cfe0=cfe[keep]

    id0=ids[keep]
    id_child0=id_child[keep]
    age0=age[keep]

    linking_length = 0.004 #4 parsec (unit here is in kpc)
    ncut           = 5 #4 star particles

    #ind, xcm, ycm, zcm, mtot, grpid, r90, r50, rmax =fof.find(x0,y0,z0, b=linking_length, mass=mass0, ncut=ncut)
    ind, xcm_llc, ycm_llc, zcm_llc, mtot, grpid, r90, r50, rmax =fof.find(x0_llc,y0_llc,z0_llc, b=linking_length, mass=mass0, ncut=ncut)

    ngroup = len(mtot)
    
    if mtot[0]<0:
        print("No clusters found at snapshot: ",snapnumber)
        snapnumber+=3
        continue
    
    export_cluster_llc={}
    for grp_index in range(ngroup):  #iterate over each group
        cluster={}
        ids_in_cluster = id0[ind[grp_index]]  #these are the star particle ids in each cluster
        id_children_in_cluster = id_child0[ind[grp_index]]
        age=age0[ind[grp_index]]
        nstar = len(ids_in_cluster)
        groupid=grpid[grp_index]
        print('------------------------------------------------------------------------------------------------------------------')
        print('grpid, nstar, xcm (kpc), ycm (kpc), zcm (kpc), mtot (msun), rmax (pc)')
        print('%s     %i     %.4f     %.4f    %.4f     %.2e     %.1f ' % (grpid[grp_index], nstar, xcm_llc[grp_index], ycm_llc[grp_index], zcm_llc[grp_index], mtot[grp_index], 1000*rmax[grp_index]))
        #print('%s     %i     %.4f     %.4f    %.4f     %.2e     %.1f ' % (grpid[grp_index], nstar, xcm[grp_index], ycm[grp_index], zcm[grp_index], mtot[grp_index], 1000*rmax[grp_index]))
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
        cluster={"cluster_groupid":groupid,"no_of_star":nstar,"id":ids_in_cluster,
        "id_children":id_children_in_cluster,"xcm_llc":xcm_llc[grp_index],"ycm_llc":ycm_llc[grp_index],
        "zcm_llc":zcm_llc[grp_index],"mtot":mtot[grp_index],"r90":r90[grp_index],"r50":r50[grp_index],
        "rmax":rmax[grp_index],"x_llc":x0_llc[ind[grp_index]],"y_llc":y0_llc[ind[grp_index]],"z_llc":z0_llc[ind[grp_index]],
        "age":age0[ind[grp_index]],"feh":feh0[ind[grp_index]],"mgh":mgh0[ind[grp_index]],
        "ch":ch0[ind[grp_index]],"nh":nh0[ind[grp_index]],"oh":oh0[ind[grp_index]],
        "neh":neh0[ind[grp_index]],"sih":sih0[ind[grp_index]],"sh":sh0[ind[grp_index]],
        "cah":cah0[ind[grp_index]],"mgfe":mgfe0[ind[grp_index]],"ofe":ofe0[ind[grp_index]],
        "sife":sife0[ind[grp_index]],"sfe":sfe0[ind[grp_index]],"cafe":cafe0[ind[grp_index]],
        "nefe":nefe0[ind[grp_index]],"nfe":nfe0[ind[grp_index]],"cfe":cfe0[ind[grp_index]]}
        export_cluster_llc.update({groupid:cluster})

    ######################################################################
    #######################################################################

    print('------------------------------------------------------------------------------------------------------------------')
    work_data_path="/home1/07428/binod/work2/data/clusters_from_multiple_snapshots/"
    path=work_data_path+"fire2_m12f_clusterdata_pkl_b4n5_llc_fixed/" #creating a path to store the data only if it does not exist
    if not os.path.exists(path):
        os.makedirs(path)

    file_name="fire2_clusters_"+simname+"_snapshot_"+str(snapnumber)+"_b4n5_llc_fixed.pkl" 

    with open(path+file_name, 'wb') as output:
        # Pickle dictionary using protocol 0.
        pickle.dump(export_cluster_llc, output)

    ########################
    #to test if the information was stored properly
    with open(path+file_name, "rb") as fp:
        import_cluster = pickle.load(fp)

    print("\n\n\nThe ids of the cluster with group id 1 is",import_cluster[1]["id"])
    ########################
    print("Total clusters present is",len(import_cluster))
    
    del rxyz0,Rxy0,x0_llc,y0_llc,z0_llc,mass0,feh0,mgh0,ch0,nh0,oh0,neh0,sih0,sh0,cah0,mgfe0,ofe0,
    sife0,cafe0,sfe0,nefe0,nfe0,cfe0,id0,id_child0,age0,keep
    
    snapnumber+=3

    ###################################
    ####################################
    #####################################
    #######################################
    #########################################