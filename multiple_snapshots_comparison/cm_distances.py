from __main__ import *
from sl_utilities import distance_functions

mass_691=part_691['star']['mass'] #mass of all stars in snapshot 691
mass_692=part_692['star']['mass'] #mass of all stars in snapshot 692
mass_693=part_693['star']['mass'] #mass of all stars in snapshot 693
mass_694=part_694['star']['mass'] #mass of all stars in snapshot 694
mass_695=part_695['star']['mass'] #mass of all stars in snapshot 695
mass_696=part_696['star']['mass'] #mass of all stars in snapshot 696


mass_692_tracked=mass_692[ind_tracked_id_692] #we stored the mass of stars that we are tracking in snapshot 691. These are the mass we need to find CM
mass_693_tracked=mass_693[ind_tracked_id_693]
mass_694_tracked=mass_694[ind_tracked_id_694]
mass_695_tracked=mass_695[ind_tracked_id_695]
mass_696_tracked=mass_696[ind_tracked_id_696]

mass_691_R7to9=mass_691[ind_R7to9] #selected the mass values that fall in R 7 to 9
mass_691_d_xyz_lessthan1=mass_691_R7to9[ind_d_xyz_lessthan1] #selected the mass values that fall in radius less than 1 (second cut) 
mass_691_young1_d_xyz=mass_691_d_xyz_lessthan1[ind_age_young1_d_xyz] #Selected the mass values that for the selected age group (7 to 8) in log10 scale

#These are the x,y, z and mass info of the stars we are interested in the Snapshot 691.
#x_young1_d_xyz 
#y_young1_d_xyz   
#z_young1_d_xyz
#mass_696_young1_d_xyz

#These are the x,y, z and mass info of the stars we are interested in the Snapshot 692.
#x_692_tracked
#y_692_tracked
#z_692_tracked
#mass_692_tracked

#Note: The masses (mass_691_young1_d_xyz and mass_692_tracked) are equal as they are the same stars !!!!!!!!!!!!!!!!!!!

###################I put this here to calculate the cm locally
#total_mass=sum(mass_691_tracked)  
#cm = sum(mass_691_tracked[i]*x_691_tracked[i] for i in range(len(x_691_tracked)))/ total_mass
#print("This is the cm I calculated locally",cm)
####################################################################################


#Calculating the center of mass using the x y and z and mass values from the snapshot 696
xcm_691=distance_functions.cm(x_young1_d_xyz,mass_691_young1_d_xyz)
ycm_691=distance_functions.cm(y_young1_d_xyz,mass_691_young1_d_xyz)
zcm_691=distance_functions.cm(z_young1_d_xyz,mass_691_young1_d_xyz)
delta_rxyz_691=distance_functions.dr(x_young1_d_xyz,y_young1_d_xyz,z_young1_d_xyz,mass_691_young1_d_xyz)
print("\n \n This gives the position of star relative to the center of mass in in snapshot 691 \n \n",delta_rxyz_691)

#Calculating the center of mass using the x y and z and mass values from the snapshot 691
xcm_692=distance_functions.cm(x_692_tracked,mass_692_tracked)
ycm_692=distance_functions.cm(y_692_tracked,mass_692_tracked)
zcm_692=distance_functions.cm(z_692_tracked,mass_692_tracked)
delta_rxyz_692=distance_functions.dr(x_692_tracked,y_692_tracked,z_692_tracked,mass_692_tracked)
print("\n \n This gives the position of star relative to the center of mass in in snapshot 691 \n \n",delta_rxyz_691)


#Now plotting figure 10 with y limits, x limits and the center of mass
fig10 = plt.figure()
ax1 = fig10.add_subplot(211)
ax1.scatter(y_692_tracked,x_692_tracked,marker=".",s=0.5)
ax1.plot(ycm_691,xcm_692,color='red',marker=".")
ax1.set_xlabel('y')
ax1.set_ylabel('x')
ax1.set_ylim(-2,2)
ax1.set_xlim(-10,10)
ax1.set_title('Snapshot 692: young stars age 7 to 8 dec within 1kpc sphere from 0,8,0')

plt.subplots_adjust(hspace=.5)

ax2=fig10.add_subplot(212)
ax2.scatter(y_young1_d_xyz,x_young1_d_xyz,marker=".",s=0.5)
ax2.plot(ycm_691,xcm_691,color='red',marker=".")
ax2.set_xlabel('y')
ax2.set_ylabel('x')
ax2.minorticks_on()
ax2.set_ylim(-2,2)
ax2.set_xlim(-10,10)
ax2.set_title('Snapshot 691: young stars age 7 to 8 dec within 1kpc sphere from 0,8,0')
#plt.tight_layout()
fig10.savefig("./plots/scaled_young1_age7to8_snapshots691_and_691_withCM.png")









