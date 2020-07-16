from __main__ import *
from sl_utilities import distance_functions

mass_696=part['star']['mass'] #mass of all stars in snapshot 696

mass_691=part_691['star']['mass'] #mass of all stars in snapshot 691


mass_691_tracked=mass_691[ind_tracked_id_691] #we stored the mass of stars that we are tracking in snapshot 691. These are the mass we need to find CM


mass_696_R7to9=mass_696[ind_R7to9] #selected the mass values that fall in R 7 to 9

mass_696_d_xyz_lessthan1=mass_696_R7to9[ind_d_xyz_lessthan1] #selected the mass values that fall in radius less than 1 (second cut) 
mass_696_young1_d_xyz=mass_696_d_xyz_lessthan1[ind_age_young1_d_xyz] #Selected the mass values that for the selected age group (7 to 8) in log10 scale

#These are the x,y, z and mass info of the stars we are interested in the Snapshot 696.
#x_young1_d_xyz 
#y_young1_d_xyz   
#z_young1_d_xyz
#mass_696_young1_d_xyz

#These are the x,y, z and mass info of the stars we are interested in the Snapshot 691.
#x_691_tracked
#y_691_tracked
#z_691_tracked
#mass_691_tracked


#Calculating the center of mass using the x y and z and mass values from the snapshot 696
xcm_696=distance_functions.cm(x_young1_d_xyz, mass_696_young1_d_xyz)
ycm_696=distance_functions.cm(y_young1_d_xyz, mass_696_young1_d_xyz)
zcm_696=distance_functions.cm(z_young1_d_xyz, mass_696_young1_d_xyz)
delta_rxyz_696=distance_functions.dr(xcm_696,ycm_696,zcm_696,mass_696_young1_d_xyz)


#Calculating the center of mass using the x y and z and mass values from the snapshot 691
xcm_691=distance_functions.cm(x_691_tracked, mass_691_tracked)
ycm_691=distance_functions.cm(y_691_tracked, mass_691_tracked)
zcm_691=distance_functions.cm(z_691_tracked, mass_691_tracked)
delta_rxyz_691=distance_functions.dr(x_691_tracked,y_691_tracked,z_691_tracked,mass_691_tracked)


#Now plotting figure 10 with y limits, x limits and the center of mass
fig10 = plt.figure()
ax1 = fig10.add_subplot(211)
ax1.scatter(y_691_tracked,x_691_tracked,marker=".",s=0.5)
ax1.plot(ycm_691,xcm_691)
ax1.set_xlabel('y')
ax1.set_ylabel('x')
ax1.set_ylim(-2,2)
ax1.set_xlim(-10,10)
ax1.set_title('Snapshot 691: young stars age 7 to 8 dec within 1kpc sphere from 0,8,0')

plt.subplots_adjust(hspace=.5)

ax2=fig10.add_subplot(212)
ax2.scatter(y_young1_d_xyz,x_young1_d_xyz,marker=".",s=0.5)
ax2.plot(ycm_696,xcm_696)
ax2.set_xlabel('y')
ax2.set_ylabel('x')
ax2.minorticks_on()
ax2.set_ylim(-2,2)
ax2.set_xlim(-10,10)
ax2.set_title('Snapshot 696: young stars age 7 to 8 dec within 1kpc sphere from 0,8,0')
#plt.tight_layout()
fig10.savefig("./plots/scaled_young1_age7to8_snapshots691_and_696_withCM.png")









