from __main__ import *
from sl_utilities import distinct_colours as dc
from sl_utilities import distance_functions
print("Now running the part to make selections from a different file !!!")

ind_R7to9=np.where((R_691>7) & (R_691<9)) #Select those stars with r greater than 7 and less than 9

x_R7to9=x_691[ind_R7to9] #find the x components only of the slection
y_R7to9=y_691[ind_R7to9] #find the y components only of the slection
z_R7to9=z_691[ind_R7to9]
age_R7to9=np.log10(age_691[ind_R7to9])
print("\n Note ages are already in log scale now !!!!! These are the ages from the snapshot (696) ",age_R7to9)
R_R7to9=R_691[ind_R7to9] #store the cylindrical radius in array Note: the values should be between 7 and 9
id_R7to9=id_691[ind_R7to9]
id_child_R7to9=id_child_691[ind_R7to9]


d_xyz=np.sqrt(x_R7to9**2+(y_R7to9-8)**2+z_R7to9**2)   #calculate radius of all stars in xyz plane form a point 
print(" \n Calculated the spherical radius for all points that lie inside the cylindrical radius 7 and 9") #inside the selection we made previously. Here the point is (0,8)
ind_d_xyz_lessthan1=np.where((d_xyz>0)&(d_xyz<1))
print("\n Selected the indices of those objects that are within 1 kpc spherical radius",ind_d_xyz_lessthan1)
d_xyz_lessthan1=d_xyz[ind_d_xyz_lessthan1] #here we stored those distances that are less than 1 in xyz plane
print("\n Total objects in that radius is",len(d_xyz_lessthan1)) #count the no. of such stars


x_d_xyz_lessthan1=x_R7to9[ind_d_xyz_lessthan1] 
y_d_xyz_lessthan1=y_R7to9[ind_d_xyz_lessthan1]
z_d_xyz_lessthan1=z_R7to9[ind_d_xyz_lessthan1]
R_d_xyz_lessthan1=R_R7to9[ind_d_xyz_lessthan1]
age_d_xyz_lessthan1=age_R7to9[ind_d_xyz_lessthan1] #ages of star in the region R7to9 that are withtin d<1 in xyz plane
id_d_xyz_lessthan1=id_R7to9[ind_d_xyz_lessthan1]
id_child_d_xyz_lessthan1=id_child_R7to9[ind_d_xyz_lessthan1]

print("These are the ids of the objects that are within a spherical radius of 1kpc from a point (0,8,0)",id_d_xyz_lessthan1)

#Now we are going to make age selections

ind_age_young1_d_xyz=np.where((age_d_xyz_lessthan1>7)&(age_d_xyz_lessthan1<8)) #select stars of the age between 7 and 8 in log10                       
x_young1_d_xyz=x_d_xyz_lessthan1[ind_age_young1_d_xyz] 
y_young1_d_xyz=y_d_xyz_lessthan1[ind_age_young1_d_xyz]   
z_young1_d_xyz=z_d_xyz_lessthan1[ind_age_young1_d_xyz] 
age_young1_d_xyz=age_d_xyz_lessthan1[ind_age_young1_d_xyz]   
#id_young1_d_xyz=id_d_xyz_lessthan1[ind_age_young1_d_xyz]
#id_child_young1_d_xyz=id_child_d_xyz_lessthan1[ind_age_young1_d_xyz]
#created a fake ID here
id_young1_d_xyz=np.array([46025169,55215865,14548425,27846424,21612608,26482714,45183866,23320885,13309616,57262468,19850157,17172767])
id_child_young1_d_xyz=np.array([0,0,0,0,0,0,0,0,1,0,0,0])
sortind=np.argsort(id_young1_d_xyz)
id_young1_d_xyz_sorted=id_young1_d_xyz[sortind]
id_child_young1_d_xyz_sorted=id_child_young1_d_xyz[sortind]

print("The total no. of stars in this cluster is",len(id_young1_d_xyz))

#Now plotting that selection
fig7 = plt.figure()
ax1 = fig7.add_subplot(211, projection='3d')
ax1.scatter(x_young1_d_xyz,y_young1_d_xyz,z_young1_d_xyz,marker=".",s=0.5)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('z')
ax1.set_zlim(-2,2)
ax1.set_title('young stars age 7 to 8 dec within 1kpc sphere from 0,8,0')

plt.subplots_adjust(hspace=.5)

ax2=fig7.add_subplot(212)
ax2.scatter(x_young1_d_xyz,y_young1_d_xyz,marker=".",s=0.5)
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.minorticks_on()
ax2.set_title('in 2D (y vs x)')
#plt.tight_layout()
fig7.savefig("./plots/verytiny_cut_1_young1_age7to8inlog10scale.png")



###Now let's find the matching ids in the next snapshot using a function
def matchids(id_current,id_child_current,id_next,id_child_next):  #this function returns the index of the ids in the next snapshot that match with the current
  ind=np.array(0)
  for i in range(len(id_current)):
    match=np.where((id_next==id_current[i])&(id_child_next==id_child_current[i]))
    print("\nMatched the id",id_next[match])
    ind=np.append(ind,match)
  ind_tracked_id_next=ind[1:len(ind)]  #The extra element in the beginning is removed by this process
  print("\nThese are the indices of the ids that matched in current snapshot\n",ind_tracked_id_next)
  return ind_tracked_id_next  

ind_tracked_id_691=matchids(id_young1_d_xyz_sorted,id_child_young1_d_xyz_sorted,id_691,id_child_691)
ind_tracked_id_692=matchids(id_young1_d_xyz_sorted,id_child_young1_d_xyz_sorted,id_692,id_child_692)
ind_tracked_id_693=matchids(id_young1_d_xyz_sorted,id_child_young1_d_xyz_sorted,id_693,id_child_693)
ind_tracked_id_694=matchids(id_young1_d_xyz_sorted,id_child_young1_d_xyz_sorted,id_694,id_child_694)
ind_tracked_id_695=matchids(id_young1_d_xyz_sorted,id_child_young1_d_xyz_sorted,id_695,id_child_695)
ind_tracked_id_696=matchids(id_young1_d_xyz_sorted,id_child_young1_d_xyz_sorted,id_696,id_child_696)

x_691_tracked=x_691[ind_tracked_id_691]
y_691_tracked=y_691[ind_tracked_id_691]
z_691_tracked=z_691[ind_tracked_id_691]


x_692_tracked=x_692[ind_tracked_id_692]
y_692_tracked=y_692[ind_tracked_id_692]
z_692_tracked=z_692[ind_tracked_id_692]
 
#########
x_693_tracked=x_693[ind_tracked_id_693]
y_693_tracked=y_693[ind_tracked_id_693]
z_693_tracked=z_693[ind_tracked_id_693]

#########
x_694_tracked=x_694[ind_tracked_id_694]
y_694_tracked=y_694[ind_tracked_id_694]
z_694_tracked=z_694[ind_tracked_id_694]

#########
x_695_tracked=x_695[ind_tracked_id_695]
y_695_tracked=y_695[ind_tracked_id_695]
z_695_tracked=z_695[ind_tracked_id_695]

#########
x_696_tracked=x_696[ind_tracked_id_696]
y_696_tracked=y_696[ind_tracked_id_696]
z_696_tracked=z_696[ind_tracked_id_696]


#detting distinct color for each star particle
colors = dc.get_distinct(len(x_691_tracked))


#we had these values from the data 
xcm_691=8.2097
ycm_691=-0.9572
ycm_691=ycm_691*-1
xcm_691=xcm_691*-1
rmax_691=11.4
rmax_kpc_691=0.0114


ymax_691=ycm_691+3*rmax_kpc_691
ymin_691=ycm_691-3*rmax_kpc_691
xmax_691=xcm_691+3*rmax_kpc_691
xmin_691=xcm_691-3*rmax_kpc_691



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


#Calculating the center of mass using the x y and z and mass values from the snapshot 692
xcm_692=distance_functions.cm(x_692_tracked,mass_692_tracked)
ycm_692=distance_functions.cm(y_692_tracked,mass_692_tracked)
zcm_692=distance_functions.cm(z_692_tracked,mass_692_tracked)
delta_rxyz_692=distance_functions.dr(x_692_tracked,y_692_tracked,z_692_tracked,mass_692_tracked)
rmax_692=distance_functions.drmax(x_692_tracked,y_692_tracked,z_692_tracked,mass_692_tracked)  
print("\n This is the rmax in snapshot 692 in kpc",rmax_692)
ymax_692=ycm_692+3*rmax_692
ymin_692=ycm_692-3*rmax_692
xmax_692=xcm_692+3*rmax_692
xmin_692=xcm_692-3*rmax_692


#Calculating the center of mass using the x y and z and mass values from the snapshot 693
xcm_693=distance_functions.cm(x_693_tracked,mass_693_tracked)
ycm_693=distance_functions.cm(y_693_tracked,mass_693_tracked)
zcm_693=distance_functions.cm(z_693_tracked,mass_693_tracked)
delta_rxyz_693=distance_functions.dr(x_693_tracked,y_693_tracked,z_693_tracked,mass_693_tracked)
rmax_693=distance_functions.drmax(x_693_tracked,y_693_tracked,z_693_tracked,mass_693_tracked)  
print("\n This is the rmax in snapshot 693 in kpc",rmax_693)
ymax_693=ycm_693+3*rmax_693
ymin_693=ycm_693-3*rmax_693
xmax_693=xcm_693+3*rmax_693
xmin_693=xcm_693-3*rmax_693

#Calculating the center of mass using the x y and z and mass values from the snapshot 694
xcm_694=distance_functions.cm(x_694_tracked,mass_694_tracked)
ycm_694=distance_functions.cm(y_694_tracked,mass_694_tracked)
zcm_694=distance_functions.cm(z_694_tracked,mass_694_tracked)
delta_rxyz_694=distance_functions.dr(x_694_tracked,y_694_tracked,z_694_tracked,mass_694_tracked)
rmax_694=distance_functions.drmax(x_694_tracked,y_694_tracked,z_694_tracked,mass_694_tracked)  
print("\n This is the rmax in snapshot 694 in kpc",rmax_694)
ymax_694=ycm_694+3*rmax_694
ymin_694=ycm_694-3*rmax_694
xmax_694=xcm_694+3*rmax_694
xmin_694=xcm_694-3*rmax_694


#Calculating the center of mass using the x y and z and mass values from the snapshot 695
xcm_695=distance_functions.cm(x_695_tracked,mass_695_tracked)
ycm_695=distance_functions.cm(y_695_tracked,mass_695_tracked)
zcm_695=distance_functions.cm(z_695_tracked,mass_695_tracked)
delta_rxyz_695=distance_functions.dr(x_695_tracked,y_695_tracked,z_695_tracked,mass_695_tracked)
rmax_695=distance_functions.drmax(x_695_tracked,y_695_tracked,z_695_tracked,mass_695_tracked)  
print("\n This is the rmax in snapshot 695 in kpc",rmax_695)
ymax_695=ycm_695+3*rmax_695
ymin_695=ycm_695-3*rmax_695
xmax_695=xcm_695+3*rmax_695
xmin_695=xcm_695-3*rmax_695

#Calculating the center of mass using the x y and z and mass values from the snapshot 696
xcm_696=distance_functions.cm(x_696_tracked,mass_696_tracked)
ycm_696=distance_functions.cm(y_696_tracked,mass_696_tracked)
zcm_696=distance_functions.cm(z_696_tracked,mass_696_tracked)
delta_rxyz_696=distance_functions.dr(x_696_tracked,y_696_tracked,z_696_tracked,mass_696_tracked)
rmax_696=distance_functions.drmax(x_696_tracked,y_696_tracked,z_696_tracked,mass_696_tracked)  
print("\n This is the rmax in snapshot 696 in kpc",rmax_696)
ymax_696=ycm_696+3*rmax_696
ymin_696=ycm_696-3*rmax_696
xmax_696=xcm_696+3*rmax_696
xmin_696=xcm_696-3*rmax_696



#Now plotting ind 2D the stars in two snapshots
fig8 = plt.figure()
#fig8.suptitle("Stars aged 7 to 8 dec within 1 from 0,8,0 \n \n")

ax1 = fig8.add_subplot(321)
for i in range(len(x_691_tracked)):
  ax1.scatter([-1*x_691_tracked[i]],[-1*y_691_tracked[i]],color=colors[i],marker=".",s=10)
ax1.plot(xcm_691,ycm_691,color='red',marker=".")
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.minorticks_on()
ax1.set_xlim(xmin_691,xmax_691)
ax1.set_ylim(ymin_691,ymax_691)
ax1.set_title('Snapshot 691')

ax2=fig8.add_subplot(322)
for i in range(len(x_692_tracked)):
  ax2.scatter([x_692_tracked[i]],[y_692_tracked[i]],color=colors[i],marker=".",s=10)
ax2.plot(xcm_692,ycm_692,color='red',marker=".")
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.minorticks_on()
ax2.set_xlim(xmin_692,xmax_692)
ax2.set_ylim(ymin_692,ymax_692)
ax2.set_title('Snapshot 692')


ax3=fig8.add_subplot(323)
for i in range(len(x_693_tracked)):
  ax3.scatter([x_693_tracked[i]],[y_693_tracked[i]],color=colors[i],marker=".",s=10)
ax3.plot(xcm_693,ycm_693,color='red',marker=".")
ax3.set_xlabel('x')
ax3.set_ylabel('y')
ax3.set_xlim(xmin_693,xmax_693)
ax3.set_ylim(ymin_693,ymax_693)
ax3.minorticks_on()
ax3.set_title('Snapshot 693')

#plt.subplots_adjust(hspace=.5)

ax4=fig8.add_subplot(324)
for i in range(len(x_694_tracked)):
  ax4.scatter([x_694_tracked[i]],[y_694_tracked[i]],color=colors[i],marker=".",s=10)
ax4.plot(xcm_694,ycm_694,color='red',marker=".")
ax4.set_xlabel('x')
ax4.set_ylabel('y')
ax4.minorticks_on()
ax4.set_xlim(xmin_694,xmax_694)
ax4.set_ylim(ymin_694,ymax_694)
ax4.set_title('Snapshot 694')


ax5=fig8.add_subplot(325)
for i in range(len(x_695_tracked)):  
  ax5.scatter([x_695_tracked[i]],[y_695_tracked[i]],color=colors[i],marker=".",s=10)
ax5.plot(xcm_695,ycm_695,color='red',marker=".")
ax5.set_xlabel('x')
ax5.set_ylabel('y')
ax5.minorticks_on()
ax5.set_xlim(xmin_695,xmax_695)
ax5.set_ylim(ymin_695,ymax_695)
ax5.set_title('Snapshot 695')

#plt.subplots_adjust(hspace=.5)

ax6=fig8.add_subplot(326)
for i in range(len(x_696_tracked)):  
  ax6.scatter([x_696_tracked[i]],[y_696_tracked[i]],color=colors[i],marker=".",s=10)
ax6.plot(xcm_696,ycm_696,color='red',marker=".")
ax6.set_xlabel('x')
ax6.set_ylabel('y')
ax6.minorticks_on()
ax6.set_xlim(xmin_696,xmax_696)
ax6.set_ylim(ymin_696,ymax_696)
ax6.set_title('Snapshot 696')

plt.tight_layout()
fig8.savefig("./plots/cluster_group_15_snapshots691to696.png")
