from __main__ import *
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
id_young1_d_xyz=id_d_xyz_lessthan1[ind_age_young1_d_xyz]
id_child_young1_d_xyz=id_child_d_xyz_lessthan1[ind_age_young1_d_xyz]

print("The total no. of young stars (7 to 8) in the sperical radius <1 from a point (0,8,0) in Snapshot 696 is",len(id_young1_d_xyz))

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
ax2.scatter(y_young1_d_xyz,x_young1_d_xyz,marker=".",s=0.5)
ax2.set_xlabel('y')
ax2.set_ylabel('x')
ax2.minorticks_on()
ax2.set_title('in 2D (x vs y)')
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


ind_tracked_id_692=matchids(id_young1_d_xyz,id_child_young1_d_xyz,id_692,id_child_692)
ind_tracked_id_693=matchids(id_young1_d_xyz,id_child_young1_d_xyz,id_693,id_child_693)
ind_tracked_id_694=matchids(id_young1_d_xyz,id_child_young1_d_xyz,id_694,id_child_694)
ind_tracked_id_695=matchids(id_young1_d_xyz,id_child_young1_d_xyz,id_695,id_child_695)
ind_tracked_id_696=matchids(id_young1_d_xyz,id_child_young1_d_xyz,id_696,id_child_696)




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



#Now plotting ind 2D the stars in two snapshots
fig8 = plt.figure()
#fig8.suptitle("Stars aged 7 to 8 dec within 1 from 0,8,0 \n \n")
ax1 = fig8.add_subplot(321)
ax1.scatter(y_young1_d_xyz,x_young1_d_xyz,marker=".",s=0.5)
ax1.set_xlabel('y')
ax1.set_ylabel('x')
ax1.minorticks_on()
ax1.set_ylim(-2,2)
ax1.set_xlim(-10,10)
ax1.set_title('Snapshot 691')

ax2=fig8.add_subplot(322)
ax2.scatter(y_692_tracked,x_692_tracked,marker=".",s=0.5)
ax2.set_xlabel('y')
ax2.set_ylabel('x')
ax2.minorticks_on()
ax2.set_ylim(-2,2)
ax2.set_xlim(-10,10)
ax2.set_title('Snapshot 692')


ax3=fig8.add_subplot(323)
ax3.scatter(y_693_tracked,x_693_tracked,marker=".",s=0.5)
ax3.set_xlabel('y')
ax3.set_ylabel('x')
ax3.set_ylim(-2,2)
ax3.set_xlim(-10,10)
ax3.minorticks_on()
ax3.set_title('Snapshot 693')

#plt.subplots_adjust(hspace=.5)

ax4=fig8.add_subplot(324)
ax4.scatter(y_694_tracked,x_694_tracked,marker=".",s=0.5)
ax4.set_xlabel('y')
ax4.set_ylabel('x')
ax4.minorticks_on()
ax4.set_ylim(-2,2)
ax4.set_xlim(-10,10)
ax4.set_title('Snapshot 694')


ax5=fig8.add_subplot(325)
ax5.scatter(y_695_tracked,x_695_tracked,marker=".",s=0.5)
ax5.set_xlabel('y')
ax5.set_ylabel('x')
ax5.minorticks_on()
ax5.set_ylim(-2,2)
ax5.set_xlim(-10,10)
ax5.set_title('Snapshot 695')

#plt.subplots_adjust(hspace=.5)

ax6=fig8.add_subplot(326)
ax6.scatter(y_696_tracked,x_696_tracked,marker=".",s=0.5)
ax6.set_xlabel('y')
ax6.set_ylabel('x')
ax6.minorticks_on()
ax6.set_ylim(-2,2)
ax6.set_xlim(-10,10)
ax6.set_title('Snapshot 696')

plt.tight_layout()
fig8.savefig("./plots/young1_age7to8_snapshots691to696.png")
