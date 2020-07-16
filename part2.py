from __main__ import *
print("Now running the part to make selections from a different file !!!")

ind_R7to9=np.where((R>7) & (R<9)) #Select those stars with r greater than 7 and less than 9

x_R7to9=x[ind_R7to9] #find the x components only of the slection
y_R7to9=y[ind_R7to9] #find the y components only of the slection
z_R7to9=z[ind_R7to9]
age_R7to9=np.log10(age[ind_R7to9])
print("Note ages are already in log scale now !!!!! These are the ages from the snapshot (696) ",age_R7to9)
R_R7to9=R[ind_R7to9] #store the cylindrical radius in array Note: the values should be between 7 and 9
id_R7to9=id[ind_R7to9]
id_child_R7to9=id_child[ind_R7to9]


d_xyz=np.sqrt(x_R7to9**2+(y_R7to9-8)**2+z_R7to9**2)   #calculate radius of all stars in xyz plane form a point 
print("Calculated the spherical radius for all points that lie inside the cylindrical radius 7 and 9") #inside the selection we made previously. Here the point is (0,8)
ind_d_xyz_lessthan1=np.where(d_xyz<1)
print("Selected the indices of those objects that are within 1kpc spherical radius")
d_xyz_lessthan1=d_xyz[ind_d_xyz_lessthan1] #here we stored those distances that are less than 1 in xyz plane
print("Total objects in that radius is",len(d_xyz_lessthan1)) #we found 971 elements that satisfy our condition


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
fig7.savefig("./plots/spherical_cut_1_young1_age7to8inlog10scale.png")



###Now let's find the matching ids in the previous snapshot

ind=np.array(0) #Creates an array ind with an element 0
for i in range(len(id_young1_d_xyz)):
  match=np.where((id_691==id_young1_d_xyz[i])&(id_child_691==id_child_young1_d_xyz[i])) #find each matching id and store its index in match
  print("Matched the id",id_691[match])
  ind=np.append(ind,match) # the index of matching id is appended. But it has one extra element in the beginning
  
print(ind)
ind_tracked_id_691=ind[1:len(ind)] #The extra element in the beginning is removed by this process
tracked_id_691=id_691[ind_tracked_id_691]

print("The number of stars in Snapshot 691 is",len(tracked_id_691))
print("Which is strange as they are not equal. I think 691 could have a little less but not more !!!!!!!!!!!! ") # There are duplicated IDss !!!
print("Total ids in Snapshot 696 was",len(id))
print("Total ids in Snapshot 691 was",len(id_691))


x_691_tracked=x_691[ind_tracked_id_691]
y_691_tracked=y_691[ind_tracked_id_691]
z_691_tracked=z_691[ind_tracked_id_691]


#Now plotting ind 2D the stars in two snapshots
fig8 = plt.figure()
ax1 = fig8.add_subplot(211)
ax1.scatter(y_691_tracked,x_691_tracked,marker=".",s=0.5)
ax1.set_xlabel('y')
ax1.set_ylabel('x')
ax1.set_title('Snapshot 691: young stars age 7 to 8 dec within 1kpc sphere from 0,8,0')

plt.subplots_adjust(hspace=.5)

ax2=fig8.add_subplot(212)
ax2.scatter(y_young1_d_xyz,x_young1_d_xyz,marker=".",s=0.5)
ax2.set_xlabel('y')
ax2.set_ylabel('x')
ax2.minorticks_on()

ax2.set_title('Snapshot 696: young stars age 7 to 8 dec within 1kpc sphere from 0,8,0')
#plt.tight_layout()
fig8.savefig("./plots/young1_age7to8_snapshots691_and_696.png")





#Now plotting figure 08 with y limits and x limits
fig9 = plt.figure()
ax1 = fig9.add_subplot(211)
ax1.scatter(y_691_tracked,x_691_tracked,marker=".",s=0.5)
ax1.set_xlabel('y')
ax1.set_ylabel('x')
ax1.set_ylim(-2,2)
ax1.set_xlim(-10,10)
ax1.set_title('Snapshot 691: young stars age 7 to 8 dec within 1kpc sphere from 0,8,0')

plt.subplots_adjust(hspace=.5)

ax2=fig9.add_subplot(212)
ax2.scatter(y_young1_d_xyz,x_young1_d_xyz,marker=".",s=0.5)
ax2.set_xlabel('y')
ax2.set_ylabel('x')
ax2.minorticks_on()
ax2.set_ylim(-2,2)
ax2.set_xlim(-10,10)
ax2.set_title('Snapshot 696: young stars age 7 to 8 dec within 1kpc sphere from 0,8,0')
#plt.tight_layout()
fig9.savefig("./plots/same_scale_young1_age7to8_snapshots691_and_696.png")