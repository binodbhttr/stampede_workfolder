from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import os
import pickle
import matplotlib
import matplotlib.colors as colors
from matplotlib import rc #to use Latex math symbols like 'phi'
import astropy
from astropy.io import ascii
import astropy
from astropy.io import ascii
import matplotlib
import pdb
from importlib import reload
import pandas as pd
import utilities as ut
from sl_utilities import distinct_colours as dc

matplotlib.rc('text', usetex=False)

############################################################################
#read in sim files and find relevant particles
############################################################################
#STAMPEDE
simname = 'm12i_res7100_mhdcv'
simtype="fire2"

datapath="./fire2_data_pkl/" #this is teh path where the data of our tracked clusters is

plot_path="./plots/" #creating a path to store the plots only if it does not exist
if not os.path.exists(plot_path):
  os.makedirs(plot_path)

data_start=586
data_end=696

initial_data=list() #list to store data of all the initial data of different cluster families 
later_data=list() #list to store data of all the later data of different cluster families 

for s in range(data_start,data_end,3):
   
    fn_clusters_initial=simtype+"_clusters_"+simname+"_snapshot_"+str(s)+".pkl" 
     
    with open(datapath+fn_clusters_initial, "rb") as input:
        cluster_data_initial= pickle.load(input)

    print("####################### Total clusters present is",len(cluster_data_initial))
    print("\n############## Keys to access the data: \n",cluster_data_initial[1].keys())

    initial_data.append(cluster_data_initial) #appending the dictionaries to the list





mgfe_stdev=[]
n_clusters=[]
star_count=[]

for i in range(len(initial_data)):
    cluster_count=1
    for j in range(len(initial_data[i])):
        n_clusters_temp=len(initial_data[i])
        mgfe_stdev_temp=np.std(initial_data[i][cluster_count]["mgfe"])
        star_count.append(len(initial_data[i][cluster_count]["mgfe"]))
        #print(len(initial_data[i][cluster_count]))
        cluster_count+=1
        mgfe_stdev.append(mgfe_stdev_temp)
        n_clusters.append(n_clusters_temp)
        
mgfe_stdev=np.array(mgfe_stdev)
n_clusters=np.array(n_clusters)
star_count=np.array(star_count)
#df=pd.DataFrame({"mgfe_stdev":mgfe_stdev}) 
#df.to_excel("mgfe_clusters_metallicities.xlsx")


print("Total no. of clusters present",len(star_count))
print("Total no. of stds collected",len(mgfe_stdev))
df1=pd.DataFrame(star_count,columns=["star_count"],index=None)
df2= pd.value_counts(df1.star_count).to_frame().reset_index()
df2.columns = ['star_count','n']
df2.to_excel("mgfe_star_count_clusters.xlsx")
a=np.array(df2.star_count)
a.sort()
print(a)

mean_stdevs=[]
for i in a:
    keep=np.where(star_count==i)
    mean_stdevs.append(np.mean(mgfe_stdev[keep]))
    
print(mean_stdevs)




binsize = 1
DistanceBin = ut.binning.DistanceBinClass([0,70], binsize, log_scale=False, dimension_number=1)
ArrStats = DistanceBin.get_statistics_profile(star_count, mgfe_stdev) #bin as a function of radius, calculate statistics for feh in those bins, weighted by mass         

#colors = dc.get_distinct(5)
plt.plot(ArrStats['min'],ArrStats['average'],'.-')

plt.ylabel('<[MG/Fe]> (dex)')
plt.xlabel('N_star')
plt.savefig("trail.png")

'''

fig3=plt.figure()
fig3.set_size_inches(8,5)
ax2=fig3.add_subplot(111)

plot_name="test_mgfe_scatter"+simname+"_all_clusters"

#ax2.set_title("Clusters from all 37 Snapshots (586,589,592,...694)")
ax2.scatter(star_count,mgfe_stdev,s=7,c="grey",alpha=0.7,label=r"${\sigma}$ each cluster")
ax2.set_xlabel(r"${\mathrm{N_{star}}}$ in cluster",fontsize=12)
ax2.set_ylabel(r"${\sigma}$ [Mg/Fe] (dex)",fontsize=12)  
ax2.scatter(a,mean_stdevs,c="red",s=10,label=r"< ${\sigma}$ > each ${\mathrm{N_{star}}}$")
global_mean=np.mean(mean_stdevs)
ax2.axhline(y=global_mean,c="blue",label=r"global average of ${\sigma}$ = %.3f"%global_mean)
ax2.text(41,0.033,"open clusters in m12i_res7100_mhdcv")
ax2.text(44,0.031,"step 586 - 694 with 3 Myr spacing")
ax2.text(54,0.029,"nmin = 5, b=10 pc",ha="left")
#ax2.text(40,0.055,r"global average of ${\sigma}$  = %.3f"%global_mean,fontsize=12,c="blue")
ax2.set(xlim=[5,70])
ax2.minorticks_on()
ax2.legend(bbox_to_anchor=(0.15,0.89), loc='center left')
plt.tight_layout()
fig3.savefig(plot_path+plot_name+".png",bbox_inches='tight',dpi=100)  
'''  