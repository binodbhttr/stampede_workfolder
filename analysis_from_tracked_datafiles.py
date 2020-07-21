#from __main__ import *
from sl_utilities import distinct_colours as dc
from sl_utilities import distance_functions
import utilities as ut

#this file is targetted to help avoid running tonns of snapshots each time to analyse the tracked data
# I am working towards loading each data from this file and generating all the kinds of plots here
# It is assumed we know the snapshot_start and snapshot_end value
snapshot_start=671
snapshot_end=696
n=snapshot_end-snapshot_start+1
count=snapshot_start
importdata={}
for i in range(n):
  file_name="export_tracked_data_snapshot"+str(count)
  path="./data/"
  importdata[count]=ut.io.file_hdf5(path+file_name)
  print("\n Loaded data from the snapshot no.",count,"located in filename:",file_name,"to importdata[",count,"]\n#####\n")
  count=count+1
  
################################################
################################################
  
total_subplots=snapshot_end-snapshot_start+1  
print("\n Now we are going to plot these snapshots !!!!!! \n####\nTotal plots we would need is",total_subplots)
#cols=int(total_subplots**0.5)
cols=3
print("Total columns we need in this plot is",cols)
rows=total_subplots//cols
rows=rows+total_subplots%cols
position = range(1,total_subplots + 1)

fig = plt.figure(figsize=(10,20))
snap=snapshot_start
for i in range(total_subplots): # add every single subplot to the figure with a for loop
    print("\n\nx of snapshot",snap,"is",importdata[snap]['x_tracked'])
    print("xcm of snapshot",snap,"is",importdata[snap]['xcm'])
    ax = fig.add_subplot(rows,cols,position[i])
    ax.scatter(np.absolute(importdata[snap]['x_tracked']),np.absolute(importdata[snap]['y_tracked']),marker=".",s=1)
    ax.plot(np.absolute(importdata[snap]['xcm']),np.absolute(importdata[snap]['ycm']),color='red',marker=".",markersize=1)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.minorticks_on()
    ax.set_xlim(np.absolute(importdata[snap]['xmin']),np.absolute(importdata[snap]['xmax']))
    ax.set_ylim(np.absolute(importdata[snap]['ymin']),np.absolute(importdata[snap]['ymax']))
    title="Snapshot " + str(snap)
    ax.set_title(title)
    snap=snap+1
      
plt.tight_layout()
plotname="./plots/from_file_snapshots_"+str(snapshot_start)+"_to_"+str(snapshot_end)
plt.savefig(plotname,dpi=600)
plt.close()
print("Plot generated and saved as filename:",plotname)