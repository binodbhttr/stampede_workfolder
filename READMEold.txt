# stampede_workfolder
#This is my work folder in STAMPEDE

#loadall_track_export.py file contains the code to load the x,y,z, age, id and id_child of the star paricles in snapshots. You can change the no. of snapshots to read by changing snapshot_start and snapshot_end value.
#This file will also track the snapshots and export the tracked data from each snapshot 
#First run this file using:
from loadall_track_export import *

#track_and_export_old.py file contains the code to assign the id and id_child of the star cluster to begin with.
#It also has code to compare the id and id_child with other snapshots and return the matching indices in those snapshots.
#It also calcualtes the center of mass and gives a plot automatically based on the no. of snapshots.
#run the file by running this line after above

exec(open("track_and_export_old.py").read())


#plots folder contains the plots generated.
#data folder contains the data exported by the file load_track_export.py

#You can read the exported data using the file snapshot_analysis_from_tracked_files.py
exec(open("snapshot_analysis_from_tracked_files.py").read())
