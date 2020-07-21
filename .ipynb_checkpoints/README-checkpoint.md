# stampede_workfolder
#This is my work folder in STAMPEDE

#load_data_using_dicts.py file contains the code to load the x,y,z, age, id and id_child of the star paricles in snapshots. You can change the no. of snapshots to read by changing snapshot_start and snapshot_end value.
#First run this file using:
from load_data_using_dicts import *

#selections_using_dicts.py file contains the code to assign the id and id_child of the star cluster to begin with.
#It also has code to compare the id and id_child with other snapshots and return the matching indices in those snapshots.
#It also calcualtes the center of mass and gives a plot automatically based on the no. of snapshots.
#run the file by running this line after above

exec(open('selections_using_dicts.py').read())exec(open('selections_using_dicts.py').read())


#plots folder contains the plots generated.
