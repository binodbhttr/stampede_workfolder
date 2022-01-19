import pickle

s=597

datapath="/home1/07428/binod/stampede_workfolder/fire2/young_stars_data/"

file_name="fire2_young_clusters__snapshot_"+str(s)+".pkl"

with open(datapath+file_name, "rb") as fp:
      import_cluster = pickle.load(fp)

