#This file is designed to export gas data of the snapshots into a file that can be used to create plots with gas in the background
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import os
import pickle
import matplotlib
import matplotlib.colors as colors
import gizmo_analysis as gizmo
import astropy
from astropy.io import ascii
#############################################################################
#constants
#############################################################################
MsunToGm = 1.99e33
KpcToCm = 3.086e21
mp = 1.67e-24
#bin_edge = 10.
bin_edge = 30.
bins = np.arange(-25,25,0.1)

############################################################################
#read in sim files and find relevant particles
############################################################################
simname = 'm12i_res7100_mhdcv'
simtype="fire2"  #this is the simtype eg. fire2, fire3, sf-fire3, sf-fire3-alpha01, sf-fire3-alpha03, sf-fire3-alpha05 
simdir = '/scratch/projects/xsede/GalaxiesOnFIRE/cr_suite/m12i_res7100/mhdcv/1Myr/fire2/'
gas_datapath="./"+simtype+"_gas_data_pkl/"  #path to store the gas data
if not os.path.exists(gas_datapath):
  os.makedirs(gas_datapath)
                
snap=650 #this is the snapshot at which the clusters were taken from using the fof algorithm I am using this to extract information from the snapshot where the clusters were first seen
snapshot_start=650  #snapshot to begin creating the figure
snapshot_end=650    #snapshot to stop at

snapshot_list=np.arange(snapshot_start,snapshot_end+1) #create a list of snapshot numbers to plot to plot eg. [596,597, ...]
time=snapshot_list-snapshot_start #time starts from zero here where t=0 is at the snapshot where we start from eg. [0,1,2,3,....]

#reading in snapshot_times.txt file to get snapshot numbers and times
# columns are: snapshot scale-factor redshift time[Gyr] time_width[Myr]
snapshot_times = simdir + 'snapshot_times.txt'
snaptime_data = astropy.io.ascii.read(snapshot_times, guess=False, comment="#")
snaptime_data = np.genfromtxt(snapshot_times, usecols=(0,3), skip_header=4, dtype=float) #the first and fourth columns are the only ones we need 
snaps = np.array(snaptime_data[:,0]) #col1 = first column saved from text file #This is a collection of all snapshot nos.
times = np.array(snaptime_data[:,1]) #col4 = second column saved #This is a collection of times equivalent to those snapshot nos.
#######################################
#######################################


tracked_gas_all_snaps={}
for i in range(len(time)):                              
    snapnumber=time[i]+snapshot_start      #this is the true snapshot number each time for eg, 596, 597 .. and son
    snaptime = times[np.where(snaps == snapnumber)][0] #time of snapshot in Gyr
    tracked_gas={} 
    part = gizmo.io.Read.read_snapshots(['all'],'snapshot_index', time[i]+snapshot_start, simulation_name=simname, simulation_directory=simdir, assign_hosts_rotation=True, assign_hosts=True)  
    t = np.max(part['star'].prop('form.time'))  
    
    rGas = part['gas'].prop('host.distance.principal.cylindrical')[:,0]
    zGas = part['gas'].prop('host.distance.principal.cylindrical')[:,1]
    
    xGas = part['gas'].prop('host.distance.principal.cartesian')[:,0]
    yGas = part['gas'].prop('host.distance.principal.cartesian')[:,1]
    zGas = part['gas'].prop('host.distance.principal.cartesian')[:,2]
    
    vxGas = part['gas'].prop('host.velocity.principal.cartesian')[:,0]
    vyGas = part['gas'].prop('host.velocity.principal.cartesian')[:,1]
    vzGas = part['gas'].prop('host.velocity.principal.cartesian')[:,2]
    
    mGas = part['gas']['mass']
    rhoGas = part['gas']['density']
    tGas = part['gas']['temperature']
    idGas = part['gas']['id']
    
    i_gas = np.where((rGas <= bin_edge) & (np.fabs(zGas) <= 1.5) & (part['gas']['density']*((MsunToGm/KpcToCm**3)/mp) >= 10.) & (tGas <= 1e4))
    
    x = xGas[i_gas]
    y = yGas[i_gas]
    z = zGas[i_gas]
    vx = vxGas[i_gas]
    vy = vyGas[i_gas]
    vz = vzGas[i_gas]
    m = mGas[i_gas]
    rho = part['gas'].prop('number.density')[i_gas]
    id = part['gas']['id'][i_gas]
    
    ###########################################################################
    #gas image (2d histogram)
    ###########################################################################
    #cold (< 10^4 K) gas in the midplane (|z| <= 1.5 kpc within bin_edge
    v =  np.where((rGas <= bin_edge) & (np.fabs(zGas) <= 1.5) & (tGas <= 1e4))
    face, xh, yh = np.histogram2d(part['gas'].prop('host.distance.principal.cartesian')[v,1][0],part['gas'].prop('host.distance.principal.cartesian')[v,0][0],bins=[bins,bins], weights=part['gas']['mass'][v])

    ###########################################################################
    tracked_gas={"snaptime":snaptime,"v":v,"face":face,"xh":xh,"yh":yh,"xGas":x,"yGas":y,"zGas":z,"vxGas":vx,"vyGas":vy,"vzGas":vz,"mGas":m}
    file_name=simtype+"_gas_data"+str(snapnumber)+".pkl"
    with open(gas_datapath+file_name, 'wb') as output:
      pickle.dump(tracked_gas, output)
    print("\n Stored the gas data for background plot from the snapshot no.",snap,"to filename:",file_name,"\n#####\n")
    '''
    This is how to read the data:
    file_name="gas_data_snapshot_"+str(snapnumber)+".pkl"
    with open(gas_datapath+file_name, "rb") as input:
      import_gasdata = pickle.load(input)
    print(import_gasdata["v"])
    '''
    #tracked_gas_all_snaps.update({snapnumber:tracked_gas})
###########################################
#storing all the gas data from all snapshots as a single file
'''
with open(gas_datapath+simtype+"_gas_data_all_snapshots.pkl", 'wb') as output:
  pickle.dump(tracked_gas_all_snaps, output) #
  
  #read it as tracked_gas_all_snaps[snap]["v"]
  '''