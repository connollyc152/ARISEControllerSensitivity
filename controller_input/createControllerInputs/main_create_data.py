import sys
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
sys.path.insert(0, '/Users/cconn/Documents/Explore_controller/controller')
import main_controller as CONTROLLER
sys.path.insert(0, '/Users/cconn/Documents/Explore_controller/controller_output/')
import plotFunctions as pF
import os.path
import gc
gc.enable()

import matplotlib.colors
lower = plt.cm.RdBu_r(np.linspace(0,.49, 49))
white = plt.cm.RdBu_r(np.ones(2)*0.5)
upper = plt.cm.RdBu_r(np.linspace(0.51, 1, 49))
colors = np.vstack((lower, white, upper))
tmap = matplotlib.colors.LinearSegmentedColormap.from_list('terrain_map_white', colors)
sys.path.insert(0, '/Users/cconn/Documents/Explore_controller/controller_input/createControllerInputs')

timesteps_n = 420
import pandas as pd
times = pd.date_range(start='01/01/2035', end='12/31/2069', freq='M')

f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/monthly/TREFHT_Basestates_SSP45.nc")
lats = f['lat']
lons = f['lon']
f.close()

def Add_Initial_dim(data):
    array = np.empty(shape = (1,192,288))
    array[0,:,:] = data
    array = xr.DataArray(array,
        dims = ['length','lat','lon'],
        coords=dict(
            length = [1],
            lat = (lats), 
            lon = (lons)
            )
        )
    return array

def Add_Metadata_3(data):
    data = xr.DataArray(data,
        dims = ['time','lat','lon'],
        coords=dict(
            time = (times), 
            lat = (lats), 
            lon = (lons)
            )
        )
    return data

def Plot_Gobal_Map(plot_data, title, levels, colorbar):
    ax=plt.axes(projection= ccrs.PlateCarree())
    data=plot_data
    data, lonsr = add_cyclic_point(data, coord=lons)
    cs=ax.contourf(lonsr, lats, data,
                transform = ccrs.PlateCarree(),extend='both', cmap = colorbar,
                levels = levels)

    ax.coastlines()
    cbar = plt.colorbar(cs,shrink=0.7,orientation='horizontal',label='Surface Air Temperature (K)', format='%.0f')#, pad=5)
    cbar.set_label("Surface Air Temperature Anomalies (C)")
    plt.title(title)
    
    plt.show()

def Open_Data_Combine(VARIABLE, file_folder_path):

    files = []
    for filename in os.listdir(file_folder_path):
        files.append(os.path.join(filename))
    files = list(np.sort(files))
    
    if ".DS_Store" in files:
        files.remove(".DS_Store")    

    # files = files[:10]
    empty_array = np.empty(shape = (len(files), timesteps_n, 192, 288))
    for i, file in enumerate(files):
        f = xr.open_dataset(file_folder_path + file)
        
        global lats, lons, times
        times = f["time"][:] #86
        lats = f['lat']
        lons = f['lon']
        vals = f[VARIABLE][:,:,:]
        
        empty_array[i, :, :, :] = vals
        
        f.close()
    
    return(empty_array)

def Open_Data_Calculated(VARIABLE, TYPE, NAME):
    f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/monthly/" + VARIABLE + "_" + TYPE + ".nc")
    vals = f[NAME]
    return(vals)

###############################
######CALCAULTE DATA###########
###############################
SAVE_DATA = False
SAVE_ANOMALIES = False
SAVE_BASESTATES = False
CALC_ENSO = True
CALC_NAO = True
CALC_SAM = True
CALC_VOLC = True

#Creates empty dictionary that will be filled with opened data
data = {}
#open either ARISE; "SSP45" or the large ensemble: "LE"
data_type = "SSP45" #"LE" #"SSP45"

#Indicate which components to include in the controller input. [0]: not used as
#an input. Code will calculate average between a range. Volcano can only either
#be on or off. Any year can be used in BASESTATES
ENSO_STATES = [0]#[-1.2, -1, 1, 1.2] #np.arange(-2, 2.2, .2)#(np.arange(-2.4, 2.7, .3))#[-2,-1.75,-1.5,-1.25,-1,-.75,-.5,-.25,-.01,0,.01,.25,.5,.75,1,1.25,1.5,1.75,2,2.25]
NAO_STATES = [0]#np.arange(-2, 2.2, .2)#[-1.2, -1, 1, 1.2]#np.arange(-2, 2.2, .2)#(np.arange(-2.4, 2.7, .3))#[-1, 0, 1]
SAM_STATES = [0]#np.arange(-2, 2.2, .2)(np.arange(-2.4, 2.7, .3))
VOLC_STATES = False
BASESTATES = [2035,2045] #np.arange(2035,2070,2)#[2040]

#Files output name
exp_name = "monthly_LE_T_AddingTO"

#opens and calculates basestates
if SAVE_BASESTATES or SAVE_ANOMALIES or SAVE_DATA:
    VARIABLE_BASESTATE = 'PSL'
    # folder_path = "/Users/cconn/Documents/CESM245_data/PSL/2035-2070_monthly/"
    folder_path = "/Users/cconn/Documents/CESM2-LE/SSP237-2035-2070_monthly/PSL/"
    data[VARIABLE_BASESTATE] = Open_Data_Combine(VARIABLE_BASESTATE, folder_path)


    import CalculateBasestate as Cbase
    Cbase.CalculateBasestate(SAVE_ANOMALIES, SAVE_BASESTATES, SAVE_DATA, data[VARIABLE_BASESTATE], VARIABLE_BASESTATE, 
                             lats, lons, times, ("_" + data_type)) 

#if needed, opens temperature anomalies associated with Pinatubo
if CALC_VOLC:
    f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/monthly/VOLC_anomalies.nc")
    VOLC_array = f["__xarray_dataarray_variable__"]
    
#opens or calcualtes the correct ENSO index time series
if CALC_ENSO:
     if data_type == "LE":
        if os.path.isfile("/Users/cconn/Documents/Explore_controller/controller_input/Data/monthly/ENSO_anomalies_LE.nc"):
            print("Found ENSO")
            f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/monthly/ENSO_anomalies_LE.nc")
            ENSO_index = f["ENSO_anomalies"]
        else:
            import CalculateENSO as calc_ENSO
            data["SST_anomalies_LE"] = Open_Data_Calculated("SST", "anomalies_LE", "SST_anom") 
            data["SST_anomalies_LE"]["time"] = times
            ENSO_index = np.array(calc_ENSO.Calc_ENSO34(data["SST_anomalies_LE"], False, 5, lats, lons))
     if data_type == "SSP45":
        if os.path.isfile("/Users/cconn/Documents/Explore_controller/controller_input/Data/monthly/ENSO_anomalies_SSP45.nc"):
            f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/monthly/ENSO_anomalies_SSP45.nc")
            ENSO_index = f["ENSO_anomalies"]
        else:
            import CalculateENSO as calc_ENSO
            data["SST_anomalies_SSP45"] = Open_Data_Calculated("SST", "anomalies_SSP45", "SST_anom") 
            data["SST_anomalies_SSP45"]["time"] = times
            ENSO_index = np.array(calc_ENSO.Calc_ENSO34(data["SST_anomalies_SSP45"], False, 5, lats, lons))
    
#opens or calcualtes the correct NAO index time series
if CALC_NAO:
    if data_type == "LE":
        if os.path.isfile("/Users/cconn/Documents/Explore_controller/controller_input/Data/monthly/NAO_anomalies_LE.nc"):
            print("Found NAO")
            f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/monthly/NAO_anomalies_LE.nc")
            NAO_index = f["NAO_anomalies"]
        else:
            import CalculateNAO as calc_NAO
            data["PSL_anomalies_LE"] = Open_Data_Calculated("PSL", "anomalies_LE", "PSL_anom")  
            data["PSL_anomalies_LE"]["time"] = times
            # calc_NAO.Calc_NAO_EOF(data["PSL_anomalies"])
            NAO_index = np.array(calc_NAO.Calc_NAO(data["PSL_anomalies_LE"], lats, lons)) 
    if data_type == "SSP45":
        if os.path.isfile("/Users/cconn/Documents/Explore_controller/controller_input/Data/monthly/NAO_anomalies_SSP45.nc"):
            f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/monthly/NAO_anomalies_SSP45.nc")
            NAO_index = f["NAO_anomalies"]
        else:
            import CalculateNAO as calc_NAO
            data["PSL_anomalies_SSP45"] = Open_Data_Calculated("PSL", "anomalies_SSP45", "PSL_anom")  
            data["PSL_anomalies_SSP45"]["time"] = times
            # calc_NAO.Calc_NAO_EOF(data["PSL_anomalies"])
            NAO_index = np.array(calc_NAO.Calc_NAO(data["PSL_anomalies_SSP45"], lats, lons)) 

#opens or calcualtes the correct SAM index time series
if CALC_SAM:
    if data_type == "LE":
        if os.path.isfile("/Users/cconn/Documents/Explore_controller/controller_input/Data/monthly/SAM_anomalies_LE.nc"):
            print("Found SAM")
            f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/monthly/SAM_anomalies_LE.nc")
            SAM_index = f["SAM_anomalies"]
        else:
            import CalculateSAM as calc_SAM
            data["PSL_LE"] = Open_Data_Calculated("PSL", "LE", "PSL_anom")  
            data["PSL_LE"]["time"] = times
            SAM_index = np.array(calc_SAM.Calc_SAM(data["PSL_LE"], lats, lons)) 
    if data_type == "SSP45":
        if os.path.isfile("/Users/cconn/Documents/Explore_controller/controller_input/Data/monthly/SAM_anomalies_SSP45.nc"):
            f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/monthly/SAM_anomalies_SSP45.nc")
            SAM_index = f["SAM_anomalies"]
        else:
            import CalculateSAM as calc_SAM
            data["PSL_SSP45"] = Open_Data_Calculated("PSL", "SSP45", "PSL_anom")  
            data["PSL_SSP45"]["time"] = times
            SAM_index = np.array(calc_SAM.Calc_SAM(data["PSL_SSP45"], lats, lons)) 


###############################
###########OPEN DATA###########
###############################
CONTINUE = True
ERROR_CHECK = False

#opens basestates and anomalies
if CONTINUE:
    f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/monthly/TREFHT_Basestates_SSP45.nc")

    lats = f['lat']
    lons = f['lon']
    basestates = f["basestate"]
    f.close()

    f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/monthly/TREFHT_Anomalies_" + data_type + ".nc")
    lats = f['lat']
    lons = f['lon']
    anomalies = f["TREFHT_anom"]
    f.close()
else:
    sys.exit(0)

if ERROR_CHECK:
    print(np.shape(anomalies))
    print(np.shape(basestates))

#creates empty arrays that will be filled will maps associated with different 
#basestates and internal variability modes. 
basestates = Add_Metadata_3(basestates)
BASESTATE_maps = np.empty(shape=(len(BASESTATES), 192, 288))
ENSO_maps = np.empty(shape=(len(ENSO_STATES[:-1]), 192, 288))
NAO_maps = np.empty(shape=(len(NAO_STATES[:-1]), 192, 288))
SAM_maps = np.empty(shape=(len(SAM_STATES[:-1]), 192, 288))

#Creates the one volcano map
VOLC_maps = np.empty(shape=(1, 192, 288))  
if VOLC_STATES:
    VOLC_maps[0,:,:] = VOLC_array
else:
    VOLC_maps.fill(0)
    
basestate_annual = basestates.groupby("time.year").mean()
for b, B in enumerate(BASESTATES):
    BASESTATE_maps[b, :, :] = basestate_annual.sel(year = B)

########SAM############
#Determines the temperature anomaly maps that fall within a specific range of
#SAM indexes and calculates the average
n_SAM_maps = []
SAM_index["time"] = times
anomalies["time"] = times
if list(SAM_STATES) != [0]:
    print("Compositing SAM")
    for n, na in enumerate(SAM_STATES[:-1]):
        n_SAM_maps_mem = []
        SAM_members = np.empty(shape=(len(SAM_index[:,0]), 192, 288))
        for member in np.arange(0,len(SAM_index[:,0]),1):
            if (len(anomalies[member,
                    (SAM_index[member] < SAM_STATES[(n + 1)]) & (SAM_index[member] > SAM_STATES[n]),0,0])) == 0:
                SAM_members[member, :, :] = np.nan
            else:
                SAM_members[member, :, :] = np.mean(anomalies[member,
                        (SAM_index[member] < SAM_STATES[(n + 1)]) & (SAM_index[member] > SAM_STATES[n]),
                        :,:], axis = 0)
                n_SAM_maps_mem.append(len(anomalies[member,
                        (SAM_index[member] < SAM_STATES[(n + 1)]) & (SAM_index[member] > SAM_STATES[n]),
                        0,0]))

        n_SAM_maps.append(sum(n_SAM_maps_mem))
        SAM_maps[n,:,:] = np.nanmean(SAM_members, axis = 0)

else:
    SAM_maps = np.empty(shape=(1, 192, 288))  
    SAM_maps.fill(0)
########NAO############
#Determines the temperature anomaly maps that fall within a specific range of
#NAO indexes and calculates the average
n_NAO_maps = []
NAO_index["time"] = times
if list(NAO_STATES) != [0]:
    print("Compositing NAO")
    for n, na in enumerate(NAO_STATES[:-1]):
        n_NAO_maps_mem = []
        NAO_members = np.empty(shape=(len(NAO_index[:,0]), 192, 288))
        for member in np.arange(0,len(NAO_index[:,0]),1):
            if (len(anomalies[member,
                    (NAO_index[member] < NAO_STATES[(n + 1)]) & (NAO_index[member] > NAO_STATES[n]),0,0])) == 0:
                NAO_members[member, :, :] = np.nan
         
            else:
                NAO_members[member, :, :] = np.mean(anomalies[member,
                        (NAO_index[member] < NAO_STATES[(n + 1)]) & (NAO_index[member] > NAO_STATES[n]),
                        :,:], axis = 0)
                n_NAO_maps_mem.append(len(anomalies[member,
                        (NAO_index[member] < NAO_STATES[(n + 1)]) & (NAO_index[member] > NAO_STATES[n]),
                        0,0]))

        n_NAO_maps.append(sum(n_NAO_maps_mem))
        NAO_maps[n,:,:] = np.nanmean(NAO_members, axis = 0)
        # Plot_Gobal_Map(NAO_maps[n,:,:], "NAO BETWEEN " + str(NAO_STATES[n]) + " and " + str(NAO_STATES[n + 1]), np.linspace(-2.5,2.5,15), tmap)
else:
    NAO_maps = np.empty(shape=(1, 192, 288))  
    NAO_maps.fill(0)

#######ENSO############
#Determines the temperature anomaly maps that fall within a specific range of
#ENSO indexes and calculates the average
n_ENSO_maps = []
ENSO_index["time"] = times
if list(ENSO_STATES) != [0]:
    print("Compositing ENSO")
    for e, en in enumerate(ENSO_STATES[:-1]):
        n_ENSO_maps_mem = []
        ENSO_members = np.empty(shape=(len(ENSO_index[:,0]), 192, 288))
        ENSO_index["time"] = times
        anomalies["time"] = times
        for member in np.arange(0,len(ENSO_index[:,0]),1):
            if (len(anomalies[member, (ENSO_index[member] < ENSO_STATES[(e + 1)]) & (ENSO_index[member] > ENSO_STATES[e]),0,0])) == 0:
                ENSO_members[member, :, :] = np.nan
         
            else:
                ENSO_members[member, :, :] = np.mean(anomalies[member,
                        (ENSO_index[member] < ENSO_STATES[(e + 1)]) & (ENSO_index[member] > ENSO_STATES[e]),
                        :,:], axis = 0)
                n_ENSO_maps_mem.append(len(anomalies[member,
                        (ENSO_index[member] < ENSO_STATES[(e + 1)]) & (ENSO_index[member] > ENSO_STATES[e]),
                        0,0]))

        n_ENSO_maps.append(sum(n_ENSO_maps_mem))
        ENSO_maps[e,:,:] = np.nanmean(ENSO_members, axis = 0)
        # Plot_Gobal_Map(ENSO_maps[e,:,:], "ENSO BETWEEN " + str(ENSO_STATES[e]) + " and " + str(ENSO_STATES[e + 1]), np.linspace(-3,3,15), tmap)
else:
    ENSO_maps = np.empty(shape=(1, 192, 288))  
    ENSO_maps.fill(0)

################
#Files headers are used in the output text file.
ADD_FILEHEADERS = ['BASESTATE_YEAR', 'ENSO_MIN', 'ENSO_MAX', 'ENSO_n', 'ENSO_q', 
                   'NAO_MIN', 'NAO_MAX', 'NAO_n', 'NAO_q', 'SAM_MIN', 'SAM_MAX', 'SAM_n', 'SAM_q', 'VOLC_q']

#Inject determines whether the maps actually pass through the controller and 
#a text file is output. PLOT_COMPOSITES determines if the six paneled plot
#is created. This is pretty slow but does plot all components.    
INJECT = True
PLOT_COMPOSITES = False
    
for b, B in enumerate(BASESTATES):
    for e, EN in enumerate(ENSO_maps):
        for n, NA in enumerate(NAO_maps):
            for s, NA in enumerate(SAM_maps):
                CONTROLLER_map = BASESTATE_maps[b,:,:] + ENSO_maps[e,:,:] + NAO_maps[n,:,:] + SAM_maps[s,:,:] + VOLC_maps[0,:,:]
                METADATA = [str(B)]
                print(str(B) + ", " + str(e))
                if list(ENSO_STATES) != [0]:
                    METADATA.append(str("{:.2f}".format(ENSO_STATES[e])))
                    METADATA.append(str("{:.2f}".format(ENSO_STATES[e + 1])))
                    METADATA.append(str(n_ENSO_maps[e]))
                    METADATA.append("ENSO")
                else:
                    METADATA.append(0)
                    METADATA.append(0)
                    METADATA.append(0)
                    METADATA.append("NA")
                if list(NAO_STATES) != [0]:
                    METADATA.append(str("{:.2f}".format(NAO_STATES[n])))
                    METADATA.append(str("{:.2f}".format(NAO_STATES[n + 1])))
                    METADATA.append(str(n_NAO_maps[n]))
                    METADATA.append("NAO")
                else:
                    METADATA.append(0)
                    METADATA.append(0)
                    METADATA.append(0)
                    METADATA.append("NA")
                if list(SAM_STATES) != [0]:
                    METADATA.append(str("{:.2f}".format(SAM_STATES[s])))
                    METADATA.append(str("{:.2f}".format(SAM_STATES[s + 1])))
                    METADATA.append(str(n_SAM_maps[s]))
                    METADATA.append("SAM")
                else:
                    METADATA.append(0)
                    METADATA.append(0)
                    METADATA.append(0)
                    METADATA.append("NA")
                    
                if VOLC_STATES:
                    METADATA.append("VOLC")
                else:
                    METADATA.append("NA")
                    
                    
                if INJECT:
                    CONTROLLER_map = Add_Initial_dim(CONTROLLER_map)
                    # Plot_Gobal_Map(CONTROLLER_map[0,:,:], str(B), np.linspace(220,320,20), 'Reds')
                    print(np.shape(CONTROLLER_map))
                    print(type(CONTROLLER_map))
                    CONTROLLER.Controller_Injection(CONTROLLER_map, exp_name, ADD_FILEHEADERS, METADATA)
                if PLOT_COMPOSITES:
                    import Plot_composits as PC
                    PC.COMPOSITES_SIX(lats, lons, ENSO_maps[e,:,:], NAO_maps[n,:,:], 
                                      VOLC_maps[0,:,:], SAM_maps[s,:,:], 
                                      BASESTATE_maps[b,:,:], CONTROLLER_map,
                                      METADATA)

