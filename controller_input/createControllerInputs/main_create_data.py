import sys
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
from netCDF4 import date2num,num2date
sys.path.insert(0, '/Users/cconn/Documents/Explore_controller/controller')
import main_controller as CONTROLLER
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

import pandas as pd
times = pd.date_range(start='01/01/2035', end='12/31/2069', freq='M')
f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/TREFHT_Basestates_SSP45.nc")
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
    ax=plt.axes(projection=ccrs.Robinson())
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
    empty_array = np.empty(shape = (len(files), 420, 192, 288))
    for i, file in enumerate(files):
        print(i)
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
    f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/" + VARIABLE + "_" + TYPE + ".nc")
    vals = f[NAME]
    return(vals)

###############################
######CALCAULTE DATA###########
###############################
SAVE_ANOMALIES = False
SAVE_BASESTATES = False
CALC_ENSO = True
CALC_NAO = True

data = {}
data_type = "SSP45" #"LE" #"SSP45"
PARTS = ["BASESTATES", "ENSO", "NAO"]

ENSO_STATES = (np.arange(-3 ,3.5, .5))#[-2,-1.75,-1.5,-1.25,-1,-.75,-.5,-.25,-.01,0,.01,.25,.5,.75,1,1.25,1.5,1.75,2,2.25]
NAO_STATES = [0]#(np.arange(-3, 3.5, .5))#[-1, 0, 1]
BASESTATES = [2035] #np.arange(2035,2070,2)#[2040]

exp_name = "data_P5_SSP45_DELETE"
# exp_name = "data_P3_DELETE"

if SAVE_BASESTATES or SAVE_ANOMALIES:
    VARIABLE_BASESTATE = 'TREFHT'
    folder_path = "/Users/cconn/Documents/CESM245_data/TREFHT/2035-2070/"
    # folder_path = "/Users/cconn/Documents/CESM2-LE/SSP237-2035-2070/PSL/"
    data[VARIABLE_BASESTATE] = Open_Data_Combine(VARIABLE_BASESTATE, folder_path)

    import CalculateBasestate as Cbase
    Cbase.CalculateBasestate(SAVE_ANOMALIES, SAVE_BASESTATES, data[VARIABLE_BASESTATE], VARIABLE_BASESTATE, 
                             lats, lons, times, "_SSP45")  

if CALC_ENSO:
     if data_type == "LE":
        if os.path.isfile("/Users/cconn/Documents/Explore_controller/controller_input/Data/ENSO_anomalies_LE.nc"):
            f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/ENSO_anomalies_LE.nc")
            ENSO_index = f["ENSO_anomalies"]
        else:
            import CalculateENSO as calc_ENSO
            data["SST_anomalies_LE"] = Open_Data_Calculated("SST", "anomalies_LE", "SST_anom") 
            data["SST_anomalies_LE"]["time"] = times
            ENSO_index = np.array(calc_ENSO.Calc_ENSO34(data["SST_anomalies_LE"], False, 5, lats, lons))
     if data_type == "SSP45":
        if os.path.isfile("/Users/cconn/Documents/Explore_controller/controller_input/Data/ENSO_anomalies_SSP45.nc"):
            f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/ENSO_anomalies_SSP45.nc")
            ENSO_index = f["ENSO_anomalies"]
        else:
            import CalculateENSO as calc_ENSO
            data["SST_anomalies_SSP45"] = Open_Data_Calculated("SST", "anomalies_SSP45", "SST_anom") 
            data["SST_anomalies_SSP45"]["time"] = times
            ENSO_index = np.array(calc_ENSO.Calc_ENSO34(data["SST_anomalies_SSP45"], False, 5, lats, lons))
    
if CALC_NAO:
    if data_type == "LE":
        if os.path.isfile("/Users/cconn/Documents/Explore_controller/controller_input/Data/NAO_anomalies_LE.nc"):
            f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/NAO_anomalies_LE.nc")
            NAO_index = f["NAO_anomalies"]
        else:
            import CalculateNAO as calc_NAO
            data["PSL_anomalies_LE"] = Open_Data_Calculated("PSL", "anomalies_LE", "PSL_anom")  
            data["PSL_anomalies_LE"]["time"] = times
            # calc_NAO.Calc_NAO_EOF(data["PSL_anomalies"])
            NAO_index = np.array(calc_NAO.Calc_NAO(data["PSL_anomalies_LE"], lats, lons)) 
    if data_type == "SSP45":
        if os.path.isfile("/Users/cconn/Documents/Explore_controller/controller_input/Data/NAO_anomalies_SSP45.nc"):
            f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/NAO_anomalies_SSP45.nc")
            NAO_index = f["NAO_anomalies"]
        else:
            import CalculateNAO as calc_NAO
            data["PSL_anomalies_SSP45"] = Open_Data_Calculated("PSL", "anomalies_SSP45", "PSL_anom")  
            data["PSL_anomalies_SSP45"]["time"] = times
            # calc_NAO.Calc_NAO_EOF(data["PSL_anomalies"])
            NAO_index = np.array(calc_NAO.Calc_NAO(data["PSL_anomalies_SSP45"], lats, lons)) 

###############################
###########OPEN DATA###########
###############################

CONTINUE = True
ERROR_CHECK = False
if CONTINUE:
    f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/TREFHT_Basestates_SSP45.nc")

    lats = f['lat']
    lons = f['lon']
    basestates = f["basestate"]
    f.close()

    f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/TREFHT_Anomalies_" + data_type + ".nc")
    lats = f['lat']
    lons = f['lon']
    anomalies = f["TREFHT_anom"]
    f.close()
else:
    sys.exit(0)

if ERROR_CHECK:
    print(np.shape(anomalies))
    print(np.shape(basestates))

basestates = Add_Metadata_3(basestates)
BASESTATE_maps = np.empty(shape=(len(BASESTATES), 192, 288))
ENSO_maps = np.empty(shape=(len(ENSO_STATES[:-1]), 192, 288))
NAO_maps = np.empty(shape=(len(NAO_STATES[:-1]), 192, 288))

basestate_annual = basestates.groupby("time.year").mean()
for b, B in enumerate(BASESTATES):
    BASESTATE_maps[b, :, :] = basestate_annual.sel(year = B)

########NAO############
n_NAO_maps = []
if list(NAO_STATES) != [0]:
    for n, na in enumerate(NAO_STATES[:-1]):
        n_NAO_maps_mem = []
        NAO_members = np.empty(shape=(10, 192, 288))
        for member in np.arange(0,10,1):
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
        # print(n_NAO_maps)
        # Plot_Gobal_Map(NAO_maps[n,:,:], "NAO BETWEEN " + str(NAO_STATES[n]) + " and " + str(NAO_STATES[n + 1]), np.linspace(-2.5,2.5,15), tmap)
else:
    NAO_maps = np.empty(shape=(1, 192, 288))  
    NAO_maps.fill(0)

#######ENSO############
n_ENSO_maps = []
if list(ENSO_STATES) != [0]:
    for e, en in enumerate(ENSO_STATES[:-1]):
        n_ENSO_maps_mem = []
        ENSO_members = np.empty(shape=(10, 192, 288))
        ENSO_index["time"] = times
        anomalies["time"] = times
        for member in np.arange(0,10,1):
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
        print(n_ENSO_maps)
        # Plot_Gobal_Map(ENSO_maps[e,:,:], "ENSO BETWEEN " + str(ENSO_STATES[e]) + " and " + str(ENSO_STATES[e + 1]), np.linspace(-3,3,15), tmap)
else:
    ENSO_maps = np.empty(shape=(1, 192, 288))  
    ENSO_maps.fill(0)
  
VOL_maps = np.empty(shape=(1, 192, 288))  
VOL_maps.fill(0)
SAM_maps = np.empty(shape=(1, 192, 288))  
SAM_maps.fill(0)
################
ADD_FILEHEADERS = ['BASESTATE_YEAR', 'ENSO_MIN', 'ENSO_MAX', 'ENSO_n', 'ENSO_q', 'NAO_MIN', 'NAO_MAX', 'NAO_n', 'NAO_q']
   
INJECT = True
PLOT_COMPOSITES = False
    
for b, B in enumerate(BASESTATES):
    for e, EN in enumerate(ENSO_maps):
        for n, NA in enumerate(NAO_maps):
            CONTROLLER_map = BASESTATE_maps[b,:,:] + ENSO_maps[e,:,:] + NAO_maps[n,:,:]
            METADATA = [str(B)]
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
            if INJECT:
                CONTROLLER_map = Add_Initial_dim(CONTROLLER_map)
                # Plot_Gobal_Map(CONTROLLER_map[0,:,:], str(B), np.linspace(220,320,20), 'Reds')
                CONTROLLER.Controller_Injection(CONTROLLER_map, exp_name, ADD_FILEHEADERS, METADATA)
            if PLOT_COMPOSITES:
                import Plot_composits as PC
                print(np.shape(CONTROLLER_map))
                PC.COMPOSITES_SIX(lats, lons, ENSO_maps[e,:,:], NAO_maps[n,:,:], 
                                  VOL_maps[0,:,:], SAM_maps[0,:,:], 
                                  BASESTATE_maps[b,:,:], CONTROLLER_map,
                                  METADATA)


##########EXTRA
# Plot_Gobal_Map(ENSO_maps[e,:,:], "ENSO BETWEEN " + str(ENSO_STATES[e]) + " and " + str(ENSO_STATES[e + 1]), np.linspace(-3,3,15), tmap)

# for b, B in enumerate(BASESTATES):
#     for e, EN in enumerate(ENSO_maps):
#         for n, NA in enumerate(NAO_maps):
#             CONTROLLER_map = BASESTATE_maps[b,:,:] + ENSO_maps[e,:,:] + NAO_maps[n,:,:]
#             METADATA = [str(B)]
#             if list(ENSO_STATES) != [0]:
#                 METADATA.append(str(ENSO_STATES[e]))
#                 METADATA.append(str(ENSO_STATES[e + 1]))
#                 METADATA.append(str(n_ENSO_maps[e]))
#                 METADATA.append("ENSO")
#             else:
#                 METADATA.append(0)
#                 METADATA.append(0)
#                 METADATA.append(0)
#                 METADATA.append("NA")
#             if list(NAO_STATES) != [0]:
#                 METADATA.append(str(NAO_STATES[n]))
#                 METADATA.append(str(NAO_STATES[n + 1]))
#                 METADATA.append(str(n_NAO_maps[n]))
#                 METADATA.append("NAO")
#             else:
#                 METADATA.append(0)
#                 METADATA.append(0)
#                 METADATA.append(0)
#                 METADATA.append("NA")
#             CONTROLLER_map = Add_Initial_dim(CONTROLLER_map)
#             # Plot_Gobal_Map(CONTROLLER_map[0,:,:], str(B), np.linspace(220,320,20), 'Reds')
#             CONTROLLER.Controller_Injection(CONTROLLER_map, exp_name, ADD_FILEHEADERS, METADATA)

