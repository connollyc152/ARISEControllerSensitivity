import sys
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
from netCDF4 import date2num,num2date
sys.path.insert(0, '/Users/cconn/Documents/Explore_controller/controller')
import main_controller as CONTROLLER

import matplotlib.colors
lower = plt.cm.RdBu_r(np.linspace(0,.49, 49))
white = plt.cm.RdBu_r(np.ones(2)*0.5)
upper = plt.cm.RdBu_r(np.linspace(0.51, 1, 49))
colors = np.vstack((lower, white, upper))
tmap = matplotlib.colors.LinearSegmentedColormap.from_list('terrain_map_white', colors)
sys.path.insert(0, '/Users/cconn/Documents/Explore_controller/controller_input/createControllerInputs')

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

def Add_Metadata_2(data):
    data = xr.DataArray(data,
        dims = ['lat','lon'],
        coords=dict(
            lat = (lats), 
            lon = (lons)
            )
        )
    return data

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

def Add_Metadata_4(data):
    data = xr.DataArray(data,
        dims = ['member','time','lat','lon'],
        coords=dict(
            member = (range(0,10)),
            time = (times), 
            lat = (lats), 
            lon = (lons)
            )
        )
    return data

def Plot_Gobal_Map(plot_data, title, levels, colorbar):
    ax=plt.axes(projection=ccrs.Robinson())
    #levels = np.arange(-7, 9, 2)
    data=plot_data
    data, lonsr = add_cyclic_point(data, coord=lons)
    cs=ax.contourf(lonsr, lats, data,
                transform = ccrs.PlateCarree(),extend='both', cmap = colorbar,
                levels = levels)

    ax.coastlines()
    cbar = plt.colorbar(cs,shrink=0.7,orientation='horizontal',label='Surface Air Temperature (K)', format='%.0f')#, pad=5)
    #cbar.set_ticks([-7,-5,-3,-1,0,1,3,5,7])
    cbar.set_label("Surface Air Temperature Anomalies (C)")
    plt.title(title)
    
    plt.show()

def Open_Data_Combine(VARIABLE):

    files = ["b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.001.cam.h0." + VARIABLE + ".2012501-206912.nc",
             "b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.002.cam.h0." + VARIABLE + ".2012501-206912.nc",
             "b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.003.cam.h0." + VARIABLE + ".2012501-206912.nc",
             "b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.004.cam.h0." + VARIABLE + ".2012501-206912.nc",
             "b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.005.cam.h0." + VARIABLE + ".2012501-206912.nc",
             "b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.006.cam.h0." + VARIABLE + ".2012501-206912.nc",
             "b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.007.cam.h0." + VARIABLE + ".2012501-206912.nc",
             "b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.008.cam.h0." + VARIABLE + ".2012501-206912.nc",
             "b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.009.cam.h0." + VARIABLE + ".2012501-206912.nc",
             "b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.0010.cam.h0." + VARIABLE + ".2012501-206912.nc",]
    
    datapath = "/Users/cconn/Documents/CESM245_data/" + VARIABLE + "_monthly/combined/"
    
    empty_array = np.empty(shape = (len(files), 660, 192, 288))
    for i, file in enumerate(files):
        f = xr.open_dataset(datapath + file)
        
        global lats, lons, times
        times = f["time"][:] #86
        lats = f['lat']
        lons = f['lon']
        vals = f[VARIABLE][:,:,:]
        
        from netCDF4 import date2num,num2date
        ts_time_units = 'hours since 1800-01-01'
        times = num2date(times, ts_time_units, 'noleap')
        
        empty_array[i, :, :, :] = vals
    
    return(empty_array)

def Open_Data_Calculated(VARIABLE, TYPE, NAME):
    f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/" + VARIABLE + "_" + TYPE + ".nc")
    vals = f[NAME]
    return(vals)

###############################
######CALCAULTE DATA###########
###############################
VARIABLE_BASESTATE = 'SST'

SAVE_ANOMALIES = False
SAVE_BASESTATES = False
CALC_ENSO = True

data = {}

if SAVE_BASESTATES or SAVE_ANOMALIES:
    data[VARIABLE_BASESTATE] = Open_Data_Combine(VARIABLE_BASESTATE)
    
    import CalculateBasestate as Cbase
    Cbase.CalculateBasestate(SAVE_ANOMALIES, SAVE_BASESTATES, data[VARIABLE_BASESTATE], VARIABLE_BASESTATE, 
                             lats, lons, times)
    
    
if CALC_ENSO:
    # FILES SST_anomlies must exist 
    data["SST_anomalies"] = Open_Data_Calculated("SST", "anomalies", "SST_anom")
    
    import CalculateENSO as calc_ENSO
    ENSO_index = np.array(calc_ENSO.Calc_ENSO34(data["SST_anomalies"], True, 5))

###############################
###########OPEN DATA###########
###############################

PARTS = ["BASESTATES", "ENSO"]

ENSO_STATES = [-2,-1.75,-1.5,-1.25,-1,-.75,-.5,-.25,0,.25,.5,.75,1,1.25,1.5,1.75,2,2.25]#[-2,-1.75,-1.5,-1.25,-1,-.75,-.5,-.25,0,.25,.5,.75,1,1.25,1.5,1.75,2,2.25]
BASESTATES = [2060, 2035, 2045]#np.arange(2035,2070,1)#[2040]

CONTINUE = True
ERROR_CHECK = False
if CONTINUE:
    f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/TREFHT_Basestates.nc")

    times = f["time"][:660] #86
    lats = f['lat']
    lons = f['lon']
    ts_time_units = 'hours since 1800-01-01'
    times = num2date(times, ts_time_units, 'noleap')
    basestates = f["basestate"]
    f.close()

    f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/TREFHT_Anomalies.nc")

    times = f["time"][:660] #86
    lats = f['lat']
    lons = f['lon']
    ts_time_units = 'hours since 1800-01-01'
    times = num2date(times, ts_time_units, 'noleap')
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

basestate_annual = basestates.groupby("time.year").mean()
for b, B in enumerate(BASESTATES):
    BASESTATE_maps[b, :, :] = basestate_annual.sel(year = B)

if ENSO_STATES != [0]:
    for e, en in enumerate(ENSO_STATES[:-1]):
        ENSO_members = np.empty(shape=(10, 192, 288))
        for member in np.arange(0,10,1):
            if (len(anomalies[member,
                    (ENSO_index[member] < ENSO_STATES[(e + 1)]) & (ENSO_index[member] > ENSO_STATES[e]),0,0])) == 0:
                ENSO_members[member, :, :] = np.nan
         
            else:
                ENSO_members[member, :, :] = np.mean(anomalies[member,
                        (ENSO_index[member] < ENSO_STATES[(e + 1)]) & (ENSO_index[member] > ENSO_STATES[e]),
                        :,:], axis = 0)

        ENSO_maps[e,:,:] = np.nanmean(ENSO_members, axis = 0)
        # Plot_Gobal_Map(ENSO_maps[e,:,:], "ENSO BETWEEN " + str(ENSO_STATES[e]) + " and " + str(ENSO_STATES[e + 1]), None, None)
else:
    ENSO_maps = np.empty(shape=(1, 192, 288))  
    ENSO_maps.fill(0)

ADD_FILEHEADERS = ['BASESTATE_YEAR', 'ENSO_MIN', 'ENSO_MAX']
    
for b, B in enumerate(BASESTATES):
    for e, EN in enumerate(ENSO_maps):
        CONTROLLER_map = BASESTATE_maps[b,:,:] + ENSO_maps[e,:,:]
        #METADATA = [str(B), 0, 0]
        METADATA = [str(B), str(ENSO_STATES[e]), str(ENSO_STATES[e + 1])]
        
        print("INJECTING")
        CONTROLLER_map = Add_Initial_dim(CONTROLLER_map)
        # Plot_Gobal_Map(CONTROLLER_map[0,:,:], str(B), np.linspace(200,320,20), 'Reds')
        CONTROLLER.Controller_Injection(CONTROLLER_map, "EXPLORING_ENSO_STATE", ADD_FILEHEADERS, METADATA)


##########EXTRA
# Plot_Gobal_Map(ENSO_maps[e,:,:], "ENSO BETWEEN " + str(ENSO_STATES[e]) + " and " + str(ENSO_STATES[e + 1]), np.linspace(-3,3,15), tmap)


