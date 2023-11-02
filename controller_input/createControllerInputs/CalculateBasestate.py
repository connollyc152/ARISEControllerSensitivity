import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.colors
from cartopy.util import add_cyclic_point
import netCDF4 as nc
import pandas as pd
from netCDF4 import date2num,num2date
import gc
gc.enable()

lower = plt.cm.RdBu_r(np.linspace(0,.49, 49))
white = plt.cm.RdBu_r(np.ones(2)*0.5)
upper = plt.cm.RdBu_r(np.linspace(0.51, 1, 49))
colors = np.vstack((lower, white, upper))
tmap = matplotlib.colors.LinearSegmentedColormap.from_list('terrain_map_white', colors)

def Add_Metadata_4(data, lats, lons, times):
    data = xr.DataArray(data,
        dims = ['member','time','lat','lon'],
        coords=dict(
            member = (range(0,len(data[:,0,0,0]))),
            time = (times), 
            lat = (lats), 
            lon = (lons)
            )
        )
    return data

def Add_Metadata_3(data, lats, lons, times):
    data = xr.DataArray(data,
        dims = ['time','lat','lon'],
        coords=dict(
            time = (times), 
            lat = (lats), 
            lon = (lons)
            )
        )
    return data


def Plot_Gobal_Map(plot_data, title, levels, colorbar, lats, lons):
    ax=plt.axes(projection=ccrs.Robinson())
    #levels = np.arange(-7, 9, 2)
    data=plot_data
    data, lonsr = add_cyclic_point(data, coord=lons)
    cs=ax.contourf(lonsr, lats, data,
                transform = ccrs.PlateCarree(),extend='both', cmap = colorbar,
                levels = levels)

    ax.coastlines()
    cbar = plt.colorbar(cs,shrink=0.7,orientation='horizontal',label='Surface Air Temperature (K)')#, pad=5)
    #cbar.set_ticks([-7,-5,-3,-1,0,1,3,5,7])
    cbar.set_label("Surface Air Temperature Anomalies (C)")
    plt.title(title)
    
    plt.show()

    
def CalculateBasestate(SAVE_ANOMALIES, SAVE_BASESTATES, SAVE_DATA, data, VARIABLE,
                       lats, lons, times, EXTRA):
    data = Add_Metadata_4(data, lats, lons, times)
    if SAVE_DATA:
        print("saving data")
        ts = nc.Dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/" + VARIABLE + "_" + EXTRA + ".nc", 'w' , format='NETCDF4')
        ts_members = ts.createDimension('members',len(data[:, 0, 0, 0]))
        ts_time = ts.createDimension('time',len(times))
        ts_lat = ts.createDimension('lat',len(lats))
        ts_lon = ts.createDimension('lon',len(lons))
         
        ts_instd = ts.createVariable(VARIABLE + '_anom','f4',('members','time','lat','lon'))
        ts_time = ts.createVariable('time','f4',('time'))
        ts_lat = ts.createVariable('lat','f4',('lat'))
        ts_lon = ts.createVariable('lon','f4',('lon'))
        
        ts_instd[:,:,:,:] = data
        ts_time[:] = times
        ts_lat[:] = lats
        ts_lon[:] = lons   
    
    print("CALCULATING BASESTATE")
    print(np.shape(data))
    n_members = len(data[:,0,0,0])
        
    if VARIABLE == "SST":
        data = data - 270

    ensemble_mean = np.mean(data, axis = 0) #[420, 192, 288]

    #polyfit of the ensemble_mean (K)
    basestate = np.empty(shape = (420, 192, 288))
    for lat in np.arange(0, len(lats), 1):
        for lon in np.arange(0, len(lons), 1):
            x = np.arange(0, len(ensemble_mean[:,0,0]), 1)
            z = np.polyfit(x,ensemble_mean[:,lat,lon], deg = 3)

            basestate[:,lat,lon] = np.polyval(z,x)
    
    #Remove the climate change trend  
    detrended = np.empty(shape = (len(data[:,0,0,0]), 420, 192, 288))
    for member in range(0, len(data[:,0,0,0])):
        for lat in np.arange(0, len(lats), 1):
            for lon in np.arange(0, len(lons), 1):
                print(member)
                detrended[member, :, lat, lon] = data[member,:,lat,lon] - basestate[:, lat, lon]
    
    del(data)
    #Remove the seasonal cycle
    data_NS = np.empty(shape = (n_members, 420, 192, 288))
    detrended = Add_Metadata_4(detrended, lats, lons, times)
    for member in range(0, n_members):
        print(member)
        climatology = detrended[member, :, :, :].groupby("time.month").mean("time")
        data_NS[member,:,:,:] = detrended[member,:,:,:].groupby("time.month") - climatology
    
    del(detrended)
    data_NS = Add_Metadata_4(data_NS, lats, lons, times)
    basestate = Add_Metadata_3(basestate, lats, lons, times)
    print(data_NS)
    print(basestate)

    if SAVE_BASESTATES:
        # from netCDF4 import date2num,num2date
        # ts_time_units = 'hours since 1800-01-01'
        # dates = date2num(times, ts_time_units, 'noleap')
        
        ts = nc.Dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/" + VARIABLE + "_Basestates" + EXTRA + ".nc", 'w' , format='NETCDF4')
        ts_time = ts.createDimension('time',len(times))
        ts_lat = ts.createDimension('lat',len(lats))
        ts_lon = ts.createDimension('lon',len(lons))
         
        ts_ensemble_mean = ts.createVariable('ensemble_mean','f4',('time','lat','lon'))
        ts_basestate = ts.createVariable('basestate','f4',('time','lat','lon'))
        ts_time = ts.createVariable('time','f4',('time'))
        ts_lat = ts.createVariable('lat','f4',('lat'))
        ts_lon = ts.createVariable('lon','f4',('lon'))
        
        # ts_ensemble_mean[:,:,:] = ensemble_mean
        ts_basestate[:,:,:] = basestate
        ts_time[:] = times
        ts_lat[:] = lats
        ts_lon[:] = lons

        print(basestate)
        print(times)

        ts.close()  
    del(basestate)
    
    if SAVE_ANOMALIES:
        print("saving anomalies")
        ts = nc.Dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/" + VARIABLE + "_anomalies" + EXTRA + ".nc", 'w' , format='NETCDF4')
        ts_members = ts.createDimension('members',len(data_NS[:, 0, 0, 0]))
        ts_time = ts.createDimension('time',len(times))
        ts_lat = ts.createDimension('lat',len(lats))
        ts_lon = ts.createDimension('lon',len(lons))
         
        ts_instd = ts.createVariable(VARIABLE + '_anom','f4',('members','time','lat','lon'))
        ts_time = ts.createVariable('time','f4',('time'))
        ts_lat = ts.createVariable('lat','f4',('lat'))
        ts_lon = ts.createVariable('lon','f4',('lon'))
        
        ts_instd[:,:,:,:] = data_NS
        ts_time[:] = times
        ts_lat[:] = lats
        ts_lon[:] = lons     
    
    


