import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from cartopy.util import add_cyclic_point
from netCDF4 import date2num,num2date
import netCDF4 as nc

import pandas as pd
times = pd.date_range(start='01/01/2035', end='12/31/2069', freq='M')

def Save_ENSO_Anom(anoms):
    #saves ENSO index
    
    ts = nc.Dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/ENSO_anomalies.nc", 'w' , format='NETCDF4')
    ts_member = ts.createDimension('member',len(anoms[:,0]))
    ts_time = ts.createDimension('time',len(times))
     
    ts_anomalies = ts.createVariable('ENSO_anomalies','f4',('member','time'))
    ts_time = ts.createVariable('time','f4',('time'))
    
    ts_anomalies[:,:] = anoms
    ts_time[:] = times

    ts.close()   

def Plot_ENSO_Map(plot_data, title, levels, colorbar):
    ax = plt.axes(projection=ccrs.Robinson(central_longitude=180))
    ax.coastlines()
    ax.gridlines()
    plot_data.plot(
        ax=ax, transform=ccrs.PlateCarree(), cmap='coolwarm'
    )
    ax.set_extent((120, 300, 10, -10))
    plt.title(title)
    plt.show()

def Plot_Gobal_Map(plot_data, title, levels, colorbar):
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
    cbar.set_label("")
    plt.title(title)
    
    plt.show()
    
def Add_Metadata_4(data, lats, lons):
    data = xr.DataArray(data,
        dims = ['member','time','lat','lon'],
        coords=dict(
            member = (range(0,len(data[:,0,0,0]))),
            time = (times), 
            # lat = (lats.sel(lat=slice(-5, 5))), 
            # lon = (lons.sel(lon=slice(190, 240)))
            lat = (lats), 
            lon = (lons)
            )
        )
    return data


def Add_Metadata_4_cut(data, lats, lons):
    data = xr.DataArray(data,
        dims = ['member','time','lat','lon'],
        coords=dict(
            member = (range(0,len(data[:,0,0,0]))),
            time = (times), 
            lat = (lats.sel(lat=slice(-5, 5))), 
            lon = (lons.sel(lon=slice(190, 240)))
            )
        )
    return data


def Calc_ENSO34(ensemble_array, STANDARDIZE, RUNMEAN, lats, lons):
    
    ensemble_array = Add_Metadata_4(ensemble_array, lats, lons)
    if STANDARDIZE:
        ensemble_array = ensemble_array.groupby("time.month") / ensemble_array.groupby("time.month").std("time")
    
    #Select the region used to caluclate ENSO index
    tos_nino34 = ensemble_array.sel(lat=slice(-5, 5), lon=slice(190, 240))
    tos_nino34 = Add_Metadata_4_cut(tos_nino34, lats, lons)
    
    #Calculated weighted average
    weights = np.cos(np.deg2rad(tos_nino34.lat))
    weights.name = "weights"
    index_nino34 = tos_nino34.weighted(weights).mean(("lat", "lon"))
    
    #complete the running mean
    if RUNMEAN > 0:
        index_nino34 = index_nino34.rolling(time=RUNMEAN, center=True).mean()
    
    Save_ENSO_Anom(index_nino34)

    return(index_nino34)




