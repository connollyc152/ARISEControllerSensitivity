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
    # from netCDF4 import date2num,num2date
    # ts_time_units = 'hours since 1800-01-01'
    # dates = date2num(times, ts_time_units, 'noleap')
    
    ts = nc.Dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/ENSO_anomalies.nc", 'w' , format='NETCDF4')
    ts_member = ts.createDimension('member',100)
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
    # print(lats)
    data = xr.DataArray(data,
        dims = ['member','time','lat','lon'],
        coords=dict(
            member = (range(0,100)),
            time = (times), 
            # lat = (lats.sel(lat=slice(-5, 5))), 
            # lon = (lons.sel(lon=slice(190, 240)))
            lat = (lats), 
            lon = (lons)
            )
        )
    return data

def Add_Metadata_4_cut(data, lats, lons):
    # print(lats)
    data = xr.DataArray(data,
        dims = ['member','time','lat','lon'],
        coords=dict(
            member = (range(0,100)),
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
    
    tos_nino34 = ensemble_array.sel(lat=slice(-5, 5), lon=slice(190, 240))
    tos_nino34 = Add_Metadata_4_cut(tos_nino34, lats, lons)
    
    weights = np.cos(np.deg2rad(tos_nino34.lat))
    weights.name = "weights"
    index_nino34 = tos_nino34.weighted(weights).mean(("lat", "lon"))
    
    if RUNMEAN > 0:
        index_nino34 = index_nino34.rolling(time=RUNMEAN, center=True).mean()
    
    Save_ENSO_Anom(index_nino34)

    return(index_nino34)

def Calc_ENSO34_oneSST(SST):
    
    tos_nino34 = SST.sel(lat=slice(-5, 5), lon=slice(190, 240))

    weights = np.cos(np.deg2rad(tos_nino34.lat))
    weights.name = "weights"
    index_nino34 = tos_nino34.weighted(weights).mean(("lat", "lon"))

    return(index_nino34)

def cal_ninoCR(sst, anom_method, index="nino3.4", standardize=False):
    # calculate anomalies
    #anom_method: 1 for ens mean; 0 for single run
    #anomaly
    
    # if anom_method == 0:
    #     sst_a = sst.groupby(“time.month”) - sst.groupby(“time.month”).mean(“time”)
    # elif anom_method ==1:
    #     sst_a = sst - sst.mean(“ens”)
    # #standardize
    # if standardize == True:
    #     sst_a = sst_a.groupby(“time.month”) / sst_a.groupby(“time.month”).std(“time”)
    #pick key region
    if index == "ONI" or index == "nino3.4":
        reg =[-5, 5, 190, 240] # 5N-5S; 170W-120W
    sst_a = sst.sel(lat=slice(reg[0], reg[1]), lon=slice(reg[2],reg[3]))
    weights = np.cos(np.deg2rad(sst_a.lat))
    weights.name = "weights"
    nino = sst_a.weighted(weights).mean(("lat", "lon"))
    # smooth
    if index == "ONI":
        nino_smooth = nino.rolling(time=3, center=True).mean()
    return nino

# f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/SST_anomalies.nc")
# vals = f["SST_anom"]

# vals = np.array(vals)
# vals = Add_Metadata_4(vals)

# Plot_Gobal_Map(vals[1,0,:,:], 'SST', None, 'Reds')

# nino = cal_ninoCR(vals[1,0,:,:], 1, standardize = True)
# print(nino)

# nino = Calc_ENSO34_oneSST(vals[1,0,:,:])
# print(nino)

# sys.exit()
