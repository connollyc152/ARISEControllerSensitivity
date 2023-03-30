import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.colors
from cartopy.util import add_cyclic_point
import netCDF4 as nc
import pandas as pd
from netCDF4 import date2num,num2date

lower = plt.cm.RdBu_r(np.linspace(0,.49, 49))
white = plt.cm.RdBu_r(np.ones(2)*0.5)
upper = plt.cm.RdBu_r(np.linspace(0.51, 1, 49))
colors = np.vstack((lower, white, upper))
tmap = matplotlib.colors.LinearSegmentedColormap.from_list('terrain_map_white', colors)

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
    cbar = plt.colorbar(cs,shrink=0.7,orientation='horizontal',label='Surface Air Temperature (K)')#, pad=5)
    #cbar.set_ticks([-7,-5,-3,-1,0,1,3,5,7])
    cbar.set_label("Surface Air Temperature Anomalies (C)")
    plt.title(title)
    
    plt.show()
    
SAVEDATA = False
SAVEDATA_BASESTATES = False 

VARIABLE = "SST"

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

ensemble_array = np.empty(shape = (len(files), 660, 192, 288))

for i, file in enumerate(files):
    f = xr.open_dataset(datapath + file)

    times = f["time"][:] #86
    lats = f['lat']
    lons = f['lon']
    vals = f[VARIABLE][:,:,:]
    
    # ts_time_units = 'hours since 1800-01-01'
    # times = date2num(times, ts_time_units, 'noleap')
    
    if VARIABLE == "SST":
        vals = vals - 270
    
    vals_mean = np.mean(vals, axis = 0)
    
    vals_anom = vals - vals_mean
    
    for lat in np.arange(0, len(lats), 1):
        for lon in np.arange(0, len(lons), 1):
                data = vals_anom[:,lat,lon]
            
                Z = np.fft.fft(data)
                Yfft = Z/np.size(data)

                # combine symmetric parts of the FFT and plot the power spectrum as a function of frequency
                freq = np.arange(0,np.size(data)/2)/float(np.size(data))

                # the factor of 2 in front is needed or the sum won't equal the total variance of X
                Ck2 = 2.*np.abs(Yfft[0:int(np.size(data)/2)+1])**2 
                
                Z[50:60] = 0.0
                Z[-60:-50] = 0.0
            
                X_filtered = np.real(np.fft.ifft(Z))    
            
                ensemble_array[i,:,lat,lon] = X_filtered + float(vals_mean[lat,lon]) 
                #ensemble_array data without seasonal cycle (K)

ensemble_mean = np.mean(ensemble_array, axis = 0)
#ensemble_mean of data with no seasonal cycles (K)

basestate = np.empty(shape = (660, 192, 288))
#polyfit of the ensemble_mean (K)

for lat in np.arange(0, len(lats), 1):
    for lon in np.arange(0, len(lons), 1):
        x = np.arange(0, len(ensemble_mean[:,0,0]), 1)
        z = np.polyfit(x,ensemble_mean[:,lat,lon], deg = 4)

        basestate[:,lat,lon] = np.polyval(z,x)
        
    plt.plot(times, ensemble_mean[:,lat,15], label = 'ensemble mean')
    plt.plot(times, basestate[:,lat,15], label = 'basestate')
    plt.title(str(lats[lat]))
    plt.show()
    
avg_wrong = np.mean(np.mean(basestate, axis = 1), axis = 1)  
plt.plot(avg_wrong)
plt.show()

ensemble_mean = np.mean(ensemble_array, axis = 0)

ensemble_anomalies = ensemble_array - ensemble_mean

ensemble_anomalies = Add_Metadata_4(ensemble_anomalies)

Plot_Gobal_Map(ensemble_mean[0,:,:], 'Ensemble mean 01-01-2035', None, 'Reds')
Plot_Gobal_Map(ensemble_array[0,0,:,:], 'Ensemble member 1 temperature 01-01-2035', None, 'Reds')
Plot_Gobal_Map(ensemble_anomalies[0,0,:,:], 'Ensemble member 1 01-01-2035 anomaly', np.arange(-7, 9, 2), tmap)
Plot_Gobal_Map(ensemble_anomalies[5,0,:,:], 'Ensemble member 5 01-01-2035 anomaly', np.arange(-7, 9, 2), tmap)

Plot_Gobal_Map(basestate[0,:,:], 'basestate 01-01-2035', None, 'Reds')

Plot_Gobal_Map(ensemble_anomalies[0,0,:,:], 'Saved data', None, 'Reds')
sys.exit()

if SAVEDATA:
    ts = nc.Dataset("/Users/cconn/Documents/Explore_controller/createControllerInputs/CalculatedData/" + VARIABLE + "_anomalies.nc", 'w' , format='NETCDF4')
    ts_members = ts.createDimension('members',len(files))
    ts_time = ts.createDimension('time',len(times))
    ts_lat = ts.createDimension('lat',len(lats))
    ts_lon = ts.createDimension('lon',len(lons))
     
    ts_instd = ts.createVariable(VARIABLE + '_anom','f4',('members','time','lat','lon'))
    
    ts_instd[:,:,:,:] = ensemble_anomalies
  
    
if SAVEDATA_BASESTATES:
    from netCDF4 import date2num,num2date
    ts = nc.Dataset("/Users/cconn/Documents/Explore_controller/createControllerInputs/CalculatedData/Basestates.nc", 'w' , format='NETCDF4')
    ts_time = ts.createDimension('time',len(times))
    # ts_time_units = 'hours since 1800-01-01'
    # dates = date2num(times, ts_time_units, 'noleap')
    ts_lat = ts.createDimension('lat',len(lats))
    ts_lon = ts.createDimension('lon',len(lons))
     
    ts_ensemble_mean = ts.createVariable('ensemble_mean','f4',('time','lat','lon'))
    ts_basestate = ts.createVariable('basestate','f4',('time','lat','lon'))
    ts_time = ts.createVariable('time','f4',('time'))
    ts_lat = ts.createVariable('lat','f4',('lat'))
    ts_lon = ts.createVariable('lon','f4',('lon'))
    
    ts_ensemble_mean[:,:,:] = ensemble_mean
    ts_basestate[:,:,:] = basestate
    ts_time[:] = times
    ts_lat[:] = lats
    ts_lon[:] = lons

    ts.close()    
    
    
    
