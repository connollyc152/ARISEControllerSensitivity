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

def Add_Metadata_4(data, lats, lons, times):
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
    
def CalculateBasestate_old(SAVE_ANOMALIES, SAVE_BASESTATES, data, VARIABLE,
                       lats, lons, times):
        
    if VARIABLE == "SST":
        data = data - 270
        
    #Calculate Ensemble mean    
    ensemble_mean = np.mean(data, axis = 0)

    #polyfit of the ensemble_mean (K)
    basestate = np.empty(shape = (660, 192, 288))
    for lat in np.arange(0, len(lats), 1):
        for lon in np.arange(0, len(lons), 1):
            x = np.arange(0, len(ensemble_mean[:,0,0]), 1)
            z = np.polyfit(x,ensemble_mean[:,lat,lon], deg = 3)

            basestate[:,lat,lon] = np.polyval(z,x)
            
        
    #Remove the climate change trend  
    detrended = np.empty(shape = (len(data[:,0,0,0]), 660, 192, 288))
    for member in range(0, len(data[:,0,0,0])):
        for lat in np.arange(0, len(lats), 1):
            for lon in np.arange(0, len(lons), 1):
        
                detrended[member, :, lat, lon] = data[member,:,lat,lon] - basestate[:, lat, lon]
        
    #Calculate anomalies from detrended data
    vals_anom = np.empty(shape = (len(data[:,0,0,0]), 660, 192, 288))
    for member in range(0, len(data[:,0,0,0])):  
        vals_mean = np.mean(detrended[member,:,:,:], axis = 0)
        
        vals_anom[member, :, :, :] = detrended[member,:,:,:] - vals_mean
    
    #Remove seasonal cycle from detrended anomaly data
    anom_noseas = np.empty(shape = (len(data[:,0,0,0]), 660, 192, 288))
    for member in range(0, len(data[:,0,0,0])):
        for lat in np.arange(0, len(lats), 1):
            for lon in np.arange(0, len(lons), 1):
                    data_temp = vals_anom[member,:,lat,lon]

                    Z = np.fft.fft(data_temp)
                    Yfft = Z/np.size(data_temp)
    
                    # combine symmetric parts of the FFT and plot the power spectrum as a function of frequency
                    freq = np.arange(0,np.size(data_temp)/2)/float(np.size(data_temp))
    
    
                    # the factor of 2 in front is needed or the sum won't equal the total variance of X
                    Ck2 = 2.*np.abs(Yfft[0:int(np.size(data_temp)/2)+1])**2 
    
                    Z[50:60] = 0.0
                    Z[-60:-50] = 0.0
                    
                    # plt.plot(Z)
                    # plt.show()
                
                    X_filtered = np.real(np.fft.ifft(Z))    
                
                    anom_noseas[member,:,lat,lon] = X_filtered #+ float(vals_mean[lat,lon]) 

def CalculateBasestate(SAVE_ANOMALIES, SAVE_BASESTATES, data, VARIABLE,
                       lats, lons, times):
    print("CALCULATING BASESTATE")
        
    if VARIABLE == "SST":
        data = data - 270
        
    #Calculate Ensemble mean    
    ensemble_mean = np.mean(data, axis = 0)

    #polyfit of the ensemble_mean (K)
    basestate = np.empty(shape = (660, 192, 288))
    for lat in np.arange(0, len(lats), 1):
        for lon in np.arange(0, len(lons), 1):
            x = np.arange(0, len(ensemble_mean[:,0,0]), 1)
            z = np.polyfit(x,ensemble_mean[:,lat,lon], deg = 3)

            basestate[:,lat,lon] = np.polyval(z,x)
            
        
    #Remove the climate change trend  
    detrended = np.empty(shape = (len(data[:,0,0,0]), 660, 192, 288))
    for member in range(0, len(data[:,0,0,0])):
        for lat in np.arange(0, len(lats), 1):
            for lon in np.arange(0, len(lons), 1):
        
                detrended[member, :, lat, lon] = data[member,:,lat,lon] - basestate[:, lat, lon]
        
    #Calculate anomalies from detrended data
    vals_anom = np.empty(shape = (len(data[:,0,0,0]), 660, 192, 288))
    for member in range(0, len(data[:,0,0,0])):  
        vals_mean = np.mean(detrended[member,:,:,:], axis = 0)
        
        vals_anom[member, :, :, :] = detrended[member,:,:,:] - vals_mean
    
    #Remove seasonal cycle from detrended anomaly data
    #anom_noseas = np.empty(shape = (len(data[:,0,0,0]), 660, 192, 288))
    
    vals_anom = Add_Metadata_4(vals_anom, lats, lons, times)
    climatology = vals_anom.groupby("time.month").mean("time")

    anom_noseas = vals_anom.groupby("time.month") - climatology
    
    # for member in range(0, len(data[:,0,0,0])):
    #     for lat in np.arange(0, len(lats), 1):
    #         for lon in np.arange(0, len(lons), 1):
    #             print(lon)
    #             plt.plot(anom_noseas[member,:,80,lon])
    #             plt.show()
                    
    #anom_noseas: anomalies detrended and seasonal cycle removed. 
    #basestate: smoothed climate chagne trend. 
    basestate = Add_Metadata_3(basestate, lats, lons, times)
    anom_noseas = Add_Metadata_4(anom_noseas, lats, lons, times)
    Plot_Gobal_Map(basestate[0,:,:], 'Bastestate 01-01-2035', None, 'Reds', lats, lons)
    Plot_Gobal_Map(anom_noseas[0,0,:,:], 'Anomalies member 1 temperature 01-01-2035', None, 'Reds', lats, lons)            
    
    if SAVE_ANOMALIES:
        from netCDF4 import date2num,num2date
        ts_time_units = 'hours since 1800-01-01'
        dates = date2num(times, ts_time_units, 'noleap')
        
        ts = nc.Dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/" + VARIABLE + "_anomalies.nc", 'w' , format='NETCDF4')
        ts_members = ts.createDimension('members',len(data))
        ts_time = ts.createDimension('time',len(times))
        ts_lat = ts.createDimension('lat',len(lats))
        ts_lon = ts.createDimension('lon',len(lons))
         
        ts_instd = ts.createVariable(VARIABLE + '_anom','f4',('members','time','lat','lon'))
        ts_time = ts.createVariable('time','f4',('time'))
        ts_lat = ts.createVariable('lat','f4',('lat'))
        ts_lon = ts.createVariable('lon','f4',('lon'))
        
        ts_instd[:,:,:,:] = anom_noseas
        ts_time[:] = dates
        ts_lat[:] = lats
        ts_lon[:] = lons
      
        
    if SAVE_BASESTATES:
        from netCDF4 import date2num,num2date
        ts_time_units = 'hours since 1800-01-01'
        dates = date2num(times, ts_time_units, 'noleap')
        
        ts = nc.Dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/" + VARIABLE + "_Basestates.nc", 'w' , format='NETCDF4')
        ts_time = ts.createDimension('time',len(times))
        ts_lat = ts.createDimension('lat',len(lats))
        ts_lon = ts.createDimension('lon',len(lons))
         
        ts_ensemble_mean = ts.createVariable('ensemble_mean','f4',('time','lat','lon'))
        ts_basestate = ts.createVariable('basestate','f4',('time','lat','lon'))
        ts_time = ts.createVariable('time','f4',('time'))
        ts_lat = ts.createVariable('lat','f4',('lat'))
        ts_lon = ts.createVariable('lon','f4',('lon'))
        
        ts_ensemble_mean[:,:,:] = ensemble_mean
        ts_basestate[:,:,:] = basestate
        ts_time[:] = dates
        ts_lat[:] = lats
        ts_lon[:] = lons

        ts.close()   
    
  
    
    
    
