import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from cartopy.util import add_cyclic_point
from netCDF4 import date2num,num2date
import numpy.linalg as LA
import netCDF4 as nc
import matplotlib.colors
lower = plt.cm.RdBu_r(np.linspace(0,.49, 49))
white = plt.cm.RdBu_r(np.ones(2)*0.5)
upper = plt.cm.RdBu_r(np.linspace(0.51, 1, 49))
colors = np.vstack((lower, white, upper))
tmap = matplotlib.colors.LinearSegmentedColormap.from_list('terrain_map_white', colors)

import pandas as pd
times = pd.date_range(start='01/01/2035', end='12/31/2069', freq='M')
n = 100

def Plot_SAM_Map(plot_data, title, levels, colorbar):
    ax = plt.axes(projection=ccrs.Robinson())#central_longitude=180))
    ax.coastlines()
    ax.gridlines()
    plot_data.plot(
        ax=ax, transform=ccrs.PlateCarree(), cmap='coolwarm'
    )
    #ax.set_extent((100, 360, 20, 80))
    plt.title(title)
    plt.show()


def nandot(X,Y):
    
    C = np.empty([np.size(X,axis=0),np.size(Y,axis=1)])
    for row in np.arange(0,np.size(X,axis=0)):
        for col in np.arange(0,np.size(Y,axis=1)):
            C[row,col] = np.nanmean(np.multiply(X[row,:],Y[:,col]))
            
    return C

    
def Add_Metadata_4(data, lats, lons):
    data = xr.DataArray(data,
        dims = ['member','time','lat','lon'],
        coords=dict(
            member = (range(0,n)),
            time = (times), 
            # lat = (lats.sel(lat=slice(-5, 5))), 
            # lon = (lons.sel(lon=slice(190, 240)))
            lat = (lats), 
            lon = (lons)
            )
        )
    return data

def Add_Metadata_2_SAM(data, lat, lon):
    data = xr.DataArray(data,
        dims = ['lat','lon'],
        coords=dict(
            lat = lat, 
            lon = lon
            )
        )
    return data


def Save_SAM_Anom(anoms):
    #Saves SAM index values
    
    ts = nc.Dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/SAM_anomalies.nc", 'w' , format='NETCDF4')
    ts_member = ts.createDimension('member',n)
    ts_time = ts.createDimension('time',len(times))
     
    ts_anomalies = ts.createVariable('SAM_anomalies','f4',('member','time'))
    ts_time = ts.createVariable('time','f4',('time'))
    
    ts_anomalies[:,:] = anoms
    ts_time[:] = times

    ts.close()
    
def Calc_SAM(ensemble_array, lats, lons):
    
    members_n = len(ensemble_array[:,0,0,0])
    ensemble_array = Add_Metadata_4(ensemble_array, lats, lons)
    
    #Selects the region used to calculate the SAM index
    SAM_region = ensemble_array.sel(lat=slice(-90, -20))
    lat_region = lats.sel(lat=slice(-90, -20))

    #Create empty array which will be filled with SAM index values
    SAM_index_members = np.empty(shape = (members_n, len(ensemble_array[0,:,0,0])))
    
    for member in np.arange(0,members_n,1):
        data = np.array(SAM_region[member, :, :, :])/100   
        data_f = data.reshape(len(data[:,0,0]), len(data[0,:,0]) * len(data[0,0,:]))  
        
        #Calculate EOF
        C = nandot(data_f, np.transpose(data_f))   
        lam, Z = LA.eig(C)
        Z = (Z - np.nanmean(Z,axis=0))/np.nanstd(Z,axis=0)
        E = np.dot(Z.T,data_f)
        
        D = nandot(Z[:,:10].T,data.reshape(data.shape[0],data.shape[1]*data.shape[2]))
        xplot = D.reshape(D.shape[0],len(data[0,:,0]),len(data[0,0,:]))[0,:,:]
        xplot = Add_Metadata_2_SAM(xplot, lat_region, lons)
        
        plt.title(str(member) + " before fix")    
        cs = plt.contourf(lons, lat_region, xplot, cmap = tmap, levels = np.linspace(-15, 15, 20))
        plt.plot(52, -80 , 'o', color = 'red')
        plt.plot(52, -40, 'o', color = 'blue')
        plt.colorbar(cs)
        plt.show()
        
        #This checks that negative EOFs represent the negative phase of the SAM
        if xplot[50, -80] <= 0 and xplot[20, -40] >=0:
            print("fix")
            Z = Z * -1
            xplot = xplot * -1
            
        plt.title(str(member) + " after fix")    
        cs = plt.contourf(lons, lat_region, xplot, cmap = tmap, levels = np.linspace(-15, 15, 20))
        plt.colorbar(cs)
        plt.show()

        SAM_index_members[member, :] = Z[:,0]
        
    plt.plot(SAM_index_members[0,:])
    plt.plot(SAM_index_members[1,:])
    plt.show()

    Save_SAM_Anom(SAM_index_members)
    return(SAM_index_members)
        
