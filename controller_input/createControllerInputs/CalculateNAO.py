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

def Plot_NAO_Map(plot_data, title, levels, colorbar):
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
    # print(lats)
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

def Add_Metadata_2_NAO(data):
    # print(lats)
    data = xr.DataArray(data,
        dims = ['lat','lon'],
        coords=dict(
            lat = (np.linspace(20,80,64)), 
            lon = (np.linspace(-90, 40, 105))
            )
        )
    return data


def Save_NAO_Anom(anoms):
    #Saves NAO index values
    
    ts = nc.Dataset("/Users/cconn/Documents/Explore_controller/controller_input/Data/NAO_anomalies.nc", 'w' , format='NETCDF4')
    ts_member = ts.createDimension('member',len(anoms[:,0]))
    ts_time = ts.createDimension('time',len(times))
     
    ts_anomalies = ts.createVariable('NAO_anomalies','f4',('member','time'))
    ts_time = ts.createVariable('time','f4',('time'))
    
    ts_anomalies[:,:] = anoms
    ts_time[:] = times

    ts.close()
    
def Calc_NAO(ensemble_array, lats, lons):
    ensemble_array = Add_Metadata_4(ensemble_array, lats, lons)
    
    #Selects the region used to calculate the NAO index
    NAO_region1 = ensemble_array.sel(lat=slice(20, 80), lon=slice(-90, 40))
    NAO_region2 = ensemble_array.sel(lat=slice(20, 80), lon=slice(270, 360))

    #Combines the two NAO regions
    NAO_region = xr.concat([NAO_region2, NAO_region1], dim="lon")
    
    NAO_region["lon"] = (np.linspace(-90, 40, 105))
    # Plot_NAO_Map(NAO_region[0,0,:,:], "Final Section", None, None)
    
    #Create empty array which will be filled with NAO index values
    NAO_index_members = np.empty(shape = (len(ensemble_array[:,0,0,0]), len(ensemble_array[0,:,0,0])))
    
    for member in np.arange(0,len(ensemble_array[:,0,0,0]),1):
        data = np.array(NAO_region[member, :, :, :])/100     
        data_f = data.reshape(len(data[:,0,0]), len(data[0,:,0]) * len(data[0,0,:]))     
        C = nandot(data_f, np.transpose(data_f))
        
        #Calculate EOF
        lam, Z = LA.eig(C)
        Z = (Z - np.nanmean(Z,axis=0))/np.nanstd(Z,axis=0)
        E = np.dot(Z.T,data_f)
        
        D = nandot(Z[:,:10].T,data.reshape(data.shape[0],data.shape[1]*data.shape[2]))
        xplot = D.reshape(D.shape[0],len(data[0,:,0]),len(data[0,0,:]))[0,:,:]
        xplot = Add_Metadata_2_NAO(xplot)
        
        #This checks that negative EOFs represent the negative phase of the NAO
        if xplot[50, 52] >= 0 and xplot[20, 52] <=0:
            print("fix")
            Z = Z * -1
            xplot = xplot * -1
            
        cs = plt.contourf(xplot, cmap = tmap, levels = np.linspace(-5, 5, 20))
        plt.colorbar(cs)
        plt.plot(52, 20, 'ro', color = 'red')
        plt.plot(52, 50, 'ro', color = 'blue')
        plt.show()

        NAO_index_members[member, :] = Z[:,0]

    Save_NAO_Anom(NAO_index_members)

    return(NAO_index_members)
        


