import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from cartopy.util import add_cyclic_point
from netCDF4 import date2num,num2date

SST_files = ["Basestates.nc"]

datapath = "/Users/cconn/Documents/Explore_controller/createControllerInputs/CalculatedData/"

f = xr.open_dataset(datapath + SST_files[0])

times = f["time"][:660] #86
times = num2date(times[:],'hours since 1800-01-01', 'noleap')
lats = f['lat']
lons = f['lon']

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
    
def Add_Metadata_4(data):
    # print(lats)
    data = xr.DataArray(data,
        dims = ['member','time','lat','lon'],
        coords=dict(
            member = (range(0,10)),
            time = (times), 
            # lat = (lats.sel(lat=slice(-5, 5))), 
            # lon = (lons.sel(lon=slice(190, 240)))
            lat = (lats), 
            lon = (lons)
            )
        )
    return data

def Add_Metadata_4_cut(data):
    # print(lats)
    data = xr.DataArray(data,
        dims = ['member','time','lat','lon'],
        coords=dict(
            member = (range(0,10)),
            time = (times), 
            lat = (lats.sel(lat=slice(-5, 5))), 
            lon = (lons.sel(lon=slice(190, 240)))
            )
        )
    return data


def Calc_ENSO34():
    
    ensemble_array = np.empty(shape = (10, 660, 192, 288))
        
    f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/createControllerInputs/CalculatedData/SST_anomalies.nc")
    vals = f["SST_anom"]
    
    vals = np.array(vals)
    
    Plot_Gobal_Map(vals[0,0,:,:], 'SST', None, 'Reds')
    
    ensemble_array[:,:,:,:] = vals
    
    ensemble_array = Add_Metadata_4(ensemble_array)
    Plot_Gobal_Map(ensemble_array[0,0,:,:], 'Ensemble member 1 SST 01-01-2035', np.arange(0,30,1), None) 
    
    tos_nino34 = ensemble_array.sel(lat=slice(-5, 5), lon=slice(190, 240))
    
    Plot_ENSO_Map(tos_nino34[0,0,:,:], 'SST anomaly over the Nino 3.4 region', np.arange(0,30,1), None)  
    tos_nino34 = Add_Metadata_4_cut(tos_nino34)
    
    #DO I NEED THIS SINCE THEY ARE ALREADY ANOMALIES?
    # gb = tos_nino34.groupby('time.month')
    # tos_nino34_anom = gb - gb.mean(dim='time')
        
    # weights = np.cos(np.deg2rad(tos_nino34_anom.lat))
    # weights.name = "weights"
    # index_nino34 = tos_nino34_anom.weighted(weights).mean(("lat", "lon"))
    
    weights = np.cos(np.deg2rad(tos_nino34.lat))
    weights.name = "weights"
    index_nino34 = tos_nino34.weighted(weights).mean(("lat", "lon"))
    
    index_nino34_rolling_mean = index_nino34 #index_nino34.rolling(time=5, center=True).mean()
    
    plt.plot(index_nino34_rolling_mean[0,:])
    plt.plot(index_nino34_rolling_mean[1,:])
    plt.show()

    return(index_nino34_rolling_mean)

def Calc_ENSO34_single(SST):
    
    tos_nino34 = SST.sel(lat=slice(-5, 5), lon=slice(190, 240))
    
    # Plot_ENSO_Map(tos_nino34[0,0,:,:], 'SST anomaly over the Nino 3.4 region', np.arange(0,30,1), None)  
    # tos_nino34 = Add_Metadata_4_cut(tos_nino34)
    
    # gb = tos_nino34.groupby('time.month')
    # tos_nino34_anom = gb - gb.mean(dim='time')
    
    
    weights = np.cos(np.deg2rad(tos_nino34.lat))
    weights.name = "weights"
    index_nino34 = tos_nino34.weighted(weights).mean(("lat", "lon"))
    
    # index_nino34_rolling_mean = index_nino34.rolling(time=5, center=True).mean()
    
    # plt.plot(index_nino34_rolling_mean[0,:])
    # plt.plot(index_nino34_rolling_mean[1,:])
    # plt.show()

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

# f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/createControllerInputs/CalculatedData/SST_anomalies.nc")
# vals = f["SST_anom"]

# vals = np.array(vals)
# vals = Add_Metadata_4(vals)

# Plot_Gobal_Map(vals[1,0,:,:], 'SST', None, 'Reds')

# nino = cal_ninoCR(vals[1,0,:,:], 1, standardize = True)
# print(nino)

# nino = Calc_ENSO34_single(vals[1,0,:,:])
# print(nino)
