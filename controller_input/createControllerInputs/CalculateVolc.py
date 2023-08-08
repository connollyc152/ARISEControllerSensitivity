#author: Emily Prewett
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.colors
import xarray as xr
import numpy as np
import os 
import sys
import pandas as pd
import scipy as sc
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
import gc
gc.enable()
import warnings
warnings.simplefilter('ignore',  RuntimeWarning)

# file_folder_path = (r'C:\Users\emily\Desktop\Python\Barnes\Project2\data')
file_folder_path = "/Users/cconn/Documents/CESM2-LE/HIST/TREFHT/"

files = []
for filename in os.listdir(file_folder_path):
    if filename.endswith('.nc') and filename.startswith('b.'):
        files.append(os.path.join(filename))
  
files = list(np.sort(files))

mean_part = np.empty(shape = (1980, 192, 288))  
mean_part.fill(0) 

count = 100

for i, file in enumerate(files[:count]):

    # f = xr.open_dataset(file_folder_path + '\\' + file)
    f = xr.open_dataset(file_folder_path + file)
    print(i)
    
    times = f["time"][:]
    lats = f['lat']
    lons = f['lon']
    data = f["TREFHT"][:,:,:]
    
    f.close()
    
    mean_part = mean_part + data

ensemble_mean = mean_part/count
climatology = ensemble_mean.groupby("time.month").mean("time")
ensemble_mean = ensemble_mean.groupby("time.month") - climatology

ensemble_mean_10before = ensemble_mean.sel(time=slice("1981-07", "1991-06"))
months_10before = times.sel(time=slice("1981-07", "1991-06"))
months_10before = months_10before.indexes['time'].to_datetimeindex()

months_before_and_after = ensemble_mean.sel(time=slice("1981-07", "1993-06"))

x = np.arange(0,120,1)
x2 = np.arange(0,144,1)

date1 = np.arange(np.datetime64('1981-06'), np.datetime64('1991-06'))
date2 = np.arange(np.datetime64('1981-06'), np.datetime64('1993-06'))

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)

vol_map = np.empty(shape = (192, 288))

for lon in np.arange(0, len(lons), 1):
    print(lon)
    for lat in np.arange(0, len(lats), 1):

        plot_10before = ensemble_mean_10before[:,lat,lon]
        res = sc.stats.linregress(x, plot_10before)
        slope, intercept, r_value, p_value, std_err = sc.stats.linregress(x, plot_10before)
        m = slope
        b = intercept
        y = (m * x) + b
        
        months_avg = months_before_and_after[:,lat,lon]
        y2 = (m * x2) + b
        
        temp_Avolc = months_avg[-24:]
        time_Avolc = date2[-24:]
        line_Avolc = y2[-24:]
        
        diff = temp_Avolc - line_Avolc
        
        point = np.mean(diff)

        vol_map[lat,lon] = float(point)

lower = plt.cm.RdBu_r(np.linspace(0,.49, 49))
white = plt.cm.RdBu_r(np.ones(2)*0.5)
upper = plt.cm.RdBu_r(np.linspace(0.51, 1, 49))
colors = np.vstack((lower, white, upper))
tmap = matplotlib.colors.LinearSegmentedColormap.from_list('terrain_map_white', colors)

# vol_map = xr.DataArray(vol_map,
#         dims = ['lat','lon'], # new column names
#         coords=dict(
#             # time = ("time", pd.date_range(start = "01-01-1850", end = "12-31-2014", freq="M")), #frequency for time (change dates)
#             lat = (lats),
#             lon = (lons)
#             )
#         )
# vol_map.to_netcdf(path=r"/Users/cconn/Desktop/VolMap100.nc", mode='w', format='NetCDF4')

def Plot_Gobal_Map(plot_data, lats, lons, title, levels, colorbar):
    plt.figure(dpi = 300, figsize = (6, 5))#, projection= ccrs.PlateCarree())
    ax=plt.axes(projection= ccrs.PlateCarree())
    data=plot_data
    data, lonsr = add_cyclic_point(data, coord=lons)
    cs=ax.contourf(lonsr, lats, data,
                transform = ccrs.PlateCarree(),extend='both', cmap = colorbar,
                levels = levels)

    ax.coastlines()
    cbar = plt.colorbar(cs,shrink=0.7,orientation='horizontal',label='Surface Air Temperature (K)', format='%.1f')#, pad=5)
    cbar.set_label("temperature anomaly (K)")
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels())    
    
    plt.title(title)

    plt.show()

Plot_Gobal_Map(vol_map, lats, lons, "Pinatubo", np.linspace(-2,2,21), tmap)