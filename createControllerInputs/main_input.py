import sys
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
from netCDF4 import date2num,num2date
sys.path.insert(0, '/Users/cconn/Documents/Explore_controller/createControllerInputs/')
import CalculateENSO as calc_ENSO
sys.path.insert(0, '/Users/cconn/Documents/Explore_controller/controller')
import main_controller as  MC
# import CalculateBasestate as BS

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

f = xr.open_dataset("/Users/cconn/Documents/Explore_controller/createControllerInputs/CalculatedData/Basestates.nc")
times = f["time"][:660] #86
times = num2date(times[:],'hours since 1800-01-01', 'noleap')

lats = f['lat']
lons = f['lon']
basestates = f['basestate']

TREFTHT_anomalies_file = "/Users/cconn/Documents/Explore_controller/createControllerInputs/CalculatedData/TREFHT_anomalies.nc"
f = xr.open_dataset(TREFTHT_anomalies_file)
TREFHT_anom = f["TREFHT_anom"]

ENSO_index = np.array(calc_ENSO.Calc_ENSO34())

print(np.shape(ENSO_index))
print(ENSO_index[1,:])

plt.plot(ENSO_index[0,:])
plt.plot(ENSO_index[1,:])
plt.show()

# ENSO_interest = [-2.25,-2,-1.75,-1.5,-1.25,-1,-.75,-.5,-.25,0,.25,.5,.75,1,1.25,1.5,1.75,2,2.25,2.5]
# for e, en in enumerate(ENSO_interest[:-1]):
#     print(ENSO_interest[e], (ENSO_interest[e + 1]))
#     member_array = np.empty(shape=(10, 192, 288))
#     for member in np.arange(0,10,1):
#         member_array[member, :, :] = np.mean(TREFHT_anom[member,
#                                     (ENSO_index[member] < ENSO_interest[(e + 1)]) & (ENSO_index[member] > ENSO_interest[e]),
#                                     :,:], axis = 0)
    
    
#     ENSO = np.mean(member_array, axis = 0)
#     Plot_Gobal_Map(ENSO, "ENSO BETWEEN " + str(ENSO_interest[e]) + " and " + str(ENSO_interest[e + 1]), None, None)
    
#     ENSO = ENSO + basestates[400,:,:]
#     base = basestates[400,:,:]
    
#     ENSO = Add_Initial_dim(ENSO)
#     MC.Controller_Injection(ENSO, "ENSO_test", "basestate + ENSO " + str(en))

for i in np.arange(0,660,30):
    base = Add_Initial_dim(basestates[i,:,:])
    MC.Controller_Injection(base, "ClimateChnage", i)
    
    ENSO_interest = [-2,-1.5,-1,-.5,0,.5,1,1.5,2,2.5]
    for e, en in enumerate(ENSO_interest[:-1]):
        print(ENSO_interest[e], (ENSO_interest[e + 1]))
        member_array = np.empty(shape=(10, 192, 288))
        for member in np.arange(0,10,1):
            member_array[member, :, :] = np.mean(TREFHT_anom[member,
                                        (ENSO_index[member] < ENSO_interest[(e + 1)]) & (ENSO_index[member] > ENSO_interest[e]),
                                        :,:], axis = 0)
        
        
        ENSO = np.mean(member_array, axis = 0)
        Plot_Gobal_Map(ENSO, "ENSO BETWEEN " + str(ENSO_interest[e]) + " and " + str(ENSO_interest[e + 1]), None, None)
        
        ENSO = ENSO + base
        
        ENSO = Add_Initial_dim(ENSO)
        MC.Controller_Injection(ENSO, "ENSO_test_trial" + str(i), "basestate + ENSO " + str(en))


