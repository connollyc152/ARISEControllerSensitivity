#! /usr/bin/env python
# Explicit feedback for climate modeling
# Copyright (C) 2020  Ben Kravitz
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Main setup script
#
# Written by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
# Last updated 11 July 2019
#
# This script provides the main information about the run.  Note that everything
# here is specific to an individual run.  Each run should have its own copy
# of the feedback suite.
#
# This script is written in native python format.  Be careful with brackets [],
# white space, and making sure everything that needs to be a string is actually
# a string by putting it in quotes ''.  All lines beginning with # are comments.
#
# Things you need to specify:

#TAKES INPUT AS KELVIN

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point

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


casepath='/Users/cconn/Documents/Explore_controller/controller'#+runname
outputpath = "/Users/cconn/Documents/Explore_controller/controller_output"
maindir=casepath
frequency='1y'

pathtocontrol=maindir+ 'PIcontrol.py'
driver_path_name = '/Users/cconn/Documents/Explore_controller/controller/driver.py'

OPEN_AND_PUSH = False

def Controller_Injection(vals, runname, ADD_FILEHEADERS, METADATA):
    runnames=runname
    casepath='/Users/cconn/Documents/Explore_controller/controller/'#+runname
    maindir=casepath
    frequency='1y'

    pathtocontrol=maindir+ 'PIcontrol.py'
    driver_path_name = '/Users/cconn/Documents/Explore_controller/controller/driver.py'
    
    data_file = "/Users/cconn/Documents/ARISE_data/TREFHT_annual/BWSSP245cmip6_005_TREFHT_201501-206412_206501-210012_annual.nc"
    f = xr.open_dataset(data_file)

    times = f["time"][80:81] #86
    lats = f['lat']
    lons = f['lon']
    
    # Plot_Gobal_Map(np.mean(vals[:,:,:], axis = 0), "Input CI", None, None, lats, lons)
    
    times =  times.values
    lats = lats.values
    lons = lons.values
    vals = vals.values
    
    f.close()
    
    local_dict = locals()
    
    execfile(driver_path_name, local_dict) 
    
    

if OPEN_AND_PUSH:
    runnames="trial"
    data_file = "/Users/cconn/Documents/ARISE_data/TREFHT_annual/BWSSP245cmip6_005_TREFHT_201501-206412_206501-210012_annual.nc"
    f = xr.open_dataset(data_file)
    data_type = 'filler'

    times = f["time"][80:81] #86
    lats = f['lat']
    lons = f['lon']
    vals = f["TREFHT"][80:81,:,:]

    if len(vals[:,0,0]) > 1:
        print("Inputing more than one map")

    Plot_Gobal_Map(np.mean(vals[:,:,:], axis = 0), "Input", None, None, lats, lons)
    
    times =  times.values
    lats = lats.values
    lons = lons.values
    vals = vals.values
    
    execfile(driver_path_name)



