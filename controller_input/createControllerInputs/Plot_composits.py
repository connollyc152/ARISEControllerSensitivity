import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
import numpy as np

import matplotlib.colors
lower = plt.cm.RdBu_r(np.linspace(0,.49, 49))
white = plt.cm.RdBu_r(np.ones(2)*0.5)
upper = plt.cm.RdBu_r(np.linspace(0.51, 1, 49))
colors = np.vstack((lower, white, upper))
tmap = matplotlib.colors.LinearSegmentedColormap.from_list('terrain_map_white', colors)

fontPath = '/Users/cconn/Library/Fonts/'
import matplotlib.font_manager as fm
for font in fm.findSystemFonts(fontPath):
     fm.fontManager.addfont(font)

font =  'Poppins'#'Comme'#'Noto Sans JP'

title_font = {'fontname':font, 'size':'25', 'color':'black', 'weight':'bold',
              'verticalalignment':'bottom'}
axis_label_font = {'fontname':font, 'size':'13', 'color':'black', 'weight':'bold'}
colorbar_font = {'fontname':font, 'size':'12', 'color':'black', 'weight':'light'}
legend_font = {'family':font, 'size':'12', 'weight':'light'}
axisTick_font = {'family':font, 'size':'10', 'color':'black', 'weight':'light'}

def COMPOSITES_SIX(lat, lon, ENSO_map, NAO_map, VOL_maps, SAM_maps, BASE_map, FINAL_map, METADATA):
    fig, axs = plt.subplots(3, 2, figsize=(15, 10), dpi = 300, subplot_kw={'projection': ccrs.PlateCarree()})
    plt.subplots_adjust(wspace= -0.2)  
    levels_anom = np.linspace(-1,1,16) #np.linspace(-2.5,2.5,16)
    ticks_temp = [-1, 0, 1]
    
    ENSO_map, lonr = add_cyclic_point(ENSO_map, coord=lon)
    cs = axs[0,0].contourf(lonr, lat, ENSO_map, cmap = tmap, levels = levels_anom, extend="both")
    axs[0,0].set_title("ENSO", loc='left', **title_font)
    axs[0,0].coastlines()
    axs[0,0].text(180,112, "index [" + str(METADATA[1]) + ":" + str(METADATA[2]) + "]", ha = "right", **colorbar_font)
    axs[0,0].text(180,96, "n: " + str(METADATA[3]), ha = "right", **colorbar_font)
    
    NAO_map, lonr = add_cyclic_point(NAO_map, coord=lon)
    cs = axs[0,1].contourf(lonr, lat, NAO_map, cmap = tmap, levels = levels_anom, extend="both")
    axs[0,1].set_title("NAO", loc='left', **title_font)
    axs[0,1].coastlines()
    axs[0,1].text(180,112, "index [" + str(METADATA[5]) + ":" + str(METADATA[6]) + "]", ha = "right", **colorbar_font)
    axs[0,1].text(180,96, "n: " + str(METADATA[7]), ha = "right", **colorbar_font)
    
    VOL_maps, lonr = add_cyclic_point(VOL_maps, coord=lon)
    cs = axs[1, 0].contourf(lonr, lat, VOL_maps, cmap = tmap, levels = levels_anom, extend="both")
    axs[1,0].set_title("Volcanoes", loc='left', **title_font)
    axs[1,0].coastlines()
    axs[1,0].text(180,112, "volcano: none" , ha = "right", **colorbar_font)
    axs[1,0].text(180,96, "n:  0", ha = "right", **colorbar_font)
    box = axs[1,0].get_position()
    box.y0 = box.y0 - .02
    box.y1 = box.y1 - 0.02
    axs[1,0].set_position(box)
    
    SAM_maps, lonr = add_cyclic_point(SAM_maps, coord=lon)
    cs = axs[1, 1].contourf(lonr, lat, SAM_maps, cmap = tmap, levels = levels_anom, extend="both")
    axs[1,1].set_title("SAM", loc='left', **title_font)
    axs[1,1].coastlines()
    axs[1,1].text(180,112, "Index: " , ha = "right", **colorbar_font)
    axs[1,1].text(180,96, "n:  0", ha = "right", **colorbar_font)
    box = axs[1,1].get_position()
    box.y0 = box.y0 - .02
    box.y1 = box.y1 - 0.02
    axs[1,1].set_position(box)
    
    cbar_ax = fig.add_axes([0.19, 0.30, 0.64, 0.03])
    cbar = fig.colorbar(cs, cax=cbar_ax, orientation="horizontal", ticks = ticks_temp)
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), **colorbar_font)
    cbar.set_label(label='Temperature Anomalies (K)', **colorbar_font)
    
    ################################# 
    levels_temp = np.linspace(230,310,10)
    
    BASE_map, lonr = add_cyclic_point(BASE_map, coord=lon)
    cs = axs[2,0].contourf(lonr, lat, BASE_map, cmap = "Reds", levels = levels_temp, extend="both")
    axs[2,0].set_title("Basestate", loc='left', **title_font)
    axs[2,0].coastlines()
    axs[2,0].text(180,112, "Year: " + METADATA[0] , ha = "right", **colorbar_font)
    axs[2,0].text(180,96, "n:  10", ha = "right", **colorbar_font)
    box = axs[2,0].get_position()
    box.y0 = box.y0 - .16
    box.y1 = box.y1 - 0.16
    axs[2,0].set_position(box)
    
    FINAL_map = FINAL_map[0,:,:]
    FINAL_map, lonr = add_cyclic_point(FINAL_map, coord=lon)
    cs = axs[2,1].contourf(lonr, lat, FINAL_map, cmap = "Reds", levels = levels_temp, extend="both")
    axs[2,1].set_title("Complete Map", loc='left', **title_font)
    axs[2,1].coastlines()
    box = axs[2,1].get_position()
    box.y0 = box.y0 - .16
    box.y1 = box.y1 - 0.16
    axs[2,1].set_position(box)
    
    cbar_ax = fig.add_axes([0.19, -0.1, 0.64, 0.03])
    cbar = fig.colorbar(cs, cax=cbar_ax, orientation="horizontal")
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), **colorbar_font)
    cbar.set_label(label='Temperature (K)', **colorbar_font)
    
    plt.show()

    # fig, axs = plt.subplots(1, 1, figsize=(15, 10), dpi = 300, subplot_kw={'projection': ccrs.PlateCarree()})
    # cs = axs.contourf(lonr, lat, (FINAL_map - BASE_map), cmap = tmap)
    # axs.coastlines()
    # plt.colorbar(cs)
    
    # plt.show()
    