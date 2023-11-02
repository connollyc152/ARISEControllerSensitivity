import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point

##########################
######PLOTTING INFO#######
##########################

loc_title = ["Injecting at 30S","Injecting at 15S","Injecting at 15N","Injecting at 30N",]

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

title_font = {'fontname':font, 'size':'15', 'color':'black', 'weight':'bold',
              'verticalalignment':'bottom'}
axis_label_font = {'fontname':font, 'size':'13', 'color':'black', 'weight':'bold'}
colorbar_font = {'fontname':font, 'size':'10', 'color':'black', 'weight':'light'}
legend_font = {'family':font, 'size':'12', 'weight':'light'}
axisTick_font = {'family':font, 'size':'10', 'color':'black', 'weight':'light'}

##########################
########FUNCTIONS#########
##########################

def Plot_Gobal_Map(plot_data, lats, lons, title, levels, colorbar):
    fig = plt.figure(dpi = 300, figsize = (6, 5))#, projection= ccrs.PlateCarree())
    ax=plt.axes(projection= ccrs.PlateCarree())
    plt.title(title, **title_font)
    data=plot_data
    data, lonsr = add_cyclic_point(data, coord=lons)
    cs=ax.contourf(lonsr, lats, data,
                transform = ccrs.PlateCarree(),extend='both', cmap = colorbar,
                levels = levels)

    ax.coastlines()
    cbar_ax = fig.add_axes([0.15, 0.2, 0.74, 0.03])

    cbar = plt.colorbar(cs,cax=cbar_ax,orientation='horizontal',label='Surface Air Temperature (K)', format='%.0f',
                        ticks = [-3, -2, -1, 0, 1, 2, 3])#, pad=5)
    cbar.set_label("(K)", **colorbar_font)
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), **colorbar_font)
    
    
    plt.show()

def PlotLines(x, y, x_label, y_label, title, color, label, line_style, alpha, lim):
    plt.axhline(0, color = 'black', linewidth=2, linestyle='--', alpha = .5)
    if len(x) == 0:
        return
    plt.plot(x, y, color = color, linewidth=3, label = label, linestyle=line_style, alpha = alpha)
    plt.title(title, pad=14, **title_font)
    plt.xlabel(x_label, **axis_label_font)
    plt.ylabel(y_label, **axis_label_font)
    plt.yticks(**axisTick_font)        
    plt.xticks(**axisTick_font)
    plt.ylim(lim[0],lim[1])
    plt.xlim(min(x), max(x))
    
def PlotColormesh(x, y, data, x_label, y_label, title):
    plt.figure(dpi = (300), figsize = (6, 5))
    cs = plt.pcolormesh(x, y, data, cmap = tmap, vmin = -100, vmax = 100)
    # for i in cs.get_array():
    #     if i <= 0:
    #         cs.set_edgecolors([1,0,0,0])
    plt.title(title, **title_font)
    cbar = plt.colorbar(cs, ticks = [-100, -50, 0, 50, 100])
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), **colorbar_font)
    cbar.set_label(label='Percent Change', **colorbar_font)
    plt.xlabel(x_label, **axis_label_font)
    plt.ylabel(y_label, **axis_label_font)
    plt.yticks(ticks = [-2, -1, 0, 1, 2], labels = ["-2", "-1", "0", "1", "2"], **axisTick_font)        
    plt.xticks(ticks = [-2, -1, 0, 1, 2], labels = ["-2", "-1", "0", "1", "2"], **axisTick_font)
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    # plt.show()
    
    y_lines = np.arange(-2,2.2,.2)
    x_lines = np.arange(-2,2.2,.2)
    print(np.shape(x_lines))
    
    recon = cs.get_array()
    recon = recon.reshape(20,20)
    
    for xi, x_grid in enumerate(recon[:-1,0]):
        for yi, y_grid in enumerate(recon[0,:]):
            if np.abs(recon[xi,yi])/recon[xi,yi] != np.abs(recon[(xi + 1),yi])/recon[(xi + 1),yi]:
                plt.hlines(y=y_lines[xi + 1], xmin=x_lines[yi], xmax=x_lines[yi + 1], color = "black", linewidth = .5)
                
    
    for yi, y_grid in enumerate(recon[0,:-1]):
        for xi, x_grid in enumerate(recon[:,0]):
            if np.abs(recon[xi,yi])/recon[xi,yi] != np.abs(recon[(xi),(yi + 1)])/recon[(xi),(yi + 1)]:
                plt.vlines(x=x_lines[yi + 1], ymin=y_lines[xi], ymax=y_lines[xi + 1], color = "black", linewidth = .5)
                
def PlotColormeshSmall(x, y, data, x_label, y_label, title):
    plt.figure(dpi = (300), figsize = (6, 5))
    cs = plt.pcolormesh(x, y, data, cmap = tmap)#, vmin = -1, vmax = 1)
    # cs = plt.pcolormesh(x, y, data, cmap = tmap, vmin = 288, vmax = 290)
    # for i in cs.get_array():
    #     if i <= 0:
    #         cs.set_edgecolors([1,0,0,0])
    plt.title(title, **title_font)
    cbar = plt.colorbar(cs)#, ticks = [-1, -.5, -.1, -.05, 0, .05, .10, .5, 1])
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), **colorbar_font)
    cbar.set_label(label='Percent Change', **colorbar_font)
    plt.xlabel(x_label, **axis_label_font)
    plt.ylabel(y_label, **axis_label_font)
    plt.yticks(ticks = [-2, -1, 0, 1, 2], labels = ["-2", "-1", "0", "1", "2"], **axisTick_font)        
    plt.xticks(ticks = [-2, -1, 0, 1, 2], labels = ["-2", "-1", "0", "1", "2"], **axisTick_font)
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    # plt.show()
    
    y_lines = np.arange(-2,2.2,.2)
    x_lines = np.arange(-2,2.2,.2)
    print(np.shape(x_lines))
    
    recon = cs.get_array()
    recon = recon.reshape(20,20)
    
    for xi, x_grid in enumerate(recon[:-1,0]):
        for yi, y_grid in enumerate(recon[0,:]):
            if np.abs(recon[xi,yi])/recon[xi,yi] != np.abs(recon[(xi + 1),yi])/recon[(xi + 1),yi]:
                plt.hlines(y=y_lines[xi + 1], xmin=x_lines[yi], xmax=x_lines[yi + 1], color = "black", linewidth = .5)
                
    
    for yi, y_grid in enumerate(recon[0,:-1]):
        for xi, x_grid in enumerate(recon[:,0]):
            if np.abs(recon[xi,yi])/recon[xi,yi] != np.abs(recon[(xi),(yi + 1)])/recon[(xi),(yi + 1)]:
                plt.vlines(x=x_lines[yi + 1], ymin=y_lines[xi], ymax=y_lines[xi + 1], color = "black", linewidth = .5)
    
    
