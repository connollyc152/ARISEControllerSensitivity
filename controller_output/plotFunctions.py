import numpy as np
import matplotlib.pyplot as plt

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

def PlotLines(x, y, x_label, y_label, title, color, label, line_style, alpha, lim):
    plt.axhline(0, color = 'black', linewidth=2, linestyle='--', alpha = .5)
    if len(x) == 0:
        return
    plt.plot(x, y, color = color, linewidth=3, label = label, linestyle=line_style, alpha = alpha)
    plt.title(title, **title_font)
    plt.xlabel(x_label, **axis_label_font)
    plt.ylabel(y_label, **axis_label_font)
    plt.yticks(**axisTick_font)        
    plt.xticks(**axisTick_font)
    plt.ylim(lim[0],lim[1])
    plt.xlim(min(x), max(x))
    
def PlotColormesh(x, y, data, x_label, y_label, title):
    plt.figure(dpi = (300), figsize = (6, 5))
    cs = plt.pcolormesh(x, y, data, cmap = tmap, vmin = -100, vmax = 100)
    plt.title(title, **title_font)
    cbar = plt.colorbar(cs, ticks = [-100, -50, 0, 50, 100])
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), **colorbar_font)
    cbar.set_label(label='Percent Change', **colorbar_font)
    plt.xlabel(x_label, **axis_label_font)
    plt.ylabel(y_label, **axis_label_font)
    plt.yticks(**axisTick_font)        
    plt.xticks(**axisTick_font)
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.show()
    
    
