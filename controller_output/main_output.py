import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '/Users/cconn/Documents/Explore_controller/controller_output/')
import plotFunctions as pF

##########################
########OPEN DATA#########
##########################

datapath = "/Users/cconn/Documents/Explore_controller/controller_output/data/ControlLog_"
datafile = "monthly_SSP45_T_LOCKED.txt"
datafile_LE = "monthly_LE_T_LOCKED.txt"

#ADD_LE determines whether the LE data is plotted alongside the ARISE data 
ADD_LE = False
#ADD_N_LABEL adds n samples to some plots
ADD_N_LABEL = False
PLOT_T_ERROR = False
#What base years do you want to plot. Data must exist in the .txt file
Base_years = [2035, 2045]
Colors = ['#667302', '#D97904', '#665A73', '#A6A26F', '#A33E5E', '#FA4BA8']
rangeV = np.arange(-2,2.2,.2)
# rangeV = np.arange(-2.4,2.7,.3)

#Creates figures from paper
PLOT_NAO_BASESTATE_TOT_N = False #Fig 2
PLOT_ENSO_BASESTATE_TOT_N = False #Fig 2
PLOT_SAM_BASESTATE_TOT_N = False #Fig 2

PLOT_NAO_BASESTATE_TOT_N_VOLC = False #Fig 4
PLOT_ENSO_BASESTATE_TOT_N_VOLC = False #Fig 4
PLOT_SAM_BASESTATE_TOT_N_VOLC = False #Fig 4

COMPARE_NAO_ENSO = True #Fig 3
COMPARE_NAO_SAM = False #Fig 3
COMPARE_ENSO_SAM = False #Fig 3

################################
########OPEN DATA SSP45#########
################################
#Opens .txt file and reads each column into a specific variable
#this section is specifically for SSP245 data
data = np.array([x.split(' ') for x in open (datapath + datafile).readlines()])[:,:-1]

header = data[0,:]
data = data[1:,:]

base_years = np.asarray(data[:,header == 'BASESTATE_YEAR'], float)[:,0]

min_ENSO = np.asarray(data[:,header == 'ENSO_MIN'], float)[:,0]
max_ENSO = np.asarray(data[:,header == 'ENSO_MAX'], float)[:,0]
avg_ENSO = np.mean([min_ENSO, max_ENSO], axis = 0)
n_ENSO = np.asarray(data[:,header == 'ENSO_n'], float)[:,0]
ENSO_q = np.array(data[:,header == 'ENSO_q'])

min_NAO = np.asarray(data[:,header == 'NAO_MIN'], float)[:,0]
max_NAO = np.asarray(data[:,header == 'NAO_MAX'], float)[:,0]
avg_NAO = np.mean([min_NAO, max_NAO], axis = 0)
n_NAO = np.asarray(data[:,header == 'NAO_n'], float)[:,0]
NAO_q = np.array(data[:,header == 'NAO_q'])

min_SAM = np.asarray(data[:,header == 'SAM_MIN'], float)[:,0]
max_SAM = np.asarray(data[:,header == 'SAM_MAX'], float)[:,0]
avg_SAM = np.mean([min_SAM, max_SAM], axis = 0)
n_SAM = np.asarray(data[:,header == 'SAM_n'], float)[:,0]
SAM_q = np.array(data[:,header == 'SAM_q'])

S30_inj = np.asarray(data[:,header == '30S(Tg)'], float)
S15_inj = np.asarray(data[:,header == '15S(Tg)'], float)
N15_inj = np.asarray(data[:,header == '15N(Tg)'], float)
N30_inj = np.asarray(data[:,header == '30N(Tg)'], float)
total_inj = S30_inj + S15_inj + N30_inj + N15_inj

T0 = np.asarray(data[:,header == 'T0'], float)[:,0]
T1 = np.asarray(data[:,header == 'T1'], float)[:,0]
T2 = np.asarray(data[:,header == 'T2'], float)[:,0]

T0_error = np.asarray(data[:,header == 'T0_error'], float)[:,0]
T1_error = np.asarray(data[:,header == 'T1_error'], float)[:,0]
T2_error = np.asarray(data[:,header == 'T2_error'], float)[:,0]

T_Total_Error = T0_error + T1_error + T2_error

VOLC_q = np.array(data[:,header == 'VOLC_q'])
#############################
########OPEN DATA LE#########
#############################
#Opens .txt file and reads each column into a specific variable
#this section is specifically for large ensemble data
data_LE = np.array([x.split(' ') for x in open (datapath + datafile_LE).readlines()])[:,:-1]

data_LE = data_LE[1:,:]

base_years_LE = np.asarray(data_LE[:,header == 'BASESTATE_YEAR'], float)[:,0]

min_ENSO_LE = np.asarray(data_LE[:,header == 'ENSO_MIN'], float)[:,0]
max_ENSO_LE = np.asarray(data_LE[:,header == 'ENSO_MAX'], float)[:,0]
avg_ENSO_LE = np.mean([min_ENSO_LE, max_ENSO_LE], axis = 0)
n_ENSO_LE = np.asarray(data_LE[:,header == 'ENSO_n'], float)[:,0]
ENSO_q_LE = np.array(data_LE[:,header == 'ENSO_q'])

min_NAO_LE = np.asarray(data_LE[:,header == 'NAO_MIN'], float)[:,0]
max_NAO_LE = np.asarray(data_LE[:,header == 'NAO_MAX'], float)[:,0]
avg_NAO_LE = np.mean([min_NAO_LE, max_NAO_LE], axis = 0)
n_NAO_LE = np.asarray(data_LE[:,header == 'NAO_n'], float)[:,0]
NAO_q_LE = np.array(data_LE[:,header == 'NAO_q'])

min_SAM_LE = np.asarray(data_LE[:,header == 'SAM_MIN'], float)[:,0]
max_SAM_LE = np.asarray(data_LE[:,header == 'SAM_MAX'], float)[:,0]
avg_SAM_LE = np.mean([min_SAM_LE, max_SAM_LE], axis = 0)
n_SAM_LE = np.asarray(data_LE[:,header == 'SAM_n'], float)[:,0]
SAM_q_LE = np.array(data_LE[:,header == 'SAM_q'])

S30_inj_LE = np.asarray(data_LE[:,header == '30S(Tg)'], float)
S15_inj_LE = np.asarray(data_LE[:,header == '15S(Tg)'], float)
N15_inj_LE = np.asarray(data_LE[:,header == '15N(Tg)'], float)
N30_inj_LE = np.asarray(data_LE[:,header == '30N(Tg)'], float)
total_inj_LE = S30_inj_LE + S15_inj_LE + N30_inj_LE + N15_inj_LE

T0_LE = np.asarray(data_LE[:,header == 'T0'], float)[:,0]
T1_LE = np.asarray(data_LE[:,header == 'T1'], float)[:,0]
T2_LE = np.asarray(data_LE[:,header == 'T2'], float)[:,0]

T0_error_LE = np.asarray(data_LE[:,header == 'T0_error'], float)[:,0]
T1_error_LE = np.asarray(data_LE[:,header == 'T1_error'], float)[:,0]
T2_error_LE = np.asarray(data_LE[:,header == 'T2_error'], float)[:,0]

VOLC_q_LE = np.array(data_LE[:,header == 'VOLC_q'])

##########################
######BASE INJ INFO#######
##########################
#Injection total when only the basestate is passed through the controller
BaseInj_2035 = np.array([n for i, n in enumerate(total_inj) if base_years[i] == 2035 and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
BaseInj_2045 = np.array([n for i, n in enumerate(total_inj) if base_years[i] == 2045 and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])

base_T0_error_2035 = np.array([n for i, n in enumerate(T0_error) if base_years[i] == 2035 and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
base_T1_error_2035 = np.array([n for i, n in enumerate(T1_error) if base_years[i] == 2035 and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
base_T2_error_2035 = np.array([n for i, n in enumerate(T2_error) if base_years[i] == 2035 and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])

base_T0_error_2045 = np.array([n for i, n in enumerate(T0_error) if base_years[i] == 2045 and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
base_T1_error_2045 = np.array([n for i, n in enumerate(T1_error) if base_years[i] == 2045 and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
base_T2_error_2045 = np.array([n for i, n in enumerate(T2_error) if base_years[i] == 2045 and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])

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
label_font = {'fontname':'Noto Sans JP', 'size':'8', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'}
#%%
####################################
###########GRIDDED PLOTS############
####################################

if COMPARE_NAO_ENSO:
    for y, year in enumerate(Base_years):
        total_inj_t = np.array([n for i, n in enumerate(total_inj) if base_years[i] == year and ENSO_q[i] == "ENSO" and NAO_q[i] == "NAO" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
        base_inj = np.array([n for i, n in enumerate(total_inj) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
        
        total_inj_t = (total_inj_t - base_inj)/base_inj * 100
        x_dim = (rangeV[:-1] + rangeV[1:])/2        
        total_inj_t_array = total_inj_t.reshape(len(x_dim), len(x_dim))
        
        title = "Impact of internal variability \non injection amounts in year " + str(year)
        pF.PlotColormesh(x_dim, x_dim, total_inj_t_array, "NAO state", "ENSO state", title) 
        if ADD_LE:
            total_inj_t = np.array([n for i, n in enumerate(total_inj_LE) if base_years_LE[i] == year and ENSO_q_LE[i] == "ENSO" and NAO_q_LE[i] == "NAO" and SAM_q_LE[i] == "NA" and VOLC_q_LE[i] == "NA"])
            
            total_inj_t = (total_inj_t - base_inj)/base_inj * 100
            x_dim = (rangeV[:-1] + rangeV[1:])/2        
            total_inj_t_array = total_inj_t.reshape(len(x_dim), len(x_dim))
            
            title = "Impact of internal variability \non injection amounts in year LE" + str(year)
            pF.PlotColormesh(x_dim, x_dim, total_inj_t_array, "NAO Index", "ENSO Index", title)   
            
    for y, year in enumerate(Base_years):
        total_inj_t = np.array([n for i, n in enumerate(T_Total_Error) if base_years[i] == year and ENSO_q[i] == "ENSO" and NAO_q[i] == "NAO" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
        base_inj = np.array([n for i, n in enumerate(T_Total_Error) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
    
        total_inj_t = total_inj_t - base_inj
        
        x_dim = (rangeV[:-1] + rangeV[1:])/2        
        total_inj_t_array = total_inj_t.reshape(len(x_dim), len(x_dim))
        
        title = "Impact of internal variability \non T0 in year " + str(year)
        pF.PlotColormeshSmall(x_dim, x_dim, total_inj_t_array, "NAO state", "ENSO state", title) 
        if ADD_LE:
            total_inj_t = np.array([n for i, n in enumerate(T0_error_LE) if base_years_LE[i] == year and ENSO_q_LE[i] == "ENSO" and NAO_q_LE[i] == "NAO" and SAM_q_LE[i] == "NA" and VOLC_q_LE[i] == "NA"])
            
            total_inj_t = (total_inj_t - base_inj)/base_inj * 100
            x_dim = (rangeV[:-1] + rangeV[1:])/2        
            total_inj_t_array = total_inj_t.reshape(len(x_dim), len(x_dim))
            
            title = "Impact of internal variability \non T0 in year LE" + str(year)
            pF.PlotColormesh(x_dim, x_dim, total_inj_t_array, "NAO Index", "ENSO Index", title)   

if COMPARE_NAO_SAM:
    for y, year in enumerate(Base_years):
        total_inj_t = np.array([n for i, n in enumerate(total_inj) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NAO" and SAM_q[i] == "SAM" and VOLC_q[i] == "NA"])
        base_inj = np.array([n for i, n in enumerate(total_inj) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
        
        total_inj_t = (total_inj_t - base_inj)/base_inj * 100
        x_dim = (rangeV[:-1] + rangeV[1:])/2        
        total_inj_t_array = total_inj_t.reshape(len(x_dim), len(x_dim))
        
        title = "Impact of internal variability \non injection amounts in year " + str(year)
        pF.PlotColormesh(x_dim, x_dim, total_inj_t_array, "SAM state", "NAO state", title)
        
        if ADD_LE:
            total_inj_t = np.array([n for i, n in enumerate(total_inj_LE) if base_years_LE[i] == year and ENSO_q_LE[i] == "NA" and NAO_q_LE[i] == "NAO" and SAM_q_LE[i] == "SAM" and VOLC_q_LE[i] == "NA"])

            total_inj_t = (total_inj_t - base_inj)/base_inj * 100
            x_dim = (rangeV[:-1] + rangeV[1:])/2        
            total_inj_t_array = total_inj_t.reshape(len(x_dim), len(x_dim))
            
            title = "Impact of internal variability \non injection amounts in year LE" + str(year)
            pF.PlotColormesh(x_dim, x_dim, total_inj_t_array, "SAM Index", "NAO Index", title)                 
  
if COMPARE_ENSO_SAM:
    for y, year in enumerate(Base_years):
        total_inj_t = np.array([n for i, n in enumerate(total_inj) if base_years[i] == year and ENSO_q[i] == "ENSO" and NAO_q[i] == "NA" and SAM_q[i] == "SAM" and VOLC_q[i] == "NA"])
        base_inj = np.array([n for i, n in enumerate(total_inj) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
        
        total_inj_t = (total_inj_t - base_inj)/base_inj * 100
        x_dim = (rangeV[:-1] + rangeV[1:])/2        
        total_inj_t_array = total_inj_t.reshape(len(x_dim), len(x_dim))
        
        title = "Impact of internal variability \non injection amounts in year " + str(year)
        pF.PlotColormesh(x_dim, x_dim, total_inj_t_array, "SAM state", "ENSO state", title)
        
        if ADD_LE:
            total_inj_t = np.array([n for i, n in enumerate(total_inj_LE) if base_years_LE[i] == year and ENSO_q_LE[i] == "ENSO" and NAO_q_LE[i] == "NA" and SAM_q_LE[i] == "SAM" and VOLC_q_LE[i] == "NA"])
            base_inj = np.array([n for i, n in enumerate(total_inj_LE) if base_years_LE[i] == year and ENSO_q_LE[i] == "NA" and NAO_q_LE[i] == "NA" and SAM_q_LE[i] == "NA" and VOLC_q_LE[i] == "NA"])
            
            total_inj_t = (total_inj_t - base_inj)/base_inj * 100
            x_dim = (rangeV[:-1] + rangeV[1:])/2        
            total_inj_t_array = total_inj_t.reshape(len(x_dim), len(x_dim))
            
            title = "Impact of internal variability \non injection amounts in year LE" + str(year)
            pF.PlotColormesh(x_dim, x_dim, total_inj_t_array, "SAM Index", "ENSO Index", title)    


####################################
############Explore Delta###########
####################################
SUBTRACT_MODES = False
if SUBTRACT_MODES:
    total_inj_ENSO2035 = np.array([n for i, n in enumerate(total_inj) if base_years[i] == 2035 and ENSO_q[i] == "ENSO" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
    total_inj_ENSO2045 = np.array([n for i, n in enumerate(total_inj) if base_years[i] == 2045 and ENSO_q[i] == "ENSO" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
    avg_ENSO_t = np.array([n for i, n in enumerate(avg_ENSO) if base_years[i] == 2035 and ENSO_q[i] == "ENSO" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"]) 

    inj_ENSO_35 = total_inj_ENSO2035 - BaseInj_2035
    inj_ENSO_45 = total_inj_ENSO2045 - BaseInj_2045
    
    lim = [-.5, .5]
    plt.figure(dpi = (300), figsize = (8, 3))
    pF.PlotLines(avg_ENSO_t, inj_ENSO_35, "ENSO", "Injection (Tg/year)", "ENSO and DELTA", None, "2035", '-', 1, lim)
    pF.PlotLines(avg_ENSO_t, inj_ENSO_45, "ENSO", "Injection (Tg/year)", "ENSO and DELTA", None, "2045", '-', 1, lim)
    plt.legend(loc = 'lower right', frameon=False, prop=legend_font)
    
    delta_ENSO = inj_ENSO_45 - inj_ENSO_35
    lim = [-1, 1]
    plt.figure(dpi = (300), figsize = (8, 3))
    pF.PlotLines(avg_ENSO_t, delta_ENSO, "ENSO", "Injection (Tg/year)", "DELTA", None, "2035", '-', 1, lim)
    
    ##########################
    
    total_inj_NAO2035 = np.array([n for i, n in enumerate(total_inj) if base_years[i] == 2035 and ENSO_q[i] == "NA" and NAO_q[i] == "NAO" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
    total_inj_NAO2045 = np.array([n for i, n in enumerate(total_inj) if base_years[i] == 2045 and ENSO_q[i] == "NA" and NAO_q[i] == "NAO" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
    avg_NAO_t = np.array([n for i, n in enumerate(avg_NAO) if base_years[i] == 2035 and ENSO_q[i] == "NA" and NAO_q[i] == "NAO" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"]) 

    inj_NAO_35 = total_inj_NAO2035 - BaseInj_2035
    inj_NAO_45 = total_inj_NAO2045 - BaseInj_2045
    
    lim = [-.5, .5]
    plt.figure(dpi = (300), figsize = (8, 3))
    pF.PlotLines(avg_NAO_t, inj_NAO_35, "NAO", "Injection (Tg/year)", "NAO and DELTA", None, "2035", '-', 1, lim)
    pF.PlotLines(avg_NAO_t, inj_NAO_45, "NAO", "Injection (Tg/year)", "NAO and DELTA", None, "2045", '-', 1, lim)
    plt.legend(loc = 'lower right', frameon=False, prop=legend_font)
    
    delta_NAO = inj_NAO_45 - inj_NAO_35
    lim = [-1, 1]
    plt.figure(dpi = (300), figsize = (8, 3))
    pF.PlotLines(avg_NAO_t, delta_NAO, "NAO", "Injection (Tg/year)", "DELTA", None, "2035", '-', 1, lim)
    
    ##########################
    
    total_inj_SAM2035 = np.array([n for i, n in enumerate(total_inj) if base_years[i] == 2035 and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "SAM" and VOLC_q[i] == "NA"])
    total_inj_SAM2045 = np.array([n for i, n in enumerate(total_inj) if base_years[i] == 2045 and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "SAM" and VOLC_q[i] == "NA"])
    avg_SAM_t = np.array([n for i, n in enumerate(avg_SAM) if base_years[i] == 2035 and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "SAM" and VOLC_q[i] == "NA"]) 

    inj_SAM_35 = total_inj_SAM2035 - BaseInj_2035
    inj_SAM_45 = total_inj_SAM2045 - BaseInj_2045
    
    lim = [-.5, .5]
    plt.figure(dpi = (300), figsize = (8, 3))
    pF.PlotLines(avg_SAM_t, inj_SAM_35, "SAM", "Injection (Tg/year)", "SAM and DELTA", None, "2035", '-', 1, lim)
    pF.PlotLines(avg_SAM_t, inj_SAM_45, "SAM", "Injection (Tg/year)", "SAM and DELTA", None, "2045", '-', 1, lim)
    plt.legend(loc = 'lower right', frameon=False, prop=legend_font)
    
    delta_SAM = inj_SAM_45 - inj_SAM_35
    lim = [-1, 1]
    plt.figure(dpi = (300), figsize = (8, 3))
    pF.PlotLines(avg_SAM_t, delta_SAM, "SAM", "Injection (Tg/year)", "DELTA", None, "2035", '-', 1, lim)

####################################
#############LINE PLOTS#############
####################################
if PLOT_ENSO_BASESTATE_TOT_N:
    lim = [-100, 100]
    plt.figure(dpi = (300), figsize = (8, 3))
    for y, year in enumerate(Base_years):
        total_inj_t = np.array([n for i, n in enumerate(total_inj) if base_years[i] == year and ENSO_q[i] == "ENSO" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
        avg_ENSO_t = np.array([n for i, n in enumerate(avg_ENSO) if base_years[i] == year and ENSO_q[i] == "ENSO" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
        n_ENSO_t = np.array([n for i, n in enumerate(n_ENSO) if base_years[i] == year and ENSO_q[i] == "ENSO" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])

        if year == 2035:
            base_inj = BaseInj_2035
        else:
            base_inj = BaseInj_2045
            
        per_change = ((total_inj_t - base_inj)/base_inj) * 100
        pF.PlotLines(avg_ENSO_t, per_change, "ENSO", "Percent Change", "Percent Change", Colors[y], year, '-', 1, lim)
           
    plt.legend(loc = 'lower right', frameon=False, prop=legend_font)
    ##################
    if ADD_N_LABEL:
        for i, N in enumerate(n_ENSO_t):
            plt.text(avg_ENSO_t[i], lim[1] + .25, str(int(n_ENSO_t[i])), horizontalalignment='center', **label_font)
    # ##################
    
    if ADD_LE:
        for y, year in enumerate(Base_years):
            if year == 2035:
                base_inj = BaseInj_2035
            else:
                base_inj = BaseInj_2045
            
            total_inj_t = np.array([n for i, n in enumerate(total_inj_LE) if base_years_LE[i] == year and ENSO_q_LE[i] == "ENSO" and NAO_q_LE[i] == "NA" and SAM_q_LE[i] == "NA" and VOLC_q_LE[i] == "NA"])
            avg_ENSO_t = np.array([n for i, n in enumerate(avg_ENSO_LE) if base_years_LE[i] == year and ENSO_q_LE[i] == "ENSO" and NAO_q_LE[i] == "NA" and SAM_q_LE[i] == "NA" and VOLC_q_LE[i] == "NA"])

            per_change = ((total_inj_t - base_inj)/base_inj) * 100
            pF.PlotLines(avg_ENSO_t, per_change, "ENSO", "Percent Change", "Percent Change", Colors[y], year, '--', .5, lim)
            
    plt.show()
    
    lim = [-.4, .4]
    if PLOT_T_ERROR:
        for y, year in enumerate(Base_years):
            
            if year == 2035:
                base_T0_error = base_T0_error_2035
                base_T1_error = base_T1_error_2035
                base_T2_error = base_T2_error_2035
            else:
                base_T0_error = base_T0_error_2045
                base_T1_error = base_T1_error_2045
                base_T2_error = base_T2_error_2045
            
            plt.figure(dpi = (300), figsize = (8, 3))
            T0_error_t = np.array([n for i, n in enumerate(T0_error) if base_years[i] == year and ENSO_q[i] == "ENSO" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
            T1_error_t = np.array([n for i, n in enumerate(T1_error) if base_years[i] == year and ENSO_q[i] == "ENSO" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
            T2_error_t = np.array([n for i, n in enumerate(T2_error) if base_years[i] == year and ENSO_q[i] == "ENSO" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
            avg_ENSO_t = np.array([n for i, n in enumerate(avg_ENSO) if base_years[i] == year and ENSO_q[i] == "ENSO" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])

            pF.PlotLines(avg_ENSO_t, (T0_error_t - base_T0_error), "ENSO", "Error", year, None, "T0", '-', 1, lim)
            pF.PlotLines(avg_ENSO_t, (T1_error_t - base_T1_error), "ENSO", "Error", year, None, "T1", '-', 1, lim)
            pF.PlotLines(avg_ENSO_t, (T2_error_t - base_T2_error), "ENSO", "Error", year, None, "T2", '-', 1, lim)

            plt.legend(loc = 'lower right', frameon=False, prop=legend_font)
            plt.show()   
    
if PLOT_NAO_BASESTATE_TOT_N:
    lim = [-50, 50]
    plt.figure(dpi = (300), figsize = (8, 3))
    for y, year in enumerate(Base_years):
        total_inj_t = np.array([n for i, n in enumerate(total_inj) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NAO" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
        avg_NAO_t = np.array([n for i, n in enumerate(avg_NAO) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NAO" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
        n_NAO_t = np.array([n for i, n in enumerate(n_NAO) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NAO" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
        
        if year == 2035:
            base_inj = BaseInj_2035
        else:
            base_inj = BaseInj_2045
        
        per_change = ((total_inj_t - base_inj)/base_inj) * 100
        pF.PlotLines(avg_NAO_t, per_change, "NAO", "Percent Change", "Percent Change", Colors[y], year, '-', 1, lim)
           
    plt.legend(loc = 'lower left', frameon=False, prop=legend_font)
    ##################
    if ADD_N_LABEL:
        for i, N in enumerate(n_NAO_t):
            plt.text(avg_NAO_t[i], lim[1] + .25, str(int(n_NAO_t[i])), horizontalalignment='center', **label_font)
    ##################
    
    if ADD_LE:
        print("ADD LE NAO")
        for y, year in enumerate(Base_years):
            total_inj_t = np.array([n for i, n in enumerate(total_inj_LE) if base_years_LE[i] == year and ENSO_q_LE[i] == "NA" and NAO_q_LE[i] == "NAO" and SAM_q_LE[i] == "NA" and VOLC_q_LE[i] == "NA"])
            avg_NAO_t = np.array([n for i, n in enumerate(avg_NAO_LE) if base_years_LE[i] == year and ENSO_q_LE[i] == "NA" and NAO_q_LE[i] == "NAO" and SAM_q_LE[i] == "NA" and VOLC_q_LE[i] == "NA"])
            
            if year == 2035:
                base_inj = BaseInj_2035
            else:
                base_inj = BaseInj_2045

            per_change = ((total_inj_t - base_inj)/base_inj) * 100
            pF.PlotLines(avg_NAO_t, per_change, "NAO", "Percent Change", "Percent Change", Colors[y], year, '--', .5, lim)
            
    plt.show()
    lim = [-.1, .1]
    if PLOT_T_ERROR:
        for y, year in enumerate(Base_years):
            plt.figure(dpi = (300), figsize = (8, 3))
            T0_error_t = np.array([n for i, n in enumerate(T0_error) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NAO" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
            T1_error_t = np.array([n for i, n in enumerate(T1_error) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NAO" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
            T2_error_t = np.array([n for i, n in enumerate(T2_error) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NAO" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
            avg_NAO_t = np.array([n for i, n in enumerate(avg_NAO) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NAO" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])

            if year == 2035:
                base_T0_error = base_T0_error_2035
                base_T1_error = base_T1_error_2035
                base_T2_error = base_T2_error_2035
            else:
                base_T0_error = base_T0_error_2045
                base_T1_error = base_T1_error_2045
                base_T2_error = base_T2_error_2045

            pF.PlotLines(avg_NAO_t, (T0_error_t - base_T0_error), "NAO", "Error", "", None, "T0", '-', 1, lim)
            pF.PlotLines(avg_NAO_t, (T1_error_t - base_T1_error), "NAO", "Error", "", None, "T1", '-', 1, lim)
            pF.PlotLines(avg_NAO_t, (T2_error_t - base_T2_error), "NAO", "Error", "", None, "T2", '-', 1, lim)
                   
            plt.legend(loc = 'lower left', frameon=False, prop=legend_font)
            plt.show()
    
if PLOT_SAM_BASESTATE_TOT_N:
    lim = [-50, 50]
    plt.figure(dpi = (300), figsize = (8, 3))
    for y, year in enumerate(Base_years):
        total_inj_t = np.array([n for i, n in enumerate(total_inj) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "SAM" and VOLC_q[i] == "NA"])
        avg_SAM_t = np.array([n for i, n in enumerate(avg_SAM) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "SAM" and VOLC_q[i] == "NA"])
        n_SAM_t = np.array([n for i, n in enumerate(n_SAM) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "SAM" and VOLC_q[i] == "NA"])
        
        if year == 2035:
            base_inj = BaseInj_2035
        else:
            base_inj = BaseInj_2045
        
        per_change = ((total_inj_t - base_inj)/base_inj) * 100
        pF.PlotLines(avg_SAM_t, per_change, "SAM", "Percent Change", "Percent Change", Colors[y], year, '-', 1, lim)
           
    plt.legend(loc = 'lower left', frameon=False, prop=legend_font)
    ##################
    if ADD_N_LABEL:
        for i, N in enumerate(n_SAM_t):
            plt.text(avg_SAM_t[i], lim[1] + .25, str(int(n_SAM_t[i])), horizontalalignment='center', **label_font)
    ##################
    
    if ADD_LE:
        print("LE")
        for y, year in enumerate(Base_years):
            total_inj_t = np.array([n for i, n in enumerate(total_inj_LE) if base_years_LE[i] == year and ENSO_q_LE[i] == "NA" and NAO_q_LE[i] == "NA" and SAM_q_LE[i] == "SAM" and VOLC_q_LE[i] == "NA"])
            avg_SAM_t = np.array([n for i, n in enumerate(avg_SAM_LE) if base_years_LE[i] == year and ENSO_q_LE[i] == "NA" and NAO_q_LE[i] == "NA" and SAM_q_LE[i] == "SAM" and VOLC_q_LE[i] == "NA"])
            
            if year == 2035:
                base_inj = BaseInj_2035
            else:
                base_inj = BaseInj_2045


            per_change = ((total_inj_t - base_inj)/base_inj) * 100
            pF.PlotLines(avg_SAM_t, per_change, "SAM", "Percent Change", "Percent Change", Colors[y], year, '--', .5, lim)

    plt.show()
    lim = [-.5, .5]
    plt.figure(dpi = (300), figsize = (8, 3))
    if PLOT_T_ERROR:
        for y, year in enumerate(Base_years):
            plt.figure(dpi = (300), figsize = (8, 3))
            T0_error_t = np.array([n for i, n in enumerate(T0_error) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "SAM" and VOLC_q[i] == "NA"])
            T1_error_t = np.array([n for i, n in enumerate(T1_error) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "SAM" and VOLC_q[i] == "NA"])
            T2_error_t = np.array([n for i, n in enumerate(T2_error) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "SAM" and VOLC_q[i] == "NA"])
            avg_SAM_t = np.array([n for i, n in enumerate(avg_SAM) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "SAM" and VOLC_q[i] == "NA"])

            if year == 2035:
                base_T0_error = base_T0_error_2035
                base_T1_error = base_T1_error_2035
                base_T2_error = base_T2_error_2035
            else:
                base_T0_error = base_T0_error_2035
                base_T1_error = base_T1_error_2035
                base_T2_error = base_T2_error_2035

            pF.PlotLines(avg_SAM_t, (T0_error_t - base_T0_error), "SAM", "Error", year, None, "T0", '-', 1, lim)
            pF.PlotLines(avg_SAM_t, (T1_error_t - base_T1_error), "SAM", "Error", year, None, "T1", '-', 1, lim)
            pF.PlotLines(avg_SAM_t, (T2_error_t - base_T2_error), "SAM", "Error", year, None, "T2", '-', 1, lim)
               
            plt.legend(loc = 'lower left', frameon=False, prop=legend_font)
            plt.show()

if PLOT_ENSO_BASESTATE_TOT_N_VOLC:
    lim = [-120, 20]
    plt.figure(dpi = (300), figsize = (8, 3))
    for y, year in enumerate(Base_years):
        total_inj_t = np.array([n for i, n in enumerate(total_inj) if base_years[i] == year and ENSO_q[i] == "ENSO" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "VOLC"])
        avg_ENSO_t = np.array([n for i, n in enumerate(avg_ENSO) if base_years[i] == year and ENSO_q[i] == "ENSO" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "VOLC"])
        base_inj = np.array([n for i, n in enumerate(total_inj) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
        n_ENSO_t = np.array([n for i, n in enumerate(n_ENSO) if base_years[i] == year and ENSO_q[i] == "ENSO" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "VOLC"])
        
        per_change = ((total_inj_t - base_inj)/base_inj) * 100
        pF.PlotLines(avg_ENSO_t, per_change, "ENSO and Pinatubo", "Percent Change", "Percent Change", Colors[y], year, '-', 1, lim)
           
    plt.legend(loc = 'lower right', frameon=False, prop=legend_font)
    ##################
    if ADD_N_LABEL:
        for i, N in enumerate(n_ENSO_t):
            plt.text(avg_ENSO_t[i], lim[1] + .25, str(int(n_ENSO_t[i])), horizontalalignment='center', **label_font)
    # ##################
    
    if ADD_LE:
        for y, year in enumerate(Base_years):
            total_inj_t = np.array([n for i, n in enumerate(total_inj_LE) if base_years_LE[i] == year and ENSO_q_LE[i] == "ENSO" and NAO_q_LE[i] == "NA" and SAM_q_LE[i] == "NA" and VOLC_q_LE[i] == "VOLC"])
            avg_ENSO_t = np.array([n for i, n in enumerate(avg_ENSO_LE) if base_years_LE[i] == year and ENSO_q_LE[i] == "ENSO" and NAO_q_LE[i] == "NA" and SAM_q_LE[i] == "NA" and VOLC_q_LE[i] == "VOLC"])
            base_inj = np.array([n for i, n in enumerate(total_inj_LE) if base_years_LE[i] == year and ENSO_q_LE[i] == "NA" and NAO_q_LE[i] == "NA" and SAM_q_LE[i] == "NA" and VOLC_q_LE[i] == "NA"])
            
            per_change = ((total_inj_t - base_inj)/base_inj) * 100
            pF.PlotLines(avg_ENSO_t, per_change, "ENSO and Pinatubo", "Percent Change", "Percent Change", Colors[y], year, '--', .5, lim)
            
    plt.show()

if PLOT_NAO_BASESTATE_TOT_N_VOLC:
    lim = [-120, 20]
    plt.figure(dpi = (300), figsize = (8, 3))
    for y, year in enumerate(Base_years):
        total_inj_t = np.array([n for i, n in enumerate(total_inj) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NAO" and SAM_q[i] == "NA" and VOLC_q[i] == "VOLC"])
        avg_NAO_t = np.array([n for i, n in enumerate(avg_NAO) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NAO" and SAM_q[i] == "NA" and VOLC_q[i] == "VOLC"])
        base_inj = np.array([n for i, n in enumerate(total_inj) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
        n_NAO_t = np.array([n for i, n in enumerate(n_NAO) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NAO" and SAM_q[i] == "NA" and VOLC_q[i] == "VOLC"])
        
        per_change = ((total_inj_t - base_inj)/base_inj) * 100
        pF.PlotLines(avg_NAO_t, per_change, "NAO and Pinatubo", "Percent Change", "Percent Change", Colors[y], year, '-', 1, lim)
           
    # plt.legend(loc = 'lower left', frameon=False, prop=legend_font)
    ##################
    if ADD_N_LABEL:
        for i, N in enumerate(n_NAO_t):
            plt.text(avg_NAO_t[i], lim[1] + .25, str(int(n_NAO_t[i])), horizontalalignment='center', **label_font)
    ##################
    
    if ADD_LE:
        print("ADD LE NAO")
        for y, year in enumerate(Base_years):
            total_inj_t = np.array([n for i, n in enumerate(total_inj_LE) if base_years_LE[i] == year and ENSO_q_LE[i] == "NA" and NAO_q_LE[i] == "NAO" and SAM_q_LE[i] == "NA" and VOLC_q_LE[i] == "VOLC"])
            avg_NAO_t = np.array([n for i, n in enumerate(avg_NAO_LE) if base_years_LE[i] == year and ENSO_q_LE[i] == "NA" and NAO_q_LE[i] == "NAO" and SAM_q_LE[i] == "NA" and VOLC_q_LE[i] == "VOLC"])
            base_inj = np.array([n for i, n in enumerate(total_inj_LE) if base_years_LE[i] == year and ENSO_q_LE[i] == "NA" and NAO_q_LE[i] == "NA" and SAM_q_LE[i] == "NA" and VOLC_q_LE[i] == "NA"])

            per_change = ((total_inj_t - base_inj)/base_inj) * 100
            pF.PlotLines(avg_NAO_t, per_change, "NAO and Pinatubo", "Percent Change", "Percent Change", Colors[y], year, '--', .5, lim)
            
    plt.show()
    
if PLOT_SAM_BASESTATE_TOT_N_VOLC:
    lim = [-120, 20]
    plt.figure(dpi = (300), figsize = (8, 3))
    for y, year in enumerate(Base_years):
        total_inj_t = np.array([n for i, n in enumerate(total_inj) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "SAM" and VOLC_q[i] == "VOLC"])
        avg_SAM_t = np.array([n for i, n in enumerate(avg_SAM) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "SAM" and VOLC_q[i] == "VOLC"])
        base_inj = np.array([n for i, n in enumerate(total_inj) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
        n_SAM_t = np.array([n for i, n in enumerate(n_SAM) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "SAM" and VOLC_q[i] == "VOLC"])
        
        per_change = ((total_inj_t - base_inj)/base_inj) * 100
        pF.PlotLines(avg_SAM_t, per_change, "SAM and Pinatubo", "Percent Change", "Percent Change", Colors[y], year, '-', 1, lim)
           
    # plt.legend(loc = 'lower left', frameon=False, prop=legend_font)
    ##################
    if ADD_N_LABEL:
        for i, N in enumerate(n_SAM_t):
            plt.text(avg_SAM_t[i], lim[1] + .25, str(int(n_SAM_t[i])), horizontalalignment='center', **label_font)
    ##################
    
    if ADD_LE:
        print("LE")
        for y, year in enumerate(Base_years):
            total_inj_t = np.array([n for i, n in enumerate(total_inj_LE) if base_years_LE[i] == year and ENSO_q_LE[i] == "NA" and NAO_q_LE[i] == "NA" and SAM_q_LE[i] == "SAM" and VOLC_q_LE[i] == "VOLC"])
            avg_SAM_t = np.array([n for i, n in enumerate(avg_SAM_LE) if base_years_LE[i] == year and ENSO_q_LE[i] == "NA" and NAO_q_LE[i] == "NA" and SAM_q_LE[i] == "SAM" and VOLC_q_LE[i] == "VOLC"])
            base_inj = np.array([n for i, n in enumerate(total_inj_LE) if base_years_LE[i] == year and ENSO_q_LE[i] == "NA" and NAO_q_LE[i] == "NA" and SAM_q_LE[i] == "NA" and VOLC_q_LE[i] == "NA"])

            per_change = ((total_inj_t - base_inj)/base_inj) * 100
            pF.PlotLines(avg_SAM_t, per_change, "SAM and Pinatubo", "Percent Change", "Percent Change", Colors[y], year, '--', .5, lim)
            
    plt.show()

VOLC_IMPACT = False
if VOLC_IMPACT:
    year = 2045
    lim = [-120, 120]
    base_inj = np.array([n for i, n in enumerate(total_inj) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "NA" and VOLC_q[i] == "NA"])
    total_inj_ENSO_VOLC = np.array([n for i, n in enumerate(total_inj) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "SAM" and VOLC_q[i] == "VOLC"])
    total_inj_ENSO = np.array([n for i, n in enumerate(total_inj) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "SAM" and VOLC_q[i] == "NA"])
    avg_ENSO_t = np.array([n for i, n in enumerate(avg_SAM) if base_years[i] == year and ENSO_q[i] == "NA" and NAO_q[i] == "NA" and SAM_q[i] == "SAM" and VOLC_q[i] == "NA"])
    
    total_inj_ENSO_VOLC = ((total_inj_ENSO_VOLC - base_inj)/base_inj) * 100
    total_inj_ENSO = ((total_inj_ENSO - base_inj)/base_inj) * 100
    
    pF.PlotLines(avg_ENSO_t, total_inj_ENSO_VOLC, "SAM and Pinatubo", "Percent Change", "Percent Change", Colors[0], "ENSO and VOLC", '-', 1, lim)
    pF.PlotLines(avg_ENSO_t, total_inj_ENSO, "SAM and Pinatubo", "Percent Change", "Percent Change", Colors[1], "ENSO", '-', 1, lim)
    pF.PlotLines(avg_ENSO_t, (total_inj_ENSO - total_inj_ENSO_VOLC), "SAM and Pinatubo", "Percent Change", "Percent Change", Colors[2], "Diff", '-', 1, lim)


