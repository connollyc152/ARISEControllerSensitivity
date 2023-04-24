import numpy as np
import matplotlib.pyplot as plt

datapath = "/Users/cconn/Documents/Explore_controller/controller_output/data/"
datafile = "ControlLog_EXPLORING_ENSO_STATE.txt"

data = np.array([x.split(' ') for x in open (datapath + datafile).readlines()])[:,:-1]

header = data[0,:]
data = data[1:,:]

PLOT_ENSO_BASESTATE_TOT = True

base_years = np.asarray(data[:,header == 'BASESTATE_YEAR'], float)[:,0]

min_ENSO = np.asarray(data[:,header == 'ENSO_MIN'], float)[:,0]
max_ENSO = np.asarray(data[:,header == 'ENSO_MAX'], float)[:,0]
avg_ENSO = np.mean([min_ENSO, max_ENSO], axis = 0)

S30_inj = np.asarray(data[:,header == '30S(Tg)'], float)
S15_inj = np.asarray(data[:,header == '15S(Tg)'], float)
N15_inj = np.asarray(data[:,header == '15N(Tg)'], float)
N30_inj = np.asarray(data[:,header == '30N(Tg)'], float)
total_inj = np.add(S30_inj, S15_inj)
total_inj = S30_inj + S15_inj + N30_inj + N15_inj

if PLOT_ENSO_BASESTATE_TOT:
    Base_years = [2060, 2045, 2035]#, 2020]
    Colors = ['#665A73', '#A6A26F', '#A33E5E']
    plt.figure(dpi = (300), figsize = (8, 3))
    for y, year in enumerate(Base_years):
        avg_ENSO_t = avg_ENSO[base_years == year]
        total_inj_t = total_inj[base_years == year]
        if len(total_inj_t) == 0:
            continue
        plt.plot(avg_ENSO_t, total_inj_t, color = Colors[y], linewidth=3)
        plt.xlabel("ENSO")
        plt.ylabel("Total Injection (Tg)")
        plt.xlim(min(avg_ENSO), max(avg_ENSO))
    plt.show()

