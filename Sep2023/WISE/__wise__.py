import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def adjusting_plotting_pars():
    plt.rcParams['font.size'] = 20
    mpl.rcParams['xtick.labelsize'] = 20
    mpl.rcParams['ytick.labelsize'] = 20 
    mpl.rcParams['axes.labelsize'] = 20
    mpl.rcParams['ytick.major.size'] = 15
    mpl.rcParams['xtick.major.size'] = 15
    # print(mpl.rcParams.keys())

def flux_to_color(s_w1, s_w2):
    return 2.5*np.log10(s_w2/s_w1) #w1 - w2

def flux_to_mag(s_w, zero_mag):
    return 2.5*np.log10(zero_mag/s_w)

def my_autopct(pct):
    return (f'{pct:.2f}%') if pct > 5 else ''