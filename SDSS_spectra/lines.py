import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from _global_ import *
from scipy.optimize import curve_fit

def func(x, a, b):
    return a*x + b

path_in = r'E:\backup\ALMA9\SDSS_GAMA.csv'

df = pd.read_csv(path_in)

fig, axs = plt.subplots(1, 2, figsize=(12, 6), tight_layout=True)

for i, item in enumerate(df['HgF_cont_g']):
    axs[0].errorbar(item, df['Hg_cont_SDSS'][i], yerr=df['Hg_cont_err_SDSS'][i], color=cd_WHAN[df['WHAN'][i]][0], alpha = 0.4, marker='.')
    axs[1].errorbar(item, df['Hg_cont_SDSS'][i], yerr=df['Hg_cont_err_SDSS'][i], color=color_dict[df['BPT'][i]][0], alpha = 0.4, marker='.')

    axs[0].set_ylabel(r'$H\gamma cont. SDSS$')
    axs[1].set_xlabel(r'$H\gamma cont. GAMA$')
    axs[0].set_xlabel(r'$H\gamma cont. GAMA$')
    
gama_x = np.array(df['HgF_cont_g'])
sdss_y = np.array(df['Hg_cont_SDSS'])

popt, pcov = curve_fit(func, gama_x, sdss_y, sigma = df['Hg_cont_err_SDSS'])

axs[0].plot(gama_x, popt[0]*gama_x + popt[1])
axs[1].plot(gama_x, popt[0]*gama_x + popt[1])

print(1/popt[0], popt[1])
plt.savefig('./FIGURES/HG_cont.pdf')

fig, axs = plt.subplots(1, 2, figsize=(12, 6), tight_layout=True)

for i, item in enumerate(df['HA_cont_g']):
    axs[0].errorbar(item, df['HA_cont_SDSS'][i], yerr=df['HA_cont_err_SDSS'][i], color=cd_WHAN[df['WHAN'][i]][0], alpha = 0.4, marker='.')
    axs[1].errorbar(item, df['HA_cont_SDSS'][i], yerr=df['HA_cont_err_SDSS'][i], color=color_dict[df['BPT'][i]][0], alpha = 0.4, marker='.')

    axs[0].set_ylabel(r'$H\alpha cont. SDSS$')
    axs[1].set_xlabel(r'$H\alpha cont. GAMA$')
    axs[0].set_xlabel(r'$H\alpha cont. GAMA$')
    
gama_x = np.array(df['HA_cont_g'])
sdss_y = np.array(df['HA_cont_SDSS'])
popt, pcov = curve_fit(func, gama_x, sdss_y, sigma=df['HA_cont_err_SDSS'])

axs[0].plot(gama_x, popt[0]*gama_x + popt[1])
axs[1].plot(gama_x, popt[0]*gama_x + popt[1])

print(1/popt[0], popt[1])
    
plt.savefig('./FIGURES/HA_cont.pdf')

fig, axs = plt.subplots(1, 2, figsize=(12, 6), tight_layout=True)

for i, item in enumerate(df['HA_g']):
    axs[0].scatter(item/df['HA_er_g'][i], df['HA_SDSS'][i]/df['HA_er_SDSS'][i], color=cd_WHAN[df['WHAN'][i]][0], alpha = 0.4, marker='.')
    axs[1].scatter(item/df['HA_er_g'][i], df['HA_SDSS'][i]/df['HA_er_SDSS'][i], color=color_dict[df['BPT'][i]][0], alpha = 0.4, marker='.')

    axs[0].set_ylabel(r'$SN H\alpha SDSS$')
    axs[1].set_xlabel(r'$SN H\alpha GAMA$')
    axs[0].set_xlabel(r'$SN H\alpha GAMA$')
    
gama_x = np.array(df['HA_g']/df['HA_er_g'])
sdss_y = np.array(df['HA_SDSS']/df['HA_er_SDSS'])
popt, pcov = curve_fit(func, gama_x, sdss_y)

axs[0].plot(gama_x, popt[0]*gama_x + popt[1])
axs[1].plot(gama_x, popt[0]*gama_x + popt[1])

axs[0].axhline(y = 0, linestyle='--')
axs[0].axvline(x = 0, linestyle='--')
axs[1].axhline(y = 0, linestyle='--')
axs[1].axvline(x = 0, linestyle='--')

print(popt[0], popt[1])
    
plt.savefig('./FIGURES/SN_HA.pdf')





fig, axs = plt.subplots(1, 2, figsize=(12, 6), tight_layout=True)

for i, item in enumerate(df['HgF_g']):
    axs[0].scatter(item/df['HgF_er_g'][i], df['Hg_SDSS'][i]/df['Hg_er_SDSS'][i], color=cd_WHAN[df['WHAN'][i]][0], alpha = 0.4, marker='.')
    axs[1].scatter(item/df['HgF_er_g'][i], df['Hg_SDSS'][i]/df['Hg_er_SDSS'][i], color=color_dict[df['BPT'][i]][0], alpha = 0.4, marker='.')

    axs[0].set_ylabel(r'$SN H\gamma SDSS$')
    axs[1].set_xlabel(r'$SN H\gamma GAMA$')
    axs[0].set_xlabel(r'$SN H\gamma GAMA$')
    
gama_x = np.array(df['HgF_g']/df['HgF_er_g'])
sdss_y = np.array(df['Hg_SDSS']/df['Hg_er_SDSS'])
popt, pcov = curve_fit(func, gama_x, sdss_y)

axs[0].plot(gama_x, popt[0]*gama_x + popt[1])
axs[1].plot(gama_x, popt[0]*gama_x + popt[1])

axs[0].axhline(y = 0, linestyle='--')
axs[0].axvline(x = 0, linestyle='--')
axs[1].axhline(y = 0, linestyle='--')
axs[1].axvline(x = 0, linestyle='--')

print(popt[0], popt[1])
    
plt.savefig('./FIGURES/SN_HG.pdf')