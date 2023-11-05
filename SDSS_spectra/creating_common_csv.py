from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as plt
from _global_ import *
import numpy as np
from scipy.optimize import curve_fit

def func(x, a, b):
    return a*x + b

fits_db_path = "E:\databases\galSpecLine-dr8.fits"

hdul = fits.open(fits_db_path)
cols = hdul[1].columns
# print(cols.info())

dataframe_sdss = hdul[1].data

#gama_path = r"E:\backup\backup_BPT\GAMA_ETG_OLA_R.csv"
gama_path = r"E:/databases/GAMA_ETG_OLA_SDSS_R.csv"

dataframe_gama = pd.read_csv(filepath_or_buffer=gama_path, usecols=['SPEC_ID', 'SURVEY', 'BPT', 'WHAN', 'HgF', 'HgF_cont', 'HgF_er', 'HA', 'HA_cont', 'HA_er', 'GALRE_i', 'GALREERR_i'])
headers = []
with open(r'E:\backup\backup_BPT\SDSS_spectra\cols.txt', 'r') as f:
    lines = f.readlines()
    lines_stripped = [line.strip() for line in lines]
    for line in lines_stripped:
        headers.append(line.split()[0])

keys = []
for item in dataframe_gama['SPEC_ID']:
    keys.append(str(item))

hgamma = []
hgamma_cont = []
hgamma_err = []
hgamma_cont_err = []
ha = []
ha_cont = []
ha_err = []
ha_cont_err = []
spec_id = []

for i, item in enumerate(dataframe_sdss['SPECOBJID']):
    if dataframe_sdss['SPECOBJID'][i] != '':
        if str(dataframe_sdss['SPECOBJID'][i]) in keys:
            spec_id.append(dataframe_sdss['SPECOBJID'][i])
            hgamma.append(dataframe_sdss['H_GAMMA_FLUX'][i])
            hgamma_err.append(dataframe_sdss['H_GAMMA_FLUX_ERR'][i])
            hgamma_cont.append(dataframe_sdss['H_GAMMA_CONT'][i])
            hgamma_cont_err.append(dataframe_sdss['H_GAMMA_CONT_ERR'][i])
            ha.append(dataframe_sdss['H_ALPHA_FLUX'][i])
            ha_err.append(dataframe_sdss['H_ALPHA_FLUX_ERR'][i])
            ha_cont.append(dataframe_sdss['H_ALPHA_CONT'][i])
            ha_cont_err.append(dataframe_sdss['H_ALPHA_CONT_ERR'][i])

print(len(spec_id))
DICT_SDSS = {
    'SPEC_ID' : spec_id,
    'Hg_SDSS' : hgamma,
    'Hg_er_SDSS' : hgamma_err,
    'Hg_cont_SDSS' : hgamma_cont,
    'Hg_cont_err_SDSS' : hgamma_cont_err,
    'HA_SDSS' : ha,
    'HA_er_SDSS' : ha_err,
    'HA_cont_SDSS' : ha_cont,
    'HA_cont_err_SDSS' : ha_cont_err
}

dataframe_sdss = pd.DataFrame.from_dict(DICT_SDSS)

dataframe_sdss['SPEC_ID']=dataframe_sdss['SPEC_ID'].astype(str)
dataframe_gama['SPEC_ID']=dataframe_gama['SPEC_ID'].astype(str)

df = pd.merge(dataframe_gama, dataframe_sdss, how='inner', on='SPEC_ID')
df.to_csv('E:/databases/SDSS_GAMA_OLA.csv', index=False)


# fig, axs = plt.subplots(1, 2, figsize=(12, 6), tight_layout=True)

# for key in data_dict.keys():
#     if len(data_dict[key]) > 6:
#         axs[0].scatter(data_dict[key][9]/data_dict[key][5], data_dict[key][8]/data_dict[key][4], color=cd_WHAN[data_dict[key][1]][0], alpha = 0.4, marker='.')
#         axs[1].scatter(data_dict[key][9]/data_dict[key][5], data_dict[key][8]/data_dict[key][4],color=color_dict[data_dict[key][0]][0], alpha = 0.4, marker='.')

# axs[0].set_ylabel(r'$H\alpha flux. SDSS / H\alpha flux. GAMA$')
# axs[1].set_xlabel(r'$H\alpha cont. SDSS / H\alpha cont. GAMA$')
# axs[0].set_xlabel(r'$H\alpha cont. SDSS / H\alpha cont. GAMA$')
# axs[0].plot([0, 1], [0, 1])
# axs[1].plot([0, 1], [0, 1])
# axs[0].set_ylim(-1.3, 1.3)
# axs[1].set_ylim(-1.3, 1.3)

# plt.savefig('./FIGURES/flux_cont.pdf')

# fig, axs = plt.subplots(1, 2, figsize=(12, 6), tight_layout=True)

# for key in data_dict.keys():
#     if len(data_dict[key]) > 6:
#         axs[0].scatter(data_dict[key][5], data_dict[key][9], color=cd_WHAN[data_dict[key][1]][0], alpha = 0.4, marker='.')
#         axs[1].scatter(data_dict[key][5], data_dict[key][9], color=color_dict[data_dict[key][0]][0], alpha = 0.4, marker='.')

# axs[0].set_ylabel(r'$H\alpha cont. SDSS$')
# axs[1].set_xlabel(r'$H\alpha cont. GAMA$')
# axs[0].set_xlabel(r'$H\alpha cont. GAMA$')

# plt.savefig('./FIGURES/HA_cont.pdf')

# fig, axs = plt.subplots(1, 2, figsize=(12, 6), tight_layout=True)

# for i, item in enumerate(HgF_cont_g):
#     axs[0].errorbar(item, hgamma_cont[i], yerr=hgamma_cont_err[i], color=cd_WHAN[SC_WHAN[i]][0], alpha = 0.4, marker='.')
#     axs[1].errorbar(item, hgamma_cont[i], yerr=hgamma_cont_err[i], color=color_dict[BPT[i]][0], alpha = 0.4, marker='.')

#     axs[0].set_ylabel(r'$H\gamma cont. SDSS$')
#     axs[1].set_xlabel(r'$H\gamma cont. GAMA$')
#     axs[0].set_xlabel(r'$H\gamma cont. GAMA$')
    
# gama_x = np.array(HgF_cont_g)
# sdss_y = np.array(hgamma_cont)

# popt, pcov = curve_fit(func, gama_x, sdss_y)

# axs[0].plot(gama_x, popt[0]*gama_x + popt[1])
# axs[1].plot(gama_x, popt[0]*gama_x + popt[1])

# print(popt[0], popt[1])
# plt.savefig('./FIGURES/HG_cont.pdf')

# fig, axs = plt.subplots(1, 2, figsize=(12, 6), tight_layout=True)

# for i, item in enumerate(HA_cont_g):
#     axs[0].errorbar(item, ha_cont[i], yerr=ha_cont_err[i], color=cd_WHAN[SC_WHAN[i]][0], alpha = 0.4, marker='.')
#     axs[1].errorbar(item, ha_cont[i], yerr=ha_cont_err[i], color=color_dict[BPT[i]][0], alpha = 0.4, marker='.')

#     axs[0].set_ylabel(r'$H\alpha cont. SDSS$')
#     axs[1].set_xlabel(r'$H\alpha cont. GAMA$')
#     axs[0].set_xlabel(r'$H\alpha cont. GAMA$')

# gama_x = np.array(HA_cont_g)
# sdss_y = np.array(ha_cont)

# popt, pcov = curve_fit(func, gama_x, sdss_y)

# axs[0].plot(gama_x, popt[0]*gama_x + popt[1])
# axs[1].plot(gama_x, popt[0]*gama_x + popt[1])

# print(popt[0], popt[1])
    
# plt.savefig('./FIGURES/HA_cont.pdf')