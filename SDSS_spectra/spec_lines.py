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

gama_path = r"E:\backup\ALMA9\GAMA_ETG_OLA.csv"

dataframe = pd.read_csv(filepath_or_buffer=gama_path, usecols=['SPEC_ID', 'SURVEY', 'BPT', 'WHAN', 'HgF', 'HgF_cont', 'HgF_er', 'HA', 'HA_cont', 'HA_er'])
headers = []
with open(r'E:\backup\ALMA9\SDSS_spectra\cols.txt', 'r') as f:
    lines = f.readlines()
    lines_stripped = [line.strip() for line in lines]
    for line in lines_stripped:
        headers.append(line.split()[0])

data_dict = {}
for i, item in enumerate(dataframe['SPEC_ID']):
    if dataframe['SURVEY'][i] == 'SDSS':
        data_dict.update({item : [dataframe['BPT'][i], dataframe['WHAN'][i], dataframe['HgF'][i], dataframe['HgF_er'][i], dataframe['HgF_cont'][i], dataframe['HA'][i], dataframe['HA_er'][i], dataframe['HA_cont'][i]]})


for i, item in enumerate(dataframe_sdss['SPECOBJID']):
    spec_id = dataframe_sdss['SPECOBJID'][i]
    hgamma = dataframe_sdss['H_GAMMA_FLUX'][i]
    hgamma_err = dataframe_sdss['H_GAMMA_FLUX_ERR'][i]
    hgamma_cont = dataframe_sdss['H_GAMMA_CONT'][i]
    hgamma_cont_err = dataframe_sdss['H_GAMMA_CONT_ERR'][i]
    ha = dataframe_sdss['H_ALPHA_FLUX'][i]
    ha_err = dataframe_sdss['H_ALPHA_FLUX_ERR'][i]
    ha_cont = dataframe_sdss['H_ALPHA_CONT'][i]
    ha_cont_err = dataframe_sdss['H_ALPHA_CONT_ERR'][i]
    if spec_id in data_dict.keys():
        data_dict[spec_id].extend([hgamma, hgamma_cont, hgamma_err, hgamma_cont_err, ha, ha_cont, ha_err, ha_cont_err])

spec_id = []
BPT = []
SC_WHAN = []
HgF_g = []
HgF_er_g = []
HgF_cont_g = []
HA_g = []
HA_er_g = []
HA_cont_g = []

hgamma = []
hgamma_cont = []
hgamma_err = []
hgamma_cont_err = []
ha = []
ha_cont = []
ha_err = []
ha_cont_err = []

for key in data_dict.keys():
    if len(data_dict[key]) > 8:
        spec_id.append(key)
        BPT.append(data_dict[key][0])
        SC_WHAN.append(data_dict[key][1])
        HgF_g.append(data_dict[key][2])
        HgF_er_g.append(data_dict[key][3])
        HgF_cont_g.append(data_dict[key][4])
        HA_g.append(data_dict[key][5])
        HA_er_g.append(data_dict[key][6])
        HA_cont_g.append(data_dict[key][7])
        hgamma.append(data_dict[key][8])
        hgamma_cont.append(data_dict[key][9])
        hgamma_err.append(data_dict[key][10])
        hgamma_cont_err.append(data_dict[key][11])
        ha.append(data_dict[key][12])
        ha_cont.append(data_dict[key][13])
        ha_err.append(data_dict[key][14])
        ha_cont_err.append(data_dict[key][15])

out_dict = {
    'SPEC_ID' : spec_id,
    'BPT' : BPT,
    'WHAN' : SC_WHAN,
    'HgF_g' : HgF_g,
    'HgF_er_g' : HgF_er_g,
    'HgF_cont_g' : HgF_cont_g,
    'Hg_SDSS' : hgamma,
    'Hg_er_SDSS' : hgamma_err,
    'Hg_cont_SDSS' : hgamma_cont,
    'Hg_cont_err_SDSS' : hgamma_cont_err,

    'HA_g' : HA_g,
    'HA_er_g' : HA_er_g,
    'HA_cont_g' : HA_cont_g,
    'HA_SDSS' : ha,
    'HA_er_SDSS' : ha_err,
    'HA_cont_SDSS' : ha_cont,
    'HA_cont_err_SDSS' : ha_cont_err,
}

df = pd.DataFrame(out_dict)
df.to_csv('SDSS_GAMA.csv', index=False)


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

fig, axs = plt.subplots(1, 2, figsize=(12, 6), tight_layout=True)

for i, item in enumerate(HgF_cont_g):
    axs[0].errorbar(item, hgamma_cont[i], yerr=hgamma_cont_err[i], color=cd_WHAN[SC_WHAN[i]][0], alpha = 0.4, marker='.')
    axs[1].errorbar(item, hgamma_cont[i], yerr=hgamma_cont_err[i], color=color_dict[BPT[i]][0], alpha = 0.4, marker='.')

    axs[0].set_ylabel(r'$H\gamma cont. SDSS$')
    axs[1].set_xlabel(r'$H\gamma cont. GAMA$')
    axs[0].set_xlabel(r'$H\gamma cont. GAMA$')
    
gama_x = np.array(HgF_cont_g)
sdss_y = np.array(hgamma_cont)

popt, pcov = curve_fit(func, gama_x, sdss_y)

axs[0].plot(gama_x, popt[0]*gama_x + popt[1])
axs[1].plot(gama_x, popt[0]*gama_x + popt[1])

print(popt[0], popt[1])
plt.savefig('./FIGURES/HG_cont.pdf')

fig, axs = plt.subplots(1, 2, figsize=(12, 6), tight_layout=True)

for i, item in enumerate(HA_cont_g):
    axs[0].errorbar(item, ha_cont[i], yerr=ha_cont_err[i], color=cd_WHAN[SC_WHAN[i]][0], alpha = 0.4, marker='.')
    axs[1].errorbar(item, ha_cont[i], yerr=ha_cont_err[i], color=color_dict[BPT[i]][0], alpha = 0.4, marker='.')

    axs[0].set_ylabel(r'$H\alpha cont. SDSS$')
    axs[1].set_xlabel(r'$H\alpha cont. GAMA$')
    axs[0].set_xlabel(r'$H\alpha cont. GAMA$')

gama_x = np.array(HA_cont_g)
sdss_y = np.array(ha_cont)

popt, pcov = curve_fit(func, gama_x, sdss_y)

axs[0].plot(gama_x, popt[0]*gama_x + popt[1])
axs[1].plot(gama_x, popt[0]*gama_x + popt[1])

print(popt[0], popt[1])
    
plt.savefig('./FIGURES/HA_cont.pdf')