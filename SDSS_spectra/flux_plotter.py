import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from _global_ import *
from scipy.optimize import curve_fit
import matplotlib as mpl 

def func(x, a, b):
    return a*x + b

#creating colormaps

class Main:
    def __init__(self, path_in):
        self.path_in = path_in
        self.df = None
    
    def reading_plotting(self):
        self.df = pd.read_csv(self.path_in)

        Main.plotter(self, 'HgF_cont', 'Hg_cont_SDSS', r'$H\gamma \; cont. \; GAMA$', [0, 150], r'$H\gamma \; cont. \; SDSS$', [0, 40], './FIGURES/GAMA_FLUXES/HG_cont_FULL.pdf')
        print('Hgamma cont. plotted')
        Main.plotter(self, 'HA_cont', 'HA_cont_SDSS', r'$H\alpha \; cont. \; GAMA$', [0, 200], r'$H\alpha \; cont. \; SDSS$', [0, 40], './FIGURES/GAMA_FLUXES/HA_cont_FULL.pdf')
        print('Halpha cont. plotted')
        Main.plotter(self, 'HgF', 'Hg_SDSS', r'$H\gamma \; flux \; GAMA$', [-50, 1000], r'$H\gamma \; flux \; SDSS$', [-50, 200], './FIGURES/GAMA_FLUXES/HG_FULL.pdf')
        print('Hgamma flux plotted')
        Main.plotter(self, 'HA', 'HA_SDSS', r'$H\alpha \; flux \; GAMA$', [-200, 5000], r'$H\alpha \; flux \; SDSS$', [-200, 1000], './FIGURES/GAMA_FLUXES/HA_FULL.pdf')
        print('Halpha flux plotted')
        Main.sn(self, 'HA', 'HA_er', 'HA_SDSS', 'HA_er_SDSS', [-50, 200], r'$SN H\alpha GAMA$', r'$SN H\alpha SDSS$', [-50, 200], './FIGURES/SN_HA_FULL.pdf')
        Main.sn(self, 'HgF', 'HgF_er', 'Hg_SDSS', 'Hg_er_SDSS', [-20, 100], r'$SN H\gamma GAMA$', r'$SN H\gamma SDSS$', [-20, 100], './FIGURES/SN_HG_FULL.pdf')
    
    def plotter(self, X, Y, Xname, Xlim, Yname, Ylim, file_name):
        norm = mpl.colors.Normalize(vmin=0,vmax=8)
        # choose a colormap
        c_m = mpl.cm.jet
        # create a ScalarMappable and initialize a data structure
        s_m = mpl.cm.ScalarMappable(cmap=c_m, norm=norm)
        s_m.set_array([])

        fig, axs = plt.subplots(1, 2, figsize=(12, 6), tight_layout=True)

        for i, item in enumerate(self.df[X]):
    #axs[0].errorbar(item, df['Hg_cont_SDSS'][i], yerr=df['Hg_cont_err_SDSS'][i], color=cd_WHAN[df['WHAN'][i]][0], alpha = 1, marker='.')
            axs[0].scatter(item, self.df[Y][i], color=cd_WHAN[self.df['WHAN'][i]][0], alpha = 1, marker='.')
            if 0 < self.df['GALRE_i'][i] < 15:
                color = s_m.to_rgba(self.df['GALRE_i'][i])
            elif self.df['GALRE_i'][i] < 0:
                color = 'grey'
            elif self.df['GALRE_i'][i] >= 15:
                color = 'black'

    #axs[1].errorbar(item, df['Hg_cont_SDSS'][i], yerr=df['Hg_cont_err_SDSS'][i], color=color, alpha = 1, marker='.')
            axs[1].scatter(item, self.df[Y][i], color=color, alpha = 1, marker='.')

        axs[0].set_ylabel(Yname)
        cbar = fig.colorbar(s_m, ax=axs[1], location='right')
        cbar.ax.set_ylabel(r'$R_{eff} \; i, \; arcsec$')

        axes = [axs[0], axs[1]]
        for ax in axes:
            ax.plot([-100, 1000], [-100, 1000], c='red')
            ax.set_xlabel(Xname)
            ax.set_xlim(Xlim[0], Xlim[1])
            ax.set_ylim(Ylim[0], Ylim[1])
        
        gama_x = np.array(self.df[X])
        sdss_y = np.array(self.df[Y])
        popt, pcov = curve_fit(func, gama_x, sdss_y)

        axs[0].plot(gama_x, popt[0]*gama_x + popt[1])
        axs[1].plot(gama_x, popt[0]*gama_x + popt[1])

        axs[0].axhline(y = 0, linestyle='--')
        axs[0].axvline(x = 0, linestyle='--')
        axs[1].axhline(y = 0, linestyle='--')
        axs[1].axvline(x = 0, linestyle='--')

        print(popt[0], popt[1])
    

        fig.savefig(file_name)

    def sn(self, X_up, X_down, Y_up, Y_down, Xlim, Xname, Yname, Ylim, file_name):
        norm = mpl.colors.Normalize(vmin=0,vmax=8)
        # choose a colormap
        c_m = mpl.cm.jet
        # create a ScalarMappable and initialize a data structure
        s_m = mpl.cm.ScalarMappable(cmap=c_m, norm=norm)
        s_m.set_array([])


        fig, axs = plt.subplots(1, 2, figsize=(12, 6), tight_layout=True)

        for i, item in enumerate(self.df[X_up]):
            axs[0].scatter(item/self.df[X_down][i], self.df[Y_up][i]/self.df[Y_down][i], color=cd_WHAN[self.df['WHAN'][i]][0], alpha = 1, marker='.')
            if 0 < self.df['GALRE_i'][i] < 15:
                color = s_m.to_rgba(self.df['GALRE_i'][i])
            elif self.df['GALRE_i'][i] < 0:
                color = 'grey'
            elif self.df['GALRE_i'][i] >= 15:
                color = 'black'        
            axs[1].scatter(item/self.df[X_down][i], self.df[Y_up][i]/self.df[Y_down][i], color=color, alpha = 1, marker='.')

        axs[0].set_ylabel(Yname)
        cbar = fig.colorbar(s_m, ax=axs[1], location='right')
        cbar.ax.set_ylabel(r'$R_{eff} \; i, \; arcsec$')

        axes = [axs[0], axs[1]]
        for ax in axes:
            ax.plot([-100, 1000], [-100, 1000], c='red')
            ax.set_xlabel(Xname)
            ax.set_xlim(Xlim[0], Xlim[1])
            ax.set_ylim(Ylim[0], Ylim[1])
            
        gama_x = np.array(self.df[X_up]/self.df[X_down])
        sdss_y = np.array(self.df[Y_up]/self.df[Y_down])
        popt, pcov = curve_fit(func, gama_x, sdss_y)

        axs[0].plot(gama_x, popt[0]*gama_x + popt[1])
        axs[1].plot(gama_x, popt[0]*gama_x + popt[1])

        axs[0].axhline(y = 0, linestyle='--')
        axs[0].axvline(x = 0, linestyle='--')
        axs[1].axhline(y = 0, linestyle='--')
        axs[1].axvline(x = 0, linestyle='--')

        print(popt[0], popt[1])
            
        plt.savefig(file_name)
    
if __name__ == '__main__':
    obj = Main('E:/databases/SDSS_GAMA_5000.csv')
    obj.reading_plotting()


# fig, axs = plt.subplots(1, 2, figsize=(12, 6), tight_layout=True)

# for i, item in enumerate(df['HgF']):
#     axs[0].scatter(item/df['HgF_er'][i], df['Hg_SDSS'][i]/df['Hg_er_SDSS'][i], color=cd_WHAN[df['WHAN'][i]][0], alpha = 1, marker='.')
#     if 0 < df['GALRE_i'][i] < 15:
#         color = s_m.to_rgba(df['GALRE_i'][i])
#     elif df['GALRE_i'][i] < 0:
#         color = 'grey'
#     elif df['GALRE_i'][i] >= 15:
#         color = 'black'
#     axs[1].scatter(item/df['HgF_er'][i], df['Hg_SDSS'][i]/df['Hg_er_SDSS'][i], color=color, alpha = 1, marker='.')

# axs[0].set_ylabel(r'$SN H\gamma SDSS$')
# axs[1].set_xlabel(r'$SN H\gamma GAMA$')
# axs[0].set_xlabel(r'$SN H\gamma GAMA$')
# cbar = fig.colorbar(s_m, ax=axs[1], location='right')
# cbar.ax.set_ylabel(r'$R_{eff} \; i, \; arcsec$')

# gama_x = np.array(df['HgF']/df['HgF_er'])
# sdss_y = np.array(df['Hg_SDSS']/df['Hg_er_SDSS'])
# popt, pcov = curve_fit(func, gama_x, sdss_y)

# axs[0].plot(gama_x, popt[0]*gama_x + popt[1])
# axs[1].plot(gama_x, popt[0]*gama_x + popt[1])

# axs[0].axhline(y = 0, linestyle='--')
# axs[0].axvline(x = 0, linestyle='--')
# axs[1].axhline(y = 0, linestyle='--')
# axs[1].axvline(x = 0, linestyle='--')

# print(popt[0], popt[1])
    
# plt.savefig('./FIGURES/SN_HG.pdf')
