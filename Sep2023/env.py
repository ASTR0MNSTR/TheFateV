import csv
import matplotlib.pyplot as plt 
import numpy as np
import astropy
from astropy.cosmology import FlatLambdaCDM
import math
# import statistics as st
import pandas as pd
# from scipy.optimize import curve_fit
from __legpars__ import *
from __stats__ import *
from __plt__ import *
from __algo__ import *
from __reader__ import *


class hp:
    def log_er(item):
        return [[abs(math.log(1-item, 10))], [abs(math.log(1+item, 10))]]
    
    def cosmic_time(z):
        cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
        return cosmo.age(z).value
    
    def ms_func(SFR, MS, z):
        SFR_teor = (0.84 - 0.026*hp.cosmic_time(z))*MS - (6.51 - 0.11*hp.cosmic_time(z))
        return SFR - SFR_teor
    
class Main(hp):
    def __init__(self, path_to_data):
        self.path_to_data = path_to_data
        self.base_dict = {}
        self.data_dict = []
        self.color_dict_BPT = color_dict_BPT
        self.color_dict_WHAN = cd_WHAN
        self.color_dict_leg = color_dict_leg

        self.list_names_BPT = ['AGN', 'UNC', 'SFG', 'NOEL']
        self.list_names_WHAN = ['sAGN', 'wAGN', 'SGF', 'ELR', 'NER', 'LLR']

        self.BMS_dict= {
            0 : ['.', 25], #bms
            1 : ['*', 25], #ms
            -1 : ['x', 25] #??????
        }
    
    def plotting_mdms_age(self):
        gs_top = plt.GridSpec(2, 2, wspace=0, hspace=0.15)
        self.fig1 = plt.figure(figsize=(12, 12))
        adjusting_plotting_pars()
        # adjusting_figure_size(12, 12, 0.8, 0.2, 0.6, 0.2)
        ax1 = self.fig1.add_subplot(gs_top[0,0])
        ax2 = self.fig1.add_subplot(gs_top[0,1], sharey=ax1)
        ax3 = self.fig1.add_subplot(gs_top[1,0])
        ax4 = self.fig1.add_subplot(gs_top[1,1], sharey=ax3)

        self.topaxes = [ax1, ax2]
        self.bottomaxes = [ax3, ax4]

        ax1.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, right=True, direction='in', labelsize=17)
        ax3.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, right=True, direction='in', labelsize=17)
        ax2.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, left=True, labelleft=False, right=True, labelright=False, direction='in', labelsize=17)
        ax4.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, left=True, labelleft=False, right=True, labelright=False, direction='in', labelsize=17)

        ax1.set_ylabel(r'$\log \mathrm{(\sigma), \: Mpc^{-2}}$')
        ax1.set_xlabel(r'$\log \mathrm{(M_{stellar} \: / \: M_\odot)}$')
        ax2.set_xlabel(r'$\log \mathrm{(M_{stellar} \: / \: M_\odot)}$')
        ax3.set_ylabel(r'$\log \mathrm{(\sigma), \: Mpc^{-2}}$')
        ax3.set_xlabel(r'$\log \mathrm{(age \: / \: yr)}$')
        ax4.set_xlabel(r'$\log \mathrm{(age \: / \: yr)}$')
        for ax in self.topaxes:  
            ax.set_xlim([9.9, 11.7])
            ax.set_ylim([-2, 3])
            ax.axhline(1, color = 'k', linestyle='-')
            ax.set_xticks(np.arange(10.0, 11.7, 0.3))
        
        for ax in self.bottomaxes:  
            ax.set_xlim([8.7, 10.1])
            ax.set_ylim([-2, 3])
            ax.axhline(1, color = 'k', linestyle='-')
            ax.set_xticks(np.arange(8.8, 10.1, 0.2))
        
        k = 0.8
        generating_annotation(ax1, 9.9 + k*(11.6 - 9.9), -2 + k*(3 - (-2)), 'BPT')    
        generating_annotation(ax2, 9.9 + k*(11.6 - 9.9), -2 + k*(3 - (-2)), 'WHAN')

        bids_mass = [[10.0, 10.25], [10.25, 10.5], [10.5, 10.75], [10.75, 11], [11, 11.25], [11.25, 11.5]]
        ax1, ax2 = Main.plotter(self, bids_mass, [ax1, ax2], 'mass_stellar_percentile50', 0)
        bids_age = [[8.8, 9.0], [9.0, 9.2], [9.2, 9.4], [9.4, 9.6], [9.6, 9.8], [9.8, 10.0]]
        ax3, ax4 = Main.plotter(self, bids_age, [ax3, ax4], 'ager_percentile50', 1)
        self.fig1.savefig('./FIGURES_IN_PAPER_DR4/SurfaceDensity.pdf', dpi=300, transparent = True, bbox_inches = 'tight', pad_inches = 0.0001)
        #plt.show()

    def plotter(self, bids, axes, xkey, legend_key):

        tot_age = []
        tot = []
        tot_up = []
        tot_down = []
        BPT_keys = []
        WHAN_keys = []
        ks = []

        df = pd.read_csv(self.path_to_data)
        clean_df = df.dropna()
        clean_df.reset_index(inplace=True, drop=True)
        for index, row in clean_df.iterrows():
                X = row[xkey]
                Y = row['SurfaceDensity']
                Y_er = row['SurfaceDensityErr']
                Y_flag = row['SurfaceDensityFlag']
                BPT = row['BPT']
                WHAN = row['WHAN']
                if Y > 0 and Y_er > 0 and Y - 2*Y_er >= 0 and Y_flag == 0:
                    Y_up = np.log10(Y + Y_er)
                    Y_down = np.log10(Y - Y_er)
                    Y = np.log10(Y)
                    axes[0].scatter(X, Y, alpha= 0.5, color = self.color_dict_BPT[BPT][0], marker='.', s = 30)
                    axes[1].scatter(X, Y, alpha= 0.5, color = self.color_dict_WHAN[WHAN][0], marker='.', s = 30)
                    k = 1
                    tot_age.append(X)
                    tot.append(Y)
                    tot_up.append(Y_up)
                    tot_down.append(Y_down)
                    BPT_keys.append(BPT)
                    WHAN_keys.append(WHAN)
                    ks.append(k)

                elif (Y > 0 and Y_er == -999.9) or (Y_flag in [1, 2]): #case 2
                    Y_up = np.log10(Y)
                    Y_down = np.log10(Y)
                    Y = np.log10(Y)
                    axes[0].arrow(X, Y, 0, -0.1, head_width=0.01, head_length=0.03, color=self.color_dict_BPT[BPT][0], alpha=0.5)
                    axes[1].arrow(X, Y, 0, -0.1, head_width=0.01, head_length=0.03, color=self.color_dict_WHAN[WHAN][0], alpha=0.5)
                    k = 0
                    tot_age.append(X)
                    tot.append(Y)
                    tot_up.append(Y_up)
                    tot_down.append(Y_down)
                    BPT_keys.append(BPT)
                    WHAN_keys.append(WHAN)
                    ks.append(k)

                    
                elif Y > 0 and Y_er > 0 and Y - 2*Y_er < 0:
                    Y_up = np.log10(Y + 2*Y_er)
                    Y_down = np.log10(Y + 2*Y_er)
                    Y = np.log10(Y + 2*Y_er)
                    axes[0].arrow(X, Y, 0, -0.1, head_width=0.01, head_length=0.03, color=self.color_dict_BPT[BPT][0], alpha=0.5)
                    axes[1].arrow(X, Y, 0, -0.1, head_width=0.01, head_length=0.03, color=self.color_dict_WHAN[WHAN][0], alpha=0.5)
                    k = -1
                    
                    tot_age.append(X)
                    tot.append(Y)
                    tot_up.append(Y_up)
                    tot_down.append(Y_down)
                    BPT_keys.append(BPT)
                    WHAN_keys.append(WHAN)
                    ks.append(k)


                
                #Y = hp.ms_func(item[y], item[x], item['Z'])

                #marker = self.color_dict[AGN][2]
        
        class_list_BPT = class_list_creator_w_err_out(tot_age, tot, tot_up, tot_down, BPT_keys, 'BPT', ks)
        class_list_WHAN = class_list_creator_w_err_out(tot_age, tot, tot_up, tot_down, WHAN_keys, 'WHAN', ks)
        
        classlist_plotter_uplim(axes[0], class_list_BPT, bids)
        classlist_plotter_uplim(axes[1], class_list_WHAN, bids)

        if legend_key == 1:
            for j, item in enumerate(class_list_BPT):
                axes[0].scatter(-99, -99, alpha = 1, color=item[4][0], marker=item[4][1], s = 150, label=list_names_BPT_1[j], edgecolors='black')
            axes[0].legend(loc=3, fontsize=15)

            for j, item in enumerate(class_list_WHAN):
                axes[1].scatter(-99, -99, alpha = 1, color=item[4][0], marker=item[4][1], s = 150, label=list_names_WHAN[j], edgecolors='black')
            axes[1].legend(loc=3, fontsize=15)
        
        return axes

if __name__ == '__main__':
    obj = Main(r'E:\databases\GAMAs4\DETG_DR4.csv')
    # obj.reading()
    # obj.matching()
    obj.plotting_mdms_age()