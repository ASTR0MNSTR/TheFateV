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
    def __init__(self, ola_file, out):
        self.ola_file = ola_file
        self.out = out
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
    
    def reading(self):
        k = 0
        t = 0
        with open(self.ola_file, 'r') as input:
            lines = input.readlines()
        lines_stripped = [line.strip() for line in lines]
        for line in lines_stripped:
            GAMAID = int(line.split()[0])
            try:
                m_s = float(line.split()[35])
                age = float(line.split()[125])
                Y = float(line.split()[1011])
                Y_err = float(line.split()[1012])
                Y_flag = float(line.split()[1013])
                self.base_dict.update({GAMAID: [m_s, age, Y, Y_err, Y_flag]})
                t += 1
            except:
                pass
        
        print(f'All galaxies with data: {t}')
   
        #print(self.base_dict)

        self.data_dict = []

        with open(self.out, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                self.data_dict.append({'GAMAID': int(row['CATAID_1']), 'AGN': row['BPT'], 'SC_WHAN' : row['WHAN']})
                
        #print(self.data_dict)

    def matching(self):
        for item in self.data_dict:
            try:
                self.base_dict[item['GAMAID']].extend([item['AGN'], item['SC_WHAN']])
            except KeyError:
                pass
            
        # print(self.base_dict)
    
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
        ax1, ax2 = Main.plotter(self, bids_mass, [ax1, ax2], 0, 0)
        bids_age = [[8.8, 9.0], [9.0, 9.2], [9.2, 9.4], [9.4, 9.6], [9.6, 9.8], [9.8, 10.0]]
        ax3, ax4 = Main.plotter(self, bids_age, [ax3, ax4], 1, 1)
        self.fig1.savefig('./FIGURES_IN_PAPER/SurfaceDensity.pdf', dpi=300, transparent = True, bbox_inches = 'tight', pad_inches = 0.0001)
        #plt.show()

    def plotter(self, bids, axes, xkey, legend_key):

        tot_age = []
        tot = []
        tot_up = []
        tot_down = []
        BPT_keys = []
        WHAN_keys = []
        ks = []

        for key in self.base_dict.keys():
            if len(self.base_dict[key]) == 7:
                X = self.base_dict[key][xkey]
                Y = self.base_dict[key][2]
                Y_er = self.base_dict[key][3]
                Y_flag = self.base_dict[key][4]
                BPT = self.base_dict[key][5]
                WHAN = self.base_dict[key][6]
                if Y > 0 and Y_er > 0 and Y - 2*Y_er >= 0 and Y_flag == 0:
                    Y_up = np.log10(Y + Y_er)
                    Y_down = np.log10(Y - Y_er)
                    Y = np.log10(Y)
                    axes[0].scatter(X, Y, alpha= 0.5, color = self.color_dict_BPT[BPT][0], marker='.', s = 30)
                    axes[1].scatter(X, Y, alpha= 0.5, color = self.color_dict_WHAN[WHAN][0], marker='.', s = 30)
                    k = 1
                elif (Y > 0 and Y_er == -999.9) or (Y_flag in [1, 2]): #case 2
                    Y_up = np.log10(Y)
                    Y_down = np.log10(Y)
                    Y = np.log10(Y)
                    axes[0].arrow(X, Y, 0, -0.1, head_width=0.01, head_length=0.03, color=self.color_dict_BPT[BPT][0], alpha=0.5)
                    axes[1].arrow(X, Y, 0, -0.1, head_width=0.01, head_length=0.03, color=self.color_dict_WHAN[WHAN][0], alpha=0.5)
                    k = 0
                    
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
                axes[0].scatter(-99, -99, alpha = 1, color=item[4][0], marker=item[4][1], s = 150, label=list_names_BPT_1[j])
            axes[0].legend(loc=3, fontsize=15)

            for j, item in enumerate(class_list_WHAN):
                axes[1].scatter(-99, -99, alpha = 1, color=item[4][0], marker=item[4][1], s = 150, label=list_names_WHAN[j])
            axes[1].legend(loc=3, fontsize=15)
        
        return axes

if __name__ == '__main__':
    obj = Main('E:\LICENSE\ProgsData\main\GAMAforOleg_1.txt', 'GAMA_ETG_OLA.csv')
    obj.reading()
    obj.matching()
    obj.plotting_mdms_age()