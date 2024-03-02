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

        self.list_names_BPT = ['AGN', 'UNC', 'SF', 'NOEL']
        self.list_names_WHAN = ['sAGN', 'wAGN', 'SF', 'ELR', 'NER', 'LLR']

        self.BMS_dict= {
            0 : ['.', 25], #bms
            1 : ['*', 25], #ms
            -1 : ['x', 25] #??????
        }
    
    def reading(self):
        with open(self.ola_file, 'r') as input:
            lines = input.readlines()
        lines_stripped = [line.strip() for line in lines]
        for line in lines_stripped:
            GAMAID = int(line.split()[0])
            try:
                X = float(line.split()[125])
                Y = float(line.split()[1011])
                Y_err = float(line.split()[1012])
                if Y_err*2 < Y and Y + Y_err > 0:
                    Y_up = np.log10(Y + Y_err)
                    Y_down = np.log10(Y - Y_err)
                    Y = np.log10(Y)
                    self.base_dict.update({GAMAID: [X, Y, Y_up, Y_down]})
            except:
                pass
   
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
        gs_top = plt.GridSpec(1, 2, wspace=0)
        self.fig1 = plt.figure(figsize=(12, 6), tight_layout=True)

        self.ax4 = self.fig1.add_subplot(gs_top[:,0])
        self.ax5 = self.fig1.add_subplot(gs_top[:,1], sharey=self.ax4)

        self.topaxes = [self.ax5, self.ax4]

        self.ax4.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, right=True, direction='in')
        self.ax5.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, left=True, labelleft=False, right=True, labelright=False, direction='in')

        self.ax4.set_ylabel(r'$log(\sigma), \: Mpc^{-2}$')
        for ax in self.topaxes:  
            ax.set_xlabel(r'$log(age/yr)$') 
            ax.set_xlim([8.7, 10.0])
            ax.set_ylim([-2, 3])

        bids = [[8.8, 9.0], [9.0, 9.2], [9.2, 9.4], [9.4, 9.6], [9.6, 9.8], [9.8, 10.0]]
        Main.plotter_mdms_BPT(self, bids)
        Main.plotter_mdms_WHAN(self, bids)
        self.fig1.savefig('./FIGURES_IN_PAPER/SurfaceDensity_age.pdf')
        #plt.show()

    def plotter_mdms_BPT(self, bids):

        tot_age = []
        tot = []
        tot_up = []
        tot_down = []
        BPT_keys = []

        for key in self.base_dict.keys():
            if len(self.base_dict[key]) > 5:
                X = self.base_dict[key][0]
                Y = self.base_dict[key][1]
                Y_up = self.base_dict[key][2]
                Y_down = self.base_dict[key][3]
                AGN = self.base_dict[key][4]
                #Y = hp.ms_func(item[y], item[x], item['Z'])
                self.ax4.scatter(X, Y, alpha= 0.3, color = self.color_dict_BPT[AGN][0], marker='.', s = 30)
                tot_age.append(X)
                tot.append(Y)
                tot_up.append(Y_up)
                tot_down.append(Y_down)
                BPT_keys.append(AGN)
                #marker = self.color_dict[AGN][2]
        
        print(len(tot))
        
        class_list = class_list_creator_w_err(tot_age, tot, tot_up, tot_down, BPT_keys, 'BPT')
        
        errs = []
        means = []
        for item in class_list:
            X_plot = []
            Y_plot = []
            err_plot = []
            err_up = []
            err_down = []
            # up_lim_end = []
            
            X, Y, err, length = monte_carlo(item[0], item[1], item[2], item[3], bids)
            # up_lim = up_lim_analysis(item[0], item[5], bids)
            means.append(Y)
            errs.append(err)
            for j in range(len(X)):
                if Y[j] != -99:
                    # axis.errorbar(X[j], Y[j], alpha = 1, xerr=0, yerr= err[j], color=item[4][0], fmt=item[4][1], ms = 12)
                    # axis.scatter(X[j], Y[j], alpha = 1, color=item[4][0], marker=item[4][1], s = 100)
                    # self.ax4.text(X[j], Y[j], length[j], c = 'red')
                    X_plot.append(X[j])
                    Y_plot.append(Y[j])
                    err_plot.append(err[j])
                    # up_lim_end.append(up_lim[j])
                    
            X_plot = np.asarray(X_plot)
            Y_plot = np.asarray(Y_plot)
            for err in err_plot:
                err_up.append(err[1][0])
                err_down.append(err[0][0])
            err_up = np.asarray(err_up)
            err_down = np.asarray(err_down)
            
            self.ax4.fill_between(X_plot, Y_plot + err_up, Y_plot - err_down, color = item[4][0], alpha = 0.17)
            self.ax4.scatter(X_plot, Y_plot, alpha = 1, color=item[4][0], marker=item[4][1], s = 100)
            # for i, elem in enumerate(up_lim_end):
            #     if elem:
            #         self.ax4.arrow(X_plot[i], Y_plot[i], 0, -0.3, width = 0.007, alpha = 1, color=item[4][0])
            self.ax4.plot(X_plot, Y_plot + err_up, alpha = 1, color=item[4][0])
            self.ax4.plot(X_plot, Y_plot - err_down, alpha = 1, color=item[4][0])
            self.ax4.plot(X_plot, Y_plot, alpha = 1, color=item[4][0], linestyle = '--')

        for j, item in enumerate(class_list):
            self.ax4.scatter(-99, -99, alpha = 1, color=item[4][0], marker=item[4][1], s = 150, label=list_names_BPT_1[j])
        self.ax4.legend(loc=3, fontsize='13')

    
    def plotter_mdms_WHAN(self, bids):

        
        tot_age = []
        tot = []
        tot_up = []
        tot_down = []
        WHAN_keys = []

        for key in self.base_dict.keys():
            if len(self.base_dict[key]) > 5:
                X = self.base_dict[key][0]
                Y = self.base_dict[key][1]
                Y_up = self.base_dict[key][2]
                Y_down = self.base_dict[key][3]
                AGN = self.base_dict[key][5]
                self.ax5.scatter(X, Y, alpha= 0.3, color = self.color_dict_WHAN[AGN][0], marker='.', s = 30)
                tot_age.append(X)
                tot.append(Y)
                tot_up.append(Y_up)
                tot_down.append(Y_down)
                #marker = self.color_dict[AGN][2]
                WHAN_keys.append(AGN)
        
        #class_list = [[yes_temp_age, yes_temp, 'midnightblue'], [UNC_temp_age, UNC_temp, 'red'], [no_temp_age, no_temp, 'mediumvioletred'], [noel_temp_age, noel_temp, 'orchid'], [llr_age, llr, 'maroon']]
        
        class_list = class_list_creator_w_err(tot_age, tot, tot_up, tot_down, WHAN_keys, 'WHAN')
        self.means = []
        self.errs = []
        self.ages = []

        errs = []
        means = []
        for item in class_list:
            X_plot = []
            Y_plot = []
            err_plot = []
            err_up = []
            err_down = []
            # up_lim_end = []
            
            X, Y, err, length = monte_carlo(item[0], item[1], item[2], item[3], bids)
            
            # up_lim = up_lim_analysis(item[0], item[5], bids)
            means.append(Y)
            errs.append(err)
            for j in range(len(X)):
                if Y[j] != -99:
                    # axis.errorbar(X[j], Y[j], alpha = 1, xerr=0, yerr= err[j], color=item[4][0], fmt=item[4][1], ms = 12)
                    # axis.scatter(X[j], Y[j], alpha = 1, color=item[4][0], marker=item[4][1], s = 100)
                    # self.ax5.text(X[j], Y[j], length[j], c = 'red')
                    X_plot.append(X[j])
                    Y_plot.append(Y[j])
                    err_plot.append(err[j])
                    # up_lim_end.append(up_lim[j])
                    
            X_plot = np.asarray(X_plot)
            Y_plot = np.asarray(Y_plot)
            for err in err_plot:
                err_up.append(err[1][0])
                err_down.append(err[0][0])
            err_up = np.asarray(err_up)
            err_down = np.asarray(err_down)
            
            self.ax5.fill_between(X_plot, Y_plot + err_up, Y_plot - err_down, color = item[4][0], alpha = 0.17)
            self.ax5.scatter(X_plot, Y_plot, alpha = 1, color=item[4][0], marker=item[4][1], s = 100)
            # for i, elem in enumerate(up_lim_end):
            #     if elem:
            #         self.ax5.arrow(X_plot[i], Y_plot[i], 0, -0.3, width = 0.007, alpha = 1, color=item[4][0])
            self.ax5.plot(X_plot, Y_plot + err_up, alpha = 1, color=item[4][0])
            self.ax5.plot(X_plot, Y_plot - err_down, alpha = 1, color=item[4][0])
            self.ax5.plot(X_plot, Y_plot, alpha = 1, color=item[4][0], linestyle = '--')

        for j, item in enumerate(class_list):
            self.ax5.scatter(-99, -99, alpha = 1, color=item[4][0], marker=item[4][1], s = 150, label=list_names_WHAN[j])
        self.ax5.legend(loc=3, fontsize='13')

if __name__ == '__main__':
    obj = Main('E:\LICENSE\ProgsData\main\GAMAforOleg_1.txt', 'GAMA_ETG_OLA.csv')
    obj.reading()
    obj.matching()
    obj.plotting_mdms_age()