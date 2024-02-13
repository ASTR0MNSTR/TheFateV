import csv
import matplotlib.pyplot as plt 
import numpy as np
import math
# import statistics as st
import pandas as pd
# from scipy.optimize import curve_fit
from __legpars__ import *
from __stats__ import *
from __plt__ import *

class hp:
    def log_er(item):
        return [[abs(math.log(1-item, 10))], [abs(math.log(1+item, 10))]]

class Main(hp):
    def __init__(self, ola_file, out):
        self.ola_file = ola_file
        self.out = out
        self.base_dict = {}
        self.data_dict = []
        self.color_dict_BPT = color_dict_BPT
        self.color_dict_WHAN = cd_WHAN
        self.color_dict_leg = color_dict_leg

        self.list_names_BPT = list_names_BPT
        self.list_names_WHAN = list_names_WHAN

        self.BMS_dict= {
            0 : ['o', 12], #bms
            1 : ['*', 12], #ms
            -1 : ['x', 15] #??????
        }
    
    def reading(self):
        with open(self.ola_file, 'r') as input:
            lines = input.readlines()
        lines_stripped = [line.strip() for line in lines]
        for line in lines_stripped:
            GAMAID = int(line.split()[0])
            X = float(line.split()[125])
            Y = float(line.split()[53])
            Y_up = float(line.split()[54])
            Y_down = float(line.split()[52])
            P100 = float(line.split()[724])
            P100_er = float(line.split()[725])
            #SFR_up = float(line.split()[96])
            #SFR_down = float(line.split()[94])
            #SFR_er = [[abs(SFR_down-SFR_0)], [abs(SFR_up-SFR_0)]]
            self.base_dict.update({GAMAID: [X, Y, Y_up, Y_down, P100, P100_er]})

        self.data_dict = []

        with open(self.out, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                self.data_dict.append({'GAMAID': int(row['CATAID_1']), 'AGN': row['BPT'], 'SC_WHAN' : row['WHAN']})

    def matching(self):
        for item in self.data_dict:
            try:
                X, Y, Y_up, Y_down, P100, P100_er = self.base_dict[item['GAMAID']]
                item.update({'X' : X, 'Y' : Y, 'Y_up' : Y_up, 'Y_down' : Y_down, 'P100' : P100, 'P100_er' : P100_er})
            except: 
                print(item['GAMAID'])
    
    def plotting_mdms_age(self):
        gs_top = plt.GridSpec(1, 2, wspace=0)
        self.fig1 = plt.figure(figsize=(12, 6), tight_layout=True)

        self.ax4 = self.fig1.add_subplot(gs_top[:,0])
        self.ax5 = self.fig1.add_subplot(gs_top[:,1], sharey=self.ax4)

        self.topaxes = [self.ax5, self.ax4]

        self.ax4.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, right=True, direction='in')
        self.ax5.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, left=True, labelleft=False, right=True, labelright=False, direction='in')

        self.ax4.set_ylabel(r'$T_{warm \; dust}, K$')
        for ax in self.topaxes:    
            ax.set_xlim([8.7, 10.1])
            ax.set_ylim([30, 61])
            ax.set_xlabel(r'$log(age/yr)$')
            #ax.set_yticks(np.arange(15, 25.9, 2))

        Main.plotter_mdms_BPT(self, 'X', 'Y', 'Y_up', 'Y_down', 'AGN', 'SC_WHAN', True, 'P100', 'P100_er')
        
        self.fig1.savefig('./FIGURES_IN_PAPER/TDW.pdf')
        #plt.show()

    def plotter_mdms_BPT(self, x, y, up, down, BPT_key, WHAN_key, bids_key, p100, p100_er):
        XX = []
        YY = []
        Y_up = []
        Y_down = []
        BPT_keys = []
        WHAN_keys = []
        for item in self.data_dict:
            X = item[x]
            Y = item[y]
            AGN = item[BPT_key]
            WHAN = item[WHAN_key]
            
            P100 = item[p100]
            P100_er = item[p100_er]
            if Y > 0 and P100 > 2*P100_er and P100_er >= 0:
                XX.append(X)
                YY.append(Y)
                Y_up.append(item[up])
                Y_down.append(item[down])
                BPT_keys.append(AGN)
                WHAN_keys.append(WHAN)
                self.ax4.scatter(X, Y, alpha= 0.4, color = self.color_dict_BPT[AGN][0], marker='.', s = self.color_dict_BPT[AGN][1])
                self.ax5.scatter(X, Y, alpha= 0.4, color = self.color_dict_WHAN[WHAN][0], marker='.', s = self.color_dict_WHAN[WHAN][1])
        
        class_list_BPT = class_list_creator_w_err(XX, YY, Y_up, Y_down, BPT_keys, 'BPT')
        class_list_WHAN = class_list_creator_w_err(XX, YY, Y_up, Y_down, WHAN_keys, 'WHAN')

        self.means = []
        self.errs = []
        self.ages = []
        if bids_key == True:
            bids = [[8.8, 9.0], [9.0, 9.2], [9.2, 9.4], [9.4, 9.6], [9.6, 9.8], [9.8, 10.0]]
            ages_const = np.arange(8.8, 10.2, 0.2)
        else:
            bids = [[10.0, 10.25], [10.25, 10.5], [10.5, 10.75], [10.75, 11], [11, 11.25], [11.25, 11.5]]
            ages_const = np.arange(7, 12, 1)

        errs = []
        means = []
        for item in class_list_BPT:
            X_plot = []
            Y_plot = []
            X, Y, err, length = monte_carlo(item[0], item[1], item[2], item[3], bids)
            means.append(Y)
            errs.append(err)
            for j in range(len(X)):
                if Y[j] != -99:
                    self.ax4.errorbar(X[j], Y[j], alpha = 1, xerr=0, yerr= err[j], color=item[4][0], fmt=item[4][1], ms = 12)
                    self.ax4.text(X[j], Y[j], length[j], c = 'red')
                    X_plot.append(X[j])
                    Y_plot.append(Y[j])
            self.ax4.plot(X_plot, Y_plot, alpha = 1, color=item[4][0])

        #for key in self.color_dict_leg.keys():
        #        self.ax4.scatter(-99, -99, alpha= 1, color = self.color_dict_leg[key][0], marker = self.color_dict_leg[key][2], s = self.color_dict_leg[key][1], label=key)
        self.ax4.set_xticks(ages_const)
        for j, item in enumerate(class_list_BPT):
            self.ax4.scatter(-99, -99, alpha = 1, color=item[4][0], marker=item[4][1], s = 150, label=list_names_BPT_1[j])
        self.ax4.legend(loc=3, fontsize='13')
        
        errs = []
        means = []
        for item in class_list_WHAN:
            X_plot = []
            Y_plot = []
            X, Y, err, length = monte_carlo(item[0], item[1], item[2], item[3], bids)
            means.append(Y)
            errs.append(err)
            for j in range(len(X)):
                if Y[j] != -99:
                    self.ax5.errorbar(X[j], Y[j], alpha = 1, xerr=0, yerr= err[j], color=item[4][0], fmt=item[4][1], ms = 12)
                    self.ax5.text(X[j], Y[j], length[j], c = 'red')
                    X_plot.append(X[j])
                    Y_plot.append(Y[j])
            self.ax5.plot(X_plot, Y_plot, alpha = 1, color=item[4][0])

        #for key in self.color_dict_leg.keys():
        #        self.ax4.scatter(-99, -99, alpha= 1, color = self.color_dict_leg[key][0], marker = self.color_dict_leg[key][2], s = self.color_dict_leg[key][1], label=key)
        self.ax5.set_xticks(ages_const)
        for j, item in enumerate(class_list_WHAN):
            self.ax5.scatter(-99, -99, alpha = 1, color=item[4][0], marker=item[4][1], s = 150, label=list_names_WHAN[j])
        self.ax5.legend(loc=3, fontsize='13')

        #self.ax4.set_xlim(8, 10.1)
        #self.ax4.set_ylim(14, 26)
        #self.fig.savefig(name_file + '.pdf')
        #plt.show()

if __name__ == '__main__':
    obj = Main('E:\LICENSE\ProgsData\main\GAMAforOleg_1.txt', 'GAMA_ETG_OLA.csv')
    obj.reading()
    obj.matching()
    obj.plotting_mdms_age()