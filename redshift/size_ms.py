import csv
import matplotlib.pyplot as plt 
import numpy as np
import math
import statistics as st
from scipy.stats import pearsonr
from scipy.stats import sem
import pandas as pd
from scipy.optimize import curve_fit
import random
from _global_ import *

#redshift -stellar mass plot with bins!  - anomaly or interesting conclusion

class hp:
    def log_er(item):
        return [[abs(math.log(1-item, 10))], [abs(math.log(1+item, 10))]]

    def temp_stats_1(x, y_mid, y_up, y_down, x_bids):
        y_values = [[], [], [], [], [], []]
        stmeaner = []
        stmean = []
        medians = [(pair[0] + pair[1])/2 for pair in x_bids]

        for j, item in enumerate(x):
            for i, pair in enumerate(x_bids):
                if item >= pair[0] and item < pair[1]:
                    y_values[i].append([y_mid[j] - y_down[j], y_up[j] - y_mid[j], y_mid[j]])
                    break
        
        length = [len(item) for item in y_values]

        for mean in y_values:
            if len(mean) <= 3:
                stmean.append(-99)
                stmeaner.append(0)
            else:
                data = []
                for i in range(100):
                    av = []
                    for pair in mean:
                        if random.random() > 0.5:
                            av.append(pair[2] + abs(random.gauss(0, pair[1])))
                        else:
                            av.append(pair[2] - abs(random.gauss(0, pair[0])))
                    data.append(st.mean(av))
                
                data.sort()
                stmean.append(data[49])
                stmeaner.append([[data[49] - data[15]], [data[83] - data[49]]])


        return medians, stmean, stmeaner, length
         
    
class Main(hp):
    def __init__(self, ola_file, out):
        self.ola_file = ola_file
        self.out = out
        self.base_dict = {}
        self.data_dict = []
        self.color_dict_BPT = color_dict
        self.color_dict_WHAN = cd_WHAN
        self.color_dict_leg = color_dict_leg

        self.list_names_BPT = ['AGN', 'UNC', 'SF', 'NOEL']
        self.list_names_WHAN = ['sAGN', 'ELR', 'SF', 'RG', 'LLR', 'wAGN']

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
            X = float(line.split()[1])
            Y = 3/float(line.split()[301])
            Y_up = 3/float(line.split()[301])
            Y_down = 3/float(line.split()[301])
            #SFR_up = float(line.split()[96])
            #SFR_down = float(line.split()[94])
            #SFR_er = [[abs(SFR_down-SFR_0)], [abs(SFR_up-SFR_0)]]
            self.base_dict.update({GAMAID: [X, Y, Y_up, Y_down]})

        self.data_dict = []

        with open(self.out, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                self.data_dict.append({'GAMAID': int(row['GAMAID']), 'AGN': row['BPT'], 'SC_WHAN' : row['WHAN'], 'BMS' : int(row['BMS'])})

    def matching(self):
        for item in self.data_dict:
            try:
                X, Y, Y_up, Y_down = self.base_dict[item['GAMAID']]
                item.update({'X' : X, 'Y' : Y, 'Y_up' : Y_up, 'Y_down' : Y_down})
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

        self.ax4.set_ylabel(r'$Aperture in R_{eff} \; i, \;$')
        for ax in self.topaxes:    
            ax.set_xlim([0, 0.35])

        #self.ax5.set_xlabel(r'$log(age/yr)$')

        Main.plotter_BPT(self, 'X', 'Y', 'Y_up', 'Y_down', 'AGN', True)
        Main.plotter_WHAN(self, 'X', 'Y', 'Y_up', 'Y_down', 'SC_WHAN', True)
        self.fig1.savefig('./FIGURES/AP_RED.pdf')
        plt.show()

    def plotter_BPT(self, x, y, up, down, AGN_key, bids):

        yes_temp = []
        UNC_temp = []
        no_temp = []
        noel_temp = []

        yes_temp_age = []
        UNC_temp_age = []
        no_temp_age = []
        noel_temp_age = []

        yes_temp_up = []
        UNC_temp_up = []
        no_temp_up = []
        noel_temp_up = []

        yes_temp_down = []
        UNC_temp_down = []
        no_temp_down = []
        noel_temp_down = []

        tot_age = []
        tot = []
        tot_up = []
        tot_down = []

        for item in self.data_dict:
            X = item[x]
            Y = item[y]
            Y_up = item[up]
            Y_down = item[down]
            AGN = item[AGN_key]
            self.ax4.scatter(X, Y, alpha= 0.4, color = self.color_dict_BPT[AGN][0], marker='.', s = self.color_dict_BPT[AGN][1])
            tot_age.append(X)
            tot.append(Y)
            tot_up.append(Y_up)
            tot_down.append(Y_down)
            #marker = self.color_dict[AGN][2]
            if AGN in ['AGNXY', 'AGNX', 'AGNY']:
                    yes_temp.append(Y)
                    yes_temp_age.append(X)
                    yes_temp_up.append(Y_up)
                    yes_temp_down.append(Y_down)
            elif AGN in ['UNCXY', 'UNCX', 'UNCY']:
                    UNC_temp.append(Y)
                    UNC_temp_age.append(X)
                    UNC_temp_up.append(Y_up)
                    UNC_temp_down.append(Y_down)
            elif AGN in ['SFXY', 'SFX', 'SFY']:
                    no_temp.append(Y)
                    no_temp_age.append(X)
                    no_temp_up.append(Y_up)
                    no_temp_down.append(Y_down)
            elif AGN == 'NOEL':
                    noel_temp.append(Y)
                    noel_temp_age.append(X)
                    noel_temp_up.append(Y_up)
                    noel_temp_down.append(Y_down)
        
        class_list = [[yes_temp_age, yes_temp, 'midnightblue', yes_temp_down, yes_temp_up, 'P'], [UNC_temp_age, UNC_temp, 'springgreen', UNC_temp_down, UNC_temp_up, 'H'], [no_temp_age, no_temp, 'mediumvioletred', no_temp_down, no_temp_up, '*'], [noel_temp_age, noel_temp, 'orchid', noel_temp_down, noel_temp_up, 'o']]

        self.means = []
        self.errs = []
        self.ages = []
        
        age_bids = [[0.0, 0.05], [0.05, 0.1], [0.1, 0.15], [0.15, 0.2], [0.2, 0.25], [0.25, 0.34]]
        ages_const = np.arange(0.0, 0.35, 0.05)

        for item in class_list:
            X_plot = []
            Y_plot = []
            X, Y, err, length = hp.temp_stats_1(item[0], item[1], item[3], item[4], age_bids)
            self.ages = X
            self.means.append(Y)
            self.errs.append(err)
            for j in range(len(X)):
                if Y[j] != -99:
                    self.ax4.errorbar(X[j], Y[j], alpha = 1, xerr=0, yerr= err[j], color=item[2], fmt=item[5], ms = 12)
                    self.ax4.text(X[j], Y[j], length[j], c = 'red')
                    X_plot.append(X[j])
                    Y_plot.append(Y[j])
            self.ax4.plot(X_plot, Y_plot, alpha = 1, color=item[2])

        for key in self.color_dict_leg.keys():
                self.ax4.scatter(-99, -99, alpha= 1, color = self.color_dict_leg[key][0], marker = self.color_dict_leg[key][2], s = self.color_dict_leg[key][1], label=key)
        self.ax4.set_xticks(ages_const)
        for j, item in enumerate(class_list):
            self.ax4.scatter(-99, -99, alpha = 1, color=item[2], marker=item[5], s = 150, label=self.list_names_BPT[j])
        self.ax4.legend()

        #self.ax4.set_xlim(8, 10.1)
        #self.ax4.set_ylim(14, 26)
        #self.fig.savefig(name_file + '.pdf')
        #plt.show()
    
    def plotter_WHAN(self, x, y, up, down, AGN_key, bids):

        
        tot_age = []
        tot = []
        tot_up = []
        tot_down = []

        yes_temp = []
        UNC_temp = []
        no_temp = []
        noel_temp = []

        yes_temp_age = []
        UNC_temp_age = []
        no_temp_age = []
        noel_temp_age = []

        yes_temp_up = []
        UNC_temp_up = []
        no_temp_up = []
        noel_temp_up = []

        yes_temp_down = []
        UNC_temp_down = []
        no_temp_down = []
        noel_temp_down = []

        wAGN_temp = []
        wAGN_temp_age = []
        wAGN_temp_down = []
        wAGN_temp_up = []

        llr_temp = []
        llr_temp_age = []
        llr_temp_down = []
        llr_temp_up = []


        for item in self.data_dict:
            X = item[x]
            Y = item[y]
            Y_up = item[up]
            Y_down = item[down]
            AGN = item[AGN_key]
            self.ax5.scatter(X, Y, alpha= 0.4, color = self.color_dict_WHAN[AGN][0], marker='.', s = self.color_dict_WHAN[AGN][1])
            tot_age.append(X)
            tot.append(Y)
            tot_up.append(Y_up)
            tot_down.append(Y_down)
            #marker = self.color_dict[AGN][2]
            if AGN in ['sAGN']:
                    yes_temp.append(Y)
                    yes_temp_age.append(X)
                    yes_temp_up.append(Y_up)
                    yes_temp_down.append(Y_down)
            elif AGN in ['ELR']:
                    UNC_temp.append(Y)
                    UNC_temp_age.append(X)
                    UNC_temp_up.append(Y_up)
                    UNC_temp_down.append(Y_down)
            elif AGN in ['SF']:
                    no_temp.append(Y)
                    no_temp_age.append(X)
                    no_temp_up.append(Y_up)
                    no_temp_down.append(Y_down)
            elif AGN in ['RG']:
                    noel_temp.append(Y)
                    noel_temp_age.append(X)
                    noel_temp_up.append(Y_up)
                    noel_temp_down.append(Y_down)
            elif AGN in ['LLR']:
                    llr_temp.append(Y)
                    llr_temp_age.append(X)
                    llr_temp_up.append(Y_up)
                    llr_temp_down.append(Y_down)
            elif AGN in ['wAGN']:
                    wAGN_temp.append(Y)
                    wAGN_temp_age.append(X)
                    wAGN_temp_up.append(Y_up)
                    wAGN_temp_down.append(Y_down)
        
        #class_list = [[yes_temp_age, yes_temp, 'midnightblue'], [UNC_temp_age, UNC_temp, 'red'], [no_temp_age, no_temp, 'mediumvioletred'], [noel_temp_age, noel_temp, 'orchid'], [llr_age, llr, 'maroon']]
        class_list = [[yes_temp_age, yes_temp, 'midnightblue', yes_temp_down, yes_temp_up, 'P'], [UNC_temp_age, UNC_temp, 'sandybrown', UNC_temp_down, UNC_temp_up, 'D'], [no_temp_age, no_temp, 'mediumvioletred', no_temp_down, no_temp_up, '*'], [noel_temp_age, noel_temp, 'chocolate', noel_temp_down, noel_temp_up, 'o'], [llr_temp_age, llr_temp, 'maroon', llr_temp_down, llr_temp_up, 'o'], [wAGN_temp_age, wAGN_temp, 'blue', wAGN_temp_down, wAGN_temp_up, 'P']]

        self.means = []
        self.errs = []
        self.ages = []
    
        age_bids = [[0.0, 0.05], [0.05, 0.1], [0.1, 0.15], [0.15, 0.2], [0.2, 0.25], [0.25, 0.34]]
        ages_const = np.arange(0.0, 0.35, 0.05)

        for item in class_list:
            X_plot = []
            Y_plot = []
            X, Y, err, length = hp.temp_stats_1(item[0], item[1], item[3], item[4], age_bids)
            self.ages = X
            self.means.append(Y)
            self.errs.append(err)
            for j in range(len(X)):
                if Y[j] != -99:
                    self.ax5.errorbar(X[j], Y[j], alpha = 1, xerr=0, yerr= err[j], color=item[2], fmt=item[5], ms = 12)
                    self.ax5.text(X[j], Y[j], length[j], c = 'red')
                    X_plot.append(X[j])
                    Y_plot.append(Y[j])
            self.ax5.plot(X_plot, Y_plot, alpha = 1, color=item[2])

        for key in self.color_dict_WHAN.keys():
            self.ax5.scatter(-99, -99, alpha= 0.7, color = self.color_dict_WHAN[key][0], marker = self.color_dict_WHAN[key][2], s = self.color_dict_WHAN[key][1], label=key)
        self.ax5.set_xticks(ages_const)
        for j, item in enumerate(class_list):
            self.ax5.scatter(-99, -99, alpha = 1, color=item[2], marker=item[5], s = 150, label=self.list_names_WHAN[j])
        self.ax5.legend()

        #self.ax4.set_xlim(8, 10.1)
        #self.ax4.set_ylim(14, 26)
        #self.fig.savefig(name_file + '.pdf')
        #plt.show()

if __name__ == '__main__':
    obj = Main('E:\LICENSE\ProgsData\main\GAMAforOleg.txt', 'GAMA_ETG_OLA.csv')
    obj.reading()
    obj.matching()
    obj.plotting_mdms_age()