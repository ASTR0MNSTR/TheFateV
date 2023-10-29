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
        
        for mean in y_values:
            if len(mean) <= 2:
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


        return medians, stmean, stmeaner
         
    
class Main(hp):
    def __init__(self, ola_file, table, out):
        self.ola_file = ola_file
        self.out = out
        self.table = table
        self.base_dict = {}
        self.data_dict = []
        self.color_dict_BPT ={
            'AGN' : ['midnightblue', 6, '.'],
            'AGNX' : ['midnightblue', 6, '.'],
            'AGNY' : ['midnightblue', 6, '.'],
            'UNC' : ['springgreen', 6, '.'],
            'UNCX' : ['springgreen', 6, '.'],
            'UNCY' : ['springgreen', 6, '.'],
            'SF' : ['mediumvioletred', 6, '.'],
            'SFX' : ['mediumvioletred', 6, '.'],
            'SFY' : ['mediumvioletred', 6, '.'],
            'NOEL' : ['orchid', 6, '.'],
            'NDA' : ['yellow', 6, '.']
        }

        self.color_dict_WHAN = {
            'UNC' : ['springgreen', 9, '.'],
            'NDA' : ['gold', 9, '.'],
            'LLR' : ['maroon', 9, '.'],
            'ELR' : ['red', 9, '.'],
            'SF' : ['mediumvioletred', 9, '.'],
            'wAGN' : ['blue', 15, '.'],
            'sAGN' : ['midnightblue', 15, '.'],
            'RG' : ['chocolate', 9, '.']
        }

        self.color_dict_leg ={
            'AGN' : ['midnightblue', 6, '.'],
            'UNC' : ['springgreen', 6, '.'],
            'SF' : ['mediumvioletred', 6, '.'],
            'NOEL' : ['orchid', 6, '.'],
            'NDA' : ['yellow', 6, '.']
        }

        self.list_names_BPT = ['AGN', 'UNC', 'SF', 'NOEL']
        self.list_names_WHAN = ['sAGN + wAGN', 'ELR + RG', 'SF', 'LLR']

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
            X = float(line.split()[35])
            if float(line.split()[696])/3631 <= 0:
                U = -99
            else:
                U = -2.5*math.log10(float(line.split()[696])/3631)

            if float(line.split()[700])/3631 <= 0:
                R = 0
            else:
                R = -2.5*math.log10(float(line.split()[700])/3631)
            Y = U - R
            Y_up = U - R + 0.02
            Y_down = U - R - 0.02
            #SFR_up = float(line.split()[96])
            #SFR_down = float(line.split()[94])
            #SFR_er = [[abs(SFR_down-SFR_0)], [abs(SFR_up-SFR_0)]]
            self.base_dict.update({GAMAID: [X, Y, Y_up, Y_down]})

        self.data_dict = []

        with open(self.out, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                self.data_dict.append({'OI': int(row['0I']), 'GAMAID': int(row['GAMAID']), 'age': float(row['age']), 'AGN': row['AGN'], 'SC_WHAN' : row['SC_WHAN'], 'BMS' : int(row['BMS'])})

    def matching(self):
        for item in self.data_dict:
            try:
                X, Y, Y_up, Y_down = self.base_dict[item['GAMAID']]
                item.update({'X' : X, 'Y' : Y, 'Y_up' : Y_up, 'Y_down' : Y_down})
            except: 
                print(item['GAMAID'])
    
    def plotting_mdms_age(self):
        gs_top = plt.GridSpec(2, 1, hspace=0)
        self.fig1 = plt.figure(figsize=(10,8))

        self.ax4 = self.fig1.add_subplot(gs_top[0,:])
        self.ax5 = self.fig1.add_subplot(gs_top[1,:], sharex=self.ax4)

        self.topaxes = [self.ax5, self.ax4]

        self.ax5.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, right=True, direction='in')
        self.ax4.tick_params(top=True, labeltop=False, bottom=True, labelbottom=False, right=True, direction='in')

        for ax in self.topaxes:    
            ax.set_ylabel(r'$u - r, \; colour$')
            ax.set_xlim([9, 12])
            ax.set_ylim([0, 6])
            ax.set_yticks(np.arange(0, 6, 0.5))

        self.ax5.set_xlabel(r'$M_s, \; dex(M_{Sun})$')

        Main.plotter_BPT(self, 'X', 'Y', 'Y_up', 'Y_down', 'AGN', False)
        Main.plotter_WHAN(self, 'X', 'Y', 'Y_up', 'Y_down', 'SC_WHAN', False)
        self.fig1.savefig('green.pdf')
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
            if Y > -20:
                Y_up = item[up]
                Y_down = item[down]
                AGN = item[AGN_key]
                self.ax4.scatter(X, Y, alpha= 0.4, color = self.color_dict_BPT[AGN][0], marker='.', s = self.color_dict_BPT[AGN][1])
                tot_age.append(X)
                tot.append(Y)
                tot_up.append(Y_up)
                tot_down.append(Y_down)
            #marker = self.color_dict[AGN][2]
                if AGN in ['AGN', 'AGNX', 'AGNY']:
                    yes_temp.append(Y)
                    yes_temp_age.append(X)
                    yes_temp_up.append(Y_up)
                    yes_temp_down.append(Y_down)
                elif AGN in ['UNC', 'UNCX', 'UNCY']:
                    UNC_temp.append(Y)
                    UNC_temp_age.append(X)
                    UNC_temp_up.append(Y_up)
                    UNC_temp_down.append(Y_down)
                elif AGN in ['SF', 'SFX', 'SFY']:
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
        if bids == True:
            age_bids = [[8.8, 9.0], [9.0, 9.2], [9.2, 9.4], [9.4, 9.6], [9.6, 9.8], [9.8, 10.0]]
            ages_const = np.arange(8, 10.2, 0.2)
        else:
            age_bids = [[10.0, 10.25], [10.25, 10.5], [10.5, 10.75], [10.75, 11], [11, 11.25], [11.25, 11.5]]
            ages_const = np.arange(9, 12, 1)

        for item in class_list:
            X_plot = []
            Y_plot = []
            X, Y, err = hp.temp_stats_1(item[0], item[1], item[3], item[4], age_bids)
            self.ages = X
            self.means.append(Y)
            self.errs.append(err)
            for j in range(len(X)):
                if Y[j] != -99:
                    self.ax4.errorbar(X[j], Y[j], alpha = 1, xerr=0, yerr= err[j], color=item[2], fmt=item[5], ms = 12)
                    X_plot.append(X[j])
                    Y_plot.append(Y[j])
            self.ax4.plot(X_plot, Y_plot, alpha = 1, color=item[2])

        x = np.arange(9, 13, 0.1)
        self.ax4.plot(x, -0.24 + 0.25*x, color='olive')
        self.ax4.plot(x, -0.75 + 0.25*x, color='olive')

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

        for item in self.data_dict:
            X = item[x]
            Y = item[y]
            if Y > -20:
                Y_up = item[up]
                Y_down = item[down]
                AGN = item[AGN_key]
                self.ax5.scatter(X, Y, alpha= 0.4, color = self.color_dict_WHAN[AGN][0], marker='.', s = self.color_dict_WHAN[AGN][1])
                tot_age.append(X)
                tot.append(Y)
                tot_up.append(Y_up)
                tot_down.append(Y_down)
                #marker = self.color_dict[AGN][2]
                if AGN in ['sAGN', 'wAGN']:
                    yes_temp.append(Y)
                    yes_temp_age.append(X)
                    yes_temp_up.append(Y_up)
                    yes_temp_down.append(Y_down)
                elif AGN in ['ELR', 'RG']:
                    UNC_temp.append(Y)
                    UNC_temp_age.append(X)
                    UNC_temp_up.append(Y_up)
                    UNC_temp_down.append(Y_down)
                elif AGN in ['SF']:
                    no_temp.append(Y)
                    no_temp_age.append(X)
                    no_temp_up.append(Y_up)
                    no_temp_down.append(Y_down)
                elif AGN in ['LLR']:
                    noel_temp.append(Y)
                    noel_temp_age.append(X)
                    noel_temp_up.append(Y_up)
                    noel_temp_down.append(Y_down)
        
        #class_list = [[yes_temp_age, yes_temp, 'midnightblue'], [UNC_temp_age, UNC_temp, 'red'], [no_temp_age, no_temp, 'mediumvioletred'], [noel_temp_age, noel_temp, 'orchid'], [llr_age, llr, 'maroon']]
        class_list = [[yes_temp_age, yes_temp, 'midnightblue', yes_temp_down, yes_temp_up, 'P'], [UNC_temp_age, UNC_temp, 'red', UNC_temp_down, UNC_temp_up, 'D'], [no_temp_age, no_temp, 'mediumvioletred', no_temp_down, no_temp_up, '*'], [noel_temp_age, noel_temp, 'maroon', noel_temp_down, noel_temp_up, 'o']]

        self.means = []
        self.errs = []
        self.ages = []
        if bids == True:
            age_bids = [[8.8, 9.0], [9.0, 9.2], [9.2, 9.4], [9.4, 9.6], [9.6, 9.8], [9.8, 10.0]]
            ages_const = np.arange(8, 10.2, 0.2)
        else:
            age_bids = [[10.0, 10.25], [10.25, 10.5], [10.5, 10.75], [10.75, 11], [11, 11.25], [11.25, 11.5]]
            ages_const = np.arange(9, 12, 1)

        for item in class_list:
            X_plot = []
            Y_plot = []
            X, Y, err = hp.temp_stats_1(item[0], item[1], item[3], item[4], age_bids)
            self.ages = X
            self.means.append(Y)
            self.errs.append(err)
            for j in range(len(X)):
                if Y[j] != -99:
                    self.ax5.errorbar(X[j], Y[j], alpha = 1, xerr=0, yerr= err[j], color=item[2], fmt=item[5], ms = 12)
                    X_plot.append(X[j])
                    Y_plot.append(Y[j])
            self.ax5.plot(X_plot, Y_plot, alpha = 1, color=item[2])

        x = np.arange(4, 13, 0.1)
        self.ax5.plot(x, -0.24 + 0.25*x, color='olive')
        self.ax5.plot(x, -0.75 + 0.25*x, color='olive')

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
    obj = Main('GAMAforOleg.txt', 'sample_out_spectra.csv', 'GAMA_ETG_OLA.csv')
    obj.reading()
    obj.matching()
    obj.plotting_mdms_age()