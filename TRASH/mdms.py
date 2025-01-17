import csv
import matplotlib.pyplot as plt 
import numpy as np
import math
# import statistics as st
import pandas as pd
# from scipy.optimize import curve_fit
from __legpars__ import *
from __stats__ import *

class hp:
    
    def fit(x, y_mid, y_up, y_down):
        n_sim = 50
        a_list= []
        b_list = []
        b_err_list = []
        for i in range(1000):
            data = []
            x_fit = []
            for j in range(len(x)):  
                pair = [y_mid[j] - y_down[j], y_up[j] - y_mid[j], y_mid[j]]
                if random.random() > 0.5:
                    data.append(pair[2] + abs(random.gauss(0, pair[1])))
                else:
                    data.append(pair[2] - abs(random.gauss(0, pair[0])))

                x_fit.append(10**x[j])

            x_fit = np.array(x_fit)
            y_fit = np.array(data)
            popt, pcov = curve_fit(lambda age, a, b: a - 0.4343*age/b, x_fit, y_fit)

            a_list.append(popt[0])
            b_list.append(popt[1])
            b_err_list.append(np.sqrt(np.diag(pcov))[1])
        
        a_list.sort()
        b_list.sort()
        a = a_list[499]
        a_err = [[a_list[499] - a_list[159]], [a_list[839] - a_list[499]]]
        b = b_list[499]
        print(b/10**9)
        b_err = [[b_list[499] - b_list[159]], [b_list[839] - b_list[499]]]
        print((b_list[499] - b_list[0]), (b_list[-1] - b_list[499]))
        print(b_err[0][0]/10**9, b_err[1][0]/10**9)
        print(b_err_list[499]/10**9)
        
        print('\n')
        return a, b
         
    
class Main(hp):
    def __init__(self, ola_file, out):
        self.ola_file = ola_file
        self.out = out
        self.base_dict = {}
        self.data_dict = []
        self.list_names_BPT = ['AGN', 'UNC', 'SF', 'NOEL']
        self.list_names_WHAN = ['sAGN', 'wAGN', 'SF','ELR', 'RG', 'LLR']

        self.color_dict_BPT = color_dict_BPT
        self.color_dict_WHAN = cd_WHAN
        self.color_dict_leg = color_dict_leg

        self.BMS_dict= {
            0 : ['o', 12], #bms
            1 : ['*', 12], #ms
            -1 : ['x', 15] #??????
        }
        self.a = 0
        self.b = 0
    
    def reading(self):
        with open(self.ola_file, 'r') as input:
            lines = input.readlines()
        lines_stripped = [line.strip() for line in lines]
        for line in lines_stripped:
            GAMAID = int(line.split()[0])
            X = float(line.split()[125])
            MS = float(line.split()[35])
            MDMS = float(line.split()[89]) - float(line.split()[35])
            
            MS_up = float(line.split()[36]) - float(line.split()[35])
            MS = float(line.split()[35])
            MS_down = float(line.split()[35]) - float(line.split()[34])

            MD_up = float(line.split()[90]) - float(line.split()[89])
            MD = float(line.split()[89])
            MD_down = float(line.split()[89]) - float(line.split()[88])

            sigma_up = (MD_up**2 + MS_down**2)**0.5
            sigma_down = (MD_down**2 + MS_up**2)**0.5

            MDMS_up = MDMS + sigma_up
            MDMS_down = MDMS - sigma_down
            #SFR_up = float(line.split()[96])
            #SFR_down = float(line.split()[94])
            #SFR_er = [[abs(SFR_down-SFR_0)], [abs(SFR_up-SFR_0)]]
            self.base_dict.update({GAMAID: [X, MDMS, MDMS_up, MDMS_down]})

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

        self.ax4.set_ylabel(r'$log(M_{d}/M_{*})$')
        for ax in self.topaxes:    
            ax.set_xlim([7.9, 10.1])
            ax.set_ylim([-5, -0.1])

        self.ax5.set_xlabel(r'$log(age)$')

        Main.plotter_BPT(self, 'X', 'Y', 'Y_up', 'Y_down', 'AGN', True)
        Main.plotter_WHAN(self, 'X', 'Y', 'Y_up', 'Y_down', 'SC_WHAN', True)
        self.fig1.savefig('./FIGURES/MDMS.pdf')
        #plt.show()

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
        if bids == True:
            age_bids = [[8.8, 9.0], [9.0, 9.2], [9.2, 9.4], [9.4, 9.6], [9.6, 9.8], [9.8, 10.0]]
            ages_const = np.arange(8.8, 10.2, 0.2)
        else:
            age_bids = [[10.0, 10.25], [10.25, 10.5], [10.5, 10.75], [10.75, 11], [11, 11.25], [11.25, 11.5]]
            ages_const = np.arange(7, 12, 1)

        for item in class_list:
            X_plot = []
            Y_plot = []
            X, Y, err, length = bootstrapper(item[0], item[1], item[3], item[4], age_bids)
            self.ages = X
            self.means.append(Y)
            self.errs.append(err)
            for j in range(len(X)):
                if Y[j] != -99:
                    self.ax4.errorbar(X[j], Y[j], alpha = 1, xerr=0, yerr= err[j], color=item[2], fmt=item[5], ms = 12)
                    #self.ax4.text(X[j], Y[j], length[j], c = 'red')
                    X_plot.append(X[j])
                    Y_plot.append(Y[j])
            self.ax4.plot(X_plot, Y_plot, alpha = 1, color=item[2])

            #a, b = hp.fit(item[0], item[1], item[3], item[4])
            #x = np.arange(7.9, 10.1, 0.01)
            #a = item[6]
            #b = item[7]
            #self.ax4.plot(x, ((-(10**x)/10**9) / b) - a, alpha = 1, color=item[2], linestyle='dashed')

            #curve_fit
        
        x = np.arange(7.9, 10.1, 0.01)
        a_all = 2.3152363967664957
        b_all = 5.276003555081468
        self.ax4.plot(x, ((-(10**x)/10**9) / b_all) - a_all, alpha = 1, color='k')

        #for key in self.color_dict_leg.keys():
        #        self.ax4.scatter(-99, -99, alpha= 1, color = self.color_dict_leg[key][0], marker = self.color_dict_leg[key][2], s = self.color_dict_leg[key][1], label=key)
        self.ax4.set_xticks(ages_const)
        for j, item in enumerate(class_list):
            self.ax4.scatter(-99, -99, alpha = 1, color=item[2], marker=item[5], s = 150, label=self.list_names_BPT[j])
        self.ax4.legend(fontsize="13")

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
        class_list = [[yes_temp_age, yes_temp, 'midnightblue', yes_temp_down, yes_temp_up, 'P'], [wAGN_temp_age, wAGN_temp, 'blue', wAGN_temp_down, wAGN_temp_up, 'P'], [no_temp_age, no_temp, 'mediumvioletred', no_temp_down, no_temp_up, '*'], [UNC_temp_age, UNC_temp, 'sandybrown', UNC_temp_down, UNC_temp_up, 'D'], [noel_temp_age, noel_temp, 'chocolate', noel_temp_down, noel_temp_up, 'o'], [llr_temp_age, llr_temp, 'maroon', llr_temp_down, llr_temp_up, 'o']]

        self.means = []
        self.errs = []
        self.ages = []
        if bids == True:
            age_bids = [[8.8, 9.0], [9.0, 9.2], [9.2, 9.4], [9.4, 9.6], [9.6, 9.8], [9.8, 10.0]]
            ages_const = np.arange(8.8, 10.2, 0.2)
        else:
            age_bids = [[10.0, 10.25], [10.25, 10.5], [10.5, 10.75], [10.75, 11], [11, 11.25], [11.25, 11.5]]
            ages_const = np.arange(7, 12, 1)

        for item in class_list:
            X_plot = []
            Y_plot = []
            X, Y, err, length = bootstrapper(item[0], item[1], item[3], item[4], age_bids)
            self.ages = X
            self.means.append(Y)
            self.errs.append(err)
            for j in range(len(X)):
                if Y[j] != -99:
                    self.ax5.errorbar(X[j], Y[j], alpha = 1, xerr=0, yerr= err[j], color=item[2], fmt=item[5], ms = 12)
                    #self.ax5.text(X[j], Y[j], length[j], c = 'red')
                    X_plot.append(X[j])
                    Y_plot.append(Y[j])
            self.ax5.plot(X_plot, Y_plot, alpha = 1, color=item[2])
            #x = np.arange(7.9, 10.1, 0.01)
            #a = item[6]
            #b = item[7]
            #self.ax5.plot(x, ((-(10**x)/10**9) / b) - a, alpha = 1, color=item[2], linestyle='dashed')

            #curve_fit
        
        x = np.arange(7.9, 10.1, 0.01)
        a_all = 2.3152363967664957
        b_all = 5.276003555081468
        self.ax5.plot(x, ((-(10**x)/10**9) / b_all) - a_all, alpha = 1, color='k')

        #for key in self.color_dict_WHAN.keys():
        #   self.ax5.scatter(-99, -99, alpha= 0.7, color = self.color_dict_WHAN[key][0], marker = self.color_dict_WHAN[key][2], s = self.color_dict_WHAN[key][1], label=key)
        self.ax5.set_xticks(ages_const)
        for j, item in enumerate(class_list):
            self.ax5.scatter(-99, -99, alpha = 1, color=item[2], marker=item[5], s = 150, label=self.list_names_WHAN[j])
        self.ax5.legend(fontsize="13")

        #self.ax4.set_xlim(8, 10.1)
        #self.ax4.set_ylim(14, 26)
        #self.fig.savefig(name_file + '.pdf')
        #plt.show()

if __name__ == '__main__':
    obj = Main(r'E:/LICENSE/ProgsData/main/GAMAforOleg.txt', r'E:/backup/backup_BPT/GAMA_ETG_OLA.csv')
    obj.reading()
    obj.matching()
    obj.plotting_mdms_age()