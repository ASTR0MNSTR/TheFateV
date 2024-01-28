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
from __algo__new import *
from __reader__ import *
    
def outflow_wAGN(SFR, LAGN, MS):
    try: 
        return 1.14*math.log10(0.52*(10**SFR) + 0.51*LAGN) - 0.41*(MS - 11)
        # return 1.14*math.log10(0.52*(10**SFR)) - 0.41*(MS - 11)
    except:
        return -99

def outflow_woAGN(SFR, MS):
    try: 
        return 1.14*math.log10(0.52*(10**SFR)) - 0.41*(MS - 11)
    except:
        return -99
    
class Main:
    def __init__(self, file):
        self.file = file 
        
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
        
        self.DataFrame = None
    
    def reading(self):
        source_path = r'E:/LICENSE/ProgsData/main/GAMAforOleg.txt'
        input_path = r'E:/backup/backup_BPT/GAMA_ETG_OLA.csv'
        output_path = r'E:/databases/Merged_Out.csv'
        merge_phys_databases(source_path, input_path, output_path)
        self.DataFrame = pd.read_csv(output_path)
        
        OUTFLOW = []
        OUTFLOW_up = []
        OUTFLOW_down = []
        
        for i in range(len(self.DataFrame['SFR_0_1Gyr_percentile50'])):
                
            if self.DataFrame['LAGN_er'][i] in ['upNon', 'upAbs']:
                OUTFLOW.append(outflow_wAGN(self.DataFrame['SFR_0_1Gyr_percentile84'][i], self.DataFrame['LAGN'][i], self.DataFrame['mass_stellar_percentile16'][i]))
                OUTFLOW_up.append(self.DataFrame['LAGN_er'][i])
                OUTFLOW_down.append(self.DataFrame['LAGN_er'][i])
            else:
                OUTFLOW.append(outflow_wAGN(self.DataFrame['SFR_0_1Gyr_percentile50'][i], self.DataFrame['LAGN'][i], self.DataFrame['mass_stellar_percentile50'][i]))
                OUTFLOW_up.append(outflow_wAGN(self.DataFrame['SFR_0_1Gyr_percentile84'][i], self.DataFrame['LAGN'][i] + float(self.DataFrame['LAGN_er'][i]), self.DataFrame['mass_stellar_percentile16'][i]))
                OUTFLOW_down.append(outflow_wAGN(self.DataFrame['SFR_0_1Gyr_percentile16'][i], self.DataFrame['LAGN'][i] - float(self.DataFrame['LAGN_er'][i]), self.DataFrame['mass_stellar_percentile84'][i]))
            
            # OUTFLOW.append(outflow_woAGN(self.DataFrame['SFR_0_1Gyr_percentile50'][i], self.DataFrame['mass_stellar_percentile50'][i]))
            # OUTFLOW_up.append(outflow_woAGN(self.DataFrame['SFR_0_1Gyr_percentile84'][i], self.DataFrame['mass_stellar_percentile16'][i]))
            # OUTFLOW_down.append(outflow_woAGN(self.DataFrame['SFR_0_1Gyr_percentile16'][i], self.DataFrame['mass_stellar_percentile84'][i]))
            
            if str(self.DataFrame['SPECID_x'][i]) == '1031459319616399360':
                print(self.DataFrame['SFR_0_1Gyr_percentile50'][i])
                print(self.DataFrame['mass_stellar_percentile50'][i])
                print(self.DataFrame['LAGN'][i])
                print(OUTFLOW_up[-1])
                print(OUTFLOW[-1])
                print(OUTFLOW_down[-1])
        
        self.DataFrame['OUTFLOW'] = OUTFLOW
        self.DataFrame['OUTFLOW_up'] = OUTFLOW_up
        self.DataFrame['OUTFLOW_down'] = OUTFLOW_down
    
    def plotting(self):
        gs_top = plt.GridSpec(1, 2, wspace=0)
        self.fig1 = plt.figure(figsize=(12, 6), tight_layout=True)

        self.ax4 = self.fig1.add_subplot(gs_top[:,0])
        self.ax5 = self.fig1.add_subplot(gs_top[:,1], sharey=self.ax4)

        self.topaxes = [self.ax5, self.ax4]

        self.ax4.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, right=True, direction='in')
        self.ax5.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, left=True, labelleft=False, right=True, labelright=False, direction='in')

        self.ax4.set_ylabel(r'$log(M_{H_2}/(M_\odot/yr))$')
        for ax in self.topaxes:    
            ax.set_xlim([8.7, 10.0])
            ax.set_ylim([-4, 10])
            ax.set_xlabel(r'$log(age/yr)$')
            #ax.set_yticks(np.arange(15, 25.9, 2))

        Main.plotter(self, 'ager_percentile50', 'OUTFLOW', 'OUTFLOW_up', 'OUTFLOW_down', 'BPT', 'WHAN', True)
        
        self.fig1.savefig('./FIGURES_IN_PAPER/OUTFLOW_SN2_w_fullDC.pdf')
        #plt.show()

    def plotter(self, x, y, up, down, BPT_key, WHAN_key, bids_key):
        for i in range(len(self.DataFrame[x])):
            AGN = self.DataFrame[BPT_key][i]
            WHAN = self.DataFrame[WHAN_key][i]
            X = self.DataFrame[x][i]
            Y = self.DataFrame[y][i]
            Y_up = self.DataFrame[up][i]
            Y_down = self.DataFrame[down][i]
            if Y_up == 'upNon':
                self.ax4.arrow(X, Y, 0, -0.1, head_width=0.01, head_length=0.03, color=self.color_dict_BPT[AGN][0], alpha=0.5)
                self.ax5.arrow(X, Y, 0, -0.1, head_width=0.01, head_length=0.03, color=self.color_dict_WHAN[WHAN][0], alpha=0.5)
            elif Y_up == 'upAbs':
                self.ax4.arrow(X, Y, 0, -0.1, head_width=0.01, head_length=0.03, color=self.color_dict_BPT[AGN][0], alpha=0.1)
                self.ax5.arrow(X, Y, 0, -0.1, head_width=0.01, head_length=0.03, color=self.color_dict_WHAN[WHAN][0], alpha=0.1)
            else:
                self.ax4.errorbar(X, Y, yerr = [[Y - float(Y_down)], [float(Y_up) - Y]], alpha = 0.5, color = self.color_dict_BPT[AGN][0], marker = '.')
                self.ax5.errorbar(X, Y, yerr = [[Y - float(Y_down)], [float(Y_up) - Y]], alpha = 0.5, color = self.color_dict_WHAN[WHAN][0], marker = '.')
                
        self.ax4.axhline(y = 0, color = 'black', linestyle='dashed')        
        self.ax5.axhline(y = 0, color = 'black', linestyle='dashed')
               
        
        # class_list_BPT = class_list_creator_w_err(XX, YY, Y_up, Y_down, BPT_keys, 'BPT')
        # class_list_WHAN = class_list_creator_w_err(XX, YY, Y_up, Y_down, WHAN_keys, 'WHAN')

        # self.means = []
        # self.errs = []
        # self.ages = []
        # if bids_key == True:
        #     bids = [[8.8, 9.0], [9.0, 9.2], [9.2, 9.4], [9.4, 9.6], [9.6, 9.8], [9.8, 10.0]]
        #     ages_const = np.arange(8.8, 10.2, 0.2)
        # else:
        #     bids = [[10.0, 10.25], [10.25, 10.5], [10.5, 10.75], [10.75, 11], [11, 11.25], [11.25, 11.5]]
        #     ages_const = np.arange(7, 12, 1)

        # errs = []
        # means = []
        # for item in class_list_BPT:
        #     X_plot = []
        #     Y_plot = []
        #     X, Y, err, length = monte_carlo(item[0], item[1], item[2], item[3], bids)
        #     means.append(Y)
        #     errs.append(err)
        #     for j in range(len(X)):
        #         if Y[j] != -99:
        #             self.ax4.errorbar(X[j], Y[j], alpha = 1, xerr=0, yerr= err[j], color=item[4][0], fmt=item[4][1], ms = 12)
        #             self.ax4.text(X[j], Y[j], length[j], c = 'red')
        #             X_plot.append(X[j])
        #             Y_plot.append(Y[j])
        #     self.ax4.plot(X_plot, Y_plot, alpha = 1, color=item[4][0])

        #for key in self.color_dict_leg.keys():
        #        self.ax4.scatter(-99, -99, alpha= 1, color = self.color_dict_leg[key][0], marker = self.color_dict_leg[key][2], s = self.color_dict_leg[key][1], label=key)
        # self.ax4.set_xticks(ages_const)
        # for j, item in enumerate(class_list_BPT):
        #     self.ax4.scatter(-99, -99, alpha = 1, color=item[4][0], marker=item[4][1], s = 150, label=list_names_BPT_1[j])
        # self.ax4.legend(loc=3, fontsize='13')
        
        # errs = []
        # means = []
        # for item in class_list_WHAN:
        #     X_plot = []
        #     Y_plot = []
        #     X, Y, err, length = monte_carlo(item[0], item[1], item[2], item[3], bids)
        #     means.append(Y)
        #     errs.append(err)
        #     for j in range(len(X)):
        #         if Y[j] != -99:
        #             self.ax5.errorbar(X[j], Y[j], alpha = 1, xerr=0, yerr= err[j], color=item[4][0], fmt=item[4][1], ms = 12)
        #             self.ax5.text(X[j], Y[j], length[j], c = 'red')
        #             X_plot.append(X[j])
        #             Y_plot.append(Y[j])
        #     self.ax5.plot(X_plot, Y_plot, alpha = 1, color=item[4][0])

        #for key in self.color_dict_leg.keys():
        #        self.ax4.scatter(-99, -99, alpha= 1, color = self.color_dict_leg[key][0], marker = self.color_dict_leg[key][2], s = self.color_dict_leg[key][1], label=key)
        # self.ax5.set_xticks(ages_const)
        # for j, item in enumerate(class_list_WHAN):
        #     self.ax5.scatter(-99, -99, alpha = 1, color=item[4][0], marker=item[4][1], s = 150, label=list_names_WHAN[j])
        # self.ax5.legend(loc=3, fontsize='13')

        #self.ax4.set_xlim(8, 10.1)
        #self.ax4.set_ylim(14, 26)
        #self.fig.savefig(name_file + '.pdf')
        plt.show()

if __name__ == '__main__':
    obj = Main('GAMA_ETG_OLA.csv')
    obj.reading()
    obj.plotting()