import csv
import matplotlib.pyplot as plt 
import numpy as np
import math
# import statistics as st
import pandas as pd
from scipy.optimize import curve_fit
from __legpars__ import *
from __stats__ import *
from __plt__ import *
from __algo__ import *
from __reader__ import *
    
def outflow_wAGN(SFR, LAGN, MS):
    try: 
        return 1.14*math.log10(0.52*(10**SFR) + 0.51*LAGN) - 0.41*math.log10((10**(MS))/(10**11))
        # return 1.14*math.log10(0.52*(10**SFR)) - 0.41*(MS - 11)
    except:
        return -99

def outflow_woAGN(SFR, MS):
    try: 
        return 1.14*math.log10(0.52*(10**SFR)) - 0.41*math.log10((10**(MS))/(10**11))
    except:
        return -99
    
def fit(x, y_mid, y_up, y_down):
    n_sim = 50
    a_list= []
    b_list = []
    b_err_list = []
    a_err_list = []
    for i in range(1001):
        data = []
        x_fit = []
        for j in range(len(x)):  
            pair = [y_mid[j] - y_down[j], y_up[j] - y_mid[j], y_mid[j]]
            if random.random() > 0.5:
                data.append(pair[2] + abs(random.gauss(0, pair[1])))
            else:
                data.append(pair[2] - abs(random.gauss(0, pair[0])))

            x_fit.append(x[j])

        x_fit = np.array(x_fit)
        y_fit = np.array(data)
        popt, pcov = curve_fit(func, x_fit, y_fit)

        a_list.append(popt[0])
        b_list.append(popt[1])
        a_err_list.append(np.sqrt(np.diag(pcov))[0])
        b_err_list.append(np.sqrt(np.diag(pcov))[1])
        
    a_list = np.asarray(a_list)    
    b_list = np.asarray(b_list)    
    a_err_list = np.asarray(a_err_list)    
    b_err_list = np.asarray(b_err_list)    
    
    a = np.median(a_list)
    a_up = np.percentile(a_list, 84)
    a_down = np.percentile(a_list, 16)
    
    b = np.median(b_list)
    b_up = np.percentile(b_list, 84)
    b_down = np.percentile(b_list, 16)

    a_i = np.where(a_list == a)[0][0]
    b_i = np.where(b_list == b)[0][0]
    
    a_err = [[math.sqrt(a_err_list[a_i]**2 + (a - a_down)**2)], [math.sqrt(a_err_list[a_i]**2 + (a_up - a)**2)]]
    b_err = [[math.sqrt(b_err_list[b_i]**2 + (b - b_down)**2)], [math.sqrt(b_err_list[b_i]**2 + (b_up - b)**2)]]
    
    return a, a_err, b, b_err 
    
def func(x, A, B):
    return np.log10(A) + (-1*(10**x)/B)*math.log10(math.exp(1))
    
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
        
        self.ax4 = None
        self.ax5 = None
        self.ax1 = None
        self.ax2 = None
    
    def reading(self):
        source_path = r'E:/LICENSE/ProgsData/main/GAMAv3_1.txt'
        input_path = r'E:/backup/backup_BPT/GAMA_ETG_OLA.csv'
        output_path = r'E:/databases/Merged_Out.csv'
        merge_phys_databases_outflow(source_path, input_path, output_path)
        self.DataFrame = pd.read_csv(output_path)
        
        self.DataFrame_or = pd.read_csv(input_path)
        
        OUTFLOW = []
        OUTFLOW_up = []
        OUTFLOW_down = []
        OUTFLOW_off = []
        OUTFLOW_up_off = []
        OUTFLOW_down_off = []
        
        for i in range(len(self.DataFrame['SFR_0_1Gyr_percentile50'])):
                
            if self.DataFrame['LAGN_er'][i] in [-1, -2]:
                OUTFLOW.append(outflow_wAGN(self.DataFrame['SFR_0_1Gyr_percentile50'][i], self.DataFrame['LAGN'][i], self.DataFrame['mass_stellar_percentile50'][i]))
                OUTFLOW_up.append(self.DataFrame['LAGN_er'][i])
                OUTFLOW_down.append(self.DataFrame['LAGN_er'][i])
            else:
                OUTFLOW.append(outflow_wAGN(self.DataFrame['SFR_0_1Gyr_percentile50'][i], self.DataFrame['LAGN'][i], self.DataFrame['mass_stellar_percentile50'][i]))
                OUTFLOW_up.append(outflow_wAGN(self.DataFrame['SFR_0_1Gyr_percentile50'][i], self.DataFrame['LAGN'][i] + float(self.DataFrame['LAGN_er'][i]), self.DataFrame['mass_stellar_percentile50'][i]))
                OUTFLOW_down.append(outflow_wAGN(self.DataFrame['SFR_0_1Gyr_percentile50'][i], self.DataFrame['LAGN'][i] - float(self.DataFrame['LAGN_er'][i]), self.DataFrame['mass_stellar_percentile50'][i]))
            
            OUTFLOW_off.append(outflow_woAGN(self.DataFrame['SFR_0_1Gyr_percentile50'][i], self.DataFrame['mass_stellar_percentile50'][i]))
            OUTFLOW_up_off.append(outflow_woAGN(self.DataFrame['SFR_0_1Gyr_percentile84'][i], self.DataFrame['mass_stellar_percentile16'][i]))
            OUTFLOW_down_off.append(outflow_woAGN(self.DataFrame['SFR_0_1Gyr_percentile16'][i], self.DataFrame['mass_stellar_percentile84'][i]))
            
            # if str(self.DataFrame['SPECID'][i]) == '584413008453724160':
            #     print(self.DataFrame['Z_1'][i])
            #     print(self.DataFrame['SFR_0_1Gyr_percentile50'][i])
            #     print(self.DataFrame['mass_stellar_percentile50'][i])
            #     print(self.DataFrame['LAGN'][i])
            #     print(OUTFLOW_up[-1])
            #     print(OUTFLOW[-1])
            #     print(OUTFLOW_down[-1])
                
        self.DataFrame['OUT_OR'] = OUTFLOW
        self.DataFrame['OUT_OR_OFF'] = OUTFLOW_off
        self.DataFrame.to_csv(output_path, index=False)
        
        self.DataFrame['OUTFLOW'] = OUTFLOW_off
        self.DataFrame['OUTFLOW_up'] = OUTFLOW_up_off
        self.DataFrame['OUTFLOW_down'] = OUTFLOW_down_off
        
        gs_top = plt.GridSpec(2, 2, wspace=0, hspace=0)
        self.fig1 = plt.figure(figsize=(12, 12))
        # adjusting_figure_size(12, 12, 1, 0.2, 0.7, 0.2)
        adjusting_plotting_pars()
        
        ax1 = self.fig1.add_subplot(gs_top[0,0])
        ax2 = self.fig1.add_subplot(gs_top[0,1], sharey=ax1)
        ax3 = self.fig1.add_subplot(gs_top[1,0])
        ax4 = self.fig1.add_subplot(gs_top[1,1], sharey=ax3)

        ax1.tick_params(top=True, labeltop=False, bottom=True, labelbottom=False, right=True, direction='in', labelsize=20)
        ax2.tick_params(top=True, labeltop=False, bottom=True, labelbottom=False, left=True, labelleft=False, right=True, labelright=False, direction='in', labelsize=20)
        ax3.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, right=True, direction='in', labelsize=20)
        ax4.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, left=True, labelleft=False, right=True, labelright=False, direction='in', labelsize=20)

        k = 0.8
        
        ax1.text(8.7 + k*(10.0-8.7), -3 + k*(4.0 - (-3.0)), 'BPT', fontsize=20)
        ax2.text(8.7 + k*(10.0-8.7), -3 + k*(4.0 - (-3.0)), 'WHAN', fontsize=20)
        
        ax1.set_ylabel(r'$\log \mathrm{(\dot{M}_{H_2} \: / \: (M_\odot/yr))}$', fontsize=20)
        ax3.set_ylabel(r'$\log \mathrm{(\dot{M}_{H_2} \: / \: (M_\odot/yr))}$', fontsize=20)
        ax3.set_xlabel(r'$\log \mathrm{(age \: / \: yr)}$', fontsize=20)
        ax4.set_xlabel(r'$\log \mathrm{(age \: / \: yr)}$', fontsize=20)
        
        topaxes = [ax1,ax2, ax3, ax4]
        for ax in topaxes:    
            ax.set_xlim([8.7, 10.1])
            ax.set_xticks(np.arange(8.8, 10.1, 0.2))
            ax.set_ylim([-3.0, 4.0])
            ax.set_yticks(np.arange(-3.0, 3.9, 1.0))

        ax1, ax2 = Main.plotting_init(self,
            {'X' : 'ager_percentile50',
             'Y' : 'OUTFLOW',
             'Y_up' : 'OUTFLOW_up',
             'Y_down' : 'OUTFLOW_down',
             'bids_chain' : [[8.8, 9.0], [9.0, 9.2], [9.2, 9.4], [9.4, 9.6], [9.6, 9.8], [9.8, 10.0]],
             'legend' : False
             },
            ax1, ax2
        )
        
        self.DataFrame['OUTFLOW'] = OUTFLOW
        self.DataFrame['OUTFLOW_up'] = OUTFLOW_up
        self.DataFrame['OUTFLOW_down'] = OUTFLOW_down   
        
        ax3, ax4 = Main.plotting_init(self,
            {'X' : 'ager_percentile50',
             'Y' : 'OUTFLOW',
             'Y_up' : 'OUTFLOW_up',
             'Y_down' : 'OUTFLOW_down',
             'bids_chain' : [[8.8, 9.0], [9.0, 9.2], [9.2, 9.4], [9.4, 9.6], [9.6, 9.8], [9.8, 10.0]],
             'legend' : True
             },
            ax3, ax4
        )
        
        self.fig1.savefig('./FIGURES_IN_PAPER/OUTFLOW.pdf', dpi=70, transparent = True, bbox_inches = 'tight', pad_inches = 0.0001)
    
    def plotting_init(self, pars, ax1, ax2):
        Main.plotter(self, [ax1, ax2], pars['X'], pars['Y'], pars['Y_up'], pars['Y_down'], 'BPT', 'WHAN', pars['bids_chain'], pars['legend'])
        return ax1, ax2


    def plotter(self, axes, x, y, up, down, BPT_key, WHAN_key, bids_chain, legend):
        XX = []
        YY = []
        Y_ups = []
        Y_downs = []
        BPT_keys = []
        WHAN_keys = []
        ks = []
        
        X_fit = []
        Y_fit = []
        Y_fit_up = []
        Y_fit_down = []
                
        for i in range(len(self.DataFrame[x])):
            AGN = self.DataFrame[BPT_key][i]
            WHAN = self.DataFrame[WHAN_key][i]
            X = self.DataFrame[x][i]
            Y = self.DataFrame[y][i]
            Y_up = self.DataFrame[up][i]
            Y_down = self.DataFrame[down][i]
            
            alpha = 0.2
            if Y_up == -2:
                axes[0].arrow(X, Y, 0, -0.1, head_width=0.01, head_length=0.03, color=self.color_dict_BPT[AGN][0], alpha=alpha)
                axes[1].arrow(X, Y, 0, -0.1, head_width=0.01, head_length=0.03, color=self.color_dict_WHAN[WHAN][0], alpha=alpha)
                k = 0
                Y_up = 0
                Y_down = 0
            elif Y_up == -1:
                axes[0].arrow(X, Y, 0, -0.1, head_width=0.01, head_length=0.03, color=self.color_dict_BPT[AGN][0], alpha=alpha)
                axes[1].arrow(X, Y, 0, -0.1, head_width=0.01, head_length=0.03, color=self.color_dict_WHAN[WHAN][0], alpha=alpha)
                k = -1
                Y_up = 0
                Y_down = 0
            else:
                # axes[0].errorbar(X, Y, yerr = [[Y - float(Y_down)], [float(Y_up) - Y]], alpha = 0.5, color = self.color_dict_BPT[AGN][0], marker = '.')
                # axes[1].errorbar(X, Y, yerr = [[Y - float(Y_down)], [float(Y_up) - Y]], alpha = 0.5, color = self.color_dict_WHAN[WHAN][0], marker = '.')
                axes[0].scatter(X, Y, alpha = alpha, color = self.color_dict_BPT[AGN][0], marker = '.')
                axes[1].scatter(X, Y, alpha = alpha, color = self.color_dict_WHAN[WHAN][0], marker = '.')
                X_fit.append(X)
                Y_fit.append(Y)
                Y_fit_up.append(Y_up)
                Y_fit_down.append(Y_down)
                k = 1
            
            if Y > -10:
                XX.append(X)
                YY.append(Y)
                Y_ups.append(Y_up)
                Y_downs.append(Y_down)
                BPT_keys.append(AGN)
                WHAN_keys.append(WHAN)
                ks.append(k)
              
        axes[0].axhline(y = 0, color = 'black', linestyle='dashed')        
        axes[1].axhline(y = 0, color = 'black', linestyle='dashed')
               
        
        class_list_BPT = class_list_creator_w_err_out(XX, YY, Y_ups, Y_downs, BPT_keys, 'BPT', ks)
        class_list_WHAN = class_list_creator_w_err_out(XX, YY, Y_ups, Y_downs, WHAN_keys, 'WHAN', ks)

        self.means = []
        self.errs = []
        self.ages = []
        

        bids = bids_chain
        
        # ages_const = np.arange(7, 12, 1)

        # axes[0].set_xticks(ages_const)
        # axes[1].set_xticks(ages_const)
        
        classlist_plotter_uplim(axes[0], class_list_BPT, bids)
        classlist_plotter_uplim(axes[1], class_list_WHAN, bids)

        if legend == True:
            for j, item in enumerate(class_list_BPT):
                axes[0].scatter(-99, -99, alpha = 1, color=item[4][0], marker=item[4][1], s = 150, label=list_names_BPT_1[j])
            axes[0].legend(loc=3, fontsize=15)

            for j, item in enumerate(class_list_WHAN):
                axes[1].scatter(-99, -99, alpha = 1, color=item[4][0], marker=item[4][1], s = 150, label=list_names_WHAN[j])
            axes[1].legend(loc=3, fontsize=15)
        
        
        # for ax in axes:
        #     ax.plot([-3.5, 3.5], [-3.5, 3.5], c = 'k', linestyle = '--')
        
        #fancy functions part#
        
        ages = np.arange(8.8, 10, 0.05)
        
        a, SE_A, b, SE_B = fit(X_fit, Y_fit, Y_fit_up, Y_fit_down)
        
        print('a : {:3e}, +: {:3e} -: {:3e}'.format(a, SE_A[0][0], SE_A[1][0]))
        print('b : {:3e}, +: {:3e} -: {:3e}'.format(b, SE_B[0][0], SE_B[1][0]))
        
        axes[0].plot(ages, np.log10(a) + (-1*(10**ages)/b)*math.log10(math.exp(1)), c='k') 
        axes[1].plot(ages, np.log10(a) + (-1*(10**ages)/b)*math.log10(math.exp(1)), c='k') 


if __name__ == '__main__':
    obj = Main('GAMA_ETG_OLA.csv')
    obj.reading()
    # obj.plotting()