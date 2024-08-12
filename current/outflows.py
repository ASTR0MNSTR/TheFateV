import csv
import matplotlib.pyplot as plt 
import numpy as np
import math
# import statistics as st
import pandas as pd
from scipy.optimize import curve_fit
from __legpars__ import *
from __stats__ import *
from __plt__2 import *     
    
def plotting_init(pars, ax1, ax2):
    plotter(pars['db'], [ax1, ax2], pars['X'], pars['Y'], pars['Y_up'], pars['Y_down'], 'BPT', 'WHAN', pars['bids_chain'], pars['legend'])
    return ax1, ax2

def plotter(DataFrame, axes, x, y, up, down, BPT_key, WHAN_key, bids_chain, legend):
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
                
        for i in range(len(DataFrame[x])):
            AGN = DataFrame[BPT_key][i]
            WHAN = DataFrame[WHAN_key][i]
            X = DataFrame[x][i]
            Y = DataFrame[y][i]
            Y_up = DataFrame[up][i]
            Y_down = DataFrame[down][i]
            
            alpha = 0.2
            if Y_up == -2:
                # axes[0].arrow(X, Y, 0, -0.1, head_width=0.01, head_length=0.03, color=color_dict_BPT[AGN][0], alpha=alpha)
                # axes[1].arrow(X, Y, 0, -0.1, head_width=0.01, head_length=0.03, color=cd_WHAN[WHAN][0], alpha=alpha)
                k = 0
                Y_up = 0
                Y_down = 0
            elif Y_up == -1:
                # axes[0].arrow(X, Y, 0, -0.1, head_width=0.01, head_length=0.03, color=color_dict_BPT[AGN][0], alpha=alpha)
                # axes[1].arrow(X, Y, 0, -0.1, head_width=0.01, head_length=0.03, color=cd_WHAN[WHAN][0], alpha=alpha)
                k = -1
                Y_up = 0
                Y_down = 0
            else:
                # axes[0].errorbar(X, Y, yerr = [[Y - float(Y_down)], [float(Y_up) - Y]], alpha = 0.5, color = color_dict_BPT[AGN][0], marker = '.')
                # axes[1].errorbar(X, Y, yerr = [[Y - float(Y_down)], [float(Y_up) - Y]], alpha = 0.5, color = cd_WHAN[WHAN][0], marker = '.')
                # axes[0].scatter(X, Y, alpha = alpha, color = color_dict_BPT[AGN][0], marker = '.', s = color_dict_BPT[AGN][1])
                # axes[1].scatter(X, Y, alpha = alpha, color = cd_WHAN[WHAN][0], marker = '.', s = cd_WHAN[WHAN][1])
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

        means = []
        errs = []
        ages = []
        
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
        