import numpy as np
import matplotlib.pyplot as plt
from __legpars__ import *
from __stats__ import *
import pandas as pd

def plotting(pars_dict):
    cols = [pars_dict['x'], pars_dict['y'], pars_dict['up'], pars_dict['down'], 'BPT', 'WHAN']
    DataFrame = pd.read_csv(pars_dict['input_path'], usecols=cols)
    gs_top = plt.GridSpec(1, 2, wspace=0)
    fig1 = plt.figure(figsize=(12, 6), tight_layout=True)

    ax4 = fig1.add_subplot(gs_top[:,0])
    ax5 = fig1.add_subplot(gs_top[:,1], sharey=ax4)

    topaxes = [ax5, ax4]

    ax4.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, right=True, direction='in')
    ax5.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, left=True, labelleft=False, right=True, labelright=False, direction='in')

    ax4.set_ylabel(pars_dict['ylabel'])
    for ax in topaxes:    
        ax.set_xlim(pars_dict['xlim'])
        ax.set_ylim(pars_dict['ylim'])
        ax.set_yticks(pars_dict['yticks'])
        ax.set_xticks(pars_dict['xticks'])
        ax.set_xlabel(pars_dict['xlabel'])
            
    bids = pars_dict['bids']
        
    ax4 = phys_plotter(ax4, DataFrame[pars_dict['x']], DataFrame[pars_dict['y']], DataFrame[pars_dict['up']], DataFrame[pars_dict['down']], DataFrame['BPT'], bids, 'BPT', True)
    ax5 = phys_plotter(ax5, DataFrame[pars_dict['x']], DataFrame[pars_dict['y']], DataFrame[pars_dict['up']], DataFrame[pars_dict['down']], DataFrame['WHAN'], bids, 'WHAN', True)

    #fig1.savefig(pars_dict['save_path'])

def phys_plotter(axis, x, y, up, down, AGN_keys, bids, WHAN_or_BPT, leg):    

    if WHAN_or_BPT == 'BPT':
        color_dict = color_dict_BPT
    else:
        color_dict = cd_WHAN
    for i, item in enumerate(y):
        X = x[i]
        Y = item
        AGN = AGN_keys[i]
        axis.scatter(X, Y, alpha=0.4, color = color_dict[AGN][0], marker='.', s = color_dict[AGN][1])
                    
    class_list = class_list_creator(x, y, up, down, AGN_keys, WHAN_or_BPT)
    
    errs = []
    means = []
    
    for item in class_list:
        X_plot = []
        Y_plot = []
        X, Y, err, length = bootstrapper(item[0], item[1], item[2], item[3], bids)
        means.append(Y)
        errs.append(err)
        for j in range(len(X)):
            if Y[j] != -99:
                axis.errorbar(X[j], Y[j], alpha = 1, xerr=0, yerr= err[j], color=item[4][0], fmt=item[4][1], ms = 12)
                axis.text(X[j], Y[j], length[j], c = 'red')
                X_plot.append(X[j])
                Y_plot.append(Y[j])
        axis.plot(X_plot, Y_plot, alpha = 1, color=item[4][0])
    
    if leg == True:
        for j, item in enumerate(class_list):
            if WHAN_or_BPT == 'BPT':
                list_names = list_names_BPT
            elif WHAN_or_BPT == 'WHAN':
                list_names = list_names_WHAN
            axis.scatter(-99, -99, alpha = 1, color=item[4][0], marker=item[4][1], s = 150, label=list_names[j])
        axis.legend(fontsize="13")
    
    return axis
        
        

def class_list_creator(x, y, up, down, AGN_keys, WHAN_or_BPT):
    if WHAN_or_BPT == 'BPT':
        keys = [['AGNXY', 'AGNX', 'AGNY'], ['UNCXY', 'UNCX', 'UNCY'], ['SFXY', 'SFX', 'SFY'], ['NOEL']]
        colors_markers = [['midnightblue', 'P'], ['springgreen', 'H'], ['mediumvioletred', '*'], ['orchid', 'o']]
    elif WHAN_or_BPT == 'WHAN':
        keys = [['sAGN'], ['wAGN'], ['SF'], ['ELR'], ['RG'], ['LLR']]
        colors_markers = [['midnightblue', 'P'], ['blue', 'P'], ['mediumvioletred', '*'], ['sandybrown', 'D'], ['chocolate', 'o'], ['maroon', 'o']]
    class_list = []
    
    for j, chain in enumerate(keys):
        X = []
        Y = []
        Y_up = []
        Y_down = []
        for i, item in enumerate(y):
            AGN = AGN_keys[i]
            if AGN in chain:
                X.append(x[i])
                Y.append(y[i])
                Y_up.append(up[i])
                Y_down.append(down[i])
        
        class_list.append([X, Y, Y_down, Y_up, colors_markers[j]])
    return class_list
        # class_list = [[yes_temp_age, yes_temp, 'midnightblue', yes_temp_down, yes_temp_up, 'P'], [UNC_temp_age, UNC_temp, 'springgreen', UNC_temp_down, UNC_temp_up, 'H'], [no_temp_age, no_temp, 'mediumvioletred', no_temp_down, no_temp_up, '*'], [noel_temp_age, noel_temp, 'orchid', noel_temp_down, noel_temp_up, 'o']]
    