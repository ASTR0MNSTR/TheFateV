import numpy as np
import matplotlib.pyplot as plt
from __legpars__ import *
from __stats__ import *
import pandas as pd
from scipy import ndimage
import matplotlib as mpl

#ROUND HISTOS

def my_level_list(data, kwarg):
    list = []
    if kwarg == 'WHAN':
        labels = ['sAGN', 'wAGN', 'UNC', 'SFG', 'ELR', 'LLR', 'NER']
    elif kwarg == 'BPT':
        labels = ['AGNXY', 'AGNX', 'UNCXY', 'UNCX', 'UNCY', 'SFGXY', 'SFGX', 'SFGY', 'NOEL']

    for i in range(len(data)):
        if (data[i]*100/np.sum(data)) > 3: #2%
            list.append(labels[i])
        else:
            list.append('')
    return list

def my_autopct_BPT(pct):
    return (f'{pct:.2f}%') if pct > 3 else ''
    
def my_autopct_WHAN(pct):
    return (f'{pct:.2f}%') if pct > 3 and round(pct, 2) != 4.96 else ''

def merging_BPT(list_obj):
    AGNXY, AGNX, UNCXY, UNCX, UNCY, SFXY, SFX, SFY, NOEL = list_obj
        
    AGNs = [AGNXY, AGNX]
    UNCs = [UNCXY, UNCX, UNCY]
    SFs = [SFXY, SFX, SFY]
    NOELs = [NOEL]
        
    arrays = [AGNs, UNCs, SFs, NOELs]
    indexes = []
    results = []
    for i, array in enumerate(arrays):
        if len(array) == 1 or array.count(0) == len(array) - 1:
            indexes.append(i)
        results.append(sum(array))

    return results
    
def merging_WHAN(list_obj):

    sAGN, wAGN, UNC, SF, ELR, LLR, RG = list_obj

    AGNs = [sAGN, wAGN]
    UNCs = [UNC]
    SFs = [SF]
    RGs = [ELR, LLR, RG]
        
    arrays = [AGNs, UNCs, SFs, RGs]
    indexes = []
    results = []
    for i, array in enumerate(arrays):
        if len(array) == 1 or array.count(0) == len(array) - 1:
            indexes.append(i)
        results.append(sum(array))

    return results

def short_WHAN_in(data):
    def my_format(pct):
        total = sum(data)
        val = int(round(pct*total/100.0))
        if data[1] == int(val) or data[2] == int(val) or pct == 100 or pct < 3:
            return ''
        else:
            return (f'{pct:.2f}%')
    return my_format
    
def short_BPT_in(data):
    def my_format(pct):
        total = sum(data)
        val = int(round(pct*total/100.0))
        if data[3] == int(val) or pct == 100 or pct < 3:
            return ''
        else:
            return (f'{pct:.2f}%')
    return my_format

#PHYSICS

def adjusting_plotting_pars():
    plt.rcParams['font.size'] = 20
    mpl.rcParams['xtick.labelsize'] = 20
    mpl.rcParams['ytick.labelsize'] = 20 
    mpl.rcParams['axes.labelsize'] = 20
    mpl.rcParams['ytick.major.size'] = 15
    mpl.rcParams['xtick.major.size'] = 15
    # print(mpl.rcParams.keys())
    
def adjusting_figure_size(figw, figh, l, r, b, t):
    # plt.subplots_adjust(left=1/figw, right=1-0.2/figw, bottom=0.7/figh, top=1-0.2/figh)
    plt.subplots_adjust(left=l/figw, right=1-r/figw, bottom=b/figh, top=1-t/figh)
    
def generating_annotation(axis, x, y, text):
    axis.text(x, y, text)

def theor_lines(axes, key):
    
    for ax in axes:
        if key == 'mdms':
            # x = np.arange(7.9, 10.1, 0.01)
            # a_all = 2.3152363967664957
            # b_all = 5.276003555081468
            # y = ((-(10**x)/10**9) / b_all) - a_all
            # ax.plot(x, y, color='k', linestyle='solid')
            
            # x = np.arange(7.9, 10.1, 0.01)
            # a_all = 2.3152363967664957
            # b_all = 5.276003555081468
            # y = ((-(10**x)/10**9) / b_all) - a_all
            # ax.plot(x, y, color='k', linestyle='solid')
            
            x = np.arange(7.9, 10.1, 0.01)
            a_all = 2.31
            b_all = 2.26 
            y = ((-(10**x)/10**9) / b_all)*np.log10(2.7182) - a_all
            ax.plot(x, y, color='k', linestyle='solid')
            
            # x = np.arange(7.9, 10.1, 0.01)
            # a_all = 2.31
            # b_all = 2.26 
            # y = ((-(10**x)/10**9) / b_all)*np.log10(2.7182) - a_all
            # ax.plot(x, y, color='k', linestyle='dashed')
            
        elif key == 'sfrsm':
            x = np.arange(6.9, 12, 0.1)
            # ax.plot(x, (0.84 - 0.026*11.636)*x - (6.51 - 0.11*11.636), color='k', linestyle='solid')
            # ax.text(11.3, (0.84 - 0.026*11.636)*11.3 - (6.5 - 0.11*11.636), 'z = 0.17', rotation=9)
            
            ax.plot(x, (0.84 - 0.026*10.594)*x - (6.51 - 0.11*10.594), color='k', linestyle='solid')
            ax.text(11.3, (0.84 - 0.026*10.594)*11.3 - (6.5 - 0.11*10.594), 'z = 0.27', rotation=9)
            
            # ax.plot(x, (0.84 - 0.026*11.4074)*x - (6.51 - 0.11*11.4074), linestyle='dashed', color='k')
            # ax.text(11.3, (0.84 - 0.026*11.4074)*11.3 - (6.5 - 0.11*11.4074), 'z = 0.19', rotation=9)
            # ax.plot(x, (0.84 - 0.026*9.8615801)*x - (6.51 - 0.11*9.8615801), linestyle='solid', color='k')
            # ax.text(11.3, (0.84 - 0.026*9.8615801)*11.3 - (6.5 - 0.11*9.8615801), 'z = 0.33', rotation=9)
            
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
        
    a_array = np.array(a_list)
    b_array = np.array(b_list)
    b_array = b_array/(10**9)
    a = np.percentile(a_array, 50)
    b = np.percentile(b_array, 50)
    print(a, '+', np.percentile(a_array, 84) - a, '-', a - np.percentile(a_array, 16))
    print(b, '+', np.percentile(b_array, 84) - b, '-', b - np.percentile(b_array, 16))
    print('\n')
    return a, b
    
            
def plotting(pars_dict):
    if 'err' in pars_dict.keys():
        cols = [pars_dict['x'], pars_dict['y'], pars_dict['err'], 'BPT', 'WHAN']
    else:
        cols = [pars_dict['x'], pars_dict['y'], pars_dict['up'], pars_dict['down'], 'BPT', 'WHAN']
    DataFrame = pd.read_csv(pars_dict['input_path'], usecols=cols)
    gs_top = plt.GridSpec(1, 2, wspace=0)
    fig1 = plt.figure(figsize=(12, 6), tight_layout=True)
    adjusting_plotting_pars()

    ax4 = fig1.add_subplot(gs_top[:,0])
    ax5 = fig1.add_subplot(gs_top[:,1], sharey=ax4)
    
    topaxes = [ax5, ax4]

    ax4.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, right=True, direction='in', labelsize=17)
    ax5.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, left=True, labelleft=False, right=True, labelright=False, direction='in', labelsize=17)
    
    ax4.set_ylabel(pars_dict['ylabel'], fontsize=17)
    for ax in topaxes:    
        ax.set_xlim(pars_dict['xlim'])
        ax.set_ylim(pars_dict['ylim'])
        ax.set_yticks(pars_dict['yticks'])
        ax.set_xticks(pars_dict['xticks'])
        ax.set_xlabel(pars_dict['xlabel'], fontsize=17)
    
    k = 0.8
    generating_annotation(ax4, pars_dict['xlim'][0] + k*(pars_dict['xlim'][1] - pars_dict['xlim'][0]), pars_dict['ylim'][0] + k*(pars_dict['ylim'][1] - pars_dict['ylim'][0]), 'BPT')    
    generating_annotation(ax5, pars_dict['xlim'][0] + k*(pars_dict['xlim'][1] - pars_dict['xlim'][0]), pars_dict['ylim'][0] + k*(pars_dict['ylim'][1] - pars_dict['ylim'][0]), 'WHAN')    
    # generating_annotation(ax5, 11.0, 1.8, '0.22 < z < 0.33')    
    bids = pars_dict['bids']
    
    try:
        theor_lines([ax4, ax5], pars_dict['theor_lines'])
    except:
        pass
    
    if 'err' in pars_dict.keys():
        # ax4 = phys_plotter(ax4, DataFrame[pars_dict['x']], DataFrame[pars_dict['y']], DataFrame[pars_dict['y']] + DataFrame[pars_dict['err']], DataFrame[pars_dict['y']] - DataFrame[pars_dict['err']], DataFrame['BPT'], bids, 'BPT', True)
        ax4 = phys_plotter(ax4, DataFrame[pars_dict['x']], DataFrame[pars_dict['y']], DataFrame[pars_dict['y']] + DataFrame[pars_dict['err']], DataFrame[pars_dict['y']] - DataFrame[pars_dict['err']], DataFrame['BPT'], bids, 'BPT', True)
        # ax5 = phys_plotter(ax5, DataFrame[pars_dict['x']], DataFrame[pars_dict['y']], DataFrame[pars_dict['y']] + DataFrame[pars_dict['err']], DataFrame[pars_dict['y']] - DataFrame[pars_dict['err']], DataFrame['WHAN'], bids, 'WHAN', True)
        ax5 = phys_plotter(ax5, DataFrame[pars_dict['x']], DataFrame[pars_dict['y']], DataFrame[pars_dict['y']] + DataFrame[pars_dict['err']], DataFrame[pars_dict['y']] - DataFrame[pars_dict['err']], DataFrame['WHAN'], bids, 'WHAN', True)
    else:
        ax4 = phys_plotter(ax4, DataFrame[pars_dict['x']], DataFrame[pars_dict['y']], DataFrame[pars_dict['up']], DataFrame[pars_dict['down']], DataFrame['BPT'], bids, 'BPT', True)
        ax5 = phys_plotter(ax5, DataFrame[pars_dict['x']], DataFrame[pars_dict['y']], DataFrame[pars_dict['up']], DataFrame[pars_dict['down']], DataFrame['WHAN'], bids, 'WHAN', True)

    fig1.savefig(pars_dict['save_path'], dpi=300, transparent = True, bbox_inches = 'tight', pad_inches = 0.0001)

def contours(x , y , lev, sigma, bin1, bin2):
        #mode{‘reflect’, ‘constant’, ‘nearest’, ‘mirror’, ‘wrap’}, optional\n”,
        x_g1d = ndimage.gaussian_filter1d(x, sigma, mode = 'mirror')
        y_g1d = ndimage.gaussian_filter1d(y, sigma, mode = 'mirror')
        H, xedges, yedges = np.histogram2d(y_g1d, x_g1d, bins=(bin1,bin2))
        extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
        levels = []
        for i in range(len(lev)):
            levels.append(np.max(H) * lev[i])
        return (H, extent, levels)
    
def contour_plotter(axis, classlist):
    for item in classlist:
        X = np.array(item[0])
        Y = np.array(item[1])
        line_widths = (2.9,2,1.)
        levels = [0.05, 0.5] # theses levels are ‘lev’ in the above def, so they define the percentage of the contours: eg.: 0.05 is the contour at 5% of the hight (enclosing 95% of the data), and 0.9 is the contours at 90% of the hight (enclosing 10% of the data)
        H = contours(X, Y , levels, sigma = 2, bin1 = 20, bin2 = 20)
        axis.contour(H[0], levels = H[2], origin='lower', colors=item[4][0], linewidths=line_widths, extent=H[1], alpha = 1)

def rainbow_plotter(axis, classlist):
    for item in classlist:
        popt, dev, ssfr50, ssfr_down, ssfr_up = width_estimation(item[0], item[1])
        X_plot = np.linspace(10.0, 11.5, 1000)
        axis.fill_between(X_plot, linear_function(X_plot, popt[0], popt[1] + dev), linear_function(X_plot, popt[0], popt[1] - dev), color = item[4][0], alpha = 0.17)
        axis.plot(X_plot, linear_function(X_plot, popt[0], popt[1] + dev), alpha = 1, color=item[4][0])
        axis.plot(X_plot, linear_function(X_plot, popt[0], popt[1] - dev), alpha = 1, color=item[4][0])
        axis.plot(X_plot, linear_function(X_plot, popt[0], popt[1]), alpha = 1, color=item[4][0], linestyle = '--')
        print(f'{item[-1]} & ${ssfr50}\pm^{ssfr_up}_{ssfr_down}$ & {popt[0]} & {popt[1]} & {dev} {chr(92)}{chr(92)}')
        
    
def classlist_plotter(axis, classlist, bids):
    errs = []
    means = []
    for item in classlist:
        X_plot = []
        Y_plot = []
        err_plot = []
        err_up = []
        err_down = []
        X, Y, err, length, res = monte_carlo(item[0], item[1], item[2], item[3], bids)
        # X, Y, err, length = monte_carlo(item[0], item[1], item[2], item[3], bids)
        means.append(Y)
        errs.append(err)
        for j in range(len(X)):
            if Y[j] != -99:
                    # axis.errorbar(X[j], Y[j], alpha = 1, xerr=0, yerr= err[j], color=item[4][0], fmt=item[4][1], ms = 12)
                    # axis.scatter(X[j], Y[j], alpha = 1, color=item[4][0], marker=item[4][1], s = 100)
                # axis.text(X[j], Y[j], length[j], c = 'red')
                X_plot.append(X[j])
                Y_plot.append(Y[j])
                err_plot.append(err[j])
                    
        X_plot = np.asarray(X_plot)
        Y_plot = np.asarray(Y_plot)
        for err in err_plot:
            err_up.append(err[1][0])
            err_down.append(err[0][0])
        err_up = np.asarray(err_up)
        err_down = np.asarray(err_down)
        
        axis.fill_between(X_plot, Y_plot + err_up, Y_plot - err_down, color = item[4][0], alpha = 0.17)
        axis.scatter(X_plot, Y_plot, alpha = 1, color=item[4][0], marker=item[4][1], s = 100, edgecolors='black')
        axis.plot(X_plot, Y_plot + err_up, alpha = 1, color=item[4][0])
        axis.plot(X_plot, Y_plot - err_down, alpha = 1, color=item[4][0])
        axis.plot(X_plot, Y_plot, alpha = 1, color=item[4][0], linestyle = '--')
        axis.text(X_plot[-1], Y_plot[-1], round(res.pvalue, 5), color=item[4][0])

def classlist_plotter_uplim(axis, classlist, bids):
    errs = []
    means = []
    pert = np.linspace(-0.05, 0.05, num=len(classlist))
    for i, item in enumerate(classlist):
        X_plot = []
        Y_plot = []
        err_plot = []
        err_up = []
        err_down = []
        up_lim_end = []
            
        X, Y, err, length, res = monte_carlo(item[0], item[1], item[2], item[3], bids)
        # X, Y, err, length = monte_carlo(item[0], item[1], item[2], item[3], bids)
        up_lim = up_lim_analysis(item[0], item[5], bids)
        means.append(Y)
        errs.append(err)
        for j in range(len(X)):
            if Y[j] != -99:
                    # axis.errorbar(X[j], Y[j], alpha = 1, xerr=0, yerr= err[j], color=item[4][0], fmt=item[4][1], ms = 12)
                    # axis.scatter(X[j], Y[j], alpha = 1, color=item[4][0], marker=item[4][1], s = 100)
                    # self.ax4.text(X[j], Y[j], length[j], c = 'red')
                X_plot.append(X[j] + pert[i])
                Y_plot.append(Y[j])
                err_plot.append(err[j])
                up_lim_end.append(up_lim[j])
                    
        X_plot = np.asarray(X_plot)
        Y_plot = np.asarray(Y_plot)
        for err in err_plot:
            err_up.append(err[1][0])
            err_down.append(err[0][0])
        err_up = np.asarray(err_up)
        err_down = np.asarray(err_down)
            
        axis.fill_between(X_plot, Y_plot + err_up, Y_plot - err_down, color = item[4][0], alpha = 0.17)
        axis.scatter(X_plot, Y_plot, alpha = 1, color=item[4][0], marker=item[4][1], s = 100, edgecolors='black')
        for i, elem in enumerate(up_lim_end):
            if elem:
                axis.arrow(X_plot[i], Y_plot[i], 0, -0.3, width = 0.007, alpha = 1, color=item[4][0])
        axis.plot(X_plot, Y_plot + err_up, alpha = 1, color=item[4][0])
        axis.plot(X_plot, Y_plot - err_down, alpha = 1, color=item[4][0])
        axis.plot(X_plot, Y_plot, alpha = 1, color=item[4][0], linestyle = '--')
        axis.text(X_plot[-1], Y_plot[-1], round(res.pvalue, 5), color=item[4][0])
        print(item[4][0], ' : ', round(res.pvalue, 5))

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
                    
    class_list = class_list_creator_w_err(x, y, up, down, AGN_keys, WHAN_or_BPT)
    
    # a, b = fit(x, y, up, down)
    
    # a = -2.257
    # b = 2.213
    # x_theor = np.arange(8.6, 10.2, 0.01)
    # axis.plot(x_theor, a - 0.4343*(10**x_theor)/(b*(10**9)), color='orange', linestyle='dashed', linewidth=3.)
    
    classlist_plotter(axis, class_list, bids)
    # rainbow_plotter(axis, class_list)
    # contour_plotter(axis, class_list)
    
    if leg == True:
        for j, item in enumerate(class_list):
            if WHAN_or_BPT == 'BPT':
                list_names = list_names_BPT_1
            elif WHAN_or_BPT == 'WHAN':
                list_names = list_names_WHAN
            axis.scatter(-99, -99, alpha = 1, color=item[4][0], marker=item[4][1], s = 150, label=list_names[j], edgecolors='black')
        axis.legend(loc=3, fontsize=15)
    
    return axis   

def class_list_creator_w_err(x, y, up, down, AGN_keys, WHAN_or_BPT):
    if WHAN_or_BPT == 'BPT':
        keys = [['AGNXY'], ['AGNX'], ['UNCXY'], ['UNCX'], ['SFGXY'], ['SFGX'], ['NOEL']]
        colors_markers = [['midnightblue', 'P'], ['dodgerblue', 'P'], ['springgreen', 'H'], ['darkgreen', 'H'], ['mediumvioletred', '*'], ['crimson', 'p'], ['silver', 'o']]
    elif WHAN_or_BPT == 'WHAN':
        # keys = [['sAGN'], ['wAGN'], ['SFG'], ['ELR'], ['NER'], ['LLR']]
        keys = [['sAGN'], ['wAGN'], ['SFG'], ['ELR'], ['NER'], ['LLR']]
        # colors_markers = [['midnightblue', 'P'], ['blue', 'P'], ['mediumvioletred', '*'], ['sandybrown', 'D'], ['chocolate', 'o'], ['maroon', 'o']]
        colors_markers = [['midnightblue', 'P'], ['blue', 'o'], ['mediumvioletred', '*'], ['sandybrown', 'D'], ['chocolate', '^'], ['maroon', 'v']]
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
        
        class_list.append([X, Y, Y_up, Y_down, colors_markers[j]])
    return class_list
        # class_list = [[yes_temp_age, yes_temp, 'midnightblue', yes_temp_down, yes_temp_up, 'P'], [UNC_temp_age, UNC_temp, 'springgreen', UNC_temp_down, UNC_temp_up, 'H'], [no_temp_age, no_temp, 'mediumvioletred', no_temp_down, no_temp_up, '*'], [noel_temp_age, noel_temp, 'orchid', noel_temp_down, noel_temp_up, 'o']]

def class_list_creator_wo_err(x, y, age, AGN_keys, WHAN_or_BPT):
    if WHAN_or_BPT == 'BPT':
        keys = [['AGNXY'], ['AGNX'], ['UNCXY'], ['UNCX'], ['SFGXY'], ['SFGX'], ['NOEL']]
        colors_markers = [['midnightblue', 'P'], ['dodgerblue', 'P'], ['springgreen', 'H'], ['darkgreen', 'H'], ['mediumvioletred', '*'], ['crimson', 'p'], ['silver', 'o']]
    elif WHAN_or_BPT == 'WHAN':
        keys = [['sAGN'], ['wAGN'], ['SFG'], ['ELR'], ['NER'], ['LLR']]
        colors_markers = [['midnightblue', 'P'], ['blue', 'o'], ['mediumvioletred', '*'], ['sandybrown', 'D'], ['chocolate', '^'], ['maroon', 'v']]
    class_list = []
    
    for j, chain in enumerate(keys):
        X = []
        Y = []
        AGE = []
        for i, item in enumerate(y):
            AGN = AGN_keys[i]
            if AGN in chain:
                X.append(x[i])
                Y.append(y[i])
                AGE.append(age[i])
                
        class_list.append([X, Y, AGE, colors_markers[j]])
    return class_list

def class_list_creator_w_err_out(x, y, up, down, AGN_keys, WHAN_or_BPT, ks):
    if WHAN_or_BPT == 'BPT':
        keys = [['AGNXY'], ['AGNX'], ['UNCXY'], ['UNCX'], ['SFGXY'], ['SFGX'], ['NOEL'], ['AGNXY', 'AGNX', 'UNCXY', 'UNCX', 'SFGXY', 'SFGX', 'NOEL']]
        colors_markers = [['midnightblue', 'P'], ['dodgerblue', 'P'], ['springgreen', 'H'], ['darkgreen', 'H'], ['mediumvioletred', '*'], ['crimson', 'p'], ['silver', 'o'], ['black', 'h']]
    elif WHAN_or_BPT == 'WHAN':
        # keys = [['sAGN'], ['wAGN'], ['SFG'], ['ELR'], ['NER'], ['LLR']]
        keys = [['sAGN'], ['wAGN'], ['SFG'], ['ELR'], ['NER'], ['LLR'], ['sAGN', 'wAGN', 'SFG', 'ELR', 'NER', 'LLR']]
        # colors_markers = [['midnightblue', 'P'], ['blue', 'P'], ['mediumvioletred', '*'], ['sandybrown', 'D'], ['chocolate', 'o'], ['maroon', 'o']]
        colors_markers = [['midnightblue', 'P'], ['blue', 'o'], ['mediumvioletred', '*'], ['sandybrown', 'D'], ['chocolate', '^'], ['maroon', 'v'], ['black', 'h']]
    class_list = []
    
    for j, chain in enumerate(keys):
        X = []
        Y = []
        Y_up = []
        Y_down = []
        KS = []
        for i, item in enumerate(y):
            AGN = AGN_keys[i]
            if AGN in chain:
                X.append(x[i])
                Y.append(y[i])
                Y_up.append(up[i])
                Y_down.append(down[i])
                KS.append(ks[i])
        
        class_list.append([X, Y, Y_up, Y_down, colors_markers[j], KS])
    return class_list

def bin_stats(pars_dict):
    
    DataFrame = pd.read_csv(pars_dict['input_path'])
    
    adjusting_plotting_pars()
    
    gs_top = plt.GridSpec(6, 1, hspace=0, wspace=0)
    fig = plt.figure(figsize=(12, 18), tight_layout=True)

    ax1 = fig.add_subplot(gs_top[0,0])
    ax2 = fig.add_subplot(gs_top[1,0])
    ax3 = fig.add_subplot(gs_top[2,0])
    ax4 = fig.add_subplot(gs_top[3,0])
    ax5 = fig.add_subplot(gs_top[4,0])
    ax6 = fig.add_subplot(gs_top[5,0])
    
    plotter_histo_BPT(ax3, [0], DataFrame[pars_dict['x']], DataFrame['BPT'], 'BMS, %', DataFrame['BMS'], pars_dict['bins'])

    plotter_histo_BPT(ax2, [1], DataFrame[pars_dict['x']], DataFrame['BPT'], 'MS, %', DataFrame['BMS'], pars_dict['bins'])

    plotter_histo_BPT(ax1, [0, 1], DataFrame[pars_dict['x']], DataFrame['BPT'], 'Total, %', DataFrame['BMS'], pars_dict['bins'])

    plotter_histo_WHAN(ax6, [0], DataFrame[pars_dict['x']], DataFrame['WHAN'], 'BMS, %', DataFrame['BMS'], pars_dict['bins'])

    plotter_histo_WHAN(ax5, [1], DataFrame[pars_dict['x']], DataFrame['WHAN'], 'MS, %', DataFrame['BMS'], pars_dict['bins'])

    plotter_histo_WHAN(ax4, [0,1], DataFrame[pars_dict['x']], DataFrame['WHAN'], 'Total, %', DataFrame['BMS'], pars_dict['bins'])

    ax2.legend(fontsize="10", loc='center right')
    ax5.legend(fontsize="13", loc='center right')

    top_axes = [ax2, ax3, ax5]
    for item in top_axes:
        item.tick_params(top=False, labeltop=False, bottom=True, labelbottom=False, right=True, direction='in')

    ax1.tick_params(top=False, labeltop=False, bottom=True, labelbottom=False, right=True, direction='in')
    ax1.text(5.5, 80, 'BPT', fontsize='15', ha='center', va='center')
    ax4.text(5.5, 80, 'WHAN', fontsize='15', ha='center', va='center')
    #ax1.xaxis.set_label_position('top')
    # ax1.set_title('BPT classification')
    ax4.tick_params(top=False, labeltop=False, bottom=True, labelbottom=False, right=True, direction='in')
    #ax4.xaxis.set_label_position('top')
    # ax4.set_title('WHAN classification')

    ax6.set_xlabel(pars_dict['xlabel'], fontsize=17)
    # ax6.set_xticks([r for r in range(6)], ['<0.05', '0.05-0.10', '0.10-0.15', '0.15-0.20', '0.20-0.25', '0.25-0.33'])
    ax6.set_xticks([r for r in range(len(pars_dict['bins_names']))], pars_dict['bins_names'], fontsize=17)

    plt.savefig(pars_dict['save_path'])

def empty():
    return [0, 0, 0, 0, 0, 0, 0]

def plotter_histo_BPT(axes, BMS_condition, age, SC_BPT, y_name, BMS, bins):
        
        UNC = empty()
        UNCX = empty()
        UNCY = empty()
        SFG = empty()
        SFX = empty()
        SFY = empty()
        AGN = empty()
        AGNX = empty()
        AGNY = empty()
        NOEL = empty()
        NDA = empty()
        SAMPLE = empty()
        UNC_perc = empty()
        UNCX_perc = empty()
        UNCY_perc = empty()
        SF_perc = empty()
        SFX_perc = empty()
        SFY_perc = empty()
        AGN_perc = empty()
        AGNX_perc = empty()
        AGNY_perc = empty()
        NOEL_perc = empty()
        NDA_perc = empty()

        for i in range(len(bins)):
            for j in range(len(age)):
                if age[j] > bins[i][0] and age[j] <= bins[i][1] and BMS[j] in BMS_condition: #below-MS=0 / MS galaxies=1 and item['BMS'] == 1:
                    SAMPLE[i] += 1
                    if SC_BPT[j] == 'SFGXY':
                        SFG[i] += 1
                    elif SC_BPT[j] == 'SFGX':
                        SFX[i] += 1
                    elif SC_BPT[j] == 'SFGY':
                        SFY[i] += 1
                    elif SC_BPT[j] == 'NOEL':
                        NOEL[i] += 1
                    elif SC_BPT[j] == 'AGNXY':
                        AGN[i] += 1
                    elif SC_BPT[j] == 'AGNX':
                        AGNX[i] += 1
                    elif SC_BPT[j] == 'AGNY':
                        AGNY[i] += 1
                    elif SC_BPT[j] == 'UNCXY':
                        UNC[i] += 1
                    elif SC_BPT[j] == 'UNCX':
                        UNCX[i] += 1
                    elif SC_BPT[j] == 'UNCY':
                        UNCY[i] += 1
                    elif SC_BPT[j] == 'NDA':
                        #NDA[i] += 1
                        SAMPLE[i] -= 1
        
        for j in range(len(SAMPLE)):
            try:
                SFY_perc[j] = (SFY[j]/SAMPLE[j])*100
                SFX_perc[j] = (SFY[j] + SFX[j])*100/SAMPLE[j]
                SF_perc[j] = (SFY[j] + SFX[j]+SFG[j])*100/SAMPLE[j]
                NDA_perc[j] = (NDA[j] + SFY[j] + SFX[j]+SFG[j])*100/SAMPLE[j]
                NOEL_perc[j] = (NDA[j] + SFY[j] + SFX[j]+SFG[j]+NOEL[j])*100/SAMPLE[j]
                UNCY_perc[j] = (NDA[j] + SFY[j] + SFX[j]+SFG[j]+NOEL[j] + UNCY[j])*100/SAMPLE[j]
                UNCX_perc[j] = UNCY_perc[j] + (UNCX[j]/SAMPLE[j])*100
                UNC_perc[j] = UNCX_perc[j] + (UNC[j]/SAMPLE[j])*100
                
                AGNY_perc[j] = UNC_perc[j] + (AGNY[j]/SAMPLE[j])*100
                AGNX_perc[j] = AGNY_perc[j] + (AGNX[j]/SAMPLE[j])*100
                AGN_perc[j] = AGNX_perc[j] + (AGN[j]/SAMPLE[j])*100
                
            except:
                pass

        print('AGN', AGN, sum(AGN)) 
        print('AGNX', AGNX, sum(AGNX))
        print('AGNY', AGNY, sum(AGNY))
        print('UNC', UNC, sum(UNC))
        print('UNCX', UNCX, sum(UNCX))
        print('UNCY', UNCY, sum(UNCY))
        print('SFGXY', SFG, sum(SFG))
        print('SFGX', SFX, sum(SFX))
        print('SFGY', SFY, sum(SFY))
        print('NOEL', NOEL, sum(NOEL))
        print('NDA', NDA, sum(NDA))
        print('TOT', SAMPLE, sum(SAMPLE))

        print("\n")

        barWidth = 0.4
        br1 = np.arange(len(SAMPLE))
 
        # Make the plot
        bar1 = axes.bar(br1, AGN_perc, color ='midnightblue', width = barWidth,
        edgecolor ='grey', label ='AGNXY')
        axes.bar(br1, AGNX_perc, color ='dodgerblue', width = barWidth, edgecolor ='grey', label ='AGNX')
        axes.bar(br1, AGNY_perc, color ='aqua', width = barWidth, edgecolor ='grey', label ='AGNY')
        axes.bar(br1, UNC_perc, color ='springgreen', width = barWidth, edgecolor ='grey', label ='UNCXY')
        axes.bar(br1, UNCX_perc, color ='darkgreen', width = barWidth, edgecolor ='grey', label ='UNCX')
        axes.bar(br1, UNCY_perc, color ='limegreen', width = barWidth, edgecolor ='grey', label ='UNCY')
        axes.bar(br1, NOEL_perc, color ='silver', width = barWidth, edgecolor ='grey', label ='NOEL')
        axes.bar(br1, SF_perc, color ='mediumvioletred', width = barWidth, edgecolor ='grey', label ='SFGXY')
        axes.bar(br1, SFX_perc, color ='crimson', width = barWidth, edgecolor ='grey', label ='SFGX')
        axes.bar(br1, SFY_perc, color ='fuchsia', width = barWidth, edgecolor ='grey', label ='SFGY')
        
        yticks = np.arange(0, 101, 25)
        axes.set_ylabel(f'{y_name}', fontsize=17)
        axes.set_yticks(yticks, yticks, fontsize=17)
        axes.set_ylim(0, 110)
        
        # Adding Xticks

        k = 0
        for rect in bar1:
            height = rect.get_height()
            axes.text(rect.get_x() + rect.get_width() / 2.0, height, str(SAMPLE[k]), ha='center', va='bottom', fontsize='13')
            k+=1

def plotter_histo_WHAN(axes, BMS_condition, age, SC_WHAN, y_name, BMS, bins):
 
    SFG = empty()
    wAGN = empty()
    sAGN = empty()
    NOEL = empty()
    ELR = empty()
    LLR = empty()
    NER = empty()
    NDA = empty()
    SAMPLE = empty()

    SF_perc = empty()
    wAGN_perc = empty()
    sAGN_perc = empty()
    NOEL_perc = empty()
    ELR_perc = empty()
    LLR_perc = empty()
    RG_perc = empty()
    NDA_perc = empty()

    for i in range(len(bins)):
        for j in range(len(age)):
            if age[j] > bins[i][0] and age[j] <= bins[i][1] and int(BMS[j]) in BMS_condition: #below-MS=0 / MS galaxies=1 and item['BMS'] == 1:
                SAMPLE[i] += 1
                if SC_WHAN[j] == 'SFG':
                    SFG[i] += 1
                elif SC_WHAN[j] == 'wAGN':
                    wAGN[i] += 1
                elif SC_WHAN[j] == 'sAGN':
                    sAGN[i] += 1
                elif SC_WHAN[j] == 'UNC':
                    NOEL[i] += 1
                elif SC_WHAN[j] == 'ELR':
                    ELR[i] += 1
                elif SC_WHAN[j] == 'LLR':
                    LLR[i] += 1
                elif SC_WHAN[j] == 'NER':
                    NER[i] += 1
                elif SC_WHAN[j] in ['NDA', 'NDA0', 'NDA1']:
                    #NDA[i] += 1
                    SAMPLE[i] -= 1
        
    print(SAMPLE)
            
    for j in range(len(SAMPLE)):
        try:
            LLR_perc[j] = (LLR[j]/SAMPLE[j])*100
            RG_perc[j] = ((LLR[j] + NER[j])/SAMPLE[j])*100
            ELR_perc[j] = (LLR[j] + NER[j] + ELR[j])*100/SAMPLE[j]
            SF_perc[j] = (LLR[j] + NER[j] + ELR[j]+SFG[j])*100/SAMPLE[j]
            NDA_perc[j] = (NDA[j] + NER[j] + LLR[j] + ELR[j]+SFG[j])*100/SAMPLE[j]
            NOEL_perc[j] = (NDA[j] + NER[j] + LLR[j] + ELR[j]+SFG[j]+NOEL[j])*100/SAMPLE[j]
            wAGN_perc[j] = (NDA[j] + NER[j] + LLR[j] + ELR[j]+SFG[j]+NOEL[j] + wAGN[j])*100/SAMPLE[j]
            sAGN_perc[j] = (NDA[j] + NER[j] + LLR[j] + ELR[j]+SFG[j]+NOEL[j] + wAGN[j] + sAGN[j])*100/SAMPLE[j]
        except:
            pass
    
    barWidth = 0.4
    br1 = np.arange(len(SAMPLE))
 
        # Make the plot
    bar1 = axes.bar(br1, sAGN_perc, color ='midnightblue', width = barWidth,
    edgecolor ='grey', label ='sAGN')
    axes.bar(br1, wAGN_perc, color ='dodgerblue', width = barWidth, edgecolor ='grey', label ='wAGN')
    axes.bar(br1, NOEL_perc, color ='springgreen', width = barWidth, edgecolor ='grey', label ='UNC')
    axes.bar(br1, SF_perc, color ='mediumvioletred', width = barWidth, edgecolor ='grey', label ='SFG')
    axes.bar(br1, ELR_perc, color ='sandybrown', width = barWidth, edgecolor ='grey', label ='ELR')
    axes.bar(br1, RG_perc, color ='chocolate', width = barWidth, edgecolor ='grey', label ='NER')
    axes.bar(br1, LLR_perc, color ='maroon', width = barWidth, edgecolor ='grey', label ='LLR')
 
        # Adding Xticks

    yticks = np.arange(0, 101, 25)
    axes.set_ylabel(f'{y_name}', fontsize=17)
    axes.set_yticks(yticks, yticks, fontsize=17)
    axes.set_ylim(0, 110)
        

    k = 0
    for rect in bar1:
        height = rect.get_height()
        axes.text(rect.get_x() + rect.get_width() / 2.0, height, str(SAMPLE[k]), ha='center', va='bottom', fontsize='13')
        k+=1
    all = [sAGN, wAGN, NOEL, SFG, ELR, NER, LLR, SAMPLE]

    for group in all:
        for i, item in enumerate(group):
            print('${}$ & '.format(item), end='')
        print(sum(group), r' \\')

    print('\n')




