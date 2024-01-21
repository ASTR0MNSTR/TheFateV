import numpy as np
import matplotlib.pyplot as plt
from __legpars__ import *
from __stats__ import *
import pandas as pd

def theor_lines(axes, key):
    
    for ax in axes:
        if key == 'mdms':
            x = np.arange(7.9, 10.1, 0.01)
            a_all = 2.3152363967664957
            b_all = 5.276003555081468
            y = ((-(10**x)/10**9) / b_all) - a_all
            ax.plot(x, y, color='k', linestyle='solid')
        elif key == 'sfrsm':
            x = np.arange(6.9, 12, 0.1)
            ax.plot(x, (0.84 - 0.026*13.323023)*x - (6.51 - 0.11*13.323023), color='k', linestyle='dotted')
            ax.plot(x, (0.84 - 0.026*11.4074)*x - (6.51 - 0.11*11.4074), linestyle='dashed', color='k')
            ax.plot(x, (0.84 - 0.026*9.8615801)*x - (6.51 - 0.11*9.8615801), linestyle='solid', color='k')
    
            
def plotting(pars_dict):
    if 'err' in pars_dict.keys():
        cols = [pars_dict['x'], pars_dict['y'], pars_dict['err'], 'BPT', 'WHAN', 'P100_flux', 'P100_fluxerr']
    else:
        cols = [pars_dict['x'], pars_dict['y'], pars_dict['up'], pars_dict['down'], 'BPT', 'WHAN', 'P100_flux', 'P100_fluxerr']
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
    
    try:
        theor_lines([ax4, ax5], pars_dict['theor_lines'])
    except:
        pass
    
    if 'err' in pars_dict.keys():
        # ax4 = phys_plotter(ax4, DataFrame[pars_dict['x']], DataFrame[pars_dict['y']], DataFrame[pars_dict['y']] + DataFrame[pars_dict['err']], DataFrame[pars_dict['y']] - DataFrame[pars_dict['err']], DataFrame['BPT'], bids, 'BPT', True)
        ax4 = phys_plotter(ax4, DataFrame[pars_dict['x']], DataFrame[pars_dict['y']], DataFrame[pars_dict['y']], DataFrame[pars_dict['y']], DataFrame['BPT'], bids, 'BPT', True)
        # ax5 = phys_plotter(ax5, DataFrame[pars_dict['x']], DataFrame[pars_dict['y']], DataFrame[pars_dict['y']] + DataFrame[pars_dict['err']], DataFrame[pars_dict['y']] - DataFrame[pars_dict['err']], DataFrame['WHAN'], bids, 'WHAN', True)
        ax5 = phys_plotter(ax5, DataFrame[pars_dict['x']], DataFrame[pars_dict['y']], DataFrame[pars_dict['y']], DataFrame[pars_dict['y']], DataFrame['WHAN'], bids, 'WHAN', True)
    else:
        ax4 = phys_plotter(ax4, DataFrame[pars_dict['x']], DataFrame[pars_dict['y']], DataFrame[pars_dict['up']], DataFrame[pars_dict['down']], DataFrame['BPT'], bids, 'BPT', True)
        ax5 = phys_plotter(ax5, DataFrame[pars_dict['x']], DataFrame[pars_dict['y']], DataFrame[pars_dict['up']], DataFrame[pars_dict['down']], DataFrame['WHAN'], bids, 'WHAN', True)

    fig1.savefig(pars_dict['save_path'])

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
    
    errs = []
    means = []
    
    for item in class_list:
        X_plot = []
        Y_plot = []
        X, Y, err, length = monte_carlo(item[0], item[1], item[2], item[3], bids)
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
                list_names = list_names_BPT_1
            elif WHAN_or_BPT == 'WHAN':
                list_names = list_names_WHAN
            axis.scatter(-99, -99, alpha = 1, color=item[4][0], marker=item[4][1], s = 150, label=list_names[j])
        axis.legend(loc=3, fontsize="13")
    
    return axis   

def class_list_creator_w_err(x, y, up, down, AGN_keys, WHAN_or_BPT):
    if WHAN_or_BPT == 'BPT':
        keys = [['AGNXY'], ['AGNX'], ['UNCXY'], ['UNCX'], ['SFXY'], ['SFX'], ['NOEL']]
        colors_markers = [['midnightblue', 'P'], ['dodgerblue', 'P'], ['springgreen', 'H'], ['darkgreen', 'H'], ['mediumvioletred', '*'], ['deeppink', '*'], ['orchid', 'o']]
    elif WHAN_or_BPT == 'WHAN':
        keys = [['sAGN'], ['wAGN'], ['SF'], ['ELR'], ['NLR'], ['LLR'], ['sAGN', 'wAGN', 'SF', 'ELR', 'NLR', 'LLR']]
        colors_markers = [['midnightblue', 'P'], ['blue', 'P'], ['mediumvioletred', '*'], ['sandybrown', 'D'], ['chocolate', 'o'], ['maroon', 'o'], ['black', 'h']]
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

def class_list_creator_wo_err(x, y, age, AGN_keys, WHAN_or_BPT):
    if WHAN_or_BPT == 'BPT':
        keys = [['AGNXY'], ['AGNX'], ['UNCXY'], ['UNCX'], ['SFXY'], ['SFX'], ['NOEL']]
        colors_markers = [['midnightblue', 'P'], ['dodgerblue', 'P'], ['springgreen', 'H'], ['darkgreen', 'H'], ['mediumvioletred', '*'], ['deeppink', '*'], ['orchid', 'o']]
    elif WHAN_or_BPT == 'WHAN':
        keys = [['sAGN'], ['wAGN'], ['SF'], ['ELR'], ['NLR'], ['LLR']]
        colors_markers = [['midnightblue', 'P'], ['blue', 'P'], ['mediumvioletred', '*'], ['sandybrown', 'D'], ['chocolate', 'o'], ['maroon', 'o']]
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

def bin_stats(pars_dict):
    
    DataFrame = pd.read_csv(pars_dict['input_path'])
    
    fig, axs = plt.subplots(6, 1, figsize=(12, 18), tight_layout=True)
    plt.subplots_adjust(wspace=0, hspace=0)
    ax1 = axs[0]
    ax2 = axs[1]
    ax3 = axs[2]
    ax4 = axs[3]
    ax5 = axs[4]
    ax6 = axs[5]
    plotter_histo_BPT(ax3, [0], DataFrame[pars_dict['x']], DataFrame['BPT'], 'BMS, %', DataFrame['BMS'], pars_dict['bins'])

    plotter_histo_BPT(ax2, [1], DataFrame[pars_dict['x']], DataFrame['BPT'], 'MS, %', DataFrame['BMS'], pars_dict['bins'])

    plotter_histo_BPT(ax1, [0, 1], DataFrame[pars_dict['x']], DataFrame['BPT'], 'Total, %', DataFrame['BMS'], pars_dict['bins'])

    plotter_histo_WHAN(ax6, [0], DataFrame[pars_dict['x']], DataFrame['WHAN'], 'BMS, %', DataFrame['BMS'], pars_dict['bins'])

    plotter_histo_WHAN(ax5, [1], DataFrame[pars_dict['x']], DataFrame['WHAN'], 'MS, %', DataFrame['BMS'], pars_dict['bins'])

    plotter_histo_WHAN(ax4, [0,1], DataFrame[pars_dict['x']], DataFrame['WHAN'], 'Total, %', DataFrame['BMS'], pars_dict['bins'])

    ax2.legend(fontsize="10", loc='lower left')
    ax5.legend(fontsize="13", loc='lower left')

    top_axes = [ax2, ax3, ax5]
    for item in top_axes:
        item.tick_params(top=False, labeltop=False, bottom=False, labelbottom=False, right=True, direction='in')

    ax1.tick_params(top=False, labeltop=False, bottom=False, labelbottom=False, right=True, direction='in')
    #ax1.xaxis.set_label_position('top')
    ax1.set_title('BPT classification')
    ax4.tick_params(top=False, labeltop=False, bottom=False, labelbottom=False, right=True, direction='in')
    #ax4.xaxis.set_label_position('top')
    ax4.set_title('WHAN classification')

    ax6.set_xlabel(pars_dict['xlabel'])
    # ax6.set_xticks([r for r in range(6)], ['<0.05', '0.05-0.10', '0.10-0.15', '0.15-0.20', '0.20-0.25', '0.25-0.33'])
    ax6.set_xticks([r for r in range(len(pars_dict['bins_names']))], pars_dict['bins_names'])

    plt.savefig(pars_dict['save_path'])

def empty():
    return [0, 0, 0, 0, 0, 0, 0]

def plotter_histo_BPT(axes, BMS_condition, age, SC_BPT, y_name, BMS, bins):
        
        UNC = empty()
        UNCX = empty()
        UNCY = empty()
        SF = empty()
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
                    if SC_BPT[j] == 'SFXY':
                        SF[i] += 1
                    elif SC_BPT[j] == 'SFX':
                        SFX[i] += 1
                    elif SC_BPT[j] == 'SFY':
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
                SF_perc[j] = (SFY[j] + SFX[j]+SF[j])*100/SAMPLE[j]
                NDA_perc[j] = (NDA[j] + SFY[j] + SFX[j]+SF[j])*100/SAMPLE[j]
                NOEL_perc[j] = (NDA[j] + SFY[j] + SFX[j]+SF[j]+NOEL[j])*100/SAMPLE[j]
                UNCY_perc[j] = (NDA[j] + SFY[j] + SFX[j]+SF[j]+NOEL[j] + UNCY[j])*100/SAMPLE[j]
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
        print('SF', SF, sum(SF))
        print('SFX', SFX, sum(SFX))
        print('SFY', SFY, sum(SFY))
        print('NOEL', NOEL, sum(NOEL))
        print('NDA', NDA, sum(NDA))
        print('TOT', SAMPLE, sum(SAMPLE))

        print("\n")

        barWidth = 0.25
        br1 = np.arange(len(SAMPLE))
 
        # Make the plot
        bar1 = axes.bar(br1, AGN_perc, color ='midnightblue', width = barWidth,
        edgecolor ='grey', label ='AGNXY')
        axes.bar(br1, AGNX_perc, color ='dodgerblue', width = barWidth, edgecolor ='grey', label ='AGNX')
        axes.bar(br1, AGNY_perc, color ='aqua', width = barWidth, edgecolor ='grey', label ='AGNY')
        axes.bar(br1, UNC_perc, color ='springgreen', width = barWidth, edgecolor ='grey', label ='UNCXY')
        axes.bar(br1, UNCX_perc, color ='darkgreen', width = barWidth, edgecolor ='grey', label ='UNCX')
        axes.bar(br1, UNCY_perc, color ='limegreen', width = barWidth, edgecolor ='grey', label ='UNCY')
        axes.bar(br1, NOEL_perc, color ='white', width = barWidth, edgecolor ='grey', label ='NOEL')
        axes.bar(br1, SF_perc, color ='mediumvioletred', width = barWidth, edgecolor ='grey', label ='SFXY')
        axes.bar(br1, SFX_perc, color ='deeppink', width = barWidth, edgecolor ='grey', label ='SFX')
        axes.bar(br1, SFY_perc, color ='fuchsia', width = barWidth, edgecolor ='grey', label ='SFY')
        
        axes.set_ylabel(f'{y_name}')
        axes.set_ylim(0, 110)
        
        # Adding Xticks

        k = 0
        for rect in bar1:
            height = rect.get_height()
            axes.text(rect.get_x() + rect.get_width() / 2.0, height, str(SAMPLE[k]), ha='center', va='bottom')
            k+=1

def plotter_histo_WHAN(axes, BMS_condition, age, SC_WHAN, y_name, BMS, bins):
 
    SF = empty()
    wAGN = empty()
    sAGN = empty()
    NOEL = empty()
    ELR = empty()
    LLR = empty()
    NLR = empty()
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
                if SC_WHAN[j] == 'SF':
                    SF[i] += 1
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
                elif SC_WHAN[j] == 'NLR':
                    NLR[i] += 1
                elif SC_WHAN[j] in ['NDA', 'NDA0', 'NDA1']:
                    #NDA[i] += 1
                    SAMPLE[i] -= 1
        
    print(SAMPLE)
            
    for j in range(len(SAMPLE)):
        try:
            LLR_perc[j] = (LLR[j]/SAMPLE[j])*100
            RG_perc[j] = ((LLR[j] + NLR[j])/SAMPLE[j])*100
            ELR_perc[j] = (LLR[j] + NLR[j] + ELR[j])*100/SAMPLE[j]
            SF_perc[j] = (LLR[j] + NLR[j] + ELR[j]+SF[j])*100/SAMPLE[j]
            NDA_perc[j] = (NDA[j] + NLR[j] + LLR[j] + ELR[j]+SF[j])*100/SAMPLE[j]
            NOEL_perc[j] = (NDA[j] + NLR[j] + LLR[j] + ELR[j]+SF[j]+NOEL[j])*100/SAMPLE[j]
            wAGN_perc[j] = (NDA[j] + NLR[j] + LLR[j] + ELR[j]+SF[j]+NOEL[j] + wAGN[j])*100/SAMPLE[j]
            sAGN_perc[j] = (NDA[j] + NLR[j] + LLR[j] + ELR[j]+SF[j]+NOEL[j] + wAGN[j] + sAGN[j])*100/SAMPLE[j]
        except:
            pass
    
    barWidth = 0.25
    br1 = np.arange(len(SAMPLE))
 
        # Make the plot
    bar1 = axes.bar(br1, sAGN_perc, color ='midnightblue', width = barWidth,
    edgecolor ='grey', label ='sAGN')
    axes.bar(br1, wAGN_perc, color ='dodgerblue', width = barWidth, edgecolor ='grey', label ='wAGN')
    axes.bar(br1, NOEL_perc, color ='springgreen', width = barWidth, edgecolor ='grey', label ='UNC')
    axes.bar(br1, SF_perc, color ='mediumvioletred', width = barWidth, edgecolor ='grey', label ='SF')
    axes.bar(br1, ELR_perc, color ='sandybrown', width = barWidth, edgecolor ='grey', label ='ELR')
    axes.bar(br1, RG_perc, color ='chocolate', width = barWidth, edgecolor ='grey', label ='NLR')
    axes.bar(br1, LLR_perc, color ='maroon', width = barWidth, edgecolor ='grey', label ='LLR')
 
        # Adding Xticks

    axes.set_ylabel(f'{y_name}')
    axes.set_ylim(0, 110)

    k = 0
    for rect in bar1:
        height = rect.get_height()
        axes.text(rect.get_x() + rect.get_width() / 2.0, height, str(SAMPLE[k]), ha='center', va='bottom')
        k+=1
    all = [sAGN, wAGN, NOEL, SF, ELR, NLR, LLR, SAMPLE]

    for group in all:
        for i, item in enumerate(group):
            print('${}$ & '.format(item), end='')
        print(sum(group), r' \\')

    print('\n')




