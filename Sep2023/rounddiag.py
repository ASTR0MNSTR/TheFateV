import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from __plt__ import *

class Main:

    def merging_WHAN(self, list_obj):

        sAGN, wAGN, UNC, SF, ELR, LLR, RG = list_obj

        AGN = sAGN + wAGN
        UNC = UNC
        SF = SF
        ret = ELR + LLR + RG

        probs = [0, 0]
        if sAGN == 0 or wAGN == 0:
            probs[0] = 1
        if ELR == ret or RG == ret or LLR == ret:
            probs[1] = 1

        return [AGN, UNC, SF, ret], probs

    def my_autopct(pct):
        return (f'{pct:.2f}%') if pct > 5 else ''
    
    def short_WHAN_in(data):
        def my_format(pct):
            total = sum(data[0])
            val = int(round(pct*total/100.0))
            if data[0][1] == int(val) or data[0][2] == int(val) or pct == 100 or pct < 5:
                return ''
            elif (data[0][0] == int(val) and data[1][0] == 1) or (data[0][3] == int(val) and data[1][1] == 1):
                return ''
            else:
                return (f'{pct:.2f}%')
        return my_format

    def __init__(self, file):
        self.file = file
        self.dataframe = None

        self.WHAN_labels = ['sAGN', 'wAGN', 'UNC', 'SFG', 'ELR', 'LLR', 'NER']
        self.WHAN_colors = ['midnightblue', 'blue', 'springgreen', 'mediumvioletred', 'sandybrown', 'maroon', 'chocolate']

        self.WHAN_colors_merged = ['royalblue', 'lime', 'hotpink', 'brown']


        self.BPT_labels = ['AGN', 'UNC', 'SFG', 'NOEL']
        self.BPT_colors = ['midnightblue', 'springgreen', 'mediumvioletred', 'silver']
    

    def reading(self):
        usecols= ['WHAN', 'BPT', 'BMS']
        self.dataframe = pd.read_csv(self.file, usecols=usecols)
    
    def plotting(self):

        fig, axs = plt.subplots(4, 3, figsize=(12, 16), tight_layout=True)
        plt.subplots_adjust(wspace=0, hspace=0)
        # adjusting_plotting_pars()

        # [0, 0]
        # [1, 0] 1 - row
        # [0, 1] 1 - col

        Main.histo(self, axs[0, 0], ['AGNX'])
        Main.histo(self, axs[0, 1], ['AGNXY'])
        size = 0.45
        keys = ['AGNXY', 'AGNX', 'AGNY']
        axs[0, 2].pie(Main.sorting_forWHAN(self, keys), radius=1, labels=self.WHAN_labels, colors=self.WHAN_colors, autopct=Main.my_autopct, wedgeprops=dict(width=size, edgecolor='w'), pctdistance=0.8, labeldistance=1.1)
        axs[0, 2].pie(Main.merging_WHAN(self, Main.sorting_forWHAN(self, keys))[0], radius=1-size, colors=self.WHAN_colors_merged, autopct=Main.short_WHAN_in(Main.merging_WHAN(self, Main.sorting_forWHAN(self, keys))), wedgeprops=dict(width=size, edgecolor='w'))
        axs[0, 2].set(aspect='equal')
        
        #axs[0, 2].legend(title='WHAN classes:', loc="best", fontsize="13")
        axs[0, 2].set(aspect='equal')

######################################

        Main.histo(self, axs[1, 0], ['UNCX'])
        Main.histo(self, axs[1, 1], ['UNCXY'])
        Main.histo(self, axs[1, 2], ['UNCXY', 'UNCX', 'UNCY'])

######################################

        Main.histo(self, axs[2, 0], ['SFGX'])
        Main.histo(self, axs[2, 1], ['SFGXY'])
        Main.histo(self, axs[2, 2], ['SFGXY', 'SFGX', 'SFGY'])

######################################

        Main.histo(self, axs[3, 2], ['NOEL'])

        axs[0, 0].set_title(r"$\log \mathrm{([NII]/H\alpha)} \: \mathrm{(X)}$")
        axs[0, 1].set_title(r'$\log \mathrm{([NII]/H\alpha)} \: & \: \log \mathrm{([OIII]/H\beta)} \: \mathrm{(XY)}$')
        axs[0, 2].set_title('All objects (X, Y, XY)')
        axs[0, 0].set_ylabel('AGN (BPT)')
        axs[1, 0].set_ylabel('UNC (BPT)')
        axs[2, 0].set_ylabel('SFG (BPT)')
        axs[3, 2].set_ylabel('NOEL (BPT)')
        
        axs[3, 0].remove()
        axs[3, 1].remove()
        fig.savefig('./FIGURES_IN_PAPER/DIAG.pdf', dpi=70, transparent = True, bbox_inches = 'tight', pad_inches = 0.0001)
        plt.show()

    def sorting_forWHAN(self, keys):
        #SF
        SF = 0
        ELR = 0
        LLR = 0
        sAGN = 0
        wAGN = 0
        UNC = 0
        RG = 0
        NDA = 0
        self.total1 = 0
        for i in range(len(self.dataframe['BPT'])):
            if self.dataframe['BPT'][i] in keys:
                self.total1 += 1
                if self.dataframe['WHAN'][i] == 'SFG':
                    SF += 1
                elif self.dataframe['WHAN'][i] == 'ELR':
                    ELR += 1
                elif self.dataframe['WHAN'][i] == 'LLR':
                    LLR += 1
                elif self.dataframe['WHAN'][i] == 'sAGN':
                    sAGN += 1
                elif self.dataframe['WHAN'][i] == 'wAGN':
                    wAGN += 1
                elif self.dataframe['WHAN'][i] == 'UNC':
                    UNC += 1
                elif self.dataframe['WHAN'][i] == 'NER':
                    RG += 1
                elif self.dataframe['WHAN'][i] == 'NDA':
                    NDA += 1
                else:
                    print('AKHRANA, ATMENA')
        
        print(keys, self.total1)
        
        return [sAGN, wAGN, UNC, SF, ELR, LLR, RG]
    
    def my_level_list(data):
        list = []
        WHAN_labels = ['sAGN', 'wAGN', 'UNC', 'SFG', 'ELR', 'LLR', 'NER']
        for i in range(len(data)):
            if (data[i]*100/np.sum(data)) > 2 : #2%
                list.append(WHAN_labels[i])
            else:
                list.append('')
        return list

    def histo(self, figure, keys):
        size = 0.45
        patches, texts, autotexts = figure.pie(Main.sorting_forWHAN(self, keys), radius=1, labels=Main.my_level_list(Main.sorting_forWHAN(self, keys)), colors=self.WHAN_colors, autopct=Main.my_autopct, wedgeprops=dict(width=size, edgecolor='w'), pctdistance=0.8, labeldistance=1.1)
        [autotext.set_color('black') for autotext in autotexts]
        autotexts[0].set_color('white')
        autotexts[1].set_color('white')
        autotexts[-2].set_color('white')
        figure.pie(Main.merging_WHAN(self, Main.sorting_forWHAN(self, keys))[0], radius=1-size, colors=self.WHAN_colors_merged, autopct=Main.short_WHAN_in(Main.merging_WHAN(self, Main.sorting_forWHAN(self, keys))), wedgeprops=dict(width=size, edgecolor='w'))
        figure.set(aspect='equal')

if __name__ == '__main__':
    obj = Main('GAMA_ETG_OLA.csv')
    obj.reading()
    obj.plotting()