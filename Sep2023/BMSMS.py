import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from __plt__ import *

class Main:

    def __init__(self, file):
        self.file = file
        self.dataframe = None
        self.WHAN_labels = ['sAGN', 'wAGN', 'UNC', 'SFG', 'ELR', 'LLR', 'NER']
        self.WHAN_colors = ['midnightblue', 'blue', 'springgreen', 'mediumvioletred', 'sandybrown', 'maroon', 'chocolate']

        self.BPT_labels = ['AGNXY', 'AGNX', 'UNCXY', 'UNCX', 'UNCY', 'SFGXY', 'SFGX', 'SFGY', 'NOEL']
        self.BPT_colors = ['midnightblue', 'dodgerblue', 'springgreen', 'darkgreen', 'limegreen', 'mediumvioletred', 'crimson', 'fuchsia', 'silver']

        self.BPT_colors_merged = ['royalblue', 'lime', 'hotpink', 'w']
        self.WHAN_colors_merged = ['royalblue', 'lime', 'hotpink', 'brown']
    
    def reading(self):
        usecols= ['WHAN', 'BPT', 'BMS']
        self.dataframe = pd.read_csv(self.file, usecols=usecols)
    
    def plotting(self):

        fig, axs = plt.subplots(2, 3, figsize=(12, 8), tight_layout=True)
        plt.subplots_adjust(wspace=0, hspace=0)
        # adjusting_plotting_pars()

        ax1 = axs[0, 0]
        ax3 = axs[0, 1]
        ax5 = axs[0, 2]

        ax2 = axs[1, 0]
        ax4 = axs[1, 1]
        ax6 = axs[1, 2]

        Main.histo(self, ax1, [0], 'BPT')
        Main.histo(self, ax3, [1], 'BPT')
        Main.histo(self, ax5, [0, 1], 'BPT')
        Main.histo(self, ax2, [0], 'WHAN')
        Main.histo(self, ax4, [1], 'WHAN')
        Main.histo(self, ax6, [0, 1], 'WHAN')
        ax1.set(title='below-MS (1298 galaxies)')
        ax3.set(title='MS (697 galaxies)')
        ax5.set(title='Total (1995 galaxies)')     
        #self.ax1.legend(title = 'BPT: ', loc=2)
        #self.ax2.legend(title = 'WHAN:', loc=2)

        fig.savefig('./FIGURES_IN_PAPER/BMSMS.pdf', dpi=70, transparent = True, bbox_inches = 'tight', pad_inches = 0.0001)
        plt.show()

    def sorting_forWHAN(self, keys):
        #SF
        SF = 0
        ELR = 0
        LLR = 0
        sAGN = 0
        wAGN = 0
        UNC = 0
        NER = 0
        NDA = 0
        self.total1 = 0
        for i in range(len(self.dataframe['BMS'])):
            if self.dataframe['BMS'][i] in keys:
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
                    NER += 1
                elif self.dataframe['WHAN'][i] == 'NDA':
                #    NDA += 1
                    self.total1 -= 1
                else:
                    print('AKHRANA, ATMENA')
        
        print(keys, self.total1)
        print([sAGN, wAGN, UNC, SF, ELR, LLR, NER])
        
        return [sAGN, wAGN, UNC, SF, ELR, LLR, NER]

    def sorting_forBPT(self, keys):
        #SF
        SFX = 0
        SFY = 0
        SF = 0
        AGN = 0
        AGNX = 0
        AGNY = 0
        UNC = 0
        UNCX = 0
        UNCY = 0
        NDA = 0
        NOEL = 0

        self.total2 = 0
        for i in range(len(self.dataframe['BMS'])):
            if self.dataframe['BMS'][i] in keys:
                self.total2 += 1
                if self.dataframe['BPT'][i] in ['SFGXY']:
                    SF += 1
                elif self.dataframe['BPT'][i] in ['SFGX']:
                    SFX += 1
                elif self.dataframe['BPT'][i] in ['SFGY']:
                    SFY += 1
                elif self.dataframe['BPT'][i] in ['AGNXY']:
                    AGN += 1
                elif self.dataframe['BPT'][i] in ['AGNX']:
                    AGNX += 1
                elif self.dataframe['BPT'][i] in ['UNCXY']:
                    UNC += 1
                elif self.dataframe['BPT'][i] in ['UNCX']:
                    UNCX += 1
                elif self.dataframe['BPT'][i] in ['UNCY']:
                    UNCY += 1
                elif self.dataframe['BPT'][i] == 'NDA':
                #   NDA += 1
                    pass
                elif self.dataframe['BPT'][i] == 'NOEL':
                    NOEL += 1
                else:
                    print('AKHRANA, ATMENA')
        
        print(keys, self.total2)
        print([AGN, AGNX, UNC, UNCX, UNCY, SF, SFX, SFY, NOEL])
        return [AGN, AGNX, UNC, UNCX, UNCY, SF, SFX, SFY, NOEL]

    
    def histo(self, figure, keys, kwarg):
        size = 0.45
        merged_WHAN = merging_WHAN(self, Main.sorting_forWHAN(self, keys))
        merged_BPT = merging_BPT(self, Main.sorting_forBPT(self, keys))
        if kwarg == 'WHAN':
            patches, texts, autotexts = figure.pie(Main.sorting_forWHAN(self, keys), radius=1, labels=my_level_list(Main.sorting_forWHAN(self, keys), 'WHAN'), colors=self.WHAN_colors, autopct=my_autopct_WHAN, wedgeprops=dict(width=size, edgecolor='w'), pctdistance=0.8, labeldistance=1.1)
            [autotext.set_color('black') for autotext in autotexts]
            autotexts[0].set_color('white')
            autotexts[1].set_color('white')
            autotexts[-2].set_color('white')
            figure.pie(merged_WHAN, radius=1-size, colors=self.WHAN_colors_merged, autopct=short_WHAN_in(merging_WHAN(Main.sorting_forWHAN(self, keys))), wedgeprops=dict(width=size, edgecolor='w'))
            figure.set(aspect='equal')
        elif kwarg == 'BPT':
            patches, texts, autotexts = figure.pie(Main.sorting_forBPT(self, keys), radius=1, labels=my_level_list(Main.sorting_forBPT(self, keys), 'BPT'), colors=self.BPT_colors, autopct=my_autopct_BPT, wedgeprops=dict(width=size, edgecolor='w'), pctdistance=0.8, labeldistance=1.1)
            [autotext.set_color('black') for autotext in autotexts]
            autotexts[0].set_color('white')
            autotexts[3].set_color('white')
            figure.pie(merged_BPT, radius=1-size, colors=self.BPT_colors_merged, autopct=short_BPT_in(merging_BPT(Main.sorting_forBPT(self, keys))), wedgeprops=dict(width=size, edgecolor='w'))
            figure.set(aspect='equal')
            
            
if __name__ == '__main__':
    obj = Main('GAMA_ETG_OLA.csv')
    obj.reading()
    obj.plotting()