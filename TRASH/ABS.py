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
        self.dataframe = pd.read_csv(self.file)
    
    def plotting(self):

        fig, axs = plt.subplots(2, 4, figsize=(12, 6), tight_layout=True)
        plt.subplots_adjust(wspace=0, hspace=0)
        # adjusting_plotting_pars()
        
        ax1 = axs[0, 0]
        ax3 = axs[0, 1]
        ax5 = axs[0, 2]
        ax7 = axs[0, 3]

        ax2 = axs[1, 0]
        ax4 = axs[1, 1]
        ax6 = axs[1, 2]
        ax8 = axs[1, 3]

        Main.histo(self, ax1, 'OIII', 'OIII_er', 'BPT')
        Main.histo(self, ax2, 'OIII', 'OIII_er', 'WHAN')
        Main.histo(self, ax3, 'HB', 'HB_er', 'BPT')
        Main.histo(self, ax4, 'HB', 'HB_er', 'WHAN')
        Main.histo(self, ax5, 'NII', 'NII_er', 'BPT')
        Main.histo(self, ax6, 'NII', 'NII_er', 'WHAN')
        Main.histo(self, ax7, 'HA', 'HA_er', 'BPT')
        Main.histo(self, ax8, 'HA', 'HA_er', 'WHAN')

        ax1.set(title=r'$\mathrm{[OIII]}$' + ', 321 gal.')
        ax3.set(title=r'$\mathrm{H\beta}$' + ', 777 gal.')
        ax5.set(title=r'$\mathrm{[NII]}$' + ', 15 gal.')         
        ax7.set(title=r'$\mathrm{H\alpha}$' + ', 326 gal.')         
        #self.ax1.legend(title = 'BPT: ', loc=2)
        #self.ax2.legend(title = 'WHAN:', loc=2)

        fig.savefig('./FIGURES_IN_PAPER/ABS.pdf')
        plt.show()

    def sorting_forWHAN(self, flux_key, flux_er_key):
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
        for i in range(len(self.dataframe['BMS'])):
            if self.dataframe[flux_key][i] < (-2)*self.dataframe[flux_er_key][i]:
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
                #    NDA += 1
                    self.total1 -= 1
                else:
                    print('AKHRANA, ATMENA')
        
        print('WHAN', flux_key, self.total1)
        
        return [sAGN, wAGN, UNC, SF, ELR, LLR, RG]

    def sorting_forBPT(self, flux_key, flux_er_key):
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
            if self.dataframe[flux_key][i] < (-2)*self.dataframe[flux_er_key][i]:
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
        
        print('BPT', flux_key, self.total2)
        return [AGN, AGNX, UNC, UNCX, UNCY, SF, SFX, SFY, NOEL]
    
    def histo(self, figure, flux_key, flux_er_key, kwarg):
        size = 0.45
        if kwarg == 'WHAN':
            patches, texts, autotexts = figure.pie(Main.sorting_forWHAN(self, flux_key, flux_er_key), radius=1, labels=my_level_list(Main.sorting_forWHAN(self, flux_key, flux_er_key), 'WHAN'), colors=self.WHAN_colors, autopct=my_autopct_WHAN, wedgeprops=dict(width=size, edgecolor='w'), pctdistance=0.8, labeldistance=1.1)
            [autotext.set_color('black') for autotext in autotexts]
            autotexts[0].set_color('white')
            autotexts[1].set_color('white')
            # autotexts[-2].set_color('white')
            # figure.pie(Main.merging_WHAN(self, Main.sorting_forWHAN(self, flux_key, flux_er_key)), radius=1-size, colors=self.WHAN_colors_merged, autopct=short_WHAN_in(Main.merging_WHAN(self, Main.sorting_forWHAN(self, flux_key, flux_er_key))), wedgeprops=dict(width=size, edgecolor='w'))
            figure.pie(merging_WHAN(self, Main.sorting_forWHAN(self, flux_key, flux_er_key)), radius=1-size, colors=self.WHAN_colors_merged, wedgeprops=dict(width=size, edgecolor='w'))
            figure.set(aspect='equal')
        elif kwarg == 'BPT':
            patches, texts, autotexts = figure.pie(Main.sorting_forBPT(self, flux_key, flux_er_key), radius=1, labels=my_level_list(Main.sorting_forBPT(self, flux_key, flux_er_key), 'BPT'), colors=self.BPT_colors, autopct=my_autopct_BPT, wedgeprops=dict(width=size, edgecolor='w'), pctdistance=0.8, labeldistance=1.1)
            [autotext.set_color('black') for autotext in autotexts]
            autotexts[0].set_color('white')
            # autotexts[3].set_color('white')
            # figure.pie(Main.merging_BPT(self, Main.sorting_forBPT(self, flux_key, flux_er_key)), radius=1-size, colors=self.BPT_colors_merged, autopct=short_BPT_in(Main.merging_BPT(self, Main.sorting_forBPT(self, flux_key, flux_er_key))), wedgeprops=dict(width=size, edgecolor='w'))
            figure.pie(merging_BPT(self, Main.sorting_forBPT(self, flux_key, flux_er_key)), radius=1-size, colors=self.BPT_colors_merged, wedgeprops=dict(width=size, edgecolor='w'))
            figure.set(aspect='equal')

if __name__ == '__main__':
    obj = Main('GAMA_ETG_OLA.csv')
    obj.reading()
    obj.plotting()