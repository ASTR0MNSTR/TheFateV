import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from __plt__ import *

def plotting_pars():
    plt.rcParams['font.size'] = 13

class Main:
    
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
        plotting_pars()
        # adjusting_plotting_pars()

        # [0, 0]
        # [1, 0] 1 - row
        # [0, 1] 1 - col

        Main.histo(self, axs[0, 0], ['AGNX'])
        Main.histo(self, axs[0, 1], ['AGNXY'])
        Main.histo(self, axs[0, 2], ['AGNXY', 'AGNX'])
        
        #axs[0, 2].legend(title='WHAN classes:', loc="best", fontsize="13")

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
        axs[0, 0].set_ylabel('AGN (BPT)', fontsize=14)
        axs[1, 0].set_ylabel('UNC (BPT)', fontsize=14)
        axs[2, 0].set_ylabel('SFG (BPT)', fontsize=14)
        axs[3, 2].set_ylabel('NOEL (BPT)', fontsize=14)
        
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
            if (data[i]*100/np.sum(data)) > 3 : #2%
                list.append(WHAN_labels[i])
            else:
                list.append('')
        return list

    def histo(self, figure, keys):
        merged_WHAN = merging_WHAN(Main.sorting_forWHAN(self, keys))
        size = 0.45
        patches, texts, autotexts = figure.pie(Main.sorting_forWHAN(self, keys), radius=1, labels=my_level_list(Main.sorting_forWHAN(self, keys), 'WHAN'), colors=self.WHAN_colors, autopct=my_autopct_WHAN, wedgeprops=dict(width=size, edgecolor='w'), pctdistance=0.75, labeldistance=1.08)
        [autotext.set_color('black') for autotext in autotexts]
        autotexts[0].set_color('white')
        autotexts[1].set_color('white')
        autotexts[-2].set_color('white')
        # figure.pie(merged_WHAN, radius=1-size, colors=self.WHAN_colors_merged, autopct=short_WHAN_in(merging_WHAN(Main.sorting_forWHAN(self, keys))), wedgeprops=dict(width=size, edgecolor='w'))
        figure.pie(merged_WHAN, radius=1-size, colors=self.WHAN_colors_merged, wedgeprops=dict(width=size, edgecolor='w'))
        figure.set(aspect='equal')
        
        
if __name__ == '__main__':
    obj = Main('GAMA_ETG_OLA.csv')
    obj.reading()
    obj.plotting()