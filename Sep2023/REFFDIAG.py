import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from __plt__ import *
from __reader__ import *

def plotting_pars():
    plt.rcParams['font.size'] = 13

class Main:

    def __init__(self, path_to_data):
        self.path_to_data = path_to_data
        self.dataframe = None
        self.WHAN_labels = ['sAGN', 'wAGN', 'UNC', 'SFG', 'ELR', 'LLR', 'NER']
        self.WHAN_colors = ['midnightblue', 'blue', 'springgreen', 'mediumvioletred', 'sandybrown', 'maroon', 'chocolate']

        self.BPT_labels = ['AGNXY', 'AGNX', 'UNCXY', 'UNCX', 'UNCY', 'SFGXY', 'SFGX', 'SFGY', 'NOEL']
        self.BPT_colors = ['midnightblue', 'dodgerblue', 'springgreen', 'darkgreen', 'limegreen', 'mediumvioletred', 'crimson', 'fuchsia', 'silver']

        self.BPT_colors_merged = ['royalblue', 'lime', 'hotpink', 'w']
        self.WHAN_colors_merged = ['royalblue', 'lime', 'hotpink', 'brown']
        
        self.indexes_BPT = []
    
    def reading(self):
        usecols= ['WHAN', 'BPT', 'in_aperture', 'Z', 'mass_stellar_percentile50']
        # source_path = r'E:/LICENSE/ProgsData/main/GAMAforOleg.txt'
        # input_path = r'E:/databases/GAMA_ETG_OLA_R_r_1.csv'
        # merge_phys_databases(source_path, input_path, output_path)
        self.dataframe = pd.read_csv(self.path_to_data, usecols=usecols)
    
    def plotting(self):

        fig, axs = plt.subplots(2, 3, figsize=(12, 8), tight_layout=True)
        plt.subplots_adjust(wspace=0, hspace=0)
        # adjusting_plotting_pars()
        plotting_pars()

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
        ax1.set(title=r'$> 1 \; R_{eff}$' + ' (1202 gal.) \n [compact galaxies]')
        ax3.set(title=r'$< 1 \; R_{eff}$' + ' (679 gal.) \n [extended galaxies]')
        ax5.set(title='Total (1881 gal.)')     
        #self.ax1.legend(title = 'BPT: ', loc=2)
        #self.ax2.legend(title = 'WHAN:', loc=2)

        fig.savefig('./FIGURES_IN_PAPER_DR4/APERTURE_DIAG.pdf', dpi=300, transparent = True, bbox_inches = 'tight', pad_inches = 0.0001)
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
        redshift = []
        for i in range(len(self.dataframe['in_aperture'])):
            if self.dataframe['in_aperture'][i] in keys:
                self.total1 += 1
                redshift.append(float(self.dataframe['Z'][i]))
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
                    print(self.dataframe['WHAN'][i])
                    
        redshift = np.array(redshift)
        
        print(keys, self.total1, np.median(redshift) - np.percentile(redshift, 16), np.median(redshift), np.percentile(redshift, 84) - np.median(redshift))
        
        return [sAGN, wAGN, UNC, SF, ELR, LLR, RG]

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
        for i in range(len(self.dataframe['in_aperture'])):
            if self.dataframe['in_aperture'][i] in keys:
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
                    self.total1 -= 1
                elif self.dataframe['BPT'][i] == 'NOEL':
                    NOEL += 1
                else:
                    print(self.dataframe['BPT'][i])
        
        print(keys, self.total2)
        return [AGN, AGNX, UNC, UNCX, UNCY, SF, SFX, SFY, NOEL]
    
    def histo(self, figure, keys, kwarg):
        size = 0.5
        merged_WHAN = merging_WHAN(Main.sorting_forWHAN(self, keys))
        merged_BPT = merging_BPT(Main.sorting_forBPT(self, keys))
        if kwarg == 'WHAN':
            patches, texts, autotexts = figure.pie(Main.sorting_forWHAN(self, keys), radius=1, labels=my_level_list(Main.sorting_forWHAN(self, keys), 'WHAN'), colors=self.WHAN_colors, autopct=my_autopct_WHAN, wedgeprops=dict(width=size, edgecolor='w'), pctdistance=0.725, labeldistance=1.05)
            [autotext.set_color('black') for autotext in autotexts]
            autotexts[0].set_color('white')
            autotexts[1].set_color('white')
            autotexts[-2].set_color('white')
            # figure.pie(merged_WHAN, radius=1-size, colors=self.WHAN_colors_merged, autopct=short_WHAN_in(merging_WHAN(Main.sorting_forWHAN(self, keys))), wedgeprops=dict(width=size, edgecolor='w'))
            figure.pie(merged_WHAN, radius=1-size, colors=self.WHAN_colors_merged, wedgeprops=dict(width=size, edgecolor='w'))
            figure.set(aspect='equal')
        elif kwarg == 'BPT':
            patches, texts, autotexts = figure.pie(Main.sorting_forBPT(self, keys), radius=1, labels=my_level_list(Main.sorting_forBPT(self, keys), 'BPT'), colors=self.BPT_colors, autopct=my_autopct_BPT, wedgeprops=dict(width=size, edgecolor='w'), pctdistance=0.725, labeldistance=1.05)
            [autotext.set_color('black') for autotext in autotexts]
            autotexts[0].set_color('white')
            autotexts[3].set_color('white')
            # figure.pie(merged_BPT, radius=1-size, colors=self.BPT_colors_merged, autopct=short_BPT_in(merging_BPT(Main.sorting_forBPT(self, keys))), wedgeprops=dict(width=size, edgecolor='w'))
            figure.pie(merged_BPT, radius=1-size, colors=self.BPT_colors_merged, wedgeprops=dict(width=size, edgecolor='w'))
            figure.set(aspect='equal')

if __name__ == '__main__':
    obj = Main(r"E:\databases\GAMAs4\DETG_DR4.csv")
    obj.reading()
    obj.plotting()