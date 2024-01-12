import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

class Main:

    def my_autopct_BPT(pct):
        return (f'{pct:.2f}%') if pct > 0.1 else ''
    
    def my_autopct_WHAN(pct):
        return (f'{pct:.2f}%') if pct > 1 else ''

    def __init__(self, file):
        self.file = file
        self.dataframe = None
        self.WHAN_labels = ['sAGN', 'wAGN', 'UNC', 'SF', 'ELR', 'LLR', 'NLR']
        self.WHAN_colors = ['midnightblue', 'blue', 'springgreen', 'mediumvioletred', 'sandybrown', 'maroon', 'chocolate']

        self.BPT_labels = ['AGNXY', 'AGNX', 'UNCXY', 'UNCX', 'UNCY', 'SFXY', 'SFX', 'SFY', 'NOEL']
        self.BPT_colors = ['midnightblue', 'dodgerblue', 'springgreen', 'darkgreen', 'limegreen', 'mediumvioletred', 'deeppink', 'fuchsia', 'white']

        self.BPT_colors_merged = ['royalblue', 'lime', 'hotpink', 'w']
        self.WHAN_colors_merged = ['royalblue', 'lime', 'hotpink', 'brown']
    
    def reading(self):
        usecols= ['WHAN', 'BPT', 'Aperture_1_Reff']
        self.dataframe = pd.read_csv(self.file, usecols=usecols)
    
    def plotting(self):

        fig, axs = plt.subplots(2, 3, figsize=(12, 8), tight_layout=True)
        plt.subplots_adjust(wspace=0, hspace=0)

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
        ax1.set(title=r'$> 1 \; R_{eff}$ (969 gal.)')
        ax3.set(title=r'$< 1 \; R_{eff}$ (652 gal.)')
        ax5.set(title='Total (1621 gal.)')     
        #self.ax1.legend(title = 'BPT: ', loc=2)
        #self.ax2.legend(title = 'WHAN:', loc=2)

        fig.savefig('./FIGURES/APERTURE_DIAG_SDSS_r.pdf')
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
        for i in range(len(self.dataframe['Aperture_1_Reff'])):
            if self.dataframe['Aperture_1_Reff'][i] in keys:
                self.total1 += 1
                if self.dataframe['WHAN'][i] == 'SF':
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
                elif self.dataframe['WHAN'][i] == 'NLR':
                    RG += 1
                elif self.dataframe['WHAN'][i] == 'NDA':
                #    NDA += 1
                    self.total1 -= 1
                else:
                    print(self.dataframe['WHAN'][i])
        
        print(keys, self.total1)
        
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
        for i in range(len(self.dataframe['Aperture_1_Reff'])):
            if self.dataframe['Aperture_1_Reff'][i] in keys:
                self.total2 += 1
                if self.dataframe['BPT'][i] in ['SFXY']:
                    SF += 1
                elif self.dataframe['BPT'][i] in ['SFX']:
                    SFX += 1
                elif self.dataframe['BPT'][i] in ['SFY']:
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
                    print(self.dataframe['BPT'][i])
        
        print(keys, self.total2)
        return [AGN, AGNX, UNC, UNCX, UNCY, SF, SFX, SFY, NOEL]
    
    def merging_BPT(self, list_obj):
        AGN = list_obj[0] + list_obj[1]
        UNC = list_obj[2] + list_obj[3] + list_obj[4]
        SF =  list_obj[5] + list_obj[6] + list_obj[7]
        NOEL = list_obj[8]

        return [AGN, UNC, SF, NOEL]
    
    def merging_WHAN(self, list_obj):

        sAGN, wAGN, UNC, SF, ELR, LLR, RG = list_obj

        AGN = sAGN + wAGN
        UNC = UNC
        SF = SF
        ret = ELR + LLR + RG

        return [AGN, UNC, SF, ret]
    
    def my_level_list(data, kwarg):
        list = []
        if kwarg == 'WHAN':
            labels = ['sAGN', 'wAGN', 'UNC', 'SF', 'ELR', 'LLR', 'NLR']
        elif kwarg == 'BPT':
            labels = ['AGNXY', 'AGNX', 'UNCXY', 'UNCX', 'UNCY', 'SFXY', 'SFX', 'SFY', 'NOEL']

        for i in range(len(data)):
            if (data[i]*100/np.sum(data)) > 0.1: #2%
                list.append(labels[i])
            else:
                list.append('')
        return list
    
    def short_WHAN_in(data):
        def my_format(pct):
            total = sum(data)
            val = int(round(pct*total/100.0))
            if data[1] == int(val) or data[2] == int(val) or pct == 100 or pct < 1:
                return ''
            else:
                return (f'{pct:.2f}%')
        return my_format
    
    def short_BPT_in(data):
        def my_format(pct):
            total = sum(data)
            val = int(round(pct*total/100.0))
            if data[3] == int(val) or pct == 100 or pct < 1:
                return ''
            else:
                return (f'{pct:.2f}%')
        return my_format
    
    def histo(self, figure, keys, kwarg):
        size = 0.45
        if kwarg == 'WHAN':
            figure.pie(Main.sorting_forWHAN(self, keys), radius=1, labels=Main.my_level_list(Main.sorting_forWHAN(self, keys), 'WHAN'), colors=self.WHAN_colors, autopct=Main.my_autopct_WHAN, wedgeprops=dict(width=size, edgecolor='w'), pctdistance=0.8, labeldistance=1.1)
            figure.pie(Main.merging_WHAN(self, Main.sorting_forWHAN(self, keys)), radius=1-size, colors=self.WHAN_colors_merged, autopct=Main.short_WHAN_in(Main.merging_WHAN(self, Main.sorting_forWHAN(self, keys))), wedgeprops=dict(width=size, edgecolor='w'))
            figure.set(aspect='equal')
        elif kwarg == 'BPT':
            figure.pie(Main.sorting_forBPT(self, keys), radius=1, labels=Main.my_level_list(Main.sorting_forBPT(self, keys), 'BPT'), colors=self.BPT_colors, autopct=Main.my_autopct_BPT, wedgeprops=dict(width=size, edgecolor='w'), pctdistance=0.8, labeldistance=1.1)
            figure.pie(Main.merging_BPT(self, Main.sorting_forBPT(self, keys)), radius=1-size, colors=self.BPT_colors_merged, autopct=Main.short_BPT_in(Main.merging_BPT(self, Main.sorting_forBPT(self, keys))), wedgeprops=dict(width=size, edgecolor='w'))
            figure.set(aspect='equal')

if __name__ == '__main__':
    obj = Main('E:/databases/GAMA_ETG_OLA_R_r_1.csv')
    obj.reading()
    obj.plotting()