import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from __algo__ import *
from __legpars__ import *
from __stats__ import *
from __plt__ import *

class Main:
    
    def __init__(self, path_to_data):
        self.path_to_data = path_to_data

        # self.radio_sources = [39495, 272962, 375530, 959586, 210108, 137037, 419441, 31076, 250557, 228288, 238593, 238593, 373280, 622622, 185507]
        
        self.dict_coord = {}
        
        self.color_dict = color_dict_BPT
        self.cd_WHAN = cd_WHAN
        self.cd_WHAN_leg = cd_WHAN_leg
        
    def getting_coords(self):
        df = pd.read_csv(self.path_to_data, sep=',')
        df.info()
        for index, row in df.iterrows():   
            AGN, X, pair_x_flags, Y, pair_y_flags, SC_WHAN, LOIII, LOIII_er, HA_ew, HA_ew_err, pair_HA = AGN_reg(row['OIIIR_FLUX'], row['OIIIR_FLUX_ERR'], row['HB_FLUX'], row['HB_FLUX_ERR'],
                                                                                                                 row['NIIR_FLUX'], row['NIIR_FLUX_ERR'], row['HA_FLUX'], row['HA_FLUX_ERR'], row['HA_EW'], row['HA_EW_ERR'], row['Z'])
            X_er = 0
            Y_er = 0
            self.dict_coord.update({row['SPECID'] : [X, X_er, pair_x_flags, Y, Y_er, pair_y_flags, HA_ew, HA_ew_err, pair_HA]})
            
    def plotting_arrows(self, ax, x, y, pair_x_flags, pair_y_flags, color, m_x, m_y, alpha):

        try:
            pair_x_flag = pair_x_flags[0]
        except:
            pair_x_flag = ''

        try:
            pair_y_flag = pair_y_flags[0]
        except:
            pair_y_flag = ''

        if pair_x_flag == 'down':
            pair_x_flag = 'left'
        elif pair_x_flag == 'up':
            pair_x_flag = 'right'

        coord_dict = {
            'down' : [0, m_y*(-0.07)],
            'left' : [m_x*(-0.07), 0],
            'up' : [0, m_y*(0.07)],
            'right' : [m_x*(0.07), 0]
        }
        try:
            ax.arrow(x, y, coord_dict[pair_x_flag][0], coord_dict[pair_x_flag][1], head_width=0.03,
            head_length=0.03, color=color, alpha=alpha)
        except:
            pass

        try:
            ax.arrow(x, y, coord_dict[pair_y_flag][0], coord_dict[pair_y_flag][1], head_width=0.03,
            head_length=0.03, color=color, alpha=alpha)
        except:
            pass

    def plotting_BPT(self):
        self.gs_top = plt.GridSpec(2, 3, wspace=0, hspace=0)
        #self.fig = plt.figure(figsize=(12, 12), tight_layout=True)
        self.fig = plt.figure(figsize=(20, 12))
        adjusting_figure_size(20, 12, l=1.2, r=1.8, b=0.6, t=0.3)
        adjusting_plotting_pars()

        self.ax4 = self.fig.add_subplot(self.gs_top[0,0])
        self.ax5 = self.fig.add_subplot(self.gs_top[0,1], sharey=self.ax4)
        self.ax_med_BPT = self.fig.add_subplot(self.gs_top[0,2], sharey=self.ax4)

        self.ax5.tick_params(top=True, labeltop=False, bottom=True, labelbottom=False, right=True, labelleft=False, left=True, direction='in', labelsize=20)
        self.ax_med_BPT.tick_params(top=True, labeltop=False, bottom=True, labelbottom=False, right=True, labelleft=False, left=True, direction='in', labelsize=20)
        self.ax4.tick_params(top=True, labeltop=False, bottom=True, labelbottom=False, right=True, labelleft=True, left=True, direction='in', labelsize=20)

        self.topaxes = [self.ax5, self.ax4, self.ax_med_BPT]
        for ax in self.topaxes:    
            # ax.set_xlabel(r'$\log{\mathrm{([NII]/H\alpha)}$')
            ax.set_yticks(np.arange(-1, 2.1, 0.5))
            ax.set_xticks(np.arange(-2, 1.2, 0.5))
            ax.set_xlim(-1.9, 1.2)
            ax.set_ylim(-1.2, 1.5) 
            X_1 = np.arange(-4, 0.4, 0.01)
            X_111 = np.arange(-4, 0, 0.01)
            X_11 = np.arange(-0.2, 1.5, 0.01)
            ax.plot(X_1, (0.61/(X_1 - 0.47)) + 1.19,
                      c='k', linewidth=3)  # Kauffmann, 2001
            ax.plot(X_111, (0.61/(X_111 - 0.05)) + 1.3,
                      c='k', linestyle='dashed', linewidth=3)  # Kewley, 2001
            # https://adsabs.harvard.edu/full/2003MNRAS.346.1055K
            # X = np.arange(-4, 0.4, 0.01)
            # ax.plot(X, (-30.787 + 1.1358*X + 0.27297*(X**2))*np.tanh(5.7409*X) - 31.093, linestyle='solid', linewidth=5, c='violet')
            # https://ui.adsabs.harvard.edu/abs/2006MNRAS.371..972S/abstract
            ax.plot(X_11, 1.01*X_11 + 0.48, c='black', linestyle='dotted', linewidth=3)
            ax.text(-1.5, 1.2, 'AGN')
            ax.text(0, -1, 'UNC', ha='center', va='center')
            ax.text(-1.5, 0.1, 'SFG')
            ax.text(0.5, -0.5, 'LINER')
            ax.set_box_aspect(1)

        self.ax4.set_ylabel(r'$\log{\mathrm{([OIII]/H\beta)}}$')

        norm = mpl.colors.Normalize(vmin=8.8,vmax=10.0)
        c_m = mpl.cm.viridis
        self.s_m = mpl.cm.ScalarMappable(cmap=c_m, norm=norm)
        self.s_m.set_array([])
        
        cmap = plt.cm.viridis  
        cmaplist = [cmap(i) for i in range(cmap.N)]

        cmap_segment = mpl.colors.LinearSegmentedColormap.from_list(
        'Custom cmap', cmaplist, cmap.N)

        # define the bins and normalize
        bounds = np.linspace(8.8, 10.0, 7)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        k = 0
        
        self.s_m_INT = mpl.cm.ScalarMappable(cmap=cmap_segment, norm=norm)
        self.s_m_INT.set_array([])
        
        # norm = mpl.colors.Normalize(vmin=8.8,vmax=10.25)
        # choose a colormap
        # c_m = mpl.cm.jet
        # create a ScalarMappable and initialize a data structure
        # self.s_m = mpl.cm.ScalarMappable(cmap=c_m, norm=norm)
        # self.s_m.set_array([])
        
        X = []
        Y = []
        AGE = []
        AGN_flags = []

        df = pd.read_csv(self.path_to_data)
        for index, row in df.iterrows(): 
            plots = self.dict_coord[row['SPECID']]
            
            AGN = row['BPT']
            SC_WHAN = row['WHAN']
            age = row['ager_percentile50']

            x = plots[0]
            x_er = plots[1]
            pair_x_flags = plots[2]
            y = plots[3]
            y_er = plots[4]
            pair_y_flags = plots[5]

            if x != -100 and x != -99 and y != -99 and y != 100 and AGN != 'NDA' and SC_WHAN != 'NDA':
                X.append(x)
                Y.append(y)
                AGE.append(age)
                AGN_flags.append(AGN)
                if len(pair_x_flags) == 0 and len(pair_y_flags) == 0:
                    try:
                        self.ax4.scatter(
                            x, y, s=5, color=self.color_dict[AGN][0], alpha=1)
                    #self.ax4.scatter(x, y, s=1.5, color=self.s_m.to_rgba(age), alpha=1)
                        self.ax5.scatter(x, y, s=5, color=self.cd_WHAN[SC_WHAN][0], alpha=1)
                        #self.ax4.scatter(x, y, s=self.cd_WHAN[SC_WHAN][1], color=self.cd_WHAN[SC_WHAN][0], marker =self.cd_WHAN[SC_WHAN][2], alpha=0.5)
                        #self.ax5.scatter(x, y, s=self.cd_WHAN[SC_WHAN][1], color=self.s_m.to_rgba(age), marker =self.cd_WHAN[SC_WHAN][2], alpha=0.5)
                        self.ax_med_BPT.scatter(x, y, s=5, color=self.s_m.to_rgba(age), alpha=0.15)
                        k += 1
                    except KeyError:
                        pass
                    # self.axes[i].scatter(plots[i][0], plots[i][1], s=1.5, color=self.s_m.to_rgba(age), alpha=1)
                    # self.axes[i].errorbar(plots[i][0], plots[i][1], xerr = plots[i][2], yerr = plots[i][3], fmt = 'o', color=self.s_m.to_rgba(age), markersize=2, alpha=0.2)
                else:
                    Main.plotting_arrows(self, self.ax5, x, y, pair_x_flags, pair_y_flags, self.cd_WHAN[SC_WHAN][0], m_x = 0.7, m_y = 0.7, alpha=1)
                    #Main.plotting_arrows(self, self.ax4, x, y, pair_x_flags, pair_y_flags, self.cd_WHAN[SC_WHAN][0], m_x = 1, m_y = 1)
                    #Main.plotting_arrows(self, self.ax5, x, y, pair_x_flags, pair_y_flags, self.s_m.to_rgba(age), m_x = 1, m_y = 1)
                    Main.plotting_arrows(self, self.ax4, x, y, pair_x_flags, pair_y_flags, self.color_dict[AGN][0], m_x = 0.7, m_y = 0.7, alpha=1)
                    Main.plotting_arrows(self, self.ax_med_BPT, x, y, pair_x_flags, pair_y_flags, self.s_m.to_rgba(age), m_x = 0.7, m_y = 0.7, alpha=0.5)
            

        class_list = class_list_creator_wo_err(X, Y, AGE, AGN_flags, 'BPT')
        
        for item in class_list:
            big_point_X, big_point_Y, big_point_age = median_position(item[0], item[1], item[2], [[8.8, 9.0], [9.0, 9.2], [9.2, 9.4], [9.4, 9.6], [9.6, 9.8], [9.8, 10.0]])
            self.ax_med_BPT.scatter(big_point_X, big_point_Y, s = 150, color = self.s_m_INT.to_rgba(big_point_age), marker=item[3][1], edgecolors='black')

        
        for key in BPT_color_plt:
            self.ax_med_BPT.scatter(-99, -99, alpha=1, color = BPT_color_plt[key][0], s = BPT_color_plt[key][1], marker= BPT_color_plt[key][2], label=key, edgecolors='black')
        
        self.ax_med_BPT.legend(loc=3)
        print('Number of points on BPT: ', k)

        #self.ax4.scatter(-99, -99, alpha= 1, color = 'midnightblue', label='AGN', s = 30, marker='o')
        #self.ax4.scatter(-99, -99, alpha= 1, color = 'springgreen', label='UNC', s = 30, marker='o')
        #self.ax4.scatter(-99, -99, alpha= 1, color = 'mediumvioletred', label='SF', s = 30, marker='o')
        #self.ax4.scatter(-99, -99, color='none', edgecolors='crimson', s=20, label='Abs. lines BPT')
        #self.ax4.scatter(-99, -99, color='none', edgecolors='black', s=20, label='Abs. lines WHAN')
        #self.ax5.scatter(-99, -99, color='none', edgecolors='crimson', s=20, label='Abs. lines BPT')
        #self.ax5.scatter(-99, -99, color='none', edgecolors='black', s=20, label='Abs. lines WHAN')

        # for key in self.cd_WHAN_leg.keys():
        #     self.ax5.scatter(-99, -99, alpha= 1, color = self.cd_WHAN_leg[key][0], marker = self.cd_WHAN_leg[key][2], s = self.cd_WHAN_leg[key][1], label=key)

        #self.ax5.legend(loc=3, fontsize="13")
        #self.ax4.legend(loc=3, fontsize="13")
        #self.fig.savefig('./FIGURES/BPT.pdf')
        #self.fig.savefig('BPT.pdf')
        # plt.show()

    
    def plotting_WHAN(self):

        self.ax6 = self.fig.add_subplot(self.gs_top[1,0])
        self.ax7 = self.fig.add_subplot(self.gs_top[1,1], sharey=self.ax6)
        self.ax_med_WHAN = self.fig.add_subplot(self.gs_top[1,2], sharey=self.ax6)

        # gs_top = plt.GridSpec(2, 1, hspace=0)
        # self.fig = plt.figure(figsize=(8,8), tight_layout=True)

        # self.ax6 = self.fig.add_subplot(gs_top[0,:])
        # self.ax7 = self.fig.add_subplot(gs_top[1,:], sharex=self.ax6)

        self.topaxes = [self.ax7, self.ax6, self.ax_med_WHAN]

        self.ax7.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, right=True, labelleft=False, left=True, direction='in', labelsize=20)
        self.ax_med_WHAN.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, right=True, labelleft=False, left=True, direction='in', labelsize=20)
        self.ax6.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, right=True, labelleft=True, left=True, direction='in', labelsize=20)
        #for ax in self.topaxes[1:]:
        #plt.setp(ax.get_xticklabels(), visible=False)

        for ax in self.topaxes:
            ax.set_xticks(np.arange(-2.0, 1.2, 0.5))
            ax.set_yticks(np.arange(-2, 3.0, 0.5))    
            ax.set_xlim([-1.9, 1.2])
            ax.set_ylim([-2, 2.7])
            ax.axhline(y = 0.47712, color = 'black', linestyle='dashed', linewidth=3)
            ax.axhline(y = -0.301, color = 'black', linestyle='dotted', linewidth=3)
            ax.text(0.8, 0, 'ELR')
            ax.text(0.8, -1.5, 'LLR')

            X_wAGN = np.arange(-0.4, 2.5, 0.01)
            ax.plot(X_wAGN, 0.77815125+X_wAGN*0, 'black', linewidth=3)
            ax.text(-1.5, 2, 'SFG')
            ax.text(0.6, 0.5, 'wAGN')
            ax.text(0.6, 2, 'sAGN')

            Y_sAGN = np.arange(0.47712, 3, 0.01)
            ax.plot(-0.4+Y_sAGN*0, Y_sAGN, 'black', linewidth=3)
            ax.set_xlabel(r"$\log \mathrm{([NII]/H\alpha)}$")
            ax.set_box_aspect(1)
        

        self.ax6.set_ylabel(r"$\log \mathrm{(EW_{\mathrm{H\alpha}})}$")
        
        k = 0
        X = []
        Y = []
        AGE = []
        AGN_flags = []

        df = pd.read_csv(self.path_to_data)
        for index, row in df.iterrows(): 
            age = row['ager_percentile50']
            plots = self.dict_coord[row['SPECID']]
            x = plots[0]
            pair_x_flags = plots[2]
            y = plots[6]
            pair_y_flags = plots[8]
            AGN = row['BPT']
            SC_WHAN = row['WHAN']

            if y >= -3 and y <= 3 and x >= -3.5 and x <= 2:
                k += 1
                
            if x != -100 and x != -99 and y != -99 and y != 100 and SC_WHAN != 'NDA' and SC_WHAN != 'NDA':
                X.append(x)
                Y.append(y)
                AGE.append(age)
                AGN_flags.append(SC_WHAN)
                if len(pair_x_flags) == 0 and len(pair_y_flags) == 0:
                    self.ax6.scatter(x, y, s=30, color=self.color_dict[AGN][0], alpha=1, marker=self.color_dict[AGN][2])
                    self.ax_med_WHAN.scatter(x, y, s=30, color=self.s_m.to_rgba(age), alpha=0.15, marker=self.color_dict[AGN][2])
                    self.ax7.scatter(x, y, s=30, color=self.cd_WHAN[SC_WHAN][0], marker = '.', alpha=1)
                    #self.ax6.scatter(x, y, s=self.cd_WHAN[SC_WHAN][1], color=self.cd_WHAN[SC_WHAN][0], marker = '.', alpha=0.5)
                    #self.ax7.scatter(x, y, s=self.cd_WHAN[SC_WHAN][1], color=self.s_m.to_rgba(age), marker = '.', alpha=0.5)
            
                else:
                    Main.plotting_arrows(self, self.ax7, x, y, pair_x_flags, pair_y_flags, self.cd_WHAN[SC_WHAN][0], m_y = 1, m_x = 0.8, alpha=1)
                    Main.plotting_arrows(self, self.ax_med_WHAN, x, y, pair_x_flags, pair_y_flags, self.s_m.to_rgba(age), m_y = 1, m_x = 0.8, alpha=0.5)
                    Main.plotting_arrows(self, self.ax6, x, y, pair_x_flags, pair_y_flags, self.color_dict[AGN][0], m_y = 1, m_x = 0.8, alpha = 1)
                    #Main.plotting_arrows(self, self.ax6, x, y, pair_x_flags, pair_y_flags, self.cd_WHAN[SC_WHAN][0], m_y = 1, m_x = 0.8)
                    
        class_list = class_list_creator_wo_err(X, Y, AGE, AGN_flags, 'WHAN')
        
        for j, item in enumerate(class_list):
            big_point_X, big_point_Y, big_point_age = median_position(item[0], item[1], item[2], [[8.8, 9.0], [9.0, 9.2], [9.2, 9.4], [9.4, 9.6], [9.6, 9.8], [9.8, 10.0]])
            self.ax_med_WHAN.scatter(big_point_X, big_point_Y, s = 150, color = self.s_m.to_rgba(big_point_age), marker=item[3][1], edgecolors='black')
        
        for key in WHAN_color_plt:
            self.ax_med_WHAN.scatter(-99, -99, alpha=1, color = WHAN_color_plt[key][0], s = WHAN_color_plt[key][1], marker= WHAN_color_plt[key][2], label=key, edgecolors='black')

        for key in self.cd_WHAN_leg.keys():
            self.ax7.scatter(-99, -99, alpha= 1, color = self.cd_WHAN_leg[key][0], marker = self.cd_WHAN_leg[key][2], s = 50, label=key)

        self.ax6.scatter(-99, -99, alpha= 1, color = 'midnightblue', label='AGNXY', s = 50)
        self.ax6.scatter(-99, -99, alpha= 1, color = 'dodgerblue', label='AGNX', s = 50)
        self.ax6.scatter(-99, -99, alpha= 1, color = 'springgreen', label='UNCXY', s = 50)
        self.ax6.scatter(-99, -99, alpha= 1, color = 'darkgreen', label='UNCX', s = 50)
        self.ax6.scatter(-99, -99, alpha= 1, color = 'mediumvioletred', label='SFGXY', s = 50)
        self.ax6.scatter(-99, -99, alpha= 1, color = 'crimson', label='SFGX', s = 50)
        self.ax6.legend(loc=3)
        self.ax7.legend(loc=3)
        self.ax_med_WHAN.legend(loc=3)
        
        
        #self.ax6.scatter(-99, -99, color='none', edgecolors='crimson', s=20, label='Abs. lines BPT')
        #self.ax6.scatter(-99, -99, color='none', edgecolors='black', s=20, label='Abs. lines WHAN')
        #self.ax7.scatter(-99, -99, color='none', edgecolors='crimson', s=20, label='Abs. lines BPT')
        #self.ax7.scatter(-99, -99, color='none', edgecolors='black', s=20, label='Abs. lines WHAN')

        # self.fig.subplots_adjust(right=0.89)
        cbar_ax = self.fig.add_axes([0.91, 0.05, 0.03, 0.925])
        self.fig.colorbar(self.s_m, cax=cbar_ax, label=r'$\log \mathrm{(age \: / \: yr)}$')
        #self.ax7.legend(loc=3, fontsize="13")
        
        self.fig.savefig('./FIGURES_IN_PAPER_DR4/BPT_WHAN.pdf', dpi=300, transparent = True, bbox_inches = 'tight', pad_inches = 0.0001)
        #self.fig.savefig('WHAN.pdf')

        # plt.show()


if __name__ == '__main__':
    obj = Main(r"E:\databases\GAMAs4\DETG_DR4.csv")
    obj.getting_coords()
    obj.plotting_BPT()
    obj.plotting_WHAN()