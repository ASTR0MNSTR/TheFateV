import csv
import matplotlib.pyplot as plt
import math
import numpy as np
import matplotlib as mpl
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

from AGN_reg import *
from _global_ import *

# version 0.5 importing dust-correction


class hp:
    def log(figure):
        return math.log(figure, 10)

    def log_er(item):
        return [[abs(math.log(1-item, 10))], [abs(math.log(1+item, 10))]]

    def atten_coef(A_V, mode):
        R_V = 4.05
        OIII_c = 2.659*(-2.156 + (1.509/0.5007) -
                        (0.198/(0.5007**2)) + (0.011/(0.5007**3))) + R_V
        NII_c = 2.659*(-1.857 + (1.040/0.6583)) + R_V
        SII_c = 2.659*(-1.857 + (1.040/0.6732)) + R_V
        OI_c = 2.659*(-1.857 + (1.040/0.6300)) + R_V
        HB_c = 2.659*(-2.156 + (1.509/0.4861) -
                      (0.198/(0.4861**2)) + (0.011/(0.4861**3))) + R_V
        HD_c = 2.659*(-2.156 + (1.509/0.4101734) -
                      (0.198/(0.4101734**2)) + (0.011/(0.4101734**3))) + R_V
        HG_c = 2.659*(-2.156 + (1.509/0.4340472) -
                      (0.198/(0.4340472**2)) + (0.011/(0.4340472**3))) + R_V
        res = [OIII_c, NII_c, SII_c, OI_c, HB_c, HD_c, HG_c]
        res_out = [item*A_V/R_V for item in res]
        res_out_out = [10**(0.4*item) for item in res_out]
        return res_out_out

    def SFR(HA, HA_er, z):
        if HA == -99999.0 or HA_er < 0:
            return -99.9, -99.9

        cosmo = FlatLambdaCDM(H0=75 * u.km / u.s / u.Mpc,
                              Tcmb0=2.725 * u.K, Om0=0.3)
        d_mpc = float(str(cosmo.luminosity_distance(z)).split(' ')[0])
        if abs(HA) < 2*HA_er:
            HA += 2*HA_er
            SFR = 7.52*10**(-10)*4*np.pi*(d_mpc**2)*HA/(1+z)
            SFR_er = '-'
            return SFR, SFR_er
        elif HA > 2*HA_er:
            SFR = 7.52*10**(-10)*4*np.pi*(d_mpc**2)*HA/(1+z)
            SFR_er = SFR*(HA_er/HA)
            return SFR, SFR_er
        else:
            return -99.9, -99.9
    
    def ew_proc(HA_ew, HA_ew_err):
        pair_HA = []
        HA_ew_or = HA_ew
        SN = 2
        if HA_ew == -99999.0 or HA_ew_err < 0:
            return -99999.0, -99999.0, pair_HA, HA_ew_or
        else:
            if HA_ew > 0 and HA_ew < SN*HA_ew_err:
                pair_HA = ['down']
                HA_ew += SN*HA_ew_err
            elif HA_ew < 0 and HA_ew + SN*HA_ew_err <= 0:
                pair_HA = ['down']
                HA_ew = SN*HA_ew_err
            elif HA_ew < 0 and HA_ew + SN*HA_ew_err > 0:
                pair_HA = ['down']
                HA_ew += SN*HA_ew_err
            elif HA_ew > 0 and HA_ew >= SN*HA_ew_err:
                pair_HA = []
            else:
                print("Forgot me!")
        
        return math.log(HA_ew, 10), HA_ew_err, pair_HA, HA_ew_or

class Main(hp):
    def __init__(self, ola_file, gama_file, filename_out):
        self.filename_out = filename_out
        self.gama_file = gama_file
        self.ola_file = ola_file

        self.gama_dict = []
        self.sample_dict = {}
        self.headers = []  # headers for spectra_file

        self.flux_er_mod = []

        self.flux_er_mod9 = []
        self.flux_er_mod52 = []

        self.sample9 = []
        self.sample52 = []

        # EXPORTING

        self.OI = []
        self.SURV = []
        self.SURV_CODE = []
        self.IS_BEST = []
        self.IS_SBEST = []
        self.GAMAID = []
        self.SFR_HA = []
        self.SFR_HA_er = []
        self.coef = []
        self.met = []
        self.col1 = []
        self.col2 = []
        self.AGE = []
        self.AGN = []

        self.HA_l = []
        self.HB_l = []
        self.OI_l = []
        self.OIII_l = []
        self.NII_l = []
        self.SII_l = []
        self.HA_EW_l = []
        self.OII_l = []

        self.HA_l_er = []
        self.HB_l_er = []
        self.OI_l_er = []
        self.OIII_l_er = []
        self.NII_l_er = []
        self.SII_l_er = []
        self.HA_EW_l_er = []
        self.OII_l_er = []

        self.BMS = []
        self.NOEL_flags = []
        self.color_dict = {}

        self.RA = []
        self.DEC = []
        self.Z = []
        self.SPEC_ID = []
        self.SM = []
        self.SC_WHAN = []

        self.HdA = []
        self.HdA_er = []
        self.HdF = []
        self.HdF_er = []
        self.HgA = []
        self.HgA_er = []
        self.HgF = []
        self.HgF_er = []

        self.count_1 = 0
        self.count_2 = 0
        self.count_3 = 0
        self.count_4 = 0
        
        self.SN_HdF = []
        self.SN_HgF = []

        self.HgF_cont = []
        self.HdF_cont = []
        self.HA_cont = []

        #self.color_dict ={
        #    'AGN' : ['midnightblue', 15, '.'],
        #    'AGNX' : ['midnightblue', 27, '>'],
        #    'AGNY' : ['midnightblue', 27, '^'],
        #    'UNC' : ['springgreen', 3, '.'],
        #    'UNCX' : ['springgreen', 9, 's'],
        #    'UNCY' : ['springgreen', 9, 'D'],
        #    'SF' : ['mediumvioletred', 3, '.'],
        #    'SFX' : ['mediumvioletred', 9, '<'],
        #    'SFY' : ['mediumvioletred', 9, 'v'],
        #    'NOEL' : ['orchid', 3, '.'],
        #     'NDA' : ['slategrey', 3, '.']
        #}
        size = 3
        self.color_dict = color_dict

        #self.cd_WHAN = {
        #    'NOEL' : ['orchid', 24, '.'],
        #    'NDA' : ['gold', 24, '.'],
        #    'LLR' : ['maroon', 24, '.'],
        #    'ELR' : ['red', 24, '.'],
        #    'SF' : ['mediumvioletred', 24, '.'],
        #    'wAGN' : ['blue', 24, '.'],
        #    'sAGN' : ['midnightblue', 24, '.'],
        #}
        self.cd_WHAN = cd_WHAN

        self.cd_WHAN_leg = cd_WHAN_leg


    def ola_reading(self):
        with open(self.ola_file, 'r') as input:
            reader = csv.DictReader(input)
            for row in reader:
                self.sample_dict.update({int(row['CATAID_1']): int(row['below=0/MS=1'])})

    def samples_get(self):
        #data_dict = []
        pass
        #with open(self.sample_file, 'r') as csvfile:
        #    reader = csv.DictReader(csvfile)
        #    for row in reader:
        #        data_dict.append(
        #            {'OI': row['OI'], 'GAMAID': row['GAMAID'], 'age': row['age']})

        #for galaxy in data_dict:
        #    if galaxy['OI'] == '':
        #        IND = 0
        #    else:
        #        IND = int(galaxy['OI'])
        #    self.sample_dict.update(
        #        {int(galaxy['GAMAID']): [float(galaxy['age']), -1, IND, -99.9]})

    def gama_reading(self):
        count = 0
        with open(self.gama_file, 'r') as input:
            while True:
                line = input.readline()
                if not line:
                    break
                count += 1
                self.flux_er_mod.append(Main.dividing(self, line))
        print(count)
        print(len(self.flux_er_mod))

    def dividing(self, line):
        line_out = line.strip()
        galaxy_pars = line_out.split()
        CATAID = int(galaxy_pars[1])
        #bms = -1
        if CATAID in self.sample_dict.keys():
            bms = self.sample_dict[CATAID]
        else:
            return 1

        Z = float(galaxy_pars[4])
        RA = float(galaxy_pars[2])
        DEC = float(galaxy_pars[3])
        SPEC_ID = galaxy_pars[0]
        HA = float(galaxy_pars[82])
        HA_er = float(galaxy_pars[81])
        HA_cont = float(galaxy_pars[78])
        HB = float(galaxy_pars[182])
        HB_er = float(galaxy_pars[181])
        OIII = float(galaxy_pars[167])
        OIII_er = float(galaxy_pars[166])  # OIIIR, 5007
        NII = float(galaxy_pars[72])
        NII_er = float(galaxy_pars[71])  # NIIR, 6583

        HdA = float(galaxy_pars[262])
        HdA_er = float(galaxy_pars[261])

        HdF = float(galaxy_pars[257])
        HdF_er = float(galaxy_pars[256])
        HdF_cont = float(galaxy_pars[253])

        HgA = float(galaxy_pars[232])
        HgA_er = float(galaxy_pars[231])

        HgF = float(galaxy_pars[227])
        HgF_er = float(galaxy_pars[226])
        HgF_cont = float(galaxy_pars[223])

        HA_EW = float(galaxy_pars[80])
        HA_EW_ERR = float(galaxy_pars[79])
        
        OII = float(galaxy_pars[267]) # OII 3727
        OII_er = float(galaxy_pars[266])

        OI = float(galaxy_pars[97]) #OI, 6300 A
        OI_er = float(galaxy_pars[96])

        SII = float(galaxy_pars[67]) #SIIB, 6716 A
        SII_er = float(galaxy_pars[66])


        HA_EW, HA_EW_ERR, pair_HA, HA_EW_OR = hp.ew_proc(HA_EW, HA_EW_ERR)

        # dust_correction module:
        
        AGN, X, X_er, pair_x_flags, Y, Y_er, pair_y_flags, SC_WHAN = AGN_reg(OIII, OIII_er, HB, HB_er, NII, NII_er, HA, HA_er, HA_EW, HA_EW_ERR, pair_HA)

        SFR_HA, SFR_HA_er = hp.SFR(HA, HA_er, Z)

        SURV = galaxy_pars[6]
        SURV_CODE = galaxy_pars[7]
        IS_BEST = galaxy_pars[8]
        IS_SBEST = galaxy_pars[9]

        FLUXES = [HA, HB, OIII, NII, SII, OI, OII]
        FLUXES_ER = [HA_er, HB_er, OIII_er, NII_er, SII_er, OI_er, OII_er]

        #par_met_0 = [[OIII, OIII_er], [HB, HB_er]]
        #par_met_1 = [[NII, NII_er], [HA, HA_er]]

        #par_mets = [hp.indexing(par) for par in [par_met_0, par_met_1]]
        #met_out = hp.metallicity(hp.val_mist([par_mets[0], par_mets[1]]))

        res_out = [X, X_er, pair_x_flags, Y, Y_er, pair_y_flags]

        kwargs = [SURV, SURV_CODE, IS_BEST, IS_SBEST, CATAID,
                  [SFR_HA, SFR_HA_er], AGN, FLUXES, FLUXES_ER, bms, [RA, DEC, Z, SPEC_ID], [X, pair_x_flags, HA_EW, pair_HA], SC_WHAN,
                  [HdA, HdA_er, HdF, HdF_er, HgA, HgA_er, HgF, HgF_er], [HA_EW_OR, HA_EW_ERR], HgF_cont, HdF_cont, HA_cont]
        res_out_out = [res_out, kwargs]
        return res_out_out

    def sorting(self):
        counter_all = 0
        counter_best = 0
        for pars in self.flux_er_mod:
            if type(pars) == list:
                counter_all += 1
                if pars[-1][2] == 'true' and pars[-1][3] == 'true' and (pars[-1][0] == 'GAMA' or pars[-1][0] == 'SDSS') and pars[-1][6] != 'NDA': #!!!!
                #if pars[-1][2] == 'true' and pars[-1][3] == 'true' and pars[-1][6] != 'NDA': #!!!!
                #if pars[-1][2] == 'true' and pars[-1][3] == 'true':
                #if pars[-3][6] != 'NDA': #!!!!
                    Main.exporting(self, pars)
                    counter_best += 1
                    self.flux_er_mod9.append(pars)
            else:
                pass

        print('All processed galaxies: ', counter_all)
        print('Galaxies with proper data: ', counter_best)
        Main.file_out(self)

    def exporting(self, pars):
        self.SURV.append(pars[-1][0])
        self.SURV_CODE.append(pars[-1][1])
        self.IS_BEST.append(pars[-1][2])
        self.IS_SBEST.append(pars[-1][3])
        self.GAMAID.append(pars[-1][4])
        self.SFR_HA.append(pars[-1][5][0])
        self.SFR_HA_er.append(pars[-1][5][1])

        # changed line
        self.AGN.append(pars[-1][6])

        self.HA_l.append(pars[-1][7][0])
        self.HB_l.append(pars[-1][7][1])
        self.OIII_l.append(pars[-1][7][2])
        self.NII_l.append(pars[-1][7][3])
        self.SII_l.append(pars[-1][7][4])
        self.OI_l.append(pars[-1][7][5])
        self.OII_l.append(pars[-1][7][6])

        self.HA_l_er.append(pars[-1][8][0])
        self.HB_l_er.append(pars[-1][8][1])
        self.OIII_l_er.append(pars[-1][8][2])
        self.NII_l_er.append(pars[-1][8][3])
        self.SII_l_er.append(pars[-1][8][4])
        self.OI_l_er.append(pars[-1][8][5])
        self.OII_l_er.append(pars[-1][8][6])

        self.BMS.append(pars[-1][9])
        self.RA.append(pars[-1][10][0])
        self.DEC.append(pars[-1][10][1])
        self.Z.append(pars[-1][10][2])
        self.SPEC_ID.append(pars[-1][10][3])
        self.SC_WHAN.append(pars[-1][12])

        self.HdA.append(pars[-1][13][0])
        self.HdA_er.append(pars[-1][13][1])
        self.HdF.append(pars[-1][13][2])
        self.HdF_er.append(pars[-1][13][3])
        self.HgA.append(pars[-1][13][4])
        self.HgA_er.append(pars[-1][13][5])
        self.HgF.append(pars[-1][13][6])
        self.HgF_er.append(pars[-1][13][7])

        self.SN_HdF.append(pars[-1][13][2]/pars[-1][13][3])
        self.SN_HgF.append(pars[-1][13][6]/pars[-1][13][7])

        self.HA_EW_l.append(pars[-1][14][0])
        self.HA_EW_l_er.append(pars[-1][14][1])

        self.HgF_cont.append(pars[-1][15])
        self.HdF_cont.append(pars[-1][16])
        self.HA_cont.append(pars[-1][17])

    def file_out(self):

        Dict = {
            'SPEC_ID' : self.SPEC_ID,
            'GAMAID': self.GAMAID,
            'SURVEY': self.SURV,
            'IS_BEST': self.IS_BEST,
            'IS_SBEST': self.IS_SBEST,
            'RA' : self.RA,
            'DEC' : self.DEC,
            'Z' : self.Z,
            'BPT': self.AGN,
            'WHAN' : self.SC_WHAN,
            'SN_HgF' : self.SN_HgF,
            'SN_HdF' : self.SN_HdF,
            'HgF' : self.HgF,
            'HgF_cont' : self.HgF_cont,
            'HgF_er' : self.HgF_er,
            'HdF' : self.HdF,
            'HdF_cont' : self.HdF_cont,
            'HdF_er' : self.HdF_er,
            'HA_EW' : self.HA_EW_l,
            'HA_EW_ERR' : self.HA_EW_l_er,
            'HA': self.HA_l,
            'HA_cont' : self.HA_cont,
            'HA_er': self.HA_l_er,
            'HB': self.HB_l,
            'HB_er': self.HB_l_er,
            'OIII': self.OIII_l,
            'OIII_er': self.OIII_l_er,
            'NII': self.NII_l,
            'NII_er': self.NII_l_er,
            'HdA' : self.HdA,
            'HdA_er' : self.HdA_er,
            'HgA' : self.HgA,
            'HgA_er' : self.HgA_er, 
            'SII' : self.SII_l,
            'SII_er' : self.SII_l_er,
            'OI' : self.OI_l,
            'OI_er' : self.OI_l_er,
            'OII' : self.OII_l,
            'OII_er' : self.OII_l_er,
            'BMS': self.BMS
        }

        df = pd.DataFrame(Dict)
        df.to_csv(self.filename_out, index=False)

    def plotting_arrows(self, ax, dict, x, y, pair_x_flags, pair_y_flags, mode, m_x, m_y):

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
            head_length=0.03, color=dict[mode][0], alpha=1)
        except:
            pass

        try:
            ax.arrow(x, y, coord_dict[pair_y_flag][0], coord_dict[pair_y_flag][1], head_width=0.03,
            head_length=0.03, color=dict[mode][0], alpha=1)
        except:
            pass

    
    def plotting_arrows_WHAN(self, ax, dict, x, y, pair_x_flags, pair_y_flags, mode, m_x, m_y):

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
            head_length=0.03, color=dict[mode][0], alpha=1)
        except:
            pass

        try:
            ax.arrow(x, y, coord_dict[pair_y_flag][0], coord_dict[pair_y_flag][1], head_width=0.03,
            head_length=0.03, color=dict[mode][0], alpha=1)
        except:
            pass


    def plotting_BPT(self):
        gs_top = plt.GridSpec(1, 2, wspace=0)
        self.fig = plt.figure(figsize=(12, 6))

        self.ax4 = self.fig.add_subplot(gs_top[:,0])
        self.ax5 = self.fig.add_subplot(gs_top[:,1], sharey=self.ax4)

        self.ax5.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, right=True, labelleft=False, left=True, direction='in')
        self.ax4.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, right=True, labelleft=True, left=True, direction='in')

        self.topaxes = [self.ax5, self.ax4]
        for ax in self.topaxes:    
            ax.set_xlabel(r'$log(N[II]/H\alpha)$')
            ax.set_xlim(-2, 0.8)
            ax.set_ylim(-1.5, 1.5) 
            X_1 = np.arange(-4, 0.4, 0.01)
            X_111 = np.arange(-4, 0, 0.01)
            X_11 = np.arange(-0.45, 0.6, 0.01)
            ax.plot(X_1, (0.61/(X_1 - 0.47)) + 1.19,
                      c='k')  # Kewley, 2001
            ax.plot(X_111, (0.61/(X_111 - 0.05)) + 1.3,
                      c='k', linestyle='dashed')  # Kauffman, 2003
            # https://adsabs.harvard.edu/full/2003MNRAS.346.1055K
            ax.plot(X_11, 1.01*X_11 + 0.48, c='k', linestyle='dotted')
            ax.text(-1, 2, 'AGN', fontweight='semibold')
            #self.ax1.text(0.1, 0.2, 'LINER')
            ax.text(0, -1.5, 'C', fontweight='semibold')
            ax.text(-1.7, -2, 'SF', fontweight='semibold')
            ax.text(1, 0.5, 'LINER', fontweight='semibold')
            ax.set_yticks(np.arange(-1.5, 1.5, 0.5))
            ax.set_xticks(np.arange(-2, 1, 0.5))

        self.ax4.set_ylabel(r'$log(O[III]/H\beta)$')

        ages = []
        for pars in self.flux_er_mod9:
            ages.append(pars[1])

        # norm = mpl.colors.Normalize(vmin=8.8,vmax=10.25)
        # choose a colormap
        # c_m = mpl.cm.jet
        # create a ScalarMappable and initialize a data structure
        # self.s_m = mpl.cm.ScalarMappable(cmap=c_m, norm=norm)
        # self.s_m.set_array([])
        k = 0

        for pars in self.flux_er_mod9:
            plots = pars[0]
            age = pars[-2]
            AGN = pars[-1][6]
            SC_WHAN = pars[-1][12]

            x = plots[0]
            x_er = plots[1]
            pair_x_flags = plots[2]
            y = plots[3]
            y_er = plots[4]
            pair_y_flags = plots[5]

            #if y >= -3 and y <= 3 and x >= -3.5 and x <= 2:
            if AGN[-1] == '!':
                AGN = AGN[:-1]
                self.ax4.scatter(x, y, color='none', edgecolors='crimson', s=20)
                self.ax5.scatter(x, y, color='none', edgecolors='crimson', s=20)
                

            if SC_WHAN[-1] == '!':
                SC_WHAN = SC_WHAN[:-1]
                self.ax4.scatter(x, y, color='none', edgecolors='black', s=50)
                self.ax5.scatter(x, y, color='none', edgecolors='black', s=50)

            if len(pair_x_flags) == 0 and len(pair_y_flags) == 0:
                try:
                    self.ax4.scatter(
                        x, y, s=1.5, color=self.color_dict[AGN][0], alpha=1)
                    self.ax5.scatter(x, y, s=self.cd_WHAN[SC_WHAN][1], color=self.cd_WHAN[SC_WHAN][0], marker =self.cd_WHAN[SC_WHAN][2], alpha=1)

                    k += 1
                except KeyError:
                    pass
                    # self.axes[i].scatter(plots[i][0], plots[i][1], s=1.5, color=self.s_m.to_rgba(age), alpha=1)
                    # self.axes[i].errorbar(plots[i][0], plots[i][1], xerr = plots[i][2], yerr = plots[i][3], fmt = 'o', color=self.s_m.to_rgba(age), markersize=2, alpha=0.2)
            else:
                Main.plotting_arrows(
                    self, self.ax5, self.cd_WHAN, x, y, pair_x_flags, pair_y_flags, SC_WHAN, m_x = 1, m_y = 1)
                Main.plotting_arrows(
                    self, self.ax4, self.color_dict, x, y, pair_x_flags, pair_y_flags, AGN, m_x = 1, m_y = 1)
            
            #self.ax1.scatter(x, y, s=self.cd_WHAN[SC_WHAN][1], color=self.cd_WHAN[SC_WHAN][0], marker =self.cd_WHAN[SC_WHAN][2], alpha=1)
                
            #if mod > 0:
                #self.ax1.text(x, y, pars[-1])
            #    self.ax1.scatter(x, y, color='black', marker="x")
            #if mod == 0:
            #    self.ax1.scatter(x, y, color='black', marker="x")


        # self.fig.subplots_adjust(right=0.85)
        # cbar_ax = self.fig.add_axes([0.85, 0.15, 0.05, 0.7])
        # for map in self.axes:
        #    self.fig.colorbar(self.s_m, cax=cbar_ax)

        print('Number of points on BPT: ', k)

        self.ax4.scatter(-99, -99, alpha= 1, color = 'midnightblue', label='AGN', s = 30, marker='.')
        self.ax4.scatter(-99, -99, alpha= 1, color = 'springgreen', label='UNC', s = 30, marker='.')
        self.ax4.scatter(-99, -99, alpha= 1, color = 'mediumvioletred', label='SF', s = 30, marker='.')
        self.ax4.scatter(-99, -99, color='none', edgecolors='crimson', s=20, label='Abs. lines BPT')
        self.ax4.scatter(-99, -99, color='none', edgecolors='black', s=20, label='Abs. lines WHAN')
        self.ax5.scatter(-99, -99, color='none', edgecolors='crimson', s=20, label='Abs. lines BPT')
        self.ax5.scatter(-99, -99, color='none', edgecolors='black', s=20, label='Abs. lines WHAN')

        self.ax4.set_box_aspect(1)
        self.ax5.set_box_aspect(1)

        for key in self.cd_WHAN_leg.keys():
            self.ax5.scatter(-99, -99, alpha= 1, color = self.cd_WHAN_leg[key][0], marker = self.cd_WHAN_leg[key][2], s = self.cd_WHAN_leg[key][1], label=key)
        self.ax5.legend(loc=3)

        self.ax4.legend(loc=3)
        self.fig.savefig('./FIGURES/BPT.pdf')
        # plt.show()

    
    def plotting_WHAN(self):
        

        gs_top = plt.GridSpec(2, 1, hspace=0)
        self.fig = plt.figure(figsize=(8,8))

        self.ax4 = self.fig.add_subplot(gs_top[0,:])
        self.ax5 = self.fig.add_subplot(gs_top[1,:], sharex=self.ax4)

        self.topaxes = [self.ax5, self.ax4]

        self.ax5.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, right=True, direction='in')
        self.ax4.tick_params(top=True, labeltop=False, bottom=True, labelbottom=False, right=True, direction='in')
        #for ax in self.topaxes[1:]:
        #plt.setp(ax.get_xticklabels(), visible=False)

        for ax in self.topaxes:    
            ax.set_ylabel(r"$log(EW_{H\alpha})$")
            ax.set_xlim([-2, 2.5])
            ax.set_ylim([-3, 3])
            ax.axhline(y = 0.47712, color = 'black', linestyle='dashed')
            ax.axhline(y = -0.301, color = 'black', linestyle='dotted')
            ax.text(1.5, 0, 'ELR')
            ax.text(1.5, -2, 'LLR')

            X_wAGN = np.arange(-0.4, 2.5, 0.01)
            ax.plot(X_wAGN, 0.77815125+X_wAGN*0, 'black')
            ax.text(-1.5, 2, 'SF')
            ax.text(1.5, 0.5, 'wAGN')
            ax.text(1.5, 2, 'sAGN')

            Y_sAGN = np.arange(0.47712, 3, 0.01)
            ax.plot(-0.4+Y_sAGN*0, Y_sAGN, 'black')
        
        self.ax5.set_yticks(np.arange(-3, 2.1, 1))

        self.ax5.set_xlabel(r"$log(N[II]/H\alpha)$")

        k = 0
        # norm = mpl.colors.Normalize(vmin=8.8,vmax=10.25)
        # choose a colormap
        # c_m = mpl.cm.jet
        # create a ScalarMappable and initialize a data structure
        # self.s_m = mpl.cm.ScalarMappable(cmap=c_m, norm=norm)
        # self.s_m.set_array([])

        for pars in self.flux_er_mod9:
            coord = pars[-1][11]
            x = coord[0]
            pair_x_flags = coord[1]
            y = coord[2]
            pair_y_flags = coord[3]
            age = pars[-2]
            AGN = pars[-1][6]
            SC_WHAN = pars[-1][12]

            if AGN[-1] == '!':
                AGN = AGN[:-1]
                self.ax4.scatter(x, y, color='none', edgecolors='crimson', s=20)
                self.ax5.scatter(x, y, color='none', edgecolors='crimson', s=20)
                

            if SC_WHAN[-1] == '!':
                SC_WHAN = SC_WHAN[:-1]
                self.ax4.scatter(x, y, color='none', edgecolors='black', s=50)
                self.ax5.scatter(x, y, color='none', edgecolors='black', s=50)

            if y >= -3 and y <= 3 and x >= -3.5 and x <= 2:
                k += 1
            if len(pair_x_flags) == 0 and len(pair_y_flags) == 0:
                self.ax5.scatter(x, y, s=self.color_dict[AGN][1], color=self.color_dict[AGN][0], alpha=1, marker=self.color_dict[AGN][2])
                self.ax4.scatter(x, y, s=self.cd_WHAN[SC_WHAN][1], color=self.cd_WHAN[SC_WHAN][0], marker = '.', alpha=1)
                #if len(pair_x_flags) != 0:
                #    Main.plotting_arrows(
                #        self, self.ax4, x, y, pair_x_flags, [], AGN)
                # self.ax4.scatter(plots[i][0], plots[i][1], color=self.s_m.to_rgba(age), marker="x")
            else:
                Main.plotting_arrows_WHAN(
                    self, self.ax4, self.cd_WHAN, x, y, pair_x_flags, pair_y_flags, SC_WHAN, m_y = 1, m_x = 0.3)
                Main.plotting_arrows_WHAN(
                    self, self.ax5, self.color_dict, x, y, pair_x_flags, pair_y_flags, AGN, m_y = 1, m_x = 0.3)
        
        print(k)

        for key in self.cd_WHAN_leg.keys():
            self.ax4.scatter(-99, -99, alpha= 1, color = self.cd_WHAN[key][0], marker = self.cd_WHAN[key][2], s = self.cd_WHAN[key][1], label=key)
        self.ax4.legend()

        self.ax5.scatter(-99, -99, alpha= 1, color = 'midnightblue', label='AGN', s=9)
        self.ax5.scatter(-99, -99, alpha= 1, color = 'springgreen', label='UNC', s=9)
        self.ax5.scatter(-99, -99, alpha= 1, color = 'mediumvioletred', label='SF', s=9)
        self.ax4.scatter(-99, -99, color='none', edgecolors='crimson', s=20, label='Abs. lines BPT')
        self.ax4.scatter(-99, -99, color='none', edgecolors='black', s=20, label='Abs. lines WHAN')
        self.ax5.scatter(-99, -99, color='none', edgecolors='crimson', s=20, label='Abs. lines BPT')
        self.ax5.scatter(-99, -99, color='none', edgecolors='black', s=20, label='Abs. lines WHAN')

        self.ax4.legend(loc=3)
        self.ax5.legend(loc=3)
        self.fig.savefig('./FIGURES/WHAN.pdf')

        # plt.show()


    # def plotting_WHAN_age(self):
        # 
        # self.fig2 = plt.figure(figsize=(8,8))
        # # self.ax4 = self.fig2.add_subplot()
        # self.ax4.set_xlabel('log[N[II]/H_A]')
        # self.ax4.set_ylabel('log(HA_EW)')

        # norm = mpl.colors.Normalize(vmin=8.8,vmax=10.25)
        # # choose a colormap
        # c_m = mpl.cm.jet
        # # create a ScalarMappable and initialize a data structure
        # self.s_m = mpl.cm.ScalarMappable(cmap=c_m, norm=norm)
        # self.s_m.set_array([])

        # # norm = mpl.colors.Normalize(vmin=8.8,vmax=10.25)
        # # choose a colormap
        # # c_m = mpl.cm.jet
        # # create a ScalarMappable and initialize a data structure
        # # self.s_m = mpl.cm.ScalarMappable(cmap=c_m, norm=norm)
        # # self.s_m.set_array([])

        # for pars in self.flux_er_mod9:
        #     coord = pars[-3][11]
        #     x = coord[0]
        #     pair_x_flags = coord[1]
        #     y = coord[2]
        #     age = pars[-2]
        #     mod = pars[-1]
        #     AGN = pars[-3][6]
        #     if y != -99999.0:
        #         self.ax4.scatter(
        #             x, y, s=6, color=self.s_m.to_rgba(age), alpha=1, marker=self.color_dict[AGN][2])
        #         #if len(pair_x_flags) != 0:
        #         #    Main.plotting_arrows(
        #         #        self, self.ax4, x, y, pair_x_flags, [], AGN)
        #         if mod > 0:
        #             self.ax4.text(x, y, pars[-1])
        #         # self.ax4.scatter(plots[i][0], plots[i][1], color=self.s_m.to_rgba(age), marker="x")
        #     # if mod == 0:
        #         # i = 3
        #         # self.ax4.scatter(plots[i][0], plots[i][1], color=self.s_m.to_rgba(age), marker="x")

        # # self.fig2.subplots_adjust(right=0.85)
        # # cbar_ax = self.fig2.add_axes([0.85, 0.15, 0.05, 0.7])
        # # self.fig2.colorbar(self.s_m, cax=cbar_ax)

        # self.fig2.subplots_adjust(right=0.85)
        # cbar_ax = self.fig2.add_axes([0.85, 0.15, 0.05, 0.7])
        # self.fig2.colorbar(self.s_m, cax=cbar_ax) 

        # self.ax4.axhline(y = 0.47712, color = 'black', linestyle='dotted')
        # self.ax4.axhline(y = -0.301, color = 'black', linestyle='dashed')
        # self.ax4.text(0, 0, 'Retired galaxies')

        # X_wAGN = np.arange(-0.4, 1, 0.01)
        # self.ax4.plot(X_wAGN, 0.77815125+X_wAGN*0, 'black')
        # self.ax4.text(-1.5, 1, 'SF')
        # self.ax4.text(1, 0.5, 'wAGN')
        # self.ax4.text(1, 2.5, 'sAGN')

        # Y_sAGN = np.arange(0.47712, 3, 0.01)
        # self.ax4.plot(-0.4+Y_sAGN*0, Y_sAGN, 'black')

        # plt.show()


if __name__ == '__main__':
    obj = Main('E:\LICENSE\ProgsData\main\Oleg_GAMA_belowMS.csv', 'E:\LICENSE\ProgsData\main\DirectSummationv05', 'GAMA_ETG_OLA.csv')
    obj.ola_reading()
    obj.samples_get()
    obj.gama_reading()
    obj.sorting()
    #obj.plotting_BPT()
    #obj.plotting_WHAN()