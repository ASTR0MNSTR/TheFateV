import csv
import matplotlib.pyplot as plt
import math
import numpy as np
import matplotlib as mpl
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

from __algo__ import *
from __legpars__ import *
from __stats__ import *
from __plt__ import *

# version 0.5 importing dust-correction

class hp:
    def log(figure):
        return math.log(figure, 10)

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
        
        return math.log(HA_ew, 10), HA_ew_err/(HA_ew*math.log(10)), pair_HA, HA_ew_or

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
        # self.SFR_HA = []
        # self.SFR_HA_er = []
        self.coef = []
        # self.met = []
        self.col1 = []
        self.col2 = []
        self.AGE = []
        self.AGN = []

        self.HA_l = []
        self.HB_l = []
        # self.OI_l = []
        self.OIII_l = []
        self.NII_l = []
        # self.SII_l = []
        self.HA_EW_l = []
        self.OII_l = []

        self.HA_l_er = []
        self.HB_l_er = []
        # self.OI_l_er = []
        self.OIII_l_er = []
        self.NII_l_er = []
        # self.SII_l_er = []
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
        
        # self.HdA = []
        # self.HdA_er = []
        # self.HdF = []
        # self.HdF_er = []
        # self.HgA = []
        # self.HgA_er = []
        # self.HgF = []
        # self.HgF_er = []

        self.count_1 = 0
        self.count_2 = 0
        self.count_3 = 0
        self.count_4 = 0
        
        # self.SN_HdF = []
        # self.SN_HgF = []

        # self.HgF_cont = []
        # self.HdF_cont = []
        # self.HA_cont = []
        self.GAMAIDs = []
        self.lines_prob = []
        self.ages = []
        
        self.LAGNs = []
        self.LAGNs_er = []
        self.k_OIII = []
        self.logLAGN = []
        self.k_HA = []
        self.E_B_V = []
        
        self.RA_m = []
        self.DEC_m = []
        self.OIII_m = []
        self.OIII_er_m = []
        self.LAGN_m = []
        self.LAGN_max_m = []
        self.out_m = []
        self.out_m_off = []
        
        self.OIII_m_out = []
        self.OIII_er_m_out = []
        self.LAGN_m_out = []
        self.LAGN_max_m_out = []
        self.out_m_out = []
        self.out_m_off_out = []

        # self.radio_sources = [39495, 272962, 375530, 959586, 210108, 137037, 419441, 31076, 250557, 228288, 238593, 238593, 373280, 622622, 185507]

        size = 3
        self.color_dict = color_dict_BPT

        self.cd_WHAN = cd_WHAN

        self.cd_WHAN_leg = cd_WHAN_leg


    def ola_reading(self):
        with open(self.ola_file, 'r') as input:
            reader = csv.DictReader(input)
            for row in reader:
                self.sample_dict.update({int(row['CATAID_1']): [int(row['below=0/MS=1']), float(row['ager_percentile50'])]})

    def samples_get(self):
        with open(r'E:\backup\backup_BPT\Sep2023\outflow4OR_1.txt', 'r') as f:
            header = f.readline()
            print(header)
            lines = f.readlines()
            lines_strip = [item.strip() for item in lines]
            for line in lines_strip:
                self.RA_m.append(float(line.split()[0]))
                self.DEC_m.append(float(line.split()[1]))
                self.OIII_m.append(float(line.split()[2]))
                self.OIII_er_m.append(float(line.split()[3]))
                self.LAGN_m.append(float(line.split()[4]))
                self.LAGN_max_m.append(float(line.split()[5]))
                self.out_m_off.append(float(line.split()[6]))
                try:
                    self.out_m.append(float(line.split()[7]))
                except:
                    self.out_m.append(-99)
        print(self.RA_m)            

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
        RA = float(galaxy_pars[2])
        DEC = float(galaxy_pars[3])
        bms = -1
        age = -1
        if CATAID in self.sample_dict.keys():
            bms, age = self.sample_dict[CATAID]
        else:
            return 1
        
        # k = 0
        OIII = float(galaxy_pars[167])
        
        # for i in range(len(self.RA_m)):
        #     if abs(self.RA_m[i] - RA) <= 0.01 and abs(self.DEC_m[i] - DEC) <= 0.01 and abs(self.OIII_m[i] - OIII) <= 0.01:
        #         k = 1
        #         print('Gotcha')
        #         OIII_m = self.OIII_m[i]
        #         OIII_er_m = self.OIII_er_m[i]
        #         LAGN_m = self.LAGN_m[i]
        #         LAGN_max = self.LAGN_max_m[i]
        #         out_m = self.out_m[i]
        #         out_m_off = self.out_m_off[i]
        #         break
        
        # if k == 0:
        #     return 1
                
        Z = float(galaxy_pars[4])

        SPEC_ID = galaxy_pars[0]
        HA = float(galaxy_pars[82])
        HA_er = float(galaxy_pars[81])
        HB = float(galaxy_pars[182])
        HB_er = float(galaxy_pars[181])
        OIII = float(galaxy_pars[167])
        OIII_er = float(galaxy_pars[166])  # OIIIR, 5007
        OIIIB = float(galaxy_pars[172])
        OIIIB_er = float(galaxy_pars[171])
        NII = float(galaxy_pars[72])
        NII_er = float(galaxy_pars[71])  # NIIR, 6583

        HA_cont = float(galaxy_pars[78])
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

        # dust_correction module:
        
        AGN, X, pair_x_flags, Y, pair_y_flags, SC_WHAN, LAGN, LAGN_er, HA_ew, HA_ew_err, pair_HA = AGN_reg(OIII, OIII_er, HB, HB_er, NII, NII_er, HA, HA_er, HA_EW, HA_EW_ERR, Z)

        
        X_er = 0
        Y_er = 0
        try:
            log_LAGN = math.log10(LAGN*(10**43))
        except ValueError:
            log_LAGN = -99999.0

        SFR_HA, SFR_HA_er = hp.SFR(HA, HA_er, Z)

        SURV = galaxy_pars[6]
        SURV_CODE = galaxy_pars[7]
        IS_BEST = galaxy_pars[8]
        IS_SBEST = galaxy_pars[9]

        FLUXES = [HA, HB, OIII, NII, SII, OI, OII]
        FLUXES_ER = [HA_er, HB_er, OIII_er, NII_er, SII_er, OI_er, OII_er]

        res_out = [X, X_er, pair_x_flags, Y, Y_er, pair_y_flags]

        kwargs = [SURV, SURV_CODE, IS_BEST, IS_SBEST, CATAID,
                  [SFR_HA, SFR_HA_er], AGN, FLUXES, FLUXES_ER, bms, [RA, DEC, Z, SPEC_ID], [X, pair_x_flags, HA_ew, pair_HA], SC_WHAN,
                  [HdA, HdA_er, HdF, HdF_er, HgA, HgA_er, HgF, HgF_er], [HA_EW, HA_EW_ERR], HgF_cont, HdF_cont, HA_cont, age, LAGN, LAGN_er]
                #   [HdA, HdA_er, HdF, HdF_er, HgA, HgA_er, HgF, HgF_er], [HA_EW_1, HA_EW_ERR_1], HgF_cont, HdF_cont, HA_cont, age, LAGN, LAGN_er, k_OIII, log_LAGN, k_HA, E_B_V, OIII_m, OIII_er_m, LAGN_m, out_m, LAGN_max, out_m_off]
        res_out_out = [res_out, kwargs]
        return res_out_out

    def sorting(self):
        counter_all = 0
        counter_best = 0
        count_OIII = 0
        for pars in self.flux_er_mod:
            if type(pars) == list:
                counter_all += 1
                if pars[-1][2] == 'true' and pars[-1][3] == 'true' and (pars[-1][0] == 'GAMA' or pars[-1][0] == 'SDSS') and pars[-1][6] != 'NDA': #!!!!
                # if pars[-1][2] == 'true' and pars[-1][3] == 'true': #!!!!
                    # if pars[-1][4] in self.radio_sources:
                    #     print(pars[-1][4], pars[-1][10][3], pars[-1][6], pars[-1][12])
                    Main.exporting(self, pars)
                    counter_best += 1
                    self.flux_er_mod9.append(pars)
                    if pars[-1][7][2] >= 3*pars[-1][8][2]:
                        count_OIII += 1
                else:
                    pass

        print('All processed galaxies: ', counter_all)
        print('Galaxies with proper data: ', counter_best)
        print('Galaxies with 3sigma OIII: ', count_OIII)
        Main.file_out(self)

    def exporting(self, pars):
        self.SURV.append(pars[-1][0])
        self.SURV_CODE.append(pars[-1][1])
        self.IS_BEST.append(pars[-1][2])
        self.IS_SBEST.append(pars[-1][3])
        self.GAMAID.append(pars[-1][4])
        # self.SFR_HA.append(pars[-1][5][0])
        # self.SFR_HA_er.append(pars[-1][5][1])

        # changed line
        self.AGN.append(pars[-1][6])

        self.HA_l.append(pars[-1][7][0])
        self.HB_l.append(pars[-1][7][1])
        self.OIII_l.append(pars[-1][7][2])
        self.NII_l.append(pars[-1][7][3])
        # self.SII_l.append(pars[-1][7][4])
        # self.OI_l.append(pars[-1][7][5])
        self.OII_l.append(pars[-1][7][6])

        self.HA_l_er.append(pars[-1][8][0])
        self.HB_l_er.append(pars[-1][8][1])
        self.OIII_l_er.append(pars[-1][8][2])
        self.NII_l_er.append(pars[-1][8][3])
        # self.SII_l_er.append(pars[-1][8][4])
        # self.OI_l_er.append(pars[-1][8][5])
        self.OII_l_er.append(pars[-1][8][6])

        self.BMS.append(pars[-1][9])
        self.RA.append(pars[-1][10][0])
        self.DEC.append(pars[-1][10][1])
        self.Z.append(pars[-1][10][2])
        self.SPEC_ID.append(pars[-1][10][3])
        self.SC_WHAN.append(pars[-1][12])

        # self.HdA.append(pars[-1][13][0])
        # self.HdA_er.append(pars[-1][13][1])
        # self.HdF.append(pars[-1][13][2])
        # self.HdF_er.append(pars[-1][13][3])
        # self.HgA.append(pars[-1][13][4])
        # self.HgA_er.append(pars[-1][13][5])
        # self.HgF.append(pars[-1][13][6])
        # self.HgF_er.append(pars[-1][13][7])

        # self.SN_HdF.append(pars[-1][13][2]/pars[-1][13][3])
        # self.SN_HgF.append(pars[-1][13][6]/pars[-1][13][7])

        self.HA_EW_l.append(pars[-1][14][0])
        self.HA_EW_l_er.append(pars[-1][14][1])

        # self.HgF_cont.append(pars[-1][15])
        # self.HdF_cont.append(pars[-1][16])
        # self.HA_cont.append(pars[-1][17])
        self.ages.append(pars[-1][18])
        self.LAGNs.append(pars[-1][19])
        self.LAGNs_er.append(pars[-1][20])
        # self.OIII_m_out.append(pars[-1][25])
        # self.OIII_er_m_out.append(pars[-1][26])
        # self.LAGN_m_out.append(pars[-1][27])
        # self.out_m_out.append(pars[-1][28])
        # self.LAGN_max_m_out.append(pars[-1][29])
        # self.out_m_off_out.append(pars[-1][30])

    def file_out(self):

        Dict = {
            'SPECID' : self.SPEC_ID,
            'CATAID_1': self.GAMAID,
            'SURVEY': self.SURV,
            'IS_BEST': self.IS_BEST,
            'IS_SBEST': self.IS_SBEST,
            'RA' : self.RA,
            'DEC' : self.DEC,
            'Z' : self.Z,
            'BPT': self.AGN,
            'WHAN' : self.SC_WHAN,
            # 'SN_HgF' : self.SN_HgF,
            # 'SN_HdF' : self.SN_HdF,
            # 'HgF' : self.HgF,
            # 'HgF_cont' : self.HgF_cont,
            # 'HgF_er' : self.HgF_er,
            # 'HdF' : self.HdF,
            # 'HdF_cont' : self.HdF_cont,
            # 'HdF_er' : self.HdF_er,
            'HA_EW' : self.HA_EW_l,
            'HA_EW_ERR' : self.HA_EW_l_er,
            'HA': self.HA_l,
            # 'HA_cont' : self.HA_cont,
            'HA_er': self.HA_l_er,
            'HB': self.HB_l,
            'HB_er': self.HB_l_er,
            'OIII': self.OIII_l,
            # 'OIII_MM' :self.OIII_m_out,
            'OIII_er': self.OIII_l_er,
            # 'OIII_er_MM' : self.OIII_er_m_out,
            'NII': self.NII_l,
            'NII_er': self.NII_l_er,
            # 'HdA' : self.HdA,
            # 'HdA_er' : self.HdA_er,
            # 'HgA' : self.HgA,
            # 'HgA_er' : self.HgA_er, 
            # 'SII' : self.SII_l,
            # 'SII_er' : self.SII_l_er,
            # 'OI' : self.OI_l,
            # 'OI_er' : self.OI_l_er,
            # 'OII' : self.OII_l,
            # 'OII_er' : self.OII_l_er,
            # 'log(LAGN)' : self.logLAGN,
            'LAGN' : self.LAGNs,
            # 'LAGN_MM_max' : self.LAGN_max_m_out,
            'LAGN_er' : self.LAGNs_er,
            'BMS': self.BMS,
            # 'LAGN_MM' : self.LAGN_m_out,
            # 'OUT_MM' : self.out_m_out,
            # 'OUT_MM_mol' : self.out_m_off_out
        }

        df = pd.DataFrame(Dict)
        df.to_csv(self.filename_out, index=False)

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
            X_11 = np.arange(-0.45, 1.5, 0.01)
            ax.plot(X_1, (0.61/(X_1 - 0.47)) + 1.19,
                      c='k')  # Kewley, 2001
            ax.plot(X_111, (0.61/(X_111 - 0.05)) + 1.3,
                      c='k', linestyle='dashed')  # Kauffman, 2003
            # https://adsabs.harvard.edu/full/2003MNRAS.346.1055K
            ax.plot(X_11, 1.01*X_11 + 0.48, c='r', linestyle='dotted')
            ax.text(-1.5, 1.2, 'AGN')
            ax.text(0, -1, 'UNC', ha='center', va='center')
            ax.text(-1.5, -0.5, 'SF')
            ax.text(0.5, -0.5, 'LINER')
            ax.set_box_aspect(1)

        self.ax4.set_ylabel(r'$\log{\mathrm{([OIII]/H\beta)}}$')

        norm = mpl.colors.Normalize(vmin=8.8,vmax=10.0)
        c_m = mpl.cm.jet
        self.s_m = mpl.cm.ScalarMappable(cmap=c_m, norm=norm)
        self.s_m.set_array([])
        
        cmap = plt.cm.jet  
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

        for pars in self.flux_er_mod9:
            plots = pars[0]
            
            AGN = pars[-1][6]
            SC_WHAN = pars[-1][12]
            age = pars[-1][18]

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
            if x != -100 and x != -99 and y != -99 and y != 100:
                X.append(x)
                Y.append(y)
                AGE.append(age)
                AGN_flags.append(AGN)
                if len(pair_x_flags) == 0 and len(pair_y_flags) == 0:
                    try:
                        self.ax4.scatter(
                            x, y, s=1.5, color=self.color_dict[AGN][0], alpha=1)
                    #self.ax4.scatter(x, y, s=1.5, color=self.s_m.to_rgba(age), alpha=1)
                        self.ax5.scatter(x, y, s=self.cd_WHAN[SC_WHAN][1], color=self.cd_WHAN[SC_WHAN][0], marker =self.cd_WHAN[SC_WHAN][2], alpha=1)
                        #self.ax4.scatter(x, y, s=self.cd_WHAN[SC_WHAN][1], color=self.cd_WHAN[SC_WHAN][0], marker =self.cd_WHAN[SC_WHAN][2], alpha=0.5)
                        #self.ax5.scatter(x, y, s=self.cd_WHAN[SC_WHAN][1], color=self.s_m.to_rgba(age), marker =self.cd_WHAN[SC_WHAN][2], alpha=0.5)
                        self.ax_med_BPT.scatter(x, y, s=1.5, color=self.s_m.to_rgba(age), alpha=0.3)
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
            self.ax_med_BPT.scatter(big_point_X, big_point_Y, s = 150, color = self.s_m_INT.to_rgba(big_point_age), marker=item[3][1])

        
        for key in BPT_color_plt:
            self.ax_med_BPT.scatter(-99, -99, alpha=1, color = BPT_color_plt[key][0], s = BPT_color_plt[key][1], marker= BPT_color_plt[key][2], label=key)
        
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
            ax.axhline(y = 0.47712, color = 'black', linestyle='dashed')
            ax.axhline(y = -0.301, color = 'black', linestyle='dotted')
            ax.text(0.8, 0, 'ELR')
            ax.text(0.8, -1.5, 'LLR')

            X_wAGN = np.arange(-0.4, 2.5, 0.01)
            ax.plot(X_wAGN, 0.77815125+X_wAGN*0, 'black')
            ax.text(-1.5, 2, 'SF')
            ax.text(0.6, 0.5, 'wAGN')
            ax.text(0.6, 2, 'sAGN')

            Y_sAGN = np.arange(0.47712, 3, 0.01)
            ax.plot(-0.4+Y_sAGN*0, Y_sAGN, 'black')
            ax.set_xlabel(r"$\log \mathrm{([NII]/H\alpha)}$")
            ax.set_box_aspect(1)
        

        self.ax6.set_ylabel(r"$\log \mathrm{(EW_{\mathrm{H\alpha}})}$")
        
        k = 0
        X = []
        Y = []
        AGE = []
        AGN_flags = []

        for pars in self.flux_er_mod9:
            age = pars[-1][18]
            coord = pars[-1][11]
            x = coord[0]
            pair_x_flags = coord[1]
            y = coord[2]
            pair_y_flags = coord[3]
            AGN = pars[-1][6]
            SC_WHAN = pars[-1][12]

            if AGN[-1] == '!':
                AGN = AGN[:-1]
                self.ax7.scatter(x, y, color='none', edgecolors='crimson', s=20)
                self.ax6.scatter(x, y, color='none', edgecolors='crimson', s=20)
                

            if SC_WHAN[-1] == '!':
                SC_WHAN = SC_WHAN[:-1]
                self.ax7.scatter(x, y, color='none', edgecolors='black', s=50)
                self.ax6.scatter(x, y, color='none', edgecolors='black', s=50)

            if y >= -3 and y <= 3 and x >= -3.5 and x <= 2:
                k += 1
                
            if x != -100 and x != -99 and y != -99 and y != 100:
                X.append(x)
                Y.append(y)
                AGE.append(age)
                AGN_flags.append(SC_WHAN)
                if len(pair_x_flags) == 0 and len(pair_y_flags) == 0:
                    self.ax6.scatter(x, y, s=self.color_dict[AGN][1], color=self.color_dict[AGN][0], alpha=1, marker=self.color_dict[AGN][2])
                    self.ax_med_WHAN.scatter(x, y, s=self.color_dict[AGN][1], color=self.s_m.to_rgba(age), alpha=0.3, marker=self.color_dict[AGN][2])
                    self.ax7.scatter(x, y, s=self.cd_WHAN[SC_WHAN][1], color=self.cd_WHAN[SC_WHAN][0], marker = '.', alpha=1)
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
            self.ax_med_WHAN.scatter(big_point_X, big_point_Y, s = 150, color = self.s_m.to_rgba(big_point_age), marker=item[3][1])
        
        for key in WHAN_color_plt:
            self.ax_med_WHAN.scatter(-99, -99, alpha=1, color = WHAN_color_plt[key][0], s = WHAN_color_plt[key][1], marker= WHAN_color_plt[key][2], label=key)

        for key in self.cd_WHAN_leg.keys():
            self.ax7.scatter(-99, -99, alpha= 1, color = self.cd_WHAN_leg[key][0], marker = self.cd_WHAN_leg[key][2], s = self.cd_WHAN_leg[key][1], label=key)

        self.ax6.scatter(-99, -99, alpha= 1, color = 'midnightblue', label='AGNXY', s = 30, marker='o')
        self.ax6.scatter(-99, -99, alpha= 1, color = 'dodgerblue', label='AGNX', s = 30, marker='o')
        self.ax6.scatter(-99, -99, alpha= 1, color = 'springgreen', label='UNCXY', s = 30, marker='o')
        self.ax6.scatter(-99, -99, alpha= 1, color = 'darkgreen', label='UNCX', s = 30, marker='o')
        self.ax6.scatter(-99, -99, alpha= 1, color = 'mediumvioletred', label='SFGXY', s = 30, marker='o')
        self.ax6.scatter(-99, -99, alpha= 1, color = 'crimson', label='SFGX', s = 30, marker='o')
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
        
        self.fig.savefig('./FIGURES_IN_PAPER/BPT_WHAN.pdf', dpi=70, transparent = True, bbox_inches = 'tight', pad_inches = 0.0001)
        #self.fig.savefig('WHAN.pdf')

        # plt.show()
        
    def plotting_EW_age(self):
        
        self.gs_top = plt.GridSpec(1, 1, wspace=0)
        #self.fig = plt.figure(figsize=(12, 12), tight_layout=True)
        self.fig = plt.figure(figsize=(12, 12))

        self.ax = self.fig.add_subplot(self.gs_top[0,0])
        
        self.ax.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, right=True, labelleft=True, left=True, direction='in')
        
        self.ax.set_xlabel(r"$log(age/yr)$")
        self.ax.set_xlim([8.7, 10.1])
        self.ax.set_ylim([-3, 2.5])
        self.ax.set_xticks(np.arange(8.8, 10.1, 0.2))
        self.ax.set_yticks(np.arange(-3, 2.6, 0.5))
        
        self.ax.set_ylabel(r"$log(EW_{H\alpha})$")
        
        for pars in self.flux_er_mod9:
            age = pars[-1][18]
            coord = pars[-1][11]
            x = age
            pair_x_flags = []
            y = coord[2]
            pair_y_flags = coord[3]
            AGN = pars[-1][6]
            SC_WHAN = pars[-1][12]
            
            if len(pair_y_flags) == 0:
                    #self.ax6.scatter(x, y, s=self.color_dict[AGN][1], color=self.s_m.to_rgba(age), alpha=1, marker=self.color_dict[AGN][2])
                    #self.ax7.scatter(x, y, s=self.cd_WHAN[SC_WHAN][1], color=self.cd_WHAN[SC_WHAN][0], marker = '.', alpha=1)
                self.ax.scatter(x, y, s=self.cd_WHAN[SC_WHAN][1], color=self.cd_WHAN[SC_WHAN][0], marker = '.', alpha=1)
                #self.ax.scatter(x, y, s=self.cd_WHAN[SC_WHAN][1], color=self.s_m.to_rgba(age), marker = '.', alpha=0.5)
            else:
                Main.plotting_arrows(self, self.ax, x, y, pair_x_flags, pair_y_flags, self.cd_WHAN[SC_WHAN][0], m_y = 1, m_x = 0.8)
        
        self.ax.set_box_aspect(1)
        self.fig.savefig('./FIGURES/EW_age.pdf')


if __name__ == '__main__':
    obj = Main('E:\LICENSE\ProgsData\main\Oleg_GAMA_belowMS.csv', 'E:\LICENSE\ProgsData\main\DirectSummationv05', r'E:\backup\backup_BPT\GAMA_ETG_OLA.csv')
    obj.ola_reading()
    # obj.samples_get()
    obj.gama_reading()
    obj.sorting()
    obj.plotting_BPT()
    obj.plotting_WHAN()
    #obj.plotting_EW_age()