import csv
import matplotlib.pyplot as plt
import math
import numpy as np
import matplotlib as mpl
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

#version 0.5 importing dust-correction

class hp:
    def log(figure):
        return math.log(figure, 10)
    
    def log_er(item):
        return [[abs(math.log(1-item, 10))], [abs(math.log(1+item, 10))]]
    
    def quality_coeff(pair):
        pars = pair[0]
        coef = 0
        for item in pars:
            if item != 0 and len(item) == 4:
                coef+=10
            elif item != 0 and len(item) == 3:
                coef+=1
        return coef
    
    def atten_coef(A_V, mode):            
        R_V = 4.05
        OIII_c = 2.659*(-2.156 + (1.509/0.5007) - (0.198/(0.5007**2)) + (0.011/(0.5007**3))) + R_V
        NII_c = 2.659*(-1.857 + (1.040/0.6583)) + R_V
        SII_c = 2.659*(-1.857 + (1.040/0.6732)) + R_V
        OI_c = 2.659*(-1.857 + (1.040/0.6300)) + R_V
        HB_c = 2.659*(-2.156 + (1.509/0.4861) - (0.198/(0.4861**2)) + (0.011/(0.4861**3))) + R_V
        res = [OIII_c, NII_c, SII_c, OI_c, HB_c]
        res_out = [item*A_V/R_V for item in res]
        res_out_out = [10**(0.4*item) for item in res_out]
        if mode > 0:
            print(f'{mode} {res_out_out[0]} {res_out_out[1]} {res_out_out[2]} {res_out_out[3]} {res_out_out[4]} {10**(0.4*A_V)}\n')
        return res_out_out
    
    def SFR(HA, HA_er, z):
        if HA == -99999.0 or HA_er < 0: 
            return -99.9, -99.9
        
        cosmo = FlatLambdaCDM(H0=75 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)
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
    
    def indexing(pair):
        k = 0
        n = 0
        for i in range(len(pair)):
            if pair[i][0] == -99999.0 or pair[i][1] < 0: #????????
                pair.append(-1)
                return pair
        
            if abs(pair[i][0]) < 2*pair[i][1]: 
                k += 1
                n = i
            elif pair[i][0] < 0:
                pair.append(-1)
                return pair
        
        if k == 2:
            pair.append(-1)
            return pair
        elif k == 1:
            pair[n][0] = pair[n][0] + 2*pair[n][1]
            pair.append(n)
            return pair
        else:
            pair.append(10)
            return pair
    
    def val_mist(quadro): #x-pair, y-pair
        for item in quadro:
            if item[2] == -1:
                return 0
        
        x = hp.log(quadro[0][0][0]/quadro[0][1][0])        
        y = hp.log(quadro[1][0][0]/quadro[1][1][0]) 

        sym = []
        if quadro[0][2] == 0:
            sym.append(0)
        if quadro[0][2] == 1:
            sym.append(1)
        if quadro[1][2] == 0:
            sym.append(2)
        if quadro[1][2] == 1:
            sym.append(3)
        if len(sym) == 1:
            sym.append(4)
        if len(sym) == 0:
            sym.append(4)
            sym.append(4)
        
        try:
            if sym[0] == 4 and sym[1] == 4:
                x_err = hp.log_er((quadro[0][0][1]/quadro[0][0][0]) + (quadro[0][1][1]/quadro[0][1][0]))
                y_err = hp.log_er((quadro[1][0][1]/quadro[1][0][0]) + (quadro[1][1][1]/quadro[1][1][0]))
                return [x, y, x_err, y_err]
            else:
                return [x, y, sym]
        except:
            return [x, y, sym]
    
#    def metallicity(fluxes):
#        if fluxes == 0:
#            return -99.9, -99.9, -99.9
#        else:
#            met = fluxes[0] + 0.264*fluxes[1] - 3.23
#            if len(fluxes) == 4:
#                err_plus = fluxes[2][1][0] + 0.264*fluxes[3][1][0]
#                err_minus = fluxes[2][0][0] + 0.264*fluxes[3][0][0]
#                return met, err_plus, err_minus
#            else:
#                signs = []
#                for item in fluxes[2]:
#                    if item == 1 or item == 3:
#                        signs.append('+')
#                    elif item == 2 or item == 0:
#                        signs.append('-')
#                    else:
#                        signs.append('0')
#                return met, signs[0], signs[1]

    def metallicity(NII, NII_er, HA, HA_er):
        if (HA != -99999.0 and HA > 2*HA_er) and (NII != -99999.0 and NII > 2*NII_er):
            met = 0.73*math.log(NII/HA, 10) - 2.88 #log(O/H)
            err_plus = 0.73*hp.log_er((NII_er/NII) + (HA_er/HA))[1][0]
            err_minus = 0.73*hp.log_er((NII_er/NII) + (HA_er/HA))[0][0]
            return met, err_plus, err_minus
        elif (HA != -99999.0 and HA + 2*HA_er > 0 and HA < 2*HA_er) and (NII != -99999.0 and NII > 2*NII_er):
            met = 0.73*math.log(NII/(HA+2*HA_er), 10) - 2.88 #log(O/H)
            signs = ['+', 0]
            return met, signs[0], signs[1]
        elif (NII != -99999.0 and NII + 2*NII_er > 0 and NII < 2*NII_er) and (HA != -99999.0 and HA > 2*HA_er):
            met = 0.73*math.log((NII+2*NII_er)/HA, 10) - 2.88 #log(O/H)
            signs = ['-', 0]
            return met, signs[0], signs[1]
        else:
            return -99.9, -99.9, -99.9

    
    def func_anal(res_out, i, AGN_count, UNC_count, SF_count, lim_divup, lim_divdown, lim_sum):
        if len(res_out[i]) == 3:
            if res_out[i][0] > lim_divdown and res_out[i][1] < lim_sum:
                if res_out[i][2][0] != 0 and res_out[i][2][1] != 0:
                    AGN_count += 1
                else:
                    UNC_count += 1
            elif res_out[i][0] > lim_divdown and res_out[i][1] > lim_sum:
                if (res_out[i][2][0] == 0 and res_out[i][2][1] == 2) or (res_out[i][2][0] == 2 and res_out[i][2][1] == 0):
                    UNC_count += 1
                else:
                    AGN_count += 1
                            
            elif res_out[i][1] > lim_sum and res_out[i][0] < lim_divdown:
                if res_out[i][2][0] != 2 and res_out[i][2][1] != 2:
                    AGN_count += 1
                else:
                    UNC_count += 1

            else:
                if res_out[i][1] > (lim_divup/(res_out[i][0] - lim_divdown)) + lim_sum:
                    if (res_out[i][2][0] == 1 and res_out[i][2][1] == 4) or (res_out[i][2][0] == 4 and res_out[i][2][1] == 1) or (res_out[i][2][0] == 3 and res_out[i][2][1] == 4) or (
                    res_out[i][2][0] == 4 and res_out[i][2][1] == 3) or (res_out[i][2][0] == 3 and res_out[i][2][1] == 1) or (res_out[i][2][0] == 1 and res_out[i][2][1] == 3):
                        AGN_count += 1
                    else:
                        UNC_count += 1
                        
                elif res_out[i][1] < (lim_divup/(res_out[i][0] - lim_divdown)) + lim_sum:
                    if (res_out[i][2][0] == 0 and res_out[i][2][1] == 4) or (res_out[i][2][0] == 4 and res_out[i][2][1] == 0) or (res_out[i][2][0] == 2 and res_out[i][2][1] == 4) or (
                    res_out[i][2][0] == 4 and res_out[i][2][1] == 2) or (res_out[i][2][0] == 0 and res_out[i][2][1] == 2) or (res_out[i][2][0] == 2 and res_out[i][2][1] == 0):
                        SF_count += 1
                    else:
                        UNC_count += 1        
                else:
                    print('Wierd!')
                    
        elif len(res_out[i]) == 4:
            if res_out[i][1] > (lim_divup/(res_out[i][0] - lim_divdown)) + lim_sum: 
                if res_out[i][1] - abs(res_out[i][3][0][0]) > (lim_divup/(res_out[i][0] - lim_divdown)) + lim_sum and res_out[i][1] > (lim_divup/(res_out[i][0] - abs(res_out[i][2][0][0]) - lim_divdown)) + lim_sum:
                    AGN_count += 1
                else:
                    UNC_count += 1
                        
            elif res_out[i][1] < (lim_divup/(res_out[i][0] - lim_divdown)) + lim_sum: 
                if res_out[i][1] + abs(res_out[i][3][1][0]) > (lim_divup/(res_out[i][0] - lim_divdown)) + lim_sum or res_out[i][0] + abs(res_out[i][2][1][0]) > lim_divdown or res_out[i][1] > (lim_divup/(res_out[i][0] + abs(res_out[i][2][1][0]) - lim_divdown)) + lim_sum:
                    UNC_count += 1
                else:
                    SF_count += 1
            else:
                print('DONT FORGET ME')
        else:
            pass
                        #check the moment with points in extreme area
        return AGN_count, UNC_count, SF_count
    
    def AGN_registration(res_out):
        AGN_count = 0
        UNC_count = 0
        SF_count = 0
        total_count = 0
        if res_out[0] != 0:
            total_count += 1
            AGN_count, UNC_count, SF_count = hp.func_anal(res_out, 0, AGN_count, UNC_count, SF_count, 0.61, 0.47, 1.19)
                #elif i == 1:
                #    AGN_count, UNC_count, SF_count = hp.func_anal(res_out, i, AGN_count, UNC_count, SF_count, 0.72, 0.32, 1.3)
                #elif i == 2: 
                #    AGN_count, UNC_count, SF_count = hp.func_anal(res_out, i, AGN_count, UNC_count, SF_count, 0.73, -0.59, 1.33)
        
        if (total_count == 0):
            AGN = 'NOEL'
        elif SF_count == total_count:
            AGN = 'NO'
        elif AGN_count == total_count:
            AGN = 'YES'
        else:
            AGN = 'UNC'
        return AGN
    
class Main(hp):
    def __init__(self, ola_file, gama_file, sample_file, num):
        self.gama_file = gama_file
        self.sample_file = sample_file
        self.ola_file = ola_file
        self.num = num

        self.gama_dict = []
        self.sample_dict = {}
        self.headers = [] #headers for spectra_file

        self.flux_er_mod = []

        self.flux_er_mod9 = []
        self.flux_er_mod52 = []

        self.sample9 = []
        self.sample52 = []

        #EXPORTING

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
        self.HB_l= []
        self.OI_l = []
        self.OIII_l = []
        self.NII_l = []
        self.SII_l = []
        self.HA_EW_l = []

        self.HA_l_er = []
        self.HB_l_er = []
        self.OI_l_er = []
        self.OIII_l_er = []
        self.NII_l_er = []
        self.SII_l_er = []
        self.HA_EW_l_er = []

    def ola_reading(self):
        with open(self.ola_file, 'r') as input:
            lines = input.readlines()
        lines_stripped = [line.strip() for line in lines]
        for line in lines_stripped:
            GAMAID = int(line.split()[0])
            AGE = float(line.split()[105])
            self.sample_dict.update({GAMAID : [AGE, -1]})
    
    def samples_get(self):
        data_dict = []
        with open(self.sample_file, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                data_dict.append({'OI' : row['OI'], 'GAMAID' : row['GAMAID'], 'age' : row['age']})
        
        for galaxy in data_dict:
            if galaxy['OI'] == '':
                IND = 0
            else:
                IND = int(galaxy['OI'])
            self.sample_dict.update({int(galaxy['GAMAID']) : [float(galaxy['age']), IND]})

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
        mod = -2
        age = -1
        try:
            age, mod = self.sample_dict[CATAID]
        except:
            return 1
        
        Z = float(galaxy_pars[4])
        HA = float(galaxy_pars[82])
        HA_er = float(galaxy_pars[81])
        HB = float(galaxy_pars[182]) #should I include errors of HA and HB
        HB_er = float(galaxy_pars[181])
        OIII = float(galaxy_pars[172])
        OIII_er = float(galaxy_pars[171]) #OIIIB, 5007
        OI = float(galaxy_pars[97])
        OI_er = float(galaxy_pars[96]) #OIB(OI), 6300
        NII = float(galaxy_pars[77])
        NII_er = float(galaxy_pars[76]) #NIIB, 6583
        SII = float(galaxy_pars[67])
        SII_er = float(galaxy_pars[66]) #SIIB, 6731, should I change/add to SIIR        

        #dust_correction module:
        balm_coef = 2.86
        f_HA = 1
        if HA_er > 0 and HB_er > 0 and HA != -99999.0 and HB != -99999.0 and HA > 2*HA_er and HB > 2*HB_er and HA/HB > 2.86:
            f_HA = ((HA/HB)/2.86)**(2.114)
            A_V = 2.5*math.log(f_HA, 10)
            HA = HA*f_HA
            HA_er = HA_er*f_HA
            coefs = hp.atten_coef(A_V, mod)
            HB = HB*coefs[-1]
            HB_er = HB_er*coefs[-1]
            if OIII != -99999.0 and OIII_er > 0:
                OIII = OIII*coefs[0]
                OIII_er = OIII_er*coefs[0]
            if NII != -99999.0 and NII_er > 0:
                NII = NII*coefs[1]
                NII_er = NII_er*coefs[1]
            if SII != -99999.0 and SII_er > 0:
                SII = SII*coefs[2]
                SII_er = SII_er*coefs[2]
            if OI != -99999.0 and OI_er > 0:
                OI = OI*coefs[3]
                OI_er = OI_er*coefs[3]

        HB_trick = False
        if (HB != -99999.0 and HB < 0 and HB + HB_er*2 < 0) and (OIII != -99999.0 and OIII > 2*OIII_er):
            HB = 0
            HB_trick = True

        SFR_HA, SFR_HA_er = hp.SFR(HA, HA_er, Z)

        HA_EW = float(galaxy_pars[80])
        HA_EW_ERR = float(galaxy_pars[79])
        SURV = galaxy_pars[6]
        SURV_CODE = galaxy_pars[7]
        IS_BEST = galaxy_pars[8]
        IS_SBEST = galaxy_pars[9]

        FLUXES = [HA, HB, OI, OIII, NII, SII, HA_EW]
        FLUXES_ER = [HA_er, HB_er, OI_er, OIII_er, NII_er, SII_er, HA_EW_ERR]

        par_y = [[OIII, OIII_er], [HB, HB_er]]
        par_y_2 = [[HA_EW, HA_EW_ERR], [1, 0]]

        par_x_0 = [[NII, NII_er], [HA, HA_er]]
        par_x_1 = [[SII, SII_er], [HA, HA_er]]
        par_x_2 = [[OI, OI_er], [HA, HA_er]]

        par_met_0 = [[NII, NII_er], [SII, SII_er]]
        par_met_1 = [[NII, NII_er], [HA, HA_er]]

        par_mets = [hp.indexing(par) for par in [par_met_0, par_met_1]]
        met_out = [hp.val_mist([par_mets[0], par_mets[1]])]

        pars = [par_x_0, par_x_1, par_x_2, par_y_2, par_y]
        pairs = [hp.indexing(par) for par in pars]

        pair_y = pairs[-1]
        pair_y_2 = pairs[-2]

        res_out = [hp.val_mist([pair, pair_y]) for pair in pairs[0:3]]
        res_out.append(hp.val_mist([pars[0], pair_y_2]))


        #AGN_registration

        if HB_trick == True:
            if (NII != -99999.0 and NII < 0 and NII + NII_er*2 < 0) or (HA != -99999.0 and HA < 0 and HA + HA_er*2 < 0):
                if math.log(OIII/(2*HB_er), 10) > 0.8:
                    AGN = 'YES'
                else: #OIII/HB < 0.8
                    AGN = 'UNC'
            else:
                AGN = hp.AGN_registration(res_out)
        elif HB == -99999.0 or OIII == -99999.0 or HA == -99999.0 or NII == -99999.0:
            AGN = 'NDA'
        else:
            AGN = hp.AGN_registration(res_out)

        #changed line
        met_out = hp.metallicity(NII, NII_er, HA, HA_er)

        kwargs = [SURV, SURV_CODE, IS_BEST, IS_SBEST, CATAID, [SFR_HA, SFR_HA_er], AGN, FLUXES, FLUXES_ER]
        res_out_out = [res_out, met_out, kwargs, age, mod]
        return res_out_out

    def sorting(self):
        for pars in self.flux_er_mod:
            if type(pars) == list:
                if pars[-1] >= -1 and pars[-3][2] == 'true' and pars[-3][3] == 'true' and (pars[-3][0] == 'GAMA' or pars[-3][0] == 'SDSS'):
                    Main.exporting(self, pars)
                    self.flux_er_mod9.append(pars)
            else:
                pass
        Main.file_out(self)
    
    def exporting(self, pars):
        self.OI.append(pars[-1])
        self.SURV.append(pars[-3][0])
        self.SURV_CODE.append(pars[-3][1])
        self.IS_BEST.append(pars[-3][2])
        self.IS_SBEST.append(pars[-3][3])
        self.GAMAID.append(pars[-3][4])
        self.SFR_HA.append(pars[-3][5][0])
        self.SFR_HA_er.append(pars[-3][5][1])
        self.coef.append(hp.quality_coeff(pars))

        #changed line
        met, col1, col2 = pars[-4]
        self.met.append(met)
        self.col1.append(col1)
        self.col2.append(col2)
        self.AGE.append(pars[-2])
        self.AGN.append(pars[-3][6])
        self.HA_l.append(pars[-3][7][0])
        self.HB_l.append(pars[-3][7][1])
        self.OI_l.append(pars[-3][7][2])
        self.OIII_l.append(pars[-3][7][3])
        self.NII_l.append(pars[-3][7][4])
        self.SII_l.append(pars[-3][7][5])
        self.HA_EW_l.append(pars[-3][7][6])

        self.HA_l_er.append(pars[-3][8][0])
        self.HB_l_er.append(pars[-3][8][1])
        self.OI_l_er.append(pars[-3][8][2])
        self.OIII_l_er.append(pars[-3][8][3])
        self.NII_l_er.append(pars[-3][8][4])
        self.SII_l_er.append(pars[-3][8][5])
        self.HA_EW_l_er.append(pars[-3][8][6])
        

    def file_out(self):
        Dict = {
            '0I' : self.OI,
            'GAMAID' : self.GAMAID,
            'SURVEY' : self.SURV,
            'SURVEY_CODE' : self.SURV_CODE,
            'IS_BEST' : self.IS_BEST,
            'IS_SBEST' : self.IS_SBEST,
            'QC' : self.coef,
            'met' : self.met,
            'met_er0' : self.col1,
            'met_er1' : self.col2,
            'SFR_HA' : self.SFR_HA,
            'SFR_HA_er' : self.SFR_HA_er,
            'age' : self.AGE,
            'AGN' : self.AGN,
            'HA' : self.HA_l,
            'HA_er' : self.HA_l_er,
            'HB' : self.HB_l,
            'HB_er' : self.HB_l_er,
            'OI' : self.OI_l,
            'OI_er' : self.OI_l_er,
            'OIII' : self.OIII_l,
            'OIII_er' : self.OIII_l_er,
            'NII' : self.NII_l,
            'NII_er' : self.NII_l_er,
            'SII' : self.SII_l,
            'SII_er' : self.SII_l_er,
            'HA_EW' : self.HA_EW_l,
            'HA_EW_er' : self.HA_EW_l_er,
        }

        df = pd.DataFrame(Dict)
        df.to_csv('out.csv', index=False)
    
    def plotting_arrows(self, ax, x, y, sym, age, mode):
        coord = [[-0.07, 0], [0.07, 0], [0, -0.07], [0, 0.07]]
        for item in sym:
            if item != 4 and mode == 'YES':
                ax.arrow(x, y, coord[item][0], coord[item][1], head_width=0.03, head_length=0.03, color='midnightblue', alpha=1)
            elif item != 4 and mode == 'UNC':
                ax.arrow(x, y, coord[item][0], coord[item][1], head_width=0.03, head_length=0.03, color='springgreen', alpha=1)
            elif item != 4 and mode == 'NO':
                ax.arrow(x, y, coord[item][0], coord[item][1], head_width=0.03, head_length=0.03, color='mediumvioletred', alpha=1)
            elif item != 4:
                ax.arrow(x, y, coord[item][0], coord[item][1], head_width=0.03, head_length=0.03, color=self.s_m.to_rgba(age), alpha=1)
            else:
                pass

    def plotting(self):

        gs_top = plt.GridSpec(1, 3, wspace=0)
        self.fig = plt.figure()

        self.ax1 = self.fig.add_subplot(gs_top[:,0])
        self.ax1.set_ylim([-3,3])
        self.ax2 = self.fig.add_subplot(gs_top[:,1], sharey=self.ax1)
        self.ax3 = self.fig.add_subplot(gs_top[:,2], sharey=self.ax1)

        self.topaxes = [self.ax1, self.ax2, self.ax3]

        for ax in self.topaxes[1:]:
            plt.setp(ax.get_yticklabels(), visible=False)

        self.ax1.set_ylabel('log[O[III]/H_B]')
        self.ax1.set_xlabel('log[N[II]/H_A]')
        self.ax2.set_xlabel('log[S[II]/H_A]')
        self.ax3.set_xlabel('log[O[I]/H_A]')

        self.axes = [self.ax1, self.ax2, self.ax3]

        ages = []
        for pars in self.flux_er_mod9:
            ages.append(pars[1])

        for item in self.axes:
            item.set_xlim([-2.5, 0.5])
        
        norm = mpl.colors.Normalize(vmin=8.8,vmax=10.25)
        # choose a colormap
        c_m = mpl.cm.jet
        # create a ScalarMappable and initialize a data structure
        self.s_m = mpl.cm.ScalarMappable(cmap=c_m, norm=norm)
        self.s_m.set_array([])

        for pars in self.flux_er_mod9:
            plots = pars[0] 
            age = pars[-2]
            mod = pars[-1]
            AGN = pars[-3][6]
            for i in range(len(plots)-1):
                if plots[i] == 0:
                    continue
                if len(plots[i]) == 4:

                    if AGN == 'YES':
                        self.axes[i].scatter(plots[i][0], plots[i][1], s=1.5, color='midnightblue', alpha=1)
                    elif AGN == 'UNC': 
                        self.axes[i].scatter(plots[i][0], plots[i][1], s=1.5, color='springgreen', alpha=1)
                    elif AGN == 'NO': 
                        self.axes[i].scatter(plots[i][0], plots[i][1], s=1.5, color='mediumvioletred', alpha=1)
                    else: 
                        self.axes[i].scatter(plots[i][0], plots[i][1], s=1.5, color=self.s_m.to_rgba(age), alpha=1) 
                    #self.axes[i].errorbar(plots[i][0], plots[i][1], xerr = plots[i][2], yerr = plots[i][3], fmt = 'o', color=self.s_m.to_rgba(age), markersize=2, alpha=0.2) 
                else:
                    Main.plotting_arrows(self, self.axes[i], plots[i][0], plots[i][1], plots[i][2], age, AGN)
                if mod > 0:
                    self.axes[i].text(plots[i][0], plots[i][1], pars[-1])
                    self.axes[i].scatter(plots[i][0], plots[i][1], color=self.s_m.to_rgba(age), marker="x")
                if mod == 0:
                    self.axes[i].scatter(plots[i][0], plots[i][1], color=self.s_m.to_rgba(age), marker="x")
                #print(plots[i])
        
        self.fig.subplots_adjust(right=0.85)
        cbar_ax = self.fig.add_axes([0.85, 0.15, 0.05, 0.7])
        for map in self.axes:
            self.fig.colorbar(self.s_m, cax=cbar_ax) 

        X_1 = np.arange(-4, 0.3, 0.01)
        X_111 = np.arange(-4, 0, 0.01)
        X_11 = np.arange(-0.45, 0.3, 0.01)
        self.ax1.plot(X_1, (0.61/(X_1 - 0.47)) + 1.19, c='crimson') #Kewley, 2001
        self.ax1.plot(X_111, (0.61/(X_111 - 0.05)) + 1.3, c='crimson', linestyle='dashed') #Kauffman, 2003
        self.ax1.plot(X_11, (2.145*(X_11+0.45)) - 0.5, c='crimson', linestyle='dotted') #https://adsabs.harvard.edu/full/2003MNRAS.346.1055K
        self.ax1.text(-1, 1.5, 'AGN/Seyfert')
        self.ax1.text(0.1, 0.2, 'LINER')
        self.ax1.text(-0.1, -1.5, 'Comp')
        self.ax1.text(-1.7, -1.5, 'SF')
        self.ax1.set_yticks(np.arange(-3, 3, 0.2))
        self.ax1.set_xticks(np.arange(-2.5, 0.5, 0.2))

        X_2 = np.arange(-4, 0.15, 0.01)
        X_22 = np.arange(-0.3, 0.5, 0.01)
        self.ax2.plot(X_2, (0.72/(X_2 - 0.32)) + 1.3, c='crimson')
        self.ax2.plot(X_22, (1.89*X_22) + 0.76, c='crimson', linestyle='dotted')
        self.ax2.text(-0.7, 1.5, 'Seyfert')
        self.ax2.text(0.1, 0.2, 'LINER')
        self.ax2.text(-1.7, -1.5, 'SF')
        self.ax2.set_xticks(np.arange(-2.5, 0.5, 0.2))

        X_3 = np.arange(-4, -0.7, 0.01)
        X_33 = np.arange(-1.1, 0.5, 0.01)
        self.ax3.plot(X_3, (0.73/(X_3 + 0.59)) + 1.33, c='crimson')
        self.ax3.plot(X_33, (1.18*X_33) + 1.3, c='crimson', linestyle='dotted')
        self.ax3.text(-1.5, 1.5, 'Seyfert')
        self.ax3.text(-0.5, -0.5, 'LINER')
        self.ax3.text(-1.7, -1.5, 'SF')
        self.ax3.set_xticks(np.arange(-2.5, 0.5, 0.2))
        plt.show()


    def plotting_2(self):
        
        print('NEW FIGURE')
        self.fig2 = plt.figure()
        self.ax4 = self.fig2.add_subplot()
        self.ax4.set_xlabel('log[N[II]/H_A]')
        self.ax4.set_ylabel('log(HA_EW)')

        self.ax4.set_xlim([-2.2, 1.2])
        self.ax4.set_ylim([-1.5, 3.2])

        norm = mpl.colors.Normalize(vmin=8.8,vmax=10.25)
        # choose a colormap
        c_m = mpl.cm.jet
        # create a ScalarMappable and initialize a data structure
        self.s_m = mpl.cm.ScalarMappable(cmap=c_m, norm=norm)
        self.s_m.set_array([])


        for pars in self.flux_er_mod9:
            plots = pars[0] 
            age = pars[-2]
            mod = pars[-1]
            AGN = pars[-3][6]
            if plots[3] == 0:
                continue
            if len(plots[3]) == 4:
                #self.ax4.errorbar(plots[3][0], plots[3][1], xerr = plots[3][2], yerr = plots[3][3], fmt = 'o', color=self.s_m.to_rgba(age), markersize=2, alpha=0.2) 
                i = 3
                if AGN == 'YES':
                    self.ax4.scatter(plots[i][0], plots[i][1], s=1.5, color='midnightblue', alpha=1)
                elif AGN == 'UNC': 
                    self.ax4.scatter(plots[i][0], plots[i][1], s=1.5, color='springgreen', alpha=1)
                elif AGN == 'NO': 
                    self.ax4.scatter(plots[i][0], plots[i][1], s=1.5, color='mediumvioletred', alpha=1)
                elif AGN == 'NOEL':
                    self.ax4.scatter(plots[i][0], plots[i][1], s=1.5, color='orchid', alpha=1)
                elif AGN == 'HBTRICK':
                    self.ax4.scatter(plots[i][0], plots[i][1], s=1.5, color='mediumblue', alpha = 1)
                else: 
                    self.ax4.scatter(plots[i][0], plots[i][1], s=1.5, color=self.s_m.to_rgba(age), alpha=1) 
            else:
                Main.plotting_arrows(self, self.ax4, plots[i][0], plots[i][1], plots[i][2], age, AGN)
            if mod > 0:
                i = 3
                self.ax4.text(plots[i][0], plots[i][1], pars[-1])
                #self.ax4.scatter(plots[i][0], plots[i][1], color=self.s_m.to_rgba(age), marker="x")
            #if mod == 0:
                #i = 3
                #self.ax4.scatter(plots[i][0], plots[i][1], color=self.s_m.to_rgba(age), marker="x")   
            
        self.fig2.subplots_adjust(right=0.85)
        cbar_ax = self.fig2.add_axes([0.85, 0.15, 0.05, 0.7])
        self.fig2.colorbar(self.s_m, cax=cbar_ax) 

        X_LLR = np.arange(-2, 1, 0.01)
        self.ax4.plot(X_LLR, -0.3+X_LLR*0, 'k--')
        self.ax4.text(-1.5, -1, 'LLR')

        self.ax4.plot(X_LLR, 0.45+X_LLR*0, 'black')
        self.ax4.text(-1.5, 0, 'ELR')

        X_wAGN = np.arange(-0.3, 1, 0.01)
        self.ax4.plot(X_wAGN, 0.8+X_wAGN*0, 'black')
        self.ax4.text(-1.5, 1.5, 'SF')
        self.ax4.text(0.5, 0.5, 'wAGN')
        self.ax4.text(0.5, 2.5, 'sAGN')

        Y_sAGN = np.arange(0.45, 3, 0.01)
        self.ax4.plot(-0.3+Y_sAGN*0, Y_sAGN, 'black')
        
        plt.show()

if __name__ == '__main__':
    obj = Main('GAMAforOleg.txt', 'DirectSummationv05', 'sample_out_Spectra.csv', 100)
    obj.ola_reading()
    obj.samples_get()
    obj.gama_reading()
    obj.sorting()
    obj.plotting()
    obj.plotting_2()
        