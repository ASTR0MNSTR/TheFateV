import csv
import matplotlib.pyplot as plt 
import numpy as np
import math

class hp:
    def log_er(item):
        return [[abs(math.log(1-item, 10))], [abs(math.log(1+item, 10))]]

class Main(hp):

    def __init__(self, ola_file, table, out):
        self.ola_file = ola_file
        self.out = out
        self.table = table
        self.base_dict = {}
        self.data_dict = []

    def reading(self):
        with open(self.ola_file, 'r') as input:
            lines = input.readlines()
        lines_stripped = [line.strip() for line in lines]
        for line in lines_stripped:
            GAMAID = int(line.split()[0])
            #SFR_0 = float(line.split()[95])
            #SFR_0 = float(line.split()[130])+float(line.split()[35])
            #SFR_0 = float(line.split()[135])+float(line.split()[35])
            SFR_0 = float(line.split()[140])+float(line.split()[35])
            MS = float(line.split()[35])
            TD_C = float(line.split()[47])
            TD_W = float(line.split()[53])
            #SFR_up = float(line.split()[96])
            #SFR_down = float(line.split()[94])
            #SFR_er = [[abs(SFR_down-SFR_0)], [abs(SFR_up-SFR_0)]]
            SFR_er = [[0], [0]]
            AGE_er = [[abs(float(line.split()[105]) - float(line.split()[104]))], [abs(float(line.split()[105]) - float(line.split()[106]))]]
            self.base_dict.update({GAMAID: [SFR_0, SFR_er, AGE_er, MS, TD_C, TD_W]})

        self.data_dict = []

        with open(self.out, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                self.data_dict.append({'OI': int(row['0I']), 'GAMAID': int(row['GAMAID']), 'met': float(
                    row['met']), 'met_er0': row['met_er0'], 'met_er1': row['met_er1'], 'SFR_HA': float(row['SFR_HA']), 'SFR_HA_er': row['SFR_HA_er'], 'age': float(row['age']), 'AGN': row['AGN']})

        mydata = []
        with open(self.table, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                mydata.append({'GAMAID': int(row['GAMAID']), 'SFR_0' : float(row['SFR']), 'SFR_0_er' : 0, 'AGE_er' : 0, 'Ms' : float(row['Ms']), 'TD_C' : float(row['Td']), 'TD_W' : -99})
        
        for item in mydata:
            self.base_dict.update({item['GAMAID'] : [item['SFR_0'], item['SFR_0_er'], item['AGE_er'], item['Ms'], item['TD_C'], item['TD_W']]})
                

    def matching(self):
        for item in self.data_dict:
            try:
                SFR_0, SFR_er, age_er, Ms, TD_C, TD_W = self.base_dict[item['GAMAID']]
                item.update({'SFR_0': SFR_0, 'SFR_0_er': SFR_er, 'age_er' : age_er, 'Ms' : Ms, 'TD_C' : TD_C, 'TD_W' : TD_W})
            except:
                print(item['GAMAID'])
    
    def arrow_plotting_1(self, ax, x, y, cmd, mode):
        coord = [[0, 0.07], [0, -0.07]]
        if mode == 'YES':
            if cmd == 'up':
                ax.arrow(x, y, coord[0][0], coord[0][1], head_width=0.03, head_length=0.03, alpha=0.2, color = 'midnightblue')
            else:
                ax.arrow(x, y, coord[1][0], coord[1][1], head_width=0.03, head_length=0.03, alpha=0.2, color = 'midnightblue')
        elif mode == 'UNC':
            if cmd == 'up':
                ax.arrow(x, y, coord[0][0], coord[0][1], head_width=0.03, head_length=0.03, alpha=0.2, color = 'springgreen')
            else:
                ax.arrow(x, y, coord[1][0], coord[1][1], head_width=0.03, head_length=0.03, alpha=0.2, color = 'springgreen')
        elif mode == 'NO':
            if cmd == 'up':
                ax.arrow(x, y, coord[0][0], coord[0][1], head_width=0.03, head_length=0.03, alpha=0.2, color = 'mediumvioletred')
            else:
                ax.arrow(x, y, coord[1][0], coord[1][1], head_width=0.03, head_length=0.03, alpha=0.2, color = 'mediumvioletred')
        elif mode == 'NOEL':
            if cmd == 'up':
                ax.arrow(x, y, coord[0][0], coord[0][1], head_width=0.03, head_length=0.03, alpha=0.2, color = 'orchid')
            else:
                ax.arrow(x, y, coord[1][0], coord[1][1], head_width=0.03, head_length=0.03, alpha=0.2, color = 'orchid')


    def arrow_plotting_2(self, ax, x, y, cmd, mode):
        coord = [[0, 0.1], [0, -0.1]]
        if mode == 'YES':
            if cmd == 'up':
                ax.arrow(x, y, coord[0][0], coord[0][1], head_width=0.03, head_length=0.03, alpha=0.2, color = 'midnightblue')
            else:
                ax.arrow(x, y, coord[1][0], coord[1][1], head_width=0.03, head_length=0.03, alpha=0.2, color = 'midnightblue')
        elif mode == 'UNC':
            if cmd == 'up':
                ax.arrow(x, y, coord[0][0], coord[0][1], head_width=0.03, head_length=0.03, alpha=0.2, color = 'springgreen')
            else:
                ax.arrow(x, y, coord[1][0], coord[1][1], head_width=0.03, head_length=0.03, alpha=0.2, color = 'springgreen')
        elif mode == 'NO':
            if cmd == 'up':
                ax.arrow(x, y, coord[0][0], coord[0][1], head_width=0.03, head_length=0.03, alpha=0.2, color = 'mediumvioletred')
            else:
                ax.arrow(x, y, coord[1][0], coord[1][1], head_width=0.03, head_length=0.03, alpha=0.2, color = 'mediumvioletred')
        elif mode == 'NOEL':
            if cmd == 'up':
                ax.arrow(x, y, coord[0][0], coord[0][1], head_width=0.03, head_length=0.03, alpha=0.2, color = 'orchid')
            else:
                ax.arrow(x, y, coord[1][0], coord[1][1], head_width=0.03, head_length=0.03, alpha=0.2, color = 'orchid')

    def plotting_agemet(self):

        self.fig1 = plt.figure()
        self.ax1 = self.fig1.add_subplot()

        self.ax1.set_ylabel('12 + log(O/H)')
        self.ax1.set_xlabel('log(age)')

        for item in self.data_dict:
            met = item['met']
            age = item['age']
            age_er = item['age_er']
            met_er0 = item['met_er0']
            met_er1 = item['met_er1']
            OI = item['OI']
            AGN = item['AGN']
            if met != -99.9:
                if met_er0 == '+' or met_er1 == '+':
                    Main.arrow_plotting_1(self, self.ax1, age, 12+met, 'up', AGN)
                elif met_er0 == '-' or met_er1 == '-':
                    Main.arrow_plotting_1(self, self.ax1, age, 12+met, 'down', AGN)
                else:
                    try:
                        if AGN == 'YES':
                            #self.ax1.errorbar(age, 12+met, xerr=age_er, yerr=[[abs(float(met_er0))], [abs(float(met_er1))]], fmt = 'o', alpha= 0.3, color = 'black')
                            self.ax1.scatter(age, 12+met, alpha= 0.3, color = 'midnightblue')
                        elif AGN == 'UNC':
                            #self.ax1.errorbar(age, 12+met, xerr=age_er, yerr=[[abs(float(met_er0))], [abs(float(met_er1))]], fmt = 'o', alpha= 0.3, color = 'tab:purple')
                            self.ax1.scatter(age, 12+met, alpha= 0.3, color = 'springgreen')
                        elif AGN == 'NO':
                            #self.ax1.errorbar(age, 12+met, xerr=age_er, yerr=[[abs(float(met_er0))], [abs(float(met_er1))]], fmt = 'o', alpha= 0.1, color = 'c')
                            self.ax1.scatter(age, 12+met, alpha= 0.3, color = 'mediumvioletred')
                        elif AGN == 'NOEL':
                            #self.ax1.errorbar(age, 12+met, xerr=age_er, yerr=[[abs(float(met_er0))], [abs(float(met_er1))]], fmt = 'o', alpha= 0.1, color = 'c')
                            self.ax1.scatter(age, 12+met, alpha= 0.3, color = 'orchid')
                    except:
                        print('Bad point')
                if int(OI) > 0:
                    self.ax1.text(age, 12+met, OI)
        x = np.arange(8.2, 10.1, 0.1)
        self.ax1.plot(x, 0*x + 8.69, color='tab:orange', linestyle='--', label='Metallicity of the Sun')

        #legend-making-stuff
        self.ax1.scatter(-99, -99, alpha= 0.3, color = 'midnightblue', label='YES')
        self.ax1.scatter(-99, -99, alpha= 0.3, color = 'springgreen', label='UNC')
        self.ax1.scatter(-99, -99, alpha= 0.3, color = 'mediumvioletred', label='NO')
        self.ax1.scatter(-99, -99, alpha= 0.3, color = 'orchid', label='NOEL')
        self.ax1.set_xlim(8.1, 10.1)
        self.ax1.set_ylim(5.9, 12.1)
        self.ax1.legend()
        print(len(self.data_dict))
        plt.show()


    def plotting_sfrs(self):

        self.fig2 = plt.figure()
        self.ax2 = self.fig2.add_subplot()

        self.ax2.set_ylabel('log(SFR_HA), M/yr')
        self.ax2.set_xlabel('log(SFR_GAMA), M/yr')
        self.ax2.set_aspect('equal', adjustable='box')

        for item in self.data_dict:
            SFR_HA_norm = item['SFR_HA']
            if SFR_HA_norm > 0:
                SFR_HA = math.log(item['SFR_HA'], 10)
                SFR_0 = item['SFR_0']
                SFR_0_er = item['SFR_0_er']
                SFR_HA_er = item['SFR_HA_er']
                OI = item['OI']
                AGN = item['AGN']
                if SFR_HA_er == '+':
                    Main.arrow_plotting_2(self, self.ax2, SFR_0, SFR_HA, 'up', AGN)
                elif SFR_HA_er == '-':
                    Main.arrow_plotting_2(self, self.ax2, SFR_0, SFR_HA, 'down', AGN)
                else:
                    if AGN == 'YES':
                        #self.ax2.errorbar(SFR_0, SFR_HA, xerr=SFR_0_er, yerr= hp.log_er(float(SFR_HA_er)/SFR_HA_norm), fmt = 'o', alpha = 0.3, color = 'black')
                        self.ax2.scatter(SFR_0, SFR_HA, alpha = 0.3, color = 'midnightblue')
                    elif AGN == 'UNC':
                        #self.ax1.errorbar(age, 12+met, xerr=age_er, yerr=[[abs(float(met_er0))], [abs(float(met_er1))]], fmt = 'o', alpha= 0.3, color = 'tab:purple')
                        self.ax2.scatter(SFR_0, SFR_HA, alpha= 0.3, color = 'springgreen')
                    elif AGN == 'NO':
                        #self.ax1.errorbar(age, 12+met, xerr=age_er, yerr=[[abs(float(met_er0))], [abs(float(met_er1))]], fmt = 'o', alpha= 0.1, color = 'c')
                        self.ax2.scatter(SFR_0, SFR_HA, alpha= 0.3, color = 'mediumvioletred')
                    elif AGN == 'NOEL':
                        #self.ax1.errorbar(age, 12+met, xerr=age_er, yerr=[[abs(float(met_er0))], [abs(float(met_er1))]], fmt = 'o', alpha= 0.1, color = 'c')
                        self.ax2.scatter(SFR_0, SFR_HA, alpha= 0.3, color = 'orchid')
                if int(OI) > 0:
                    self.ax2.text(SFR_0, SFR_HA, OI)
            
        #rowlands, young stars correction
        #dust correction

        self.ax2.plot([-4, 4], [-4, 4], linestyle='dashed', color='green')

        self.ax2.scatter(-99, -99, alpha= 0.3, color = 'midnightblue', label='YES')
        self.ax2.scatter(-99, -99, alpha= 0.3, color = 'springgreen', label='UNC')
        self.ax2.scatter(-99, -99, alpha= 0.3, color = 'mediumvioletred', label='NO')
        self.ax2.scatter(-99, -99, alpha= 0.3, color = 'orchid', label='NOEL')
        self.ax2.set_xlim(-4.2, 4.2)
        self.ax2.set_ylim(-4.2, 4.2)
        self.ax2.legend()
        
        plt.show()
    
    def plotting_smmet(self):

        self.fig3 = plt.figure()
        self.ax3 = self.fig3.add_subplot()

        self.ax3.set_ylabel('met')
        self.ax3.set_xlabel('Ms')
        MET = []
        MS = []
        for item in self.data_dict:
            met = item['met']
            Ms = item['Ms']
            met_er0 = item['met_er0']
            met_er1 = item['met_er1']
            OI = item['OI']
            AGN = item['AGN']
            if 12+met > 0 and Ms > 0:
                if met_er0 == '+' or met_er1 == '+':
                    Main.arrow_plotting_1(self, self.ax3, Ms, 12+met, 'up', AGN)
                elif met_er0 == '-' or met_er1 == '-':
                    Main.arrow_plotting_1(self, self.ax3, Ms, 12+met, 'down', AGN)
                else:
                    try:
                        if AGN == 'YES':
                            #self.ax1.errorbar(age, 12+met, xerr=age_er, yerr=[[abs(float(met_er0))], [abs(float(met_er1))]], fmt = 'o', alpha= 0.3, color = 'black')
                            self.ax3.scatter(Ms, 12+met, alpha= 0.6, color = 'midnightblue')
                        elif AGN == 'UNC':
                            #self.ax1.errorbar(age, 12+met, xerr=age_er, yerr=[[abs(float(met_er0))], [abs(float(met_er1))]], fmt = 'o', alpha= 0.3, color = 'tab:purple')
                            self.ax3.scatter(Ms, 12+met, alpha= 0.6, color = 'springgreen')
                        elif AGN == 'NO':
                            #self.ax1.errorbar(age, 12+met, xerr=age_er, yerr=[[abs(float(met_er0))], [abs(float(met_er1))]], fmt = 'o', alpha= 0.1, color = 'c')
                            self.ax3.scatter(Ms, 12+met, alpha= 0.6, color = 'mediumvioletred')
                            MET.append(12+met)
                            MS.append(Ms)
                        elif AGN == 'NOEL':
                            #self.ax1.errorbar(age, 12+met, xerr=age_er, yerr=[[abs(float(met_er0))], [abs(float(met_er1))]], fmt = 'o', alpha= 0.1, color = 'c')
                            self.ax3.scatter(Ms, 12+met, alpha= 0.6, color = 'orchid')
                            MET.append(12+met)
                            MS.append(Ms)
                    except:
                        print('Bad point in sm-met')

                if int(OI) > 0:
                    self.ax3.text(Ms, 12+met, OI)

        #fitting part
        MET_ar = np.array(MET)
        MS_ar = np.array(MS)

        z = np.polyfit(MS_ar, MET_ar, 3)

        coefs = z.tolist()
        print(coefs)

        x = np.arange(8.5, 11.5, 0.01)
        
        self.ax3.plot(x, -8.771 + (4.15003*x) - (0.322156*(x**2)) + (0.00818179*(x**3)), color = 'green', linestyle = 'dashed', label = ' Denicol´o et al. (2002)')
        self.ax3.plot(x, coefs[3] + (coefs[2]*x) + (coefs[1]*(x**2)) + coefs[0]*(x**3), color = 'tab:blue', linestyle = 'dashed', label = 'My fit to NO and NOEL')
        self.ax3.plot(x, 0*x + 8.69, color='tab:orange', linestyle='--', label='Metallicity of the Sun')

        self.ax3.scatter(-99, -99, alpha= 0.6, color = 'midnightblue', label='YES')
        self.ax3.scatter(-99, -99, alpha= 0.6, color = 'springgreen', label='UNC')
        self.ax3.scatter(-99, -99, alpha= 0.2, color = 'mediumvioletred', label='NO')
        self.ax3.scatter(-99, -99, alpha= 0.2, color = 'orchid', label='NOEL')
        self.ax3.set_xlim(8, 12)
        self.ax3.set_ylim(7, 10)
        self.ax3.legend()
        
        plt.show()

    def plotting_tdage(self):

        self.fig4 = plt.figure()
        self.ax4 = self.fig4.add_subplot()

        self.ax4.set_ylabel('Td')
        self.ax4.set_xlabel('age')
        for item in self.data_dict:
            TD_C = item['TD_C']
            age = item['age']
            OI = item['OI']
            AGN = item['AGN']
            try:
                        if AGN == 'YES':
                            #self.ax1.errorbar(age, 12+met, xerr=age_er, yerr=[[abs(float(met_er0))], [abs(float(met_er1))]], fmt = 'o', alpha= 0.3, color = 'black')
                            self.ax4.scatter(age, TD_C, alpha= 0.6, color = 'midnightblue')
                        elif AGN == 'UNC':
                            #self.ax1.errorbar(age, 12+met, xerr=age_er, yerr=[[abs(float(met_er0))], [abs(float(met_er1))]], fmt = 'o', alpha= 0.3, color = 'tab:purple')
                            self.ax4.scatter(age, TD_C, alpha= 0.6, color = 'springgreen')
                        elif AGN == 'NO':
                            #self.ax1.errorbar(age, 12+met, xerr=age_er, yerr=[[abs(float(met_er0))], [abs(float(met_er1))]], fmt = 'o', alpha= 0.1, color = 'c')
                            self.ax4.scatter(age, TD_C, alpha= 0.6, color = 'mediumvioletred')
                        elif AGN == 'NOEL':
                            #self.ax1.errorbar(age, 12+met, xerr=age_er, yerr=[[abs(float(met_er0))], [abs(float(met_er1))]], fmt = 'o', alpha= 0.1, color = 'c')
                            self.ax4.scatter(age, TD_C, alpha= 0.6, color = 'orchid')
            except:
                        print('Bad point in TDAGE')

            if int(OI) > 0:
                    self.ax4.text(age, TD_C, OI)

        #self.ax4.scatter(-99, -99, alpha= 0.6, color = 'midnightblue', label='YES')
        #self.ax4.scatter(-99, -99, alpha= 0.6, color = 'springgreen', label='UNC')
        #self.ax4.scatter(-99, -99, alpha= 0.2, color = 'mediumvioletred', label='NO')
        #self.ax4.scatter(-99, -99, alpha= 0.2, color = 'orchid', label='NOEL')
        #self.ax4.legend()
        
        plt.show()

if __name__ == '__main__':
    obj = Main('GAMAforOleg.txt', 'sample_out_spectra.csv', 'out.csv')
    obj.reading()
    obj.matching()
    #obj.plotting_agemet()
    #bj.plotting_sfrs()
    #obj.plotting_smmet()
    obj.plotting_tdage()
