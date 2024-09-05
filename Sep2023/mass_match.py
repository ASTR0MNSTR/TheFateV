import csv
import matplotlib.pyplot as plt 
import numpy as np
import math
import statistics as st
from scipy.stats import pearsonr
from scipy.stats import sem
import pandas as pd
from scipy.optimize import curve_fit
import random
from __legpars__ import *


class Main:

    def __init__(self, SDSS_from_GAMA, GAMAv3):
        self.SDSS_from_GAMA = SDSS_from_GAMA
        self.GAMAv3 = GAMAv3

        self.reff_dict = {}
        self.spec_ids = []
        self.df = None
    
    def reading_files(self):
        self.df = pd.read_csv(self.SDSS_from_GAMA)
        for item in self.df['SPECID']:
            self.spec_ids.append(str(item))
        print(type(self.spec_ids[3]))

        count = 0
        with open(self.GAMAv3, 'r') as f:
            header = f.readline()
            while True:
                line = f.readline()
                if not line:
                    break
                count += 1
                Main.processing_lines(self, line)

    def processing_lines(self, line):
        line = line.strip()
        SPEC_ID = line.split()[734]
        SURVEY = line.split()[740]
        if SURVEY == 'SDSS' or SURVEY == 'GAMA':
            if SPEC_ID in self.spec_ids:
                Y = float(line.split()[266])
                ER = float(line.split()[273])

                if ER == -9999.0 or 2*ER > Y:
                    Y = -99999.0
                    ER = -99999.0
                    trust = -1
                else:
                    if SURVEY == 'GAMA':
                        if (Y - ER) > 2:
                            trust = 1
                        elif (Y + ER) < 2:
                            trust = 0
                        else:
                            trust = -1 
                    elif SURVEY == 'SDSS':
                        if (Y - ER) > 3:
                            trust = 1
                        elif (Y + ER) < 3:
                            trust = 0
                        else:
                            trust = -1 
                
                self.reff_dict.update({SPEC_ID: [Y, ER, trust]})

    def appending_csv(self):
        spec_ids_out = []
        reff = []
        reff_err = []
        trust_val = []
        for key in self.reff_dict.keys():
            spec_ids_out.append(key)
            reff.append(self.reff_dict[key][0])
            reff_err.append(self.reff_dict[key][1])
            trust_val.append(self.reff_dict[key][2])

        Dict_for_dataframe ={
            'SPECID' : spec_ids_out,
            'GALRE_r' : reff,
            'GALREERR_r' : reff_err,
            'Aperture_1_Reff': trust_val
        }

        df2 = pd.DataFrame.from_dict(Dict_for_dataframe)
        
        csv_out = pd.merge(self.df, df2, how="inner", on="SPECID")
        csv_out.to_csv('E:/databases/GAMA_ETG_OLA_R_r_1.csv', index=False)

if __name__ == '__main__':
    obj = Main(r'E:\backup\backup_BPT\GAMA_ETG_OLA.csv', r'E:\LICENSE\ProgsData\main\GAMAv3.txt')
    obj.reading_files()
    obj.appending_csv() 