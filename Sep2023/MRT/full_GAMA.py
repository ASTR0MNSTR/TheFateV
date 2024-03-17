import numpy as np
import pandas as pd
from astropy.cosmology import Planck13 as cosmo
from __algo__ import *

def MS_flagging(SFR_array, MS_array, Z_array):
    bms = []
    for i in range(len(Z_array)):
        if Z_array[i] != '""' and MS_array[i] != '""' and SFR_array[i] != '""': 
            time = cosmo.age(float(Z_array[i])).value
            SFR_MS = (0.84 - 0.026*time)*float(MS_array[i]) - (6.51 - 0.11*time)
            delta_SFR = float(SFR_array[i]) - SFR_MS
        
            if delta_SFR < -0.2:
                bms.append(0)
            else:
                bms.append(1)
        else:
            bms.append(-1)
    
    bms = np.array(bms)
    return bms

def spectra_processing(OIII, OIII_er, HB, HB_er, NII, NII_er, HA, HA_er, HA_EW, HA_EW_er, Z):
    lines = [OIII, OIII_er, HB, HB_er, NII, NII_er, HA, HA_er, HA_EW, HA_EW_er, Z]
    if '""' not in lines:
        lines_float = [float(item) for item in lines]
        if -99999.0 not in lines_float:
            OIII, OIII_er, HB, HB_er, NII, NII_er, HA, HA_er, HA_EW, HA_EW_er, Z = lines_float
            BPT, X, pair_x_flags, Y, pair_y_flags, WHAN, LAGN, LAGN_er, pair_HA = AGN_reg(OIII, OIII_er, HB, HB_er, NII, NII_er, HA, HA_er, HA_EW, HA_EW_er, Z)
        else:
            BPT = 'NDA'
            WHAN = 'NDA'
            LAGN = -99999.0
            LAGN_er = -99999.0
    else:
        BPT = 'PROB'
        WHAN = 'PROB'
        LAGN = -99999.0
        LAGN_er = -99999.0
        
    return BPT, WHAN, LAGN, LAGN_er  

def outflow_wAGN(SFR, LAGN, MS):
    try: 
        return 1.14*math.log10(0.52*(10**float(SFR)) + 0.51*float(LAGN)) - 0.41*math.log10((10**(float(MS)))/(10**11))
        # return 1.14*math.log10(0.52*(10**SFR)) - 0.41*(MS - 11)
    except:
        return -99

def outflow_woAGN(SFR, MS):
    try: 
        return 1.14*math.log10(0.52*(10**float(SFR))) - 0.41*math.log10((10**(float(MS)))/(10**11))
    except:
        return -99              

class Main:
    
    def __init__(self, path):
        self.path = path
        self.out_path = r'E:/databases/GAMA_MACHINE.csv'
        self.final_out_path = r'E:/databases/MRT.csv'
        
    def extract(self):
        required_cols = [
            'CATAID_2b',
            'RA_1b',
            'DEC_1b',
            'Z_2',
            'NIIR_FLUX',
            'NIIR_FLUX_ERR',
            'OIIIR_FLUX',
            'OIIIR_FLUX_ERR',
            'HA_FLUX',
            'HA_FLUX_ERR',
            'HA_EW',
            'HA_EW_ERR',
            'HB_FLUX',
            'HB_FLUX_ERR',
            'SFR_0_1Gyr_percentile16',
            'SFR_0_1Gyr_percentile50',
            'SFR_0_1Gyr_percentile84',
            'mass_stellar_percentile16',
            'mass_stellar_percentile50',
            'mass_stellar_percentile84',
        ]
        
        DataFrame = pd.read_csv(self.path, sep="\s+", engine='python', usecols=required_cols)
        print(DataFrame.shape)
        DataFrame.to_csv(self.out_path, index=False)

    def read(self):
        DataFrame = pd.read_csv(self.out_path, sep=',')
        MS_flag = MS_flagging(DataFrame['SFR_0_1Gyr_percentile50'], DataFrame['mass_stellar_percentile50'], DataFrame['Z_2'])
        BPTs = []
        WHANs = []
        LAGNs = []
        LAGN_ers = []
        OUTFLOW = []
        OUTFLOW_up = []
        OUTFLOW_down = []
        OUTFLOW_off = []
        OUTFLOW_up_off = []
        OUTFLOW_down_off = []
        
        for i in range(len(DataFrame['OIIIR_FLUX'])):
            BPT, WHAN, LAGN, LAGN_er = spectra_processing(DataFrame['OIIIR_FLUX'][i], DataFrame['OIIIR_FLUX_ERR'][i], DataFrame['HB_FLUX'][i], DataFrame['HB_FLUX_ERR'][i],
                                                          DataFrame['NIIR_FLUX'][i], DataFrame['NIIR_FLUX_ERR'][i], DataFrame['HA_FLUX'][i], DataFrame['HA_FLUX_ERR'][i],
                                                          DataFrame['HA_EW'][i], DataFrame['HA_EW_ERR'][i], DataFrame['Z_2'][i])
                
            if LAGN_er in ['upNon', 'upAbs']:
                OUTFLOW.append(outflow_wAGN(DataFrame['SFR_0_1Gyr_percentile50'][i], LAGN, DataFrame['mass_stellar_percentile50'][i]))
                OUTFLOW_up.append(LAGN_er)
                OUTFLOW_down.append(LAGN_er)
            else:
                OUTFLOW.append(outflow_wAGN(DataFrame['SFR_0_1Gyr_percentile50'][i], LAGN, DataFrame['mass_stellar_percentile50'][i]))
                OUTFLOW_up.append(outflow_wAGN(DataFrame['SFR_0_1Gyr_percentile50'][i], LAGN + LAGN_er, DataFrame['mass_stellar_percentile50'][i]))
                OUTFLOW_down.append(outflow_wAGN(DataFrame['SFR_0_1Gyr_percentile50'][i], LAGN - LAGN_er, DataFrame['mass_stellar_percentile50'][i]))
            
            OUTFLOW_off.append(outflow_woAGN(DataFrame['SFR_0_1Gyr_percentile50'][i], DataFrame['mass_stellar_percentile50'][i]))
            OUTFLOW_up_off.append(outflow_woAGN(DataFrame['SFR_0_1Gyr_percentile84'][i], DataFrame['mass_stellar_percentile16'][i]))
            OUTFLOW_down_off.append(outflow_woAGN(DataFrame['SFR_0_1Gyr_percentile16'][i], DataFrame['mass_stellar_percentile84'][i]))
            
            BPTs.append(BPT)
            WHANs.append(WHAN)
            LAGNs.append(LAGN)
            LAGN_ers.append(LAGN_er)
        
        
        NewFrameDict = {
            'CATAID' : DataFrame['CATAID_2b'],
            'RA' : DataFrame['RA_1b'],
            'DEC' : DataFrame['DEC_1b'],
            'Z' : DataFrame['Z_2'],
            'bMS/MS' : MS_flag,
            'BPT' : BPTs,
            'WHAN' : WHANs, 
            'outflow_agn_on_percentile16' : OUTFLOW_down,
            'outflow_agn_on_percentile50' : OUTFLOW, 
            'outflow_agn_on_percentile84' : OUTFLOW_up, 
            'outflow_agn_off_percentile16' : OUTFLOW_down_off, 
            'outflow_agn_off_percentile50' :  OUTFLOW_off,
            'outflow_agn_off_percentile84' : OUTFLOW_up_off 
        }
        
        NewDataFrame = pd.DataFrame(NewFrameDict)
        NewDataFrame.to_csv(self.final_out_path, index=False)
        
if __name__ == '__main__':
    obj = Main(r'E:/LICENSE/ProgsData/main/GAMAv3.txt')
    # obj.extract()
    obj.read()