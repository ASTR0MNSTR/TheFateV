import numpy as np
import pandas as pd
from astropy.cosmology import Planck13 as cosmo
from __algo__ import *
from __plt__ import *

def aperture_flagging(REFF, REFF_err, SURVEY):
    apers = []
    for i, surv in enumerate(SURVEY):
        ER = REFF_err[i]
        Y = REFF[i]
        if ER == -9999.0 or 2*ER > Y:
            Y = -99999.0
            ER = -99999.0
            trust = -1
        else:
            if surv == 'GAMA':
                if (Y - ER) > 2:
                    trust = 1
                elif (Y + ER) < 2:
                    trust = 0
                else:
                    trust = -1 
            elif surv == 'SDSS':
                if (Y - ER) > 3:
                    trust = 1
                elif (Y + ER) < 3:
                    trust = 0
                else:
                    trust = -1
        apers.append(trust)
    return apers
        

def MS_flagging(SFR_array, MS_array, Z_array):
    bms = []
    for i in range(len(Z_array)):
        if Z_array[i] != -99999.0 and MS_array[i] != -99999.0 and SFR_array[i] != -99999.0: 
            time = cosmo.age(float(Z_array[i])).value
            SFR_MS = (0.84 - 0.026*time)*float(MS_array[i]) - (6.51 - 0.11*time)
            delta_SFR = float(SFR_array[i]) - SFR_MS
            bms.append(delta_SFR)
        else:
            bms.append(-99999.0)
    
    bms = np.array(bms)
    return bms

def morphology_flagging(index_array, index_err_array, Z_array):
    morph = []
    for i in range(len(index_array)):
        if 0.01 < Z_array[i] < 0.32:
            if index_array[i] > 2*index_err_array[i] and index_array[i] > 0 and index_err_array[i] >= 0:
                if index_array[i] + 2*index_err_array[i] < 2.5:
                    morph.append('S')
                elif index_array[i] - 2*index_err_array[i] > 4:
                    morph.append('E')
                else:
                    morph.append('U')
            else:
                morph.append('U')  
        else:
            morph.append('U')
            
    return morph
    

def spectra_processing(OIII, OIII_er, HB, HB_er, NII, NII_er, HA, HA_er, HA_EW, HA_EW_er, Z):
    lines = [OIII, OIII_er, HB, HB_er, NII, NII_er, HA, HA_er, HA_EW, HA_EW_er, Z]
    if '""' not in lines:
        lines_float = [float(item) for item in lines]
        if -99999.0 not in lines_float:
            OIII, OIII_er, HB, HB_er, NII, NII_er, HA, HA_er, HA_EW, HA_EW_er, Z = lines_float
            BPT, X, pair_x_flags, Y, pair_y_flags, WHAN, LOIII, LOIII_er, HA_ew, HA_ew_err, pair_HA = AGN_reg(OIII, OIII_er, HB, HB_er, NII, NII_er, HA, HA_er, HA_EW, HA_EW_er, Z)
        else:
            BPT = 'NDA'
            WHAN = 'NDA'
            LOIII = -99999.0
            LOIII_er = -99999.0
    else:
        BPT = 'PROB'
        WHAN = 'PROB'
        LOIII = -99999.0
        LOIII_er = -99999.0
        
    return BPT, WHAN, LOIII, LOIII_er  

def outflow_wAGN(SFR, LOIII, MS):
    return 1.14*math.log10(0.52*(10**SFR) + 0.51*(10**(LOIII + np.log10(3500) - 43))) - 0.41*math.log10((10**(MS))/(10**11))

def outflow_woAGN(SFR, MS):
    return 1.14*math.log10(0.52*(10**(float(SFR)))) - 0.41*math.log10((10**(float(MS)))/(10**11)) 

class Main:
    
    def __init__(self, path):
        self.path = path
        self.out_path = r'E:/databases/GAMA_MACHINE.csv'
        self.final_out_path = r'E:/databases/MRT/MRT_old.csv'
        
    def completely_different_extraction(self):
        DirectSummation_path = r"E:\databases\GAMAs4\DirectSummationv05"
        DS_cols ={
            'SPECID' : 0,
            'CATAID' : 1,
            'RA' : 2,
            'DEC' : 3,
            'Z' : 4,
            'SURVEY' : 6,
            'IS_BEST' : 8,
            'IS_SBEST' : 9,
            'NIIR_FLUX_ERR': 71,
            'NIIR_FLUX' : 72,
            'HA_EW_ERR' : 79,
            'HA_EW' : 80,
            'HA_FLUX_ERR' : 81,
            'HA_FLUX' : 82,
            'OIIIR_FLUX_ERR' : 166,
            'OIIIR_FLUX' : 167,
            'HB_FLUX_ERR' : 181,
            'HB_FLUX' : 182
        }
        DirectSummation = pd.read_csv(DirectSummation_path, sep=r"\s+", engine='python', usecols=DS_cols.values(), names=DS_cols.keys())
        DirectSummation = DirectSummation[DirectSummation.IS_BEST == True]
        DirectSummation = DirectSummation[DirectSummation.IS_SBEST == True]
        print(DirectSummation.shape)
        # print(DirectSummation)
        
        MagPhys_path = r'E:\databases\GAMAs4\MagPhysv06'
        MP_cols = {
            'CATAID' : 0,
            'mass_stellar_percentile16' : 34,
            'mass_stellar_percentile50' : 35,
            'mass_stellar_percentile84' : 36,
            'SFR_0_1Gyr_percentile16' : 94,
            'SFR_0_1Gyr_percentile50' : 95,
            'SFR_0_1Gyr_percentile84' : 96
        }
        
        MagPhys = pd.read_csv(MagPhys_path, sep=r"\s+", engine='python', usecols=MP_cols.values(), names=MP_cols.keys())
        print(MagPhys.shape)
        
        Sersic_path = r'E:\databases\GAMAs4\SersicCatSDSSv09'
        S_cols = {
            'CATAID' : 0,
            'GALINDEX_r' : 86,
            'GALINDEXERR_r' : 93
        }
        Sersic = pd.read_csv(Sersic_path, sep=r"\s+", engine='python', usecols=S_cols.values(), names=S_cols.keys())
        print(Sersic.shape)
        # print(Sersic)
        
        Result1 = pd.merge(DirectSummation, MagPhys, how="left", on='CATAID')
        print(Result1.shape)
        FinalDataFrame = pd.merge(Result1, Sersic, how="left", on='CATAID')
        print(FinalDataFrame.shape)
        FinalDataFrame.fillna(-99999.0, inplace=True)
        FinalDataFrame.RA = FinalDataFrame.RA.round(6)
        FinalDataFrame.DEC = FinalDataFrame.DEC.round(6)
        FinalDataFrame.Z = FinalDataFrame.Z.round(8)
        FinalDataFrame.to_csv(self.out_path, index=False)
        
        FinalDataFrame.info()
        print('Finished merge!')

    def read_process(self):
        DataFrame = pd.read_csv(self.out_path, sep=',')
        MS_flag = MS_flagging(DataFrame['SFR_0_1Gyr_percentile50'], DataFrame['mass_stellar_percentile50'], DataFrame['Z'])
        print('deltaMS calculated!')
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
            BPT, WHAN, LOIII, LOIII_er = spectra_processing(DataFrame['OIIIR_FLUX'][i], DataFrame['OIIIR_FLUX_ERR'][i], DataFrame['HB_FLUX'][i], DataFrame['HB_FLUX_ERR'][i],
                                                          DataFrame['NIIR_FLUX'][i], DataFrame['NIIR_FLUX_ERR'][i], DataFrame['HA_FLUX'][i], DataFrame['HA_FLUX_ERR'][i],
                                                          DataFrame['HA_EW'][i], DataFrame['HA_EW_ERR'][i], DataFrame['Z'][i])
            
            if LOIII == -99999.0 or DataFrame['SFR_0_1Gyr_percentile50'][i] == -99999.0 or DataFrame['mass_stellar_percentile50'][i] == -99999.0:
                OUTFLOW_up.append(-99999.0)
                OUTFLOW.append(-99999.0)
                OUTFLOW_down.append(-99999.0)
            else:    
                if LOIII_er < 0:
                    OUTFLOW_up.append(-100000.0)
                    OUTFLOW.append(outflow_wAGN(DataFrame['SFR_0_1Gyr_percentile50'][i], LOIII, DataFrame['mass_stellar_percentile50'][i]))
                    OUTFLOW_down.append(-100000.0)
                else:
                    OUTFLOW_up.append(outflow_wAGN(DataFrame['SFR_0_1Gyr_percentile50'][i], np.log10(10**LOIII + 10**LOIII_er), DataFrame['mass_stellar_percentile50'][i]))
                    OUTFLOW.append(outflow_wAGN(DataFrame['SFR_0_1Gyr_percentile50'][i], LOIII, DataFrame['mass_stellar_percentile50'][i]))
                    OUTFLOW_down.append(outflow_wAGN(DataFrame['SFR_0_1Gyr_percentile50'][i], np.log10(10**LOIII - 10**LOIII_er), DataFrame['mass_stellar_percentile50'][i]))
            
            if DataFrame['SFR_0_1Gyr_percentile50'][i] == -99999.0 or DataFrame['mass_stellar_percentile50'][i] == -99999.0:
                OUTFLOW_up_off.append(-99999.0)
                OUTFLOW_off.append(-99999.0)
                OUTFLOW_down_off.append(-99999.0) 
            else:               
                OUTFLOW_up_off.append(outflow_woAGN(DataFrame['SFR_0_1Gyr_percentile84'][i], DataFrame['mass_stellar_percentile16'][i]))
                OUTFLOW_off.append(outflow_woAGN(DataFrame['SFR_0_1Gyr_percentile50'][i], DataFrame['mass_stellar_percentile50'][i]))
                OUTFLOW_down_off.append(outflow_woAGN(DataFrame['SFR_0_1Gyr_percentile16'][i], DataFrame['mass_stellar_percentile84'][i]))
            
            BPTs.append(BPT)
            WHANs.append(WHAN)
            LAGNs.append(LOIII)
            LAGN_ers.append(LOIII_er)
        
        print('Spectra processing finished!')
        
        NewFrameDict = {
            'SPECID' : DataFrame['SPECID'],
            'CATAID' : DataFrame['CATAID'],
            'RA' : DataFrame['RA'],
            'DEC' : DataFrame['DEC'],
            'Z' : DataFrame['Z'],
            'SURVEY' : DataFrame['SURVEY'],
            'deltaMS' : MS_flag,
            'BPT' : BPTs,
            'WHAN' : WHANs, 
            'GALINDEX_r' : DataFrame['GALINDEX_r'],
            'GALINDEXERR_r' : DataFrame['GALINDEXERR_r'],
            'out_on_16' : OUTFLOW_down,
            'out_on_50' : OUTFLOW, 
            'out_on_84' : OUTFLOW_up, 
            'out_off_16' : OUTFLOW_down_off, 
            'out_off_50' :  OUTFLOW_off,
            'out_off_84' : OUTFLOW_up_off 
        }
        
        NewDataFrame = pd.DataFrame(NewFrameDict)
        
        #rounding part
        # NewDataFrame = NewDataFrame.round(decimals=8)
        # NewDataFrame.RA = NewDataFrame.RA.round(5)
        # NewDataFrame.DEC = NewDataFrame.DEC.round(5)
        # NewDataFrame.Z = NewDataFrame.Z.round(8)
        NewDataFrame.deltaMS = NewDataFrame.deltaMS.round(3)
        NewDataFrame.GALINDEX_r = NewDataFrame.GALINDEX_r.round(4)
        NewDataFrame.GALINDEXERR_r = NewDataFrame.GALINDEXERR_r.round(4)
        NewDataFrame.out_on_16 = NewDataFrame.out_on_16.round(3)
        NewDataFrame.out_on_50 = NewDataFrame.out_on_50.round(3)
        NewDataFrame.out_on_84 = NewDataFrame.out_on_84.round(3)
        NewDataFrame.out_off_16 = NewDataFrame.out_off_16.round(3)
        NewDataFrame.out_off_50 = NewDataFrame.out_off_50.round(3)
        NewDataFrame.out_off_84 = NewDataFrame.out_off_84.round(3)
        
        NewDataFrame.to_csv(self.final_out_path, index=False)
        
if __name__ == '__main__':
    obj = Main(r'E:/LICENSE/ProgsData/main/GAMAv3.txt')
    obj.completely_different_extraction()
    obj.read_process()