import numpy as np
import pandas as pd
from astropy.cosmology import Planck13 as cosmo
from __algo__ import *
from __plt__ import *

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

def spectra_processing(OIII, OIII_er, HB, HB_er, NII, NII_er, HA, HA_er, HA_EW, HA_EW_er, Z):
    lines = [OIII, OIII_er, HB, HB_er, NII, NII_er, HA, HA_er, HA_EW, HA_EW_er, Z]
    if '""' not in lines:
        lines_float = [float(item) for item in lines]
        if -99999.0 not in lines_float:
            OIII, OIII_er, HB, HB_er, NII, NII_er, HA, HA_er, HA_EW, HA_EW_er, Z = lines_float
            BPT, X, pair_x_flags, Y, pair_y_flags, WHAN, LAGN, LAGN_er, HA_ew, HA_ew_err, pair_HA = AGN_reg(OIII, OIII_er, HB, HB_er, NII, NII_er, HA, HA_er, HA_EW, HA_EW_er, Z)
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
        return 1.14*math.log10(0.52*(10**(float(SFR))) + 0.51*float(LAGN)) - 0.41*math.log10((10**(float(MS)))/(10**11))
        # return 1.14*math.log10(0.52*(10**SFR)) - 0.41*(MS - 11)
    except:
        return -99999.0

def outflow_woAGN(SFR, MS):
    try: 
        return 1.14*math.log10(0.52*(10**(float(SFR)))) - 0.41*math.log10((10**(float(MS)))/(10**11))
    except:
        return -99999.0         

class Main:
    
    def __init__(self, path):
        self.path = path
        self.out_path = r'E:/databases/GAMA_MACHINE.csv'
        self.final_out_path = r'E:/databases/MRT.csv'
        
    def completely_different_extraction(self):
        DirectSummation_path = r'E:\LICENSE\ProgsData\main\DirectSummationv05'
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
        DirectSummation = pd.read_csv(DirectSummation_path, sep="\s+", engine='python', usecols=DS_cols.values(), names=DS_cols.keys())
        DirectSummation = DirectSummation[DirectSummation.IS_BEST == True]
        DirectSummation = DirectSummation[DirectSummation.IS_SBEST == True]
        # print(DirectSummation)
        
        MagPhys_path = r'E:\databases\MagPhys\MagPhys'
        MP_cols = {
            'CATAID' : 0,
            'mass_stellar_percentile16' : 34,
            'mass_stellar_percentile50' : 35,
            'mass_stellar_percentile84' : 36,
            'SFR_0_1Gyr_percentile16' : 94,
            'SFR_0_1Gyr_percentile50' : 95,
            'SFR_0_1Gyr_percentile84' : 96
        }
        
        MagPhys = pd.read_csv(MagPhys_path, sep="\s+", engine='python', usecols=MP_cols.values(), names=MP_cols.keys())
        print(MagPhys)
        
        Sersic_path = r'E:\databases\SersicCatSDSS'
        S_cols = {
            'CATAID' : 0,
            'GALINDEX_r' : 86,
            'GALINDEXERR_r' : 93
        }
        Sersic = pd.read_csv(Sersic_path, sep="\s+", engine='python', usecols=S_cols.values(), names=S_cols.keys())
        # print(Sersic)
        
        Result1 = pd.merge(DirectSummation, MagPhys, how="left", on='CATAID')
        FinalDataFrame = pd.merge(Result1, Sersic, how="left", on='CATAID')
        FinalDataFrame.fillna(-99999.0, inplace=True)
        FinalDataFrame.to_csv(self.out_path, index=False)
        print('Finished merge!')
        
    # def extract(self):
    #     # required_cols = [
    #     #     'CATAID_2b',
    #     #     'RA_1b',
    #     #     'DEC_1b',
    #     #     'Z_2',
    #     #     'NIIR_FLUX',
    #     #     'NIIR_FLUX_ERR',
    #     #     'OIIIR_FLUX',
    #     #     'OIIIR_FLUX_ERR',
    #     #     'HA_FLUX',
    #     #     'HA_FLUX_ERR',
    #     #     'HA_EW',
    #     #     'HA_EW_ERR',
    #     #     'HB_FLUX',
    #     #     'HB_FLUX_ERR',
    #     #     'SFR_0_1Gyr_percentile16',
    #     #     'SFR_0_1Gyr_percentile50',
    #     #     'SFR_0_1Gyr_percentile84',
    #     #     'mass_stellar_percentile16',
    #     #     'mass_stellar_percentile50',
    #     #     'mass_stellar_percentile84',
    #     #     'GALINDEX_r',
    #     #     'GALINDEXERR_r'
    #     # ]
        
    #     required_cols_2 = [
    #         'CATAID_2b',
    #         'RA_1b',
    #         'DEC_1b',
    #         'Z_2',
    #         'NIIR_FLUX',
    #         'NIIR_FLUX_ERR',
    #         'OIIIR_FLUX',
    #         'OIIIR_FLUX_ERR',
    #         'HA_FLUX',
    #         'HA_FLUX_ERR',
    #         'HA_EW',
    #         'HA_EW_ERR',
    #         'HB_FLUX',
    #         'HB_FLUX_ERR',
    #         'mass_stellar_percentile50',
    #         'SFR_0_1Gyr_percentile16',
    #         'SFR_0_1Gyr_percentile50',
    #         'SFR_0_1Gyr_percentile84',
    #         'GALINDEX_r',
    #         'GALINDEXERR_r'
    #     ]
        
    #     DataFrame = pd.read_csv(self.path, sep="\s+", engine='python', usecols=required_cols_2)
    #     DataFrame = DataFrame.drop(DataFrame[DataFrame['CATAID_2b'] == '""'].index)
    #     print(DataFrame.shape)
    #     DataFrame.to_csv(self.out_path, index=False)
    
    # def reader(self):
    #     DataFrame = pd.read_csv(self.out_path, sep=',')
    #     BPTs = []
    #     WHANs = []
    #     for i in range(len(DataFrame['OIIIR_FLUX'])):
    #         if float(DataFrame['GALINDEX_r'][i]) + 2*float(DataFrame['GALINDEXERR_r'][i]) < 2.5 and float(DataFrame['GALINDEXERR_r'][i]) > 0 and float(DataFrame['GALINDEX_r'][i]) > 0 and float(DataFrame['Z_2'][i]) < 0.33 and float(DataFrame['Z_2'][i]) > 0.26:
    #             BPT, WHAN, LAGN, LAGN_er = spectra_processing(DataFrame['OIIIR_FLUX'][i], DataFrame['OIIIR_FLUX_ERR'][i], DataFrame['HB_FLUX'][i], DataFrame['HB_FLUX_ERR'][i],
    #                                                       DataFrame['NIIR_FLUX'][i], DataFrame['NIIR_FLUX_ERR'][i], DataFrame['HA_FLUX'][i], DataFrame['HA_FLUX_ERR'][i],
    #                                                       DataFrame['HA_EW'][i], DataFrame['HA_EW_ERR'][i], DataFrame['Z'][i])
    #             BPTs.append(BPT)
    #             WHANs.append(WHAN)
    #         else:
    #             DataFrame = DataFrame.drop(i)
            
    #     DataFrame['BPT'] = BPTs
    #     DataFrame['WHAN'] = WHANs
        
    #     DataFrame = DataFrame.drop(DataFrame[DataFrame['BPT'] == 'NDA'].index)
    #     DataFrame = DataFrame.drop(DataFrame[DataFrame['WHAN'] == 'NDA'].index)
    #     print(DataFrame.shape)
    #     DataFrame.to_csv(self.out_path, index=False)
    
    # def plotter(self):
    #     bids_mass_plt = [[10.0, 10.25], [10.25, 10.5], [10.5, 10.75], [10.75, 11], [11, 11.25], [11.25, 11.5]]
        
    #     plotting({
    #         'input_path' : self.out_path,
    #         'x' : 'mass_stellar_percentile50',
    #         'y' : 'SFR_0_1Gyr_percentile50',
    #         'up' : 'SFR_0_1Gyr_percentile84',
    #         'down' : 'SFR_0_1Gyr_percentile16',
    #         'xlim' : [9.9, 11.6],
    #         'ylim' : [-3.1, 2.3],
    #         'xticks' : np.arange(10.0, 11.6, 0.25),
    #         'yticks' : np.arange(-3, 2.2, 1),
    #         'xlabel' : r'$log(M_s / M_{\odot})$',
    #         'ylabel' : r'$log(SFR) / M_{\odot} yr^{-1})$',
    #         'bids': bids_mass_plt,
    #         'save_path' : r'SFRSM_LTG_033.pdf',
    #         'theor_lines' : 'sfrsm'
    #     })

    def read(self):
        DataFrame = pd.read_csv(self.out_path, sep=',')
        MS_flag = MS_flagging(DataFrame['SFR_0_1Gyr_percentile50'], DataFrame['mass_stellar_percentile50'], DataFrame['Z'])
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
                                                          DataFrame['HA_EW'][i], DataFrame['HA_EW_ERR'][i], DataFrame['Z'][i])
            
            out = outflow_wAGN(DataFrame['SFR_0_1Gyr_percentile50'][i], LAGN, DataFrame['mass_stellar_percentile50'][i])
            OUTFLOW.append(out)
            if out == -99999.0:
                OUTFLOW_up.append(-99999.0)
                OUTFLOW_down.append(-99999.0)
            else:    
                if LAGN_er in [-1, -2]:
                    OUTFLOW_up.append(-100000.0)
                    OUTFLOW_down.append(-100000.0)
                else:
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
            'CATAID' : DataFrame['CATAID'],
            'RA' : DataFrame['RA'],
            'DEC' : DataFrame['DEC'],
            'Z' : DataFrame['Z'],
            'bMS/MS' : MS_flag,
            'BPT' : BPTs,
            'WHAN' : WHANs, 
            'GALINDEX_r' : DataFrame['GALINDEX_r'],
            'GALINDEXERR_r' : DataFrame['GALINDEXERR_r'],
            'outflow_agn_on_percentile16' : OUTFLOW_down,
            'outflow_agn_on_percentile50' : OUTFLOW, 
            'outflow_agn_on_percentile84' : OUTFLOW_up, 
            'outflow_agn_off_percentile16' : OUTFLOW_down_off, 
            'outflow_agn_off_percentile50' :  OUTFLOW_off,
            'outflow_agn_off_percentile84' : OUTFLOW_up_off 
        }
        
        NewDataFrame = pd.DataFrame(NewFrameDict)
        NewDataFrame = NewDataFrame.round(decimals=4)
        NewDataFrame.to_csv(self.final_out_path, index=False)
        
if __name__ == '__main__':
    obj = Main(r'E:/LICENSE/ProgsData/main/GAMAv3.txt')
    obj.completely_different_extraction()
    obj.read()