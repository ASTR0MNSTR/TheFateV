import math
from Sep2023.AGN_reg_abs import *

def ew_proc(HA_ew, HA_ew_err):
    pair_HA = []
    HA_ew_or = HA_ew
    if HA_ew == -99999.0 or HA_ew_err < 0:
        return -99999.0, -99999.0, pair_HA
    else:
        if HA_ew > 0 and HA_ew < 2*HA_ew_err:
            pair_HA = ['down']
            HA_ew += 2*HA_ew_err
        elif HA_ew < 0 and HA_ew + 2*HA_ew_err <= 0:
            pair_HA = ['down']
            HA_ew = 2*HA_ew_err
        elif HA_ew < 0 and HA_ew + 2*HA_ew_err > 0:
            pair_HA = ['down']
            HA_ew += 2*HA_ew_err
        elif HA_ew > 0 and HA_ew >= 2*HA_ew_err:
            pair_HA = []
        else:
            print("Forgot me!")
        
    return math.log(HA_ew, 10), HA_ew_err, pair_HA


def processor(line):
        line_out = line.strip()
        galaxy_pars = line_out.split()
        try:
            HA = float(galaxy_pars[816])
            HA_er = float(galaxy_pars[815])
            HB = float(galaxy_pars[916])  # should I include errors of HA and HB
            HB_er = float(galaxy_pars[915])
            OIII = float(galaxy_pars[901])
            OIII_er = float(galaxy_pars[900])  # OIIIB, 5007
            NII = float(galaxy_pars[806])
            NII_er = float(galaxy_pars[805])  # NIIB, 6583
            HA_EW = float(galaxy_pars[814])
            HA_EW_ERR = float(galaxy_pars[813])
        except:
            if len(galaxy_pars) != 1042:
                print(len(galaxy_pars))
            return 'NDA', 'NDA'
        
        HA_EW, HA_EW_ERR, pair_HA = ew_proc(HA_EW, HA_EW_ERR)

        AGN, X, X_er, pair_x_flags, Y, Y_er, pair_y_flags, SC_WHAN = AGN_reg(OIII, OIII_er, HB, HB_er, NII, NII_er, HA, HA_er, HA_EW, HA_EW_ERR, pair_HA)

        return AGN, SC_WHAN

def processor_2(line):
    line_out = line.strip()
    galaxy_pars = line_out.split()
    try:
        HA = float(galaxy_pars[816])
        HA_er = float(galaxy_pars[815])
        HB = float(galaxy_pars[916])  # should I include errors of HA and HB
        HB_er = float(galaxy_pars[915])
        OIII = float(galaxy_pars[906])
        OIII_er = float(galaxy_pars[905])  # OIIIB, 5007
        NII = float(galaxy_pars[811])
        NII_er = float(galaxy_pars[810])  # NIIB, 6583
        GAMAID = int(galaxy_pars[735])
        SURVEY = galaxy_pars[740]
        ISBEST = galaxy_pars[742]
        ISSBEST = galaxy_pars[743]
        Z = float(galaxy_pars[738])
    except:
        if len(galaxy_pars) != 1042:
            print(len(galaxy_pars))
        return 0
    SNR = 2
    if HA > SNR*HA_er and HB > SNR*HB_er and OIII > SNR*OIII_er and NII > SNR*NII_er and (SURVEY == 'GAMA' or SURVEY == 'SDSS') and ISBEST == 'true' and ISSBEST == 'true' and Z < 0.3 and HA != -99999.0 and HB != -99999.0 and OIII != -99999.0 and NII != -99999.0:
        return [GAMAID, HA, HB, OIII, NII]
    else:
        return 0
    
