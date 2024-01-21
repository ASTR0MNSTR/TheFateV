import math

def log_er(item):
    return [[abs(math.log(1-item, 10))], [abs(math.log(1+item, 10))]]

def atten_coef(A_V):
    R_V = 4.05
    OIII_c = 2.659*(-2.156 + (1.509/0.5007) -
                        (0.198/(0.5007**2)) + (0.011/(0.5007**3))) + R_V
    NII_c = 2.659*(-1.857 + (1.040/0.6583)) + R_V
    SII_c = 2.659*(-1.857 + (1.040/0.6732)) + R_V
    OI_c = 2.659*(-1.857 + (1.040/0.6300)) + R_V
    HB_c = 2.659*(-2.156 + (1.509/0.4861) -
                      (0.198/(0.4861**2)) + (0.011/(0.4861**3))) + R_V
    #HD_c = 2.659*(-2.156 + (1.509/0.4101734) -
    #                  (0.198/(0.4101734**2)) + (0.011/(0.4101734**3))) + R_V
    #    HG_c = 2.659*(-2.156 + (1.509/0.4340472) -
    #                  (0.198/(0.4340472**2)) + (0.011/(0.4340472**3))) + R_V
    OII_c = 2.659*(-2.156 + (1.509/0.3727) -
                      (0.198/(0.3727**2)) + (0.011/(0.3727**3))) + R_V
    res = [OIII_c, NII_c, SII_c, OI_c, HB_c, OII_c]
    res_out = [item*A_V/R_V for item in res]
    res_out_out = [10**(0.4*item) for item in res_out]
    return res_out_out

def dust_correction(HA, HA_er, HB, HB_er):
    SN = 2
    if HA_er > 0 and HB_er > 0 and HA != -99999.0 and HB != -99999.0 and HA > SN*HA_er and HB > SN*HB_er and HA/HB > 2.86:
        f_HA = ((HA/HB)/2.86)**(2.114)
        A_V = 2.5*math.log(f_HA, 10)
        coefs = atten_coef(A_V)
        coefs.append(f_HA)
    else:
        coefs = [1, 1, 1, 1, 1, 1, 1]
    return coefs

def line_flagging(pair, spec_coefs):
    SN = 2
    abs = 0
    pair_flags = []
    flags_dict = {
        0 : 'down',
        1 : 'up'
    }
    m = 0
    for j, item in enumerate(pair):
        if item[0] == -99999.0 or item[1] < 0:
            return -100, -100, [-99, -99], 0
        elif item[0] < 0 and item[1] > 0 and SN*item[1] + item[0] < 0: #abs
            pair_flags.append(flags_dict[j])
            m += 1
            item[0] = SN*item[1]*spec_coefs[j]
            abs += 1
        elif item[0] <= SN*item[1]: #non-det
            pair_flags.append(flags_dict[j])
            m += 1
            #item[0] += SN*item[1]*spec_coefs[j]
            item[0] += SN*item[1]
            item[0] *= spec_coefs[j]
            item[1] *= spec_coefs[j]
        elif item[0] > SN*item[1]: #det
            item[0] *= spec_coefs[j]
        else:
            print("SMTH WRONG")
    
    if m == 2:
        return -99, -99, [-99, -99], 0
    else:
        res = math.log(pair[0][0]/pair[1][0], 10)

    if m == 0:
        res_er = log_er((pair[0][0]*pair[1][1] + pair[0][1]*pair[1][0])/(pair[1][0]*pair[0][0]))
    else:
        res_er = -99

    return res, res_er, pair_flags, abs

def double_line(res, res_er, pair_flags, limit_SF, limit_AGN, flag):
    AGN = 'PROB'
    if res_er != -99:
        if res - res_er[0][0] > limit_AGN:
            AGN = 'AGN'
        elif res + res_er[1][0] < limit_SF:
            AGN = 'SF'
        else:
            AGN = 'UNC'
    else:
        if pair_flags[0] == 'down' and res < limit_SF:
            AGN = 'SF'
        elif pair_flags[0] == 'up' and res > limit_AGN:
            AGN = 'AGN'
        else:
            AGN = 'UNC'
    AGN += flag
    return AGN

def def_four_lines(X, X_er, Y, Y_er):
    if Y + Y_er[1][0] < (0.61 / ((X + X_er[1][0]) - 0.05)) + 1.3 and X < 0.05:
        AGN = 'SFXY'
    elif Y - Y_er[0][0] > (0.61 / ((X - X_er[0][0]) - 0.47)) + 1.19 and X < 0.47:
        AGN = 'AGNXY'
    elif X > 0.47:
        AGN = 'AGNXY'
    else:
        AGN = 'UNCXY'
    
    return AGN

def four_lines(X, pair_x_flags, Y, pair_y_flags):

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
     
    if X > 0.47 and Y < 1.19:
        if pair_x_flag == 'left':
            AGN = 'UNCXY'
        else:
            AGN = 'AGNXY'
    elif X < 0.47 and Y > 1.3:
        if pair_y_flag == 'down':
            AGN = 'UNCXY'
        else:
            AGN = 'AGNXY' 
    elif X > 0.47 and Y > 1.19:
        if pair_y_flag == 'down' and pair_x_flag == 'left':
            AGN = 'UNCXY'
        else:
            AGN = 'AGNXY'
    else:
        if Y > (0.61 / (X - 0.47)) + 1.19:
            if pair_y_flag == 'down' or pair_x_flag == 'left':
                AGN = 'UNCXY'
            else:
                AGN = 'AGNXY'
        
        elif Y < (0.61 / (X - 0.05)) + 1.3:
            if pair_y_flag == 'up' or pair_x_flag == 'right':
                AGN = 'UNCXY'
            else:
                AGN = 'SFXY'
                
        else:
            AGN = 'UNCXY'

    return AGN

def metallicity(X, X_er, pair_x_flags, Y, Y_er, pair_y_flags):
    #met = 8.73 - 0.32*(Y - X)
    met = 9.37 + 2.03*X + 1.26*(X**2) + 0.32*(X**3)
        #err_plus = 0.32*fluxes[2][1][0] + 0.32*fluxes[3][1][0]
        #err_minus = 0.32*fluxes[2][0][0] + 0.32*fluxes[3][0][0]
    return met, 0, 0

def WHAN(X, X_er, pair_x_flags, HA_ew, HA_ew_err, pair_HA):
    if HA_ew == -99999.0:
        return 'NDA'
    else:
        if HA_ew <= -0.307:
            return 'LLR'
        elif HA_ew < 0.47712 and len(pair_HA) == 0:
            return 'ELR'
        elif HA_ew < 0.47712 and len(pair_HA) != 0:
            return 'RG'
        else:
            if X == -99:
                return 'LLR'
            elif X < -0.4 and 'up' not in pair_x_flags:
                return 'SF'
            elif HA_ew < 0.77815125 and len(pair_x_flags) == 0:
                return 'wAGN'
            elif HA_ew >= 0.77815125 and len(pair_x_flags) == 0:
                return 'sAGN'
            else:
                return 'UNC'
            
def AGN_reg(OIII, OIII_er, HB, HB_er, NII, NII_er, HA, HA_er, HA_ew, HA_ew_err, pair_HA):

    #for dust correction
    coefs = dust_correction(HA, HA_er, HB, HB_er)

    #coefs = [1, 1, 1, 1, 1, 1, 1]
    #coefs = [OIII_c, NII_c, SII_c, OI_c, HB_c, OII_c, HA_c]
    
    Y, Y_er, pair_y_flags, abs_Y = line_flagging([[OIII, OIII_er], [HB, HB_er]], [coefs[0], coefs[4]])
    X, X_er, pair_x_flags, abs_X = line_flagging([[NII, NII_er], [HA, HA_er]], [coefs[1], coefs[6]])

    if X == -100 or Y == -100:
        AGN = 'NDA'
        SC_WHAN = 'NDA'
    elif X == -99 and Y == -99:
        AGN = 'NOEL'
        SC_WHAN = WHAN(X, X_er, pair_x_flags, HA_ew, HA_ew_err, pair_HA)
    elif X == -99 and Y != -99:
        AGN = double_line(Y, Y_er, pair_y_flags, 1.3, 1.19, 'Y')
        SC_WHAN = WHAN(X, X_er, pair_x_flags, HA_ew, HA_ew_err, pair_HA)
    elif X != 99 and Y == -99:
        AGN = double_line(X, X_er, pair_x_flags, 0.05, 0.47, 'X')
        SC_WHAN = WHAN(X, X_er, pair_x_flags, HA_ew, HA_ew_err, pair_HA)
    elif len(pair_x_flags) == 0 and len(pair_y_flags) == 0:
        AGN = def_four_lines(X, X_er, Y, Y_er)
        SC_WHAN = WHAN(X, X_er, pair_x_flags, HA_ew, HA_ew_err, pair_HA)
    else:
        SC_WHAN = WHAN(X, X_er, pair_x_flags, HA_ew, HA_ew_err, pair_HA)
        AGN = four_lines(X, pair_x_flags, Y, pair_y_flags)
    
    #if abs_X > 0 or abs_Y > 0:
    #    AGN += '!'
    #    SC_WHAN += '!'

    return AGN, X, pair_x_flags, Y, pair_y_flags, SC_WHAN
