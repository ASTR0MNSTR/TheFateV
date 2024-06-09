import numpy as np

def flux_to_color(s_w1, s_w2):
    return 2.5*np.log10(s_w2/s_w1) #w1 - w2

def flux_to_mag(s_w, zero_mag):
    return 2.5*np.log10(zero_mag/s_w)