import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u 
from matplotlib import gridspec
import pandas as pd
from astropy.constants import c

def extract_sdss(path_to_file):
    hdul = fits.open(path_to_file)
    hdu = hdul[0]
    flux = hdu.data[0,:]/0.85
    i = np.linspace(int(hdu.header['CRPIX1']), int(hdu.header['NAXIS1']), int(hdu.header['NAXIS1']))
    wavelength = 10**(float(hdu.header['CRVAL1']) + float(hdu.header['CDELT1'])*i)
    return wavelength, flux

def extract_gama(path_to_file):
    hdul = fits.open(path_to_file)
    hdu = hdul[0]
    flux = hdu.data[0,:]
    wavelength = np.linspace(float(hdu.header['WMIN']), float(hdu.header['WMAX']), int(hdu.header['NAXIS1']))
    return wavelength, flux

def gaussian(x, c1, mu1, sigma1):
    res =  c1 * np.exp( - (x - mu1)**2.0 / (2.0 * sigma1**2.0) )
    return res

def continuum(x, a, b, wave):
    res = a*wave*(x/c.to('km/s').value + 1) + b
    return res

# def multi_gaussian(x, *pars):
#     gs = 0
#     for i in range(int(len(pars) / 3)):
#         gs += gaussian(x, pars[3*i], pars[3*i+1], pars[3*i+2])
#     return gs

def velocity_profile_plotter(ax, flux, wavelength, line_observed, Z, pars_cont, pars_gauss):
    res_min = np.abs(wavelength - line_observed - 300)
    res_max = np.abs(wavelength - line_observed + 300)
    i_max = np.where(res_min == np.min(res_min))[0][0]
    i_min = np.where(res_max == np.min(res_max))[0][0]
    x = wavelength
    vel = (x[i_min:i_max] - line_observed)*c.to('km/s').value/x[i_min:i_max]
    Y_DATA_FLUX = flux[i_min:i_max]
    ax.plot(vel, Y_DATA_FLUX, linewidth=1, c='black')
    # ax.plot(vel, continuum(vel, pars_cont[0], pars_cont[1], line_observed), color='red', linestyle='dashed')
    model_gaussian_top = continuum(vel, pars_cont[0], pars_cont[1], line_observed)
    model_gaussian_bottom = continuum(vel, pars_cont[0], pars_cont[1], line_observed)
    for item in pars_gauss:
        if -99999.0 not in item:
            # print('FOUND FUCKER')
            color = item[-1]
            cont = continuum(vel, pars_cont[0], pars_cont[1], line_observed)
            gaussian_top = gaussian(vel, (item[0]+item[1]), item[4]*c.to('km/s').value/line_observed, (item[2]+item[3])*c.to('km/s').value/(line_observed))
            gaussian_bottom = gaussian(vel, (item[0]-item[1]), item[4]*c.to('km/s').value/line_observed, (item[2]-item[3])*c.to('km/s').value/(line_observed))
            ax.plot(vel, continuum(vel, pars_cont[0], pars_cont[1], line_observed) + gaussian(vel, item[0], item[4]*c.to('km/s').value/line_observed, item[2]*c.to('km/s').value/(line_observed)), color=color)
            ax.fill_between(vel, cont+gaussian_top, cont+gaussian_bottom, color=color)
            model_gaussian_top += gaussian_top
            model_gaussian_bottom += gaussian_bottom
    
    ax.fill_between(vel, model_gaussian_top, model_gaussian_bottom, color='r')
    ax.set_xlabel('Velocity, km/s')
    ax.set_xlim(-5500, 5500)
    # ax.set_ylabel('Flux, 10^-17 erg cm-2 s-1 AA-1')
    return ax

def universal_plotter(path_to_file, path_to_save, data_row):
    if data_row['SPECID'][0][0] == 'G':
        wavelength, flux = extract_gama(path_to_file)
    else:
        wavelength, flux = extract_sdss(path_to_file)
    Z = data_row['Z'][0]
        
    fig = plt.figure(figsize=(16,6), tight_layout=True)
    spec = gridspec.GridSpec(ncols=3, nrows=1,
                         width_ratios=[2, 1, 1], wspace=0.15,
                         hspace=0)
    
    ax_obs = fig.add_subplot(spec[0])
    ax_obs.set_title(data_row['SPECID'][0])
    ax_HB = fig.add_subplot(spec[1])
    ax_HA = fig.add_subplot(spec[2])
    
    ax_obs.plot(wavelength, flux, linewidth=1, c='black')
    ax_obs.set_xlabel('Wavelength, AA')
    ax_obs.set_ylabel('Flux, 10^-17 erg cm-2 s-1 AA-1')
    
    pars_cont_HB = [data_row['HB_GRAD'][0], data_row['HB_CONT'][0]]
    pars_gauss_HB = [
        [data_row['AMP_HB'][0], data_row['AMP_HB_ERR'][0], data_row['SIG_HB'][0], data_row['SIG_HB_ERR'][0], data_row['POS_HB'][0], data_row['POS_HB_ERR'][0], 'lime'],
        [data_row['AMP_HB_B'][0], data_row['AMP_HB_B_ERR'][0], data_row['SIG_HB_B'][0], data_row['SIG_HB_B_ERR'][0], data_row['POS_HB_B'][0], data_row['POS_HB_B_ERR'][0], 'darkorange'],
        [data_row['AMP_OIIIR'][0], data_row['AMP_OIIIR_ERR'][0], data_row['SIG_OIIIR'][0], data_row['SIG_OIIIR_ERR'][0], data_row['POS_OIIIR'][0], data_row['POS_OIIIR_ERR'][0], 'blue'],
        [data_row['AMP_OIIIB'][0], data_row['AMP_OIIIB_ERR'][0], data_row['SIG_OIIIB'][0], data_row['SIG_OIIIB_ERR'][0], data_row['POS_OIIIB'][0], data_row['POS_OIIIB_ERR'][0], 'violet']
    ]

    ax_HB = velocity_profile_plotter(ax_HB, flux, wavelength, data_row['HB_CEN'][0], Z, pars_cont_HB, pars_gauss_HB)
    ax_HB.set_title(r'$\mathrm{H \beta}$')
    
    pars_cont_HA = [data_row['HA_GRAD'][0], data_row['HA_CONT'][0]]
    pars_gauss_HA = [
        [data_row['AMP_HA'][0], data_row['AMP_HA_ERR'][0], data_row['SIG_HA'][0], data_row['SIG_HA_ERR'][0], data_row['POS_HA'][0], data_row['POS_HA_ERR'][0], 'lime'],
        [data_row['AMP_HA_B'][0], data_row['AMP_HA_B_ERR'][0], data_row['SIG_HA_B'][0], data_row['SIG_HA_B_ERR'][0], data_row['POS_HA_B'][0], data_row['POS_HA_B_ERR'][0], 'darkorange'],
        [data_row['AMP_NIIR'][0], data_row['AMP_NIIR_ERR'][0], data_row['SIG_NIIR'][0], data_row['SIG_NIIR_ERR'][0], data_row['POS_NIIR'][0], data_row['POS_NIIR_ERR'][0], 'blue'],
        [data_row['AMP_NIIB'][0], data_row['AMP_NIIB_ERR'][0], data_row['SIG_NIIB'][0], data_row['SIG_NIIB_ERR'][0], data_row['POS_NIIB'][0], data_row['POS_NIIB_ERR'][0], 'violet']
    ]
      
    ax_HA = velocity_profile_plotter(ax_HA, flux, wavelength, data_row['HA_CEN'][0], Z, pars_cont_HA, pars_gauss_HA)
    ax_HA.set_title(r'$\mathrm{H \alpha}$')
    
    # fig.suptitle(f'SPECID = {specid} \n CATAID = {cataid}')
    
    fig.savefig(path_to_save, dpi=70, bbox_inches = 'tight', pad_inches = 0.0001)