import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u 
from matplotlib import gridspec
import pandas as pd
from astropy.constants import c

def plotter_extractor_sdss(path_to_file, output_path, wave):
    hdul = fits.open(path_to_file)
    data1 = hdul[1].data
 
    wavelength = 10 ** data1.field('loglam') * u.Unit('AA')
    flux = data1.field('flux') * 10**-17 * u.Unit('erg cm-2 s-1 AA-1')
    model = data1.field('model') * 10**-17 * u.Unit('erg cm-2 s-1 AA-1')
    
    fig = plt.figure(figsize=(12,12), tight_layout=True)
    spec = gridspec.GridSpec(ncols=2, nrows=2,
                         width_ratios=[2, 1], wspace=0.1,
                         hspace=0, height_ratios=[1, 1])
    # fig.suptitle(path_to_file.split('S\\')[-1])
    
    ax0 = fig.add_subplot(spec[0])
    ax1 = fig.add_subplot(spec[2])
 
    ax0.plot(wavelength, flux, linewidth=1, c='black')
    ax1.plot(wavelength, model, linewidth=1, c='r')
    
    data3 = hdul[3].data
    
    line_names = data3.field('LINENAME')
    line_waves = data3.field('LINEWAVE')
    line_z = data3.field('LINEZ')
    line_area = data3.field('LINEAREA')
 
    df = pd.DataFrame(
    {'name': [i for i in line_names],
     'lam_rest': [i for i in line_waves],
     'z': [i for i in line_z],
     'area': [i for i in line_area]
     }
    )
 
    df = df[abs(df['area']) != 0]
    df.reset_index(inplace=True)
    df['lam_obs'] = df['lam_rest'] * (1 + df['z'])
    label_y = np.random.uniform(
    low=model.min().value*1.1,
    high=model.max().value*0.9,
    size=len(df)
    )
    for ax in [ax0, ax1]:
        for i in range(len(df)):
            ax.axvline(x=df['lam_obs'].iloc[i],
               color='r',
               alpha=0.3,
               label=df['name'].iloc[i],
               ls='--',
               lw=0.7)
     
            ax.text(x=df['lam_obs'].iloc[i] ,
                y=label_y[i],
                s=df['name'].iloc[i],
                fontsize='small',
                rotation=90)
 
    ax1.set_xlabel('Wavelength')
    ax0.set_ylabel('Flux (raw)')
    ax1.set_ylabel('Flux (fit)')
    
    Halpha_observed = wave * (1 + df['z'][0]) 
    wavelength = 10 ** data1.field('loglam')
    res_min = np.abs(wavelength - Halpha_observed - 200)
    res_max = np.abs(wavelength - Halpha_observed + 200)
    i_max = np.where(res_min == np.min(res_min))[0][0]
    i_min = np.where(res_max == np.min(res_max))[0][0]
    
    x = wavelength
    vel = (x[i_min:i_max] - Halpha_observed)*c.to('km/s').value/x[i_min:i_max]
    Y_DATA_FLUX = flux[i_min:i_max]
    Y_DATA_MODEL = model[i_min:i_max]
    
    ax2 = fig.add_subplot(spec[1])
    ax2.plot(vel, Y_DATA_FLUX, color='black')
    ax3 = fig.add_subplot(spec[3])
    ax3.plot(vel, Y_DATA_MODEL, color='r')
    ax3.set_xlabel('Velocity, km/s')
    
    fig.savefig(output_path, dpi=300, bbox_inches = 'tight', pad_inches = 0.0001)
    
def plotter_extractor_gama(path_to_file, output_path, wave):
    hdul = fits.open(path_to_file)
    hdu = hdul[0]
    flux = hdu.data[0,:]
    wavelength = np.linspace(float(hdu.header['WMIN']), float(hdu.header['WMAX']), int(hdu.header['NAXIS1']))
    
    fig = plt.figure(figsize=(12,6), tight_layout=True)
    spec = gridspec.GridSpec(ncols=2, nrows=1,
                         width_ratios=[2, 1], wspace=0.1,
                         hspace=0)
    
    ax0 = fig.add_subplot(spec[0])
    ax0.plot(wavelength, flux, linewidth=1, c='black')
    
    ax0.set_xlabel('Wavelength')
    ax0.set_ylabel('Flux (raw)')
    
    Halpha_observed = wave * (1 + float(hdu.header['Z'])) 
    res_min = np.abs(wavelength - Halpha_observed - 200)
    res_max = np.abs(wavelength - Halpha_observed + 200)
    i_max = np.where(res_min == np.min(res_min))[0][0]
    i_min = np.where(res_max == np.min(res_max))[0][0]
    
    x = wavelength
    vel = (x[i_min:i_max] - Halpha_observed)*c.to('km/s').value/x[i_min:i_max]
    Y_DATA_FLUX = flux[i_min:i_max]
    
    ax2 = fig.add_subplot(spec[1])
    ax2.plot(vel, Y_DATA_FLUX, color='black')
    ax2.set_xlabel('Velocity, km/s')
    fig.savefig(output_path, dpi=300, bbox_inches = 'tight', pad_inches = 0.0001)