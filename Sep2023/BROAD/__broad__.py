from astropy.constants import c
from astropy import units as u

def wavelength_to_velocity_corr(waves, constant):
    return (waves*c.to('km/s').value /constant)