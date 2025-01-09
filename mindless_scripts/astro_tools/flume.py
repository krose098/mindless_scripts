from astropy import units as u
from astropy.constants import c
import numpy as np
import click
 
@click.command()
@click.option("-g","--galactic", is_flag=True, show_default=True,default=False, help="If set, assumes distance is in pc instead of Mpc")
@click.option("-s","--specific", is_flag=True, show_default=True,default=False, help="If set, skips total luminosity calculation")
@click.argument("flux", type=float)
@click.argument("dist", type=float)
@click.argument("bandwidth", type=float)
def main(flux,dist,bandwidth,galactic,specific):
    """
    Convert radio flux to luminosity using astropy.
    
    Parameters:
    flux_mJy (float): Flux in millijanskys (mJy)
    distance_Mpc (float): Distance in megaparsecs (Mpc)
    bandwidth_MHz (float): Observing bandwidth in megahertz (MHz)
    
    Returns:
    prints luminosity per Hz and total luminosity (if specific flag is not set)
    """
    # Convert inputs to astropy quantities
    flux = flux * u.mJy
    bandwidth = bandwidth * u.MHz
    if galactic == False:
        distance = dist * u.Mpc
    else:
        distance = dist * u.pc
    # Convert flux to CGS units (erg/s/cm^2/Hz)
    flux_cgs = flux.to(u.erg / (u.s * u.cm**2 * u.Hz))
    
    # Convert distance to cm
    distance_cm = distance.to(u.cm)
    
    luminosity_per_Hz = (4 * np.pi * distance_cm**2 * flux_cgs).to(u.erg / u.s / u.Hz)
    print(f"Luminosity per Hz: {luminosity_per_Hz:.3e}")

    if specific == False:
        total_luminosity = (luminosity_per_Hz * bandwidth.to(u.Hz)).to(u.erg / u.s)
        print(f"Total Luminosity: {total_luminosity:.3e}")
    
    return

if __name__ == "__main__":
	main()
