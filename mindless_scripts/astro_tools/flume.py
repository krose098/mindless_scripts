from astropy import units as u
import numpy as np
import click
from astropy.coordinates import Distance
from astropy.cosmology import Planck18, WMAP9

def redshift_to_distance(z,cosmology=WMAP9):
    """
    Convert redshift to distance using astropy.
    
    Parameters:
    z (float): Redshift
    cosmology (astropy.cosmology): Cosmology type for Hubble constant -Default: WMAP9
    
    Returns:
    distance (astropy.units.quantity.Quantity): Distance in Mpc
    """
    d = Distance(unit=u.pc, z = z,cosmology=cosmology).to(u.Mpc)
    return d.value

def flume(flux,dist,bandwidth,galactic,specific):
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
    if galactic is False:
        distance = dist * u.Mpc
    else:
        distance = dist * u.pc
    
    luminosity_per_Hz = (4 * np.pi * distance**2 * flux).to(u.erg / u.s / u.Hz)
    print(f"Luminosity per Hz: {luminosity_per_Hz:.3e}")

    if specific is False:
        total_luminosity = (luminosity_per_Hz * bandwidth.to(u.Hz)).to(u.erg / u.s)
        print(f"Total Luminosity: {total_luminosity:.3e}")
    
    return

@click.command()
@click.option("-g","--galactic", is_flag=True, show_default=True,default=False, help="If set, assumes distance is in pc instead of Mpc")
@click.option("-z","--redshift", is_flag=True, show_default=True,default=False, help="If set, assumes redshift given instead of distance")
@click.option("-s","--specific", is_flag=True, show_default=True,default=False, help="If set, skips total luminosity calculation")
@click.argument("flux", type=float)
@click.argument("dist", type=float)
@click.argument("bandwidth", type=float)
def main(flux,dist,bandwidth,galactic,redshift,specific):
    if redshift is True:
        dist = redshift_to_distance(dist)
    flume(flux,dist,bandwidth,galactic,specific)
    return
if __name__ == "__main__":
	main()
