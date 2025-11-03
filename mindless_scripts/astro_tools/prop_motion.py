

import astropy.units as u
from astropy.coordinates import SkyCoord, Distance
from astropy.time import Time
import warnings
import click

from mindless_scripts.astro_tools.unicoord import unicoord

def prop_motion(coords, pm_ra, pm_dec, old_time, new_time, distance=None, degrees=False):
    """
    Propagate coordinates to a new time given proper motion and parallax.

    Parameters:
    coords : SkyCoord
        Initial coordinates.
    pm_ra : Quantity
        Proper motion in right ascension (with cos(dec) factor).
    pm_dec : Quantity
        Proper motion in declination.
    distance : Quantity, optional
        Distance to the object in parsecs.
    old_time : Time
        Original observation time.
    new_time : Time
        New time to propagate to.
    Returns:
    SkyCoord
        New coordinates at the specified time.
    """
    warnings.filterwarnings("ignore", category=UserWarning)
    pos_old= SkyCoord(
        ra=coords.ra.deg*u.deg,
        dec=coords.dec.deg*u.deg,
        frame='icrs',
        distance=distance,
        pm_ra_cosdec=pm_ra*u.mas/u.yr,
        pm_dec=pm_dec*u.mas/u.yr,
        obstime=Time(old_time))
    t_new = Time(new_time) 
    pos_new = pos_old.apply_space_motion(t_new)
    separation = pos_old.separation(pos_new).arcsecond
    print(f"Separation: {separation:.3f} arcsec")  

    if degrees is True:
        print(f"Old Position ({pos_old.obstime.isot.split('.')[0]}): RA: {pos_old.ra.deg:.6f} deg, DEC: {pos_old.dec.deg:.6f} deg")
        print(f"New Position ({new_time}): RA: {pos_new.ra.deg:.6f} deg, DEC: {pos_new.dec.deg:.6f} deg")
    else:
        print(f"Old Position ({pos_old.obstime.isot.split('.')[0]}): RA: {pos_old.ra.to_string(unit=u.hour, sep=':', pad=True)}, DEC: {pos_old.dec.to_string(unit=u.deg, sep=':', pad=True, alwayssign=True)}")
        print(f"New Position ({new_time}): RA: {pos_new.ra.to_string(unit=u.hour, sep=':', pad=True)}, DEC: {pos_new.dec.to_string(unit=u.deg, sep=':', pad=True, alwayssign=True)}")
    return pos_new

@click.command()
@click.option("-p","--parallax", is_flag=True, show_default=True,default=False, help="If set, calculate distance from parallax")
@click.option("-g","--galactic", is_flag=True, show_default=True,default=False, help="If set, assumes galactic coordinates")
@click.option("-d","--degrees", is_flag=True, show_default=True,default=False, help="If set, prints coordinates in degrees")
@click.argument("coords", nargs=1, type=str)
@click.argument("pm_ra", nargs=1, type=float, default=0)
@click.argument("pm_dec", nargs=1, type=float, default=0)
@click.argument("old_time", nargs=1, type=str, default='J2000')
@click.argument("new_time", nargs=1, type=str, default='2025-01-01T12:00:00')
@click.argument("distance", nargs=1, type=float, default=None)
def main(coords, pm_ra, pm_dec, old_time, new_time, distance, parallax,galactic,degrees):
    pos_eq, __ = unicoord(coords,galactic,display=False)
    if parallax is True:
        distance=Distance(parallax = distance*u.mas)
    else:
        distance=Distance(distance*u.pc)
    prop_motion(pos_eq,pm_ra, pm_dec, old_time, new_time, distance,degrees )
    return

if __name__ == "__main__":
	main()