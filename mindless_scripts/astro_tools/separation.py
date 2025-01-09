from astropy.coordinates import SkyCoord
import astropy.units as u
import click
 
@click.command()
@click.argument("coords1", nargs=2, type=float)
@click.argument("coords2", nargs=2, type=float)
def main(coords1,coords2):  
    """
    Convert angular separation between two sets of coordinates.
    
    Parameters:
    coords1 (tuple): Tuple of (ra, dec) in degrees
    coords2 (tuple): Tuple of (ra, dec) in degrees
    Returns:
    angular separation in arcseconds
    """
    c1= SkyCoord(ra=coords1[0]*u.deg,dec=coords1[1]*u.deg)
    c2= SkyCoord(ra=coords2[0]*u.deg,dec=coords2[1]*u.deg)
    print(c1.separation(c2).arcsec.round(3),'"')
    return

if __name__ == "__main__":
	main()
