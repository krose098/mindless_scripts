from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd
import numpy as np

##TODO clickify this code
## TODO import unicoord from reduction_tools and use it here instead of copy-pasting

def unicoord(coords,galactic,display):
	"""
	Universal coordinate converter for equatorial and galactic coordinates.

	Parameters:	
	coords (string): Source coordinates in any standard format
	galactic (bool): If set, assumes galactic coordinates
	display (bool): If set, display all coordinates
	
	Returns:
	ICRS and Galactic coordinates as SkyCoord objects
	"""
	ra,dec = coords.split(" ")
	raunit = "hourangle" if ":" in ra or "h" in ra else "deg"
	if galactic is False:
		pos_eq = SkyCoord(ra=ra, dec=dec, unit=(raunit, "deg"))
		pos_gal = pos_eq.galactic
	else:
		pos_gal = SkyCoord(l=ra, b=dec, unit=("deg", "deg"),frame='galactic')
		pos_eq = pos_gal.transform_to('icrs')
	if display is True:
		print(29*"=","\nEquatorial Coordinates:")
		print(pos_eq.to_string(style="decimal", precision=6))
		print(pos_eq.to_string(style="hmsdms", precision=3))
		print(29*"=","\nGalactic Coordinates:")
		print(pos_gal.to_string(style="decimal", precision=6))
		print(29*"=")
	return pos_eq, pos_gal

def coord2name(ra, dec, catalogue_name: str = "ASKAP"):
    """
    Build IAU-style source names from RA/Dec in degrees.

    Parameters
    ----------
    ra : float, array-like, or pandas.Series
        Right Ascension in degrees.
    dec : float, array-like, or pandas.Series
        Declination in degrees.
    catalogue_name : str, optional
        Prefix catalogue name, default "ASKAP".

    Returns
    -------
    pandas.Series
        Source names like 'ASKAP JHHMMSS.SDDMMSS'.
    """
    # Normalise inputs to 1D numpy arrays
    if np.isscalar(ra):
        ra_arr = np.array([ra], dtype=float)
        dec_arr = np.array([dec], dtype=float)
        scalar_input = True
    else:
        ra_arr = np.asarray(ra, dtype=float)
        dec_arr = np.asarray(dec, dtype=float)
        scalar_input = False

    coords = SkyCoord(ra=ra_arr * u.deg, dec=dec_arr * u.deg, frame="icrs")

    # Use built-in hmsdms string then reformat
    hmsdms = coords.to_string('hmsdms', sep=':', precision=1, pad=True)
    print(hmsdms)
    names = []
    for s in hmsdms:
        ra_str, dec_str = s.split()
        # RA: HH:MM:SS.SS -> HHMMSS.S
        hh, mm, ss = ra_str.split(':')
        ra_compact = f"{hh}{mm}{ss[:4]}"  # truncate to nearst tenth of a second 
        # Dec: ±DD:MM:SS.SS -> ±DDMMSS
        sign = '+' if dec_str.strip()[0] == '+' else '-'
        dd, dm, ds = dec_str.replace('+','').replace('-','').split(':')
        dec_compact = f"{sign}{dd}{dm}{ds[:2]}"  # truncate to nearest second 
        names.append(f"{catalogue_name} J{ra_compact}{dec_compact}")

    result = pd.Series(names)
    if scalar_input:
        return result.iloc[0]
    return result