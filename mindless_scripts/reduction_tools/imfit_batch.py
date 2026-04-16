from casatasks import imfit, imstat
from pathlib import Path
import numpy as np
from astropy.coordinates import SkyCoord
import click
import csv
import sys
from io import StringIO
import subprocess
import re

def extract_fits_header(fits_path: str) -> dict:
    """
    Extract specific header values from a FITS file using the fitshdr command.
    
    Args:
        fits_path: Path to the FITS file
        
    Returns:
        Dictionary with keys: DATE-OBS, PROJECT, CRVAL3, SBID, DURATION
    """
    keys = ["DATE-OBS", "PROJECT", "CRVAL3", "SBID", "DURATION"]
    
    try:
        result = subprocess.run(
            ["fitshdr", fits_path],
            capture_output=True, text=True, check=True
        )
        header_text = result.stdout
    except FileNotFoundError:
        # Fallback to fitsheader (astropy) if fitshdr not available
        result = subprocess.run(
            ["fitsheader", fits_path],
            capture_output=True, text=True, check=True
        )
        header_text = result.stdout
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Failed to read FITS header: {e.stderr}") from e

    values = {}
    for key in keys:
        pattern = rf"^{key}\s*=\s*([^/\n]+)"
        match = re.search(pattern, header_text, re.MULTILINE)
        if match:
            raw = match.group(1).strip().strip("'").strip()
            # Try numeric conversion
            try:
                values[key] = float(raw) if '.' in raw or 'e' in raw.lower() else int(raw)
            except ValueError:
                values[key] = raw
        else:
            values[key] = None

    # Apply transformations
    if values.get("CRVAL3") is not None:
        values["CRVAL3"] = values["CRVAL3"] / 1e9

    if values.get("DURATION") is not None:
        values["DURATION"] = values["DURATION"] / 60

    return values

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

def flux_extractor_to_dict(path, coords, stokes, radius, rms_radius, scale):
    """
    Extract flux and positional fit statistics from a .fits image and return them as a dictionary.

    This function runs a Gaussian fit (via CASA's `imfit`) on a specified circular
    region around given sky coordinates in an image, and computes associated flux
    and positional uncertainties. It also estimates the local RMS noise from an
    annular region around the source (via CASA's `imstat`), and combines this with
    the fitting error and a user-provided fractional calibration uncertainty.

    Parameters
    ----------
    path : str
        Path to the CASA image (e.g., a FITS or CASA image cube) to analyze.
    coords : astropy.coordinates.SkyCoord
        Sky coordinates of the target position used to center the fit and RMS
        regions. Must be in an equatorial frame with accessible `.ra` and `.dec`.
    stokes : str
        Stokes parameter of the image to analyze (e.g., 'I', 'Q', 'U', 'V').
    radius : float
        Radius of the circular fit region in arcseconds, centered on `coords`.
    rms_radius : float
        Outer radius of the annulus in arcseconds used to compute the RMS noise.
        The inner radius of the annulus is set to `radius`.
    scale : float
        Fractional calibration (or systematic) uncertainty to be applied to the
        peak flux. The corresponding absolute error term is `scale * peak_flux`.

    Returns
    -------
    dict
        A dictionary containing:
        
        - 'file' : str
            The input image path.
        - 'ra_input' : float
            Input right ascension of `coords` in degrees.
        - 'dec_input' : float
            Input declination of `coords` in degrees.
        - 'rms_mJy' : float or None
            RMS noise in the annulus region, in mJy.
        - 'peak_flux_mJy' : float or None
            Fitted peak flux density from `imfit`, in mJy.
        - 'peak_err_mJy' : float or None
            Formal peak flux uncertainty from `imfit`, in mJy.
        - 'combined_err_mJy' : float or None
            Quadrature combination of fit error, RMS noise, and `scale * peak_flux`,
            in mJy.
        - 'ra_fit_deg' : float or None
            Fitted right ascension in degrees.
        - 'ra_err_deg' : float or None
            Uncertainty on the fitted right ascension in degrees.
        - 'dec_fit_deg' : float or None
            Fitted declination in degrees.
        - 'dec_err_deg' : float or None
            Uncertainty on the fitted declination in degrees.
        - 'fit_failed' : bool
            True if `imfit` failed or did not converge, False otherwise.

    Notes
    -----
    This function assumes that CASA tasks `imfit` and `imstat` are available in the
    execution environment. If the fit fails or does not converge, flux- and
    position-related fields in the returned dictionary will remain `None`, and
    `fit_failed` will be set to True, but the RMS estimate will still be returned
    if possible.
    """

    result = {
        'file': path,
        'ra_input': coords.ra.deg,
        'dec_input': coords.dec.deg,
        'rms_mJy': None,
        'peak_flux_mJy': None,
        'peak_err_mJy': None,
        'combined_err_mJy': None,
        'ra_fit_deg': None,
        'ra_err_deg': None,
        'dec_fit_deg': None,
        'dec_err_deg': None,
        'fit_failed': False,
    }

    hmsdms = coords.to_string('hmsdms', sep=':', precision=1, pad=True)
    ra_str, dec_str = hmsdms.split()
    hh, mm, ss = ra_str.split(':')
    ra_casa = f"{hh}h{mm}m{ss}s"
    sign = '+' if dec_str.strip()[0] == '+' else '-'
    dd, dm, ds = dec_str.replace('+','').replace('-','').split(':')
    dec_casa = f"{sign}{dd}.{dm}.{ds[:4]}"
    region_temp = f"{ra_casa},{dec_casa}"
    region = 'circle[[{}], {}arcsec]'.format(region_temp, radius)
    annulus_region = 'annulus[[{}], [{}arcsec, {}arcsec]]'.format(region_temp, radius, rms_radius)

    source_name = "J"+region_temp[0:2]+region_temp[3:5]
    name = f"{source_name}_{stokes}_{radius}_arcsec"
    summary = 'fit_summary_{}'.format(name)
    logfile = 'fit_log_{}'.format(name)

    try:
        fit = imfit(
            imagename=path,
            region=region,
            chans='',
            stokes="",
            summary=summary,
            logfile=logfile,
        )
        fit_failed = not fit.get('converged', [False])[0]
    except Exception as e:
        fit = None
        fit_failed = True

    stats = imstat(imagename=path, stokes=stokes, region=annulus_region)
    rms = stats["rms"] * 10**3
    result['rms_mJy'] = round(rms[0], 4)
    result['fit_failed'] = fit_failed

    if not fit_failed and fit is not None:
        temp = fit['results']['component0']['peak']
        peak_flux = temp['value'] * 10**3
        peak_err = temp['error'] * 10**3
        combined_err = (peak_err**2 + rms**2 + (scale * peak_flux)**2)**(0.5)

        fit_coords = fit['results']['component0']['shape']['direction']
        ra = fit_coords['m0']['value'] * 180 / np.pi % 360
        ra_err = fit_coords['error']['longitude']['value'] / 3600
        dec = fit_coords['m1']['value'] * 180 / np.pi
        dec_err = fit_coords['error']['latitude']['value'] / 3600

        result.update({
            'peak_flux_mJy': round(peak_flux, 4),
            'peak_err_mJy': round(peak_err, 4),
            'combined_err_mJy': round(combined_err[0], 4),
            'ra_fit_deg': round(ra, 6),
            'ra_err_deg': round(ra_err, 6),
            'dec_fit_deg': round(dec, 6),
            'dec_err_deg': round(dec_err, 6),
        })

    return result


@click.group()
def cli():
    pass


@cli.command()
@click.option("-g", "--galactic", is_flag=True, default=False, help="If set, assumes galactic coordinates")
@click.argument("path", nargs=1, type=str)
@click.argument("coords", nargs=1, type=str)
@click.argument("stokes", nargs=1, type=str, default="I")
@click.argument("radius", nargs=1, type=float, default=5)
@click.argument("rms_radius", nargs=1, type=float, default=30)
@click.argument("scale", nargs=1, type=float, default=0.1)
def single(path, coords, stokes, radius, rms_radius, scale, galactic):
    pos_eq, __ = unicoord(coords, galactic, display=False)
    flux_extractor(path, pos_eq, stokes, radius, rms_radius, scale)


@cli.command()
@click.option("-g", "--galactic", is_flag=True, default=False, help="If set, assumes galactic coordinates")
@click.option("-o", "--output", default="results.csv", show_default=True, help="Output CSV file")
@click.argument("paths", nargs=-1, type=str, required=True)
@click.argument("coords", nargs=1, type=str)
@click.argument("stokes", nargs=1, type=str, default="I")
@click.argument("radius", nargs=1, type=float, default=5)
@click.argument("rms_radius", nargs=1, type=float, default=30)
@click.argument("scale", nargs=1, type=float, default=0.1)
def batch(paths, coords, stokes, radius, rms_radius, scale, galactic, output):
    pos_eq, __ = unicoord(coords, galactic, display=False)
    rows = []
    for path in paths:
        print(f"Processing {path}...")
        row = flux_extractor_to_dict(path, pos_eq, stokes, radius, rms_radius, scale)
        rows.append(row)

    with open(output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nSaved {len(rows)} rows to {output}")


if __name__ == "__main__":
    cli()