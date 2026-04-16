from casatasks import imfit, imstat
from pathlib import Path
import numpy as np
from astropy.coordinates import SkyCoord
import click

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

def flux_extractor(path, coords, stokes, radius, rms_radius, scale):
        
    hmsdms = coords.to_string('hmsdms', sep=':', precision=1, pad=True)
    print(hmsdms)
    ra_str, dec_str = hmsdms.split()
    # RA: HH:MM:SS.SS -> HHhMMmSSs.SS
    hh, mm, ss = ra_str.split(':')
    ra_casa = f"{hh}h{mm}m{ss}s"
    # Dec: ±DD:MM:SS.SS -> ±DD.MM.SS.SS
    sign = '+' if dec_str.strip()[0] == '+' else '-'
    dd, dm, ds = dec_str.replace('+','').replace('-','').split(':')
    dec_casa = f"{sign}{dd}.{dm}.{ds}"  
    region_temp = f"{ra_casa},{dec_casa}"
    region = 'circle[[{}], {}arcsec]'.format(region_temp,radius)
    print(region_temp)
    source_name = "J"+region_temp[0:2]+region_temp[3:5]
    name = f"{source_name}_{stokes}_{radius}_arcsec"

    annulus_region = 'annulus[[{}], [{}arcsec, {}arcsec]]'.format(region_temp,radius,rms_radius) 

    imagename=path #TODO save the input parameters and present them as an optional starting point for the next run

    summary='fit_summary_{}'.format(name)
    logfile='fit_log_{}'.format(name)

    try:
        fit = imfit(
            imagename=imagename,
            region=region,
            chans='',
            stokes="",
            summary=summary,
            logfile=logfile,
        )
        fit_failed = not fit.get('converged', [False])[0]
    except Exception as e:
        print(f"imfit failed with error: {e}")
        fit = None
        fit_failed = True

    stats=imstat(imagename=imagename,stokes=stokes, region=annulus_region)

    print('\n'+42*'-')
    print('For {}:'.format(name))

    rms=stats["rms"]*10**3 #rms in mJy

    print("RMS = {} mJy".format(round(rms[0],4)))
    print(42*'-')

    if fit is None or fit_failed:
        print("imfit failed, unable to extract flux and centroid information.")
    else:
        temp=fit['results']['component0']['peak']#['value']
        peak_flux = temp['value']*10**3 #peak flux in mJy
        peak_err = temp['error']*10**3 #peak flux error in mJy

        print('Peak Flux = {} +/- {} mJy'.format(round(peak_flux,4),round(peak_err,4)))
        combined_err = (peak_err**2+rms**2+(scale*peak_flux)**2)**(0.5)
        print("Combined Error = {} mJy".format(round(combined_err[0],4)))

        coords = fit['results']['component0']['shape']['direction']
        ra = coords['m0']['value']*180/np.pi % 360 #conversion from radians to degrees
        ra_err = coords['error']['longitude']['value']/3600 #conversion from arcsec to degrees
        dec = coords['m1']['value']*180/np.pi #conversion from radians to degrees
        dec_err = coords['error']['latitude']['value']/3600 #conversion from arcsec to degrees

        print("Centroid:\nRA = {:.6f} +/- {:.6f} deg \nDec = {:.6f} +/- {:.6f} deg ".format(ra,ra_err,dec,dec_err))
        print(42*'-'+'\n')

@click.command()
@click.option("-g","--galactic", is_flag=True, show_default=True,default=False, help="If set, assumes galactic coordinates")
@click.argument("path", nargs=1, type=str)
@click.argument("coords", nargs=1, type=str)
@click.argument("stokes", nargs=1, type=str, default="I")
@click.argument("radius", nargs=1, type=float, default=5)
@click.argument("rms_radius", nargs=1, type=float, default=30)
@click.argument("scale", nargs=1, type=float, default=0.1)
def main(path, coords ,galactic,stokes, radius, rms_radius, scale):
    pos_eq, __ = unicoord(coords,galactic,display=False)
    flux_extractor(path, pos_eq, stokes, radius, rms_radius, scale)  # noqa: F821
    return

if __name__ == "__main__":
	main()