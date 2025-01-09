

from astropy.coordinates import SkyCoord
import click
 
@click.command()
@click.option("-g","--galactic", is_flag=True, show_default=True,default=False, help="If set, assumes galactic coordinates")
@click.argument("coords", nargs=1, type=str)
def main(coords,galactic):
	ra,dec = coords.split(" ")
	raunit = "hourangle" if ":" in ra or "h" in ra else "deg"
	if galactic is False:
		pos_eq = SkyCoord(ra=ra, dec=dec, unit=(raunit, "deg"))
		pos_gal = pos_eq.galactic
	else:
		pos_gal = SkyCoord(l=ra, b=dec, unit=("deg", "deg"),frame='galactic')
		pos_eq = pos_gal.transform_to('icrs')
	print(29*"=","\nEquatorial Coordinates:")
	print(pos_eq.to_string(style="decimal", precision=6))
	print(pos_eq.to_string(style="hmsdms", precision=3))
	print(29*"=","\nGalactic Coordinates:")
	print(pos_gal.to_string(style="decimal", precision=6))
	print(29*"=")

	return 

if __name__ == "__main__":
	main()


