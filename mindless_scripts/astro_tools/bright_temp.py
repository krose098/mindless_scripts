from astropy import units as u
import astropy.constants as const
import numpy as np
import click
#TODO Add uncertainty calculation

@click.command()
@click.option("-j","--jupiter", is_flag=True, show_default=True,default=False, help="If set, skips total luminosity calculation")
@click.argument("flux", type=float)
@click.argument("dist", type=float)
@click.argument("freq", type=float)
@click.argument("radius", type=float)
def main(flux, dist, freq, radius, jupiter):

    flux = flux * u.mJy
    dist = dist * u.pc
    freq = freq * u.MHz
    if jupiter == False:
        radius = radius * u.R_sun
    else:
        radius = radius * u.R_jup

    omega = (np.pi * np.arctan(radius / dist)**2).value

    tb = (flux * const.c**2 / (2 * const.k_B * freq**2 * omega)).decompose()

    # tb_err = tb * np.sqrt((flux_err / flux)**2 + 2 * (dist_err / dist)**2 + 
    #                       2 * (radius_err / radius)**2 + 2 * (freq_err / freq)**2)

    # return tb.to(u.K), tb_err.to(u.K)

    print(f"Brightness Temperature: {tb:.3e}")
    return tb


if __name__ == "__main__":
	main()
