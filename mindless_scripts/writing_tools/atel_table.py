
from pathlib import Path
import pandas as pd
from datetime import datetime, timedelta
from mindless_scripts.astro_tools.unicoord import unicoord
import click
 

def mjd_to_date(mjd_list):
    formatted_dates = []
    for mjd in mjd_list:
        # Convert MJD to Julian Date (JD)
        jd = mjd + 2400000.5
        
        # Convert JD to a datetime object
        jd_offset = jd - 1721425.5  # Julian date of 0001-01-01
        date = datetime(1, 1, 1) + timedelta(days=jd_offset)
        
        # Extract year, month, day, and fractional day
        year = date.year
        month = date.strftime("%b")  # Get the first three letters of the month
        day_fractional = date.day + (date.hour / 24.0) + (date.minute / 1440.0) + (date.second / 86400.0)
        
        # Format the date as required
        formatted_date = f"{year} {month} {day_fractional:.2f}"
        formatted_dates.append(formatted_date)
    
    return formatted_dates


@click.command()
# @click.option('-lm"--long", is_flag=True, help="Extended output")
@click.argument("csv", type=Path)
def main(csv):
    
    df = pd.read_csv(str(csv))
    flux = df['flux_peak']
    err = df['flux_err_quad']
    dates = df['obs_start']
    freqs = df['freq']
    dur = df['obs_length']/60.0
    ra = df['ra_deg_cont']
    dec = df['dec_deg_cont']
    
    formatted_dates = mjd_to_date(dates)
    for idx,date in enumerate(formatted_dates):
        print("{}: {:.2f} +/- {:.2f} mJy/beam {:.0f} MHz {:.0f} min".format(date, flux[idx], err[idx], round(freqs[idx]), round(dur[idx])))
    for idx,date in enumerate(formatted_dates):
         coords = "{} {}".format(ra[idx], dec[idx])
         pos_eq, pos_gal = unicoord(coords, galactic=False,display=False)
         print("(R.A., Dec) = ({} deg) = = ({})".format(pos_eq.to_string(style="decimal", precision=4), pos_eq.to_string(style="hmsdms", precision=1)))
    
    # print("Mean position (R.A., Dec) = ({})".format(pos_eq.to_string(style="decimal", precision=4)))
    # if l == True:
    #     # Print the full DataFrame
    #     print(df.to_string(index=False))
    # return

if __name__ == "__main__":
	main()
