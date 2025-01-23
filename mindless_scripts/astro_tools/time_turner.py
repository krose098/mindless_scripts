from astropy.time import Time
import click
 
def time_to_mjd(time_raw,num):
    '''
    Convert a time string or numerical value from various formats into Modified Julian Date (MJD) format.
    '''
    if num is True:
        time_raw = float(time_raw)
    
    if isinstance(time_raw, str):
        if 'T' in time_raw:
            print(f"Detected ISOT format: {time_raw}")
            time = Time(time_raw, format='isot', scale='utc')
        else:
            print(f"Detected ISO format: {time_raw}")
            time = Time(time_raw, format='iso', scale='utc')
    elif isinstance(time_raw, float):
        if 5e4 < time_raw < 1e5:
            print(f"Detected Modified Julian Date (MJD) format: {time_raw}")
            time = Time(time_raw, format='mjd', scale='utc')
        elif 2.4e6 < time_raw < 2.6e6:
            print(f"Detected Julian Date (JD) format: {time_raw}")
            time = Time(time_raw, format='jd', scale='utc')
        elif 1e9 < time_raw < 2e9:
            print(f"Detected Unix timestamp format: {time_raw}")
            time = Time(time_raw, format='unix', scale='utc')
        else:
            raise ValueError(f"Could not identify numerical time format for value: {time_raw}")
    else:
        raise ValueError(f"Could not identify time format for value: {time_raw}")
    print("MJD format: {:.3f}".format(time.mjd))
    return time.mjd

@click.command()
@click.argument("time", nargs=1, type=str)
@click.option("-n","--num", is_flag=True, show_default=True,default=False, help="If set, assumes time is a numerical value")
def main(time,num):
    time_to_mjd(time,num)
    return
if __name__ == "__main__":
	main()