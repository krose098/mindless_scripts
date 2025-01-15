from pathlib import Path
from astropy.time import Time
from astropy import coordinates as coord, units as u
from scipy.signal import find_peaks, peak_widths
from mindless_scripts.astro_tools.unicoord import unicoord
import pandas as pd
import click
import matplotlib.pyplot as plt


def process_astro_data(csv, observatory_name, format='isot', coords=None,galactic=False):
    """
    Process astronomical time series data from a pandas DataFrame.
    
    Parameters:
    -----------
    df : pandas.DataFrameå
        Input DataFrame containing time, flux, and flux error columns
    observatory_name : str
        Name of the observatory (must be recognized by astropy)
        
    Returns:
    --------
    pandas.DataFrame
        Processed DataFrame with standardized time column
    """
    time_col = None
    flux_col = None
    flux_err_col = None
    df = pd.read_csv(str(csv))
    
    for col in df.columns:
        col_lower = col.lower()
        if 'time' in col_lower or 'date' in col_lower:
            time_col = col
        elif 'flux' in col_lower and 'err' not in col_lower:
            flux_col = col
        elif 'err' in col_lower and 'flux' in col_lower:
            flux_err_col = col
    
    if not all([time_col, flux_col, flux_err_col]):
        missing = []
        if not time_col: missing.append('time')
        if not flux_col: missing.append('flux')
        if not flux_err_col: missing.append('flux error')
        raise ValueError(f"Could not find columns for: {', '.join(missing)}")
    
    if observatory_name.lower() == "atca":
        observatory_location = ('149.5646d', '-30.3139d','237m')
    else:
        try:
            observatory_location = coord.EarthLocation.of_site(observatory_name)
        except Exception as e:
            raise ValueError(f"Invalid observatory name: {observatory_name}. Error: {str(e)}")
  
    sample_time = df[time_col].dropna().iloc[0]
    result_df = pd.DataFrame()
    
    try:
        if pd.api.types.is_datetime64_any_dtype(df[time_col]):
            print(f"Detected pandas datetime format in time column")
            # Convert pandas datetime series to numpy array of datetime64
            datetime_array = df[time_col].values
            times = Time(datetime_array, format='datetime64', scale='utc',
                        location=observatory_location)
        else:
            if isinstance(sample_time, str):
                if 'T' in sample_time:
                    print(f"Detected ISOT format in time column: {sample_time}")
                    times = Time(df[time_col].values, format='isot', scale='utc',
                               location=observatory_location)
                else:
                    print(f"Detected datetime string format in time column: {sample_time}")
                    datetime_objects = pd.to_datetime(df[time_col]).values
                    times = Time(datetime_objects, format='datetime64', scale='utc',
                               location=observatory_location)
            elif isinstance(sample_time, (int, float)):
                if 5e4 < sample_time < 1e5:
                    print(f"Detected Modified Julian Date (MJD) format: {sample_time}")
                    times = Time(df[time_col].values, format='mjd', scale='utc',
                               location=observatory_location)
                elif 2.4e6 < sample_time < 2.6e6:
                    print(f"Detected Julian Date (JD) format: {sample_time}")
                    times = Time(df[time_col].values, format='jd', scale='utc',
                               location=observatory_location)
                elif 1e9 < sample_time < 2e9:
                    print(f"Detected Unix timestamp format: {sample_time}")
                    times = Time(df[time_col].values, format='unix', scale='utc',
                               location=observatory_location)
                else:
                    raise ValueError(f"Could not identify numerical time format for value: {sample_time}")
            else:
                raise ValueError(f"Could not identify time format for value: {sample_time}")

        result_df['flux'] = df[flux_col]
        result_df['flux_err'] = df[flux_err_col]

        if format == 'mjd':
            result_df['time'] = times.mjd
        elif format == 'isot':
            result_df['time'] = times.isot
        elif format == 'tdb':
            pos_eq, pos_gal = unicoord(coords,galactic,display=False)
            result_df['time'] = times.tdb.mjd+times.light_travel_time(pos_eq).value

    except Exception as e:
        raise ValueError(f"Error converting times: {str(e)}")
    
    return result_df

def peak_finder(df,thresh):
    """
    Find peaks in a time series DataFrame.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input DataFrame containing time, flux, and flux error columns
    thresh : float
        Threshold for peak detection
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing peak times and fluxes
    """

    x = df.flux
    x_err = df.flux_err
    t_raw = df.time
    t = (t_raw - t_raw[0])*24
    tmin=t.min()
    tmax=t.max()
    peak_thresh = thresh*x.std()
    n=tmax/len(t)

    peaks, _ = find_peaks(x, height=1.0*peak_thresh)
    fwhm_width, fwhm_height, fwhm_left, fwhm_right  = peak_widths(x, peaks, rel_height=0.5)

    t_left = []
    t_right = []
    t_width = []
    y_height = []
    for idx, peak in enumerate(peaks):
        t_left.append(t[peak]/peak*fwhm_left[idx])
        t_right.append(t[peak]/peak*fwhm_right[idx])
        t_width.append(t[peak]/peak*fwhm_width[idx])
        y_height.append(x[peak]/peak*fwhm_height[idx])


    fig = plt.figure(figsize=(10,7))
    ax1 = fig.add_subplot(1, 1, 1) 
    ax1.plot(t,x)#, label="Stokes I (ATCA 5)")
    ax1.plot(t[peaks], x[peaks]*1.0, marker="x", markersize=15, linestyle='',label = 'Peak',alpha=0.4)
    ax1.errorbar(t[peaks], x[peaks], yerr=x_err[peaks], fmt='none', ecolor='gray', capsize=4,alpha=0.4)
    ax1.hlines(y = 0, xmin = tmin, xmax = tmax, linestyles=":", color='gray')
    for i in range(thresh):
        if i % 2 != 0 and i>2:
            ax1.hlines(y = i*peak_thresh/thresh, xmin = tmin, xmax = tmax, linestyles="--", 
                       color='C{}'.format(i),label = '{}σ'.format(i),alpha=0.3)
    ax1.set_xlabel(r"$\Delta$t [hours]",fontsize = 18)
    ax1.set_ylabel(r'Flux Density [mJy beam $^{-1}$]',rotation = 90,fontsize=18)
    ax1.tick_params(axis='both', which='major', labelsize=20)
    ax1.set_xlim(tmin,tmax)
    ax1.legend(fontsize=14,loc="upper left",ncol=2)


    fig.tight_layout(pad=1.5)
    fig.savefig('temp_fig.pdf',dpi=300)

    print(46*"*")
    print("The {}σ peaks are found at:".format(thresh))
    print(46*"=")
    print("Date [MJD]| Δt [hr] | Flux Density [mJy/beam]|")
    print(46*"-")
    for idx, peak in enumerate(peaks):
        print("{:.3f} | {:.3f}   | {:.3f} +/- {:.3f}\t     |".format(t_raw[idx],t[peak], x[peak], x_err[peak]))
    print(46*"=")
    print("with a median peak flux of {:.3f} +/- {:.3f} mJy/beam".format(x[peaks].median(),x_err[peaks].median()))
    print(46*"*")
    return


@click.command()
@click.argument("observatory_name", type=str)
@click.argument("csv", type=Path)
@click.argument("format", type=str,default='mjd')
@click.argument("thresh", type=int,default=4)
@click.option("-c","coords", nargs=1, type=str,default=None, help='Coordinates for light travel time correction')
def main(observatory_name,csv,format,coords,thresh):
    df = process_astro_data(csv,observatory_name,format,coords)
    peak_finder(df,thresh)
if __name__ == "__main__":
	main()

