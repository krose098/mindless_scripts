
#TODO: Add twin axis for AEDT time (WIP)

import click
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import warnings
from mindless_scripts.astro_tools.unicoord import unicoord
from astropy.coordinates import EarthLocation, SkyCoord
from pytz import timezone
from astropy.time import Time, TimezoneInfo

from astroplan import Observer
from astroplan import FixedTarget
from astroplan.plots import plot_altitude
from astroplan.plots import dark_style_sheet

def astro_planner(name,coords,date,time,horizon,plot=True):
    ''' Takes in five parameters as strings: The source name, right ascension in the format of XXhYYmZZs, 
    declination in the format of XXdYYmZZs, observation date in the format YYYY-MM-DD, 
    and the 24-hour time in the format HH:MM:SS'''

    warnings.filterwarnings("ignore", category=UserWarning)

    ### Telescope Definitions ###
    longitude = '149d33m00.500s'
    latitude = '-30d18m46.385s'
    elevation = 237 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)
    observer = Observer(name='ATCA',
                location=location,
                timezone=timezone('Australia/Sydney'),
                description="Australia Telescope Compact Array, Narrabri")

    ### Source Definitions ###
    calibrator_coordinates = SkyCoord('19h39m25.026s', '-63d42m45.63s', frame='icrs')
    calibrator = FixedTarget(name='1934-638', coord=calibrator_coordinates)
    calibrator_coordinates2 = SkyCoord('08h25m26.869s', '-50d10m38.49s', frame='icrs')
    calibrator2 = FixedTarget(name='0823-500', coord=calibrator_coordinates2)
    target_coordinates = coords #SkyCoord(ra,dec, frame='icrs')
    target = FixedTarget(name=name, coord=target_coordinates)

    # phase_cal_coordinates = SkyCoord('16h36m55.375s', '-41d02m00.518s', frame='icrs')
    # phase_cal = FixedTarget(name='J1636-4101', coord=phase_cal_coordinates)
    # phase_cal2_coordinates = SkyCoord('17h17m38.597s', '-39d48m52.586s', frame='icrs')
    # phase_cal2 = FixedTarget(name='J1714-397', coord=phase_cal2_coordinates)


    ### Visibility Time Calculations ###
    # observe_time = Time('2023-02-01 12:00:00')
    datetime=date+" "+time
    observe_time = Time(datetime)

    source_rise=observer.target_rise_time(observe_time,target,horizon=horizon*u.deg)
    source_set=observer.target_set_time(observe_time,target,horizon=horizon*u.deg)
    utc_plus_eleven_hours = TimezoneInfo(utc_offset=11*u.hour)
    source_rise_aedt=str(source_rise.to_datetime(timezone=utc_plus_eleven_hours)).split('+')[0]
    source_set_aedt=str(source_set.to_datetime(timezone=utc_plus_eleven_hours)).split('+')[0]

    ### Visibility Printout ###
    print('The target is above the horizon for {:.1f} hours\n'.format(24*(source_set-source_rise).value))
    print('|  TZ\t|\tRise \t\t  |\t\tSet\t\t|') 
    print(65*"=")
    print('| UTC \t| {} | {}\t|'.format(source_rise.iso,source_set.iso)) 
    print(65*"-")
    print('| AEDT \t| {} | {} \t|'.format(source_rise_aedt[:-3],source_set_aedt[:-3])) 

    observe_time = observe_time + np.linspace(-7, 8, 70)*u.hour
    ### Plotting ###
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(2, 1, 1)
    standard_style = {'linestyle': '-', 'marker': None, 'alpha': 0.75, 'linewidth': 2}
    plot_altitude(target, observer, observe_time, ax=ax, airmass_yaxis=False,min_altitude=horizon,style_kwargs=standard_style,style_sheet=dark_style_sheet)
    # plot_altitude(phase_cal, observer, observe_time, airmass_yaxis=False,min_altitude=horizon)
    # plot_altitude(phase_cal2, observer, observe_time, airmass_yaxis=False,min_altitude=horizon)
    plot_altitude(calibrator, observer, observe_time, ax=ax, airmass_yaxis=False,min_altitude=horizon,style_kwargs=standard_style)
    plot_altitude(calibrator2, observer, observe_time, ax=ax, airmass_yaxis=False,min_altitude=horizon,style_kwargs=standard_style)
    
    # observe_time_aedt = observe_time + 11*u.hour+ np.linspace(-7, 8, 70)*u.hour
    # ax2 = ax.twiny()
    # ax2.plot(observe_time_aedt, np.ones_like(observe_time_aedt), ' ')
    ax.legend()
    fig.tight_layout()
    if plot:
        fig.savefig('obs_plan_{}_{}.pdf'.format(name,datetime),dpi=300)
        plt.show()
    return


@click.command()
@click.option("-p","--plot", is_flag=True, show_default=True,default=False, help="If set, save and display the plot")
@click.option("-g","--galactic", is_flag=True, show_default=True,default=False, help="If set, assumes galactic coordinates")
@click.argument("name", nargs=1, type=str)
@click.argument("coords", nargs=1, type=str)
@click.argument("date", nargs=1, type=str)
@click.argument("time", nargs=1, type=str)
@click.argument("horizon", nargs=1, type=float, default=15)
def main(name,coords,date,time,horizon,plot,galactic):
    pos_eq, __ = unicoord(coords,galactic,display=False)
    astro_planner(name,pos_eq,date,time,horizon,plot)
    return

if __name__ == "__main__":
	main()
     
#Original example usage: astro_planner('J1745','17h45m08.90s', '-50d51m49.9s','2024-12-22','23:00:00',horizon=15)

#New example usage: 'SN 2024abfo' '03h57m25.622s -46d11m07.67s' '2025-03-03' '05:00:00' 12.5