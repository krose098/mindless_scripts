

#TODO: Add a function to generalise for all coordinate inputs
#TODO: Add saveplot option (set to default True)
#TODO: Add relevant structure to make it a click script
#TODO: Add twin axis for AEDT time

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord
from pytz import timezone
from astropy.time import Time, TimezoneInfo

from astroplan import Observer
from astroplan import FixedTarget
from astroplan.plots import plot_altitude

import warnings

def astro_planner(name,ra,dec,date,time,horizon):
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
    target_coordinates = SkyCoord(ra,dec, frame='icrs')
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
    standard_style = {'linestyle': '-', 'marker': None, 'alpha': 0.5, 'linewidth': 2}
    plot_altitude(target, observer, observe_time, airmass_yaxis=False,min_altitude=horizon,style_kwargs=standard_style)
    # plot_altitude(phase_cal, observer, observe_time, airmass_yaxis=False,min_altitude=horizon)
    # plot_altitude(phase_cal2, observer, observe_time, airmass_yaxis=False,min_altitude=horizon)
    plot_altitude(calibrator, observer, observe_time, airmass_yaxis=False,min_altitude=horizon,style_kwargs=standard_style)
    plot_altitude(calibrator2, observer, observe_time, airmass_yaxis=False,min_altitude=horizon,style_kwargs=standard_style)

    plt.legend()
    plt.tight_layout()
    plt.show()
    return


#example usage: astro_planner('J1745','17h45m08.90s', '-50d51m49.9s','2024-12-22','23:00:00',horizon=15)
