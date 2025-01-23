
from pathlib import Path
from astropy.time import Time
from astropy import coordinates as coord, units as u
from mindless_scripts.astro_tools.unicoord import unicoord
from mindless_scripts.astro_tools.time_turner import time_to_mjd
from scipy.signal import find_peaks
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def combApprox(tmin,tmax,zero_peak_time,T,N=8,Fs=1000):
    t = np.arange(tmin, tmax, 1/Fs)
    sigSum = np.ones_like(t)
    for n in range(1,N+1):
        part = np.cos(2*np.pi*n*(t-zero_peak_time)/T)
        sigSum = sigSum + part
    return t, sigSum/sigSum.max()

def toa_pred(t_obs,peak_times_obs,t_pred_start,duration,observatory_name,coords,galactic,T,err_T):
    if observatory_name.lower() == "atca":
        observatory_location = ('149.5646d', '-30.3139d','237m')
    else:
        try:
            observatory_location = coord.EarthLocation.of_site(observatory_name)
        except Exception as e:
            raise ValueError(f"Invalid observatory name: {observatory_name}. Error: {str(e)}")
        
    pos_eq, pos_gal = unicoord(coords,galactic,display=False)

    T=T/24 
    tmin1 = t_obs.min()
    tmax1 = t_obs.max()
    t_pred_start_mjd = time_to_mjd(t_pred_start,num=False)
    t_pred_start_mjd = Time(t_pred_start_mjd, format='mjd', scale='utc', location=observatory_location)
    t_pred_start_bary = t_pred_start_mjd.tdb.mjd + t_pred_start_mjd.light_travel_time(pos_eq) 
    tmin2 = t_pred_start_bary.value
    tmax2 = tmin2 + duration/24

    time, signal = combApprox(tmin1,tmax2,peak_times_obs[0],T)

    tmax2_rel = (tmax2-tmin1)*24
    tmax1_rel = (tmax1-tmin1)*24
    tmin2_rel = (tmin2-tmin1)*24
    peak_obs_rel = (peak_times_obs[0]-tmin1)*24
    time_rel, signal_rel = combApprox(0,tmax2_rel,peak_obs_rel,T*24)

    peaks_pred, _ = find_peaks(signal, height=0.80)
    peak_times_full = time[peaks_pred]
    peak_times_pred = peak_times_full[(peak_times_full >= tmin2) & (peak_times_full <= tmax2)]

    time_pred = Time(peak_times_pred, format='mjd', scale='tdb', location=observatory_location)
    time_pred_ltt = time_pred - time_pred.light_travel_time(pos_eq).value
    time_pred_mjd = time_pred_ltt.utc.mjd 
    time_pred_iso = time_pred_ltt.utc.iso

    # if corrected is True:
    #     time_pred = Time(peak_times_pred, format='mjd', scale='tdb', location=observatory_location)
    #     time_pred_ltt = time_pred - time_pred.light_travel_time(pos_eq)
    #     time_pred_mjd = time_pred_ltt.tdb.mjd 
    #     time_pred_iso = time_pred_ltt.tdb.iso

    # pos_eq, pos_gal = unicoord(coords,galactic,display=False)
    # time_pred_ltt = time_pred_bary+time_pred.light_travel_time(pos_eq).value

    time_peak_hrs = (time_pred_mjd - tmin2)*24
    print(60*"*")
    print("The predicted peaks will occur at:")
    print(48*"=")
    print("Date [MJD] |\t   Date [ISO]\t     | Δt [hr] |")
    print(48*"-")
    for idx, peak in enumerate(time_pred_mjd):
        print("{:.4f} | {} | {:.3f}   |".format(peak, time_pred_iso[idx], time_peak_hrs[idx]))
    print(48*"=")
    print("with a total of {} peaks during the {}hr observation.".format(len(time_pred_mjd),duration))
    print(60*"*")

    fig = plt.figure(figsize=(10, 8))
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)
    ax3 = ax1.twiny()
    ax4 = ax2.twiny()
    ax3.plot(time_rel, signal_rel, 'g-', lw=1, zorder=-1, alpha=.0)
    ax4.plot(time_rel, signal_rel, 'g-', lw=1, zorder=-1, alpha=.0)
    ax1.plot(time, signal, 'r-', lw=1, zorder=-1, alpha=.3)
    ax2.plot(time, signal, 'r-', lw=1, zorder=-1, alpha=.3)
    for idx,peak_time in enumerate(peak_times_obs):
        ax1.axvline(x=peak_time,lw=2,alpha=0.7)
        ax1.axvspan(peak_time-T*err_T,peak_time+T*err_T, facecolor='cornflowerblue', alpha=.15)
    for idx,peak_time in enumerate(peak_times_pred):
        ax2.axvline(x=peak_time,lw=2, linestyle='--',alpha=0.7)
        ax2.axvspan(peak_time-T*err_T,peak_time+T*err_T, facecolor='cornflowerblue', alpha=.15)

    ax1.set_xlim(tmin1,tmax1)
    ax2.set_xlim(tmin2,tmax2)
    ax3.set_xlim(0,tmax1_rel)
    ax4.set_xlim(0,tmax2_rel-tmin2_rel)
    ax1.set_xlabel('Barycentric Date [MJD]',fontsize=14)
    ax1.set_ylabel('Flux [norm.]',fontsize=14)
    ax2.set_xlabel('Barycentric Date [MJD]',fontsize=14)
    ax2.set_ylabel('Flux [norm.]',fontsize=14)
    ax3.set_xlabel(r'$\Delta$t [hr]',fontsize=14)
    ax4.set_xlabel(r'$\Delta$t [hr]',fontsize=14)
    ax1.set_title('Observed Peaks',fontsize=17,fontdict={'style': 'italic'})
    ax2.set_title('Predicted Peaks',fontsize=17,fontdict={'style': 'italic'})

    ax3.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax3.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax4.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax4.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    for ax in [ax1,ax2,ax3,ax4]:
        ax.tick_params(axis='x',which='major',direction='in',length=7,labelsize=14)
        ax.tick_params(axis='x',which='minor',direction='in',length=5,labelsize=14)
        ax.tick_params(axis='y',which='major',direction='out',length=7,labelsize=12)
        ax.tick_params(axis='y',which='minor',direction='out',length=5,labelsize=12)

    fig.tight_layout(pad=2.5)
    fig.savefig('temp_fig_toa.pdf',dpi=300)    
    return 

def main(t_obs,peak_times_obs,t_pred_start,duration,observatory_name,coords,galactic,T,err_T):
    toa_pred(t_obs,peak_times_obs,t_pred_start,duration,observatory_name,coords,galactic,T,err_T)
    return
if __name__ == "__main__":
	main()