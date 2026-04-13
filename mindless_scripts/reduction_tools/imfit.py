from casatasks import imfit, imstat
from pathlib import Path
import numpy as np

import click
 
@click.command()
@click.argument("ms", type=Path)



path = input("Enter the .fits path: ")
region_temp =input("Enter the source coordinates (AAhBBmCCs,XX.YY.ZZ format): ") #TODO generalise this
radius = input("Enter a radius for the peak flux region: ")
region = 'circle[[{}], {}arcsec]'.format(region_temp,radius)

source_name = "J"+region_temp[0:2]+region_temp[3:5]
stokes = input("Choose the stokes parameter:")
obs_name = input("Enter the observation identifier:")

name = source_name+"_"+obs_name+"_"+stokes

region_temp =input("Enter the empty coordinates (AAhBBmCCs,XX.YY.ZZ format): ") #TODO generalise this
radius = input("Enter a radius for the RMS test region: ")
clean_region = 'circle[[{}], {}arcsec]'.format(region_temp,radius)

imagename=path #TODO save the input parameters and present them as an optional starting point for the next run

summary='fit_summary_{}'.format(name)
logfile='fit_log_{}'.format(name)

fit = imfit(imagename=imagename,region=region, chans='',stokes="",summary=summary,
            logfile=logfile)#,residual=name+'_residual')
stats=imstat(imagename=imagename,stokes=stokes, region=clean_region) #TODO make option for peak from imstat if imfit fails

print('\n'+42*'-')
print('For {}:'.format(name))

rms=stats["rms"]*10**3 #rms in mJy

print("RMS = {} mJy".format(round(rms[0],4)))
print(42*'-')

temp=fit['results']['component0']['peak']#['value']
peak_flux = temp['value']*10**3 #peak flux in mJy
peak_err = temp['error']*10**3 #peak flux error in mJy

print('Peak Flux = {} +/- {} mJy'.format(round(peak_flux,4),round(peak_err,4)))
scale=0.10
combined_err = (peak_err**2+rms**2+(scale*peak_flux)**2)**(0.5)
print("Combined Error = {} mJy".format(round(combined_err[0],4)))

coords = fit['results']['component0']['shape']['direction']
ra = coords['m0']['value']*180/np.pi + 360 #conversion from radians to degrees
ra_err = coords['error']['longitude']['value']/3600 #conversion from arcsec to degrees
dec = coords['m1']['value']*180/np.pi #conversion from radians to degrees
dec_err = coords['error']['latitude']['value']/3600 #conversion from arcsec to degrees

print("Centroid:\nRA = {:.6f} +/- {:.6f} deg \nDec = {:.6f} +/- {:.6f} deg ".format(ra,ra_err,dec,dec_err))
print(42*'-'+'\n')

def main(ms):
	flagdata(
		vis=ms, 
		mode="tfcrop",
		)
	pass

if __name__ == "__main__":
	main()