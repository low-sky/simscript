#----------------------------------------------------------
# Python Routine to produce FITS files from RADMC output
#
# two fits files are produced: 1) PPV, and 2) the mom0 map
#     units of the PPV pixels are 'erg/s/cm^2/Hz/ster', 
#     and that of the mom0 map is in K km/s (e.g. Wco)
#
# Created June 7, 2013
#
# toK -- Binary flag to write out a PPV data cube in K
#----------------------------------------------------------

import matplotlib
matplotlib.use('Agg')
import numpy as np
import readRADMC as RADMC
import astropy.io.fits as pyfits
import math
import sys

pi = math.pi
def makefits (fitsfile='image.fits', dpc = None, mom0fitsfile='mom0.fits', lambda0=2600.7576, toK = True):
	Boltz = 1.380700e-16
	ckms = 2.9979e5
	a = RADMC.readimage()

	lambdanum=a.nrfr
	xnum=a.nx
	ynum=a.ny
	xsize=a.sizepix_x
	ysize=a.sizepix_y

	a.image[a.image==0.0]=np.nan

#	ind = 1
#	dpc = 1

	x1 = a.x[0]
	y1 = a.y[0]

	firstlambda = a.wavelength[0]

	if (lambdanum > 1):
		deltalam = float(a.wavelength[1]) - float(a.wavelength[0])
		v0 = ckms * (1 - lambda0/a.wavelength[0])
		v1 = ckms * (1 - lambda0/a.wavelength[1])
		delta = v1 - v0
	else:
		deltalam = a.wavelength[0]

	if(toK == True):
		scalefac = (0.5*(a.wavelength[ind]*1e-4)**2)/Boltz # Uses a single wavelength for the sake of memory
	else:
		scalefac = 1

	hdu = pyfits.PrimaryHDU(data = (scalefac*a.image))
	hd = hdu.header
	
	if(dpc != None):
		if (dpc == 0):
			print('0 pc to object?')
		hd['CDELT1'] = -xsize/(dpc*3.086e18)*(2*pi/360)
		hd['CRPIX1'] = xnum/2
		hd['CRVAL1'] = 180e0
		hd['CTYPE1'] = 'RA---CAR'
		hd['CDELT2'] = ysize/(dpc*3.086e18)*(2*pi/360)
		hd['CRVAL2'] = 0.0
		hd['CTYPE2'] = 'DEC--CAR'
		hd['EQUINOX'] = 2000.0
		hd['radesys'] = 'FK5'
		hd['TELESCOP'] = 'RADMC3D'
		hd['object'] = 'SIM@ %.0f PC'%dpc
	else:
		hd['CDELT1'] = xsize
		hd['CRPIX1'] = 1.0
		hd['CRVAL1'] = x1
		hd['CTYPE1'] = 'position - cm'
		hd['CDELT2'] = ysize
		hd['CRPIX2'] = 1.0
		hd['CRVAL2'] = y1
		hd['CTYPE2'] = 'position - cm'
	if (toK == True):
		lam0 = 2.99792458e8/110.2013542798e9*1e6
		vaxis = 2.99792458e8*(lam0-a.wavelength)/lam0
		dv = np.median(vaxis-Shift(vaxis,1)) 
		ind = np.argmin(abs(vaxis))
		hd['CDELT3'] = dv
		hd['CRPIX3'] = ind + 1
		hd['CRVAL3'] = vaxis[ind+1]
		hd['CTYPE3'] = 'V0PT'
		hd.update('BUNIT','K','Tmb')
	else:
		scalefac = 1
		hd['CDELT3'] = deltalam
		hd['CRPIX3'] = 1.0
		hd['CRVAL3'] = firstlambda
		hd['CTYPE3'] = 'Wavelength - micron'
		hd['BUNIT'] = 'erg/s/cm^2/Hz/ster'
#		print(type(a.image))
	hdu.writeto(fitsfile,clobber=True)
#	mom = pyfits.PrimaryHDU(data = mom0fitsfile)
#	momhd = mom.header
#	if len(mom0fitsfile) > 0 and lambdanum > 1:
#		mom0arr = np.zeros((xnum,ynum), 'f')
#		for i in range(xnum - 1):
#			for j in range(ynum - 1):
#				for k in range(lambdanum - 1):
#					#print(mom0arr[1,1,1])
#					#print(a.image[1,1])
#					mom0arr[i,j] = mom0arr[i,j] + 0.5 * (a.wavelength(k)*1.0e-4)**(2*a.image[i,j,k]*deltav/Boltz)

#		hd['BUNIT'] = 'K km/s'
#		momhd = hd
#		mom.writeto(mom0arr)
def Shift(l,n):
	l = list(l)
	n = int(n)
	n = n % len(l)
	n = -n
	shifted = l[n:] + l[:n]
	return shifted
#print(sys.argv)
#makefits(sys.argv[1][1:-1],sys.argv[2])
#makefits('image.out')
