#!/usr/bin/env python

from yt.mods import *
import glob
import astropy.io.fits as pyfits
import subprocess
import commands
import sys
import os
import tempfile
import re
import numpy as np
import array 
import readRADMC as RADMC

def FileSetup(targetdir,timestep=21,face=0,level=0,ppdir='/home/eros/simscript/postproc/'):

    homedir = os.getenv('HOME')
    if not os.access(targetdir,os.R_OK):
        workdir = os.getenv('PWD')
    else:
        workdir=''

    execfile(ppdir+'problem_setup.py')
    execfile(ppdir+'makefits.py')
    execfile(ppdir+'readRADMC.py')

    name = targetdir.replace('/','')

    tempdir = os.getenv('PWD')+'/tmp_'+name
    print(commands.getoutput('rm -rf '+tempdir))
    os.mkdir(tempdir)

    print('Building files with prefix %s' %(name))
    num = timestep
    path = workdir+'/'+targetdir

    fn = path+"/DD%04i/DD%04i" % (num, num) # parameter file to load

    initfile = open(fn,'r')
    inittxt = initfile.readlines()
    initfile.close()
    EOSSound = float((([line for line in inittxt if "EOSSoundSpeed" in line][0]).split())[2])
    Mu = float((([line for line in inittxt if "Mu " in line][0]).split())[2])
    Tk = EOSSound**2*Mu*1.67e-24/1.38e-16

    print(inittxt)
    pf = load(fn) # load data
    # This is the resolution we will extract at

    DIMS = (pf.domain_dimensions)*(2**level)
    print('Extracting a to hierarchy level %i given total grid size %i' % (level,DIMS[0]))

    # Now, we construct an object that describes the data region and structure we
    # want

    cube = pf.h.covering_grid(int(level), # The level we are willing to extract to; higher
                              # levels than this will not contribute to the data!
                              # Now we set our spatial extent...
                              left_edge=[0.0, 0.0, 0.0],
                              right_edge=[1.0, 1.0, 1.0],
                              # How many dimensions along each axis
                              dims=DIMS,
                              # And any fields to preload (this is optional!)
                              fields=["Density"])
    
    FlatFileName = tempdir+'/%s_flatrho_%04i_%02i.fits' %(name, num,face)
    
    pyfits.writeto(FlatFileName, cube["Density"],clobber=True)
    pyfits.append(FlatFileName, cube["x-velocity"])
    pyfits.append(FlatFileName, cube["y-velocity"])
    pyfits.append(FlatFileName, cube["z-velocity"])
    
    command = 'cp '+ppdir+'/dustkappa_silicate.inp '+tempdir
    print(command)
    commands.getoutput(command)
    command = 'cp '+ppdir+'/molecule_13co.inp '+tempdir
    commands.getoutput(command)
    command = 'cp '+ppdir+'/camera_wavelength_micron.inp '+tempdir
    commands.getoutput(command)

    OutputDict = {'GasTemp':Tk,
                  'FileName':FlatFileName,
                  'TempDir':tempdir,
                  'GridSize':DIMS[0]}

    return(OutputDict)
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

pi = np.pi
def MakeFits (fitsfile='image.fits', dpc = None,
              mom0fitsfile='mom0.fits', lambda0=2600.7576, toK = True):

    Boltz = 1.380700e-16
    ckms = 2.9979e5
    a = RADMC.readimage(filename='image.out')
    
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
            scalefac = (0.5*(firstlambda*1e-4)**2)/Boltz # Uses a single wavelength for the sake of memory
    else:
        scalefac = 1
            
    hdu = pyfits.PrimaryHDU(data = (scalefac*a.image))
    hd = hdu.header
            
    if(dpc != None):
        if (dpc != 0):
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
            print('0 pc to object?')
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
            hd['CTYPE3'] = 'VOPT'
            hd.update('BUNIT','K','Tmb')
	else:
            scalefac = 1
            hd['CDELT3'] = deltalam
            hd['CRPIX3'] = 1.0
            hd['CRVAL3'] = firstlambda
            hd['CTYPE3'] = 'Wavelength - micron'
            hd['BUNIT'] = 'erg/s/cm^2/Hz/ster'
            #		print(type(a.image))
    print("Writing Output to %s" % fitsfile)
    hdu.writeto(fitsfile,clobber=True)

def Shift(l,n):
	l = list(l)
	n = int(n)
	n = n % len(l)
	n = -n
	shifted = l[n:] + l[:n]
	return shifted

def ProblemSetup(file, face = 0, dust_temp = 10, gas_temp = 10, SetGasTemp = False):

	densityscale = 1.0
	tempscale = 1.0

	# Microturbulence? in cm s^-1

	inclmicrotub = 1
	turbvels = 0.05e5
	nphot = 10000

	# Specify constant dust temperature, so only line transfer can be done?
	dust_temp = 1
	const_dusttemp=dust_temp

	# Dust to gas ratio

	dust_gas_ratio=1.0/100.0

	# Mode of computing population levels:

	lines_mode = 3
	maxJbariter = 1000
	lvgconvergeto = 0.01

	#NOTE TO SELF: If there is an error in the output it is probably from this part of the code
	rho = pyfits.getdata(file, 0)*(densityscale)
	rho = rho.transpose(Shift([0,1,2],face))
	vx = (pyfits.getdata(file, 1)*(tempscale**(0.5))).transpose(Shift([0,1,2],face))
	vy = (pyfits.getdata(file, 2)*(tempscale**(0.5))).transpose(Shift([0,1,2],face))
	vz = (pyfits.getdata(file, 3)*(tempscale**(0.5))).transpose(Shift([0,1,2],face))

	n1 = long ((len(rho[1])))
	n2 = long ((len(rho[2])))
	n3 = long ((len(rho[3])))

	PC = 3.0857200e+18
	#NOTE TO SELF: PC parsec. From natconst? In cm?	
	sizex = 10*PC*(densityscale**(-0.5))*(tempscale**(0.5))
	sizey = 10*PC*(densityscale**(-0.5))*(tempscale**(0.5))
	sizez = 10*PC*(densityscale**(-0.5))*(tempscale**(0.5))
	#NOTE TO SELF: DindgenEdit only exists in idl

	x = DindgenEdit(n1+1,sizex/n1, -sizex/2.)
	y = DindgenEdit(n2+1,sizey/n2, -sizey/2.)
	z = DindgenEdit(n3+1,sizez/n3, -sizez/2.)

	temp = TDarr([n1, n2, n3],gas_temp)

	# To number Density
	# Use Orion conversion: 2.33 * AMU

	nden = rho / float(2.33 * 1.67e-24)

	# Extract 13CO Density
	coden = nden * float(5.625e-5 / 77.0) # 1e-4

	# Make grid
	# amr_grid has positions of cell walls, so need to include nx+1... etc positions
	n1lims=n1-1
	n2lims=n2-1
	n3lims=n3-1
	n1lim=n1
	n2lim=n2
	n3lim=n3
	x1lim = x
	x2lim = y
	x3lim = z
	#NOTE TO SELF: Changed 1: to 0:
	v1lim = vx[1:n1lim, 1:n2lim, 1:n3lim]
	v2lim = vy[1:n1lim, 1:n2lim, 1:n3lim]
	v3lim = vz[1:n1lim, 1:n2lim, 1:n3lim]
	templim=temp[0:n1lim][0:n2lim][0:n3lim]
	H2denlim=nden[1:n1lim, 1:n2lim, 1:n3lim]
	codenlim=coden[1:n1lim, 1:n2lim, 1:n3lim]
	dustdenlim = rho[1:n1lim, 1:n2lim, 1:n3lim] * dust_gas_ratio
	#
	# Write the wavelength_micron.inp file
	#

	lambda1 = float(0.1)
	lambda2 = float(7.0)
	lambda3 = float(25.)
	lambda4 = float(10000.000)
	n12     = 20
	n23     = 100
	n34     = 30
	lam12   = lambda1 * np.power((lambda2/lambda1), DindgenEdit(n12,1/(float(1)*n12)))
	lam23   = lambda2 * np.power((lambda3/lambda2), DindgenEdit(n23,1/(float(1)*n23)))
	lam34   = lambda3 * np.power((lambda4/lambda3),
                                     DindgenEdit(n34,1/(float(1)*(float(n34-1)))))
	lam     = list(lam12)
	lam.extend(lam23)
	lam.extend(lam34)
	nlam    = len(lam)

	#
	# Write the wavelength file
	#

	temp_file = open('wavelength_micron.inp', 'w')
	temp_file.write("%s\n"%nlam)
	for ilam in range(0,nlam):
		temp_file.write("%s\n"%lam[ilam])
	temp_file.close()

	#
	# Write the grid file
	#

	temp_file = open('amr_grid.inp', 'w')
	temp_file.write("%s\n"%1) # iformat
	temp_file.write("%s\n"%0) # AMR grid style (0 = regular grid, no AMR)
	temp_file.write("%s\n"%0) # Coordinate system
	temp_file.write("%s\n"%0) # Gridinfo
	temp_file.write("%s\t%s\t%s\n"%(1,1,1)) # Include x, y, z coordinate
	for item in [n1lim-1,n2lim-1,n3lim-1]:
		temp_file.write("%s\t"%item) # Size of grid
	temp_file.write("\n")
	for i in range(0,n1lim):
		temp_file.write("%s\n"%x1lim[i])
	for i in range(0,n2lim):
		temp_file.write("%s\n"%x2lim[i])
	for i in range(0,n3lim):
		temp_file.write("%s\n"%x3lim[i])
	temp_file.close()

	#
	# Write the density file
	#

	temp_file = open('dust_density.inp', 'w')
	temp_file.write("%s\n"%1) # Format number
	temp_file.write("%s\n"%(n1lims*n2lims*n3lims)) # Nr of cells
	temp_file.write("%s\n"%1)
	for iy in range(0,n2lim-1):
		for ix in range(0,n1lim-1):
			for iz in range(0,n3lim-1):
				temp_file.write("%s\n"%dustdenlim[iz,ix,iy])
	temp_file.close()

	#
	# Write the molecule number density file
	#

	comom0arr = array.array('f',[n1lim,n2lim])
	dz = x3lim[1]-x3lim[0]
	temp_file = open('numberdens_13co.inp', 'w')
	temp_file.write("%s\n"%1) # Format number
	temp_file.write("%s\n"%(n1lims*n2lims*n3lims)) # Nr of cells
	for iy in range(0,n2lim-1):
		for ix in range(0,n1lim-1):
			for iz in range(0,n3lim-1):
				temp_file.write("%s\n"%codenlim[iz,ix,iy])
	temp_file.close()

	#
	# Write the collisional partner number density field
	#

	#NOTE TO SELF: H2mom0arr not used anywhere but these two lines
	#H2mom0arr = fltarr(n1lim,n2lim)
	#H2mom0arr(*,*) = 0.
	temp_file = open('numberdens_h2.inp', 'w')
	temp_file.write("%s\n"%1) # Format number
	temp_file.write("%s\n"%(n1lims*n2lims*n3lims)) # Nr of cells
	for iy in range(0,n2lim-1):
		for ix in range(0,n1lim-1):
			for iz in range(0,n3lim-1):
				temp_file.write("%s\n"%H2denlim[iz,ix,iy]) # for units of cc^-3
	temp_file.close()

	#
	# Write the gas velocity field
	#

	temp_file = open('gas_velocity.inp', 'w')
	temp_file.write("%s\n"%1)
	temp_file.write("%s\n"%(n1lims*n2lims*n3lims))
	for iy in range(0,n2lim-1):
		for ix in range(0,n1lim-1):
			for iz in range(0,n3lim-1):
				temp_file.write("%s\t%s\t%s\n" % (v1lim[iz,ix,iy],
                                                                  v2lim[iz,ix,iy],
                                                                  v3lim[iz,ix,iy]))
	temp_file.close()

	# Write the microturbulence file
	#

	if (inclmicrotub == 1):
		temp_file = open('microturbulence.inp', 'w')
		temp_file.write("%s\n"%1) # Format number
		temp_file.write("%s\n"%(n1lims*n2lims*n3lims)) # Nr of cells
		for iy in range(0,n2lim-1):
			for ix in range(0,n1lim-1):
				for iz in range(0,n3lim-1):
					temp_file.write("%s\n"%turbvels)
		temp_file.close()

	#
	# Dust opacity control file
	#

	temp_file = open('dustopac.inp', 'w')
	temp_file.write("2               Format number of this file \n")
	temp_file.write("1               Nr of dust species\n")
	temp_file.write("============================================================================\n")
	temp_file.write("1               Way in which this dust species is read\n")
	temp_file.write("0               0=Thermal grain\n")
	temp_file.write("silicate        Extension of name of dustkappa_***.inp file\n")
	temp_file.write("----------------------------------------------------------------------------\n")
	temp_file.close()

	#
	# Write the lines.inp control file
	#

	temp_file = open('lines.inp', 'w')
	temp_file.write('2 \n')
	temp_file.write('1 \n')
	temp_file.write('13co    leiden    0    0  1 \n')
	temp_file.write('h2 \n')
	temp_file.close()

	#
	# Write the radmc.inp control file
	#

	temp_file = open('radmc.inp', 'w')
	temp_file.write('nphot = %s \n' % nphot)
	temp_file.write('scattering_mode_max = 0\n')
        if not(SetGasTemp):
            temp_file.write('tgas_eq_tdust   = 1\n')

	temp_file.write('lines_mode = %s\n' % lines_mode)
	temp_file.write('lines_nonlte_maxiter = %s\n' % maxJbariter)
	temp_file.close()
	if(SetGasTemp):
            temp_file = open('gas_temperature.inp', 'w')
            temp_file.write("%s\n"%1) # Format number
            temp_file.write("%s\n"%(n1lims*n2lims*n3lims)) # Nr of cells)
            for iy in range(0,n2lim-1):
                for ix in range(0,n1lim-1):
                    for iz in range(0,n3lim-1):
                        temp_file.write("%s\n" % templim[iz][ix][iy])
                        temp_file.close()

	if(dust_temp==1):
		temp_file = open('dust_temperature.dat', 'w')
		temp_file.write("%s\n"%1) # Format number
		temp_file.write("%s\n"%(n1lims*n2lims*n3lims)) # Nr of cells
		temp_file.write("%s\n"%1) # Nr of species
		for iy in range(0,n2lim-1):
			for ix in range(0,n1lim-1):
				for iz in range(0,n3lim-1):
					temp_file.write("%s\n"%const_dusttemp)
		temp_file.close()

def DindgenEdit(n, multi = 1, offset = 0):
	return_val = array.array('f')
	for i in range(0,n):
		return_val.append(i*multi + offset)
	return return_val

def TDarr(arg, offset = 0.):
	#Three-Dimensional Array
	return_val = array.array('f')
	return_val = [offset]*arg[0]
	for i in range(arg[0]):
		return_val[i] = [offset]*arg[1]
		for j in range(arg[1]):
			return_val[i][j] = [offset]*arg[2]
	return return_val
