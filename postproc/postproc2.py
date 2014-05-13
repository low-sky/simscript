#!/usr/bin/env python

from yt.mods import *
import glob
import astropy.io.fits as pyfits
import subprocess
import sys
import os
import tempfile
import re

homedir = os.getenv('HOME')
workdir = os.getenv('PWD')

ppdir = homedir+'/postproc/'
execfile(ppdir+'problem_setup.py')
execfile(ppdir+'makefits.py')
execfile(ppdir+'readRADMC.py')

#tempdir = tempfile.mkdtemp()
tempdir = workdir

targetdir = sys.argv[1]
timestep = float(sys.argv[2])
face = float(sys.argv[3])
level = float(sys.argv[4])  

# we will use pyfits to write to our output file.
#import h5py

#name = sys.argv[1]
#num = float(sys.argv[2])
#path=sys.argv[4]

name = targetdir.replace('/','')
print(name)
num = timestep
path = workdir+targetdir

print(workdir)
fn = path+"/DD%04i/DD%04i" % (num, num) # parameter file to load

pf = load(fn) # load data
# This is the resolution we will extract at

DIMS = (pf.domain_dimensions)*(2**level)
print(level,DIMS)
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

FlatFileName = tempdir+'/%s_flatrho_%04i.fits' %(name, num)

pyfits.writeto(FlatFileName, cube["Density"])
pyfits.append(FlatFileName, cube["x-velocity"])
pyfits.append(FlatFileName, cube["y-velocity"])
pyfits.append(FlatFileName, cube["z-velocity"])

#arg = 'caleb_flatrho_0025.fits/1/10'
#


subprocess.call('cp '+ppdir+'/dustkappa_silicate.inp '+tempdir)
subprocess.call('cp '+ppdir+'/molecule_13co.inp '+tempdir)
subprocess.call('cp '+ppdir+'/camera_wavelength_micron.inp '+tempdir)

Problem_setup(FlatFileName, face = face, dust_temp = 10.0)
command = ppdir+'radmc3d image npix '+str(DIMS[0])+' iline 1 widthkms 10 linenlam 500 loadlambda fluxcons inclline linelist nostar writepop doppcatch sizepc 10'
subprocess.call(command)
makefits()


