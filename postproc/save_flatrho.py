import sys
from yt.mods import *
import astropy.io.fits as pyfits
import string

# we will use pyfits to write to our output file.
#import h5py

name = sys.argv[1]
num = float(sys.argv[2])
path=sys.argv[4]



fn = path+"/DD%04i/DD%04i" % (num, num) # parameter file to load

print(fn)

pf = load(fn) # load data

# This is the resolution we will extract at
DIMS = int(sys.argv[3])
# Now, we construct an object that describes the data region and structure we
# want
cube = pf.h.covering_grid(0, # The level we are willing to extract to; higher
                          left_edge=[0.0, 0.0, 0.0],
                          right_edge=[1.0, 1.0, 1.0],
                          dims=[DIMS,DIMS,DIMS],
                          fields=["Density"])

#print cube

pyfits.writeto('%s_flatrho_%04i.fits' %(name, num), cube["Density"])
pyfits.append('%s_flatrho_%04i.fits' %(name, num), cube["x-velocity"])
pyfits.append('%s_flatrho_%04i.fits' %(name, num), cube["y-velocity"])
pyfits.append('%s_flatrho_%04i.fits' %(name, num), cube["z-velocity"])
