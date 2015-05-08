#/usr/bin/env python
#import matplotlib
#matplotlib.use('pdf')
import yt
import numpy as np
from astropy.table import Table
import string
#from itertools import cycle

import glob
import os

# color=["red","green","blue","black","orange","cyan","gray"]
# lines = ["-","--",":"]
# linecycle = cycle(lines)
# colorcycle=cycle(color)

RootDir = '/lustre/home/eros/SimSuite8/'
DirList = glob.glob(RootDir+'Design*')+glob.glob(RootDir+'Fiducial*')
NameList = np.array(['DD0021','DD0022','DD0023','DD0024','DD0025','DD0026','DD0027','DD0028','DD0029','DD0030'])
NameList = NameList[::-1]

t = Table(names=('Name','Mass','SurfDens','SigmaV'),dtypes=('S80','f8','f8','f8'))

for ctr,DirName in enumerate(DirList):
    for Name in NameList:
        FileName = DirName+'/'+Name+'/'+Name
        if os.path.exists(FileName):
            t.add_row()
            print(FileName)
            ds = yt.load(FileName)
            ad = ds.all_data()
            mlist = ad.quantities.total_mass()
            try:
                mass = mlist[0]+mlist[1]
            except:
                mass = mlist
            mass = mass.convert_to_units('Msun')
            sigx,vx0 = ad.quantities.weighted_variance('x-velocity','Density')
            sigy,vy0 = ad.quantities.weighted_variance('y-velocity','Density')
            sigz,vz0 = ad.quantities.weighted_variance('z-velocity','Density')
            sig3d = ((sigx.convert_to_units('cm/s'))**2+\
                     (sigy.convert_to_units('cm/s'))**2+\
                     (sigz.convert_to_units('cm/s'))**2)**0.5
            t[-1]['Name'] = FileName
            t[-1]['Mass'] = mass.value
            t[-1]['SigmaV'] = sig3d.value/1e5
            t[-1]['SurfDens'] = mass.value / (100.)
t.write('/lustre/home/eros/simprops.fits')
