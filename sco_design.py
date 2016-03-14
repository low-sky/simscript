import numpy as np
import matplotlib.pyplot as plt
import os
import astropy.units as u
import astropy.constants as const
import subprocess
import shutil
from astropy.table import Table

values = np.loadtxt('design_51run_5fac.txt')

MinLogMass = 2
MaxLogMass = 4
Mass = 1e1**(values[:,0]*(MaxLogMass-MinLogMass)+MinLogMass)

# Design parameters
logVPmin = np.log10(0.3)
logVPmax = np.log10(3)
logVPvals = logVPmin+(logVPmax-logVPmin)*values[:,1]
VPvals = 1e1**logVPvals

logZmin = np.log10(0.003)
logZmax = np.log10(2)
logZvals = logZmin + (logZmax-logZmin)*values[:,2]
Zvals = 1e1**(logZvals)

logGmin = -1
logGmax = 3
logGvals = logGmin + (logGmax-logGmin)*values[:,3]

lognmin = 1.5
lognmax = 4
lognvals = lognmin + (lognmax-lognmin)*values[:,4]
nvals = 1e1**lognvals

t = Table([Mass,VPvals,Zvals,logGvals,nvals],names=('Mass','VirialParameter','Z','log10_G0','VolumeDensity'))
t.write('sco_design.csv')

