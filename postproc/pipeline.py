import commands
import sys
import postproc as pp
import os
import shutil
targetdir = sys.argv[1]
timestep = float(sys.argv[2])
face = float(sys.argv[3])
level = float(sys.argv[4])

ppdir = '/home/eros/simscript/postproc/'
D = pp.FileSetup(targetdir,face=face,level=level,timestep=timestep,ppdir=ppdir)
os.chdir(D['TempDir'])
pp.ProblemSetup(D['FileName'], face = face, dust_temp = D['GasTemp'])

command = ppdir+'radmc3d image npix '+str(int(D['GridSize']))+' iline 1 widthkms 10 linenlam 500 loadlambda fluxcons inclline linelist nostar writepop doppcatch sizepc 10'
print(command)
result = commands.getoutput(command)
print(result)

pp.MakeFits(fitsfile=D['FileName']+'_radmc.fits',dpc = 260.0,toK=True)
shutil.move(D['FileName']+'_radmc.fits','/home/eros/runs/fitsfiles/.')

