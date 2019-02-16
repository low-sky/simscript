import commands
import sys
import postproc_yt as pp
import os
import shutil

targetdir = sys.argv[1]
timestep = float(sys.argv[2])
face = float(sys.argv[3])
level = float(sys.argv[4])

ppdir = os.getenv('PPDIR')
outdir = os.getenv('PPOUTDIR')

D = pp.FileSetup(targetdir, face=face, level=level,
                 timestep=timestep, ppdir=ppdir)

#pp.ProblemSetup(D['FileName'], face = face, dust_temp = D['GasTemp'])

os.chdir(D['TempDir'])

command = ppdir + 'radmc3d image npix ' + \
    str(int(D['GridSize'])) + \
    ' iline 1 widthkms 10 linenlam 500 loadlambda fluxcons inclline linelist nostar writepop doppcatch sizepc 10 norefine'
print(command)
result = commands.getoutput(command)
print(result)

save_name = os.path.join(outdir, D['FileName'][17:-5] + '_radmc.fits')

pp.MakeFits(fitsfile=save_name, dpc=260.0, toK=True)

shutil.move(save_name, outdir)
shutil.rmtree(D['TempDir'])
