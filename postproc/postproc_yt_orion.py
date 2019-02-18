import yt

try:
    # New yt analysis package
    from yt_astro_analysis.radmc3d_export.api import RadMC3DWriter
except ImportError:
    from yt.analysis_modules.radmc3d_export.api import RadMC3DWriter

import astropy.io.fits as pyfits
import subprocess
import os
import numpy as np
import readRADMC as RADMC


def FileSetup(targetdir, data_file,
              face=0, level=0,
              ppdir='/home/eros/code/simscript/postproc/',
              Tk=10):

    # homedir = os.getenv('HOME')
    if not os.access(targetdir, os.R_OK):
        workdir = os.getenv('PWD')
    else:
        workdir = ''

#    execfile(ppdir+'problem_setup.py')
#    execfile(ppdir+'makefits.py')
#    execfile(ppdir+'readRADMC.py')

    # name = os.path.join(targetdir,
    name = '{0}{1}_{2}'.format(data_file.split("/")[-1].rstrip('hdf5'),
                               face, level)

    tempdir = os.path.join(workdir, 'tmp_' + name + '/')
    try:
        print(subprocess.call('rm -rf ' + tempdir))
    except FileNotFoundError:
        print("No previous temp files found.")
    os.mkdir(tempdir)

    print('Building files with prefix %s' % (name))
    # num = timestep
    # path = os.path.join(workdir, targetdir)

    # fn = path + "/DD%04i/DD%04i" % (num, num)  # parameter file to load

    # initfile = open(fn, 'r')
    # inittxt = initfile.readlines()
    # initfile.close()
    # EOSSound = float(
    #     (([line for line in inittxt if "EOSSoundSpeed" in line][0]).split())[2])
    # Mu = float((([line for line in inittxt if "Mu " in line][0]).split())[2])
    # Tk = EOSSound**2 * Mu * 1.67e-24 / 1.38e-16

    print(ppdir)
    print(tempdir)

    DIMS = YTexport(data_file, tempdir)
    WriteAuxFiles(tempdir)
    command = 'cp ' + ppdir + '/dustkappa_silicate.inp ' + tempdir
    print(command)
    subprocess.call(command)
    command = 'cp ' + ppdir + '/molecule_13co.inp ' + tempdir
    subprocess.call(command)
    command = 'cp ' + ppdir + '/camera_wavelength_micron.inp ' + tempdir
    command = 'cp ' + ppdir + '/wavelength_micron.inp ' + tempdir
    subprocess.call(command)

    FlatFileName = tempdir + '/{}_flatrho.fits'.format(name)

    OutputDict = {'GasTemp': Tk,
                  'TempDir': tempdir,
                  'FileName': FlatFileName,
                  'GridSize': DIMS[0]}

    return OutputDict


def WriteAuxFiles(TempDir, nPhot=10000, lines_mode=3,
                  maxJbariter=1000, lvgconvergeto=0.01):
    temp_file = open(TempDir + 'dustopac.inp', 'w')
    temp_file.write("2               Format number of this file \n")
    temp_file.write("1               Nr of dust species\n")
    temp_file.write(
        "============================================================================\n")
    temp_file.write("1               Way in which this dust species is read\n")
    temp_file.write("0               0=Thermal grain\n")
    temp_file.write(
        "silicate        Extension of name of dustkappa_***.inp file\n")
    temp_file.write(
        "----------------------------------------------------------------------------\n")
    temp_file.close()

    #
    # Write the lines.inp control file
    #

    temp_file = open(TempDir + 'lines.inp', 'w')
    temp_file.write('2 \n')
    temp_file.write('1 \n')
    temp_file.write('13co    leiden    0    0  1 \n')
    temp_file.write('h2 \n')
    temp_file.close()

    #
    # Write the radmc.inp control file
    #

    temp_file = open(TempDir + 'radmc.inp', 'w')
    temp_file.write('nphot = %s \n' % nPhot)
    temp_file.write('scattering_mode_max = 0\n')

    temp_file.write('lines_mode = %s\n' % lines_mode)
    temp_file.write('lines_nonlte_maxiter = %s\n' % maxJbariter)
    temp_file.close()


def YTexport(filename, TempDir):

    origdir = os.getenv('PWD')
    os.chdir(TempDir)

    x_co = 5.695e-5 / 77.0
    mu_h = yt.YTQuantity(4.00e-24, 'g')

    def _NumberDensityCO(field, data):
        return (x_co / mu_h) * data["density"]

    def _NumberDensityH2(field, data):
        return (1 / mu_h) * data["density"]
    dust_to_gas = 0.01

    def _DustDensity(field, data):
        return dust_to_gas * data["density"]

    def _MicroTurbulence(field, data):
        return (3 * yt.physical_constants.boltzmann_constant_cgs *
                data["temperature"] / mu_h)**0.5
                # put in 0.01 to see if Doppler is working

    def _vx(field, data):
        return (data["x-velocity"])

    def _vy(field, data):
        return (data["y-velocity"])

    def _vz(field, data):
        return (data["z-velocity"])

    # Need to define a temperature field
    def _temp(field, data):
        return (data['ones'] * 10 * yt.units.K)

    ds = yt.load(filename)

    ds.add_field(("gas", "dust_density"),
                 function=_DustDensity, units="g/cm**3")
    ds.add_field(("gas", "number_density_CO"),
                 function=_NumberDensityCO, units="cm**-3")
    ds.add_field(("gas", "number_density"), function=_NumberDensityH2,
                 units="cm**-3", force_override=True)
    ds.add_field(("gas", "temperature"), function=_temp,
                 units="K", force_override=True)
    ds.add_field(("gas", "microturbulence"),
                 function=_MicroTurbulence, units="cm/s", force_override=True)
    ds.add_field(("gas", "x-velocity-cgs"), function=_vx,
                 units="cm/s", force_override=True)
    ds.add_field(("gas", "y-velocity-cgs"), function=_vy,
                 units="cm/s", force_override=True)
    ds.add_field(("gas", "z-velocity-cgs"), function=_vz,
                 units="cm/s", force_override=True)

    writer = RadMC3DWriter(ds, 0)  # Use only root grid extraction
    writer.write_amr_grid()
    writer.write_dust_file(("gas", "dust_density"), "dust_density.inp")
    writer.write_line_file(("gas", "number_density_CO"), "numberdens_13co.inp")
    writer.write_line_file(("gas", "number_density"), "numberdens_h2.inp")
    writer.write_line_file(("gas", "microturbulence"), "microturbulence.inp")
    writer.write_line_file(("gas", "temperature"), "gas_temperature.inp")
    writer.write_dust_file(("gas", "temperature"), "dust_temperature.dat")
    velocity_fields = ["x-velocity-cgs", "y-velocity-cgs", "z-velocity-cgs"]
    writer.write_line_file(velocity_fields, "gas_velocity.inp")

    # filelist = ["dust_density.inp","numberdens_13co.inp",
    #             "numberdens_h2.inp", "microturbulence.inp",
    #             "gas_temperature.inp", "dust_temperature.dat",
    #             "gas_velocity.inp","amr_grid.inp"]
    # for thisfile in filelist:
    #     command = 'mv '+thisfile+' '+TempDir+'/'
    #     subprocess.call(command)

    os.chdir(origdir)
    return(ds.domain_dimensions)


def MakeFits(fitsfile='image.fits', dpc=None,
             mom0fitsfile='mom0.fits', lambda0=2600.7576, toK=True):

    Boltz = 1.380700e-16
    ckms = 2.9979e5
    a = RADMC.readimage(filename='image.out')

    lambdanum = a.nrfr
    xnum = a.nx
    ynum = a.ny
    xsize = a.sizepix_x
    ysize = a.sizepix_y

    a.image[a.image == 0.0] = np.nan

    #   ind = 1
    #   dpc = 1

    x1 = a.x[0]
    y1 = a.y[0]

    firstlambda = a.wavelength[0]

    if (lambdanum > 1):
        deltalam = float(a.wavelength[1]) - float(a.wavelength[0])
        v0 = ckms * (1 - lambda0 / a.wavelength[0])
        v1 = ckms * (1 - lambda0 / a.wavelength[1])
        delta = v1 - v0
    else:
        deltalam = a.wavelength[0]

    if toK:
        # Uses a single wavelength for the sake of memory
        scalefac = (0.5 * (firstlambda * 1e-4)**2) / Boltz
    else:
        scalefac = 1

    hdu = pyfits.PrimaryHDU(data=(scalefac * a.image))
    hd = hdu.header

    if dpc is not None:
        if (dpc != 0):
            hd['CDELT1'] = -xsize / (dpc * 3.086e18) * (2 * np.pi / 360)
            hd['CRPIX1'] = xnum / 2
            hd['CRVAL1'] = 180e0
            hd['CTYPE1'] = 'RA---CAR'
            hd['CDELT2'] = ysize / (dpc * 3.086e18) * (2 * np.pi / 360)
            hd['CRPIX2'] = ynum / 2
            hd['CRVAL2'] = 0.0
            hd['CTYPE2'] = 'DEC--CAR'
            hd['EQUINOX'] = 2000.0
            hd['radesys'] = 'FK5'
            hd['TELESCOP'] = 'RADMC3D'
            hd['object'] = 'SIM@ %.0f PC' % dpc
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
        if toK:
            lam0 = 2.99792458e8 / 110.2013542798e9 * 1e6
            vaxis = 2.99792458e8 * (lam0 - a.wavelength) / lam0
            dv = np.median(vaxis - np.roll(vaxis, 1))
            ind = np.argmin(abs(vaxis))
            hd['CDELT3'] = dv
            hd['CRPIX3'] = ind + 1
            hd['CRVAL3'] = vaxis[ind + 1]
            hd['CTYPE3'] = 'VOPT'
            hd.update('BUNIT', 'K', 'Tmb')
        else:
            scalefac = 1
            hd['CDELT3'] = deltalam
            hd['CRPIX3'] = 1.0
            hd['CRVAL3'] = firstlambda
            hd['CTYPE3'] = 'Wavelength - micron'
            hd['BUNIT'] = 'erg/s/cm^2/Hz/ster'
            #       print(type(a.image))
    print("Writing Output to %s" % fitsfile)
    hdu.writeto(fitsfile, clobber=True)
