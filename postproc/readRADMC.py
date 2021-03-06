#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#             PYTHON PACKAGE TO READ AND PLOT RESULTS FROM RADMC-3D
#            Based on the included readradmc.pro
#           Written by: Kevin Sooley
#                 July 22, 2012
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.cm as cm
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#                        ROUTINES FOR IMAGES
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
#                     READ THE RECTANGULAR TELESCOPE IMAGE
#---------------------------------------------------------------------------


class Object(object):
    pass

def readimage(filename=None, imagefile=None):
    if filename is None:
        filename = 'image.out'
    f = open(filename)
    lines = f.readlines()
    # iformat = int(lines[0])
    nx, ny = lines[1].split()
    nf = int(lines[2])
    sizepix_x, sizepix_y = lines[3].split()
    sizepix_x = float(sizepix_x)
    sizepix_y = float(sizepix_y)
    wavelength = np.zeros(int(nf))
    for i in range(len(wavelength)):
        wavelength[i] = lines[i + 4]

    image_lines = lines[(len(wavelength)+4):]
    image = np.array([np.float64(line.strip()) for line in image_lines
                      if not line.isspace()])
    image = image.reshape(int(nf), int(nx), int(ny)) 
    print(image.shape)
    xi = ((np.arange(int(nx)) + 0.5) / (int(nx)) - 0.5) * sizepix_x * int(nx)
    yi = ((np.arange(int(ny)) + 0.5) / (int(ny)) - 0.5) * sizepix_y * int(ny)

    # Compute the flux
    flux = np.sum(image)
    pc = 3.0857200e+18
    flux = flux / pc**2

    returnVar = Object()
    setattr(returnVar, 'nx', int(nx))
    setattr(returnVar, 'ny', int(ny))
    setattr(returnVar, 'nrfr', int(nf))
    setattr(returnVar, 'sizepix_x', sizepix_x)
    setattr(returnVar, 'sizepix_y', sizepix_y)
    setattr(returnVar, 'image', image)
    setattr(returnVar, 'flux', flux)
    setattr(returnVar, 'x', xi)
    setattr(returnVar, 'y', yi)
    setattr(returnVar, 'wavelength', wavelength)
    print('return')
    return returnVar


#---------------------------------------------------------------------------
#               read and plot the spectrum
#---------------------------------------------------------------------------

def B_lambda(T, wavelength):
    h = 6.62606957e-27
    c = 299792458
    kb = 1.3806488e-16
    return (2.0 * h * c**2 / wavelength**5) * (1.0 / (np.exp(h * c / (wavelength * kb * T)) - 1.0))


def readspectrum(filename=None, imagefile=None):
    if filename is None:
        filename = 'spectrum.out'
    if imagefile is None:
        imagefile = 'spectrum.png'
    f = open(filename)
    lines = f.readlines()
    nlam = int(lines[1])
    wavelength = np.zeros(nlam)
    flux = np.zeros(nlam)
    for i in range(nlam):
        wavelength[i], flux[i] = lines[i + 3].split()
    wave = np.logspace(-8, 1, num=3000)
    plt.xlim(0.3e-6, 1e-6)
    # plt.ylim(10e-3,2)
    plt.plot(wavelength * 10**-6, flux / flux[0], 'r-')
    plt.plot(wave, B_lambda(18000, wave) / B_lambda(18000, 0.3e-6), 'b-')
    plt.xlabel(r'$\lambda (\mu m)$')
    plt.ylabel(r'Flux')
    plt.savefig(imagefile)
#---------------------------------------------------------------------------
#   Let's get fancy and call radmc3d from inside python
#---------------------------------------------------------------------------
# if __name__ == "__main__":
    #readimage('image_U.out', 'image_U.png')
    #readimage('image_B.out', 'image_B.png')
    #readimage('image_V.out', 'image_V.png')
    #readimage('image_R.out', 'image_R.png')
    #readimage('image_I.out', 'image_I.png')
    # readspectrum()
#   readimage()
