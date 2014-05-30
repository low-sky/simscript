import h5py
import numpy as np
import numpy.fft
from math import *
from optparse import OptionParser
import sys

######################################
# Generates random turbulent fields
# Modified from ORION2 version (SO)
####################################


def init_perturbations(n, kmin, kmax, dtype):
    kx = np.zeros(n, dtype=dtype)
    ky = np.zeros(n, dtype=dtype)
    kz = np.zeros(n, dtype=dtype)
    # perform fft k-ordering convention shifts
    for j in range(0,n[1]):
        for k in range(0,n[2]):
            kx[:,j,k] = n[0]*np.fft.fftfreq(n[0])
    for i in range(0,n[0]):
        for k in range(0,n[2]):
            ky[i,:,k] = n[1]*np.fft.fftfreq(n[1])
    for i in range(0,n[0]):
        for j in range(0,n[1]):
            kz[i,j,:] = n[2]*np.fft.fftfreq(n[2])
            
    kx = np.array(kx, dtype=dtype)
    ky = np.array(ky, dtype=dtype)
    kz = np.array(kz, dtype=dtype)
    k = np.sqrt(np.array(kx**2+ky**2+kz**2, dtype=dtype))

    # only use the positive frequencies
    inds = np.where(np.logical_and(k**2 >= kmin**2, k**2 < (kmax+1)**2))
    nr = len(inds[0])

    phasex = np.zeros(n, dtype=dtype)
    phasex[inds] = 2.*pi*np.random.uniform(size=nr)
    fx = np.zeros(n, dtype=dtype)
    fx[inds] = np.random.normal(size=nr)
    
    phasey = np.zeros(n, dtype=dtype)
    phasey[inds] = 2.*pi*np.random.uniform(size=nr)
    fy = np.zeros(n, dtype=dtype)
    fy[inds] = np.random.normal(size=nr)
    
    phasez = np.zeros(n, dtype=dtype)
    phasez[inds] = 2.*pi*np.random.uniform(size=nr)
    fz = np.zeros(n, dtype=dtype)
    fz[inds] = np.random.normal(size=nr)

    # rescale perturbation amplitude so that low number statistics
    # at low k do not throw off the desired power law scaling.
    for i in range(kmin, kmax+1):
        slice_inds = np.where(np.logical_and(k >= i, k < i+1))
        rescale = sqrt(np.sum(np.abs(fx[slice_inds])**2 + np.abs(fy[slice_inds])**2 + np.abs(fz[slice_inds])**2))
        fx[slice_inds] = fx[slice_inds]/rescale
        fy[slice_inds] = fy[slice_inds]/rescale
        fz[slice_inds] = fz[slice_inds]/rescale

    # set the power law behavior
    # wave number bins
    fx[inds] = fx[inds]*k[inds]**-(0.5*alpha)
    fy[inds] = fy[inds]*k[inds]**-(0.5*alpha)
    fz[inds] = fz[inds]*k[inds]**-(0.5*alpha)

    # add in phases
    fx = np.cos(phasex)*fx + 1j*np.sin(phasex)*fx
    fy = np.cos(phasey)*fy + 1j*np.sin(phasey)*fy
    fz = np.cos(phasez)*fz + 1j*np.sin(phasez)*fz

    return fx, fy, fz, kx, ky, kz


def normalize(fx, fy, fz):
    norm = np.sqrt(np.sum(fx**2 + fy**2 + fz**2)/np.product(n))
    fx = fx/norm
    fy = fy/norm
    fz = fz/norm
    return fx, fy, fz


def make_perturbations(n, kmin, kmax, f_solenoidal):
    fx, fy, fz, kx, ky, kz = init_perturbations(n, kmin, kmax, dtype)
    if f_solenoidal != None:
        k2 = kx**2+ky**2+kz**2
        # solenoidal part
        fxs = 0.; fys =0.; fzs = 0.
        if f_solenoidal != 0.0:
            fxs = fx - kx*(kx*fx+ky*fy+kz*fz)/np.maximum(k2,1e-16)
            fys = fy - ky*(kx*fx+ky*fy+kz*fz)/np.maximum(k2,1e-16)
            fzs = fz - kz*(kx*fx+ky*fy+kz*fz)/np.maximum(k2,1e-16)
            ind = np.where(k2 == 0)
            fxs[ind] = 0.; fys[ind] = 0.; fzs[ind] = 0.
        # compressive part
        # get a different random cube for the compressive part
        # so that we can target the RMS solenoidal fraction,
        # instead of setting a constant solenoidal fraction everywhere.
        fx, fy, fz, kx, ky, kz = init_perturbations(n, kmin, kmax, dtype)
        fxc = 0.; fyc =0.; fzc = 0.
        if f_solenoidal != 1.0:
            fxc = kx*(kx*fx+ky*fy+kz*fz)/np.maximum(k2,1e-16)
            fyc = ky*(kx*fx+ky*fy+kz*fz)/np.maximum(k2,1e-16)
            fzc = kz*(kx*fx+ky*fy+kz*fz)/np.maximum(k2,1e-16)
            ind = np.where(k2 == 0)
            fxc[ind] = 0.; fyc[ind] = 0.; fzc[ind] = 0.
    
        # back to real space
        pertx = np.real(np.fft.ifftn(f_solenoidal*fxs + (1.-f_solenoidal)*fxc))
        perty = np.real(np.fft.ifftn(f_solenoidal*fys + (1.-f_solenoidal)*fyc))
        pertz = np.real(np.fft.ifftn(f_solenoidal*fzs + (1.-f_solenoidal)*fzc))
    else:
        # just convert to real space
        pertx = np.real(np.fft.ifftn(fx))
        perty = np.real(np.fft.ifftn(fy))
        pertz = np.real(np.fft.ifftn(fz))
    
    # subtract off COM (assuming uniform density)
    pertx = pertx-np.average(pertx)
    perty = perty-np.average(perty)
    pertz = pertz-np.average(pertz)
    # scale RMS of perturbation cube to unity
    pertx, perty, pertz = normalize(pertx, perty, pertz)
    return pertx, perty, pertz


def get_erot_ke_ratio(pertx, perty, pertz):
    x, y, z = np.mgrid[0:n[0], 0:n[1], 0:n[2]]
    x = x - (n[0]-1)/2.
    y = y - (n[1]-1)/2.
    z = z - (n[2]-1)/2.
    r2 = x**2+y**2+z**2
    erot_ke_ratio = (np.sum(y*pertz-z*perty)**2 +
                     np.sum(z*pertx-x*pertz)**2 +
                     np.sum(x*perty-y*pertx)**2)/(np.sum(r2)*np.product(n)) 
    return erot_ke_ratio


def plot_spectrum1D(pertx, perty, pertz):
    # plot the 1D power to check the scaling.
    fx = np.abs(np.fft.fftn(pertx))
    fy = np.abs(np.fft.fftn(perty))
    fz = np.abs(np.fft.fftn(pertz))
    fx = np.abs(fx)
    fy = np.abs(fy)
    fz = np.abs(fz)
    kx = np.zeros(n, dtype=dtype)
    ky = np.zeros(n, dtype=dtype)
    kz = np.zeros(n, dtype=dtype)
    # perform fft k-ordering convention shifts
    for j in range(0,n[1]):
        for k in range(0,n[2]):
            kx[:,j,k] = n[0]*np.fft.fftfreq(n[0])
    for i in range(0,n[0]):
        for k in range(0,n[2]):
            ky[i,:,k] = n[1]*np.fft.fftfreq(n[1])
    for i in range(0,n[0]):
        for j in range(0,n[1]):
            kz[i,j,:] = n[2]*np.fft.fftfreq(n[2])
    k = np.sqrt(np.array(kx**2+ky**2+kz**2,dtype=dtype))
    k1d = []
    power = []
    for i in range(kmin,kmax+1):
        slice_inds = np.where(np.logical_and(k >= i, k < i+1))
        k1d.append(i+0.5)
        power.append(np.sum(fx[slice_inds]**2 + fy[slice_inds]**2 + fz[slice_inds]**2))
        print i,power[-1]
    import matplotlib.pyplot as plt
    plt.loglog(k1d, power)
    plt.show()


###################
# input parameters, read from command line
###################
parser = OptionParser()
parser.add_option('--kmin', dest='kmin', 
                  help='minimum wavenumber.', 
                  default=-1)
parser.add_option('--kmax', dest='kmax', 
                  help='maximum wavenumber.', 
                  default=-1)
parser.add_option('--size', dest='size', 
                  help='size of each direction of data cube.  default=256', 
                  default=256)
parser.add_option('--alpha', dest='alpha', 
                  help='negative of power law slope.  (Power ~ k^-alpha) '+
                  'supersonic turbulence is near alpha=2.  '+
                  'driving over a narrow band of two modes is often done with alpha=0', 
                  default = None)
parser.add_option('--seed', dest='seed', 
                  help='seed for random # generation.  default=0', 
                  default = 0)
parser.add_option('--f_solenoidal', dest='f_solenoidal', 
                  help='volume RMS fraction of solenoidal componet of the perturbations relative to the total.  ' + 
                       'If --f_solenoidal=None, the motions are purely random.  For low wave numbers ' +
                       'the relative imporance of solenoidal to compressive may be sensitve to the ' +
                       'choice of radom seed.  It has been suggested (Federrath 2008) that ' +
                       'f_solenoidal=2/3 is the most natural driving mode and this is currently' + 
                       'the suggested best-practice.', 
                  default = -1)

(options, args) = parser.parse_args()

# size of the data domain
n = [int(options.size), int(options.size), int(options.size)]
# range of perturbation length scale in units of the smallest side of the domain
kmin = int(options.kmin)
kmax = int(options.kmax)
print('Size: {0}, kmin: {1}, kmax:{2}'.format(n,kmin,kmax))
if kmin > kmax or kmin < 0 or kmax < 0:
    print "kmin must be < kmax, with kmin > 0, kmax > 0.  See --help."
    sys.exit(0)
if kmax > floor(np.min(n))/2:
    print "kmax must be <= floor(size/2).  See --help."
    sys.exit(0)
f_solenoidal = options.f_solenoidal
if f_solenoidal == "None" or f_solenoidal == "none":
    f_solenoidal = None
else:
    f_solenoidal = float(options.f_solenoidal)
    if f_solenoidal > 1. or f_solenoidal < 0.:
        print "You must choose f_solenoidal.  See --help."
        sys.exit(0)
alpha = options.alpha
if alpha==None:
    print "You must choose a power law slope, alpha.  See --help."
    sys.exit(0)
alpha = float(options.alpha)
if alpha < 0.:
    print "alpha is less than zero. Thats probably not what ou want.  See --help."
    sys.exit(0)
seed = int(options.seed)
# data precision
dtype = np.float64
# ratio of solenoidal to compressive components
if options.f_solenoidal=="None" or options.f_solenoidal==None:
    f_solenoidal = None
else:
    f_solenoidal = min(max(float(options.f_solenoidal), 0.), 1.)

###################
# begin computation
###################

np.random.seed(seed=seed)
pertx, perty, pertz = make_perturbations(n, kmin, kmax, f_solenoidal)
erot_ke_ratio = get_erot_ke_ratio(pertx, perty, pertz)
print "erot_ke_ratio = ", erot_ke_ratio

#plot_spectrum1D(pertx,perty,pertz)

#divV_rms = 0
#for i in range(1,n[0]-1):
#    for j in range(1,n[1]-1):
#        for k in range(1,n[2]-1):
#            divV_rms += ((pertx[i+1,j,k]-pertx[i-1,j,k])+(perty[i,j+1,k]-perty[i,j-1,k])+(pertz[i,j,k+1]-pertz[i,j,k-1]))**2
#divV_rms = sqrt(divV_rms/np.product(n))
#print divV_rms


# hdf5 output
f = h5py.File('RandomField1', 'w')
ds = f['/'].create_dataset('RandomField1', n, dtype=np.float)
ds[:] = pertx
ds.attrs['kmin'] = kmin
ds.attrs['kmax'] = kmax
ds.attrs['alpha'] = alpha
if f_solenoidal!=None: ds.attrs['f_solenoidal'] = f_solenoidal
ds.attrs['erot_ke_ratio'] = erot_ke_ratio
ds.attrs['seed'] = seed
ds.attrs['Component_Rank']=1 # Must be equal to Npart = 1
ds.attrs['Rank']=3
ds.attrs['Component_Size']=options.size 
ds.attrs['Dimensions'] = n
f.close()

f = h5py.File('RandomField2', 'w')
ds = f['/'].create_dataset('RandomField2', n, dtype=np.float)
ds[:] = perty
ds.attrs['kmin'] = kmin
ds.attrs['kmax'] = kmax
ds.attrs['alpha'] = alpha
if f_solenoidal!=None: ds.attrs['f_solenoidal'] = f_solenoidal
ds.attrs['erot_ke_ratio'] = erot_ke_ratio
ds.attrs['seed'] = seed
ds.attrs['Component_Rank']=1 # Must be equal to Npart = 1
ds.attrs['Rank']=3
ds.attrs['Component_Size']=options.size 
ds.attrs['Dimensions'] = n
f.close()

f = h5py.File('RandomField3', 'w')
ds = f['/'].create_dataset('RandomField3', n, dtype=np.float)
ds[:] = pertz
ds.attrs['kmin'] = kmin
ds.attrs['kmax'] = kmax
ds.attrs['alpha'] = alpha
if f_solenoidal!=None: ds.attrs['f_solenoidal'] = f_solenoidal
ds.attrs['erot_ke_ratio'] = erot_ke_ratio
ds.attrs['seed'] = seed
ds.attrs['Component_Rank']=1 # Must be equal to Npart = 1
ds.attrs['Rank']=3
ds.attrs['Component_Size']=options.size 
ds.attrs['Dimensions'] = n
f.close()
