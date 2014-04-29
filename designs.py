import numpy as np
import matplotlib.pyplot as plt
import os
import astropy.units as u
import astropy.constants as const
import subprocess
import shutil
from astropy.table import Table

values = np.loadtxt('gpRuns24N4D.txt')

rootdir = './SimSuite/'
rootname = 'Design'

GenerateFields = False

BoxSize = 10*u.pc

# Design parameters
Tmin = 8
Tmax = 30

bmin = 0
bmax = 1

logBmin = -7.0
logBmax = -5.0

MachMin = 5
MachMax = 20
np.random.seed(8675309)
seeds = np.random.randint(long(2)**24,size=values.shape[0])

# Fixed parameters for this simulation. 
kMin = 2*np.ones(values.shape[0]).astype('int')
kMax = 4*np.ones(values.shape[0]).astype('int')
RootGridSize = 128

Tvals = (Tmin + (Tmax-Tmin)*values[:,0])*u.K
bvals = bmin + (bmax-bmin)*values[:,1]
logBvals = logBmin + (logBmax-logBmin)*values[:,2]
MachVals = MachMin + (MachMax-MachMin)*values[:,3]
Bvals = (1e1**logBvals)*u.G


params = Table([Tvals,bvals,Bvals,MachVals,kMin,kMax,seeds],\
                    names=('Kinetic Temperature','Solenoidal Fraction','Magnetic Field','Mach Number','kMin','kMax','Seed'))

soundspeed = (const.k_B*Tvals/(2.33*const.m_n))**(0.5)
tcross = (BoxSize/(MachVals*soundspeed)).to(u.s)

params['Index'] = np.arange(values.shape[0]).astype('int')
params['Sound Speed']=soundspeed.to(u.cm/u.s)
params['Crossing Time']=tcross.to(u.s)
params['PL Index']=2*np.ones(values.shape[0])
params['Box Size']=RootGridSize*np.ones(values.shape[0]).astype('int')


for idx,bval in enumerate(bvals):

    dirname = rootdir+rootname+str(params['Index'][idx])+'/'
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    # Generate a driving field into the RandomField1,2,3 etc.
    if GenerateFields:
        callstring = "perturbation_enzo.py --size={5} --alpha={4} --kmin={0} --kmax={1} --f_solenoidal={2:.3f} --seed={3}".\
            format(params[idx]['kMin'],
                   params[idx]['kMax'],
                   params[idx]['Solenoidal Fraction'],
                   params[idx]['Seed'],
                   params[idx]['PL Index'],
                   params[idx]['Box Size'])
        callstring = 'python '+os.getcwd()+'/'+callstring
        print(callstring)
        subprocess.call(callstring,shell=True)
        shutil.move(os.getcwd()+'/RandomField1',dirname)
        shutil.move(os.getcwd()+'/RandomField2',dirname)
        shutil.move(os.getcwd()+'/RandomField3',dirname)

    shutil.copy(os.getcwd()+'/templates/MHDstart.pbs',dirname)
    shutil.copy(os.getcwd()+'/templates/MHDrestart.pbs',dirname)
    shutil.copy(os.getcwd()+'/templates/DrivenTurbulenceCTMHD',dirname)
    with open(dirname+'DrivenTurbulenceCTMHD','a') as template:
        template.write('TopGridDimensions = {0} {1} {2}\n'.\
                           format(params[idx]['Box Size'],
                                  params[idx]['Box Size'],
                                  params[idx]['Box Size']))
        template.write('InitialBfield = {0:6e}\n'.\
                           format(params[idx]['Magnetic Field']/\
                                      np.sqrt(4*np.pi)))
                       # This 4pi is for CTMHD units
        template.write('MachNumber = {0:4f}\n'.\
                           format(params[idx]['Mach Number']))
        template.write('RandomSeed = {0}\n'.\
                           format(params[idx]['Seed']))
        template.write('EOSSoundSpeed = {0:6e}\n'.\
                           format(params[idx]['Sound Speed']))
        template.write('SoundVelocity = {0:6e}\n'.\
                           format(params[idx]['Sound Speed']))
        template.write('TimeUnits = {0:6e}\n'.\
                           format(params[idx]['Crossing Time']))
        template.close()
    print(dirname)
params.write('parameter_table.ascii',format='ascii.fixed_width')
