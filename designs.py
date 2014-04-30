import numpy as np
import matplotlib.pyplot as plt
import os
import astropy.units as u
import astropy.constants as const
import subprocess
import shutil
from astropy.table import Table

values = np.loadtxt('gpRuns24N4D.txt')

rootdir = os.path.expanduser("~")+'/SimSuite/'
rootname = 'Design'

GenerateFields = False

# Domain Definition
# Fixed parameters for this simulation.

BoxSize = 10*u.pc
Tvals = 10.0*u.K*np.ones(values.shape[0])
kMin = 2*np.ones(values.shape[0]).astype('int')
kMax = 4*np.ones(values.shape[0]).astype('int')
RootGridSize = 128

# Design parameters
logVPmin = -1.0
logVPmax = +1.0

bmin = 0
bmax = 1

logbetamin = -1.0
logbetamax = 1.0

MachMin = 5
MachMax = 20

np.random.seed(8675309)
seeds = np.random.randint(long(2)**24,size=values.shape[0])

logVPvals = logVPmin+(logVPmax-logVPmin)*values[:,0]
logbetavals = logbetamin + (logbetamax-logbetamin)*values[:,2]

#fundamental design params
bvals = bmin + (bmax-bmin)*values[:,1]
MachVals = MachMin + (MachMax-MachMin)*values[:,3]
betavals = (1e1**logbetavals)
VPvals = (1e1**logVPvals)

#Derived parameters
SoundSpeed = ((const.k_B*Tvals/(2.33*const.m_n))**(0.5)).to(u.cm/u.s)
density = ((5*MachVals**2*SoundSpeed**2)/
           (6*const.G*VPvals*BoxSize**2)).to(u.g/u.cm**3)
tcross = (BoxSize/(MachVals*SoundSpeed)).to(u.s)
Bvals = ((8*np.pi*density*SoundSpeed**2/betavals)**(0.5)).value*(u.G)

params = Table([Tvals,bvals,Bvals,MachVals,kMin,kMax,seeds],\
                    names=('Kinetic Temperature','Solenoidal Fraction','Magnetic Field','Mach Number','kMin','kMax','Seed'))


params['Virial Parameter'] = VPvals
params['Plasma Beta'] = betavals
params['Index'] = np.arange(values.shape[0]).astype('int')
params['Sound Speed']=SoundSpeed.to(u.cm/u.s)
params['Crossing Time']=tcross.to(u.s)
params['PL Index']=2*np.ones(values.shape[0])
params['Box Size']=RootGridSize*np.ones(values.shape[0]).astype('int')
params['Density']=density

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
        template.write('Density = {0:6e}\n'.\
                           format(params[idx]['Density']))
        template.close()
    print(dirname)
params.write('parameter_table.ascii',format='ascii.fixed_width')
