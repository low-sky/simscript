import numpy as np
import matplotlib.pyplot as plt
import os
import astropy.units as u
import astropy.constants as const
import subprocess
import shutil
from astropy.table import Table


# Start with the seed
np.random.seed(649810806)


rootdir = os.path.expanduser("~")+'/SimSuite10/'
rootname = 'Design'
fid_rootname = 'Fiducial'
GenerateFields = True
AppendPermutations = 1
Fiducials = True
NFiducials = 5

# Load Design; permute if requested.
values = np.loadtxt('xd5n25.txt')
if AppendPermutations>0:
     originals = np.copy(values)
     for counter in np.arange(AppendPermutations):
          newvalues = np.copy(originals)
          np.random.shuffle(newvalues.T)
          values = np.concatenate((values,newvalues),axis=0)
          

# Domain Definition
# Fixed parameters for this simulation.

BoxSize = 10*u.pc
Tvals = 10.0*u.K*np.ones(values.shape[0])
kMin = (np.ceil(4*values[:,4])).astype(np.int)
# Trap the bottom case.
kMin[kMin==0]=1
kMax = 2*kMin

RootGridSize = 256

# Design parameters
logVPmin = np.log10(0.5)
logVPmax = logVPmin+np.log10(10)

bmin = 0.25
bmax = 0.75

logbetamin = -0.3
logbetamax = 1.0

MachMin = 5
MachMax = 15

seeds = np.random.randint(long(2)**24,size=values.shape[0])

logVPvals = logVPmin+(logVPmax-logVPmin)*values[:,0]
logbetavals = logbetamin + (logbetamax-logbetamin)*values[:,2]

#fundamental design params
bvals = bmin + (bmax-bmin)*values[:,1]
MachVals = MachMin + (MachMax-MachMin)*values[:,3]
betavals = (1e1**logbetavals)
VPvals = (1e1**logVPvals)
kminVals = kMin
kmaxVals = kMax
print(kminVals)
print(kmaxVals)

#Derived parameters
SoundSpeed = ((const.k_B*Tvals/(2.33*const.m_n))**(0.5)).to(u.cm/u.s)
density = ((5*MachVals**2*SoundSpeed**2)/
           (2*const.G*VPvals*BoxSize**2)).to(u.g/u.cm**3)
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
params['kMin'] = kminVals.astype(np.int)
params['kMax'] = kmaxVals.astype(np.int)

for idx,bval in enumerate(bvals):
     dirname = rootdir+rootname+str(params['Index'][idx]).zfill(2)+'/'
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
         template.write('RandomForcingMachNumber = {0:4f}\n'.\
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
params.write('d10_parameter_table.ascii',format='ascii.fixed_width')
params.write('d10_parameter_table.csv',format='ascii',delimiter=',')

if Fiducials:
    fidval = 0.5*np.ones(NFiducials)
    seeds = np.random.randint(long(2)**24,size=(NFiducials))
    logVPvals = logVPmin+(logVPmax-logVPmin)*fidval
    logbetavals = logbetamin + (logbetamax-logbetamin)*fidval
    Tvals = 10*np.ones(NFiducials)*u.K
    kMin = 2*np.ones(NFiducials).astype('int')
    kMax = 4*np.ones(NFiducials).astype('int')

    bvals = bmin + (bmax-bmin)*fidval
    MachVals = MachMin + (MachMax-MachMin)*fidval
    betavals = (1e1**logbetavals)
    VPvals = (1e1**logVPvals)

    SoundSpeed = ((const.k_B*Tvals/(2.33*const.m_n))**(0.5)).to(u.cm/u.s)
    density = ((5*MachVals**2*SoundSpeed**2)/
               (6*const.G*VPvals*BoxSize**2)).to(u.g/u.cm**3)
    tcross = (BoxSize/(MachVals*SoundSpeed)).to(u.s)
    Bvals = ((8*np.pi*density*SoundSpeed**2/betavals)**(0.5)).value*(u.G)
    kminVals = 2*np.ones(len(betavals))
    kmaxVals = 4*np.ones(len(betavals))

    fparams = Table([Tvals,bvals,Bvals,MachVals,kMin,kMax,seeds],\
                   names=('Kinetic Temperature','Solenoidal Fraction',
                          'Magnetic Field','Mach Number',
                          'kMin','kMax','Seed'))
    

    fparams['Virial Parameter'] = VPvals
    fparams['Plasma Beta'] = betavals
    fparams['Index'] = np.arange(NFiducials).astype('int')
    fparams['Sound Speed']=SoundSpeed.to(u.cm/u.s)
    fparams['Crossing Time']=tcross.to(u.s)
    fparams['PL Index']=2*np.ones(NFiducials)
    fparams['Box Size']=RootGridSize*np.ones(NFiducials).astype('int')
    fparams['Density']=density
    fparams['kMin'] = kminVals
    fparams['kMax'] = kmaxVals

    for idx,bval in enumerate(bvals):

        dirname = rootdir+fid_rootname+str(fparams['Index'][idx]).zfill(2)+'/'
        if not os.path.isdir(dirname):
            os.makedirs(dirname)
            # Generate a driving field into the RandomField1,2,3 etc.
            if GenerateFields:
                callstring = "perturbation_enzo.py --size={5} --alpha={4} --kmin={0} --kmax={1} --f_solenoidal={2:.3f} --seed={3}".\
                             format(fparams[idx]['kMin'],
                                    fparams[idx]['kMax'],
                                    fparams[idx]['Solenoidal Fraction'],
                                    fparams[idx]['Seed'],
                                    fparams[idx]['PL Index'],
                                    fparams[idx]['Box Size'])
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
                               format(fparams[idx]['Box Size'],
                                      fparams[idx]['Box Size'],
                                      fparams[idx]['Box Size']))
                template.write('InitialBfield = {0:6e}\n'.\
                               format(fparams[idx]['Magnetic Field']/\
                                      np.sqrt(4*np.pi)))
                # This 4pi is for CTMHD units
                template.write('MachNumber = {0:4f}\n'.\
                               format(fparams[idx]['Mach Number']))
                template.write('RandomForcingMachNumber = {0:4f}\n'.\
                               format(params[idx]['Mach Number']))
                template.write('RandomSeed = {0}\n'.\
                               format(fparams[idx]['Seed']))
                template.write('EOSSoundSpeed = {0:6e}\n'.\
                               format(fparams[idx]['Sound Speed']))
                template.write('SoundVelocity = {0:6e}\n'.\
                               format(fparams[idx]['Sound Speed']))
                template.write('TimeUnits = {0:6e}\n'.\
                               format(fparams[idx]['Crossing Time']))
                template.write('Density = {0:6e}\n'.\
                               format(fparams[idx]['Density']))
                template.close()
                print(dirname)
    fparams.write('d10_fiducial_table.ascii',format='ascii.fixed_width')
    fparams.write('d10_fiducial_table.csv',format='ascii',delimiter=',')
