#!/bin/env python
import os
import shutil
import sys

if sys.argv[1] is None:
    try:
        simdir = os.listdir('/global/scratch/eros/runs/SimSuite3/')
    except OSError:
        simdir = os.listdir('.')
else:
    simdir = sys.argv[1]

simnames = [name for name in os.listdir(simdir) if os.path.isdir(simdir+name)]
steplist = list((21,22,23,24,25,26,27,28,29,30))
facelist = list((0,1,2))
un = os.uname()
for nom in simnames:
    if sys.argv[2] is None:
        if os.access('/global/scratch/eros/runs/',os.W_OK):
            target_dir = '/global/scratch/eros/runs/'
        else:
            target_dir ='./'
    else:
        target_dir = sys.argv[2]
    target_name = 'Radmc_{0}.pbs'.format(nom)
    shutil.copy(os.getcwd()+'/templates/Radmc_template_'+un[1]+'.pbs',target_dir+target_name)
    with open(target_dir+target_name,'a') as template:
        for face in facelist:
            for step in steplist:
                template.write('python $PPDIR/pipeline.py '+simdir+'/'+nom+' {0} {1} 0\n'.format(step,face))
    template.close()




    
