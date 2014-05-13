#!/bin/env python
import os
import shutil

simnames = os.listdir('/global/scratch/eros/runs/SimSuite3/')
steplist = list((21,22,23,24,25,26,27,28,29,30))
facelist = list((0,1,2))
for nom in simnames:
    target_dir = '/global/scratch/eros/runs/'
    target_name = 'Radmc_{0}.pbs'.format(nom)
    shutil.copy(os.getcwd()+'/templates/Radmc_template.pbs',target_dir+target_name)
    with open(target_dir+target_name,'a') as template:
        for face in facelist:
            for step in steplist:
                template.write('python $PPDIR/pipeline.py SimSuite3/'+nom+' {0} {1} 0\n'.format(step,face))
    template.close()




    
