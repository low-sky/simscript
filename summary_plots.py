#!/bin/env python
import os
import shutil
import sys

from yt.mods import *

stepnames = [name for name in os.listdir(os.getcwd()) if os.path.isdir(name)]
stepnames2 = [name for name in stepnames if 'DD002' in name]
stepnames3 = [name for name in stepnames if 'DD003' in name]
stepnames = stepnames2+stepnames3
designname = (os.getcwd().split('/'))[-1]
sides = ('x','y','z')
quants = ('Density','TotalEnergy')
dx = {'x':'DrivingField2','y':'DrivingField1','z':'DrivingField1'}
dy = {'x':'DrivingField3','y':'DrivingField3','z':'DrivingField2'}
if not os.path.isdir('Figures'):
    os.mkdir('Figures')
for sim in stepnames:
    pf = load(sim+'/'+sim)
    for face in sides:
        for q in quants:
            if q=='Density':
                p = ProjectionPlot(pf,face,q,weight_field=None,
                                   width=(10,'pc'))
                p.annotate_quiver(dx[face],dy[face],factor=24)
            else:
                p = ProjectionPlot(pf,face,q,weight_field=None,width=(10,'pc'))
                p.annotate_magnetic_field(factor=24)
            p.set_cmap(field=q,cmap='Greys')
            if pf.get_parameter('NumberOfParticles',type=int)>0:
                print('Found {0} particles'.format(pf.get_parameter('NumberOfParticles',type=int)))
                p.annotate_particles(4,p_size=24,col='blue')   
            p.annotate_title(designname+' '+sim+' '+face)
            p.save('Figures/'+designname+'.'+sim+'.'+q+'.'+face+'.png')




            
