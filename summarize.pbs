#!/bin/bash

#PBS -S /bin/bash
#PBS -l pmem=2000mb
#PBS -l procs=1
#PBS -l walltime=24:00:00
#PBS -m bea
#PBS -M erosolow@ualberta.ca
#PBS -l epilogue=/home/eros/bin/epilogue.sh

export EHOME=/lustre/home/eros/
export PPDIR=$EHOME/code/simscript/postproc/
export PPOUTDIR=/lustre/home/eros/fitsfiles/
export PYTHONPATH=$PYTHONPATH+':'+$PPDIR
source $EHOME/yt-x86_64/bin/activate
cd $PBS_O_WORKDIR
echo "Current working directory is `pwd`"

cd /lustre/home/eros/SimSuite8/Design0
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design1
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design2
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design3
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design4
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design5
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design6
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design7
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design8
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design9
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design10
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design11
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design12
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design13
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design14
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design15
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design16
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design17
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design18
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design19
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design20
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design21
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design22
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design23
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design24
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design25
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design26
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design27
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design28
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design29
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design30
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Design31
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Fiducial0
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Fiducial1
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Fiducial2
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Fiducial3
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/Fiducial4
python $EHOME/code/simscript/summary_plots.py
cd /lustre/home/eros/SimSuite8/
tar cvfz FIGURES.tar.gz `find . | grep Figures`

