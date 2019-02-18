#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --job-name=radmc-%A-%a
#SBATCH --output=radmc-%A-%a.out
#SBATCH --array=1,3,45,6,7,8,11,12,14,15

export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE

module restore my_default

source $HOME/.bashrc
source $HOME/preload.bash

module load NiaEnv
module load gcc/7.3.0
module load openmpi/3.1.1
module load hdf5/1.10.2

export PPDIR=$HOME/code/simscript/postproc/
export DATADIR=/home/e/eros/eros/project/
export PPOUTDIR=$HOME/scratch/Orion_sims/
# export PYTHONPATH=$PYTHONPATH+':'+$PPDIR
# source $HOME/yt-x86_64/bin/activate

sim_name="Dsgn0$SLURM_ARRAY_TASK_ID"

cd $PPOUTDIR
echo "Current working directory is `pwd`"

if [ ! -d $sim_name ]; then
    mkdir $sim_name
fi

cd $sim_name

pids=

# Loop through the data files for each simulation since
# the timesteps are not homogeneous.
for data_file in $DATADIR/$sim_name/data.0*.hdf5; do
    $HOME/anaconda3/bin/python $HOME/code/simscript/postproc/pipeline_orion.py $PPOUTDIR/$sim_name $data_file 0 0 &
    pids+=" $!"
done

wait $pids || { echo "There was an error" >&2; exit 1; }

echo "All jobs exited."

# targetdir = sys.argv[1]
# timestep = float(sys.argv[2])
# face = float(sys.argv[3])
# level = float(sys.argv[4])

