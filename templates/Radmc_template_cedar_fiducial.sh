#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=512000M
#SBATCH --job-name=radmc-%A-%a
#SBATCH --output=radmc-%A-%a.out
#SBATCH --array=0-5

export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE

module restore my_default

source $HOME/.bashrc
source $HOME/preload.bash

module load gcc/7.3.0
module load openmpi/2.1.1
module load hdf5/1.8.20

export PPDIR=$HOME/code/simscript/postproc/
export DATADIR=$HOME/scratch/Orion_sims/
export PPOUTDIR=$HOME/scratch/Orion_sims/radmc_13co
# export PYTHONPATH=$PYTHONPATH+':'+$PPDIR
# source $HOME/yt-x86_64/bin/activate

sim_name="Fiducial0$SLURM_ARRAY_TASK_ID"

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
    $HOME/anaconda3/bin/python  $HOME/code/simscript/postproc/pipeline_orion.py $PPOUTDIR/$sim_name $data_file 0 0 &
    pids+=" $!"
done

wait $pids || { echo "There was an error" >&2; exit 1; }

echo "All jobs exited."

