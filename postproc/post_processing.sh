#!/bin/bash
#source /home/cward/yt-x86_64/bin/activate
if [ $1 == "--help" ]; then
	
	echo "1) enzo data directory"
	echo "2) time step"
	echo "3) cube dimension e.g. 128"
	echo "4) is face number 0,1,2"
	#echo "5) Temp"
else
PPDIR=$HOME/runs/postproc/ # directory where post_processing script is stored
current_dir=`pwd`
input_dir=$1 # enzo data directory
t_s=$2  # time step
dim=$3 # cube dimension e.g. 128
face=$4 # $4 is face number 0,1,2
#T=$5  # Temp


# Time step formatting: (kind of arbitrary these few steps are to make exsiting scripts work)
len_t_s=${#t_s}
pro_set_zeros=$((4-$len_t_s)) # zeros for pading problem_setup.pro script input fits file
makefits_zeros=$((2-$len_t_s)) # zeros for pading makefits.pro script input fits file

zeros_p=`perl -e 'print "0" x '$pro_set_zeros` #zero padding for problem_setup
zeros_m=`perl -e 'print "0" x '$makefits_zeros` #zero padding for makefits

# Extract fits file name for later processing (e.g. if input directory is ~miayan/runs/enzo/run_128_5hd_12_20 then name is run_128_5hd_12_20 )

name=( `echo $input_dir | tr '/' ' '` )
name=${name[$((${#name[@]}-1))]}

# Extracttemp parameter from run name
T=( `echo $name | tr '_' ' '` )
T=${T[$((${#T[@]}-1))]}

# Create local self contained directory in order to avoid conflicts read/write conflicts during parallel execution of post processing

mkdir ${name}_$t_s
cd ${name}_$t_s

# Files needed to ran RADMC-3D

cp $PPDIR/dustkappa_silicate.inp .
cp $PPDIR/molecule_13co.inp .
cp $PPDIR/camera_wavelength_micron.inp .

# save_flatrho.py 
# -----------------
# input:
# 1) out_file_name 2) time_step 3) cube dimension 4) input directory where enzo data is found
#source /home/localad/yt/yt-x86_64/bin/activate
source $HOME/yt-x86_64/bin/activate
echo $t_s $dim $input_dir
python $PPDIR/save_flatrho.py ${name}_face$face $t_s $dim $input_dir

# problem_setup.pro / problem_setup.py
# -------------------------------------
# input:
# arg='name of the flatrho file produced by save_flatrho.py/face/simulation temperature'
 
arg="'${name}_face${face}_flatrho_$zeros_p$t_s.fits/$face/$T'"
python $PPDIR/problem_setup.py $arg
#idl -quiet -e "problem_setup,"$arg


<<COMMENT
Example radmc3d ran:
______________________
radmc3d image npix 255 loadlambda fluxcons inclline linelist nostar writepop doppcatch sizepc 10

doppcatch - turns on doppler catching (this interpolates between velocity jumps in the input data to create a smoother line); set in problem_setup
sizepc - the size of the output image
npix - number of pixels in the output image
nostar - no sources
writepop - outputs the population stats

COMMENT

#$PPDIR/radmc3d image npix $((dim-1)) iline 1 widthkms 10 linenlam 500 loadlambda fluxcons inclline linelist nostar writepop doppcatch sizepc 10

# makefits.pro / makefits.py
# ---------------------------
# input:
# 8) output name for final fits file
arg1="'$PPDIR/output/${name}_$zeros_m${t_s}_face$face.fits'"
python $PPDIR/makefits.py $arg1 dpc=260 #toK='toK'
#idl -quiet -e "makefits, fits ="$arg1", dpc = 260, /toK"

#cleanup
cd ..
rm -r ${name}_$t_s

fi
