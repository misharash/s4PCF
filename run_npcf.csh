#!/bin/csh

##################### DOCUMENTATION #####################
### Shell script for running the C++ s4PCF NPCF-estimator function on a data and data-random catalog, then combining the outputs, including edge-correction (Michael Rashkovetskyi, adapted from Oliver Philcox, 2021).
#
# This can be run either from the terminal or as a SLURM script (using the below parameters).
# The code should be compiled (with the relevant options, i.e. N-bins, ell-max and 4PCF/5PCF/6PCF) before this script is run. The isotropic 2PCF and 3PCF will always be computed.
# The script should be run from the code directory
# This is adapted from a similar script by Daniel Eisenstein.
# In the input directory, we expect compressed .gz files labelled {root}.data.gz, {ranroot}.ran.{IJ}.gz where {root} is a user-set name, and {IJ} indexes the random catalogs, from 0 - 31.
# We expect the summed weights to be the same for the data and each random catalog, but the random weights should be negative
# This script will compute the D^N counts, the (D-R)^N counts for 32 random subsets, and the R^N counts for one subset (should be sufficient).
# If the connected flag is set (and -DDISCONNECTED added to the makefile) we compute also the Gaussian contribution to the 4PCF.
# The output will be a set of .zeta_{N}pcf.txt files in the specified directory as well as a .tgz compressed directory of other intermediary outputs
# It is important to check the errlog file in the output directory to make sure all went well!
# Note that performing edge-correction is slow for the 5PCF and 6PCF since 9j symbols must be computed. Furthermore, the output multipoles are only accurate up to (ORDER-1), i.e. to compute an accurate edge-corrected spectrum with ell=5, we must compute (D-R) and R counts up to ell=6.
#
# NB: If needed, we could access a task ID by SLURM_ARRAY_TASK_ID, if we're running with SLURM
##########################################################

#SBATCH -n 16 # cpus
#SBATCH -N 1 # tasks
#SBATCH -t 0-02:59:59 # time
#SBATCH --mem-per-cpu=1GB
#SBATCH -o /home/ophilcox/out/boss4pcfSall_run.%A.out         # File to which STDOUT will be written (make sure the directory exists!)
#SBATCH -e /home/ophilcox/out/boss4pcfSall_run.%A.err         # File to which STDERR will be written
#SBATCH --mail-type=END,FAIL        # Type of email notification
#SBATCH --mail-user=ophilcox@princeton.edu # Email to which notifications will be sent

##################### INPUT PARAMETERS ###################

# Main inputs
set periodic = 0 # whether to run with periodic boundary conditions (should also be set in Makefile)
set rmin_short = 0 # minimum short radius in Mpc/h
set rmax_short = 30 # maximum short radius in Mpc/h
set rmin_long = 60 # minimum long radius in Mpc/h
set rmax_long = 240 # maximum long radius in Mpc/h
set rmin_cf = 0 # minimum 2PCF radius in Mpc/h
set rmax_cf = 300 # maximum 2PCF radius in Mpc/h

# Other inputs
set scale = 1 # rescaling for co-ordinates
set ngrid = 50 # grid-size for accelerating pair count
set boxsize = 1000 # only used if periodic=1

# File directories
set root = qpm_galaxies # root for data filenames
set ranroot = qpm_randoms # root for random filenames
set dir = `pwd`
set in = "$dir/in" # input directory (see above for required contents)
set out = "$dir/out" # output file directory
set tmp = "$dir/tmp" # temporary directory for intermediate file storage for this run (ideally somewhere with fast I/O)

# Load some python environment with numpy and sympy installed
# module load anaconda3
# conda activate ptenv

##########################################################

# Set number of threads (no SLURM)
set OMP_NUM_THREADS = 4

# Set number of threads (with SLURM)
# setenv OMP_NUM_THREADS $SLURM_NPROCS

# Define command to run the C++ code
set code = ./s4PCF

if ($periodic) then
  set command = "$code -rmax_short $rmax_short -rmin_short $rmin_short -rmax_long $rmax_long -rmin_long $rmin_long -rmax_cf $rmax_cf -rmin_cf $rmin_cf -ngrid $ngrid -scale $scale -boxsize $boxsize"
else
  set command = "$code -rmax_short $rmax_short -rmin_short $rmin_short -rmax_long $rmax_long -rmin_long $rmin_long -rmax_cf $rmax_cf -rmin_cf $rmin_cf -ngrid $ngrid -scale $scale"
endif

# Create a temporary directory for saving
/bin/rm -rf $tmp       # Delete, just in case we have crud from a previous run.
mkdir $tmp

# Copy this script in for posterity
cp run_npcf.csh $tmp

# Create output directory
if (!(-e $out)) then
    mkdir $out
endif

# Create an output file for errors
set errfile = errlog
set errlog = $out/$errfile
set tmpout = $tmp
rm -f $errlog
date > $errlog
echo Executing $0 >> $errlog
echo $command >> $errlog
echo $OMP_NUM_THREADS >> $errlog

# Filename for saved multipoles (a big file)
set multfile = $tmp/$root.mult

# Extra the data into our temporary ramdisk
gunzip -c $in/$root.data.gz > $tmp/$root.data

# Find number of galaxies (needed later for R^N periodic counts)
set Ngal = `cat $tmp/$root.data | wc -l`
set Ngal = `expr $Ngal + 1`

#### Compute D^N NPCF counts
# Note that we save the a_lm multipoles from the data here
echo Starting Computation
# echo "Starting D^N" >> $errlog
# date >> $errlog
# ($command -in $tmp/$root.data -save $multfile -outstr $root.data > $tmpout/$root.d.out) >>& $errlog
# save & load not implemented yet

# Remove the output - we don't use it
# rm output/$root.data_?pc*.txt

# echo "Done with D^N"

### Compute R^N NPCF counts
# We just use one R catalog for this and invert it such that the galaxies are positively weighted
gunzip -c $in/$ranroot.ran.00.gz > $tmp/$root.ran.00

echo "Starting R^N" >> $errlog
date >> $errlog
($command -in $tmp/$root.ran.00 -outstr $root.r -invert > $tmpout/$root.r.out) >>& $errlog
# Copy the output into the temporary directory
mv output/$root.r_*pc*.txt $tmpout/

echo "Done with R^N"

# Now make D-R for each of 49 random catalogs
foreach n ( 00 01 02 03 04 05 06 07 08 09 \
	    10 11 12 13 14 15 16 17 18 19 \
	    20 21 22 23 24 25 26 27 28 29 \
	    30 31 32 33 34 35 36 37 38 39 \
      40 41 42 43 44 45 46 47 48)

    # First copy the randoms and add the data
    /bin/cp -f $tmp/$root.data $tmp/$root.ran.$n
    gunzip -c $in/$ranroot.ran.$n.gz >> $tmp/$root.ran.$n

    ### Compute the (D-R)^N counts
    # This uses the loaded data multipoles from the D^N step
    # Note that we balance the weights here to ensure that Sum(D-R) = 0 exactly
    echo "Starting D-R $n" >> $errlog
    date >> $errlog
    # ($command -in $tmp/$root.ran.$n -load $multfile -outstr $root.n$n -balance > $tmpout/$root.n$n.out) >>& $errlog
    # save & load not implemented yet so instead
    ($command -in $tmp/$root.ran.$n -outstr $root.n$n -balance > $tmpout/$root.n$n.out) >>& $errlog
    # Copy the output into the temporary directory
    mv output/$root.n${n}_*pc*.txt $tmpout/

    # Remove the random catalog
    /bin/rm -f $tmp/$root.ran.$n
    echo Done with D-R $n

end    # foreach D-R loop

### Now need to combine the files to get the full NPCF estimate
# We do this in Python, and perform edge-correction unless the periodic flag is not set
if ($periodic) then
  echo Combining files together without performing edge-corrections (using analytic R^N counts)
  # this script not adapted yet
  python python/combine_files_periodic.py $tmpout/$root $Ngal $boxsize $rmin_short $rmax_short >>& $errlog
else
  echo Combining files together and performing edge-corrections
  python python/combine_files.py $tmpout/$root >>& $errlog
endif

# Do some cleanup
rm $tmp/$root.data # $multfile

# Now move the output files into the output directory.
# Compress all the auxilliary files and copy
echo Finished with computation.  Placing results into $out/
echo Finished with computation.  Placing results into $out/ >> $errlog
date >> $errlog
pushd $tmpout > /dev/null
echo >> $errlog
/bin/ls -l >> $errlog
/bin/cp $errlog .
tar cfz $root.tgz $root.*.out $root.*pc*.txt $errfile run_npcf.csh
popd > /dev/null
/bin/mv $tmpout/$root.tgz $tmpout/$root.zeta_*pcf.txt $out/

# Destroy ramdisk
/bin/rm -rf $tmp
