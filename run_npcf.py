##################### DOCUMENTATION #####################
### Python script for running the C++ s4PCF NPCF-estimator function on a data and data-random catalog, then combining the outputs, including edge-correction (Michael Rashkovetskyi, adapted from Oliver Philcox, 2021).
#
# The code should be compiled (with the relevant options) before this script is run.
# The script should be run from the code directory
# This is adapted from a similar script by Daniel Eisenstein.
# In the input directory, we expect one or many plain text files {datafilenames} and one random catalog {randomfilename}; all in xyzw (4 columns, 3D coords and weight, space or tab-separated) format
# Data and/or random files can have extra columns (e.g. jackknife region numbers), those will be just ignored
# We do NOT expect the summed weights to be the same for the data and each random catalog, the random weights must NOT be negative
# This script will compute the D^N counts, the (D-R)^N counts for 32 random subsets, and the R^N counts for one subset (should be sufficient).
# The output will be a set of .zeta_{N}pcf.txt files in the specified directory as well as a .tgz compressed directory of other intermediary outputs
# It is important to check the errlog file in the output directory to make sure all went well!
##########################################################

import os
from datetime import datetime
import numpy as np

##################### INPUT PARAMETERS ###################

# Main inputs
periodic = 0 # whether to run with periodic boundary conditions (should also be set in Makefile)
rmin = 0 # minimum radius in Mpc/h
rmax = 30 # maximum radius in Mpc/h
rmin_long = 60 # minimum long side radius in Mpc/h
rmax_long = 120 # maximum long side radius in Mpc/h

# Other inputs
scale = 1 # rescaling for co-ordinates
ngrid = 50 # grid size for accelerating pair count
boxsize = 1000 # only used if periodic=1

# File names and directories
datafilenames = [f"qpm_galaxies/{i+1:04d}.xyz" for i in range(50)] # data filenames
randomfilename = "qpm_randoms_50x.xyzwj" # random filename
outroot = "qpm_galaxies" # base name for outputs
workdir = os.getcwd()
indir = os.path.join(workdir, "gaussian_s4PCF") # input directory (see above for required contents)
outdir = os.path.join(workdir, "out") # output file directory
tmpdir = os.path.join(workdir, "tmp") # temporary directory for intermediate file storage for this run (ideally somewhere with fast I/O)
scriptname = "run_npcf.py"

do_full = 1 # whether do full computation or only combine existing intermediate files

##########################################################

# Set OpenMP number of threads
OMP_NUM_THREADS = 30
os.environ["OMP_NUM_THREADS"] = str(OMP_NUM_THREADS)

# Define command to run the C++ code
code = "./s4PCF"

command = f"{code} -rmax {rmax} -rmin {rmin} -rmax_long {rmax_long} -rmin_long {rmin_long} -ngrid {ngrid} -scale {scale}"
if periodic:
  command += f" -boxsize {boxsize}"

# Create a temporary directory for saving
if do_full:
  os.system(f"rm -rf {tmpdir}") # Delete, just in case we have crud from a previous run.
os.makedirs(tmpdir, exist_ok=1)

# Copy this script in for posterity
os.system(f"cp {scriptname} {os.path.normpath(tmpdir)}")

# Create output directory
os.makedirs(outdir, exist_ok=1)

# Create an output file for errors
errfilename = "errlog"
errlog = os.path.join(outdir, errfilename)

def print_and_log(s):
  print(s)
  print_log(s)
print_log = lambda l: os.system(f"echo {l} >> {errlog}")

print_and_log(datetime.now())
print_log(f"Executing {__file__}")
print_log(command)
print_log(OMP_NUM_THREADS)

print("Starting Computation")

# Find number of lines (points) in each file
count_lines = lambda fname: sum(1 for _ in open(os.path.join(indir, fname)))

Ngals = [count_lines(fname) for fname in datafilenames]
print_and_log(f"Number of points in data file(s): {Ngals}")
Nrandoms = count_lines(randomfilename)
print_and_log(f"Number of points in random file: {Nrandoms}")

Ngal_avg = sum(Ngals)/len(Ngals)
Nparts = int(np.ceil(Nrandoms/Ngal_avg))
print_and_log(f"Divide randoms to {Nparts} part(s)")

if do_full:
  ### Divide randoms to desired number of parts, randomly
  # skip if only need to combine
  # first decide indices
  print("Creating random indices")
  random_indices = np.arange(Nrandoms)
  if Nparts>1:
    np.random.shuffle(random_indices)
  random_parts_indices = [random_indices[i::Nparts] for i in range(Nparts)]
  print("Created random indices")
  # now read contents
  print_and_log("Reading random file")
  print_and_log(datetime.now())
  random_contents = np.loadtxt(os.path.join(indir, randomfilename), usecols=range(4))
  # use only X, Y, Z coords and weights (4 first columns), consistently with data reading
  print_and_log("Read random file")
  random_contents[:, 3] *= -1 # negate the weights

### Now for each of random parts
for i in range(Nparts*do_full): # skip if do_full is false
  # Compute R^N NPCF counts
  print_and_log(f"Starting R[{i}]^N")
  print_and_log(datetime.now())
  outstr = f"{outroot}.r{i}"
  filename = os.path.join(tmpdir, outstr)
  random_content = random_contents[random_parts_indices[i]] # select part of randoms
  np.savetxt(filename, random_content) # save part to file
  # run code, forward output to separate file
  os.system(f"{command} -in {filename} -outstr {outstr} -invert >> {os.path.join(tmpdir, outstr)}.out")
  os.remove(filename) # clean up
  os.system(f"mv output/{outstr}_?pc*.txt {os.path.normpath(tmpdir)}/") # move output into the temporary dir
  print_and_log(f"Done with R[{i}]^N")

  ### Now for each data file
  for j, datafilename in enumerate(datafilenames):
    # Compute (D-R)^N NCPF counts
    print_and_log(f"Starting (D[{j}]-R[{i}])^N")
    print_and_log(datetime.now())
    outstr = f"{outroot}.{j}.n{i}"
    filename = os.path.join(tmpdir, outstr)
    data_content = np.loadtxt(os.path.join(indir, datafilename), usecols=range(4)) # read data
    # use only X, Y, Z coords and weights (4 first columns), consistently with randoms reading
    np.savetxt(filename, np.concatenate((data_content, random_content))) # save data and random part to file
    # run code, forward output to separate file
    os.system(f"{command} -in {filename} -outstr {outstr} -balance >> {os.path.join(tmpdir, outstr)}.out")
    os.remove(filename) # clean up
    os.system(f"mv output/{outstr}_?pc*.txt {os.path.normpath(tmpdir)}/") # move output into the temporary dir
    print_and_log(f"Done with (D[{j}]-R[{i}])^N")
  ### end for each data file
### end for each random part

### Now need to combine the files to get the full NPCF estimate
# We do this in another python script, and perform edge-correction unless the periodic flag is not set
if periodic:
  print("Combining files together without performing edge-corrections (using analytic R^N counts)")
  # The script doesn't exist yet
  os.system(f"python python/combine_files_periodic_new.py {os.path.join(tmpdir, outroot)} {len(datafilenames)} {Nparts} | tee -a {errlog}")
else:
  print("Combining files together and performing edge-corrections")
  # run script, print output to stdout AND append to errlog
  os.system(f"python python/combine_files_new.py {os.path.join(tmpdir, outroot)} {len(datafilenames)} {Nparts} | tee -a {errlog}")

print_and_log(f"Finished with computation. Placing results into {outdir}/")
print_log(datetime.now())
os.chdir(tmpdir)
os.system(f"ls -l >> {errlog}") # list tmpdir contents
# finished writing log, now copy it to tmpdir, we are in it now
os.system(f"cp {errlog} ./")
# Compress all the auxilliary files and copy
os.system(f"tar cfz {outroot}.tgz {outroot}.*.out {outroot}.*pc*.txt {errfilename} {scriptname}")
os.chdir(workdir)
# Now move the output files into the output directory.
os.system(f"mv {os.path.join(tmpdir, outroot)}.tgz {os.path.join(tmpdir, outroot)}.zeta_*pcf.txt {os.path.normpath(outdir)}/")

# Destroy temporary dir
os.system(f"rm -rf {tmpdir}")
