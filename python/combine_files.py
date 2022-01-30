### combine_files.py (Michael Rashkovetskyi, adapted from Oliver Philcox, 2021-2022)
# This reads in a set of (data-random) and (random) particle counts and uses them to construct the N-point functions, including edge-correction
# It is designed to be used with the run_npcf.csh script
# Currently 2PCF, 3PCF and 4PCF are supported.
# The output is saved to the working directory with the same format as the NPCF counts, with the filename ...zeta_{N}pcf.txt

import sys, os
import numpy as np

## First read-in the input file string from the command line
if len(sys.argv)!=2:
    raise Exception("Need to specify the input files!")
else:
    inputs = str(sys.argv[1])

print("Reading in files starting with %s\n"%inputs)

# Decide which N we're using
Ns = []
for N in range(10):
    R_file = inputs+'.r_%dpcf.txt'%N
    if os.path.exists(R_file):
        Ns.append(N)

if len(Ns)==0:
    raise Exception("No files found with input string %s"%inputs)

for N in Ns:
    # First load in R piece
    R_file = inputs+'.r_%dpcf.txt'%N
    if N not in [2, 3, 4]:
        raise Exception("%dPCF not yet configured" % N)

    # Now load in D-R pieces and average
    countsN_all = []
    for i in range(100):
        DmR_file = inputs+'.n%s_%dpcf.txt'%(str(i).zfill(2),N)
        if not os.path.exists(DmR_file): continue
        # Extract counts
        countsN_all.append(np.loadtxt(DmR_file)[N-1:]) # skipping rows with radial bins
    countsN_all = np.asarray(countsN_all)
    N_files = len(countsN_all)
    countsN = np.mean(countsN_all,axis=0)
    # could use this next line to compute std from finite number of randoms!
    #countsNsig = np.std(countsN_all,axis=0)/np.sqrt(N_files)

    # Now compute edge-correction equations

    # isotropic 2/3/4PCF are all easy and similar
    zeta = countsN/countsR

    # Now save the output to file, copying the first few lines from the N files
    zeta_file = inputs+'.zeta_%dpcf.txt'%N
    rfile = open(R_file,"r")
    zfile = open(zeta_file,"w")
    # copy comments and bins from first random file
    lnc = 0 # counter of lines that are not comments
    for line in rfile:
        if lnc >= N-1: break # only need N-1 data lines, i.e. radial bins, can terminate loop afterwards
        if line[0] != "#": lnc += 1
        zfile.write(line)
    rfile.close()
    # write NPCF
    for a in range(len(zeta)):
        for b in range(len(zeta[0])):
            zfile.write("%.8e\t" % zeta[a, b])
        zfile.write("\n")
    zfile.close()

    print("Computed %dPCF using %d (data-random) files, saving to %s\n"%(N,N_files,zeta_file))
