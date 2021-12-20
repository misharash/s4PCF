### combine_files.py (Michael Rashkovetskyi, adapted from Oliver Philcox, 2021)
# This reads in a set of (data-random) and (random) particle counts and uses them to construct the N-point functions, including edge-correction
# It is designed to be used with the run_npcf.csh script
# Currently 2PCF, 3PCF, 4PCF, 5PCF and 6PCF are supported. The code will try to load the edge correction coupling matrices from the coupling_matrices/ directory, and recompute them if not (using multithreading to speed up the 9j manipulations)
# This can handle both odd and even parity NPCFs.
# The output is saved to the working directory with the same format as the NPCF counts, with the filename ...zeta_{N}pcf.txt

import sys, os
import subprocess
import numpy as np
import multiprocessing
from sympy.physics.wigner import wigner_3j, wigner_9j

## First read-in the input file string from the command line
if len(sys.argv)!=3:
    raise Exception("Need to specify the input files and N_threads!")
else:
    inputs = str(sys.argv[1])
    threads = int(sys.argv[2])

print("Reading in files starting with %s\n"%inputs)

def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

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
    if N==2:
        countsR = np.loadtxt(R_file,skiprows=5)
    elif N==3:
        countsR = np.loadtxt(R_file,skiprows=6)
    elif N==4:
        countsR = np.loadtxt(R_file,skiprows=9)
    else:
        raise Exception("%dPCF not yet configured"%N)
        
    # Extract radial bins
    if N==2:
        bin1 = np.loadtxt(R_file,skiprows=4,max_rows=1)
    elif N==3:
        bin1,bin2 = np.loadtxt(R_file,skiprows=4,max_rows=2)
    elif N==4:
        bin1,bin2,bin3 = np.loadtxt(R_file,skiprows=6,max_rows=3)

    # Now load in D-R pieces and average
    countsN_all = []
    total_DmR = 0
    for i in range(100):
        DmR_file = inputs+'.n%s_%dpcf.txt'%(str(i).zfill(2),N)
        if not os.path.exists(DmR_file): continue
        # Extract counts
        if N==2:
            countsN_all.append(np.loadtxt(DmR_file,skiprows=5))
        elif N==3:
            countsN_all.append(np.loadtxt(DmR_file,skiprows=6)) # skipping rows with radial bins and ell
        elif N==4:
            countsN_all.append(np.loadtxt(DmR_file,skiprows=9))
    countsN_all = np.asarray(countsN_all)
    N_files = len(countsN_all)
    countsN = np.mean(countsN_all,axis=0)
    # could use this next line to compute std from finite number of randoms!
    #countsNsig = np.std(countsN_all,axis=0)/np.sqrt(N_files)

    # Now compute edge-correction equations

    if N==2:
        # isotropic 2PCF is easy!
        zeta = countsN/countsR

        # Now save the output to file, copying the first few lines from the N files
        zeta_file = inputs+'.zeta_%dpcf.txt'%N
        rfile = open(R_file,"r")
        zfile = open(zeta_file,"w")
        for l,line in enumerate(rfile):
            if l>=5: continue
            zfile.write(line)
        for a in range(len(zeta)):
            zfile.write("%.8e\t"%zeta[a])
        zfile.close()

    if N==3:
        # isotropic 3PCF should be easy too
        zeta = countsN/countsR

        # Now save the output to file, copying the first few lines from the N files
        zeta_file = inputs+'.zeta_%dpcf.txt'%N
        rfile = open(R_file,"r")
        zfile = open(zeta_file,"w")
        for l,line in enumerate(rfile):
            if l>=6: continue
            zfile.write(line)
        for a in range(len(zeta)):
            zfile.write("%.8e\t"%zeta[a])
        zfile.close()

    if N==4:
        # isotropic 4PCF should be easy as well
        zeta = countsN/countsR

        # Now save the output to file, copying the first few lines from the N files
        zeta_file = inputs+'.zeta_%dpcf.txt'%N
        rfile = open(R_file,"r")
        zfile = open(zeta_file,"w")
        for l,line in enumerate(rfile):
            if l>=9: continue
            zfile.write(line)
        for a in range(len(zeta)):
            zfile.write("%.8e\t"%zeta[a])
        zfile.close()

    print("Computed %dPCF using %d (data-random) files, saving to %s\n"%(N,N_files,zeta_file))
