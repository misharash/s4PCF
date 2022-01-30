### combine_files_new.py (Michael Rashkovetskyi, adapted from Oliver Philcox, 2022)
# This reads in a set of (data-random) and (random) particle counts and uses them to construct the N-point functions, including edge-correction
# It is designed to be used with the run_npcf.py script
# Currently 2PCF, 3PCF and 4PCF are supported.
# The output is saved to the working directory with the same format as the NPCF counts, with the filename ...zeta_{N}pcf.txt

import sys, os
import numpy as np

## First read-in the input file string from the command line
if len(sys.argv)!=4:
    raise Exception("Need to specify the input root (base name), number of data and of randoms!")
else:
    inroot = str(sys.argv[1])
    Ndata = int(sys.argv[2])
    Nrandoms = int(sys.argv[3])

# Decide which N we're using
Ns = []
for N in range(10):
    R_file = inroot+'.r0_%dpcf.txt'%N
    if os.path.exists(R_file):
        Ns.append(N)

if len(Ns)==0:
    raise Exception("No files found with input string %s.r0"%inroot)

for N in Ns:
    # First load in R piece
    R_file = inroot+'.r0_%dpcf.txt'%N
    # Extract radial bins
    if N not in [2, 3, 4]:
        raise Exception("%dPCF not yet configured" % N)
    
    countsR_all = []
    for i in range(Nrandoms):
        R_file = inroot+'.r%d_%dpcf.txt'%(i, N)
        # Extract counts
        countsR_all.append(np.loadtxt(R_file)[N-1:]) # skipping rows with radial bins
    countsR_all = np.asarray(countsR_all)
    countsR = np.mean(countsR_all,axis=0)

    R_file = inroot+'.r0_%dpcf.txt'%N

    countsN_alldata = []
    for j in range(Ndata+1): # extra iteration to do average among all data
        if j<Ndata:
            # Now load in D-R pieces and average
            countsN_all = []
            for i in range(Nrandoms):
                DmR_file = inroot+'.%d.n%d_%dpcf.txt'%(j, i, N)
                # Extract counts
                countsN_all.append(np.loadtxt(DmR_file)[N-1:]) # skipping rows with radial bins
            countsN_all = np.asarray(countsN_all)
            countsN = np.mean(countsN_all,axis=0)
            countsN_alldata.append(countsN)
        else: # do average among all data in the last iteration
            countsN_alldata = np.asarray(countsN_alldata)
            countsN = np.mean(countsN_alldata,axis=0)
        
        # set output file name - format does not depend on N
        zeta_file = inroot+'.%d.zeta_%dpcf.txt'%(j, N)
        if j==Ndata: # if we average over data
            zeta_file = inroot+'.zeta_%dpcf.txt'%N

        # Now compute edge-correction equations

        # isotropic 2/3/4PCF are all easy and similar
        zeta = countsN/countsR

        # Now save the output to file, copying the first few lines from the N files
        rfile = open(R_file,"r")
        zfile = open(zeta_file,"w")
        # copy comments and bins from first random file
        lnc = 0 # counter of lines that are not comments
        for line in rfile:
            if lnc >= N-1: continue
            if line[0] != "#": lnc += 1
            zfile.write(line)
        rfile.close()
        # write NPCF
        for a in range(len(zeta)):
            for b in range(len(zeta[0])):
                zfile.write("%.8e\t" % zeta[a, b])
            zfile.write("\n")
        zfile.close()

        if j<Ndata:
            print("Computed %dPCF using %d (random-random and data-random) files, saving to %s\n"%(N,Nrandoms,zeta_file))
        else:
            print("Averaged %dPCF using %d (data) files, saving to %s\n"%(N,Ndata,zeta_file))
