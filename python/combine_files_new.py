### combine_files.py (Michael Rashkovetskyi, adapted from Oliver Philcox, 2021)
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

def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

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
    if N==2:
        bin1 = np.loadtxt(R_file,skiprows=4,max_rows=1)
    elif N==3:
        bin1,bin2 = np.loadtxt(R_file,skiprows=4,max_rows=2)
    elif N==4:
        bin1,bin2,bin3 = np.loadtxt(R_file,skiprows=6,max_rows=3)
    else:
        raise Exception("%dPCF not yet configured"%N)
    
    countsR_all = []
    for i in range(Nrandoms):
        R_file = inroot+'.r%d_%dpcf.txt'%(i, N)
        # Extract counts
        if N==2:
            countsR_all.append(np.loadtxt(R_file,skiprows=5))
        elif N==3:
            countsR_all.append(np.loadtxt(R_file,skiprows=6)) # skipping rows with radial bins and ell
        elif N==4:
            countsR_all.append(np.loadtxt(R_file,skiprows=9))
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
                if N==2:
                    countsN_all.append(np.loadtxt(DmR_file,skiprows=5))
                elif N==3:
                    countsN_all.append(np.loadtxt(DmR_file,skiprows=6)) # skipping rows with radial bins and ell
                elif N==4:
                    countsN_all.append(np.loadtxt(DmR_file,skiprows=9))
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

        if N==2:
            # isotropic 2PCF is easy!
            zeta = countsN/countsR

            # Now save the output to file, copying the first few lines from the N files
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
            rfile = open(R_file,"r")
            zfile = open(zeta_file,"w")
            for l,line in enumerate(rfile):
                if l>=9: continue
                zfile.write(line)
            for a in range(len(zeta)):
                zfile.write("%.8e\t"%zeta[a])
            zfile.close()

        if j<Ndata:
            print("Computed %dPCF using %d (random-random and data-random) files, saving to %s\n"%(N,Nrandoms,zeta_file))
        else:
            print("Averaged %dPCF using %d (data) files, saving to %s\n"%(N,Ndata,zeta_file))
