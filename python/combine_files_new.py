### combine_files_new.py (Michael Rashkovetskyi, adapted from Oliver Philcox, 2022)
# Handles both periodic and aperiodic cases
# In aperiodic case, this reads in a set of (data-random) and (random) particle counts and uses them to construct the N-point functions, including edge-correction
# In periodic case, this reads in only a set of (data-random) particle counts, computes analytic random counts and uses them to construct the N-point functions, without edge-correction
# It is designed to be used with the run_npcf.py script
# Currently fine2PCF, 2PCF, 3PCF and 4PCF are supported.
# The output is saved to the working directory with the same format as the NPCF counts, with the filename ...zeta_{N}pcf.txt

import sys, os
import numpy as np

## First read-in the input file string from the command line
if len(sys.argv)!=5 and len(sys.argv)!=12:
    raise Exception("Need to specify the periodicity (0 or 1), input root (base name), number of data and of randoms!")
else:
    periodic = int(sys.argv[1])
    inroot = str(sys.argv[2])
    Ndata = int(sys.argv[3])
    Nrandoms = int(sys.argv[4])

if not periodic and len(sys.argv)!=5:
    raise Exception("In non-periodic case, only need to specify the periodicity (0), input root (base name), number of data and of randoms!")

if periodic and len(sys.argv)!=12:
    raise Exception("In periodic case, need to specify the periodicity (1), input root (base name), number of data, number of randoms, box size, minimal and maximal short radii, min and max long radii, min and max fine2PCF radii!")
else:
    boxsize = float(sys.argv[5])
    rmin_short = float(sys.argv[6])
    rmax_short = float(sys.argv[7])
    rmin_long = float(sys.argv[8])
    rmax_long = float(sys.argv[9])
    rmin_cf = float(sys.argv[10])
    rmax_cf = float(sys.argv[11])
    # number of galaxies is not needed, since counts are normalized by (sum of positive weights)^-N
    radius = lambda bin, maxbin, rmin, rmax: rmin + 1.*bin*(rmax-rmin)/(maxbin+1)
    bin_volume = lambda bin, maxbin, rmin, rmax: boxsize**-3*4.*np.pi/3.*(radius(bin+1,maxbin,rmin,rmax)**3.-radius(bin,maxbin,rmin,rmax)**3.)
    bin_volume_short = lambda bin, maxbin: bin_volume(bin, maxbin, rmin_short, rmax_short)
    bin_volume_long = lambda bin, maxbin: bin_volume(bin, maxbin, rmin_long, rmax_long)
    bin_volume_cf = lambda bin, maxbin: bin_volume(bin, maxbin, rmin_cf, rmax_cf)

# Decide which N we're using
Ns = []
for N in ["2", "3", "4", "fine2"]:
    DmR_file = inroot+'.0.n0_%spcf.txt'%N # DmR file exists in any case
    if os.path.exists(DmR_file):
        Ns.append(N)

if len(Ns)==0:
    raise Exception("No files found with input string %s.0.n0"%inroot)

for N in Ns:
    n = int(N[-1]) # same as N for 2,3,4; 2 for fine2
    
    if periodic:
        DmR_file = inroot+'.0.n0_%spcf.txt'%N
        bins = np.loadtxt(DmR_file, max_rows=n-1)
        if N == "fine2":
            bins = np.arange(len(bins))[:, None] # in this case "bins" is midpoints, need to override with bin numbers and add second dim
            n_mu = len(np.loadtxt(DmR_file, skiprows=1, max_rows=1)) # count mu bins
            countsR = 0.5 * bin_volume_cf(bins, np.max(bins)) / n_mu # factor of 1/2 because we count only i > j pairs for fine2pcf
        elif N == "2":
            countsR = bin_volume_short(bins, np.max(bins))
        elif N == "3":
            countsR = bin_volume_short(bins[0], np.max(bins)) * bin_volume_short(bins[1], np.max(bins))
        elif N == "4":
            countsR = bin_volume_long(bins[0], np.max(bins[0])) * bin_volume_short(bins[1], np.max(bins[1:])) * bin_volume_short(bins[2], np.max(bins[1:]))
        else: # should never reach this
            raise Exception(f"Unrecoginized N: {N}")
    else:
        # First load in R piece
        countsR_all = []
        for i in range(Nrandoms):
            R_file = inroot+'.r%d_%spcf.txt'%(i, N)
            # Extract counts
            countsR_all.append(np.loadtxt(R_file, skiprows=(len(N)>1))[n-1:]) # skipping rows with radial bins, skip 1 more row for fine2
        countsR_all = np.asarray(countsR_all)
        countsR = np.mean(countsR_all,axis=0)

    R_file = inroot+'.0.n0_%spcf.txt'%N # technically it's DmR_file, but that will be overridden. R_file won't exist in periodic case. This will be used as source of bin data for output.

    countsN_alldata = []
    for j in range(Ndata+1): # extra iteration to do average among all data
        if j<Ndata:
            # Now load in D-R pieces and average
            countsN_all = []
            for i in range(Nrandoms):
                DmR_file = inroot+'.%d.n%d_%spcf.txt'%(j, i, N)
                # Extract counts
                countsN_all.append(np.loadtxt(DmR_file, skiprows=(len(N)>1))[n-1:]) # skipping rows with radial bins, skip 1 more row for fine2
            countsN_all = np.asarray(countsN_all)
            countsN = np.mean(countsN_all,axis=0)
            countsN_alldata.append(countsN)
        else: # do average among all data in the last iteration
            countsN_alldata = np.asarray(countsN_alldata)
            countsN = np.mean(countsN_alldata,axis=0)
        
        # set output file name - format does not depend on N
        zeta_file = inroot+'.%d.zeta_%spcf.txt'%(j, N)
        if j==Ndata: # if we average over data
            zeta_file = inroot+'.zeta_%spcf.txt'%N

        # Now compute edge-correction equations

        # isotropic 2/3/4PCF are all easy and similar
        zeta = countsN/countsR

        # Now save the output to file, copying the first few lines from the N files
        rfile = open(R_file,"r")
        zfile = open(zeta_file,"w")
        # copy comments and bins from first random file
        lnc = 0 # counter of lines that are not comments
        for line in rfile:
            if lnc >= n-1+(len(N)>1): break # only need N-1 data lines (2 for fine2PCF), i.e. radial bins, can terminate loop afterwards
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
            if periodic:
                print("Computed %sPCF using %d (data-random) files, saving to %s\n"%(N,Nrandoms,zeta_file))
            else:
                print("Computed %sPCF using %d (random-random and data-random) files, saving to %s\n"%(N,Nrandoms,zeta_file))
        else:
            print("Averaged %sPCF using %d (data) files, saving to %s\n"%(N,Ndata,zeta_file))
